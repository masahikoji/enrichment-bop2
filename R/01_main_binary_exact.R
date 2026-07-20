# Globally calibrated enrichment BOP2 design
# Exact calibration and exact operating characteristics for the globally
# calibrated enrichment BOP2 design with post-enrichment biomarker-positive
# monitoring. No Monte Carlo simulation is performed. Computationally
# intensive exact-enumeration steps are parallelized on macOS/Linux.
#
# The core analysis uses base R. CSV files are always produced. A consolidated
# XLSX workbook is also produced when either the 'writexl' or 'openxlsx' package
# is installed. Set environment variable BOP2_QUICK=1 for a short code check.

rm(list = ls())
options(stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
## 1. User settings
## -----------------------------------------------------------------------------

p0 <- 0.20
p1 <- 0.40
alpha <- 0.10
interim_n <- c(20L, 30L)
N_A <- 40L
N_pos <- 40L
n_pos_min <- 10L
post_enrichment_n <- c(20L, 30L)

prior_a <- 0.5
prior_b <- 0.5
prevalence_grid <- c(0.4, 0.5, 0.6, 0.7, 0.8)

theta_pos_grid <- c(p1, p1 + 0.10, p1 + 0.20)
theta_neg_grid <- c(p0, p0 - 0.10)

lambda_grid <- seq(0.005, 1.000, by = 0.005)

# gamma is a positive shape parameter rather than a probability and therefore
# need not be bounded by one. Values above one allow a cutoff that remains
# relatively close to one early in the trial and decreases more sharply later.
# The upper limit of eight is a numerical search bound; a warning is issued
# below if the selected design reaches this limit.
gamma_grid <- seq(0.050, 8.000, by = 0.050)

# Repository root and output directory. run_main.R and run_all.R set
# ENRICHMENT_BOP2_ROOT automatically. Direct execution from the repository root
# is also supported.
project_dir <- path.expand(Sys.getenv(
  "ENRICHMENT_BOP2_ROOT",
  unset = normalizePath(getwd(), mustWork = TRUE)
))
result_root <- file.path(project_dir, "results", "main")

if (!dir.exists(file.path(project_dir, "R"))) {
  stop(
    "Repository root was not found: ", project_dir,
    "\nRun this analysis from the repository root or set ENRICHMENT_BOP2_ROOT."
  )
}
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

# Parallel processing is used for the two computationally intensive,
# embarrassingly parallel steps: construction of biomarker-positive kernels
# across candidate boundaries and construction of all-comer state spaces across
# operating-characteristic scenarios. On macOS/Linux, forked workers are used.
# Windows falls back to sequential execution because this program is intended
# primarily for the specified macOS project environment.
#
# By default, use one fewer than the number of detected physical cores, capped
# at 8 to limit memory use. Override this choice, for example with 10 workers,
# before starting R by setting:
#   Sys.setenv(BOP2_N_CORES = 10)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- parallel::detectCores(logical = TRUE)
}
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
default_n_cores <- max(1L, min(8L, detected_cores - 1L))
requested_n_cores <- suppressWarnings(
  as.integer(Sys.getenv("BOP2_N_CORES", unset = ""))
)
n_cores <- if (!is.na(requested_n_cores) && requested_n_cores >= 1L) {
  min(requested_n_cores, detected_cores)
} else {
  default_n_cores
}
output_dir <- result_root

quick_test <- identical(Sys.getenv("BOP2_QUICK"), "1")
if (quick_test) {
  message("BOP2_QUICK=1: running a reduced code check, not the manuscript analysis.")
  lambda_grid <- seq(0.05, 1.00, by = 0.05)
  gamma_grid <- seq(0.25, 2.00, by = 0.25)
  prevalence_grid <- c(0.4, 0.8)
  theta_pos_grid <- c(0.4, 0.6)
  theta_neg_grid <- c(0.2, 0.1)
  n_cores <- 1L
  output_dir <- file.path(result_root, "quick_test")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

parallel_backend <- if (n_cores > 1L && .Platform$OS.type != "windows") {
  "parallel::mclapply (forked)"
} else {
  "sequential lapply"
}
message(
  "Parallel configuration: detected physical cores = ", detected_cores,
  "; workers used = ", n_cores,
  "; backend = ", parallel_backend
)

# Remove legacy Monte Carlo outputs so that old files are not mistaken for
# results from the current exact-only analysis.
legacy_output_files <- c(
  "monte_carlo_operating_characteristics.csv",
  "exact_vs_monte_carlo.csv",
  "simulation_scenarios.csv",
  "enrichment_bop2_results.xlsx"
)
unlink(file.path(output_dir, legacy_output_files), force = TRUE)

stopifnot(
  length(interim_n) == 2L,
  all(interim_n > 0L),
  all(diff(interim_n) > 0L),
  max(interim_n) < min(N_A, N_pos),
  n_pos_min > 0L,
  n_pos_min < N_pos,
  all(post_enrichment_n > 0L),
  all(diff(post_enrichment_n) > 0L),
  all(post_enrichment_n < N_pos),
  all(theta_pos_grid >= 0 & theta_pos_grid <= 1),
  all(theta_neg_grid >= 0 & theta_neg_grid <= 1)
)

## -----------------------------------------------------------------------------
## 2. Utility functions
## -----------------------------------------------------------------------------

parallel_lapply <- function(X, FUN, ..., n_workers = n_cores) {
  n_workers <- max(1L, min(as.integer(n_workers), length(X)))
  if (length(X) == 0L) {
    return(list())
  }
  if (n_workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      X,
      FUN,
      ...,
      mc.cores = n_workers,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(X, FUN, ...)
  }
}

write_results_workbook <- function(sheets, path) {
  # CSV files are the primary reproducible outputs. This helper additionally
  # combines the main results into one XLSX workbook when a supported package
  # is available. It never installs packages automatically.
  sheets <- lapply(sheets, function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(x) <- NULL
    x
  })

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path = path)
    return(TRUE)
  }

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    for (sheet_name in names(sheets)) {
      openxlsx::addWorksheet(wb, sheetName = sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = sheets[[sheet_name]])
      openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
      openxlsx::setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(sheets[[sheet_name]])), widths = "auto")
    }
    openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
    return(TRUE)
  }

  warning(
    "Neither 'writexl' nor 'openxlsx' is installed. CSV files were created, ",
    "but the consolidated XLSX workbook was not. Install one of these packages ",
    "and rerun the program to create the XLSX file."
  )
  FALSE
}

prob_binom_greater <- function(cutoff, size, prob) {
  # P{Binomial(size, prob) > cutoff}; vectorized in cutoff.
  out <- numeric(length(cutoff))
  out[cutoff < 0] <- 1
  out[cutoff >= size] <- 0
  idx <- cutoff >= 0 & cutoff < size
  if (any(idx)) {
    out[idx] <- pbinom(cutoff[idx], size = size, prob = prob,
                       lower.tail = FALSE)
  }
  out
}

posterior_futility_boundary <- function(sample_size, max_sample_size,
                                        lambda, gamma,
                                        p_ref = p0,
                                        a = prior_a, b = prior_b) {
  cutoff <- 1 - lambda * (sample_size / max_sample_size)^gamma
  x <- 0:sample_size
  post_prob <- pbeta(p_ref, shape1 = a + x,
                     shape2 = b + sample_size - x)
  eligible <- which(post_prob > cutoff)
  if (length(eligible) == 0L) -1L else max(eligible) - 1L
}

make_boundary_candidates <- function(sample_sizes, max_sample_size,
                                     lambda_values, gamma_values,
                                     prefix) {
  grid <- expand.grid(
    lambda = lambda_values,
    gamma = gamma_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  boundary_matrix <- matrix(
    0L,
    nrow = nrow(grid),
    ncol = length(sample_sizes),
    dimnames = list(NULL, paste0("b_", sample_sizes))
  )

  for (j in seq_along(sample_sizes)) {
    s <- sample_sizes[j]
    x <- 0:s
    post_prob <- pbeta(p0, shape1 = prior_a + x,
                       shape2 = prior_b + s - x)
    cutoffs <- 1 - grid$lambda * (s / max_sample_size)^grid$gamma

    # The posterior probability decreases monotonically with x, so the number
    # of values exceeding the cutoff minus one is the integer futility boundary.
    boundary_matrix[, j] <- vapply(
      cutoffs,
      function(cc) as.integer(sum(post_prob > cc) - 1L),
      integer(1L)
    )
  }

  keep <- !duplicated(boundary_matrix)
  out <- cbind(
    data.frame(
      candidate_id = seq_len(sum(keep)),
      component = prefix,
      lambda = grid$lambda[keep],
      gamma = grid$gamma[keep],
      stringsAsFactors = FALSE
    ),
    as.data.frame(boundary_matrix[keep, , drop = FALSE])
  )
  rownames(out) <- NULL
  out
}

make_joint_distribution <- function(n, pi, theta_pos, theta_neg) {
  # Distribution of (n_pos, x_pos, x_neg) after n all-comer patients.
  rows <- vector("list", (n + 1L) * (n + 2L) * (n + 3L) / 6L)
  idx <- 1L
  for (n_pos in 0:n) {
    p_n_pos <- dbinom(n_pos, size = n, prob = pi)
    for (x_pos in 0:n_pos) {
      p_x_pos <- dbinom(x_pos, size = n_pos, prob = theta_pos)
      n_neg <- n - n_pos
      for (x_neg in 0:n_neg) {
        rows[[idx]] <- c(
          n_pos = n_pos,
          x_pos = x_pos,
          x_neg = x_neg,
          prob = p_n_pos * p_x_pos *
            dbinom(x_neg, size = n_neg, prob = theta_neg)
        )
        idx <- idx + 1L
      }
    }
  }
  out <- as.data.frame(do.call(rbind, rows[seq_len(idx - 1L)]))
  out$n_pos <- as.integer(out$n_pos)
  out$x_pos <- as.integer(out$x_pos)
  out$x_neg <- as.integer(out$x_neg)
  out$x_all <- out$x_pos + out$x_neg
  out
}

convolve_distributions <- function(dist1, dist2, total_n) {
  # Convolution of two independent distributions of (n_pos, x_pos, x_neg).
  dim1 <- total_n + 1L
  dim2 <- total_n + 1L
  out_prob <- numeric(dim1 * dim2 * dim2)

  for (i in seq_len(nrow(dist1))) {
    n_pos <- dist1$n_pos[i] + dist2$n_pos
    x_pos <- dist1$x_pos[i] + dist2$x_pos
    x_neg <- dist1$x_neg[i] + dist2$x_neg

    # R arrays are column-major. All indices below are zero-based before +1.
    linear_index <- (n_pos + 1L) +
      dim1 * x_pos +
      dim1 * dim2 * x_neg

    out_prob[linear_index] <- out_prob[linear_index] +
      dist1$prob[i] * dist2$prob
  }

  nonzero <- which(out_prob > 0)
  idx0 <- nonzero - 1L
  n_pos <- idx0 %% dim1
  tmp <- idx0 %/% dim1
  x_pos <- tmp %% dim2
  x_neg <- tmp %/% dim2

  data.frame(
    n_pos = as.integer(n_pos),
    x_pos = as.integer(x_pos),
    x_neg = as.integer(x_neg),
    x_all = as.integer(x_pos + x_neg),
    prob = out_prob[nonzero],
    stringsAsFactors = FALSE
  )
}

state_id <- function(n_pos, x_pos) {
  # Unique ID for 0 <= x_pos <= n_pos <= 30.
  as.integer(n_pos * (n_pos + 1L) / 2L + x_pos + 1L)
}

max_interim_pos <- max(interim_n)
positive_states <- do.call(
  rbind,
  lapply(0:max_interim_pos, function(n) {
    data.frame(
      n_pos = rep.int(n, n + 1L),
      x_pos = 0:n,
      stringsAsFactors = FALSE
    )
  })
)
positive_states$state_id <- state_id(positive_states$n_pos,
                                     positive_states$x_pos)
n_positive_states <- nrow(positive_states)
stopifnot(max(positive_states$state_id) == n_positive_states)

aggregate_crossing_states <- function(dist) {
  out <- numeric(n_positive_states)
  if (nrow(dist) == 0L) return(out)
  ids <- state_id(dist$n_pos, dist$x_pos)
  agg <- rowsum(dist$prob, group = ids, reorder = FALSE)
  out[as.integer(rownames(agg))] <- agg[, 1L]
  out
}

build_allcomer_components <- function(scenario, A_candidates) {
  pi <- scenario$pi
  theta_pos <- scenario$theta_pos
  theta_neg <- scenario$theta_neg
  theta_all <- pi * theta_pos + (1 - pi) * theta_neg

  n1 <- interim_n[1L]
  n2 <- interim_n[2L]
  d2 <- n2 - n1
  final_increment <- N_A - n2

  dist_n1 <- make_joint_distribution(n1, pi, theta_pos, theta_neg)
  dist_increment <- make_joint_distribution(d2, pi, theta_pos, theta_neg)

  b1_values <- sort(unique(A_candidates[[paste0("b_", n1)]]))
  cache <- vector("list", length(b1_values))
  names(cache) <- as.character(b1_values)

  for (b1 in b1_values) {
    cross1 <- dist_n1[dist_n1$x_all <= b1, , drop = FALSE]
    survive1 <- dist_n1[dist_n1$x_all > b1, , drop = FALSE]
    dist_n2_survivors <- if (nrow(survive1) == 0L) {
      data.frame(
        n_pos = integer(), x_pos = integer(), x_neg = integer(),
        x_all = integer(), prob = numeric()
      )
    } else {
      convolve_distributions(survive1, dist_increment, total_n = n2)
    }
    cache[[as.character(b1)]] <- list(
      cross1_vector = aggregate_crossing_states(cross1),
      dist_n2_survivors = dist_n2_survivors
    )
  }

  nA <- nrow(A_candidates)
  W1 <- matrix(0, nrow = nA, ncol = n_positive_states)
  W2 <- matrix(0, nrow = nA, ncol = n_positive_states)
  prn_all <- numeric(nA)
  no_cross_prob <- numeric(nA)
  no_cross_expected_n_pos <- numeric(nA)

  for (a_idx in seq_len(nA)) {
    b1 <- A_candidates[[paste0("b_", n1)]][a_idx]
    b2 <- A_candidates[[paste0("b_", n2)]][a_idx]
    bF <- A_candidates[[paste0("b_", N_A)]][a_idx]

    cached <- cache[[as.character(b1)]]
    W1[a_idx, ] <- cached$cross1_vector
    dist2 <- cached$dist_n2_survivors

    if (nrow(dist2) > 0L) {
      cross2 <- dist2[dist2$x_all <= b2, , drop = FALSE]
      final_path <- dist2[dist2$x_all > b2, , drop = FALSE]
      W2[a_idx, ] <- aggregate_crossing_states(cross2)

      if (nrow(final_path) > 0L) {
        q_reject <- prob_binom_greater(
          cutoff = bF - final_path$x_all,
          size = final_increment,
          prob = theta_all
        )
        prn_all[a_idx] <- sum(final_path$prob * q_reject)
        no_cross_prob[a_idx] <- sum(final_path$prob)
        no_cross_expected_n_pos[a_idx] <- sum(
          final_path$prob * (final_path$n_pos + final_increment * pi)
        )
      }
    }
  }

  total_probability <- rowSums(W1) + rowSums(W2) + no_cross_prob
  if (max(abs(total_probability - 1)) > 1e-9) {
    stop("Exact path probabilities do not sum to one for scenario ",
         scenario$scenario_id)
  }

  list(
    scenario = scenario,
    W1 = W1,
    W2 = W2,
    prn_all = prn_all,
    no_cross_prob = no_cross_prob,
    no_cross_expected_n_pos = no_cross_expected_n_pos
  )
}

make_positive_boundary_vector <- function(P_candidates, P_index) {
  bP <- rep(NA_integer_, N_pos + 1L)
  bP[P_sample_sizes + 1L] <- as.integer(unlist(
    P_candidates[P_index, paste0("b_", P_sample_sizes), drop = FALSE],
    use.names = FALSE
  ))
  bP
}

compute_positive_branch_candidate <- function(bP, theta_pos) {
  # For every state at which the all-comer futility boundary is crossed,
  # calculate the complete biomarker-positive-only branch. When fewer than
  # n_pos_min biomarker-positive patients are available, the calculation first
  # bridges to n_pos_min and performs the initial subgroup assessment there.
  n_states <- nrow(positive_states)
  metric_names <- c(
    "success", "final_no_claim", "initial_futility", "enter",
    "stop20", "stop30", "expected_additional",
    "analyzed20", "analyzed30", "bridged"
  )
  out_mat <- matrix(
    0, nrow = n_states, ncol = length(metric_names),
    dimnames = list(NULL, metric_names)
  )

  active_memo <- new.env(hash = TRUE, parent = emptyenv())

  recurse_active <- function(n_current, x_current) {
    key <- paste0(n_current, ":", x_current)
    if (exists(key, envir = active_memo, inherits = FALSE)) {
      return(get(key, envir = active_memo, inherits = FALSE))
    }

    future_looks <- c(
      post_enrichment_n[post_enrichment_n > n_current],
      N_pos
    )
    next_look <- min(future_looks)
    increment <- next_look - n_current
    u <- 0:increment
    prob_u <- dbinom(u, size = increment, prob = theta_pos)
    x_next <- x_current + u

    ans <- c(
      success = 0,
      final_no_claim = 0,
      stop20 = 0,
      stop30 = 0,
      expected_additional = increment,
      analyzed20 = 0,
      analyzed30 = 0
    )

    if (next_look == N_pos) {
      claim <- x_next > bP[N_pos + 1L]
      ans["success"] <- sum(prob_u[claim])
      ans["final_no_claim"] <- sum(prob_u[!claim])
    } else {
      pass <- x_next > bP[next_look + 1L]
      if (next_look == 20L) {
        ans["analyzed20"] <- 1
        ans["stop20"] <- sum(prob_u[!pass])
      }
      if (next_look == 30L) {
        ans["analyzed30"] <- 1
        ans["stop30"] <- sum(prob_u[!pass])
      }

      pass_index <- which(pass)
      if (length(pass_index) > 0L) {
        for (k in pass_index) {
          child <- recurse_active(next_look, x_next[k])
          weight <- prob_u[k]
          ans["success"] <- ans["success"] + weight * child["success"]
          ans["final_no_claim"] <- ans["final_no_claim"] +
            weight * child["final_no_claim"]
          ans["stop20"] <- ans["stop20"] + weight * child["stop20"]
          ans["stop30"] <- ans["stop30"] + weight * child["stop30"]
          ans["expected_additional"] <- ans["expected_additional"] +
            weight * child["expected_additional"]
          ans["analyzed20"] <- ans["analyzed20"] +
            weight * child["analyzed20"]
          ans["analyzed30"] <- ans["analyzed30"] +
            weight * child["analyzed30"]
        }
      }
    }

    partition <- ans["success"] + ans["final_no_claim"] +
      ans["stop20"] + ans["stop30"]
    if (abs(partition - 1) > 1e-10) {
      stop("Active biomarker-positive path probabilities do not sum to one.")
    }

    assign(key, ans, envir = active_memo)
    ans
  }

  for (r in seq_len(n_states)) {
    m <- positive_states$n_pos[r]
    x <- positive_states$x_pos[r]

    if (m < n_pos_min) {
      d0 <- n_pos_min - m
      u0 <- 0:d0
      prob0 <- dbinom(u0, size = d0, prob = theta_pos)
      x0 <- x + u0
      pass0 <- x0 > bP[n_pos_min + 1L]

      out_mat[r, "bridged"] <- 1
      out_mat[r, "initial_futility"] <- sum(prob0[!pass0])
      out_mat[r, "enter"] <- sum(prob0[pass0])
      out_mat[r, "expected_additional"] <- d0

      pass_index <- which(pass0)
      if (length(pass_index) > 0L) {
        for (k in pass_index) {
          child <- recurse_active(n_pos_min, x0[k])
          weight <- prob0[k]
          for (nm in c(
            "success", "final_no_claim", "stop20", "stop30",
            "analyzed20", "analyzed30"
          )) {
            out_mat[r, nm] <- out_mat[r, nm] + weight * child[nm]
          }
          out_mat[r, "expected_additional"] <-
            out_mat[r, "expected_additional"] +
            weight * child["expected_additional"]
        }
      }
    } else {
      pass0 <- x > bP[m + 1L]
      out_mat[r, "initial_futility"] <- as.numeric(!pass0)
      out_mat[r, "enter"] <- as.numeric(pass0)

      # When the first assessment itself occurs at a scheduled sample size, it
      # counts as that analysis and is not repeated later.
      if (m == 20L) out_mat[r, "analyzed20"] <- 1
      if (m == 30L) out_mat[r, "analyzed30"] <- 1

      if (pass0) {
        child <- recurse_active(m, x)
        for (nm in c(
          "success", "final_no_claim", "stop20", "stop30",
          "analyzed20", "analyzed30"
        )) {
          out_mat[r, nm] <- out_mat[r, nm] + child[nm]
        }
        out_mat[r, "expected_additional"] <- child["expected_additional"]
      }
    }
  }

  initial_partition <- out_mat[, "initial_futility"] + out_mat[, "enter"]
  if (max(abs(initial_partition - 1)) > 1e-10) {
    stop("Initial biomarker-positive assessment probabilities do not sum to one.")
  }
  terminal_partition <- out_mat[, "success"] +
    out_mat[, "final_no_claim"] + out_mat[, "initial_futility"] +
    out_mat[, "stop20"] + out_mat[, "stop30"]
  if (max(abs(terminal_partition - 1)) > 1e-10) {
    stop("Biomarker-positive branch probabilities do not sum to one.")
  }
  if (any(out_mat[, "stop20"] - out_mat[, "analyzed20"] > 1e-10) ||
      any(out_mat[, "stop30"] - out_mat[, "analyzed30"] > 1e-10)) {
    stop("A stopping probability exceeds its corresponding analysis probability.")
  }

  as.data.frame(out_mat, stringsAsFactors = FALSE)
}

make_positive_kernels <- function(P_candidates, theta_pos_values) {
  nP <- nrow(P_candidates)
  success <- vector("list", length(theta_pos_values))
  names(success) <- format(theta_pos_values, trim = TRUE)

  for (theta_pos in theta_pos_values) {
    message("  constructing positive-branch success kernels for theta_pos = ",
            theta_pos)
    candidate_results <- parallel_lapply(
      seq_len(nP),
      function(p_idx) {
        bP <- make_positive_boundary_vector(P_candidates, p_idx)
        compute_positive_branch_candidate(
          bP = bP,
          theta_pos = theta_pos
        )$success
      },
      n_workers = n_cores
    )
    success[[format(theta_pos, trim = TRUE)]] <- do.call(
      cbind,
      candidate_results
    )
  }

  list(
    success = success,
    sample_sizes = P_sample_sizes
  )
}

positive_detail_cache <- new.env(hash = TRUE, parent = emptyenv())

get_positive_detail <- function(P_candidates, P_index, theta_pos) {
  key <- paste(P_index, format(theta_pos, digits = 17), sep = ":")
  if (!exists(key, envir = positive_detail_cache, inherits = FALSE)) {
    bP <- make_positive_boundary_vector(P_candidates, P_index)
    detail <- compute_positive_branch_candidate(
      bP = bP,
      theta_pos = theta_pos
    )
    assign(key, detail, envir = positive_detail_cache)
  }
  get(key, envir = positive_detail_cache, inherits = FALSE)
}

calibration_prn_matrices <- function(component, positive_kernels) {
  W <- component$W1 + component$W2
  scenario <- component$scenario
  theta_key <- format(scenario$theta_pos, trim = TRUE)
  success_kernel <- positive_kernels$success[[theta_key]]
  if (is.null(success_kernel)) {
    stop("No positive success kernel for theta_pos = ", scenario$theta_pos)
  }

  prn_positive <- W %*% success_kernel
  prn_any <- sweep(prn_positive, MARGIN = 1L,
                   STATS = component$prn_all, FUN = "+")

  list(
    PRN_any = prn_any,
    PRN_all = component$prn_all,
    PRN_positive = prn_positive
  )
}

exact_metrics_selected <- function(component, positive_kernels,
                                   A_index, P_index) {
  w1 <- component$W1[A_index, ]
  w2 <- component$W2[A_index, ]
  w <- w1 + w2
  scenario <- component$scenario

  theta_key <- format(scenario$theta_pos, trim = TRUE)
  success_vector <- positive_kernels$success[[theta_key]][, P_index]
  detail <- get_positive_detail(
    P_candidates = P_candidates,
    P_index = P_index,
    theta_pos = scenario$theta_pos
  )
  if (max(abs(success_vector - detail$success)) > 1e-10) {
    stop("Cached and detailed biomarker-positive success probabilities differ.")
  }

  n_state <- positive_states$n_pos
  prn_all <- component$prn_all[A_index]
  prn_positive <- sum(w * detail$success)

  first_assessment <- sum(w)
  bridged <- sum(w * detail$bridged)
  pet_positive_initial <- sum(w * detail$initial_futility)
  enter <- sum(w * detail$enter)
  pet_positive_look20 <- sum(w * detail$stop20)
  pet_positive_look30 <- sum(w * detail$stop30)
  pet_positive <- pet_positive_initial + pet_positive_look20 +
    pet_positive_look30
  analyzed20 <- sum(w * detail$analyzed20)
  analyzed30 <- sum(w * detail$analyzed30)

  if (abs(first_assessment - pet_positive_initial - enter) > 1e-9) {
    stop("First biomarker-positive assessment probabilities do not sum correctly.")
  }
  terminal_partition <- sum(
    w * (detail$success + detail$final_no_claim +
           detail$initial_futility + detail$stop20 + detail$stop30)
  )
  if (abs(terminal_partition - first_assessment) > 1e-9) {
    stop("Biomarker-positive branch terminal probabilities do not sum correctly.")
  }
  if (bridged - first_assessment > 1e-9 ||
      pet_positive_look20 - analyzed20 > 1e-9 ||
      pet_positive_look30 - analyzed30 > 1e-9) {
    stop("Detailed biomarker-positive path metrics are internally inconsistent.")
  }

  mean_total_sample_size <-
    component$no_cross_prob[A_index] * N_A +
    sum(w1) * interim_n[1L] +
    sum(w2) * interim_n[2L] +
    sum(w * detail$expected_additional)

  mean_positive_sample_size <-
    component$no_cross_expected_n_pos[A_index] +
    sum(w * n_state) +
    sum(w * detail$expected_additional)

  c(
    PRN_any = prn_all + prn_positive,
    PRN_all = prn_all,
    PRN_positive = prn_positive,
    PET_any = pet_positive,
    PET_positive = pet_positive,
    Allcomer_futility = first_assessment,
    First_positive_assessment = first_assessment,
    Bridge_to_minimum = bridged,
    PET_positive_initial = pet_positive_initial,
    Enter_enrichment = enter,
    Enter_enrichment_look1 = sum(w1 * detail$enter),
    Enter_enrichment_look2 = sum(w2 * detail$enter),
    Analyzed_positive_look20 = analyzed20,
    PET_positive_look20 = pet_positive_look20,
    Analyzed_positive_look30 = analyzed30,
    PET_positive_look30 = pet_positive_look30,
    Mean_total_sample_size = mean_total_sample_size,
    Mean_positive_sample_size = mean_positive_sample_size
  )
}

select_candidate_pair <- function(feasible, avg_power, min_power,
                                  type1_metric, label) {
  feasible_idx <- which(feasible)
  if (length(feasible_idx) == 0L) {
    stop("No feasible candidate pair for ", label)
  }

  tol <- 1e-12
  best <- feasible_idx
  target <- max(avg_power[best])
  best <- best[abs(avg_power[best] - target) <= tol]

  target <- max(min_power[best])
  best <- best[abs(min_power[best] - target) <= tol]

  target <- max(type1_metric[best])
  best <- best[abs(type1_metric[best] - target) <= tol]

  chosen <- best[1L]
  ij <- arrayInd(chosen, dim(feasible))
  list(
    linear_index = chosen,
    A_index = ij[1L],
    P_index = ij[2L],
    average_power = avg_power[chosen],
    minimum_power = min_power[chosen],
    type1_metric = type1_metric[chosen]
  )
}

## -----------------------------------------------------------------------------
## 3. Candidate boundary tables
## -----------------------------------------------------------------------------

message("Generating unique all-comer boundary candidates...")
A_sample_sizes <- c(interim_n, N_A)
A_candidates <- make_boundary_candidates(
  sample_sizes = A_sample_sizes,
  max_sample_size = N_A,
  lambda_values = lambda_grid,
  gamma_values = gamma_grid,
  prefix = "all-comer"
)

message("Generating unique biomarker-positive boundary candidates...")
P_sample_sizes <- sort(unique(c(
  n_pos_min:max_interim_pos, post_enrichment_n, N_pos
)))
P_candidates <- make_boundary_candidates(
  sample_sizes = P_sample_sizes,
  max_sample_size = N_pos,
  lambda_values = lambda_grid,
  gamma_values = gamma_grid,
  prefix = "biomarker-positive"
)

message("Unique boundary candidates: all-comer = ", nrow(A_candidates),
        "; biomarker-positive = ", nrow(P_candidates),
        "; candidate pairs = ",
        format(nrow(A_candidates) * nrow(P_candidates), big.mark = ","))

write.csv(A_candidates,
          file.path(output_dir, "allcomer_boundary_candidates.csv"),
          row.names = FALSE)
write.csv(P_candidates,
          file.path(output_dir, "positive_boundary_candidates.csv"),
          row.names = FALSE)
candidate_counts <- data.frame(
  component = c("all-comer", "biomarker-positive", "candidate pairs"),
  count = c(
    nrow(A_candidates),
    nrow(P_candidates),
    nrow(A_candidates) * nrow(P_candidates)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  candidate_counts,
  file.path(output_dir, "candidate_counts.csv"),
  row.names = FALSE
)

## -----------------------------------------------------------------------------
## 4. Scenarios and exact path components
## -----------------------------------------------------------------------------

null_scenarios <- data.frame(
  scenario_type = "Null",
  pi = prevalence_grid,
  theta_pos = p0,
  theta_neg = p0,
  stringsAsFactors = FALSE
)

alternative_scenarios <- expand.grid(
  pi = prevalence_grid,
  theta_pos = theta_pos_grid,
  theta_neg = theta_neg_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
alternative_scenarios <- alternative_scenarios[
  order(alternative_scenarios$pi,
        alternative_scenarios$theta_pos,
        -alternative_scenarios$theta_neg),
]
alternative_scenarios$scenario_type <- "Alternative"

# As in the original BOP2 optimization, one prespecified alternative
# configuration is used to select the design. The remaining response
# configurations are used only to evaluate the operating characteristics of
# that fixed design and do not trigger recalibration.
alternative_scenarios$used_for_calibration <-
  abs(alternative_scenarios$theta_pos - p1) < 1e-12 &
  abs(alternative_scenarios$theta_neg - p0) < 1e-12
null_scenarios$used_for_calibration <- FALSE

alternative_scenarios <- alternative_scenarios[
  , c(
    "scenario_type", "pi", "theta_pos", "theta_neg",
    "used_for_calibration"
  )
]
null_scenarios <- null_scenarios[
  , c(
    "scenario_type", "pi", "theta_pos", "theta_neg",
    "used_for_calibration"
  )
]

all_scenarios <- rbind(null_scenarios, alternative_scenarios)
all_scenarios$scenario_id <- seq_len(nrow(all_scenarios))
all_scenarios$scenario_label <- sprintf(
  "%s_pi%.1f_pos%.1f_neg%.1f",
  all_scenarios$scenario_type,
  all_scenarios$pi,
  all_scenarios$theta_pos,
  all_scenarios$theta_neg
)

write.csv(all_scenarios,
          file.path(output_dir, "operating_characteristic_scenarios.csv"),
          row.names = FALSE)

positive_theta_values <- sort(unique(all_scenarios$theta_pos))
message("Constructing biomarker-positive decision kernels...")
positive_kernels <- make_positive_kernels(P_candidates,
                                          positive_theta_values)

message("Constructing exact all-comer path components for ",
        nrow(all_scenarios), " scenarios...")
scenario_rows <- split(all_scenarios, seq_len(nrow(all_scenarios)))
components <- parallel_lapply(
  scenario_rows,
  function(s) build_allcomer_components(s[1L, ], A_candidates),
  n_workers = n_cores
)

## -----------------------------------------------------------------------------
## 5. Calibration of the proposed and comparator designs
## -----------------------------------------------------------------------------

nA <- nrow(A_candidates)
nP <- nrow(P_candidates)

max_global_type1 <- matrix(-Inf, nrow = nA, ncol = nP)
max_all_type1 <- rep(-Inf, nA)
max_positive_type1 <- matrix(-Inf, nrow = nA, ncol = nP)

avg_power <- matrix(0, nrow = nA, ncol = nP)
min_power <- matrix(Inf, nrow = nA, ncol = nP)

n_calibration_alt <- sum(alternative_scenarios$used_for_calibration)
if (n_calibration_alt != length(prevalence_grid)) {
  stop(
    "The single calibration alternative must contain one scenario for each ",
    "prevalence value."
  )
}

message("Evaluating exact null and alternative operating characteristics...")
for (s_idx in seq_along(components)) {
  metrics <- calibration_prn_matrices(components[[s_idx]], positive_kernels)
  scenario_type <- all_scenarios$scenario_type[s_idx]

  if (scenario_type == "Null") {
    max_global_type1 <- pmax(max_global_type1, metrics$PRN_any)
    max_all_type1 <- pmax(max_all_type1, metrics$PRN_all)
    max_positive_type1 <- pmax(max_positive_type1, metrics$PRN_positive)
  } else if (isTRUE(all_scenarios$used_for_calibration[s_idx])) {
    avg_power <- avg_power + metrics$PRN_any / n_calibration_alt
    min_power <- pmin(min_power, metrics$PRN_any)
  }

  if (s_idx %% 5L == 0L || s_idx == length(components)) {
    message("  completed ", s_idx, " / ", length(components), " scenarios")
  }
}

feasible_global <- max_global_type1 <= alpha + 1e-12

feasible_componentwise <-
  matrix(max_all_type1 <= alpha + 1e-12, nrow = nA, ncol = nP) &
  max_positive_type1 <= alpha + 1e-12

selected_global <- select_candidate_pair(
  feasible = feasible_global,
  avg_power = avg_power,
  min_power = min_power,
  type1_metric = max_global_type1,
  label = "globally calibrated design"
)
selected_componentwise <- select_candidate_pair(
  feasible = feasible_componentwise,
  avg_power = avg_power,
  min_power = min_power,
  type1_metric = pmax(matrix(max_all_type1, nrow = nA, ncol = nP),
                      max_positive_type1),
  label = "componentwise-calibrated comparator"
)

make_selected_row <- function(method, selected) {
  a <- A_candidates[selected$A_index, ]
  p <- P_candidates[selected$P_index, ]
  data.frame(
    Method = method,
    A_candidate_id = a$candidate_id,
    P_candidate_id = p$candidate_id,
    lambda_A = a$lambda,
    gamma_A = a$gamma,
    lambda_positive = p$lambda,
    gamma_positive = p$gamma,
    max_global_type1 = max_global_type1[selected$A_index,
                                         selected$P_index],
    max_allcomer_type1 = max_all_type1[selected$A_index],
    max_positive_type1 = max_positive_type1[selected$A_index,
                                             selected$P_index],
    average_power = selected$average_power,
    minimum_power = selected$minimum_power,
    stringsAsFactors = FALSE
  )
}

selected_designs <- rbind(
  make_selected_row("Globally calibrated", selected_global),
  make_selected_row("Componentwise calibrated", selected_componentwise)
)

if (any(abs(selected_designs$gamma_A - max(gamma_grid)) < 1e-12) ||
    any(abs(selected_designs$gamma_positive - max(gamma_grid)) < 1e-12)) {
  warning(
    "At least one selected gamma value is at the upper search limit. ",
    "Consider extending gamma_grid to verify that the selected boundary table ",
    "is not constrained by the numerical bound."
  )
}
write.csv(selected_designs,
          file.path(output_dir, "selected_designs.csv"),
          row.names = FALSE)
print(selected_designs)

saveRDS(
  list(
    settings = list(
      p0 = p0, p1 = p1, alpha = alpha,
      interim_n = interim_n, N_A = N_A, N_pos = N_pos,
      n_pos_min = n_pos_min, post_enrichment_n = post_enrichment_n,
      prior_a = prior_a, prior_b = prior_b,
      prevalence_grid = prevalence_grid,
      theta_pos_grid = theta_pos_grid,
      theta_neg_grid = theta_neg_grid,
      lambda_grid = lambda_grid,
      gamma_grid = gamma_grid
    ),
    A_candidates = A_candidates,
    P_candidates = P_candidates,
    selected_global = selected_global,
    selected_componentwise = selected_componentwise,
    selected_designs = selected_designs
  ),
  file.path(output_dir, "calibration_objects.rds")
)

## -----------------------------------------------------------------------------
## 6. Exact operating characteristics for selected designs
## -----------------------------------------------------------------------------

selected_index <- list(
  "Globally calibrated" = c(selected_global$A_index,
                              selected_global$P_index),
  "Componentwise calibrated" = c(selected_componentwise$A_index,
                                   selected_componentwise$P_index)
)

exact_rows <- vector("list", length(components) * length(selected_index))
out_idx <- 1L
for (s_idx in seq_along(components)) {
  for (method in names(selected_index)) {
    ij <- selected_index[[method]]
    vals <- exact_metrics_selected(
      component = components[[s_idx]],
      positive_kernels = positive_kernels,
      A_index = ij[1L],
      P_index = ij[2L]
    )
    exact_rows[[out_idx]] <- cbind(
      all_scenarios[s_idx, ],
      Method = method,
      as.data.frame(as.list(vals), stringsAsFactors = FALSE)
    )
    out_idx <- out_idx + 1L
  }
}
exact_results <- do.call(rbind, exact_rows)
rownames(exact_results) <- NULL
write.csv(exact_results,
          file.path(output_dir, "exact_operating_characteristics.csv"),
          row.names = FALSE)

## -----------------------------------------------------------------------------
## 7. Table-I-style boundary table for the proposed design only
## -----------------------------------------------------------------------------

selected_A <- A_candidates[selected_global$A_index, ]
selected_P <- P_candidates[selected_global$P_index, ]

allcomer_boundary_table <- data.frame(
  Population = "All-comer",
  Rule = "Futile/no efficacy claim if number of responses <=",
  t(as.integer(unlist(selected_A[1L, paste0("b_", A_sample_sizes)], use.names = FALSE))),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(allcomer_boundary_table)[-(1:2)] <- as.character(A_sample_sizes)

positive_boundary_table <- data.frame(
  Population = "Biomarker-positive",
  Rule = "Futile/no efficacy claim if number of responses <=",
  t(as.integer(unlist(selected_P[1L, paste0("b_", P_sample_sizes)], use.names = FALSE))),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(positive_boundary_table)[-(1:2)] <- as.character(P_sample_sizes)

write.csv(allcomer_boundary_table,
          file.path(output_dir, "proposed_boundary_table_allcomer.csv"),
          row.names = FALSE)
write.csv(positive_boundary_table,
          file.path(output_dir, "proposed_boundary_table_positive.csv"),
          row.names = FALSE)

latex_vector <- function(x) paste(x, collapse = " & ")
latex_boundary_vector <- function(x) {
  paste(ifelse(x < 0L, "--", as.character(x)), collapse = " & ")
}

pos_block1 <- n_pos_min:20L
pos_block2 <- c(21L:30L, N_pos)
get_p_boundary <- function(samples) {
  as.integer(unlist(selected_P[1L, paste0("b_", samples)], use.names = FALSE))
}

latex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Decision boundaries of the proposed globally calibrated enrichment BOP2 design.}",
  "\\label{tab:proposed_boundaries}",
  "\\begin{tabular}{lccc}",
  "\\hline",
  paste0("All-comer sample size & ", latex_vector(A_sample_sizes), " \\\\"),
  paste0("Futile if number of responses $\\le$ & ",
         latex_boundary_vector(as.integer(unlist(selected_A[1L,
           paste0("b_", A_sample_sizes)], use.names = FALSE))), " \\\\"),
  "\\hline",
  "\\end{tabular}",
  "\\vspace{0.6em}",
  "",
  paste0("\\begin{tabular}{l", paste(rep("c", length(pos_block1)),
                                       collapse = ""), "}"),
  "\\hline",
  paste0("Biomarker-positive sample size & ",
         latex_vector(pos_block1), " \\\\"),
  paste0("Futile if number of responses $\\le$ & ",
         latex_boundary_vector(get_p_boundary(pos_block1)), " \\\\"),
  "\\hline",
  "\\end{tabular}",
  "\\vspace{0.6em}",
  "",
  paste0("\\begin{tabular}{l", paste(rep("c", length(pos_block2)),
                                       collapse = ""), "}"),
  "\\hline",
  paste0("Biomarker-positive sample size & ",
         latex_vector(pos_block2), " \\\\"),
  paste0("Futile if number of responses $\\le$ & ",
         latex_boundary_vector(get_p_boundary(pos_block2)), " \\\\"),
  "\\hline",
  "\\end{tabular}",
  "",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\footnotesize The all-comer boundaries apply at the interim analyses of 20 and 30 patients and at the final analysis of 40 patients. Biomarker-positive boundaries apply at the first biomarker-positive assessment and, when applicable, at subsequent analyses of 20 and 30 cumulative biomarker-positive patients and at the final analysis of 40 biomarker-positive patients. At an interim analysis, satisfying the listed rule results in stopping for futility. At the final analysis, satisfying the listed rule results in no efficacy claim.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(latex_lines,
           file.path(output_dir, "proposed_boundary_table.tex"))

## -----------------------------------------------------------------------------
## 8. Exact result summaries and session information
## -----------------------------------------------------------------------------

type1_exact <- exact_results[
  exact_results$scenario_type == "Null",
  ,
  drop = FALSE
]
power_exact <- exact_results[
  exact_results$scenario_type == "Alternative",
  ,
  drop = FALSE
]

write.csv(
  type1_exact,
  file.path(output_dir, "exact_type1_error.csv"),
  row.names = FALSE
)
write.csv(
  power_exact,
  file.path(output_dir, "exact_power.csv"),
  row.names = FALSE
)

settings_table <- data.frame(
  Parameter = c(
    "Analysis type", "p0", "p1", "alpha",
    "calibration alternative",
    "all-comer interim looks", "post-enrichment positive looks",
    "N_A", "N_positive", "n_positive_min", "prior_a", "prior_b",
    "prevalence grid", "theta_positive grid", "theta_negative grid",
    "lambda grid", "gamma grid", "detected physical cores",
    "number of workers used", "parallel backend"
  ),
  Value = c(
    "Exact recursive enumeration", p0, p1, alpha,
    paste0("theta_positive = ", p1, "; theta_negative = ", p0),
    paste(interim_n, collapse = ", "),
    paste(post_enrichment_n, collapse = ", "),
    N_A, N_pos, n_pos_min, prior_a, prior_b,
    paste(prevalence_grid, collapse = ", "),
    paste(theta_pos_grid, collapse = ", "),
    paste(theta_neg_grid, collapse = ", "),
    paste0(min(lambda_grid), " to ", max(lambda_grid),
           " by ", unique(round(diff(lambda_grid), 12))),
    paste0(min(gamma_grid), " to ", max(gamma_grid),
           " by ", unique(round(diff(gamma_grid), 12))),
    detected_cores, n_cores, parallel_backend
  ),
  stringsAsFactors = FALSE
)
write.csv(
  settings_table,
  file.path(output_dir, "analysis_settings.csv"),
  row.names = FALSE
)

xlsx_path <- file.path(output_dir, "enrichment_bop2_exact_results.xlsx")
xlsx_created <- write_results_workbook(
  sheets = list(
    Settings = settings_table,
    Selected_designs = selected_designs,
    Scenarios = all_scenarios,
    Type1_exact = type1_exact,
    Power_exact = power_exact,
    Exact_OC = exact_results,
    Candidate_counts = candidate_counts,
    Boundary_allcomer = allcomer_boundary_table,
    Boundary_positive = positive_boundary_table
  ),
  path = xlsx_path
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("Exact analysis completed. Results were written to: ",
        normalizePath(output_dir, mustWork = FALSE))
if (xlsx_created) {
  message("Consolidated exact-results workbook: ",
          normalizePath(xlsx_path, mustWork = FALSE))
} else {
  message("All exact results are available as CSV files in the result directory.")
}
