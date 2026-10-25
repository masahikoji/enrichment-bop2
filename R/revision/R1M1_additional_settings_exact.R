# Additional benchmark analysis for R1-M1 and R1-M3 (binary endpoint).
# Requested settings: (p0,p1)=(0.10,0.30) and (0.30,0.50).
# Original (0.20,0.40) is included by default as a regression control.
#
# Numerical engine: supplied FINAL CORRECTED no-bridging R program.
# This is EXACT finite-state enumeration, not a Monte Carlo simulation.
# Boundaries are re-calibrated for each (p0,p1), then held fixed for evaluation.
# Calibration remains under the POINT complete null and the given prevalence
# grid. No strong-FWER or error-matched-comparator extension is introduced here.
#
# Usage:
#   Rscript R1M1_additional_settings_exact.R
# or in R/RStudio:
#   source("R1M1_additional_settings_exact.R")
#
# Optional environment variables:
#   BOP2_N_CORES=8              parallel workers (default: capped at 8)
#   BOP2_INCLUDE_BASELINE=0     run only the two added settings
#   BOP2_R1M1_OUTPUT_DIR=...    NEW/nonexistent output directory
#   BOP2_QUICK=1               coarse-grid smoke test, NOT manuscript results
#   BOP2_WRITE_XLSX=1          optional workbooks if a writer is installed
#   BOP2_PROPOSED_GAMMA_MAX=12 optional numerical search sensitivity (default 8)
# Base R is sufficient for all scientific calculations and CSV outputs.
# No files outside the newly created output directory are modified.
# Original uploaded program SHA256: 3b6fab6e727f13c20ff2085023fcc64901c884bec7eec1f0b76a95faab50de0a

run_bop2_setting <- function(p0_input, p1_input, output_root) {
  stopifnot(is.numeric(p0_input), length(p0_input)==1L, is.finite(p0_input),
    is.numeric(p1_input), length(p1_input)==1L, is.finite(p1_input),
    p0_input >= 0.10, p1_input > p0_input, p1_input + 0.20 <= 1)
# Analysis objects and caches are local; existing user objects are not cleared.

## -----------------------------------------------------------------------------
## 1. User settings
## -----------------------------------------------------------------------------

p0 <- p0_input
p1 <- p1_input
alpha <- 0.10
alpha_lower <- 0.09

interim_n <- c(20L, 30L)
N_A <- 40L
N_pos <- 40L
n_pos_min <- 10L
post_enrichment_n <- c(20L, 30L)

prior_a <- 0.5
prior_b <- 0.5
prevalence_grid <- c(0.4, 0.5, 0.6, 0.7, 0.8)

theta_pos_grid <- round(p1 + c(0, 0.10, 0.20), 10)
theta_neg_grid <- round(p0 + c(0, -0.10), 10)

# The proposed design uses the grid prespecified in the manuscript.
proposed_lambda_grid <- seq(0.005, 1.000, by = 0.005)
gamma_max <- suppressWarnings(as.numeric(Sys.getenv(
  "BOP2_PROPOSED_GAMMA_MAX", unset = "8"
)))
if (length(gamma_max) != 1L || !is.finite(gamma_max) || gamma_max < 0.05 ||
    abs(gamma_max / 0.05 - round(gamma_max / 0.05)) > 1e-8) {
  stop("BOP2_PROPOSED_GAMMA_MAX must be a positive multiple of 0.05.")
}
proposed_gamma_grid <- seq(0.050, gamma_max, by = 0.050)

# The conventional BOP2 comparator uses the finer grid employed in the
# accompanying BOP2 reproduction program. Its selection rule is standalone
# BOP2 calibration and is completely separate from selection of the proposed
# design.
bop2_lambda_grid <- seq(0.001, 1.000, by = 0.0025)
bop2_gamma_grid <- seq(0.010, 10.000, by = 0.025)

# A separate output directory is used for each (p0,p1) setting.
project_dir <- output_root
case_id <- sprintf("p0_%.2f_p1_%.2f", p0, p1)
result_root <- file.path(output_root, case_id)
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
  alpha_lower <- 0
  proposed_lambda_grid <- seq(0.05, 1.00, by = 0.05)
  proposed_gamma_grid <- seq(0.25, 2.00, by = 0.25)
  bop2_lambda_grid <- proposed_lambda_grid
  bop2_gamma_grid <- proposed_gamma_grid
  prevalence_grid <- c(0.4, 0.8)
  theta_pos_grid <- round(p1 + c(0, 0.20), 10)
  theta_neg_grid <- round(p0 + c(0, -0.10), 10)
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

# No existing result files are deleted by this script.
stopifnot(
  alpha_lower >= 0,
  alpha_lower <= alpha,
  length(interim_n) == 2L,
  all(interim_n > 0L),
  all(diff(interim_n) > 0L),
  max(interim_n) < N_A,
  N_A == N_pos,
  n_pos_min <= min(interim_n),
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
  if (!identical(Sys.getenv("BOP2_WRITE_XLSX", unset = "0"), "1")) return(FALSE)
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

compute_positive_path_candidate <- function(bP, theta_pos, enter_vector) {
  # For each state at which the all-comer boundary is crossed, calculate the
  # probability of eventual biomarker-positive success and the probabilities
  # of stopping at the post-enrichment looks. The initial enrichment assessment
  # is represented by enter_vector and is not counted again as a scheduled
  # post-enrichment look when n_pos is exactly 20 or 30.
  n_states <- nrow(positive_states)
  success <- numeric(n_states)
  final_no_claim <- numeric(n_states)
  stop20 <- numeric(n_states)
  stop30 <- numeric(n_states)
  expected_additional <- numeric(n_states)
  reach20 <- numeric(n_states)
  reach30 <- numeric(n_states)

  memo <- new.env(hash = TRUE, parent = emptyenv())

  recurse <- function(n_current, x_current) {
    key <- paste0(n_current, ":", x_current)
    if (exists(key, envir = memo, inherits = FALSE)) {
      return(get(key, envir = memo, inherits = FALSE))
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

    out <- c(
      success = 0,
      final_no_claim = 0,
      stop20 = 0,
      stop30 = 0,
      expected_additional = increment,
      reach20 = 0,
      reach30 = 0
    )

    if (next_look == N_pos) {
      claim <- x_next > bP[N_pos + 1L]
      out["success"] <- sum(prob_u[claim])
      out["final_no_claim"] <- sum(prob_u[!claim])
    } else {
      pass <- x_next > bP[next_look + 1L]
      if (next_look == 20L) {
        out["reach20"] <- 1
        out["stop20"] <- sum(prob_u[!pass])
      }
      if (next_look == 30L) {
        out["reach30"] <- 1
        out["stop30"] <- sum(prob_u[!pass])
      }

      pass_index <- which(pass)
      if (length(pass_index) > 0L) {
        for (k in pass_index) {
          child <- recurse(next_look, x_next[k])
          weight <- prob_u[k]
          out["success"] <- out["success"] + weight * child["success"]
          out["final_no_claim"] <- out["final_no_claim"] +
            weight * child["final_no_claim"]
          out["stop20"] <- out["stop20"] + weight * child["stop20"]
          out["stop30"] <- out["stop30"] + weight * child["stop30"]
          out["expected_additional"] <- out["expected_additional"] +
            weight * child["expected_additional"]
          out["reach20"] <- out["reach20"] + weight * child["reach20"]
          out["reach30"] <- out["reach30"] + weight * child["reach30"]
        }
      }
    }

    assign(key, out, envir = memo)
    out
  }

  entered_states <- which(enter_vector > 0)
  for (r in entered_states) {
    ans <- recurse(positive_states$n_pos[r], positive_states$x_pos[r])
    success[r] <- ans["success"]
    final_no_claim[r] <- ans["final_no_claim"]
    stop20[r] <- ans["stop20"]
    stop30[r] <- ans["stop30"]
    expected_additional[r] <- ans["expected_additional"]
    reach20[r] <- ans["reach20"]
    reach30[r] <- ans["reach30"]
  }

  if (length(entered_states) > 0L) {
    total <- success + final_no_claim + stop20 + stop30
    if (max(abs(total[entered_states] - 1)) > 1e-10) {
      stop("Biomarker-positive path probabilities do not sum to one.")
    }
  }

  list(
    success = success,
    final_no_claim = final_no_claim,
    stop20 = stop20,
    stop30 = stop30,
    expected_additional = expected_additional,
    reach20 = reach20,
    reach30 = reach30
  )
}

make_positive_kernels <- function(P_candidates, theta_pos_values) {
  boundary_columns <- paste0("b_", P_sample_sizes)
  B <- as.matrix(P_candidates[, boundary_columns, drop = FALSE])
  nP <- nrow(P_candidates)

  boundary_by_state <- matrix(Inf, nrow = n_positive_states, ncol = nP)
  assessed <- positive_states$n_pos >= n_pos_min
  for (n in n_pos_min:max_interim_pos) {
    state_rows <- which(positive_states$n_pos == n)
    boundary_col <- match(paste0("b_", n), boundary_columns)
    boundary_by_state[state_rows, ] <- matrix(
      B[, boundary_col],
      nrow = length(state_rows),
      ncol = nP,
      byrow = TRUE
    )
  }

  x_matrix <- matrix(
    positive_states$x_pos,
    nrow = n_positive_states,
    ncol = nP
  )
  enter <- assessed & (x_matrix > boundary_by_state)
  positive_futility_initial <- assessed & !enter

  success <- vector("list", length(theta_pos_values))
  names(success) <- format(theta_pos_values, trim = TRUE)

  for (theta_pos in theta_pos_values) {
    message("  constructing positive-path success kernels for theta_pos = ",
            theta_pos)
    candidate_results <- parallel_lapply(
      seq_len(nP),
      function(p_idx) {
        bP <- make_positive_boundary_vector(P_candidates, p_idx)
        compute_positive_path_candidate(
          bP = bP,
          theta_pos = theta_pos,
          enter_vector = enter[, p_idx]
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
    enter = enter * 1,
    positive_futility_initial = positive_futility_initial * 1,
    success = success,
    sample_sizes = P_sample_sizes
  )
}

positive_detail_cache <- new.env(hash = TRUE, parent = emptyenv())

get_positive_detail <- function(P_candidates, positive_kernels,
                                P_index, theta_pos) {
  key <- paste(P_index, format(theta_pos, digits = 17), sep = ":")
  if (!exists(key, envir = positive_detail_cache, inherits = FALSE)) {
    bP <- make_positive_boundary_vector(P_candidates, P_index)
    detail <- compute_positive_path_candidate(
      bP = bP,
      theta_pos = theta_pos,
      enter_vector = positive_kernels$enter[, P_index]
    )
    assign(key, detail, envir = positive_detail_cache)
  }
  get(key, envir = positive_detail_cache, inherits = FALSE)
}

make_entry_specific_positive_kernels <- function(
    P_candidates,
    positive_library,
    theta_pos_values
) {
  n_states <- nrow(positive_states)
  enter <- numeric(n_states)
  positive_futility_initial <- numeric(n_states)

  for (entry_n in n_pos_min:max_interim_pos) {
    library_row <- positive_library[
      positive_library$entry_n == entry_n,
      ,
      drop = FALSE
    ]
    if (nrow(library_row) != 1L) {
      stop("Expected one positive BOP2 design for entry_n = ", entry_n)
    }
    P_index <- library_row$P_index
    bP <- make_positive_boundary_vector(P_candidates, P_index)
    state_rows <- which(positive_states$n_pos == entry_n)
    pass <- positive_states$x_pos[state_rows] > bP[entry_n + 1L]
    enter[state_rows] <- as.numeric(pass)
    positive_futility_initial[state_rows] <- as.numeric(!pass)
  }

  success <- vector("list", length(theta_pos_values))
  details <- vector("list", length(theta_pos_values))
  names(success) <- format(theta_pos_values, trim = TRUE)
  names(details) <- format(theta_pos_values, trim = TRUE)

  for (theta_pos in theta_pos_values) {
    message(
      "  constructing entry-specific conventional BOP2 kernels for theta_pos = ",
      theta_pos
    )
    combined <- list(
      success = numeric(n_states),
      final_no_claim = numeric(n_states),
      stop20 = numeric(n_states),
      stop30 = numeric(n_states),
      expected_additional = numeric(n_states),
      reach20 = numeric(n_states),
      reach30 = numeric(n_states)
    )

    for (entry_n in n_pos_min:max_interim_pos) {
      library_row <- positive_library[
        positive_library$entry_n == entry_n,
        ,
        drop = FALSE
      ]
      P_index <- library_row$P_index
      bP <- make_positive_boundary_vector(P_candidates, P_index)
      state_rows <- which(positive_states$n_pos == entry_n)
      entry_vector <- numeric(n_states)
      entry_vector[state_rows] <- enter[state_rows]

      detail <- compute_positive_path_candidate(
        bP = bP,
        theta_pos = theta_pos,
        enter_vector = entry_vector
      )

      for (field in names(combined)) {
        combined[[field]][state_rows] <- detail[[field]][state_rows]
      }
    }

    success[[format(theta_pos, trim = TRUE)]] <- combined$success
    details[[format(theta_pos, trim = TRUE)]] <- combined
  }

  list(
    enter = enter,
    positive_futility_initial = positive_futility_initial,
    success = success,
    details = details,
    positive_library = positive_library
  )
}

exact_metrics_entry_specific_positive <- function(
    component,
    comparator_positive_kernels,
    A_index
) {
  w1 <- component$W1[A_index, ]
  w2 <- component$W2[A_index, ]
  w <- w1 + w2
  scenario <- component$scenario
  theta_key <- format(scenario$theta_pos, trim = TRUE)

  success_vector <- comparator_positive_kernels$success[[theta_key]]
  detail <- comparator_positive_kernels$details[[theta_key]]
  enter_vector <- comparator_positive_kernels$enter
  initial_futility_vector <-
    comparator_positive_kernels$positive_futility_initial

  if (is.null(success_vector) || is.null(detail)) {
    stop("No entry-specific positive BOP2 kernel for theta_pos = ",
         scenario$theta_pos)
  }

  low_state <- positive_states$n_pos < n_pos_min
  assessed_state <- !low_state
  n_state <- positive_states$n_pos

  prn_all <- component$prn_all[A_index]
  prn_positive <- sum(w * success_vector)
  pet_low <- sum(w[low_state])
  pet_positive_initial <- sum(w * initial_futility_vector)
  pet_positive_look20 <- sum(w * detail$stop20)
  pet_positive_look30 <- sum(w * detail$stop30)
  pet_positive <- pet_positive_initial + pet_positive_look20 +
    pet_positive_look30
  enter <- sum(w * enter_vector)

  if (abs(pet_low + pet_positive_initial + enter - sum(w)) > 1e-9) {
    stop("Initial comparator enrichment-path probabilities do not sum correctly.")
  }
  post_partition <- sum(
    w * (detail$success + detail$final_no_claim +
           detail$stop20 + detail$stop30)
  )
  if (abs(post_partition - enter) > 1e-9) {
    stop("Comparator post-enrichment path probabilities do not sum correctly.")
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
    PET_any = pet_low + pet_positive,
    PET_all = sum(w),
    PET_positive = pet_positive,
    PET_positive_initial = pet_positive_initial,
    PET_positive_look20 = pet_positive_look20,
    PET_positive_look30 = pet_positive_look30,
    PET_low_biomarker_count = pet_low,
    Reassess_positive = sum(w[assessed_state]),
    Enter_enrichment = enter,
    Enter_enrichment_look1 = sum(w1 * enter_vector),
    Enter_enrichment_look2 = sum(w2 * enter_vector),
    Reach_positive_look20 = sum(w * detail$reach20),
    Reach_positive_look30 = sum(w * detail$reach30),
    Mean_total_sample_size = mean_total_sample_size,
    Mean_positive_sample_size = mean_positive_sample_size
  )
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
  enter_vector <- positive_kernels$enter[, P_index]
  initial_futility_vector <-
    positive_kernels$positive_futility_initial[, P_index]
  detail <- get_positive_detail(
    P_candidates = P_candidates,
    positive_kernels = positive_kernels,
    P_index = P_index,
    theta_pos = scenario$theta_pos
  )
  if (max(abs(success_vector - detail$success)) > 1e-10) {
    stop("Cached and detailed biomarker-positive success probabilities differ.")
  }

  low_state <- positive_states$n_pos < n_pos_min
  assessed_state <- !low_state
  n_state <- positive_states$n_pos

  prn_all <- component$prn_all[A_index]
  prn_positive <- sum(w * success_vector)
  pet_low <- sum(w[low_state])
  pet_positive_initial <- sum(w * initial_futility_vector)
  pet_positive_look20 <- sum(w * detail$stop20)
  pet_positive_look30 <- sum(w * detail$stop30)
  pet_positive <- pet_positive_initial + pet_positive_look20 +
    pet_positive_look30
  enter <- sum(w * enter_vector)

  if (abs(pet_low + pet_positive_initial + enter - sum(w)) > 1e-9) {
    stop("Initial enrichment-path probabilities do not sum correctly.")
  }
  post_partition <- sum(
    w * (detail$success + detail$final_no_claim +
           detail$stop20 + detail$stop30)
  )
  if (abs(post_partition - enter) > 1e-9) {
    stop("Post-enrichment path probabilities do not sum correctly.")
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
    PET_any = pet_low + pet_positive,
    PET_all = sum(w),
    PET_positive = pet_positive,
    PET_positive_initial = pet_positive_initial,
    PET_positive_look20 = pet_positive_look20,
    PET_positive_look30 = pet_positive_look30,
    PET_low_biomarker_count = pet_low,
    Reassess_positive = sum(w[assessed_state]),
    Enter_enrichment = enter,
    Enter_enrichment_look1 = sum(w1 * enter_vector),
    Enter_enrichment_look2 = sum(w2 * enter_vector),
    Reach_positive_look20 = sum(w * detail$reach20),
    Reach_positive_look30 = sum(w * detail$reach30),
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

  # Selection rule for the proposed method, as specified in the manuscript:
  # maximize average PRN-any under the calibration alternative, then maximize
  # the minimum PRN-any over prevalence, and use the larger max global type I
  # error only as the final tie breaker. This is intentionally different from
  # the standalone conventional-BOP2 selection rule below.
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

# Exact rejection probability of a conventional single-population BOP2 design
# with a fixed sequence of analyses.
standalone_bop2_rejection <- function(boundaries, looks, theta) {
  active <- c(1)
  previous_n <- 0L

  for (j in seq_along(looks)) {
    n_now <- looks[j]
    increment_n <- n_now - previous_n
    increment_pmf <- dbinom(0:increment_n, increment_n, theta)
    distribution_now <- numeric(n_now + 1L)

    for (x0 in 0:previous_n) {
      if (active[x0 + 1L] == 0) next
      distribution_now[(x0 + 0:increment_n) + 1L] <-
        distribution_now[(x0 + 0:increment_n) + 1L] +
        active[x0 + 1L] * increment_pmf
    }

    distribution_now[(0:n_now) <= boundaries[j]] <- 0
    active <- distribution_now
    previous_n <- n_now
  }

  sum(active)
}

select_standalone_bop2 <- function(candidates, looks, boundary_columns,
                                   theta_null, theta_alt, label) {
  type1 <- numeric(nrow(candidates))
  power <- numeric(nrow(candidates))

  for (i in seq_len(nrow(candidates))) {
    boundaries <- as.integer(unlist(
      candidates[i, boundary_columns, drop = FALSE],
      use.names = FALSE
    ))
    type1[i] <- standalone_bop2_rejection(boundaries, looks, theta_null)
    power[i] <- standalone_bop2_rejection(boundaries, looks, theta_alt)
  }

  feasible <- which(type1 <= alpha + 1e-12)
  if (length(feasible) == 0L) {
    stop("No feasible standalone BOP2 candidate for ", label)
  }

  # Zhou et al. (2017): select the attainable type I error closest to, and not
  # above, alpha. Among designs attaining that type I error, maximize power.
  tol <- 1e-12
  closest_type1 <- max(type1[feasible])
  best <- feasible[abs(type1[feasible] - closest_type1) <= tol]
  maximum_power <- max(power[best])
  best <- best[abs(power[best] - maximum_power) <= tol]
  chosen <- best[1L]

  list(
    index = chosen,
    type1 = type1[chosen],
    power = power[chosen],
    all_type1 = type1,
    all_power = power,
    looks = looks
  )
}

build_entry_specific_positive_bop2_library <- function(P_candidates) {
  rows <- vector("list", max_interim_pos - n_pos_min + 1L)
  row_index <- 1L

  for (entry_n in n_pos_min:max_interim_pos) {
    looks <- sort(unique(c(
      entry_n,
      post_enrichment_n[post_enrichment_n > entry_n],
      N_pos
    )))
    selected <- select_standalone_bop2(
      candidates = P_candidates,
      looks = looks,
      boundary_columns = paste0("b_", looks),
      theta_null = p0,
      theta_alt = p1,
      label = paste0("positive BOP2 entered at n = ", entry_n)
    )
    candidate <- P_candidates[selected$index, , drop = FALSE]

    rows[[row_index]] <- data.frame(
      entry_n = entry_n,
      analysis_schedule = paste(looks, collapse = ","),
      P_index = selected$index,
      candidate_id = candidate$candidate_id,
      lambda = candidate$lambda,
      gamma = candidate$gamma,
      standalone_type1 = selected$type1,
      standalone_power = selected$power,
      boundary_at_entry = as.integer(candidate[[paste0("b_", entry_n)]]),
      boundary_20 = if (20L %in% looks) as.integer(candidate[["b_20"]]) else NA_integer_,
      boundary_30 = if (30L %in% looks) as.integer(candidate[["b_30"]]) else NA_integer_,
      boundary_40 = as.integer(candidate[["b_40"]]),
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }

  do.call(rbind, rows)
}

## -----------------------------------------------------------------------------
## 3. Candidate boundary tables
## -----------------------------------------------------------------------------

message("Generating proposed-design boundary candidates...")
A_sample_sizes <- c(interim_n, N_A)
P_sample_sizes <- sort(unique(c(
  n_pos_min:max_interim_pos, post_enrichment_n, N_pos
)))

A_candidates <- make_boundary_candidates(
  sample_sizes = A_sample_sizes,
  max_sample_size = N_A,
  lambda_values = proposed_lambda_grid,
  gamma_values = proposed_gamma_grid,
  prefix = "proposed all-comer"
)
P_candidates <- make_boundary_candidates(
  sample_sizes = P_sample_sizes,
  max_sample_size = N_pos,
  lambda_values = proposed_lambda_grid,
  gamma_values = proposed_gamma_grid,
  prefix = "proposed biomarker-positive"
)

message("Generating standalone conventional-BOP2 boundary candidates...")
A_bop2_candidates <- make_boundary_candidates(
  sample_sizes = A_sample_sizes,
  max_sample_size = N_A,
  lambda_values = bop2_lambda_grid,
  gamma_values = bop2_gamma_grid,
  prefix = "standalone all-comer BOP2"
)
P_bop2_candidates <- make_boundary_candidates(
  sample_sizes = P_sample_sizes,
  max_sample_size = N_pos,
  lambda_values = bop2_lambda_grid,
  gamma_values = bop2_gamma_grid,
  prefix = "standalone positive BOP2"
)

message(
  "Unique proposed candidates: all-comer = ", nrow(A_candidates),
  "; positive = ", nrow(P_candidates),
  "; pairs = ",
  format(nrow(A_candidates) * nrow(P_candidates), big.mark = ",")
)
message(
  "Unique standalone-BOP2 candidates: all-comer = ",
  nrow(A_bop2_candidates), "; positive = ", nrow(P_bop2_candidates)
)

write.csv(A_candidates,
          file.path(output_dir, "proposed_allcomer_boundary_candidates.csv"),
          row.names = FALSE)
write.csv(P_candidates,
          file.path(output_dir, "proposed_positive_boundary_candidates.csv"),
          row.names = FALSE)
write.csv(A_bop2_candidates,
          file.path(output_dir, "bop2_allcomer_boundary_candidates.csv"),
          row.names = FALSE)
write.csv(P_bop2_candidates,
          file.path(output_dir, "bop2_positive_boundary_candidates.csv"),
          row.names = FALSE)

candidate_counts <- data.frame(
  component = c(
    "proposed all-comer",
    "proposed biomarker-positive",
    "proposed candidate pairs",
    "standalone all-comer BOP2",
    "standalone positive BOP2"
  ),
  count = c(
    nrow(A_candidates),
    nrow(P_candidates),
    nrow(A_candidates) * nrow(P_candidates),
    nrow(A_bop2_candidates),
    nrow(P_bop2_candidates)
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
preferred_global <- feasible_global &
  max_global_type1 >= alpha_lower - 1e-12

if (!any(feasible_global)) {
  stop("No candidate satisfies the global type I error upper constraint.")
}
selected_from_near_nominal_band <- any(preferred_global)
selection_set <- if (selected_from_near_nominal_band) {
  preferred_global
} else {
  warning("No candidate in the near-nominal band; using all candidates below alpha.")
  feasible_global
}

selected_global <- select_candidate_pair(
  feasible = selection_set,
  avg_power = avg_power,
  min_power = min_power,
  type1_metric = max_global_type1,
  label = "globally calibrated design with a near-nominal preference"
)

message(
  "Selected proposed design; near-nominal band used = ",
  selected_from_near_nominal_band,
  "; maximum point-null global type I error = ", selected_global$type1_metric
)

proposed_calibration_diagnostics <- data.frame(
  quantity = c(
    "number of candidate pairs",
    "number satisfying max global type I <= alpha",
    "number satisfying alpha_lower <= max global type I <= alpha",
    "selected max global type I",
    "selected average calibration PRN-any",
    "selected minimum calibration PRN-any"
  ),
  value = c(
    length(max_global_type1),
    sum(feasible_global),
    sum(preferred_global),
    selected_global$type1_metric,
    selected_global$average_power,
    selected_global$minimum_power
  ),
  stringsAsFactors = FALSE
)
write.csv(
  proposed_calibration_diagnostics,
  file.path(output_dir, "proposed_calibration_diagnostics.csv"),
  row.names = FALSE
)

# Conventional BOP2 comparator.
# 1. The all-comer BOP2 is calibrated independently for analyses at 20, 30,
#    and 40 patients.
# 2. The positive-population BOP2 is calibrated independently for every
#    attainable first-assessment size m. Its analysis schedule is
#    {m, future scheduled looks among 20 and 30, 40}.
# 3. The selected boundaries are then embedded in the enrichment procedure.
# No criterion involving the combined adaptive procedure is used in either
# standalone BOP2 calibration.
standalone_A <- select_standalone_bop2(
  candidates = A_bop2_candidates,
  looks = A_sample_sizes,
  boundary_columns = paste0("b_", A_sample_sizes),
  theta_null = p0,
  theta_alt = p1,
  label = "all-comer BOP2"
)

componentwise_allcomer_boundary <- data.frame(
  analysis_n = A_sample_sizes,
  futility_or_no_claim_boundary = as.integer(unlist(
    A_bop2_candidates[
      standalone_A$index,
      paste0("b_", A_sample_sizes),
      drop = FALSE
    ],
    use.names = FALSE
  )),
  is_final = A_sample_sizes == N_A,
  standalone_type1 = standalone_A$type1,
  standalone_power = standalone_A$power,
  stringsAsFactors = FALSE
)
write.csv(
  componentwise_allcomer_boundary,
  file.path(output_dir, "componentwise_allcomer_bop2_boundary.csv"),
  row.names = FALSE
)

positive_bop2_library <-
  build_entry_specific_positive_bop2_library(P_bop2_candidates)
write.csv(
  positive_bop2_library,
  file.path(output_dir, "componentwise_positive_bop2_library.csv"),
  row.names = FALSE
)

comparator_positive_kernels <- make_entry_specific_positive_kernels(
  P_candidates = P_bop2_candidates,
  positive_library = positive_bop2_library,
  theta_pos_values = positive_theta_values
)

selected_componentwise <- list(
  A_index = standalone_A$index,
  standalone_A_type1 = standalone_A$type1,
  standalone_A_power = standalone_A$power
)

selected_A_global <- A_candidates[selected_global$A_index, , drop = FALSE]
selected_P_global <- P_candidates[selected_global$P_index, , drop = FALSE]
selected_A_comparator <-
  A_bop2_candidates[standalone_A$index, , drop = FALSE]

message("Constructing exact all-comer path components for the selected ",
        "standalone all-comer BOP2 comparator...")
comparator_components <- parallel_lapply(
  scenario_rows,
  function(s) build_allcomer_components(s[1L, ], selected_A_comparator),
  n_workers = n_cores
)

selected_designs <- rbind(
  data.frame(
    Method = "Globally calibrated",
    A_candidate_id = selected_A_global$candidate_id,
    P_candidate_id = selected_P_global$candidate_id,
    Positive_component = "Single jointly calibrated boundary function",
    lambda_A = selected_A_global$lambda,
    gamma_A = selected_A_global$gamma,
    lambda_positive = selected_P_global$lambda,
    gamma_positive = selected_P_global$gamma,
    max_global_type1 = selected_global$type1_metric,
    max_allcomer_type1 = max_all_type1[selected_global$A_index],
    max_positive_type1 =
      max_positive_type1[selected_global$A_index, selected_global$P_index],
    average_power = selected_global$average_power,
    minimum_power = selected_global$minimum_power,
    standalone_allcomer_type1 = NA_real_,
    standalone_allcomer_power = NA_real_,
    alpha_lower = alpha_lower,
    selected_from_near_nominal_band = selected_from_near_nominal_band,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Method = "Independent conventional BOP2 components",
    A_candidate_id = selected_A_comparator$candidate_id,
    P_candidate_id = NA_integer_,
    Positive_component = "Entry-specific conventional BOP2 library",
    lambda_A = selected_A_comparator$lambda,
    gamma_A = selected_A_comparator$gamma,
    lambda_positive = NA_real_,
    gamma_positive = NA_real_,
    max_global_type1 = NA_real_,
    max_allcomer_type1 = standalone_A$type1,
    max_positive_type1 = NA_real_,
    average_power = NA_real_,
    minimum_power = NA_real_,
    standalone_allcomer_type1 = standalone_A$type1,
    standalone_allcomer_power = standalone_A$power,
    alpha_lower = NA_real_,
    selected_from_near_nominal_band = NA,
    stringsAsFactors = FALSE
  )
)

if (any(abs(selected_designs$gamma_A[1L] -
            max(proposed_gamma_grid)) < 1e-12, na.rm = TRUE) ||
    any(abs(selected_designs$gamma_positive[1L] -
            max(proposed_gamma_grid)) < 1e-12, na.rm = TRUE) ||
    any(abs(selected_designs$gamma_A[2L] -
            max(bop2_gamma_grid)) < 1e-12, na.rm = TRUE) ||
    any(abs(positive_bop2_library$gamma -
            max(bop2_gamma_grid)) < 1e-12, na.rm = TRUE)) {
  warning(
    "At least one selected gamma value is at its upper search limit. ",
    "Consider extending the corresponding gamma grid before manuscript use."
  )
}

write.csv(
  selected_designs,
  file.path(output_dir, "selected_designs.csv"),
  row.names = FALSE
)
print(selected_designs)

# The calibration snapshot is saved after comparator evaluation below.

## -----------------------------------------------------------------------------
## 6. Exact operating characteristics for selected designs
## -----------------------------------------------------------------------------

exact_rows <- vector("list", length(components) * 2L)
out_idx <- 1L
for (s_idx in seq_along(components)) {
  vals_proposed <- exact_metrics_selected(
    component = components[[s_idx]],
    positive_kernels = positive_kernels,
    A_index = selected_global$A_index,
    P_index = selected_global$P_index
  )
  exact_rows[[out_idx]] <- cbind(
    all_scenarios[s_idx, ],
    Method = "Globally calibrated",
    as.data.frame(as.list(vals_proposed), stringsAsFactors = FALSE)
  )
  out_idx <- out_idx + 1L

  vals_comparator <- exact_metrics_entry_specific_positive(
    component = comparator_components[[s_idx]],
    comparator_positive_kernels = comparator_positive_kernels,
    A_index = 1L
  )
  exact_rows[[out_idx]] <- cbind(
    all_scenarios[s_idx, ],
    Method = "Independent conventional BOP2 components",
    as.data.frame(as.list(vals_comparator), stringsAsFactors = FALSE)
  )
  out_idx <- out_idx + 1L
}

exact_results <- do.call(rbind, exact_rows)
rownames(exact_results) <- NULL
exact_results$p0 <- p0
exact_results$p1 <- p1
exact_results$theta_all <- with(exact_results,
  pi * theta_pos + (1 - pi) * theta_neg)
exact_results$H0_all_true <- exact_results$theta_all <= p0 + 1e-12
exact_results$H0_positive_true <- exact_results$theta_pos <= p0 + 1e-12
exact_results$FWER_at_scenario <- with(exact_results,
  as.numeric(H0_all_true) * PRN_all + as.numeric(H0_positive_true) * PRN_positive)
exact_results$P_valid_claim <- with(exact_results,
  as.numeric(!H0_all_true) * PRN_all + as.numeric(!H0_positive_true) * PRN_positive)
stopifnot(max(abs(exact_results$PRN_any - exact_results$FWER_at_scenario -
                  exact_results$P_valid_claim)) < 1e-10)
probability_columns <- c("PRN_any", "PRN_all", "PRN_positive", "PET_any",
  "PET_positive", "PET_low_biomarker_count", "Enter_enrichment")
stopifnot(all(as.matrix(exact_results[, probability_columns]) >= -1e-10),
          all(as.matrix(exact_results[, probability_columns]) <= 1 + 1e-10),
          all(exact_results$Mean_total_sample_size >= 20 - 1e-8),
          all(exact_results$Mean_total_sample_size <= 60 + 1e-8))


# Numerical consistency checks for the mutually exclusive rejection paths.
if (max(abs(
  exact_results$PRN_any -
    exact_results$PRN_all -
    exact_results$PRN_positive
)) > 1e-10) {
  stop("PRN-any does not equal PRN-all + PRN-positive.")
}

proposed_null_check <- exact_results[
  exact_results$Method == "Globally calibrated" &
    exact_results$scenario_type == "Null",
  ,
  drop = FALSE
]
if (max(proposed_null_check$PRN_any) > alpha + 1e-10) {
  stop("The proposed design exceeds the point-null global type I error constraint.")
}
if (selected_from_near_nominal_band &&
    max(proposed_null_check$PRN_any) < alpha_lower - 1e-10) {
  stop("Selected design is inconsistent with the near-nominal band flag.")
}

write.csv(
  exact_results,
  file.path(output_dir, "exact_operating_characteristics.csv"),
  row.names = FALSE
)

# Fill the combined comparator operating characteristics only after its two
# standalone components have been embedded in the adaptive enrichment design.
comparator_null <- exact_results[
  exact_results$Method == "Independent conventional BOP2 components" &
    exact_results$scenario_type == "Null",
  ,
  drop = FALSE
]
comparator_calibration_alt <- exact_results[
  exact_results$Method == "Independent conventional BOP2 components" &
    exact_results$used_for_calibration,
  ,
  drop = FALSE
]

selected_designs$max_global_type1[2L] <- max(comparator_null$PRN_any)
selected_designs$max_positive_type1[2L] <- max(comparator_null$PRN_positive)
selected_designs$average_power[2L] <- mean(comparator_calibration_alt$PRN_any)
selected_designs$minimum_power[2L] <- min(comparator_calibration_alt$PRN_any)
write.csv(
  selected_designs,
  file.path(output_dir, "selected_designs.csv"),
  row.names = FALSE
)

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
  "\\caption{Futility boundaries of the proposed globally calibrated enrichment BOP2 design.}",
  paste0("\\label{tab:proposed_boundaries_", gsub("[^A-Za-z0-9]", "", case_id), "}"),
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
  "\\footnotesize The all-comer boundaries apply at 20 and 30 patients and at the final look of 40 patients. The biomarker-positive boundaries from 10 to 30 patients apply when enrichment is first considered; the boundaries at 20 and 30 patients also apply at post-enrichment interim looks when those looks occur after enrichment begins. The boundary at 40 patients is the final biomarker-positive boundary. At a final look, crossing the boundary means that efficacy is not declared.",
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
    "Analysis type", "branching rule", "p0", "p1", "alpha", "alpha_lower",
    "calibration alternative",
    "all-comer interim looks", "post-enrichment positive looks",
    "N_A", "N_positive", "n_positive_min", "prior_a", "prior_b",
    "prevalence grid", "theta_positive grid", "theta_negative grid",
    "proposed lambda grid", "proposed gamma grid",
    "standalone BOP2 lambda grid", "standalone BOP2 gamma grid",
    "detected physical cores",
    "number of workers used", "parallel backend"
  ),
  Value = c(
    "Exact recursive enumeration", "No bridging", p0, p1, alpha, alpha_lower,
    paste0("theta_positive = ", p1, "; theta_negative = ", p0),
    paste(interim_n, collapse = ", "),
    paste(post_enrichment_n, collapse = ", "),
    N_A, N_pos, n_pos_min, prior_a, prior_b,
    paste(prevalence_grid, collapse = ", "),
    paste(theta_pos_grid, collapse = ", "),
    paste(theta_neg_grid, collapse = ", "),
    paste0(min(proposed_lambda_grid), " to ", max(proposed_lambda_grid),
           " by ", unique(round(diff(proposed_lambda_grid), 12))),
    paste0(min(proposed_gamma_grid), " to ", max(proposed_gamma_grid),
           " by ", unique(round(diff(proposed_gamma_grid), 12))),
    paste0(min(bop2_lambda_grid), " to ", max(bop2_lambda_grid),
           " by ", unique(round(diff(bop2_lambda_grid), 12))),
    paste0(min(bop2_gamma_grid), " to ", max(bop2_gamma_grid),
           " by ", unique(round(diff(bop2_gamma_grid), 12))),
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
    Boundary_positive = positive_boundary_table,
    Proposed_candidate_counts = candidate_counts,
    Proposed_calibration = proposed_calibration_diagnostics,
    Comparator_allcomer_BOP2 = componentwise_allcomer_boundary,
    Comparator_positive_BOP2 = positive_bop2_library
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

selected_designs$p0 <- p0
selected_designs$p1 <- p1
write.csv(selected_designs, file.path(output_dir, "selected_designs.csv"), row.names = FALSE)
saveRDS(
  list(
    settings = list(
      p0 = p0, p1 = p1, alpha = alpha, alpha_lower = alpha_lower,
      interim_n = interim_n, N_A = N_A, N_pos = N_pos,
      n_pos_min = n_pos_min, post_enrichment_n = post_enrichment_n,
      prior_a = prior_a, prior_b = prior_b,
      prevalence_grid = prevalence_grid,
      theta_pos_grid = theta_pos_grid,
      theta_neg_grid = theta_neg_grid,
      proposed_lambda_grid = proposed_lambda_grid,
      proposed_gamma_grid = proposed_gamma_grid,
      bop2_lambda_grid = bop2_lambda_grid,
      bop2_gamma_grid = bop2_gamma_grid
    ),
    A_candidates = A_candidates,
    P_candidates = P_candidates,
    A_bop2_candidates = A_bop2_candidates,
    P_bop2_candidates = P_bop2_candidates,
    selected_global = selected_global,
    selected_componentwise = selected_componentwise,
    positive_bop2_library = positive_bop2_library,
    selected_designs = selected_designs
  ),
  file.path(output_dir, "calibration_objects.rds")
)

# Regression check against the supplied manuscript results, not against a new
# Monte Carlo simulation. Skipped only in the explicitly reduced quick mode.
if (!quick_test && abs(gamma_max - 8) < 1e-12 &&
    abs(p0 - 0.20) < 1e-12 && abs(p1 - 0.40) < 1e-12) {
  expected <- c(0.0869797547519798, 0.0959043271860930, 0.0989454278278690,
                0.0965970582464720, 0.0910772077224320)
  got <- proposed_null_check[order(proposed_null_check$pi), ]
  stopifnot(nrow(got) == 5L, max(abs(got$PRN_any - expected)) < 1e-10,
    identical(as.integer(unlist(selected_A[1, c("b_20", "b_30", "b_40")])), c(4L, 8L, 11L)))
  message("Baseline regression check: PASS.")
}

# Report boundary behavior that may require clinical interpretation.
positive_entry_boundaries <- as.integer(unlist(
  selected_P[1L, paste0("b_", n_pos_min:max_interim_pos)]
))
no_futility_n <- (n_pos_min:max_interim_pos)[positive_entry_boundaries < 0L]
near_gamma_limit <- selected_A$gamma >= 0.95 * max(proposed_gamma_grid) ||
                    selected_P$gamma >= 0.95 * max(proposed_gamma_grid)
if (near_gamma_limit) {
  warning("A selected proposed gamma is close to the numerical search limit. ",
    "Consider a separate run with BOP2_PROPOSED_GAMMA_MAX=12.")
}
if (length(no_futility_n) > 0L) {
  warning("Positive-path boundary -1 occurs at n = ",
    paste(no_futility_n, collapse=","),
    ". It means no response-count futility stopping at those sample sizes, ",
    "not a missing value. The minimum positive count remains required.")
}
writeLines(c(
  paste0("Proposed gamma grid maximum: ", max(proposed_gamma_grid)),
  paste0("Selected gamma_A: ", selected_A$gamma),
  paste0("Selected gamma_positive: ", selected_P$gamma),
  paste0("Near gamma grid limit: ", near_gamma_limit),
  paste0("Positive n with boundary -1: ", paste(no_futility_n, collapse=",")),
  "Boundary -1 means no futility stopping; it does not mean missing data.",
  "No extra zero-response clinical restriction was added to the main design.",
  "Exact enumeration eliminates Monte Carlo error, not grid-search limitations."),
  file.path(output_dir, "boundary_interpretation.txt"))

# Independent single-population calculation validates all-comer PRN for every
# evaluation scenario without the three-dimensional state enumeration.
for (method in unique(exact_results$Method)) {
  d <- exact_results[exact_results$Method == method, ]
  if (method == "Globally calibrated") {
    ba <- as.integer(unlist(selected_A[1L, paste0("b_", A_sample_sizes)]))
  } else {
    ba <- componentwise_allcomer_boundary$futility_or_no_claim_boundary
  }
  check <- vapply(d$theta_all, function(th) {
    standalone_bop2_rejection(ba, A_sample_sizes, th)
  }, numeric(1))
  stopifnot(max(abs(check - d$PRN_all)) < 1e-10)
}
writeLines(c("PASS: probability partitions and bounds.",
             "PASS: all-comer PRN agrees with independent one-population recursion.",
             paste0("Near-nominal band used: ", selected_from_near_nominal_band),
             "Scope: prespecified point complete null over the prevalence grid.",
             "This analysis does not establish strong FWER control."),
           file.path(output_dir, "validation_checks.txt"))
invisible(list(results = exact_results, selected = selected_designs,
  proposed_null_max = max(proposed_null_check$PRN_any), output_dir = output_dir))

}

run_r1m1 <- function() {
  old_options <- options(stringsAsFactors = FALSE)
  on.exit(options(old_options), add = TRUE)
  include_baseline <- !identical(Sys.getenv("BOP2_INCLUDE_BASELINE", unset="1"), "0")
  default_dir <- paste0("R1M1_results_", format(Sys.time(), "%Y%m%d_%H%M%S"),
                        "_", Sys.getpid())
  output_root <- path.expand(Sys.getenv("BOP2_R1M1_OUTPUT_DIR", unset=default_dir))
  if (file.exists(output_root)) {
    stop("Output path already exists. Choose a NEW directory: ", output_root)
  }
  if (!dir.create(output_root, recursive=TRUE)) stop("Cannot create output directory.")
  output_root <- normalizePath(output_root, mustWork=TRUE)
  settings <- data.frame(p0=c(0.10, 0.30), p1=c(0.30, 0.50))
  if (include_baseline) settings <- rbind(data.frame(p0=0.20, p1=0.40), settings)
  write.csv(settings, file.path(output_root, "benchmark_settings.csv"), row.names=FALSE)
  started <- Sys.time()
  runs <- vector("list", nrow(settings))
  for (i in seq_len(nrow(settings))) {
    message("\n=== R1-M1 setting ", i, "/", nrow(settings), ": p0=",
      settings$p0[i], ", p1=", settings$p1[i], " ===")
    runs[[i]] <- run_bop2_setting(settings$p0[i], settings$p1[i], output_root)
    invisible(gc())
  }
  combined <- do.call(rbind, lapply(runs, function(x) x$results))
  selected <- do.call(rbind, lapply(runs, function(x) x$selected))
  rownames(combined) <- NULL
  rownames(selected) <- NULL
  write.csv(combined, file.path(output_root, "ALL_exact_operating_characteristics.csv"),
            row.names=FALSE)
  write.csv(combined[combined$scenario_type=="Null", ],
            file.path(output_root, "ALL_exact_type1_error.csv"), row.names=FALSE)
  write.csv(combined[combined$scenario_type=="Alternative", ],
            file.path(output_root, "ALL_exact_alternatives.csv"), row.names=FALSE)
  write.csv(selected, file.path(output_root, "ALL_selected_designs.csv"), row.names=FALSE)

  summary_rows <- list()
  for (i in seq_len(nrow(settings))) {
    for (method in unique(combined$Method)) {
      d <- combined[abs(combined$p0-settings$p0[i])<1e-12 &
                    abs(combined$p1-settings$p1[i])<1e-12 & combined$Method==method, ]
      dn <- d[d$scenario_type=="Null", ]
      da <- d[d$scenario_type=="Alternative", ]
      dc <- d[d$used_for_calibration, ]
      summary_rows[[length(summary_rows)+1L]] <- data.frame(
        p0=settings$p0[i], p1=settings$p1[i], Method=method,
        min_point_null_PRN_any=min(dn$PRN_any), max_point_null_PRN_any=max(dn$PRN_any),
        mean_working_alternative_PRN_any=mean(dc$PRN_any),
        min_working_alternative_PRN_any=min(dc$PRN_any),
        min_evaluation_PRN_any=min(da$PRN_any), max_evaluation_PRN_any=max(da$PRN_any),
        min_null_mean_N=min(dn$Mean_total_sample_size),
        max_null_mean_N=max(dn$Mean_total_sample_size))
    }
  }
  summary <- do.call(rbind, summary_rows)
  write.csv(summary, file.path(output_root, "ALL_summary.csv"), row.names=FALSE)
  capture.output(sessionInfo(), file=file.path(output_root, "sessionInfo.txt"))
  writeLines(c("Analysis completed.", paste0("Started: ", started),
    paste0("Finished: ", Sys.time()),
    paste0("Quick mode: ", identical(Sys.getenv("BOP2_QUICK"), "1")),
    "Probabilities in CSV files are on the 0-1 scale; multiply by 100 for percent.",
    "No Monte Carlo standard errors apply to exact finite-state probabilities.",
    "Finite candidate grids are used; no claim of a continuous global optimum.",
    "FWER_at_scenario is descriptive and is not a strong-control guarantee.",
    "Independent BOP2 comparator retains the manuscript's unadjusted construction."),
    file.path(output_root, "RUN_COMPLETED.txt"))
  print(summary)
  message("\nAll R1-M1 results: ", output_root)
  invisible(list(output_dir=output_root, summary=summary))
}

run_r1m1()
