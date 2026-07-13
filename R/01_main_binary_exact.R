# Binary-endpoint analysis for the main manuscript.
# Uses exact recursive enumeration and relative repository paths.

options(stringsAsFactors = FALSE)

# User settings

p0 <- 0.20
p1 <- 0.40
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

theta_pos_grid <- c(p1, p1 + 0.10, p1 + 0.20)
theta_neg_grid <- c(p0, p0 - 0.10)

lambda_grid <- seq(0.005, 1.000, by = 0.005)

# gamma is a positive shape parameter rather than a probability and therefore
# need not be bounded by one. Values above one allow the cutoff to remain
# relatively close to one early in the trial and decrease more sharply later.
# The upper limit of eight is a numerical search bound.
gamma_grid <- seq(0.050, 8.000, by = 0.050)

repo_root <- normalizePath(
  Sys.getenv("ENRICHMENT_BOP2_ROOT", unset = getwd()),
  mustWork = FALSE
)
result_root <- file.path(repo_root, "results", "main")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

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

legacy_output_files <- c(
  "monte_carlo_operating_characteristics.csv",
  "exact_vs_monte_carlo.csv",
  "simulation_scenarios.csv",
  "enrichment_bop2_results.xlsx"
)
unlink(file.path(output_dir, legacy_output_files), force = TRUE)

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

# Utility functions

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
  dim1 <- total_n + 1L
  dim2 <- total_n + 1L
  out_prob <- numeric(dim1 * dim2 * dim2)

  for (i in seq_len(nrow(dist1))) {
    n_pos <- dist1$n_pos[i] + dist2$n_pos
    x_pos <- dist1$x_pos[i] + dist2$x_pos
    x_neg <- dist1$x_neg[i] + dist2$x_neg

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

# Candidate boundary tables

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

# Scenarios and exact path components

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

# A single prespecified alternative configuration is used for calibration.
# Other response configurations are used only to evaluate the operating
# characteristics of the selected fixed design.
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

# Calibration of the proposed and comparator designs

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
use_alpha_lower_band <- any(preferred_global)
global_selection_set <- if (use_alpha_lower_band) {
  preferred_global
} else {
  feasible_global
}

feasible_componentwise <-
  matrix(max_all_type1 <= alpha + 1e-12, nrow = nA, ncol = nP) &
  max_positive_type1 <= alpha + 1e-12

if (use_alpha_lower_band) {
  message(
    "Selecting the proposed design among candidates with exact maximum global ",
    "type I error in [", alpha_lower, ", ", alpha, "]."
  )
} else {
  warning(
    "No proposed-design candidate had exact maximum global type I error in [",
    alpha_lower, ", ", alpha, "]. Selection used all candidates not exceeding ",
    alpha, "."
  )
}

selected_global <- select_candidate_pair(
  feasible = global_selection_set,
  avg_power = avg_power,
  min_power = min_power,
  type1_metric = max_global_type1,
  label = "globally calibrated design"
)
selected_global$used_alpha_lower_band <- use_alpha_lower_band

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
    alpha_lower = if (method == "Globally calibrated") alpha_lower else NA_real_,
    selected_from_alpha_lower_band = if (method == "Globally calibrated") {
      use_alpha_lower_band
    } else {
      NA
    },
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
      p0 = p0, p1 = p1, alpha = alpha, alpha_lower = alpha_lower,
      calibration_alternative = c(theta_pos = p1, theta_neg = p0),
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

# Exact operating characteristics for selected designs

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

# Table-I-style boundary table for the proposed design only

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

# Exact result summaries and session information

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
    "Analysis type", "p0", "p1", "alpha", "alpha_lower",
    "calibration alternative",
    "all-comer interim looks", "post-enrichment positive looks",
    "N_A", "N_positive", "n_positive_min", "prior_a", "prior_b",
    "prevalence grid", "theta_positive grid", "theta_negative grid",
    "lambda grid", "gamma grid", "detected physical cores",
    "number of workers used", "parallel backend"
  ),
  Value = c(
    "Exact recursive enumeration", p0, p1, alpha, alpha_lower,
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
