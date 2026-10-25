# Simulation-based type I error evaluation for complex categorical endpoints.
#
# Corrected design features:
#   * NO BRIDGING: if an all-comer futility boundary is crossed with fewer
#     than 10 biomarker-positive patients, the trial stops.
#   * Proposed design: the complete adaptive procedure is globally calibrated.
#   * Comparator: a conventional all-comer BOP2 and an entry-specific library
#     of conventional positive-population BOP2 designs are calibrated
#     independently before being combined. Each standalone BOP2 follows the
#     Zhou, Lee, and Yuan (2017) rule: among candidates with type I error not
#     exceeding alpha, first choose the type I error closest to alpha and then
#     maximize power.
#
# Results are written to:
#   results/supplement
#
# Environment variables:
#   BOP2_COMPLEX_QUICK=1       reduced code check
#   BOP2_COMPLEX_N_CALIB       calibration replicates (default 50000)
#   BOP2_COMPLEX_N_VALID       validation replicates (default 100000)
#   BOP2_N_CORES               parallel workers for standalone schedules

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## -----------------------------------------------------------------------------
## 1. User settings
## -----------------------------------------------------------------------------

alpha <- 0.10
alpha_lower <- 0.09

interim_n <- c(20L, 30L)
N_A <- 40L
N_positive <- 40L
n_positive_min <- 10L
post_enrichment_n <- c(20L, 30L)
prevalence_grid <- c(0.4, 0.5, 0.6, 0.7, 0.8)

# The supplemental complex-endpoint analysis is a type I error demonstration.
# As in the previous program, gamma is fixed for the proposed design and the
# two lambda values are calibrated over the complete adaptive procedure.
proposed_gamma_A <- 1
proposed_gamma_positive <- 1
proposed_lambda_grid <- seq(0.005, 1.000, by = 0.005)

# Conventional BOP2 components use the full two-parameter threshold family.
# The fine grid matches the corrected binary comparator program.
bop2_lambda_grid <- seq(0.001, 1.000, by = 0.0025)
bop2_gamma_grid <- seq(0.010, 10.000, by = 0.025)

# Symmetric Dirichlet prior with total prior sample size one.
dirichlet_prior <- rep(0.25, 4L)

n_calibration <- as.integer(Sys.getenv(
  "BOP2_COMPLEX_N_CALIB", unset = "50000"
))
n_validation <- as.integer(Sys.getenv(
  "BOP2_COMPLEX_N_VALID", unset = "100000"
))
base_seed <- 20260722L

# The proposed-design feasibility criterion uses a one-sided 95% Wilson upper
# confidence bound. Standalone conventional BOP2 selection follows the original
# BOP2 rule using the Monte Carlo point estimate, with independent validation
# reported afterward.
proposed_use_wilson_ucb <- TRUE

endpoint_settings <- list(
  nested = list(
    endpoint_id = "nested",
    endpoint_label = "Nested efficacy",
    categories = c("CR", "PR", "SD", "PD"),
    theta0 = c(0.15, 0.15, 0.30, 0.40),
    theta1 = c(0.25, 0.25, 0.25, 0.25),
    rule_label = "CR futility AND CR/PR futility"
  ),
  coprimary = list(
    endpoint_id = "coprimary",
    endpoint_label = "Co-primary efficacy",
    categories = c(
      "OR_and_EFS6", "OR_and_no_EFS6",
      "no_OR_and_EFS6", "no_OR_and_no_EFS6"
    ),
    theta0 = c(0.05, 0.05, 0.15, 0.75),
    theta1 = c(0.15, 0.15, 0.20, 0.50),
    rule_label = "ORR futility AND EFS6 futility"
  ),
  efficacy_toxicity = list(
    endpoint_id = "efficacy_toxicity",
    endpoint_label = "Efficacy and toxicity",
    categories = c(
      "toxicity_and_OR", "no_toxicity_and_OR",
      "toxicity_and_no_OR", "no_toxicity_and_no_OR"
    ),
    theta0 = c(0.15, 0.30, 0.15, 0.40),
    theta1 = c(0.18, 0.42, 0.02, 0.38),
    rule_label = "efficacy futility OR excessive toxicity"
  )
)

for (setting in endpoint_settings) {
  if (length(setting$theta0) != 4L || length(setting$theta1) != 4L ||
      any(setting$theta0 < 0) || any(setting$theta1 < 0) ||
      abs(sum(setting$theta0) - 1) > 1e-12 ||
      abs(sum(setting$theta1) - 1) > 1e-12) {
    stop("Invalid endpoint probability vector for ", setting$endpoint_id)
  }
}

repo_root <- normalizePath(
  Sys.getenv("ENRICHMENT_BOP2_ROOT", unset = getwd()),
  mustWork = FALSE
)
result_dir <- file.path(repo_root, "results", "supplement")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
default_cores <- max(1L, min(4L, detected_cores - 1L))
requested_cores <- suppressWarnings(as.integer(Sys.getenv(
  "BOP2_N_CORES", unset = as.character(default_cores)
)))
if (is.na(requested_cores) || requested_cores < 1L) {
  requested_cores <- default_cores
}
n_cores <- max(1L, min(requested_cores, detected_cores))

quick_test <- identical(Sys.getenv("BOP2_COMPLEX_QUICK"), "1")
if (quick_test) {
  message(
    "BOP2_COMPLEX_QUICK=1: reduced code check only; not manuscript results."
  )
  prevalence_grid <- c(0.4, 0.8)
  proposed_lambda_grid <- seq(0.05, 1.00, by = 0.05)
  bop2_lambda_grid <- seq(0.05, 1.00, by = 0.05)
  bop2_gamma_grid <- seq(0.25, 2.00, by = 0.25)
  n_calibration <- 2000L
  n_validation <- 5000L
  n_cores <- 1L
  result_dir <- file.path(result_dir, "quick_test")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(
  alpha_lower >= 0,
  alpha_lower <= alpha,
  N_A == N_positive,
  all(interim_n < N_A),
  n_positive_min <= min(interim_n),
  all(post_enrichment_n < N_positive),
  all(dirichlet_prior > 0),
  abs(sum(dirichlet_prior) - 1) < 1e-12,
  n_calibration > 0L,
  n_validation > 0L
)

message(
  "Complex-endpoint simulation: calibration = ", n_calibration,
  "; validation = ", n_validation,
  "; schedule workers = ", n_cores
)

## -----------------------------------------------------------------------------
## 2. General utilities
## -----------------------------------------------------------------------------

parallel_lapply <- function(X, FUN, ..., mc_cores = n_cores) {
  if (length(X) == 0L) {
    return(list())
  }
  mc_cores <- max(1L, min(mc_cores, length(X)))
  if (mc_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      X, FUN, ...,
      mc.cores = mc_cores,
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
      openxlsx::setColWidths(
        wb, sheet = sheet_name,
        cols = seq_len(ncol(sheets[[sheet_name]])), widths = "auto"
      )
    }
    openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
    return(TRUE)
  }
  warning(
    "Neither writexl nor openxlsx is installed. CSV files were created, ",
    "but no consolidated XLSX workbook was written."
  )
  FALSE
}

wilson_upper <- function(p, n, confidence = 0.95) {
  z <- qnorm(confidence)
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z / denominator * sqrt(
    p * (1 - p) / n + z^2 / (4 * n^2)
  )
  pmin(1, center + half)
}

wilson_interval <- function(p, n, confidence = 0.95) {
  z <- qnorm(1 - (1 - confidence) / 2)
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z / denominator * sqrt(
    p * (1 - p) / n + z^2 / (4 * n^2)
  )
  c(
    lower = max(0, center - half),
    upper = min(1, center + half)
  )
}

rmat_categorical <- function(n_row, n_col, prob) {
  breaks <- cumsum(prob)
  breaks <- breaks[-length(breaks)]
  u <- runif(n_row * n_col)
  matrix(
    findInterval(u, breaks) + 1L,
    nrow = n_row,
    ncol = n_col
  )
}

row_cumulative_indicator <- function(category_matrix, category) {
  out <- (category_matrix == category) * 1L
  if (ncol(out) > 1L) {
    for (j in 2:ncol(out)) {
      out[, j] <- out[, j] + out[, j - 1L]
    }
  }
  out
}

make_cumulative_counts <- function(Y) {
  lapply(1:4, function(k) row_cumulative_indicator(Y, k))
}

counts_from_cumulative <- function(cumulative, n) {
  out <- cbind(
    cumulative[[1L]][, n], cumulative[[2L]][, n],
    cumulative[[3L]][, n], cumulative[[4L]][, n]
  )
  colnames(out) <- paste0("x", 1:4)
  out
}

counts_at_look <- function(Y, look, Z = NULL) {
  idx <- seq_len(look)
  out <- matrix(0L, nrow = nrow(Y), ncol = 4L)
  if (is.null(Z)) {
    for (k in 1:4) {
      out[, k] <- rowSums(Y[, idx, drop = FALSE] == k)
    }
  } else {
    z_sub <- Z[, idx, drop = FALSE] == 1L
    for (k in 1:4) {
      out[, k] <- rowSums((Y[, idx, drop = FALSE] == k) & z_sub)
    }
  }
  colnames(out) <- paste0("x", 1:4)
  out
}

endpoint_event_counts <- function(counts, endpoint_id) {
  counts <- as.matrix(counts)
  if (endpoint_id == "nested") {
    return(list(
      event1 = counts[, 1L],
      event2 = counts[, 1L] + counts[, 2L]
    ))
  }
  if (endpoint_id == "coprimary") {
    return(list(
      event1 = counts[, 1L] + counts[, 2L],
      event2 = counts[, 1L] + counts[, 3L]
    ))
  }
  if (endpoint_id == "efficacy_toxicity") {
    return(list(
      event1 = counts[, 1L] + counts[, 2L],
      event2 = counts[, 1L] + counts[, 3L]
    ))
  }
  stop("Unknown endpoint_id: ", endpoint_id)
}

## -----------------------------------------------------------------------------
## 3. Posterior rules and integer boundary generation
## -----------------------------------------------------------------------------

posterior_stop_score <- function(counts, endpoint_id) {
  counts <- as.matrix(counts)
  n <- rowSums(counts)
  a <- dirichlet_prior

  if (endpoint_id == "nested") {
    p1 <- pbeta(
      0.15,
      shape1 = a[1L] + counts[, 1L],
      shape2 = sum(a[-1L]) + n - counts[, 1L]
    )
    crpr <- counts[, 1L] + counts[, 2L]
    p2 <- pbeta(
      0.30,
      shape1 = a[1L] + a[2L] + crpr,
      shape2 = a[3L] + a[4L] + n - crpr
    )
    return(pmin(p1, p2))
  }

  if (endpoint_id == "coprimary") {
    or_count <- counts[, 1L] + counts[, 2L]
    efs_count <- counts[, 1L] + counts[, 3L]
    p1 <- pbeta(
      0.10,
      shape1 = a[1L] + a[2L] + or_count,
      shape2 = a[3L] + a[4L] + n - or_count
    )
    p2 <- pbeta(
      0.20,
      shape1 = a[1L] + a[3L] + efs_count,
      shape2 = a[2L] + a[4L] + n - efs_count
    )
    return(pmin(p1, p2))
  }

  if (endpoint_id == "efficacy_toxicity") {
    response_count <- counts[, 1L] + counts[, 2L]
    toxicity_count <- counts[, 1L] + counts[, 3L]
    p_eff <- pbeta(
      0.45,
      shape1 = a[1L] + a[2L] + response_count,
      shape2 = a[3L] + a[4L] + n - response_count
    )
    p_tox <- pbeta(
      0.30,
      shape1 = a[1L] + a[3L] + toxicity_count,
      shape2 = a[2L] + a[4L] + n - toxicity_count,
      lower.tail = FALSE
    )
    return(pmax(p_eff, p_tox))
  }

  stop("Unknown endpoint_id: ", endpoint_id)
}

lower_boundary_from_cutoff <- function(
    n, phi, prior_event, prior_other, cutoff
) {
  x <- 0:n
  posterior_probability <- pbeta(
    phi,
    shape1 = prior_event + x,
    shape2 = prior_other + n - x
  )
  # posterior_probability is decreasing in x. The number of values not larger
  # than cutoff is obtained by findInterval; strict posterior > cutoff gives
  # the boundary below.
  as.integer(n - findInterval(cutoff, sort(posterior_probability)))
}

upper_boundary_from_cutoff <- function(
    n, phi, prior_event, prior_other, cutoff
) {
  x <- 0:n
  posterior_probability <- pbeta(
    phi,
    shape1 = prior_event + x,
    shape2 = prior_other + n - x,
    lower.tail = FALSE
  )
  # posterior_probability is increasing in x. The minimum x satisfying the
  # strict posterior > cutoff condition equals the number of values <= cutoff.
  as.integer(findInterval(cutoff, posterior_probability))
}

calculate_boundary_vectors <- function(setting, n, lambda, gamma) {
  cutoff <- 1 - lambda * (n / N_A)^gamma
  a <- dirichlet_prior

  if (setting$endpoint_id == "nested") {
    b1 <- lower_boundary_from_cutoff(
      n, 0.15, a[1L], sum(a[-1L]), cutoff
    )
    b2 <- lower_boundary_from_cutoff(
      n, 0.30, a[1L] + a[2L], a[3L] + a[4L], cutoff
    )
    no_stop <- b1 < 0L | b2 < 0L
    b1[no_stop] <- -1L
    b2[no_stop] <- -1L
    return(list(boundary1 = b1, boundary2 = b2))
  }

  if (setting$endpoint_id == "coprimary") {
    b1 <- lower_boundary_from_cutoff(
      n, 0.10, a[1L] + a[2L], a[3L] + a[4L], cutoff
    )
    b2 <- lower_boundary_from_cutoff(
      n, 0.20, a[1L] + a[3L], a[2L] + a[4L], cutoff
    )
    no_stop <- b1 < 0L | b2 < 0L
    b1[no_stop] <- -1L
    b2[no_stop] <- -1L
    return(list(boundary1 = b1, boundary2 = b2))
  }

  if (setting$endpoint_id == "efficacy_toxicity") {
    b1 <- lower_boundary_from_cutoff(
      n, 0.45, a[1L] + a[2L], a[3L] + a[4L], cutoff
    )
    b2 <- upper_boundary_from_cutoff(
      n, 0.30, a[1L] + a[3L], a[2L] + a[4L], cutoff
    )
    return(list(boundary1 = b1, boundary2 = b2))
  }

  stop("Unknown endpoint: ", setting$endpoint_id)
}

stop_rule_scalar <- function(counts, setting, boundary1, boundary2, n) {
  events <- endpoint_event_counts(counts, setting$endpoint_id)
  if (setting$endpoint_id %in% c("nested", "coprimary")) {
    if (boundary1 < 0L || boundary2 < 0L) {
      return(rep(FALSE, nrow(counts)))
    }
    return(events$event1 <= boundary1 & events$event2 <= boundary2)
  }

  response_stop <- if (boundary1 < 0L) {
    rep(FALSE, nrow(counts))
  } else {
    events$event1 <= boundary1
  }
  toxicity_stop <- if (boundary2 > n) {
    rep(FALSE, nrow(counts))
  } else {
    events$event2 >= boundary2
  }
  response_stop | toxicity_stop
}

stop_rule_matrix <- function(counts, setting, boundary1, boundary2, n) {
  events <- endpoint_event_counts(counts, setting$endpoint_id)
  if (setting$endpoint_id %in% c("nested", "coprimary")) {
    stop1 <- outer(events$event1, boundary1, FUN = "<=")
    stop2 <- outer(events$event2, boundary2, FUN = "<=")
    valid <- boundary1 >= 0L & boundary2 >= 0L
    out <- stop1 & stop2
    if (any(!valid)) {
      out[, !valid] <- FALSE
    }
    return(out)
  }

  stop1 <- outer(events$event1, boundary1, FUN = "<=")
  stop2 <- outer(events$event2, boundary2, FUN = ">=")
  valid1 <- boundary1 >= 0L
  valid2 <- boundary2 <= n
  if (any(!valid1)) {
    stop1[, !valid1] <- FALSE
  }
  if (any(!valid2)) {
    stop2[, !valid2] <- FALSE
  }
  stop1 | stop2
}

make_master_bop2_candidates <- function(setting) {
  grid <- expand.grid(
    lambda = bop2_lambda_grid,
    gamma = bop2_gamma_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  sample_sizes <- c(n_positive_min:30L, N_positive)

  for (n in sample_sizes) {
    b <- calculate_boundary_vectors(
      setting, n, grid$lambda, grid$gamma
    )
    grid[[paste0("b1_", n)]] <- b$boundary1
    grid[[paste0("b2_", n)]] <- b$boundary2
  }
  grid$candidate_id <- seq_len(nrow(grid))
  grid
}

candidates_for_schedule <- function(master, schedule) {
  signature_columns <- as.vector(rbind(
    paste0("b1_", schedule), paste0("b2_", schedule)
  ))
  signature <- do.call(
    paste,
    c(master[signature_columns], sep = "|")
  )
  out <- master[!duplicated(signature), , drop = FALSE]
  rownames(out) <- NULL
  out
}

## -----------------------------------------------------------------------------
## 4. Standalone conventional BOP2 calibration
## -----------------------------------------------------------------------------

standalone_claim_probabilities <- function(
    cumulative, candidates, schedule, setting, chunk_size = 250L
) {
  n_rep <- nrow(cumulative[[1L]])
  n_candidate <- nrow(candidates)
  out <- numeric(n_candidate)
  starts <- seq.int(1L, n_candidate, by = chunk_size)

  for (start in starts) {
    end <- min(n_candidate, start + chunk_size - 1L)
    idx <- start:end
    active <- matrix(TRUE, nrow = n_rep, ncol = length(idx))

    for (n in schedule) {
      counts <- counts_from_cumulative(cumulative, n)
      stop_now <- stop_rule_matrix(
        counts = counts,
        setting = setting,
        boundary1 = candidates[[paste0("b1_", n)]][idx],
        boundary2 = candidates[[paste0("b2_", n)]][idx],
        n = n
      )
      active <- active & !stop_now
    }
    out[idx] <- colMeans(active)
  }
  out
}

select_standalone_bop2 <- function(
    candidates, type1, power, schedule, label
) {
  feasible <- which(type1 <= alpha + 1e-12)
  if (length(feasible) == 0L) {
    stop("No standalone BOP2 candidate controls type I error for ", label)
  }

  # Original BOP2 ordering: closest type I error not above alpha, then maximum
  # power. Because candidates are unique boundary rules, the final tie is only
  # a deterministic choice among equivalent rules.
  best_type1 <- max(type1[feasible])
  best <- feasible[abs(type1[feasible] - best_type1) <= 1e-12]
  best_power <- max(power[best])
  best <- best[abs(power[best] - best_power) <= 1e-12]
  selected_index <- best[1L]

  list(
    index = selected_index,
    candidate = candidates[selected_index, , drop = FALSE],
    type1 = type1[selected_index],
    power = power[selected_index],
    type1_upper95 = wilson_upper(type1[selected_index], n_calibration),
    schedule = schedule
  )
}

calibrate_standalone_components <- function(setting, endpoint_number) {
  message("  Generating conventional BOP2 candidates for ", setting$endpoint_label)
  master <- make_master_bop2_candidates(setting)

  set.seed(base_seed + endpoint_number * 1000000L + 1100L)
  Y0 <- rmat_categorical(n_calibration, N_positive, setting$theta0)
  cumulative0 <- make_cumulative_counts(Y0)
  rm(Y0)
  invisible(gc())

  set.seed(base_seed + endpoint_number * 1000000L + 1200L)
  Y1 <- rmat_categorical(n_calibration, N_positive, setting$theta1)
  cumulative1 <- make_cumulative_counts(Y1)
  rm(Y1)
  invisible(gc())

  calibrate_schedule <- function(schedule, label) {
    candidates <- candidates_for_schedule(master, schedule)
    type1 <- standalone_claim_probabilities(
      cumulative0, candidates, schedule, setting
    )
    power <- standalone_claim_probabilities(
      cumulative1, candidates, schedule, setting
    )
    selected <- select_standalone_bop2(
      candidates, type1, power, schedule, label
    )
    selected$n_unique_candidates <- nrow(candidates)
    selected
  }

  allcomer <- calibrate_schedule(
    c(interim_n, N_A),
    paste0(setting$endpoint_label, " all-comer BOP2")
  )

  entry_values <- n_positive_min:30L
  positive_list <- parallel_lapply(
    entry_values,
    function(entry_n) {
      schedule <- sort(unique(c(
        entry_n,
        post_enrichment_n[post_enrichment_n > entry_n],
        N_positive
      )))
      selected <- calibrate_schedule(
        schedule,
        paste0(setting$endpoint_label, " positive BOP2 entry ", entry_n)
      )
      candidate <- selected$candidate
      data.frame(
        Endpoint = setting$endpoint_label,
        Endpoint_ID = setting$endpoint_id,
        entry_n = entry_n,
        analysis_schedule = paste(schedule, collapse = ","),
        candidate_id = candidate$candidate_id,
        lambda = candidate$lambda,
        gamma = candidate$gamma,
        standalone_type1 = selected$type1,
        standalone_type1_upper95 = selected$type1_upper95,
        standalone_power = selected$power,
        n_unique_candidates = selected$n_unique_candidates,
        boundary1_at_entry = candidate[[paste0("b1_", entry_n)]],
        boundary2_at_entry = candidate[[paste0("b2_", entry_n)]],
        boundary1_20 = if (20L %in% schedule) candidate$b1_20 else NA_integer_,
        boundary2_20 = if (20L %in% schedule) candidate$b2_20 else NA_integer_,
        boundary1_30 = if (30L %in% schedule) candidate$b1_30 else NA_integer_,
        boundary2_30 = if (30L %in% schedule) candidate$b2_30 else NA_integer_,
        boundary1_40 = candidate$b1_40,
        boundary2_40 = candidate$b2_40,
        stringsAsFactors = FALSE
      )
    },
    mc_cores = n_cores
  )
  positive_library <- do.call(rbind, positive_list)
  rownames(positive_library) <- NULL
  positive_library <- positive_library[order(positive_library$entry_n), ]

  all_candidate <- allcomer$candidate
  allcomer_table <- data.frame(
    Endpoint = setting$endpoint_label,
    Endpoint_ID = setting$endpoint_id,
    analysis_schedule = paste(c(interim_n, N_A), collapse = ","),
    candidate_id = all_candidate$candidate_id,
    lambda = all_candidate$lambda,
    gamma = all_candidate$gamma,
    standalone_type1 = allcomer$type1,
    standalone_type1_upper95 = allcomer$type1_upper95,
    standalone_power = allcomer$power,
    n_unique_candidates = allcomer$n_unique_candidates,
    boundary1_20 = all_candidate$b1_20,
    boundary2_20 = all_candidate$b2_20,
    boundary1_30 = all_candidate$b1_30,
    boundary2_30 = all_candidate$b2_30,
    boundary1_40 = all_candidate$b1_40,
    boundary2_40 = all_candidate$b2_40,
    stringsAsFactors = FALSE
  )

  rm(master, cumulative0, cumulative1)
  invisible(gc())

  list(
    allcomer = allcomer_table,
    positive_library = positive_library
  )
}

## -----------------------------------------------------------------------------
## 5. Potential-data simulation for the adaptive procedure
## -----------------------------------------------------------------------------

simulate_potential_data <- function(n_rep, pi, theta0, endpoint_id, seed) {
  set.seed(seed)
  Z <- matrix(
    rbinom(n_rep * N_A, size = 1L, prob = pi),
    nrow = n_rep, ncol = N_A
  )
  Y <- rmat_categorical(n_rep, N_A, theta0)
  Y_extra_positive <- rmat_categorical(n_rep, N_positive, theta0)

  list(
    endpoint_id = endpoint_id,
    all_counts = list(
      n20 = counts_at_look(Y, 20L),
      n30 = counts_at_look(Y, 30L),
      n40 = counts_at_look(Y, 40L)
    ),
    positive_counts = list(
      n20 = counts_at_look(Y, 20L, Z),
      n30 = counts_at_look(Y, 30L, Z)
    ),
    n_positive = list(
      n20 = rowSums(Z[, 1:20, drop = FALSE]),
      n30 = rowSums(Z[, 1:30, drop = FALSE])
    ),
    extra_cumulative = make_cumulative_counts(Y_extra_positive)
  )
}

positive_counts_at_target <- function(
    entry_counts, entry_n, target_n, extra_cumulative
) {
  out <- entry_counts
  d <- target_n - entry_n
  use_extra <- which(d > 0L)
  if (length(use_extra) > 0L) {
    for (k in 1:4) {
      out[use_extra, k] <- out[use_extra, k] +
        extra_cumulative[[k]][cbind(use_extra, d[use_extra])]
    }
  }
  out
}

## -----------------------------------------------------------------------------
## 6. Proposed-design calibration over the complete adaptive procedure
## -----------------------------------------------------------------------------

build_proposed_allcomer_matrices <- function(data, lambdas) {
  score20 <- posterior_stop_score(data$all_counts$n20, data$endpoint_id)
  score30 <- posterior_stop_score(data$all_counts$n30, data$endpoint_id)
  score40 <- posterior_stop_score(data$all_counts$n40, data$endpoint_id)

  C20 <- 1 - lambdas * (20 / N_A)^proposed_gamma_A
  C30 <- 1 - lambdas * (30 / N_A)^proposed_gamma_A
  C40 <- 1 - lambdas

  stop20 <- outer(score20, C20, FUN = ">")
  stop30_raw <- outer(score30, C30, FUN = ">")
  stop30 <- (!stop20) & stop30_raw
  reject_all <- (!stop20) & (!stop30_raw) &
    outer(score40, C40, FUN = "<=")

  list(cross_look1 = stop20, cross_look2 = stop30, reject_all = reject_all)
}

build_proposed_positive_success <- function(data, entry_look, lambdas) {
  if (entry_look == 1L) {
    entry_counts <- data$positive_counts$n20
    entry_n <- data$n_positive$n20
  } else {
    entry_counts <- data$positive_counts$n30
    entry_n <- data$n_positive$n30
  }

  score_entry <- posterior_stop_score(entry_counts, data$endpoint_id)
  C_entry <- 1 - outer(
    (entry_n / N_positive)^proposed_gamma_positive,
    lambdas
  )

  # No bridging: entry_n < n_positive_min can never proceed.
  success <- outer(
    entry_n >= n_positive_min,
    rep(TRUE, length(lambdas))
  ) & (outer(score_entry, rep(1, length(lambdas))) <= C_entry)

  for (look in post_enrichment_n) {
    required <- entry_n < look
    if (!any(required)) {
      next
    }
    counts_look <- positive_counts_at_target(
      entry_counts, entry_n, look, data$extra_cumulative
    )
    score_look <- posterior_stop_score(counts_look, data$endpoint_id)
    C_look <- 1 - lambdas * (look / N_positive)^proposed_gamma_positive
    pass_look <- outer(score_look, C_look, FUN = "<=")
    required_matrix <- outer(required, rep(TRUE, length(lambdas)))
    success <- success & ((!required_matrix) | pass_look)
  }

  counts_final <- positive_counts_at_target(
    entry_counts, entry_n, N_positive, data$extra_cumulative
  )
  score_final <- posterior_stop_score(counts_final, data$endpoint_id)
  C_final <- 1 - lambdas
  success & outer(score_final, C_final, FUN = "<=")
}

candidate_probabilities_for_prevalence <- function(
    setting, pi, n_rep, seed
) {
  data <- simulate_potential_data(
    n_rep, pi, setting$theta0, setting$endpoint_id, seed
  )
  all_matrices <- build_proposed_allcomer_matrices(
    data, proposed_lambda_grid
  )
  success1 <- build_proposed_positive_success(
    data, 1L, proposed_lambda_grid
  )
  success2 <- build_proposed_positive_success(
    data, 2L, proposed_lambda_grid
  )

  storage.mode(all_matrices$cross_look1) <- "double"
  storage.mode(all_matrices$cross_look2) <- "double"
  storage.mode(all_matrices$reject_all) <- "double"
  storage.mode(success1) <- "double"
  storage.mode(success2) <- "double"

  prn_all <- colMeans(all_matrices$reject_all)
  prn_positive <- (
    crossprod(all_matrices$cross_look1, success1) +
      crossprod(all_matrices$cross_look2, success2)
  ) / n_rep
  prn_any <- sweep(prn_positive, 1L, prn_all, FUN = "+")

  rm(data, all_matrices, success1, success2)
  invisible(gc())

  list(prn_all = prn_all, prn_positive = prn_positive, prn_any = prn_any)
}

calibrate_proposed <- function(setting, endpoint_number) {
  n_lambda <- length(proposed_lambda_grid)
  max_est <- matrix(-Inf, n_lambda, n_lambda)
  min_est <- matrix(Inf, n_lambda, n_lambda)
  sum_est <- matrix(0, n_lambda, n_lambda)
  max_ucb <- matrix(-Inf, n_lambda, n_lambda)
  max_all <- rep(-Inf, n_lambda)
  max_positive <- matrix(-Inf, n_lambda, n_lambda)

  for (p_idx in seq_along(prevalence_grid)) {
    pi <- prevalence_grid[p_idx]
    message("  Proposed calibration: ", setting$endpoint_label, ", pi = ", pi)
    probs <- candidate_probabilities_for_prevalence(
      setting = setting,
      pi = pi,
      n_rep = n_calibration,
      seed = base_seed + endpoint_number * 1000000L + p_idx * 10000L + 1L
    )
    max_est <- pmax(max_est, probs$prn_any)
    min_est <- pmin(min_est, probs$prn_any)
    sum_est <- sum_est + probs$prn_any
    max_ucb <- pmax(max_ucb, wilson_upper(probs$prn_any, n_calibration))
    max_all <- pmax(max_all, probs$prn_all)
    max_positive <- pmax(max_positive, probs$prn_positive)
    rm(probs)
    invisible(gc())
  }

  feasible <- if (proposed_use_wilson_ucb) {
    max_ucb <= alpha + 1e-12
  } else {
    max_est <= alpha + 1e-12
  }
  if (!any(feasible)) {
    stop("No globally feasible proposed candidate for ", setting$endpoint_label)
  }

  preferred <- feasible & max_est >= alpha_lower - 1e-12
  selection_set <- if (any(preferred)) preferred else feasible
  best <- which(selection_set)
  target <- max(max_est[best])
  best <- best[abs(max_est[best] - target) <= 1e-12]
  target <- max(min_est[best])
  best <- best[abs(min_est[best] - target) <= 1e-12]
  target <- max(sum_est[best] / length(prevalence_grid))
  best <- best[
    abs(sum_est[best] / length(prevalence_grid) - target) <= 1e-12
  ]
  ij <- arrayInd(best[1L], dim(max_est))

  list(
    A_index = ij[1L],
    P_index = ij[2L],
    lambda_A = proposed_lambda_grid[ij[1L]],
    lambda_positive = proposed_lambda_grid[ij[2L]],
    gamma_A = proposed_gamma_A,
    gamma_positive = proposed_gamma_positive,
    calibration_max_PRN_any = max_est[ij[1L], ij[2L]],
    calibration_max_PRN_any_upper95 = max_ucb[ij[1L], ij[2L]],
    calibration_min_PRN_any = min_est[ij[1L], ij[2L]],
    calibration_mean_PRN_any =
      sum_est[ij[1L], ij[2L]] / length(prevalence_grid),
    calibration_max_PRN_all = max_all[ij[1L]],
    calibration_max_PRN_positive = max_positive[ij[1L], ij[2L]],
    selected_from_alpha_band = any(preferred)
  )
}

## -----------------------------------------------------------------------------
## 7. Evaluation of selected proposed and comparator procedures
## -----------------------------------------------------------------------------

evaluate_proposed <- function(data, lambda_A, lambda_positive) {
  all_matrices <- build_proposed_allcomer_matrices(data, lambda_A)
  success1 <- build_proposed_positive_success(data, 1L, lambda_positive)[, 1L]
  success2 <- build_proposed_positive_success(data, 2L, lambda_positive)[, 1L]

  reject_all <- all_matrices$reject_all[, 1L]
  reject_positive <-
    (all_matrices$cross_look1[, 1L] & success1) |
    (all_matrices$cross_look2[, 1L] & success2)
  reject_any <- reject_all | reject_positive

  c(
    PRN_any = mean(reject_any),
    PRN_all = mean(reject_all),
    PRN_positive = mean(reject_positive)
  )
}

positive_counts_group_target <- function(
    entry_counts, row_index, entry_n, target_n, extra_cumulative
) {
  out <- entry_counts[row_index, , drop = FALSE]
  d <- target_n - entry_n
  if (d > 0L) {
    for (k in 1:4) {
      out[, k] <- out[, k] + extra_cumulative[[k]][row_index, d]
    }
  }
  out
}

evaluate_comparator_positive_path <- function(
    data, setting, crossing, entry_look, positive_library
) {
  n_rep <- length(crossing)
  reject_positive <- rep(FALSE, n_rep)

  if (entry_look == 1L) {
    entry_counts <- data$positive_counts$n20
    entry_n_vector <- data$n_positive$n20
  } else {
    entry_counts <- data$positive_counts$n30
    entry_n_vector <- data$n_positive$n30
  }

  for (entry_n in n_positive_min:30L) {
    idx <- which(crossing & entry_n_vector == entry_n)
    if (length(idx) == 0L) {
      next
    }
    design <- positive_library[
      positive_library$entry_n == entry_n,
      , drop = FALSE
    ]
    if (nrow(design) != 1L) {
      stop("Missing positive BOP2 library row for entry_n = ", entry_n)
    }

    counts_now <- entry_counts[idx, , drop = FALSE]
    stopped <- stop_rule_scalar(
      counts_now, setting,
      design$boundary1_at_entry,
      design$boundary2_at_entry,
      entry_n
    )
    active_local <- !stopped

    future_looks <- sort(unique(c(
      post_enrichment_n[post_enrichment_n > entry_n],
      N_positive
    )))
    for (look in future_looks) {
      active_positions <- which(active_local)
      if (length(active_positions) == 0L) {
        break
      }
      active_rows <- idx[active_positions]
      counts_look <- positive_counts_group_target(
        entry_counts = entry_counts,
        row_index = active_rows,
        entry_n = entry_n,
        target_n = look,
        extra_cumulative = data$extra_cumulative
      )
      b1 <- if (look == 20L) {
        design$boundary1_20
      } else if (look == 30L) {
        design$boundary1_30
      } else {
        design$boundary1_40
      }
      b2 <- if (look == 20L) {
        design$boundary2_20
      } else if (look == 30L) {
        design$boundary2_30
      } else {
        design$boundary2_40
      }
      stopped_look <- stop_rule_scalar(
        counts_look, setting, b1, b2, look
      )
      active_local[active_positions[stopped_look]] <- FALSE
    }
    reject_positive[idx[active_local]] <- TRUE
  }

  reject_positive
}

evaluate_comparator <- function(
    data, setting, allcomer_design, positive_library
) {
  stop20 <- stop_rule_scalar(
    data$all_counts$n20, setting,
    allcomer_design$boundary1_20,
    allcomer_design$boundary2_20, 20L
  )
  stop30_raw <- stop_rule_scalar(
    data$all_counts$n30, setting,
    allcomer_design$boundary1_30,
    allcomer_design$boundary2_30, 30L
  )
  stop30 <- !stop20 & stop30_raw
  stop40 <- stop_rule_scalar(
    data$all_counts$n40, setting,
    allcomer_design$boundary1_40,
    allcomer_design$boundary2_40, 40L
  )
  reject_all <- !stop20 & !stop30_raw & !stop40

  reject_pos1 <- evaluate_comparator_positive_path(
    data, setting, stop20, 1L, positive_library
  )
  reject_pos2 <- evaluate_comparator_positive_path(
    data, setting, stop30, 2L, positive_library
  )
  reject_positive <- reject_pos1 | reject_pos2
  reject_any <- reject_all | reject_positive

  c(
    PRN_any = mean(reject_any),
    PRN_all = mean(reject_all),
    PRN_positive = mean(reject_positive)
  )
}

## -----------------------------------------------------------------------------
## 8. Boundary output helpers
## -----------------------------------------------------------------------------

make_proposed_boundary_rows <- function(setting, proposed) {
  rows <- list()
  idx <- 1L
  for (population in c("All-comer", "Biomarker-positive")) {
    sample_sizes <- if (population == "All-comer") {
      c(interim_n, N_A)
    } else {
      c(n_positive_min:30L, N_positive)
    }
    lambda <- if (population == "All-comer") {
      proposed$lambda_A
    } else {
      proposed$lambda_positive
    }
    gamma <- if (population == "All-comer") {
      proposed$gamma_A
    } else {
      proposed$gamma_positive
    }

    for (n in sample_sizes) {
      b <- calculate_boundary_vectors(setting, n, lambda, gamma)
      boundary1_output <- b$boundary1
      boundary2_output <- b$boundary2
      if (setting$endpoint_id %in% c("nested", "coprimary")) {
        if (boundary1_output < 0L || boundary2_output < 0L) {
          boundary1_output <- NA_integer_
          boundary2_output <- NA_integer_
        }
      } else {
        if (boundary1_output < 0L) boundary1_output <- NA_integer_
        if (boundary2_output > n) boundary2_output <- NA_integer_
      }
      rows[[idx]] <- data.frame(
        Endpoint = setting$endpoint_label,
        Endpoint_ID = setting$endpoint_id,
        Population = population,
        n = n,
        boundary1 = boundary1_output,
        boundary2 = boundary2_output,
        lambda = lambda,
        gamma = gamma,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

make_comparator_boundary_rows <- function(setting, standalone) {
  normalize_output_boundary <- function(boundary1, boundary2, n) {
    if (setting$endpoint_id %in% c("nested", "coprimary")) {
      if (boundary1 < 0L || boundary2 < 0L) {
        return(c(boundary1 = NA_integer_, boundary2 = NA_integer_))
      }
    } else {
      if (boundary1 < 0L) boundary1 <- NA_integer_
      if (boundary2 > n) boundary2 <- NA_integer_
    }
    c(boundary1 = boundary1, boundary2 = boundary2)
  }

  all_d <- standalone$allcomer
  all_n <- c(20L, 30L, 40L)
  all_raw1 <- c(all_d$boundary1_20, all_d$boundary1_30, all_d$boundary1_40)
  all_raw2 <- c(all_d$boundary2_20, all_d$boundary2_30, all_d$boundary2_40)
  all_normalized <- t(vapply(seq_along(all_n), function(i) {
    normalize_output_boundary(all_raw1[i], all_raw2[i], all_n[i])
  }, numeric(2L)))
  all_rows <- data.frame(
    Endpoint = setting$endpoint_label,
    Endpoint_ID = setting$endpoint_id,
    Population = "All-comer",
    entry_n = NA_integer_,
    n = all_n,
    boundary1 = all_normalized[, "boundary1"],
    boundary2 = all_normalized[, "boundary2"],
    lambda = all_d$lambda,
    gamma = all_d$gamma,
    stringsAsFactors = FALSE
  )

  positive_rows <- list()
  out_idx <- 1L
  for (i in seq_len(nrow(standalone$positive_library))) {
    d <- standalone$positive_library[i, , drop = FALSE]
    schedule <- as.integer(strsplit(d$analysis_schedule, ",", fixed = TRUE)[[1L]])
    for (n in schedule) {
      if (n == d$entry_n) {
        b1 <- d$boundary1_at_entry
        b2 <- d$boundary2_at_entry
      } else if (n == 20L) {
        b1 <- d$boundary1_20
        b2 <- d$boundary2_20
      } else if (n == 30L) {
        b1 <- d$boundary1_30
        b2 <- d$boundary2_30
      } else {
        b1 <- d$boundary1_40
        b2 <- d$boundary2_40
      }
      normalized <- normalize_output_boundary(b1, b2, n)
      positive_rows[[out_idx]] <- data.frame(
        Endpoint = setting$endpoint_label,
        Endpoint_ID = setting$endpoint_id,
        Population = "Biomarker-positive",
        entry_n = d$entry_n,
        n = n,
        boundary1 = normalized["boundary1"],
        boundary2 = normalized["boundary2"],
        lambda = d$lambda,
        gamma = d$gamma,
        stringsAsFactors = FALSE
      )
      out_idx <- out_idx + 1L
    }
  }
  rbind(all_rows, do.call(rbind, positive_rows))
}

## -----------------------------------------------------------------------------
## 9. Endpoint-by-endpoint calibration and independent validation
## -----------------------------------------------------------------------------

all_selected_rows <- list()
all_type1_rows <- list()
all_summary_rows <- list()
all_standalone_all <- list()
all_positive_libraries <- list()
all_proposed_boundaries <- list()
all_comparator_boundaries <- list()

for (endpoint_number in seq_along(endpoint_settings)) {
  setting <- endpoint_settings[[endpoint_number]]
  message("\nEndpoint: ", setting$endpoint_label)

  proposed <- calibrate_proposed(setting, endpoint_number)
  standalone <- calibrate_standalone_components(setting, endpoint_number)

  validation_rows <- list()
  out_idx <- 1L
  for (p_idx in seq_along(prevalence_grid)) {
    pi <- prevalence_grid[p_idx]
    message("  Independent validation: pi = ", pi)
    data <- simulate_potential_data(
      n_rep = n_validation,
      pi = pi,
      theta0 = setting$theta0,
      endpoint_id = setting$endpoint_id,
      seed = base_seed + endpoint_number * 1000000L + p_idx * 10000L + 5000L
    )

    proposed_est <- evaluate_proposed(
      data, proposed$lambda_A, proposed$lambda_positive
    )
    comparator_est <- evaluate_comparator(
      data, setting, standalone$allcomer, standalone$positive_library
    )

    for (method_index in 1:2) {
      method <- c(
        "Globally calibrated",
        "Independent conventional BOP2 components"
      )[method_index]
      estimates <- if (method_index == 1L) proposed_est else comparator_est
      ci <- wilson_interval(estimates["PRN_any"], n_validation)
      validation_rows[[out_idx]] <- data.frame(
        Endpoint = setting$endpoint_label,
        Endpoint_ID = setting$endpoint_id,
        Method = method,
        pi = pi,
        Monte_Carlo_replicates = n_validation,
        PRN_any = unname(estimates["PRN_any"]),
        PRN_any_MCSE = sqrt(
          estimates["PRN_any"] * (1 - estimates["PRN_any"]) / n_validation
        ),
        PRN_any_CI_lower = ci["lower"],
        PRN_any_CI_upper = ci["upper"],
        PRN_all = unname(estimates["PRN_all"]),
        PRN_positive = unname(estimates["PRN_positive"]),
        stringsAsFactors = FALSE
      )
      out_idx <- out_idx + 1L
    }
    rm(data)
    invisible(gc())
  }
  endpoint_type1 <- do.call(rbind, validation_rows)
  rownames(endpoint_type1) <- NULL

  endpoint_summary <- do.call(
    rbind,
    lapply(split(endpoint_type1, endpoint_type1$Method), function(d) {
      worst <- d[which.max(d$PRN_any), , drop = FALSE]
      data.frame(
        Endpoint = setting$endpoint_label,
        Endpoint_ID = setting$endpoint_id,
        Method = d$Method[1L],
        Max_PRN_any = max(d$PRN_any),
        Prevalence_at_max = worst$pi,
        Max_PRN_any_CI_upper = worst$PRN_any_CI_upper,
        Max_PRN_all = max(d$PRN_all),
        Max_PRN_positive = max(d$PRN_positive),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(endpoint_summary) <- NULL

  proposed_summary <- endpoint_summary[
    endpoint_summary$Method == "Globally calibrated", , drop = FALSE
  ]
  comparator_summary <- endpoint_summary[
    endpoint_summary$Method ==
      "Independent conventional BOP2 components", , drop = FALSE
  ]

  selected_rows <- rbind(
    data.frame(
      Endpoint = setting$endpoint_label,
      Endpoint_ID = setting$endpoint_id,
      Method = "Globally calibrated",
      Positive_component = "Single globally calibrated boundary function",
      lambda_A = proposed$lambda_A,
      gamma_A = proposed$gamma_A,
      lambda_positive = proposed$lambda_positive,
      gamma_positive = proposed$gamma_positive,
      calibration_max_PRN_any = proposed$calibration_max_PRN_any,
      calibration_max_PRN_any_upper95 =
        proposed$calibration_max_PRN_any_upper95,
      validation_max_PRN_any = proposed_summary$Max_PRN_any,
      validation_max_PRN_all = proposed_summary$Max_PRN_all,
      validation_max_PRN_positive = proposed_summary$Max_PRN_positive,
      standalone_allcomer_type1 = NA_real_,
      standalone_allcomer_power = NA_real_,
      positive_standalone_type1_min = NA_real_,
      positive_standalone_type1_max = NA_real_,
      positive_standalone_power_min = NA_real_,
      positive_standalone_power_max = NA_real_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Endpoint = setting$endpoint_label,
      Endpoint_ID = setting$endpoint_id,
      Method = "Independent conventional BOP2 components",
      Positive_component = "Entry-specific conventional BOP2 library",
      lambda_A = standalone$allcomer$lambda,
      gamma_A = standalone$allcomer$gamma,
      lambda_positive = NA_real_,
      gamma_positive = NA_real_,
      calibration_max_PRN_any = NA_real_,
      calibration_max_PRN_any_upper95 = NA_real_,
      validation_max_PRN_any = comparator_summary$Max_PRN_any,
      validation_max_PRN_all = comparator_summary$Max_PRN_all,
      validation_max_PRN_positive = comparator_summary$Max_PRN_positive,
      standalone_allcomer_type1 = standalone$allcomer$standalone_type1,
      standalone_allcomer_power = standalone$allcomer$standalone_power,
      positive_standalone_type1_min =
        min(standalone$positive_library$standalone_type1),
      positive_standalone_type1_max =
        max(standalone$positive_library$standalone_type1),
      positive_standalone_power_min =
        min(standalone$positive_library$standalone_power),
      positive_standalone_power_max =
        max(standalone$positive_library$standalone_power),
      stringsAsFactors = FALSE
    )
  )

  all_selected_rows[[endpoint_number]] <- selected_rows
  all_type1_rows[[endpoint_number]] <- endpoint_type1
  all_summary_rows[[endpoint_number]] <- endpoint_summary
  all_standalone_all[[endpoint_number]] <- standalone$allcomer
  all_positive_libraries[[endpoint_number]] <- standalone$positive_library
  all_proposed_boundaries[[endpoint_number]] <-
    make_proposed_boundary_rows(setting, proposed)
  all_comparator_boundaries[[endpoint_number]] <-
    make_comparator_boundary_rows(setting, standalone)
}

selected_designs <- do.call(rbind, all_selected_rows)
type1_results <- do.call(rbind, all_type1_rows)
type1_summary <- do.call(rbind, all_summary_rows)
standalone_allcomer <- do.call(rbind, all_standalone_all)
standalone_positive_library <- do.call(rbind, all_positive_libraries)
proposed_boundaries <- do.call(rbind, all_proposed_boundaries)
comparator_boundaries <- do.call(rbind, all_comparator_boundaries)

rownames(selected_designs) <- NULL
rownames(type1_results) <- NULL
rownames(type1_summary) <- NULL
rownames(standalone_allcomer) <- NULL
rownames(standalone_positive_library) <- NULL
rownames(proposed_boundaries) <- NULL
rownames(comparator_boundaries) <- NULL

# Numerical identity check for mutually exclusive claim paths.
if (max(abs(
  type1_results$PRN_any - type1_results$PRN_all -
    type1_results$PRN_positive
)) > 5 / n_validation) {
  stop("PRN-any is inconsistent with PRN-all + PRN-positive.")
}

## -----------------------------------------------------------------------------
## 10. Output files
## -----------------------------------------------------------------------------

settings_table <- data.frame(
  Parameter = c(
    "analysis type", "branching rule", "alpha", "alpha_lower",
    "all-comer looks", "positive scheduled looks", "N_A", "N_positive",
    "n_positive_min", "prevalence grid", "Dirichlet prior",
    "proposed gamma_A", "proposed gamma_positive",
    "proposed lambda grid", "standalone BOP2 lambda grid",
    "standalone BOP2 gamma grid", "calibration replicates",
    "independent validation replicates",
    "proposed feasibility criterion", "standalone BOP2 selection rule",
    "base seed", "parallel workers"
  ),
  Value = c(
    "Monte Carlo simulation", "No bridging", alpha, alpha_lower,
    paste(c(interim_n, N_A), collapse = ", "),
    paste(c(post_enrichment_n, N_positive), collapse = ", "),
    N_A, N_positive, n_positive_min,
    paste(prevalence_grid, collapse = ", "),
    paste(dirichlet_prior, collapse = ", "),
    proposed_gamma_A, proposed_gamma_positive,
    paste0(
      min(proposed_lambda_grid), " to ", max(proposed_lambda_grid), " by ",
      if (length(proposed_lambda_grid) > 1L) {
        proposed_lambda_grid[2L] - proposed_lambda_grid[1L]
      } else NA
    ),
    paste0(
      min(bop2_lambda_grid), " to ", max(bop2_lambda_grid), " by ",
      if (length(bop2_lambda_grid) > 1L) {
        bop2_lambda_grid[2L] - bop2_lambda_grid[1L]
      } else NA
    ),
    paste0(
      min(bop2_gamma_grid), " to ", max(bop2_gamma_grid), " by ",
      if (length(bop2_gamma_grid) > 1L) {
        bop2_gamma_grid[2L] - bop2_gamma_grid[1L]
      } else NA
    ),
    n_calibration, n_validation,
    if (proposed_use_wilson_ucb) {
      "one-sided 95% Wilson upper bound <= alpha"
    } else {
      "Monte Carlo estimate <= alpha"
    },
    "closest estimated type I error <= alpha, then maximum power",
    base_seed, n_cores
  ),
  stringsAsFactors = FALSE
)

scenario_rows <- lapply(endpoint_settings, function(setting) {
  data.frame(
    Endpoint_ID = setting$endpoint_id,
    Endpoint = setting$endpoint_label,
    Categories = paste(setting$categories, collapse = "; "),
    theta0 = paste(setting$theta0, collapse = ", "),
    theta1_for_standalone_BOP2 = paste(setting$theta1, collapse = ", "),
    Monitoring_rule = setting$rule_label,
    stringsAsFactors = FALSE
  )
})
endpoint_scenarios <- do.call(rbind, scenario_rows)

write.csv(
  settings_table,
  file.path(result_dir, "complex_endpoint_type1_settings.csv"),
  row.names = FALSE
)
write.csv(
  endpoint_scenarios,
  file.path(result_dir, "complex_endpoint_scenarios.csv"),
  row.names = FALSE
)
write.csv(
  selected_designs,
  file.path(result_dir, "complex_endpoint_selected_designs.csv"),
  row.names = FALSE
)
write.csv(
  type1_results,
  file.path(result_dir, "complex_endpoint_type1_by_prevalence.csv"),
  row.names = FALSE
)
write.csv(
  type1_summary,
  file.path(result_dir, "complex_endpoint_type1_summary.csv"),
  row.names = FALSE
)
write.csv(
  standalone_allcomer,
  file.path(result_dir, "complex_endpoint_standalone_allcomer_bop2.csv"),
  row.names = FALSE
)
write.csv(
  standalone_positive_library,
  file.path(result_dir, "complex_endpoint_standalone_positive_bop2_library.csv"),
  row.names = FALSE
)
write.csv(
  proposed_boundaries,
  file.path(result_dir, "complex_endpoint_proposed_boundaries.csv"),
  row.names = FALSE
)
write.csv(
  comparator_boundaries,
  file.path(result_dir, "complex_endpoint_comparator_boundaries.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    settings = settings_table,
    scenarios = endpoint_scenarios,
    selected_designs = selected_designs,
    type1_results = type1_results,
    type1_summary = type1_summary,
    standalone_allcomer = standalone_allcomer,
    standalone_positive_library = standalone_positive_library,
    proposed_boundaries = proposed_boundaries,
    comparator_boundaries = comparator_boundaries
  ),
  file.path(result_dir, "complex_endpoint_type1_results.rds")
)

xlsx_path <- file.path(result_dir, "complex_endpoint_type1_results.xlsx")
xlsx_created <- write_results_workbook(
  sheets = list(
    Settings = settings_table,
    Scenarios = endpoint_scenarios,
    Selected_designs = selected_designs,
    Type1_by_prevalence = type1_results,
    Type1_summary = type1_summary,
    Standalone_allcomer_BOP2 = standalone_allcomer,
    Positive_BOP2_library = standalone_positive_library,
    Proposed_boundaries = proposed_boundaries,
    Comparator_boundaries = comparator_boundaries
  ),
  path = xlsx_path
)

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo.txt")
)

message("\nSelected designs:")
print(selected_designs)
message("\nIndependent validation summary:")
print(type1_summary)
message("\nCompleted. Results were written to: ", result_dir)
if (xlsx_created) {
  message("Consolidated workbook: ", xlsx_path)
}
