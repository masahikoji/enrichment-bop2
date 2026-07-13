# Complex-endpoint type I error analysis for the supplementary material.
# Uses simulation-based calibration and relative repository paths.

options(stringsAsFactors = FALSE, scipen = 999)

# User settings

alpha <- 0.10

interim_n <- c(20L, 30L)
N_A <- 40L
N_pos <- 40L
n_pos_min <- 10L
post_enrichment_n <- c(20L, 30L)

prevalence_grid <- c(0.4, 0.5, 0.6, 0.7, 0.8)

gamma_A <- 1
gamma_pos <- 1
lambda_grid <- seq(0.005, 1.000, by = 0.005)

dirichlet_prior <- rep(0.25, 4L)

n_calibration <- as.integer(Sys.getenv(
  "BOP2_COMPLEX_N_CALIB", unset = "50000"
))
n_validation <- as.integer(Sys.getenv(
  "BOP2_COMPLEX_N_VALID", unset = "100000"
))
base_seed <- 20260712L

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
default_cores <- max(1L, min(3L, detected_cores - 1L))
requested_cores <- suppressWarnings(as.integer(Sys.getenv(
  "BOP2_N_CORES", unset = as.character(default_cores)
)))
if (is.na(requested_cores) || requested_cores < 1L) {
  requested_cores <- default_cores
}
n_cores <- max(1L, min(3L, requested_cores))

quick_test <- identical(Sys.getenv("BOP2_COMPLEX_QUICK"), "1")
if (quick_test) {
  message(
    "BOP2_COMPLEX_QUICK=1: running a reduced code check, ",
    "not the manuscript analysis."
  )
  prevalence_grid <- c(0.4, 0.8)
  lambda_grid <- seq(0.05, 1.00, by = 0.05)
  n_calibration <- 2000L
  n_validation <- 5000L
  n_cores <- 1L
  result_dir <- file.path(result_dir, "quick_test")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(
  length(interim_n) == 2L,
  all(diff(interim_n) > 0L),
  max(interim_n) < N_A,
  N_A == N_pos,
  n_pos_min <= min(interim_n),
  all(post_enrichment_n < N_pos),
  all(dirichlet_prior > 0),
  abs(sum(dirichlet_prior) - 1) < 1e-12,
  n_calibration > 0L,
  n_validation > 0L
)

message(
  "Parallel configuration: physical cores = ", detected_cores,
  "; endpoint-level workers = ", n_cores,
  "; backend = ",
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    "parallel::mclapply"
  } else {
    "sequential lapply"
  }
)

# Complex-endpoint global-null settings

endpoint_settings <- list(
  nested = list(
    endpoint_id = "nested",
    endpoint_label = "Nested efficacy",
    categories = c("CR", "PR", "SD", "PD"),
    theta0 = c(0.15, 0.15, 0.30, 0.40),
    rule_label = paste0(
      "Stop if both CR and CR/PR posterior futility conditions are met"
    )
  ),
  coprimary = list(
    endpoint_id = "coprimary",
    endpoint_label = "Co-primary efficacy",
    categories = c(
      "OR_and_EFS6", "OR_and_no_EFS6",
      "no_OR_and_EFS6", "no_OR_and_no_EFS6"
    ),
    theta0 = c(0.05, 0.05, 0.15, 0.75),
    rule_label = paste0(
      "Stop if both ORR and EFS6 posterior futility conditions are met"
    )
  ),
  efficacy_toxicity = list(
    endpoint_id = "efficacy_toxicity",
    endpoint_label = "Efficacy and toxicity",
    categories = c(
      "toxicity_and_OR", "no_toxicity_and_OR",
      "toxicity_and_no_OR", "no_toxicity_and_no_OR"
    ),
    theta0 = c(0.15, 0.30, 0.15, 0.40),
    rule_label = paste0(
      "Stop if efficacy is futile or toxicity is excessive"
    )
  )
)

for (setting in endpoint_settings) {
  if (length(setting$theta0) != 4L ||
      any(setting$theta0 < 0) ||
      abs(sum(setting$theta0) - 1) > 1e-12) {
    stop("Invalid null category-probability vector for ", setting$endpoint_id)
  }
}

# Utility functions

parallel_lapply <- function(X, FUN, ..., mc_cores = n_cores) {
  if (mc_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      X, FUN, ...,
      mc.cores = mc_cores,
      mc.preschedule = TRUE
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
      openxlsx::writeData(
        wb, sheet = sheet_name, x = sheets[[sheet_name]]
      )
      openxlsx::freezePane(
        wb, sheet = sheet_name, firstRow = TRUE
      )
      openxlsx::setColWidths(
        wb,
        sheet = sheet_name,
        cols = seq_len(ncol(sheets[[sheet_name]])),
        widths = "auto"
      )
    }
    openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
    return(TRUE)
  }

  warning(
    "Neither 'writexl' nor 'openxlsx' is installed. ",
    "CSV files were created, but the XLSX workbook was not."
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
  cbind(
    lower = pmax(0, center - half),
    upper = pmin(1, center + half)
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
      out[, k] <- rowSums(
        (Y[, idx, drop = FALSE] == k) & z_sub
      )
    }
  }
  colnames(out) <- paste0("x", 1:4)
  out
}

posterior_stop_score <- function(counts, endpoint_id) {
  counts <- as.matrix(counts)
  if (ncol(counts) != 4L) {
    stop("counts must have four columns.")
  }

  n <- rowSums(counts)
  a <- dirichlet_prior

  if (endpoint_id == "nested") {
    shape_cr_1 <- a[1L] + counts[, 1L]
    shape_cr_2 <- sum(a[-1L]) + n - counts[, 1L]
    p_cr_futile <- pbeta(
      0.15, shape1 = shape_cr_1, shape2 = shape_cr_2
    )

    crpr <- counts[, 1L] + counts[, 2L]
    shape_crpr_1 <- a[1L] + a[2L] + crpr
    shape_crpr_2 <- a[3L] + a[4L] + n - crpr
    p_crpr_futile <- pbeta(
      0.30, shape1 = shape_crpr_1, shape2 = shape_crpr_2
    )

    return(pmin(p_cr_futile, p_crpr_futile))
  }

  if (endpoint_id == "coprimary") {
    or_count <- counts[, 1L] + counts[, 2L]
    efs_count <- counts[, 1L] + counts[, 3L]

    p_or_futile <- pbeta(
      0.10,
      shape1 = a[1L] + a[2L] + or_count,
      shape2 = a[3L] + a[4L] + n - or_count
    )
    p_efs_futile <- pbeta(
      0.20,
      shape1 = a[1L] + a[3L] + efs_count,
      shape2 = a[2L] + a[4L] + n - efs_count
    )

    return(pmin(p_or_futile, p_efs_futile))
  }

  if (endpoint_id == "efficacy_toxicity") {
    response_count <- counts[, 1L] + counts[, 2L]
    toxicity_count <- counts[, 1L] + counts[, 3L]

    p_eff_futile <- pbeta(
      0.45,
      shape1 = a[1L] + a[2L] + response_count,
      shape2 = a[3L] + a[4L] + n - response_count
    )
    p_tox_excessive <- pbeta(
      0.30,
      shape1 = a[1L] + a[3L] + toxicity_count,
      shape2 = a[2L] + a[4L] + n - toxicity_count,
      lower.tail = FALSE
    )

    return(pmax(p_eff_futile, p_tox_excessive))
  }

  stop("Unknown endpoint_id: ", endpoint_id)
}

simulate_potential_data <- function(
    n_rep, pi, theta0, endpoint_id, seed
) {
  set.seed(seed)

  Z <- matrix(
    rbinom(n_rep * N_A, size = 1L, prob = pi),
    nrow = n_rep,
    ncol = N_A
  )
  Y <- rmat_categorical(n_rep, N_A, theta0)
  Y_extra_positive <- rmat_categorical(n_rep, N_pos, theta0)

  all_counts <- list(
    n20 = counts_at_look(Y, interim_n[1L]),
    n30 = counts_at_look(Y, interim_n[2L]),
    n40 = counts_at_look(Y, N_A)
  )
  positive_counts <- list(
    n20 = counts_at_look(Y, interim_n[1L], Z),
    n30 = counts_at_look(Y, interim_n[2L], Z)
  )
  n_positive <- list(
    n20 = rowSums(Z[, seq_len(interim_n[1L]), drop = FALSE]),
    n30 = rowSums(Z[, seq_len(interim_n[2L]), drop = FALSE])
  )

  extra_cumulative <- lapply(
    1:4,
    function(k) row_cumulative_indicator(Y_extra_positive, k)
  )

  list(
    endpoint_id = endpoint_id,
    all_counts = all_counts,
    positive_counts = positive_counts,
    n_positive = n_positive,
    extra_cumulative = extra_cumulative
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
        extra_cumulative[[k]][
          cbind(use_extra, d[use_extra])
        ]
    }
  }
  out
}

build_allcomer_candidate_matrices <- function(data, lambdas) {
  score20 <- posterior_stop_score(
    data$all_counts$n20, data$endpoint_id
  )
  score30 <- posterior_stop_score(
    data$all_counts$n30, data$endpoint_id
  )
  score40 <- posterior_stop_score(
    data$all_counts$n40, data$endpoint_id
  )

  C20 <- 1 - lambdas * (interim_n[1L] / N_A)^gamma_A
  C30 <- 1 - lambdas * (interim_n[2L] / N_A)^gamma_A
  C40 <- 1 - lambdas

  stop20 <- outer(score20, C20, FUN = ">")
  stop30_raw <- outer(score30, C30, FUN = ">")
  stop30 <- (!stop20) & stop30_raw
  reject_all <- (!stop20) & (!stop30_raw) &
    outer(score40, C40, FUN = "<=")

  list(
    cross_look1 = stop20,
    cross_look2 = stop30,
    reject_all = reject_all
  )
}

build_positive_success_matrix <- function(
    data, entry_look, lambdas
) {
  if (entry_look == 1L) {
    entry_counts <- data$positive_counts$n20
    entry_n <- data$n_positive$n20
  } else if (entry_look == 2L) {
    entry_counts <- data$positive_counts$n30
    entry_n <- data$n_positive$n30
  } else {
    stop("entry_look must be 1 or 2.")
  }

  score_entry <- posterior_stop_score(
    entry_counts, data$endpoint_id
  )
  C_entry <- 1 - outer(
    (entry_n / N_pos)^gamma_pos,
    lambdas
  )

  success <- outer(entry_n >= n_pos_min, rep(TRUE, length(lambdas))) &
    (outer(score_entry, rep(1, length(lambdas))) <= C_entry)

  for (look in post_enrichment_n) {
    required <- entry_n < look
    if (!any(required)) {
      next
    }

    counts_look <- positive_counts_at_target(
      entry_counts = entry_counts,
      entry_n = entry_n,
      target_n = look,
      extra_cumulative = data$extra_cumulative
    )
    score_look <- posterior_stop_score(
      counts_look, data$endpoint_id
    )
    C_look <- 1 - lambdas * (look / N_pos)^gamma_pos
    pass_look <- outer(score_look, C_look, FUN = "<=")

    required_matrix <- outer(
      required, rep(TRUE, length(lambdas))
    )
    success <- success & ((!required_matrix) | pass_look)
  }

  counts_final <- positive_counts_at_target(
    entry_counts = entry_counts,
    entry_n = entry_n,
    target_n = N_pos,
    extra_cumulative = data$extra_cumulative
  )
  score_final <- posterior_stop_score(
    counts_final, data$endpoint_id
  )
  C_final <- 1 - lambdas
  success <- success &
    outer(score_final, C_final, FUN = "<=")

  success
}

candidate_probabilities_for_prevalence <- function(
    setting, pi, n_rep, seed
) {
  data <- simulate_potential_data(
    n_rep = n_rep,
    pi = pi,
    theta0 = setting$theta0,
    endpoint_id = setting$endpoint_id,
    seed = seed
  )

  all_matrices <- build_allcomer_candidate_matrices(
    data, lambda_grid
  )
  success1 <- build_positive_success_matrix(
    data, entry_look = 1L, lambdas = lambda_grid
  )
  success2 <- build_positive_success_matrix(
    data, entry_look = 2L, lambdas = lambda_grid
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

  prn_any <- sweep(
    prn_positive,
    MARGIN = 1L,
    STATS = prn_all,
    FUN = "+"
  )

  rm(data, all_matrices, success1, success2)
  invisible(gc())

  list(
    pi = pi,
    prn_all = prn_all,
    prn_positive = prn_positive,
    prn_any = prn_any
  )
}

select_proposed_pair <- function(calibration) {
  feasible <- calibration$max_global_ucb <= alpha + 1e-12
  if (!any(feasible)) {
    stop(
      "No globally feasible candidate. Increase lambda_grid or ",
      "increase the calibration replicate count."
    )
  }

  # No alternative configuration is used in this supplementary type I error
  # demonstration. Select the feasible candidate with the largest estimated
  # maximum global type I error, then use the minimum and mean values as ties.
  score1 <- calibration$max_global_estimate
  best <- which(feasible)
  tol <- 1e-12

  target <- max(score1[best])
  best <- best[abs(score1[best] - target) <= tol]

  target <- max(calibration$min_global_estimate[best])
  best <- best[
    abs(calibration$min_global_estimate[best] - target) <= tol
  ]

  target <- max(calibration$mean_global_estimate[best])
  best <- best[
    abs(calibration$mean_global_estimate[best] - target) <= tol
  ]

  ij <- arrayInd(best[1L], dim(calibration$max_global_estimate))
  list(
    A_index = ij[1L],
    P_index = ij[2L],
    max_global_estimate =
      calibration$max_global_estimate[ij[1L], ij[2L]],
    max_global_ucb =
      calibration$max_global_ucb[ij[1L], ij[2L]]
  )
}

select_componentwise_pair <- function(calibration) {
  n_lambda <- length(lambda_grid)
  all_feasible <- calibration$max_all_ucb <= alpha + 1e-12
  positive_feasible <-
    calibration$max_positive_ucb <= alpha + 1e-12

  feasible <- matrix(
    all_feasible,
    nrow = n_lambda,
    ncol = n_lambda
  ) & positive_feasible

  if (!any(feasible)) {
    stop(
      "No componentwise-feasible candidate. Increase lambda_grid or ",
      "increase the calibration replicate count."
    )
  }

  all_estimate_matrix <- matrix(
    calibration$max_all_estimate,
    nrow = n_lambda,
    ncol = n_lambda
  )

  # Select the feasible pair whose two marginal type I error estimates are
  # closest to (alpha, alpha) in squared Euclidean distance. No constraint is
  # imposed on their union for the componentwise comparator.
  distance <- (alpha - all_estimate_matrix)^2 +
    (alpha - calibration$max_positive_estimate)^2
  best <- which(feasible)
  tol <- 1e-12

  target <- min(distance[best])
  best <- best[abs(distance[best] - target) <= tol]

  target <- max(calibration$max_global_estimate[best])
  best <- best[
    abs(calibration$max_global_estimate[best] - target) <= tol
  ]

  ij <- arrayInd(best[1L], dim(calibration$max_global_estimate))
  list(
    A_index = ij[1L],
    P_index = ij[2L],
    max_global_estimate =
      calibration$max_global_estimate[ij[1L], ij[2L]],
    max_global_ucb =
      calibration$max_global_ucb[ij[1L], ij[2L]]
  )
}

calibrate_endpoint <- function(setting, endpoint_number) {
  message("Calibrating endpoint: ", setting$endpoint_label)

  n_lambda <- length(lambda_grid)
  max_global_estimate <- matrix(
    -Inf, nrow = n_lambda, ncol = n_lambda
  )
  min_global_estimate <- matrix(
    Inf, nrow = n_lambda, ncol = n_lambda
  )
  sum_global_estimate <- matrix(
    0, nrow = n_lambda, ncol = n_lambda
  )
  max_global_ucb <- matrix(
    -Inf, nrow = n_lambda, ncol = n_lambda
  )

  max_positive_estimate <- matrix(
    -Inf, nrow = n_lambda, ncol = n_lambda
  )
  max_positive_ucb <- matrix(
    -Inf, nrow = n_lambda, ncol = n_lambda
  )
  max_all_estimate <- rep(-Inf, n_lambda)
  max_all_ucb <- rep(-Inf, n_lambda)

  for (p_idx in seq_along(prevalence_grid)) {
    pi <- prevalence_grid[p_idx]
    message(
      "  ", setting$endpoint_label,
      ": calibration prevalence ", pi
    )

    probabilities <- candidate_probabilities_for_prevalence(
      setting = setting,
      pi = pi,
      n_rep = n_calibration,
      seed = base_seed +
        endpoint_number * 1000000L +
        p_idx * 10000L + 1L
    )

    global_ucb <- wilson_upper(
      probabilities$prn_any, n_calibration
    )
    positive_ucb <- wilson_upper(
      probabilities$prn_positive, n_calibration
    )
    all_ucb <- wilson_upper(
      probabilities$prn_all, n_calibration
    )

    max_global_estimate <- pmax(
      max_global_estimate, probabilities$prn_any
    )
    min_global_estimate <- pmin(
      min_global_estimate, probabilities$prn_any
    )
    sum_global_estimate <- sum_global_estimate +
      probabilities$prn_any
    max_global_ucb <- pmax(max_global_ucb, global_ucb)

    max_positive_estimate <- pmax(
      max_positive_estimate,
      probabilities$prn_positive
    )
    max_positive_ucb <- pmax(
      max_positive_ucb, positive_ucb
    )
    max_all_estimate <- pmax(
      max_all_estimate, probabilities$prn_all
    )
    max_all_ucb <- pmax(max_all_ucb, all_ucb)

    rm(probabilities, global_ucb, positive_ucb, all_ucb)
    invisible(gc())
  }

  calibration <- list(
    max_global_estimate = max_global_estimate,
    min_global_estimate = min_global_estimate,
    mean_global_estimate =
      sum_global_estimate / length(prevalence_grid),
    max_global_ucb = max_global_ucb,
    max_positive_estimate = max_positive_estimate,
    max_positive_ucb = max_positive_ucb,
    max_all_estimate = max_all_estimate,
    max_all_ucb = max_all_ucb
  )

  proposed <- select_proposed_pair(calibration)
  componentwise <- select_componentwise_pair(calibration)

  selected <- data.frame(
    Endpoint = setting$endpoint_label,
    Method = c(
      "Globally calibrated",
      "Componentwise calibrated"
    ),
    lambda_A = c(
      lambda_grid[proposed$A_index],
      lambda_grid[componentwise$A_index]
    ),
    gamma_A = gamma_A,
    lambda_positive = c(
      lambda_grid[proposed$P_index],
      lambda_grid[componentwise$P_index]
    ),
    gamma_positive = gamma_pos,
    calibration_max_PRN_any = c(
      proposed$max_global_estimate,
      componentwise$max_global_estimate
    ),
    calibration_max_PRN_any_upper95 = c(
      proposed$max_global_ucb,
      componentwise$max_global_ucb
    ),
    calibration_max_PRN_all = c(
      max_all_estimate[proposed$A_index],
      max_all_estimate[componentwise$A_index]
    ),
    calibration_max_PRN_positive = c(
      max_positive_estimate[
        proposed$A_index, proposed$P_index
      ],
      max_positive_estimate[
        componentwise$A_index, componentwise$P_index
      ]
    ),
    stringsAsFactors = FALSE
  )

  list(
    setting = setting,
    endpoint_number = endpoint_number,
    selected = selected,
    proposed = proposed,
    componentwise = componentwise
  )
}

evaluate_selected_method <- function(
    data, lambda_A_value, lambda_pos_value
) {
  all_matrices <- build_allcomer_candidate_matrices(
    data, lambda_A_value
  )
  success1 <- build_positive_success_matrix(
    data, entry_look = 1L, lambdas = lambda_pos_value
  )[, 1L]
  success2 <- build_positive_success_matrix(
    data, entry_look = 2L, lambdas = lambda_pos_value
  )[, 1L]

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

validate_endpoint <- function(calibration_result) {
  setting <- calibration_result$setting
  endpoint_number <- calibration_result$endpoint_number
  selected <- calibration_result$selected
  rows <- vector(
    "list",
    length(prevalence_grid) * nrow(selected)
  )
  out_idx <- 1L

  message("Validating endpoint: ", setting$endpoint_label)

  for (p_idx in seq_along(prevalence_grid)) {
    pi <- prevalence_grid[p_idx]
    data <- simulate_potential_data(
      n_rep = n_validation,
      pi = pi,
      theta0 = setting$theta0,
      endpoint_id = setting$endpoint_id,
      seed = base_seed +
        endpoint_number * 1000000L +
        p_idx * 10000L + 5000L
    )

    for (m_idx in seq_len(nrow(selected))) {
      estimates <- evaluate_selected_method(
        data = data,
        lambda_A_value = selected$lambda_A[m_idx],
        lambda_pos_value =
          selected$lambda_positive[m_idx]
      )
      ci_any <- wilson_interval(
        estimates["PRN_any"], n_validation
      )

      rows[[out_idx]] <- data.frame(
        Endpoint = setting$endpoint_label,
        Endpoint_ID = setting$endpoint_id,
        Method = selected$Method[m_idx],
        pi = pi,
        Monte_Carlo_replicates = n_validation,
        PRN_any = unname(estimates["PRN_any"]),
        PRN_any_MCSE = sqrt(
          estimates["PRN_any"] *
            (1 - estimates["PRN_any"]) /
            n_validation
        ),
        PRN_any_CI_lower = ci_any[1L, "lower"],
        PRN_any_CI_upper = ci_any[1L, "upper"],
        PRN_all = unname(estimates["PRN_all"]),
        PRN_positive =
          unname(estimates["PRN_positive"]),
        stringsAsFactors = FALSE
      )
      out_idx <- out_idx + 1L
    }

    rm(data)
    invisible(gc())
  }

  do.call(rbind, rows)
}

# Calibration and independent validation

endpoint_jobs <- Map(
  function(setting, number) {
    list(setting = setting, number = number)
  },
  endpoint_settings,
  seq_along(endpoint_settings)
)

calibration_results <- parallel_lapply(
  endpoint_jobs,
  function(job) {
    calibrate_endpoint(
      setting = job$setting,
      endpoint_number = job$number
    )
  },
  mc_cores = n_cores
)

selected_designs <- do.call(
  rbind,
  lapply(calibration_results, function(x) x$selected)
)
rownames(selected_designs) <- NULL

validation_results <- parallel_lapply(
  calibration_results,
  validate_endpoint,
  mc_cores = n_cores
)
type1_results <- do.call(rbind, validation_results)
rownames(type1_results) <- NULL

type1_summary <- do.call(
  rbind,
  lapply(
    split(
      type1_results,
      interaction(
        type1_results$Endpoint,
        type1_results$Method,
        drop = TRUE
      )
    ),
    function(d) {
      worst <- d[which.max(d$PRN_any), , drop = FALSE]
      data.frame(
        Endpoint = d$Endpoint[1L],
        Method = d$Method[1L],
        Max_PRN_any = max(d$PRN_any),
        Prevalence_at_max = worst$pi[1L],
        Max_PRN_any_CI_upper =
          worst$PRN_any_CI_upper[1L],
        Max_PRN_all = max(d$PRN_all),
        Max_PRN_positive = max(d$PRN_positive),
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(type1_summary) <- NULL
type1_summary <- type1_summary[
  order(type1_summary$Endpoint, type1_summary$Method),
  ,
  drop = FALSE
]

# Output files

scenario_rows <- lapply(endpoint_settings, function(setting) {
  data.frame(
    Endpoint_ID = setting$endpoint_id,
    Endpoint = setting$endpoint_label,
    Categories = paste(setting$categories, collapse = "; "),
    theta0 = paste(setting$theta0, collapse = ", "),
    Monitoring_rule = setting$rule_label,
    stringsAsFactors = FALSE
  )
})
null_scenarios <- do.call(rbind, scenario_rows)
rownames(null_scenarios) <- NULL

settings_table <- data.frame(
  Parameter = c(
    "alpha",
    "all-comer interim looks",
    "post-enrichment positive looks",
    "N_A",
    "N_positive",
    "n_positive_min",
    "prevalence grid",
    "Dirichlet prior",
    "gamma_A",
    "gamma_positive",
    "lambda grid",
    "calibration replicates",
    "independent validation replicates",
    "calibration confidence bound",
    "base seed",
    "parallel workers"
  ),
  Value = c(
    alpha,
    paste(interim_n, collapse = ", "),
    paste(post_enrichment_n, collapse = ", "),
    N_A,
    N_pos,
    n_pos_min,
    paste(prevalence_grid, collapse = ", "),
    paste(dirichlet_prior, collapse = ", "),
    gamma_A,
    gamma_pos,
    paste0(
      min(lambda_grid), " to ", max(lambda_grid),
      " by ",
      if (length(lambda_grid) > 1L) {
        lambda_grid[2L] - lambda_grid[1L]
      } else {
        NA
      }
    ),
    n_calibration,
    n_validation,
    "One-sided 95% Wilson upper bound",
    base_seed,
    n_cores
  ),
  stringsAsFactors = FALSE
)

write.csv(
  settings_table,
  file.path(result_dir, "complex_endpoint_type1_settings.csv"),
  row.names = FALSE
)
write.csv(
  null_scenarios,
  file.path(result_dir, "complex_endpoint_null_scenarios.csv"),
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

saveRDS(
  list(
    settings = settings_table,
    null_scenarios = null_scenarios,
    selected_designs = selected_designs,
    type1_results = type1_results,
    type1_summary = type1_summary
  ),
  file.path(result_dir, "complex_endpoint_type1_results.rds")
)

xlsx_path <- file.path(
  result_dir, "complex_endpoint_type1_results.xlsx"
)
xlsx_created <- write_results_workbook(
  sheets = list(
    Settings = settings_table,
    Null_scenarios = null_scenarios,
    Selected_designs = selected_designs,
    Type1_by_prevalence = type1_results,
    Type1_summary = type1_summary
  ),
  path = xlsx_path
)

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo.txt")
)

message("Selected designs:")
print(selected_designs)
message("Independent validation summary:")
print(type1_summary)
message(
  "Completed. Results were written to: ",
  normalizePath(result_dir, mustWork = FALSE)
)
if (xlsx_created) {
  message(
    "XLSX workbook: ",
    normalizePath(xlsx_path, mustWork = FALSE)
  )
}
