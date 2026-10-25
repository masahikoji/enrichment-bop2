# Exact calibration and timing program for the main binary simulation setting
#
# Purpose
#   1. Jointly calibrate the all-comer and biomarker-positive BOP2 boundaries.
#   2. Calculate exact global type I error and working-alternative rejection
#      probabilities over the prespecified biomarker-prevalence set.
#   3. Record wall-clock computational time and hardware/session information for
#      reporting in the Supplementary Material.
#
# The exact enumeration and selection logic follows the binary-endpoint program
# globally_calibrated_enrichment_bop2_exact_no_bridging_final.R, with the main binary-analysis settings below. No Monte Carlo simulation is used.

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## -----------------------------------------------------------------------------
## 1. Main binary-analysis settings
## -----------------------------------------------------------------------------

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
prevalence_grid <- c(0.40, 0.50, 0.60, 0.70, 0.80)

# Same proposed-design search grid as in the main binary analysis.
proposed_lambda_grid <- seq(0.005, 1.000, by = 0.005)
proposed_gamma_grid <- seq(0.050, 8.000, by = 0.050)

# No additional clinical boundary constraint is imposed in the main binary analysis.
require_zero_response_futility <- FALSE

output_dir <- path.expand(Sys.getenv(
  "BOP2_MAIN_TIMING_OUTPUT_DIR",
  unset = getwd()
))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parallel configuration. Record this explicitly because elapsed time depends on
# the number of workers. Set BOP2_N_CORES before running to override the default.
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- parallel::detectCores(logical = TRUE)
}
if (is.na(detected_cores) || detected_cores < 1L) detected_cores <- 1L

default_n_cores <- max(1L, min(8L, detected_cores - 1L))
requested_n_cores <- suppressWarnings(as.integer(Sys.getenv(
  "BOP2_N_CORES", unset = ""
)))
n_cores <- if (!is.na(requested_n_cores) && requested_n_cores >= 1L) {
  min(requested_n_cores, detected_cores)
} else {
  default_n_cores
}

parallel_backend <- if (n_cores > 1L && .Platform$OS.type != "windows") {
  "parallel::mclapply (forked)"
} else {
  "sequential lapply"
}

stopifnot(
  alpha_lower >= 0,
  alpha_lower <= alpha,
  N_A == N_pos,
  all(interim_n > 0L),
  all(diff(interim_n) > 0L),
  max(interim_n) < N_A,
  n_pos_min <= min(interim_n),
  all(post_enrichment_n > 0L),
  all(diff(post_enrichment_n) > 0L),
  all(post_enrichment_n < N_pos)
)

parallel_lapply <- function(X, FUN, ..., n_workers = n_cores) {
  if (length(X) == 0L) return(list())
  n_workers <- max(1L, min(as.integer(n_workers), length(X)))
  if (n_workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      X, FUN, ...,
      mc.cores = n_workers,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(X, FUN, ...)
  }
}

elapsed_now <- function() unname(proc.time()["elapsed"])

program_start <- elapsed_now()
section_times <- list()

## -----------------------------------------------------------------------------
## 2. Boundary generation
## -----------------------------------------------------------------------------

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
    post_prob <- pbeta(
      p0,
      shape1 = prior_a + x,
      shape2 = prior_b + s - x
    )
    cutoffs <- 1 - grid$lambda * (s / max_sample_size)^grid$gamma

    # post_prob decreases monotonically with x. Therefore the number of x values
    # satisfying posterior probability > cutoff, minus one, is the integer
    # futility boundary. A value of -1 denotes no stopping boundary.
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

candidate_start <- elapsed_now()

max_interim_pos <- max(interim_n)
A_sample_sizes <- c(interim_n, N_A)
P_sample_sizes <- sort(unique(c(
  n_pos_min:max_interim_pos,
  post_enrichment_n,
  N_pos
)))

A_candidates <- make_boundary_candidates(
  sample_sizes = A_sample_sizes,
  max_sample_size = N_A,
  lambda_values = proposed_lambda_grid,
  gamma_values = proposed_gamma_grid,
  prefix = "proposed all-comer"
)

P_candidates_unconstrained <- make_boundary_candidates(
  sample_sizes = P_sample_sizes,
  max_sample_size = N_pos,
  lambda_values = proposed_lambda_grid,
  gamma_values = proposed_gamma_grid,
  prefix = "proposed biomarker-positive"
)

P_candidates <- P_candidates_unconstrained
if (require_zero_response_futility) {
  constrained_columns <- paste0("b_", n_pos_min:max_interim_pos)
  keep_positive <- apply(
    P_candidates[, constrained_columns, drop = FALSE],
    1L,
    function(z) all(as.integer(z) >= 0L)
  )
  P_candidates <- P_candidates[keep_positive, , drop = FALSE]
  rownames(P_candidates) <- NULL
  if (nrow(P_candidates) == 0L) {
    stop("No biomarker-positive candidate satisfies b_+(n) >= 0 for n=10,...,30.")
  }
}

section_times$candidate_generation_seconds <- elapsed_now() - candidate_start

## -----------------------------------------------------------------------------
## 3. Exact all-comer state-space calculations
## -----------------------------------------------------------------------------

prob_binom_greater <- function(cutoff, size, prob) {
  out <- numeric(length(cutoff))
  out[cutoff < 0] <- 1
  out[cutoff >= size] <- 0
  idx <- cutoff >= 0 & cutoff < size
  if (any(idx)) {
    out[idx] <- pbinom(
      cutoff[idx],
      size = size,
      prob = prob,
      lower.tail = FALSE
    )
  }
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
positive_states$state_id <- state_id(
  positive_states$n_pos,
  positive_states$x_pos
)
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
      }
    }
  }

  total_probability <- rowSums(W1) + rowSums(W2) + no_cross_prob
  if (max(abs(total_probability - 1)) > 1e-9) {
    stop("Exact all-comer path probabilities do not sum to one.")
  }

  list(
    scenario = scenario,
    W1 = W1,
    W2 = W2,
    prn_all = prn_all
  )
}

## -----------------------------------------------------------------------------
## 4. Exact biomarker-positive path kernels
## -----------------------------------------------------------------------------

make_positive_boundary_vector <- function(P_candidates, P_index) {
  bP <- rep(NA_integer_, N_pos + 1L)
  bP[P_sample_sizes + 1L] <- as.integer(unlist(
    P_candidates[
      P_index,
      paste0("b_", P_sample_sizes),
      drop = FALSE
    ],
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
      stop30 = 0
    )

    if (next_look == N_pos) {
      claim <- x_next > bP[N_pos + 1L]
      out["success"] <- sum(prob_u[claim])
      out["final_no_claim"] <- sum(prob_u[!claim])
    } else {
      pass <- x_next > bP[next_look + 1L]
      if (next_look == 20L) out["stop20"] <- sum(prob_u[!pass])
      if (next_look == 30L) out["stop30"] <- sum(prob_u[!pass])

      pass_index <- which(pass)
      if (length(pass_index) > 0L) {
        for (k in pass_index) {
          child <- recurse(next_look, x_next[k])
          weight <- prob_u[k]
          out <- out + weight * child
        }
      }
    }

    assign(key, out, envir = memo)
    out
  }

  entered_states <- which(enter_vector > 0)
  for (r in entered_states) {
    ans <- recurse(
      positive_states$n_pos[r],
      positive_states$x_pos[r]
    )
    success[r] <- ans["success"]
    final_no_claim[r] <- ans["final_no_claim"]
    stop20[r] <- ans["stop20"]
    stop30[r] <- ans["stop30"]
  }

  if (length(entered_states) > 0L) {
    total <- success + final_no_claim + stop20 + stop30
    if (max(abs(total[entered_states] - 1)) > 1e-10) {
      stop("Biomarker-positive path probabilities do not sum to one.")
    }
  }

  list(success = success)
}

make_positive_kernels <- function(P_candidates, theta_pos_values) {
  boundary_columns <- paste0("b_", P_sample_sizes)
  B <- as.matrix(P_candidates[, boundary_columns, drop = FALSE])
  nP <- nrow(P_candidates)

  boundary_by_state <- matrix(
    Inf,
    nrow = n_positive_states,
    ncol = nP
  )
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

  success <- vector("list", length(theta_pos_values))
  names(success) <- format(theta_pos_values, trim = TRUE)

  for (theta_pos in theta_pos_values) {
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

  list(enter = enter * 1, success = success)
}

positive_start <- elapsed_now()
positive_kernels <- make_positive_kernels(
  P_candidates,
  theta_pos_values = c(p0, p1)
)
section_times$positive_kernel_seconds <- elapsed_now() - positive_start

## -----------------------------------------------------------------------------
## 5. Exact global calibration
## -----------------------------------------------------------------------------

null_scenarios <- data.frame(
  scenario_type = "Null",
  pi = prevalence_grid,
  theta_pos = p0,
  theta_neg = p0,
  stringsAsFactors = FALSE
)

alternative_scenarios <- data.frame(
  scenario_type = "Alternative",
  pi = prevalence_grid,
  theta_pos = p1,
  theta_neg = p0,
  stringsAsFactors = FALSE
)

all_scenarios <- rbind(null_scenarios, alternative_scenarios)
all_scenarios$scenario_id <- seq_len(nrow(all_scenarios))

allcomer_start <- elapsed_now()
scenario_rows <- split(all_scenarios, seq_len(nrow(all_scenarios)))
components <- parallel_lapply(
  scenario_rows,
  function(s) build_allcomer_components(s[1L, ], A_candidates),
  n_workers = n_cores
)
section_times$allcomer_component_seconds <- elapsed_now() - allcomer_start

calibration_prn_matrices <- function(component, positive_kernels) {
  W <- component$W1 + component$W2
  scenario <- component$scenario
  theta_key <- format(scenario$theta_pos, trim = TRUE)
  success_kernel <- positive_kernels$success[[theta_key]]
  if (is.null(success_kernel)) {
    stop("Missing biomarker-positive success kernel for theta = ", scenario$theta_pos)
  }

  prn_positive <- W %*% success_kernel
  prn_any <- sweep(
    prn_positive,
    MARGIN = 1L,
    STATS = component$prn_all,
    FUN = "+"
  )

  list(
    PRN_any = prn_any,
    PRN_all = component$prn_all,
    PRN_positive = prn_positive
  )
}

select_candidate_pair <- function(feasible, avg_power, min_power,
                                  type1_metric) {
  feasible_idx <- which(feasible)
  if (length(feasible_idx) == 0L) {
    stop("No feasible candidate pair in the prespecified near-nominal band.")
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
    A_index = ij[1L],
    P_index = ij[2L],
    average_power = avg_power[chosen],
    minimum_power = min_power[chosen],
    max_global_type1 = type1_metric[chosen]
  )
}

calibration_start <- elapsed_now()

nA <- nrow(A_candidates)
nP <- nrow(P_candidates)
max_global_type1 <- matrix(-Inf, nrow = nA, ncol = nP)
avg_power <- matrix(0, nrow = nA, ncol = nP)
min_power <- matrix(Inf, nrow = nA, ncol = nP)

for (s_idx in seq_along(components)) {
  metrics <- calibration_prn_matrices(
    components[[s_idx]],
    positive_kernels
  )
  if (all_scenarios$scenario_type[s_idx] == "Null") {
    max_global_type1 <- pmax(max_global_type1, metrics$PRN_any)
  } else {
    avg_power <- avg_power + metrics$PRN_any / length(prevalence_grid)
    min_power <- pmin(min_power, metrics$PRN_any)
  }
}

feasible_global <- max_global_type1 <= alpha + 1e-12
preferred_global <- feasible_global &
  max_global_type1 >= alpha_lower - 1e-12

selected <- select_candidate_pair(
  feasible = preferred_global,
  avg_power = avg_power,
  min_power = min_power,
  type1_metric = max_global_type1
)

section_times$joint_selection_seconds <- elapsed_now() - calibration_start

selected_A <- A_candidates[selected$A_index, , drop = FALSE]
selected_P <- P_candidates[selected$P_index, , drop = FALSE]

## -----------------------------------------------------------------------------
## 6. Exact selected-design operating characteristics
## -----------------------------------------------------------------------------

selected_oc <- vector("list", nrow(all_scenarios))
for (s_idx in seq_along(components)) {
  component <- components[[s_idx]]
  scenario <- component$scenario
  theta_key <- format(scenario$theta_pos, trim = TRUE)
  W <- component$W1[selected$A_index, ] +
    component$W2[selected$A_index, ]
  success <- positive_kernels$success[[theta_key]][, selected$P_index]

  prn_all <- component$prn_all[selected$A_index]
  prn_positive <- sum(W * success)

  selected_oc[[s_idx]] <- data.frame(
    scenario_type = scenario$scenario_type,
    pi = scenario$pi,
    theta_pos = scenario$theta_pos,
    theta_neg = scenario$theta_neg,
    PRN_any = prn_all + prn_positive,
    PRN_all = prn_all,
    PRN_positive = prn_positive,
    stringsAsFactors = FALSE
  )
}
selected_oc <- do.call(rbind, selected_oc)

## -----------------------------------------------------------------------------
## 7. Generic boundary and timing outputs for the main binary analysis
## -----------------------------------------------------------------------------

get_boundary <- function(row, n) {
  as.integer(row[[paste0("b_", n)]])
}

allcomer_boundary_table <- data.frame(
  Population = "All-comer",
  n = A_sample_sizes,
  Boundary = vapply(A_sample_sizes, function(n) get_boundary(selected_A, n), integer(1L)),
  stringsAsFactors = FALSE
)

positive_boundary_table <- data.frame(
  Population = "Biomarker-positive",
  n = P_sample_sizes,
  Boundary = vapply(P_sample_sizes, function(n) get_boundary(selected_P, n), integer(1L)),
  stringsAsFactors = FALSE
)

calibration_summary <- data.frame(
  quantity = c(
    "all-comer boundary candidates",
    "positive boundary candidates",
    "candidate pairs evaluated",
    "candidate pairs with max global type I <= alpha",
    "candidate pairs in [alpha_lower, alpha]",
    "selected lambda_A",
    "selected gamma_A",
    "selected lambda_positive",
    "selected gamma_positive",
    "selected max global type I",
    "selected average working-alternative PRN-any",
    "selected minimum working-alternative PRN-any"
  ),
  value = c(
    nrow(A_candidates),
    nrow(P_candidates),
    nrow(A_candidates) * nrow(P_candidates),
    sum(feasible_global),
    sum(preferred_global),
    selected_A$lambda,
    selected_A$gamma,
    selected_P$lambda,
    selected_P$gamma,
    selected$max_global_type1,
    selected$average_power,
    selected$minimum_power
  ),
  stringsAsFactors = FALSE
)

section_times$total_calibration_seconds <- elapsed_now() - program_start

timing <- data.frame(
  step = names(section_times),
  elapsed_seconds = as.numeric(unlist(section_times)),
  stringsAsFactors = FALSE
)

hardware <- data.frame(
  item = c(
    "date_time",
    "R_version",
    "platform",
    "OS_type",
    "detected_physical_cores",
    "workers_used",
    "parallel_backend"
  ),
  value = c(
    as.character(Sys.time()),
    R.version.string,
    R.version$platform,
    .Platform$OS.type,
    detected_cores,
    n_cores,
    parallel_backend
  ),
  stringsAsFactors = FALSE
)

if (Sys.info()[["sysname"]] == "Darwin") {
  hw_text <- tryCatch(
    system("system_profiler SPHardwareDataType", intern = TRUE),
    error = function(e) character()
  )
  wanted <- grep(
    "Chip:|Processor Name:|Total Number of Cores:|Memory:",
    hw_text,
    value = TRUE
  )
  if (length(wanted) > 0L) {
    hardware <- rbind(
      hardware,
      data.frame(
        item = paste0("mac_hardware_", seq_along(wanted)),
        value = trimws(wanted),
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(allcomer_boundary_table,
          file.path(output_dir, "main_binary_allcomer_boundaries.csv"),
          row.names = FALSE)
write.csv(positive_boundary_table,
          file.path(output_dir, "main_binary_positive_boundaries.csv"),
          row.names = FALSE)
write.csv(selected_oc,
          file.path(output_dir, "main_binary_exact_calibration_oc.csv"),
          row.names = FALSE)
write.csv(calibration_summary,
          file.path(output_dir, "main_binary_calibration_summary.csv"),
          row.names = FALSE)
write.csv(timing,
          file.path(output_dir, "main_binary_calibration_timing.csv"),
          row.names = FALSE)
write.csv(hardware,
          file.path(output_dir, "main_binary_hardware_session.csv"),
          row.names = FALSE)

cat("\nSelected all-comer parameters and boundaries:\n")
print(selected_A[, c("lambda", "gamma", paste0("b_", A_sample_sizes))])
cat("\nSelected biomarker-positive parameters and boundaries:\n")
print(selected_P[, c("lambda", "gamma", paste0("b_", P_sample_sizes))])
cat("\nExact calibration operating characteristics:\n")
print(selected_oc)
cat("\nTiming:\n")
print(timing)
cat("\nHardware/session:\n")
print(hardware)

cat("\nSupplementary-material timing sentence template:\n")
cat(
  "Calibration of the binary decision table required [TOTAL] seconds of ",
  "wall-clock time on [HARDWARE], using R [VERSION] with [CORES] worker(s). ",
  "The calculation used exact finite-state enumeration and therefore involved ",
  "no Monte Carlo error.\n",
  sep = ""
)

cat("\nOutput directory: ", output_dir, "\n", sep = "")
