# R1-M2: reversed biomarker-response ordering (theta_positive < theta_negative)
# For the SMMR globally calibrated enrichment BOP2 manuscript.
#
# Exact finite-state evaluation, NOT Monte Carlo simulation.
# No boundary recalibration; the published (p0,p1)=(0.20,0.40) boundaries
# are fixed throughout. Both the proposed and the independent conventional
# BOP2 comparator are evaluated under the identical no-bridging procedure.
#
# Only base R and the recommended 'parallel' package are used.
# Run from the desired project folder:
#   Sys.setenv(BOP2_N_CORES = "16")
#   source("R1M2_reverse_order_exact.R")
# Optional: BOP2_R1M2_OUTPUT_DIR specifies a NEW results directory.

## 1. Fixed design specifications

p0 <- 0.20
p1 <- 0.40
prevalence <- c(0.4, 0.5, 0.6, 0.7, 0.8)
all_looks <- c(20L, 30L, 40L)
positive_looks <- c(20L, 30L, 40L)
min_positive <- 10L

# Fixed proposed-design boundaries from the original manuscript.
# b(n) is the largest response count indicating futility or no claim.
# Positive boundaries for entry m=10,...,30 and final n=40.
proposed_A <- c(`20` = 4L, `30` = 8L, `40` = 11L)
proposed_P <- rep(NA_integer_, 41L)  # element n+1 is boundary at sample size n
proposed_P[11:13] <- 1L              # n=10,...,12
proposed_P[14:17] <- 2L              # n=13,...,16
proposed_P[18:20] <- 3L              # n=17,...,19
proposed_P[21:24] <- 4L              # n=20,...,23
proposed_P[25:27] <- 5L              # n=24,...,26
proposed_P[28:30] <- 6L              # n=27,...,29
proposed_P[31] <- 7L                 # n=30
proposed_P[41] <- 11L                # n=40

# Conventional BOP2 comparator, independently calibrated before embedding in
# the same adaptive procedure. This is the ORIGINAL comparator, not an
# error-matched or multiplicity-adjusted comparator. The positive boundaries
# depend on the first-assessment sample size, m=10,...,30.
comparator_A <- c(`20` = 5L, `30` = 7L, `40` = 10L)
comparator_library <- data.frame(
  entry_n = 10:30,
  b_entry = c(2,2,3,3,3,4,3,4,4,4,5,-1,5,6,0,6,6,7,7,7,8),
  b_20 = c(5,4,5,5,4,5,5,4,5,5,5,-1,5,5,-1,5,5,5,5,5,5),
  b_30 = c(7,7,7,7,7,7,7,7,7,7,7,2,8,7,2,8,8,8,8,8,8),
  b_40 = c(9,10,9,9,10,9,10,10,10,10,10,11,10,10,11,10,10,10,10,10,10)
)

stopifnot(
  all(names(proposed_A) == c("20", "30", "40")),
  all(names(comparator_A) == c("20", "30", "40")),
  nrow(comparator_library) == 21L,
  identical(comparator_library$entry_n, 10:30),
  all(is.finite(proposed_P[c(11:31, 41)]))
)

make_proposed_P <- function(entry_n) proposed_P

make_comparator_P <- function(entry_n) {
  z <- comparator_library[comparator_library$entry_n == entry_n, , drop = FALSE]
  if (nrow(z) != 1L) stop("Missing comparator entry sample size: ", entry_n)
  b <- rep(NA_integer_, 41L)
  b[entry_n + 1L] <- z$b_entry
  if (entry_n < 20L) b[21L] <- z$b_20
  if (entry_n < 30L) b[31L] <- z$b_30
  b[41L] <- z$b_40
  b
}

designs <- list(
  Proposed = list(A = proposed_A, getP = make_proposed_P),
  Independent_BOP2 = list(A = comparator_A, getP = make_comparator_P)
)

## 2. Exact conditional success on the enrichment branch

# The argument (m,x) denotes the biomarker-positive sample size and number
# of responses at the FIRST assessment following all-comer futility.
# The first assessment is performed exactly once. Previously passed scheduled
# positive-population looks are skipped. If m<10, enrichment is impossible.
# q(n,x) is the success probability AFTER the first assessment is passed.
make_positive_success_kernel <- function(design, theta_positive) {
  n_states <- 31L * 32L / 2L   # states (m,x) with 0<=x<=m<=30
  kernel <- numeric(n_states)

  state_id <- function(m, x) as.integer(m * (m + 1L) / 2L + x + 1L)

  for (entry_n in min_positive:30L) {
    b <- design$getP(entry_n)
    if (is.na(b[entry_n + 1L]) || is.na(b[41L])) {
      stop("Missing positive boundary for entry n=", entry_n)
    }

    memo <- new.env(hash = TRUE, parent = emptyenv())
    active_success <- function(n_current, x_current) {
      key <- paste(n_current, x_current, sep = ":")
      if (exists(key, envir = memo, inherits = FALSE)) {
        return(get(key, envir = memo, inherits = FALSE))
      }
      next_n <- min(positive_looks[positive_looks > n_current])
      d <- next_n - n_current
      u <- 0:d
      pu <- dbinom(u, size = d, prob = theta_positive)
      x_next <- x_current + u

      if (next_n == 40L) {
        ans <- sum(pu[x_next > b[41L]])
      } else {
        cutoff <- b[next_n + 1L]
        if (is.na(cutoff)) stop("Missing scheduled biomarker-positive boundary")
        pass <- which(x_next > cutoff)
        if (length(pass) == 0L) {
          ans <- 0
        } else {
          descendants <- vapply(
            x_next[pass],
            function(x) active_success(next_n, x),
            numeric(1L)
          )
          ans <- sum(pu[pass] * descendants)
        }
      }
      assign(key, ans, envir = memo)
      ans
    }

    for (x in 0:entry_n) {
      if (x > b[entry_n + 1L]) {
        kernel[state_id(entry_n, x)] <- active_success(entry_n, x)
      }
    }
  }
  if (any(!is.finite(kernel)) || any(kernel < -1e-12 | kernel > 1 + 1e-12)) {
    stop("Invalid positive success probabilities")
  }
  kernel
}

## 3. Exact all-comer state probabilities and claim probabilities

# Rows represent (m, x_pos, x_neg) among n all-comer patients. These
# probabilities correspond to four independent multinomial categories:
# positive responder/nonresponder and negative responder/nonresponder.
make_joint_states <- function(n, pi, theta_positive, theta_negative) {
  blocks <- lapply(0:n, function(m) {
    n_neg <- n - m
    data.frame(
      m = rep.int(as.integer(m), (m + 1L) * (n_neg + 1L)),
      x_pos = rep(0:m, each = n_neg + 1L),
      x_neg = rep(0:n_neg, times = m + 1L)
    )
  })
  st <- do.call(rbind, blocks)
  st$x_all <- st$x_pos + st$x_neg
  st$prob <- dbinom(st$m, size = n, prob = pi) *
    dbinom(st$x_pos, size = st$m, prob = theta_positive) *
    dbinom(st$x_neg, size = n - st$m, prob = theta_negative)
  if (abs(sum(st$prob) - 1) > 1e-11) {
    stop("Initial joint-state probabilities do not sum to one (n=", n, ")")
  }
  st
}

# Hypergeometric identity: conditional on the responses among 30 all-comer
# patients, the number of responses in the FIRST 20 patients has a
# Hypergeometric(total=30, successes=x_all, draw=20) distribution.
# This exactly accounts for the selection induced by first-look survival.
# It is equivalent to the cohort-by-cohort convolution in the submitted code.
make_allcomer_crossing_weights <- function(st20, st30, A, theta_all) {
  b20 <- unname(A["20"])
  b30 <- unname(A["30"])
  b40 <- unname(A["40"])
  if (anyNA(c(b20, b30, b40))) stop("Missing all-comer boundary")

  cross20 <- st20$x_all <= b20
  # phyper(q,m,n,k,lower.tail=FALSE) = P(X>q)
  p_survive20_given30 <- phyper(
    q = b20,
    m = st30$x_all,
    n = 30L - st30$x_all,
    k = 20L,
    lower.tail = FALSE
  )
  at30_after_surviving20 <- st30$prob * p_survive20_given30
  cross30 <- st30$x_all <= b30
  survive30 <- !cross30

  weights <- numeric(31L * 32L / 2L)
  add_weights <- function(st, mass) {
    if (length(mass) == 0L) return(invisible(NULL))
    ids <- as.integer(st$m * (st$m + 1L) / 2L + st$x_pos + 1L)
    groups <- rowsum(mass, group = ids, reorder = FALSE)
    ix <- as.integer(rownames(groups))
    weights[ix] <<- weights[ix] + groups[, 1L]
    invisible(NULL)
  }
  add_weights(st20[cross20, , drop = FALSE], st20$prob[cross20])
  add_weights(st30[cross30, , drop = FALSE], at30_after_surviving20[cross30])

  p_reach_final <- sum(at30_after_surviving20[survive30])
  p_all_claim <- sum(
    at30_after_surviving20[survive30] *
      pbinom(b40 - st30$x_all[survive30],
             size = 10L, prob = theta_all, lower.tail = FALSE)
  )

  if (abs(sum(weights) + p_reach_final - 1) > 2e-10) {
    stop("All-comer route probabilities fail to sum to one")
  }
  list(weights = weights, PRN_all = p_all_claim)
}

## 4. Prespecified scenarios

scenario_specs <- data.frame(
  scenario = c("Point global null", "Working alternative",
               "Reverse: positive below null", "Reverse: positive at null"),
  theta_positive = c(0.20, 0.40, 0.10, 0.20),
  theta_negative = c(0.20, 0.20, 0.40, 0.40),
  stringsAsFactors = FALSE
)
scenarios <- merge(
  scenario_specs,
  data.frame(pi = prevalence),
  by = NULL,
  sort = FALSE
)
scenarios <- scenarios[order(match(scenarios$scenario,
                                   scenario_specs$scenario), scenarios$pi), ]
rownames(scenarios) <- NULL

# Precompute the biomarker-positive success kernels only once for each
# method and distinct theta_positive; reuse across all prevalence values.
positive_kernels <- list()
for (method in names(designs)) {
  positive_kernels[[method]] <- list()
  for (tp in sort(unique(scenarios$theta_positive))) {
    key <- sprintf("%.2f", tp)
    message("Precomputing positive kernel: ", method, ", theta_positive=", key)
    positive_kernels[[method]][[key]] <-
      make_positive_success_kernel(designs[[method]], tp)
  }
}

one_scenario <- function(i) {
  s <- scenarios[i, , drop = FALSE]
  theta_all <- s$pi * s$theta_positive + (1 - s$pi) * s$theta_negative
  st20 <- make_joint_states(20L, s$pi, s$theta_positive, s$theta_negative)
  st30 <- make_joint_states(30L, s$pi, s$theta_positive, s$theta_negative)
  rows <- lapply(names(designs), function(method) {
    component <- make_allcomer_crossing_weights(
      st20, st30, designs[[method]]$A, theta_all
    )
    key <- sprintf("%.2f", s$theta_positive)
    p_pos <- sum(component$weights * positive_kernels[[method]][[key]])
    p_all <- component$PRN_all
    p_any <- p_all + p_pos  # mutually exclusive trial branches

    # These are scenario-specific correctness and error probabilities, NOT
    # a claim of strong error control over the entire composite null space.
    A_is_null <- theta_all <= p0 + 1e-12
    P_is_null <- s$theta_positive <= p0 + 1e-12
    p_false <- (if (A_is_null) p_all else 0) +
      (if (P_is_null) p_pos else 0)
    p_correct <- (if (!A_is_null) p_all else 0) +
      (if (!P_is_null) p_pos else 0)

    if (p_any < -1e-12 || p_any > 1 + 1e-10 ||
        abs(p_any - p_false - p_correct) > 1e-10) {
      stop("Invalid PRN partition for scenario ", i, ", method ", method)
    }
    data.frame(
      scenario = s$scenario,
      pi = s$pi,
      theta_positive = s$theta_positive,
      theta_negative = s$theta_negative,
      theta_all = theta_all,
      method = method,
      PRN_any = p_any,
      PRN_all = p_all,
      PRN_positive = p_pos,
      all_comer_null_true = A_is_null,
      positive_null_true = P_is_null,
      pointwise_FWER = p_false,
      probability_valid_claim = p_correct,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

## 5. Run in parallel, validate against the submitted results, and save CSVs

n_detect <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (length(n_detect) != 1L || is.na(n_detect) || n_detect < 1L) n_detect <- 1L
n_requested <- suppressWarnings(as.integer(Sys.getenv("BOP2_N_CORES", "16")))
if (length(n_requested) != 1L || is.na(n_requested) || n_requested < 1L) {
  n_requested <- 16L
}
n_workers <- min(n_requested, n_detect, nrow(scenarios))
if (.Platform$OS.type == "windows") n_workers <- 1L
message("Detected physical cores: ", n_detect, "; scenario workers: ", n_workers)
message("Evaluating ", nrow(scenarios), " scenarios x ", length(designs),
        " methods by exact recursion...")

results_list <- if (n_workers > 1L) {
  parallel::mclapply(seq_len(nrow(scenarios)), one_scenario,
                     mc.cores = n_workers, mc.preschedule = TRUE,
                     mc.set.seed = FALSE)
} else {
  lapply(seq_len(nrow(scenarios)), one_scenario)
}
if (any(vapply(results_list, inherits, logical(1), what = "try-error"))) {
  stop("At least one scenario failed. Check earlier messages.")
}
results <- do.call(rbind, results_list)
rownames(results) <- NULL

# Strong regression checks: exact values of the originally submitted design.
check_reference <- function(scenario_name, pi_value, method, expected, tol = 2e-9) {
  z <- results[results$scenario == scenario_name &
                 abs(results$pi - pi_value) < 1e-12 &
                 results$method == method, , drop = FALSE]
  if (nrow(z) != 1L || abs(z$PRN_any - expected) > tol) {
    stop("Baseline regression check FAILED for ", method, ", ",
         scenario_name, ", pi=", pi_value, ". Expected ", expected,
         "; got ", if (nrow(z) == 1L) z$PRN_any else "missing")
  }
}
check_reference("Point global null", 0.6, "Proposed", 0.09894542782786911)
check_reference("Point global null", 0.6, "Independent_BOP2", 0.12633818863934865)
check_reference("Working alternative", 0.4, "Proposed", 0.5711563063010007)
check_reference("Working alternative", 0.4, "Independent_BOP2", 0.4844503298531837)
message("Baseline regression checks: ALL PASSED")

output_parent <- Sys.getenv("BOP2_R1M2_OUTPUT_DIR", unset = "")
if (!nzchar(output_parent)) {
  output_parent <- file.path(
    getwd(),
    paste0("R1M2_results_", format(Sys.time(), "%Y%m%d_%H%M%S"),
           "_", Sys.getpid())
  )
}
if (dir.exists(output_parent) || file.exists(output_parent)) {
  stop("Output path already exists; choose a new output directory: ", output_parent)
}
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)

reverse_rows <- grepl("^Reverse:", results$scenario)
reverse_results <- results[reverse_rows, , drop = FALSE]

paper <- reverse_results[, c("scenario", "pi", "theta_positive",
                             "theta_negative", "theta_all", "method")]
paper$PRN_any_pct <- round(reverse_results$PRN_any * 100, 4)
paper$PRN_all_pct <- round(reverse_results$PRN_all * 100, 4)
paper$PRN_positive_pct <- round(reverse_results$PRN_positive * 100, 4)
paper$pointwise_FWER_pct <- round(reverse_results$pointwise_FWER * 100, 4)

write.csv(results, file.path(output_parent, "R1M2_all_exact.csv"),
          row.names = FALSE)
write.csv(reverse_results, file.path(output_parent, "R1M2_reverse_only_exact.csv"),
          row.names = FALSE)
write.csv(paper, file.path(output_parent, "R1M2_reverse_table_percent.csv"),
          row.names = FALSE)
write.csv(comparator_library,
          file.path(output_parent, "fixed_comparator_positive_library.csv"),
          row.names = FALSE)

message("\n=== R1-M2 reverse-order scenarios: results (%) ===")
print(paper, row.names = FALSE, digits = 6)
message("\nAll results saved to: ", normalizePath(output_parent))
message("These are exact evaluation results for FIXED boundaries, not new calibrations.")
message("Pointwise FWER is scenario-specific; it is NOT a strong-control guarantee.")
