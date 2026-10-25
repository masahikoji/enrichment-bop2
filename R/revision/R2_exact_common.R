# Shared exact finite-state evaluator for Reviewer 2 Comments 2--4.
# Source this file from the same directory as the two analysis scripts.
# Base R only. No Monte Carlo simulation.
# SMMR manuscript: globally calibrated BOP2 adaptive enrichment.

R2_p0 <- 0.20
R2_p1 <- 0.40
R2_prevalences <- c(0.4, 0.5, 0.6, 0.7, 0.8)
R2_alpha <- 0.10
R2_all_looks <- c(20L, 30L, 40L)
R2_positive_looks <- c(20L, 30L, 40L)
R2_min_positive <- 10L
R2_final_n <- 40L

R2_proposed_A <- c(`20` = 4L, `30` = 8L, `40` = 11L)
R2_proposed_P <- rep(NA_integer_, 41L) # n+1 stores boundary at n
R2_proposed_P[11:13] <- 1L
R2_proposed_P[14:17] <- 2L
R2_proposed_P[18:20] <- 3L
R2_proposed_P[21:24] <- 4L
R2_proposed_P[25:27] <- 5L
R2_proposed_P[28:30] <- 6L
R2_proposed_P[31] <- 7L
R2_proposed_P[41] <- 11L

# From the submitted independent (unadjusted) BOP2 calibration at alpha=0.10.
R2_independent_A <- c(`20` = 5L, `30` = 7L, `40` = 10L)
R2_independent_library <- data.frame(
  entry_n = 10:30,
  b_entry = c(2,2,3,3,3,4,3,4,4,4,5,-1,5,6,0,6,6,7,7,7,8),
  b_20 = c(5,4,5,5,4,5,5,4,5,5,5,-1,5,5,-1,5,5,5,5,5,5),
  b_30 = c(7,7,7,7,7,7,7,7,7,7,7,2,8,7,2,8,8,8,8,8,8),
  b_40 = c(9,10,9,9,10,9,10,10,10,10,10,11,10,10,11,10,10,10,10,10,10)
)
stopifnot(identical(R2_independent_library$entry_n, 10:30))
# Regression guard: entry m=17 must use b20=4 (not 5).
stopifnot(R2_independent_library$b_20[R2_independent_library$entry_n == 17L] == 4L)

R2_make_design <- function(A, library = NULL, shared_P = NULL, name = "Unnamed") {
  stopifnot(length(A) == 3L, all(names(A) == c("20", "30", "40")))
  if (is.null(library) == is.null(shared_P)) stop("Specify exactly one of library or shared_P")
  if (!is.null(library)) {
    stopifnot(nrow(library) == 21L, identical(as.integer(library$entry_n), 10:30))
  }
  list(name = name, A = A, library = library, shared_P = shared_P)
}
R2_proposed_design <- R2_make_design(R2_proposed_A, shared_P = R2_proposed_P, name = "Proposed")
R2_independent_design <- R2_make_design(R2_independent_A, library = R2_independent_library, name = "Independent_BOP2_unadjusted")

# Probability that a standalone BOP2 reaches its final analysis and rejects.
R2_standalone_prob <- function(boundaries, looks, theta) {
  active <- 1
  old_n <- 0L
  for (j in seq_along(looks)) {
    n <- looks[j]
    d <- n - old_n
    pmf <- dbinom(0:d, d, theta)
    next_active <- numeric(n + 1L)
    for (u in 0:d) {
      next_active[(u + 1L):(u + old_n + 1L)] <- next_active[(u + 1L):(u + old_n + 1L)] + active * pmf[u + 1L]
    }
    next_active[(0:n) <= boundaries[j]] <- 0
    active <- next_active
    old_n <- n
  }
  sum(active)
}

# Evaluate all candidate boundary tables at once for standalone null and power.
R2_standalone_prob_matrix <- function(B, looks, theta) {
  B <- as.matrix(B)
  nc <- nrow(B)
  active <- matrix(1, nrow = nc, ncol = 1L)
  old_n <- 0L
  for (j in seq_along(looks)) {
    n <- as.integer(looks[j])
    d <- n - old_n
    pmf <- dbinom(0:d, d, theta)
    next_active <- matrix(0, nrow = nc, ncol = n + 1L)
    for (u in 0:d) {
      cols <- (u + 1L):(u + old_n + 1L)
      next_active[, cols] <- next_active[, cols, drop = FALSE] + active * pmf[u + 1L]
    }
    for (x in 0:n) {
      next_active[, x + 1L] <- next_active[, x + 1L] * as.numeric(x > B[, j])
    }
    active <- next_active
    old_n <- n
  }
  rowSums(active)
}

# First assessment at m is performed exactly once, and passed interim looks are skipped.
R2_positive_success_at_entry <- function(entry_n, cutoffs, theta) {
  looks <- c(entry_n, c(20L, 30L, 40L)[c(20L, 30L, 40L) > entry_n])
  stopifnot(length(cutoffs) == length(looks), tail(looks, 1L) == 40L)
  q <- as.numeric((0:40) > tail(cutoffs, 1L))
  for (j in rev(seq_len(length(looks) - 1L))) {
    n <- looks[j]
    n_next <- looks[j + 1L]
    d <- n_next - n
    pmf <- dbinom(0:d, size = d, prob = theta)
    q_next <- vapply(0:n, function(x) sum(q[(x + 0:d) + 1L] * pmf), numeric(1L))
    q_next[(0:n) <= cutoffs[j]] <- 0
    q <- q_next
  }
  stopifnot(length(q) == entry_n + 1L, all(is.finite(q)))
  q
}

R2_positive_kernel <- function(design, theta) {
  out <- numeric(31L * 32L / 2L) # (m,x) with m=0,...,30 and x=0,...,m
  for (m in 10:30) {
    looks <- c(m, c(20L, 30L, 40L)[c(20L, 30L, 40L) > m])
    if (!is.null(design$shared_P)) {
      cutoffs <- as.integer(design$shared_P[looks + 1L])
    } else {
      row <- design$library[design$library$entry_n == m, , drop = FALSE]
      if (nrow(row) != 1L) stop("Entry-specific boundary missing at m=", m)
      cutoffs <- vapply(looks, function(n) {
        if (n == m) return(as.integer(row$b_entry))
        if (n == 20L) return(as.integer(row$b_20))
        if (n == 30L) return(as.integer(row$b_30))
        as.integer(row$b_40)
      }, integer(1L))
    }
    if (anyNA(cutoffs)) stop("Missing biomarker-positive cutoff for entry=", m)
    ids <- as.integer(m * (m + 1L) / 2L + (0:m) + 1L)
    out[ids] <- R2_positive_success_at_entry(m, cutoffs, theta)
  }
  out
}

R2_joint_states <- function(n, pi, tp, tn) {
  pieces <- lapply(0:n, function(m) {
    neg_n <- n - m
    data.frame(m = rep.int(m, (m+1L) * (neg_n+1L)),
               x_pos = rep(0:m, each = neg_n+1L),
               x_neg = rep(0:neg_n, times = m+1L))
  })
  st <- do.call(rbind, pieces)
  st$x_all <- st$x_pos + st$x_neg
  st$prob <- dbinom(st$m, n, pi) * dbinom(st$x_pos, st$m, tp) * dbinom(st$x_neg, n-st$m, tn)
  if (abs(sum(st$prob)-1) > 1e-11) stop("Joint state probabilities do not sum to one")
  st
}

# Hypergeometric exchangeability identity accounts for the fact that the 30-patient
# state must have survived the first all-comer futility assessment at n=20.
R2_all_components <- function(A, st20, st30, theta_all) {
  b20 <- as.integer(A[["20"]])
  b30 <- as.integer(A[["30"]])
  b40 <- as.integer(A[["40"]])
  cross20 <- st20$x_all <= b20
  survive20_given30 <- phyper(b20, m = st30$x_all, n = 30L - st30$x_all, k = 20L, lower.tail = FALSE)
  mass30 <- st30$prob * survive20_given30
  cross30 <- st30$x_all <= b30
  final30 <- !cross30
  weights <- numeric(31L * 32L / 2L)
  add_mass <- function(state, probabilities) {
    if (!nrow(state)) return(invisible(NULL))
    key <- as.integer(state$m * (state$m + 1L) / 2L + state$x_pos + 1L)
    agg <- rowsum(probabilities, group = key, reorder = FALSE)
    ids <- as.integer(rownames(agg))
    weights[ids] <<- weights[ids] + as.numeric(agg[, 1L])
    invisible(NULL)
  }
  add_mass(st20[cross20, , drop = FALSE], st20$prob[cross20])
  add_mass(st30[cross30, , drop = FALSE], mass30[cross30])
  p_all <- sum(mass30[final30] * pbinom(b40 - st30$x_all[final30], size = 10L, prob = theta_all, lower.tail = FALSE))
  p_all_final <- sum(mass30[final30])
  if (abs(sum(weights) + p_all_final - 1) > 2e-10) stop("Invalid all-comer branch mass partition")
  list(weights = weights, PRN_all = p_all)
}

R2_classification <- function(theta_all, tp, tn, p0=R2_p0) {
  list(null_A = (theta_all <= p0 + 1e-12),
       null_positive = (tp <= p0 + 1e-12),
       null_negative = (tn <= p0 + 1e-12),
       subgroup_only = (tp > p0 + 1e-12 && theta_all <= p0 + 1e-12),
       enrichment_preferred = (tp > p0 + 1e-12 && tn <= p0 + 1e-12))
}

R2_evaluate_scenario <- function(pi, tp, tn, design, kernel = NULL, st20 = NULL, st30 = NULL) {
  ta <- pi*tp+(1-pi)*tn
  if (is.null(kernel)) kernel <- R2_positive_kernel(design, tp)
  if (is.null(st20)) st20 <- R2_joint_states(20L, pi, tp, tn)
  if (is.null(st30)) st30 <- R2_joint_states(30L, pi, tp, tn)
  comp <- R2_all_components(design$A, st20, st30, ta)
  pA <- as.numeric(comp$PRN_all)
  pP <- as.numeric(sum(comp$weights * kernel))
  pAny <- pA + pP
  cl <- R2_classification(ta, tp, tn)
  pFalseA <- if (cl$null_A) pA else 0
  pFalseP <- if (cl$null_positive) pP else 0
  pValidA <- if (!cl$null_A) pA else 0
  pValidP <- if (!cl$null_positive) pP else 0
  stopifnot(is.finite(pAny), pAny >= -1e-12, pAny <= 1+1e-10,
            abs(pFalseA+pFalseP+pValidA+pValidP-pAny) < 1e-10)
  data.frame(method = design$name, pi=pi, theta_positive=tp, theta_negative=tn, theta_all=ta,
             null_allcomer_true=cl$null_A, null_positive_true=cl$null_positive,
             null_negative_true=cl$null_negative, subgroup_only=cl$subgroup_only,
             positive_target_preferred=cl$enrichment_preferred,
             PRN_any=pAny, PRN_all=pA, PRN_positive=pP,
             P_false_allcomer=pFalseA, P_false_positive=pFalseP,
             pointwise_FWER=pFalseA+pFalseP,
             P_valid_allcomer=pValidA, P_valid_positive=pValidP,
             P_at_least_one_valid_claim=pValidA+pValidP,
             P_correct_positive_target=if (cl$enrichment_preferred) pP else NA_real_,
             P_potentially_overbroad_allcomer=if (cl$enrichment_preferred) pA else NA_real_,
             stringsAsFactors = FALSE)
}

R2_scenarios <- function() {
  original <- expand.grid(pi=R2_prevalences, theta_positive=c(.4,.5,.6), theta_negative=c(.2,.1),
                          KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  original$scenario <- "Original_30_alternatives"
  subgroup <- data.frame(pi=c(.4,.4), theta_positive=c(.3,.4), theta_negative=c(.05,.05),
                         scenario="Subgroup_only", stringsAsFactors=FALSE)
  null <- data.frame(pi=R2_prevalences, theta_positive=.2, theta_negative=.2,
                     scenario="Point_global_null", stringsAsFactors=FALSE)
  reverse <- expand.grid(pi=R2_prevalences, theta_positive=c(.1,.2), theta_negative=.4,
                         KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  reverse$scenario <- "Reverse_order_R1M2"
  rbind(null, original, subgroup, reverse)
}

R2_new_output_dir <- function(prefix) {
  output_root <- path.expand(Sys.getenv("BOP2_REVISION_RESULTS_ROOT", unset=getwd()))
  dir.create(output_root, recursive=TRUE, showWarnings=FALSE)
  p <- file.path(output_root, paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  if (file.exists(p) || dir.exists(p)) stop("Output directory already exists: ", p)
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
  normalizePath(p, winslash="/", mustWork=TRUE)
}

R2_verify_original <- function(tolerance=3e-9) {
  scenarios <- data.frame(pi=c(.6,.6,.4,.4), tp=c(.2,.2,.4,.4), tn=c(.2,.2,.2,.2))
  methods <- c("Proposed","Independent_BOP2_unadjusted","Proposed","Independent_BOP2_unadjusted")
  expect <- c(.09894542782786911,.12633818863934865,.5711563063010007,.4844503298531837)
  for (i in seq_len(nrow(scenarios))) {
    d <- if (methods[i] == "Proposed") R2_proposed_design else R2_independent_design
    r <- R2_evaluate_scenario(scenarios$pi[i],scenarios$tp[i],scenarios$tn[i],d)
    if (abs(r$PRN_any - expect[i]) > tolerance) {
      stop("Baseline regression test failed: ",methods[i], ", pi=",scenarios$pi[i],
           ", expected=",expect[i],", obtained=",r$PRN_any)
    }
  }
  message("Original submitted calibration and comparator regression checks: PASSED")
  invisible(TRUE)
}
