# Reproduce the conventional BOP2 binary-endpoint design and calculate
# BOP2 boundaries for the enrichment-trial analysis schedules.
#
# Selection rule follows the numerical-study convention of Zhou, Lee,
# and Yuan (2017, Section 4; Table I):
#   1. retain designs with exact type I error <= alpha;
#   2. among them, retain the attainable type I error closest to alpha;
#   3. among ties, maximize exact power under p1;
#   4. among remaining ties, minimize EN under p0 and then choose the first.
#
# The script creates:
#   bop2_reference_reproduction.csv
#   bop2_allcomer_current_schedule.csv
#   bop2_positive_entry_specific_library.csv
#   bop2_positive_entry_specific_boundaries.csv

options(stringsAsFactors = FALSE, scipen = 999)

## -----------------------------------------------------------------------------
## 1. User settings
## -----------------------------------------------------------------------------

p0 <- 0.20
p1 <- 0.40
alpha <- 0.10
N <- 40L
prior_a <- 0.5
prior_b <- 0.5

# Current enrichment-trial schedules.
allcomer_looks <- c(20L, 30L, 40L)
positive_first_assessment_sizes <- 10L:30L
post_enrichment_looks <- c(20L, 30L)

# Fine grid similar to the user-provided BOP2 reference program.
lambda_grid <- seq(0.001, 1.000, by = 0.0025)
gamma_grid <- seq(0.010, 10.000, by = 0.025)

# Run from the repository root. Outputs are generated below results/revision/
# (ignored by Git). Set ENRICHMENT_BOP2_PROJECT_DIR to override.
default_output_dir <- file.path(getwd(), "results", "revision", "bop2_validation")
output_dir <- path.expand(Sys.getenv(
  "ENRICHMENT_BOP2_PROJECT_DIR",
  unset = default_output_dir
))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 2. Exact single-population BOP2 calculations
## -----------------------------------------------------------------------------

conv_binom <- function(dist, increment_n, probability) {
  increment <- dbinom(0:increment_n, size = increment_n, prob = probability)
  out <- numeric(length(dist) + increment_n)
  for (i in seq_along(dist)) {
    out[i:(i + increment_n)] <-
      out[i:(i + increment_n)] + dist[i] * increment
  }
  out
}

calc_bop2_oc_exact <- function(boundary, looks, probability) {
  stopifnot(
    length(boundary) == length(looks),
    all(diff(looks) > 0L)
  )

  active_distribution <- c(1)
  previous_n <- 0L
  probability_early_termination <- 0
  expected_sample_size <- 0

  for (j in seq_along(looks)) {
    current_n <- looks[j]
    increment_n <- current_n - previous_n
    distribution_now <- conv_binom(
      active_distribution,
      increment_n = increment_n,
      probability = probability
    )

    response_counts <- 0:current_n
    stop_now <- response_counts <= boundary[j]

    if (j < length(looks)) {
      stop_probability <- sum(distribution_now[stop_now])
      probability_early_termination <-
        probability_early_termination + stop_probability
      expected_sample_size <-
        expected_sample_size + current_n * stop_probability
      distribution_now[stop_now] <- 0
      active_distribution <- distribution_now
    } else {
      rejection_probability <- sum(distribution_now[!stop_now])
      expected_sample_size <-
        expected_sample_size + current_n * sum(distribution_now)
    }

    previous_n <- current_n
  }

  c(
    PRN = rejection_probability,
    PET = probability_early_termination,
    EN = expected_sample_size
  )
}

make_unique_bop2_candidates <- function(
    looks,
    N,
    p0,
    prior_a,
    prior_b,
    lambda_grid,
    gamma_grid
) {
  stopifnot(tail(looks, 1L) == N)

  posterior_probabilities <- lapply(looks, function(n) {
    x <- 0:n
    pbeta(
      p0,
      shape1 = prior_a + x,
      shape2 = prior_b + n - x
    )
  })

  grid <- expand.grid(
    lambda = lambda_grid,
    gamma = gamma_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  boundary_matrix <- matrix(
    0L,
    nrow = nrow(grid),
    ncol = length(looks),
    dimnames = list(NULL, paste0("b_", looks))
  )

  for (j in seq_along(looks)) {
    n <- looks[j]
    cutoffs <- 1 - grid$lambda * (n / N)^grid$gamma
    post_prob <- posterior_probabilities[[j]]
    boundary_matrix[, j] <- vapply(
      cutoffs,
      function(cutoff) as.integer(sum(post_prob > cutoff) - 1L),
      integer(1L)
    )
  }

  keep <- !duplicated(boundary_matrix)
  out <- cbind(
    data.frame(
      candidate_id = seq_len(sum(keep)),
      lambda = grid$lambda[keep],
      gamma = grid$gamma[keep],
      stringsAsFactors = FALSE
    ),
    as.data.frame(boundary_matrix[keep, , drop = FALSE])
  )
  rownames(out) <- NULL
  out
}

select_bop2_design <- function(
    looks,
    p0,
    p1,
    alpha,
    N,
    prior_a,
    prior_b,
    lambda_grid,
    gamma_grid,
    tolerance = 1e-12
) {
  candidates <- make_unique_bop2_candidates(
    looks = looks,
    N = N,
    p0 = p0,
    prior_a = prior_a,
    prior_b = prior_b,
    lambda_grid = lambda_grid,
    gamma_grid = gamma_grid
  )

  boundary_columns <- paste0("b_", looks)
  type1 <- numeric(nrow(candidates))
  power <- numeric(nrow(candidates))
  en_null <- numeric(nrow(candidates))
  pet_null <- numeric(nrow(candidates))

  for (i in seq_len(nrow(candidates))) {
    boundary <- as.integer(unlist(
      candidates[i, boundary_columns, drop = FALSE],
      use.names = FALSE
    ))
    oc0 <- calc_bop2_oc_exact(boundary, looks, p0)
    oc1 <- calc_bop2_oc_exact(boundary, looks, p1)
    type1[i] <- oc0["PRN"]
    power[i] <- oc1["PRN"]
    en_null[i] <- oc0["EN"]
    pet_null[i] <- oc0["PET"]
  }

  candidates$type1 <- type1
  candidates$power <- power
  candidates$PET_null <- pet_null
  candidates$EN_null <- en_null

  feasible <- which(type1 <= alpha + tolerance)
  if (length(feasible) == 0L) {
    stop("No BOP2 candidate satisfies type I error <= alpha for looks: ",
         paste(looks, collapse = ", "))
  }

  # Published numerical-study convention (Zhou et al., Section 4): choose
  # type I error closest to and not above nominal; use power to break ties.
  # This is lexicographic, not unconstrained power maximization.
  closest_type1 <- max(type1[feasible])
  best <- feasible[abs(type1[feasible] - closest_type1) <= tolerance]

  maximum_power <- max(power[best])
  best <- best[abs(power[best] - maximum_power) <= tolerance]

  minimum_en <- min(en_null[best])
  best <- best[abs(en_null[best] - minimum_en) <= tolerance]

  chosen <- best[1L]
  boundary <- as.integer(unlist(
    candidates[chosen, boundary_columns, drop = FALSE],
    use.names = FALSE
  ))
  names(boundary) <- as.character(looks)

  list(
    looks = looks,
    boundary = boundary,
    selected = candidates[chosen, , drop = FALSE],
    candidates = candidates
  )
}

summarize_design <- function(design, label) {
  data.frame(
    Design = label,
    Looks = paste(design$looks, collapse = ","),
    Boundaries = paste(design$boundary, collapse = ","),
    lambda = design$selected$lambda,
    gamma = design$selected$gamma,
    exact_type1 = design$selected$type1,
    exact_power = design$selected$power,
    PET_null = design$selected$PET_null,
    EN_null = design$selected$EN_null,
    stringsAsFactors = FALSE
  )
}

## -----------------------------------------------------------------------------
## 3. Reproduce the published binary BOP2 example
## -----------------------------------------------------------------------------

reference_looks <- c(10L, 15L, 20L, 25L, 30L, 35L, 40L)
reference_design <- select_bop2_design(
  looks = reference_looks,
  p0 = p0,
  p1 = p1,
  alpha = alpha,
  N = N,
  prior_a = prior_a,
  prior_b = prior_b,
  lambda_grid = lambda_grid,
  gamma_grid = gamma_grid
)

reference_summary <- summarize_design(
  reference_design,
  "Published BOP2 binary example"
)
write.csv(
  reference_summary,
  file.path(output_dir, "bop2_reference_reproduction.csv"),
  row.names = FALSE
)

published_boundary <- c(1L, 2L, 4L, 5L, 7L, 9L, 10L)
if (!identical(unname(reference_design$boundary), published_boundary)) {
  warning(
    "The current numerical grid did not reproduce the published boundary ",
    paste(published_boundary, collapse = ","),
    ". Observed: ", paste(reference_design$boundary, collapse = ","),
    ". Enlarge or refine lambda_grid/gamma_grid before manuscript use."
  )
}

## -----------------------------------------------------------------------------
## 4. Current all-comer BOP2 schedule
## -----------------------------------------------------------------------------

allcomer_design <- select_bop2_design(
  looks = allcomer_looks,
  p0 = p0,
  p1 = p1,
  alpha = alpha,
  N = N,
  prior_a = prior_a,
  prior_b = prior_b,
  lambda_grid = lambda_grid,
  gamma_grid = gamma_grid
)

write.csv(
  summarize_design(allcomer_design, "All-comer BOP2"),
  file.path(output_dir, "bop2_allcomer_current_schedule.csv"),
  row.names = FALSE
)

## -----------------------------------------------------------------------------
## 5. Entry-specific positive-population BOP2 library
## -----------------------------------------------------------------------------

positive_designs <- vector(
  "list",
  length(positive_first_assessment_sizes)
)
names(positive_designs) <- as.character(positive_first_assessment_sizes)
positive_summary_rows <- vector(
  "list",
  length(positive_first_assessment_sizes)
)
positive_boundary_rows <- list()
boundary_row_index <- 1L

for (i in seq_along(positive_first_assessment_sizes)) {
  first_n <- positive_first_assessment_sizes[i]
  looks <- sort(unique(c(
    first_n,
    post_enrichment_looks[post_enrichment_looks > first_n],
    N
  )))

  design <- select_bop2_design(
    looks = looks,
    p0 = p0,
    p1 = p1,
    alpha = alpha,
    N = N,
    prior_a = prior_a,
    prior_b = prior_b,
    lambda_grid = lambda_grid,
    gamma_grid = gamma_grid
  )
  positive_designs[[as.character(first_n)]] <- design

  positive_summary_rows[[i]] <- data.frame(
    first_assessment_n = first_n,
    analysis_schedule = paste(looks, collapse = ","),
    candidate_id = design$selected$candidate_id,
    lambda = design$selected$lambda,
    gamma = design$selected$gamma,
    exact_type1 = design$selected$type1,
    exact_power = design$selected$power,
    PET_null = design$selected$PET_null,
    EN_null = design$selected$EN_null,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(looks)) {
    positive_boundary_rows[[boundary_row_index]] <- data.frame(
      first_assessment_n = first_n,
      analysis_n = looks[j],
      futility_or_no_claim_boundary = design$boundary[j],
      is_final = looks[j] == N,
      stringsAsFactors = FALSE
    )
    boundary_row_index <- boundary_row_index + 1L
  }
}

positive_summary <- do.call(rbind, positive_summary_rows)
positive_boundaries <- do.call(rbind, positive_boundary_rows)

write.csv(
  positive_summary,
  file.path(output_dir, "bop2_positive_entry_specific_library.csv"),
  row.names = FALSE
)
write.csv(
  positive_boundaries,
  file.path(output_dir, "bop2_positive_entry_specific_boundaries.csv"),
  row.names = FALSE
)

cat("\nPublished-example reproduction:\n")
print(reference_summary)
cat("\nCurrent all-comer BOP2:\n")
print(summarize_design(allcomer_design, "All-comer BOP2"))
cat("\nEntry-specific positive BOP2 library:\n")
print(positive_summary)
cat("\nOutput directory: ", output_dir, "\n", sep = "")
