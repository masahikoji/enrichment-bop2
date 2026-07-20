# Create LaTeX tables for the complex-endpoint type I error evaluation.
#
# Run through run_supplement.R from the repository root. Required inputs are
# read from results/supplement and LaTeX tables are written to
# tables/supplement.

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## -----------------------------------------------------------------------------
## 1. Paths and fixed design settings
## -----------------------------------------------------------------------------

project_dir <- path.expand(Sys.getenv(
  "ENRICHMENT_BOP2_ROOT",
  unset = normalizePath(getwd(), mustWork = TRUE)
))
quick_test <- identical(Sys.getenv("BOP2_COMPLEX_QUICK"), "1")
result_dir <- file.path(project_dir, "results", "supplement")
table_dir <- file.path(project_dir, "tables", "supplement")
if (quick_test) {
  result_dir <- file.path(result_dir, "quick_test")
  table_dir <- file.path(table_dir, "quick_test")
}
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

selected_file <- file.path(
  result_dir, "complex_endpoint_selected_designs.csv"
)
type1_file <- file.path(
  result_dir, "complex_endpoint_type1_by_prevalence.csv"
)

required_files <- c(selected_file, type1_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "The following required result files were not found:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\nRun the complex-endpoint simulation program first."
  )
}

selected <- read.csv(selected_file, check.names = FALSE)
type1 <- read.csv(type1_file, check.names = FALSE)

required_selected_columns <- c(
  "Endpoint", "Method",
  "lambda_A", "gamma_A",
  "lambda_positive", "gamma_positive"
)
required_type1_columns <- c(
  "Endpoint", "Method", "pi",
  "Monte_Carlo_replicates",
  "PRN_any", "PRN_all", "PRN_positive"
)

if (!all(required_selected_columns %in% names(selected))) {
  stop(
    "The selected-design file is missing columns: ",
    paste(
      setdiff(required_selected_columns, names(selected)),
      collapse = ", "
    )
  )
}
if (!all(required_type1_columns %in% names(type1))) {
  stop(
    "The type I error file is missing columns: ",
    paste(
      setdiff(required_type1_columns, names(type1)),
      collapse = ", "
    )
  )
}

N_A <- 40L
N_positive <- 40L
allcomer_looks <- c(20L, 30L, 40L)
positive_looks <- c(10:30, 40L)
dirichlet_prior <- rep(0.25, 4L)

endpoint_order <- c(
  "Nested efficacy",
  "Co-primary efficacy",
  "Efficacy and toxicity"
)

## -----------------------------------------------------------------------------
## 2. Boundary calculation
## -----------------------------------------------------------------------------

lower_tail_boundary <- function(
    n, phi, prior_event, prior_other, cutoff
) {
  x <- 0:n
  posterior_probability <- pbeta(
    phi,
    shape1 = prior_event + x,
    shape2 = prior_other + n - x
  )
  eligible <- x[posterior_probability > cutoff]
  if (length(eligible) == 0L) {
    return(-1L)
  }
  max(eligible)
}

upper_tail_boundary <- function(
    n, phi, prior_event, prior_other, cutoff
) {
  x <- 0:n
  posterior_probability <- pbeta(
    phi,
    shape1 = prior_event + x,
    shape2 = prior_other + n - x,
    lower.tail = FALSE
  )
  eligible <- x[posterior_probability > cutoff]
  if (length(eligible) == 0L) {
    return(n + 1L)
  }
  min(eligible)
}

calculate_endpoint_boundary <- function(
    endpoint, n, lambda, gamma, maximum_n
) {
  cutoff <- 1 - lambda * (n / maximum_n)^gamma
  a <- dirichlet_prior

  if (endpoint == "Nested efficacy") {
    cr_boundary <- lower_tail_boundary(
      n = n,
      phi = 0.15,
      prior_event = a[1L],
      prior_other = sum(a[-1L]),
      cutoff = cutoff
    )
    crpr_boundary <- lower_tail_boundary(
      n = n,
      phi = 0.30,
      prior_event = a[1L] + a[2L],
      prior_other = a[3L] + a[4L],
      cutoff = cutoff
    )

    # The rule uses AND. If either component cannot be crossed, the combined
    # stopping rule cannot be crossed.
    if (cr_boundary < 0L || crpr_boundary < 0L) {
      return(c(boundary1 = NA_integer_, boundary2 = NA_integer_))
    }
    return(c(
      boundary1 = cr_boundary,
      boundary2 = crpr_boundary
    ))
  }

  if (endpoint == "Co-primary efficacy") {
    or_boundary <- lower_tail_boundary(
      n = n,
      phi = 0.10,
      prior_event = a[1L] + a[2L],
      prior_other = a[3L] + a[4L],
      cutoff = cutoff
    )
    efs_boundary <- lower_tail_boundary(
      n = n,
      phi = 0.20,
      prior_event = a[1L] + a[3L],
      prior_other = a[2L] + a[4L],
      cutoff = cutoff
    )

    # The rule uses AND.
    if (or_boundary < 0L || efs_boundary < 0L) {
      return(c(boundary1 = NA_integer_, boundary2 = NA_integer_))
    }
    return(c(
      boundary1 = or_boundary,
      boundary2 = efs_boundary
    ))
  }

  if (endpoint == "Efficacy and toxicity") {
    response_boundary <- lower_tail_boundary(
      n = n,
      phi = 0.45,
      prior_event = a[1L] + a[2L],
      prior_other = a[3L] + a[4L],
      cutoff = cutoff
    )
    toxicity_boundary <- upper_tail_boundary(
      n = n,
      phi = 0.30,
      prior_event = a[1L] + a[3L],
      prior_other = a[2L] + a[4L],
      cutoff = cutoff
    )

    if (response_boundary < 0L) {
      response_boundary <- NA_integer_
    }
    if (toxicity_boundary > n) {
      toxicity_boundary <- NA_integer_
    }
    return(c(
      boundary1 = response_boundary,
      boundary2 = toxicity_boundary
    ))
  }

  stop("Unknown endpoint: ", endpoint)
}

boundary_key <- function(boundary1, boundary2) {
  paste(
    ifelse(is.na(boundary1), "NA", boundary1),
    ifelse(is.na(boundary2), "NA", boundary2),
    sep = "|"
  )
}

compact_consecutive_rows <- function(x) {
  if (nrow(x) == 0L) {
    return(x)
  }

  x <- x[order(x$n), , drop = FALSE]
  groups <- integer(nrow(x))
  groups[1L] <- 1L

  if (nrow(x) > 1L) {
    for (i in 2:nrow(x)) {
      same_boundary <- x$key[i] == x$key[i - 1L]
      consecutive <- x$n[i] == x$n[i - 1L] + 1L
      groups[i] <- groups[i - 1L] +
        as.integer(!(same_boundary && consecutive))
    }
  }

  pieces <- split(x, groups)
  out <- lapply(pieces, function(d) {
    data.frame(
      Endpoint = d$Endpoint[1L],
      Population = d$Population[1L],
      n_start = min(d$n),
      n_end = max(d$n),
      boundary1 = d$boundary1[1L],
      boundary2 = d$boundary2[1L],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

format_patient_range <- function(n_start, n_end) {
  if (n_start == n_end) {
    return(as.character(n_start))
  }
  paste0(n_start, "--", n_end)
}

format_stop_rule <- function(endpoint, boundary1, boundary2) {
  if (endpoint == "Nested efficacy") {
    if (is.na(boundary1) || is.na(boundary2)) {
      return("No stopping boundary")
    }
    return(paste0(
      "\\# CR $\\le ", boundary1,
      "$ and \\# CR/PR $\\le ", boundary2, "$"
    ))
  }

  if (endpoint == "Co-primary efficacy") {
    if (is.na(boundary1) || is.na(boundary2)) {
      return("No stopping boundary")
    }
    return(paste0(
      "\\# ORR $\\le ", boundary1,
      "$ and \\# EFS6 $\\le ", boundary2, "$"
    ))
  }

  if (endpoint == "Efficacy and toxicity") {
    parts <- character()
    if (!is.na(boundary1)) {
      parts <- c(
        parts,
        paste0("\\# responses $\\le ", boundary1, "$")
      )
    }
    if (!is.na(boundary2)) {
      parts <- c(
        parts,
        paste0("\\# toxicities $\\ge ", boundary2, "$")
      )
    }
    if (length(parts) == 0L) {
      return("No stopping boundary")
    }
    return(paste(parts, collapse = " or "))
  }

  stop("Unknown endpoint: ", endpoint)
}

proposed_selected <- selected[
  selected$Method == "Globally calibrated",
  ,
  drop = FALSE
]

missing_endpoints <- setdiff(
  endpoint_order, proposed_selected$Endpoint
)
if (length(missing_endpoints) > 0L) {
  stop(
    "No globally calibrated selected design was found for: ",
    paste(missing_endpoints, collapse = ", ")
  )
}

boundary_rows <- list()
boundary_index <- 1L

for (endpoint in endpoint_order) {
  design <- proposed_selected[
    proposed_selected$Endpoint == endpoint,
    ,
    drop = FALSE
  ]
  if (nrow(design) != 1L) {
    stop(
      "Expected exactly one globally calibrated design for ",
      endpoint, "."
    )
  }

  all_rows <- lapply(allcomer_looks, function(n) {
    boundary <- calculate_endpoint_boundary(
      endpoint = endpoint,
      n = n,
      lambda = design$lambda_A,
      gamma = design$gamma_A,
      maximum_n = N_A
    )
    data.frame(
      Endpoint = endpoint,
      Population = "All-comer",
      n = n,
      boundary1 = unname(boundary[1L]),
      boundary2 = unname(boundary[2L]),
      stringsAsFactors = FALSE
    )
  })
  all_rows <- do.call(rbind, all_rows)
  all_rows$key <- boundary_key(
    all_rows$boundary1, all_rows$boundary2
  )

  positive_rows <- lapply(positive_looks, function(n) {
    boundary <- calculate_endpoint_boundary(
      endpoint = endpoint,
      n = n,
      lambda = design$lambda_positive,
      gamma = design$gamma_positive,
      maximum_n = N_positive
    )
    data.frame(
      Endpoint = endpoint,
      Population = "Biomarker-positive",
      n = n,
      boundary1 = unname(boundary[1L]),
      boundary2 = unname(boundary[2L]),
      stringsAsFactors = FALSE
    )
  })
  positive_rows <- do.call(rbind, positive_rows)
  positive_rows$key <- boundary_key(
    positive_rows$boundary1, positive_rows$boundary2
  )

  boundary_rows[[boundary_index]] <-
    compact_consecutive_rows(all_rows)
  boundary_index <- boundary_index + 1L
  boundary_rows[[boundary_index]] <-
    compact_consecutive_rows(positive_rows)
  boundary_index <- boundary_index + 1L
}

compact_boundaries <- do.call(rbind, boundary_rows)
rownames(compact_boundaries) <- NULL

## -----------------------------------------------------------------------------
## 3. Compact boundary LaTeX table
## -----------------------------------------------------------------------------

boundary_tex_rows <- character()
previous_endpoint <- ""
previous_population <- ""

for (i in seq_len(nrow(compact_boundaries))) {
  d <- compact_boundaries[i, , drop = FALSE]

  endpoint_cell <- if (d$Endpoint != previous_endpoint) {
    d$Endpoint
  } else {
    ""
  }

  population_cell <- if (
    d$Endpoint != previous_endpoint ||
      d$Population != previous_population
  ) {
    d$Population
  } else {
    ""
  }

  boundary_tex_rows <- c(
    boundary_tex_rows,
    paste0(
      endpoint_cell, " & ",
      population_cell, " & ",
      format_patient_range(d$n_start, d$n_end), " & ",
      format_stop_rule(
        d$Endpoint, d$boundary1, d$boundary2
      ),
      " \\\\"
    )
  )

  next_is_new_endpoint <- if (i < nrow(compact_boundaries)) {
    compact_boundaries$Endpoint[i + 1L] != d$Endpoint
  } else {
    FALSE
  }
  if (next_is_new_endpoint) {
    boundary_tex_rows <- c(
      boundary_tex_rows, "\\addlinespace"
    )
  }

  previous_endpoint <- d$Endpoint
  previous_population <- d$Population
}

boundary_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  paste0(
    "\\caption{Compact presentation of the decision boundaries ",
    "for the globally calibrated design with complex categorical ",
    "endpoints.}"
  ),
  "\\label{tab:complex_endpoint_boundaries_compact}",
  "\\small",
  "\\begin{tabularx}{\\textwidth}{llcX}",
  "\\toprule",
  paste0(
    "Endpoint & Population & Number evaluated & ",
    "Futility or no-claim rule \\\\"
  ),
  "\\midrule",
  boundary_tex_rows,
  "\\bottomrule",
  "\\end{tabularx}",
  "\\begin{minipage}{0.98\\textwidth}",
  "\\footnotesize",
  paste0(
    "The all-comer boundaries apply at the interim analyses of 20 and 30 ",
    "patients and at the final analysis of 40 patients. Biomarker-positive ",
    "boundaries apply at the first biomarker-positive assessment and, when ",
    "applicable, at subsequent analyses of 20 and 30 cumulative ",
    "biomarker-positive patients and at the final analysis of 40 ",
    "biomarker-positive patients. At an interim analysis, satisfying the ",
    "listed rule results in stopping for futility or unacceptability. At the ",
    "final analysis, satisfying the listed rule results in no efficacy claim."
  ),
  "\\end{minipage}",
  "\\end{table}"
)

boundary_output <- file.path(
  table_dir, "table_complex_endpoint_boundaries_compact.tex"
)
writeLines(boundary_tex, boundary_output, useBytes = TRUE)

## -----------------------------------------------------------------------------
## 4. Type I error LaTeX table
## -----------------------------------------------------------------------------

type1$Endpoint <- factor(
  type1$Endpoint,
  levels = endpoint_order
)
type1 <- type1[
  order(type1$Endpoint, type1$pi, type1$Method),
  ,
  drop = FALSE
]
type1$Endpoint <- as.character(type1$Endpoint)

proposed_label <- "Globally calibrated"
componentwise_label <- "Componentwise calibrated"

if (!all(c(proposed_label, componentwise_label) %in% type1$Method)) {
  stop(
    "Both globally calibrated and componentwise calibrated ",
    "simulation results are required."
  )
}

fmt_pct <- function(x) {
  formatC(100 * x, format = "f", digits = 1L)
}
fmt_pi <- function(x) {
  formatC(x, format = "f", digits = 1L)
}

type1_tex_rows <- character()

for (endpoint in endpoint_order) {
  d_endpoint <- type1[
    type1$Endpoint == endpoint,
    ,
    drop = FALSE
  ]
  pi_values <- sort(unique(d_endpoint$pi))

  for (j in seq_along(pi_values)) {
    pi_value <- pi_values[j]

    proposed <- d_endpoint[
      d_endpoint$Method == proposed_label &
        abs(d_endpoint$pi - pi_value) < 1e-12,
      ,
      drop = FALSE
    ]
    componentwise <- d_endpoint[
      d_endpoint$Method == componentwise_label &
        abs(d_endpoint$pi - pi_value) < 1e-12,
      ,
      drop = FALSE
    ]

    if (nrow(proposed) != 1L || nrow(componentwise) != 1L) {
      stop(
        "Incomplete type I error result for endpoint = ",
        endpoint, ", pi = ", pi_value
      )
    }

    endpoint_cell <- if (j == 1L) endpoint else ""

    type1_tex_rows <- c(
      type1_tex_rows,
      paste0(
        endpoint_cell, " & ", fmt_pi(pi_value),
        " & \\textbf{", fmt_pct(proposed$PRN_any), "}",
        " & ", fmt_pct(proposed$PRN_all),
        " & ", fmt_pct(proposed$PRN_positive),
        " & \\textbf{", fmt_pct(componentwise$PRN_any), "}",
        " & ", fmt_pct(componentwise$PRN_all),
        " & ", fmt_pct(componentwise$PRN_positive),
        " \\\\"
      )
    )
  }

  if (endpoint != tail(endpoint_order, 1L)) {
    type1_tex_rows <- c(
      type1_tex_rows, "\\addlinespace"
    )
  }
}

replicate_values <- unique(type1$Monte_Carlo_replicates)
replicate_note <- if (length(replicate_values) == 1L) {
  paste0(
    format(replicate_values, big.mark = ","),
    " independent validation replicates per endpoint, method, ",
    "and prevalence"
  )
} else {
  "the independent validation replicates reported in the result file"
}

type1_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  paste0(
    "\\caption{Simulation-based global and component-specific type I ",
    "error rates for complex categorical endpoints.}"
  ),
  "\\label{tab:complex_endpoint_type1}",
  "\\small",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lcrrrrrr}",
  "\\toprule",
  paste0(
    "Endpoint & $\\pi$ & \\multicolumn{3}{c}{Proposed} & ",
    "\\multicolumn{3}{c}{Componentwise} \\\\"
  ),
  paste0(
    " & & \\textbf{PRN-any} & PRN-all & PRN-positive & ",
    "\\textbf{PRN-any} & PRN-all & PRN-positive \\\\"
  ),
  "\\midrule",
  type1_tex_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\begin{minipage}{0.98\\textwidth}",
  "\\footnotesize",
  paste0(
    "Values are percentages based on ", replicate_note, ". ",
    "Boldface indicates PRN-any, the primary global type I error ",
    "measure. The proposed design is calibrated using the probability of ",
    "an efficacy claim in either population, whereas the componentwise ",
    "design calibrates the all-comer and biomarker-positive claim ",
    "probabilities separately."
  ),
  "\\end{minipage}",
  "\\end{table}"
)

type1_output <- file.path(
  table_dir, "table_complex_endpoint_type1.tex"
)
writeLines(type1_tex, type1_output, useBytes = TRUE)

## -----------------------------------------------------------------------------
## 5. Input file and audit CSV
## -----------------------------------------------------------------------------

input_output <- file.path(
  table_dir, "supplemental_complex_endpoint_tables_input.tex"
)
writeLines(
  c(
    "\\input{table_complex_endpoint_boundaries_compact.tex}",
    "\\input{table_complex_endpoint_type1.tex}"
  ),
  input_output,
  useBytes = TRUE
)

write.csv(
  compact_boundaries,
  file.path(
    result_dir, "complex_endpoint_boundaries_compact.csv"
  ),
  row.names = FALSE
)

output_paths <- file.path(
  table_dir,
  c(
    "table_complex_endpoint_boundaries_compact.tex",
    "table_complex_endpoint_type1.tex",
    "supplemental_complex_endpoint_tables_input.tex"
  )
)

message("Created:")
for (output_path in output_paths) {
  if (file.exists(output_path)) {
    message("  ", normalizePath(output_path, mustWork = FALSE))
  } else {
    warning("Expected output file was not created: ", output_path)
  }
}
