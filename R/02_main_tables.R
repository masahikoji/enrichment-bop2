# Create LaTeX tables from the corrected exact binary-endpoint results.
#
# This version is for the NO-BRIDGING design and the corrected comparator:
#   * a conventional all-comer BOP2 calibrated independently;
#   * an entry-specific library of conventional biomarker-positive BOP2 designs;
#   * the two standalone components are embedded in the adaptive enrichment
#     procedure without constraining PRN-any.
#
# Repository inputs:
#   results/main
#
# Repository outputs:
#   tables/main

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## -----------------------------------------------------------------------------
## 1. Paths
## -----------------------------------------------------------------------------

repo_root <- normalizePath(
  Sys.getenv("ENRICHMENT_BOP2_ROOT", unset = getwd()),
  mustWork = FALSE
)
result_dir <- file.path(repo_root, "results", "main")
table_dir <- file.path(repo_root, "tables", "main")
if (identical(Sys.getenv("BOP2_QUICK"), "1")) {
  result_dir <- file.path(result_dir, "quick_test")
  table_dir <- file.path(table_dir, "quick_test")
}
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(result_dir)) {
  stop("Result directory does not exist: ", result_dir)
}

message("Reading corrected exact results from: ", result_dir)
message("Writing TeX tables to: ", table_dir)

## -----------------------------------------------------------------------------
## 2. Utility functions
## -----------------------------------------------------------------------------

read_required_csv <- function(filename) {
  path <- file.path(result_dir, filename)
  if (!file.exists(path)) {
    stop("Required result file is missing: ", path)
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tex <- function(lines, filename) {
  path <- file.path(table_dir, filename)
  writeLines(lines, con = path, useBytes = TRUE)
  message("  wrote ", filename)
  invisible(path)
}

fmt_prob <- function(x, digits = 1L) {
  ifelse(is.na(x), "--", formatC(100 * x, format = "f", digits = digits))
}

fmt_prob_bold <- function(x, digits = 1L) {
  value <- fmt_prob(x, digits)
  ifelse(value == "--", value, paste0("\\textbf{", value, "}"))
}

fmt_num <- function(x, digits = 1L) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

fmt_param <- function(x, digits = 3L) {
  ifelse(
    is.na(x),
    "--",
    sub("0+$", "", sub("\\.$", "", formatC(x, format = "f", digits = digits)))
  )
}

fmt_integer_or_dash <- function(x) {
  ifelse(is.na(x), "--", as.character(as.integer(x)))
}

method_display <- function(x) {
  out <- x
  out[out == "Globally calibrated"] <- "Proposed"
  out[out == "Independent conventional BOP2 components"] <-
    "Independent BOP2 components"
  out
}

get_method_row <- function(df, method, pi_value = NULL,
                           theta_pos_value = NULL,
                           theta_neg_value = NULL) {
  keep <- df$Method == method
  if (!is.null(pi_value)) {
    keep <- keep & abs(df$pi - pi_value) < 1e-12
  }
  if (!is.null(theta_pos_value)) {
    keep <- keep & abs(df$theta_pos - theta_pos_value) < 1e-12
  }
  if (!is.null(theta_neg_value)) {
    keep <- keep & abs(df$theta_neg - theta_neg_value) < 1e-12
  }
  out <- df[keep, , drop = FALSE]
  if (nrow(out) != 1L) {
    stop(
      "Expected one row but found ", nrow(out),
      " for method = ", method
    )
  }
  out
}

make_table_environment <- function(caption, label, column_spec,
                                   header_lines, body_lines,
                                   note_lines = character(),
                                   resize = FALSE,
                                   font_command = "\\small") {
  tabular_lines <- c(
    paste0("\\begin{tabular}{", column_spec, "}"),
    "\\toprule",
    header_lines,
    "\\midrule",
    body_lines,
    "\\bottomrule",
    "\\end{tabular}"
  )

  if (resize) {
    tabular_lines <- c(
      "\\resizebox{\\textwidth}{!}{%",
      tabular_lines,
      "}"
    )
  }

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    font_command,
    tabular_lines
  )

  if (length(note_lines) > 0L) {
    lines <- c(
      lines,
      "\\begin{minipage}{0.98\\textwidth}",
      "\\footnotesize",
      note_lines,
      "\\end{minipage}"
    )
  }

  c(lines, "\\end{table}")
}

## -----------------------------------------------------------------------------
## 3. Read corrected exact results
## -----------------------------------------------------------------------------

settings <- read_required_csv("analysis_settings.csv")
selected <- read_required_csv("selected_designs.csv")
exact <- read_required_csv("exact_operating_characteristics.csv")
boundary_all <- read_required_csv("proposed_boundary_table_allcomer.csv")
boundary_pos <- read_required_csv("proposed_boundary_table_positive.csv")
comparator_all <- read_required_csv("comparator_allcomer_bop2_boundary.csv")
comparator_pos <- read_required_csv("comparator_positive_bop2_library.csv")

proposed_method <- "Globally calibrated"
comparator_method <- "Independent conventional BOP2 components"
required_methods <- c(proposed_method, comparator_method)

if (!all(required_methods %in% unique(exact$Method))) {
  stop(
    "Expected methods were not found. Required: ",
    paste(required_methods, collapse = "; ")
  )
}

required_exact_columns <- c(
  "scenario_type", "pi", "theta_pos", "theta_neg",
  "used_for_calibration", "Method",
  "PRN_any", "PRN_all", "PRN_positive",
  "PET_any", "PET_all", "PET_positive",
  "PET_positive_initial", "PET_positive_look20",
  "PET_positive_look30", "PET_low_biomarker_count",
  "Reassess_positive", "Enter_enrichment",
  "Enter_enrichment_look1", "Enter_enrichment_look2",
  "Reach_positive_look20", "Reach_positive_look30",
  "Mean_total_sample_size", "Mean_positive_sample_size"
)
missing_exact <- setdiff(required_exact_columns, names(exact))
if (length(missing_exact) > 0L) {
  stop(
    "The exact result file is missing columns: ",
    paste(missing_exact, collapse = ", ")
  )
}

required_pos_library_columns <- c(
  "entry_n", "analysis_schedule", "lambda", "gamma",
  "standalone_type1", "standalone_power",
  "boundary_at_entry", "boundary_20", "boundary_30", "boundary_40"
)
missing_pos_library <- setdiff(required_pos_library_columns, names(comparator_pos))
if (length(missing_pos_library) > 0L) {
  stop(
    "The positive BOP2 library is missing columns: ",
    paste(missing_pos_library, collapse = ", ")
  )
}

# Basic consistency checks for the no-bridging path decomposition.
tol <- 1e-8
if (max(abs(exact$PRN_any - exact$PRN_all - exact$PRN_positive)) > tol) {
  stop("PRN-any does not equal PRN-all + PRN-positive.")
}
if (max(abs(
  exact$PET_all - exact$PET_low_biomarker_count - exact$Reassess_positive
)) > tol) {
  stop(
    "All-comer futility does not partition into low biomarker count and ",
    "immediate positive reassessment."
  )
}
if (max(abs(
  exact$Reassess_positive - exact$PET_positive_initial -
    exact$Enter_enrichment
)) > tol) {
  stop(
    "Positive reassessment does not partition into initial futility and ",
    "entry into enrichment."
  )
}

## -----------------------------------------------------------------------------
## 4. Calibration-summary table
## -----------------------------------------------------------------------------

selected <- selected[match(required_methods, selected$Method), , drop = FALSE]
selected$Method_display <- method_display(selected$Method)

positive_specification <- ifelse(
  selected$Method == proposed_method,
  paste0(
    "$(", fmt_param(selected$lambda_positive, 3), ",", 
    fmt_param(selected$gamma_positive, 3), ")$"
  ),
  "Entry-specific BOP2 library"
)
allcomer_specification <- paste0(
  "$(", fmt_param(selected$lambda_A, 3), ",",
  fmt_param(selected$gamma_A, 3), ")$"
)

selected_body <- vapply(seq_len(nrow(selected)), function(i) {
  paste0(
    selected$Method_display[i], " & ",
    allcomer_specification[i], " & ",
    positive_specification[i], " & ",
    fmt_prob(selected$max_global_type1[i]), " & ",
    fmt_prob(selected$max_allcomer_type1[i]), " & ",
    fmt_prob(selected$max_positive_type1[i]), " & ",
    fmt_prob(selected$average_power[i]), " & ",
    fmt_prob(selected$minimum_power[i]), " \\\\"
  )
}, character(1L))

selected_tex <- make_table_environment(
  caption = paste0(
    "Selected proposed design and independently calibrated conventional ",
    "BOP2 comparator."
  ),
  label = "tab:selected_designs",
  column_spec = "lllrrrrr",
  header_lines = paste0(
    "Method & All-comer $(\\lambda,\\gamma)$ & Positive component & ",
    "Max. PRN-any (\\%) & Max. PRN-all (\\%) & ",
    "Max. PRN-positive (\\%) & Average power (\\%) & ",
    "Minimum power (\\%) \\\\"
  ),
  body_lines = selected_body,
  note_lines = paste0(
    "For the proposed design, average and minimum power are calculated over ",
    "the prespecified prevalence values under ",
    "$(\\theta_+,\\theta_-)=(p_1,p_0)$. The comparator all-comer BOP2 and ",
    "each entry-specific positive-population BOP2 were calibrated as ",
    "standalone conventional BOP2 designs before being embedded in the ",
    "adaptive procedure. No constraint was imposed on their union."
  ),
  resize = TRUE
)
write_tex(selected_tex, "table_selected_designs.tex")

## -----------------------------------------------------------------------------
## 5. Global type I error table
## -----------------------------------------------------------------------------

null_df <- exact[exact$scenario_type == "Null", , drop = FALSE]
pi_values <- sort(unique(null_df$pi))

type1_body <- vapply(pi_values, function(pi_value) {
  prop <- get_method_row(null_df, proposed_method, pi_value)
  comp <- get_method_row(null_df, comparator_method, pi_value)
  paste0(
    fmt_param(pi_value, 1), " & ",
    fmt_prob_bold(prop$PRN_any), " & ",
    fmt_prob(prop$PRN_all), " & ",
    fmt_prob(prop$PRN_positive), " & ",
    fmt_prob_bold(comp$PRN_any), " & ",
    fmt_prob(comp$PRN_all), " & ",
    fmt_prob(comp$PRN_positive), " \\\\"
  )
}, character(1L))

type1_tex <- make_table_environment(
  caption = paste0(
    "Exact global and path-specific type I error rates under the global null ",
    "hypothesis."
  ),
  label = "tab:type1_exact",
  column_spec = "crrrrrr",
  header_lines = c(
    paste0(
      "$\\pi$ & \\multicolumn{3}{c}{Proposed} & ",
      "\\multicolumn{3}{c}{Independent BOP2 components} \\\\"
    ),
    paste0(
      " & \\textbf{PRN-any} & PRN-all & PRN-positive & ",
      "\\textbf{PRN-any} & PRN-all & PRN-positive \\\\"
    )
  ),
  body_lines = type1_body,
  note_lines = c(
    "Values are percentages obtained by exact recursive enumeration. Boldface indicates PRN-any, the global probability of an efficacy claim in either population.",
    "The independent comparator calibrates each conventional BOP2 component under its standalone analysis schedule and does not control the union probability."
  ),
  resize = TRUE
)
write_tex(type1_tex, "table_type1_exact.tex")

## -----------------------------------------------------------------------------
## 6. Power-comparison and proposed-design OC tables
## -----------------------------------------------------------------------------

alt <- exact[exact$scenario_type == "Alternative", , drop = FALSE]
theta_pos_values <- sort(unique(alt$theta_pos))
power_files <- character()

for (theta_pos_value in theta_pos_values) {
  sub <- alt[abs(alt$theta_pos - theta_pos_value) < 1e-12, , drop = FALSE]
  key <- unique(sub[, c("pi", "theta_neg")])
  key <- key[order(-key$theta_neg, key$pi), , drop = FALSE]

  comparison_body <- vapply(seq_len(nrow(key)), function(i) {
    pi_value <- key$pi[i]
    theta_neg_value <- key$theta_neg[i]
    prop <- get_method_row(
      sub, proposed_method, pi_value, theta_pos_value, theta_neg_value
    )
    comp <- get_method_row(
      sub, comparator_method, pi_value, theta_pos_value, theta_neg_value
    )
    paste0(
      fmt_param(theta_neg_value, 1), " & ",
      fmt_param(pi_value, 1), " & ",
      fmt_prob_bold(prop$PRN_any), " & ",
      fmt_prob(prop$PRN_all), " & ",
      fmt_prob(prop$PRN_positive), " & ",
      fmt_prob_bold(comp$PRN_any), " & ",
      fmt_prob(comp$PRN_all), " & ",
      fmt_prob(comp$PRN_positive), " \\\\"
    )
  }, character(1L))

  theta_tag <- gsub("\\.", "p", format(theta_pos_value, nsmall = 1))
  comparison_filename <- paste0(
    "table_power_comparison_exact_theta_pos_", theta_tag, ".tex"
  )

  comparison_tex <- make_table_environment(
    caption = paste0(
      "Exact efficacy-claim probabilities for $\\theta_+=",
      fmt_param(theta_pos_value, 1), "$."
    ),
    label = paste0("tab:power_comparison_exact_", theta_tag),
    column_spec = "ccrrrrrr",
    header_lines = c(
      paste0(
        "$\\theta_-$ & $\\pi$ & \\multicolumn{3}{c}{Proposed} & ",
        "\\multicolumn{3}{c}{Independent BOP2 components} \\\\"
      ),
      paste0(
        " & & \\textbf{PRN-any} & PRN-all & PRN-positive & ",
        "\\textbf{PRN-any} & PRN-all & PRN-positive \\\\"
      )
    ),
    body_lines = comparison_body,
    note_lines = c(
      "Values are percentages. Boldface indicates PRN-any, the probability of an efficacy claim in either population.",
      "The comparator values are obtained after its independently calibrated standalone BOP2 components are embedded in the no-bridging adaptive enrichment procedure."
    ),
    resize = TRUE
  )
  write_tex(comparison_tex, comparison_filename)
  power_files <- c(power_files, comparison_filename)

  prop_sub <- sub[sub$Method == proposed_method, , drop = FALSE]
  prop_sub <- prop_sub[order(-prop_sub$theta_neg, prop_sub$pi), , drop = FALSE]

  oc_body <- vapply(seq_len(nrow(prop_sub)), function(i) {
    r <- prop_sub[i, , drop = FALSE]
    paste0(
      fmt_param(r$theta_neg, 1), " & ",
      fmt_param(r$pi, 1), " & ",
      fmt_prob(r$PRN_any), " & ",
      fmt_prob(r$PRN_all), " & ",
      fmt_prob(r$PRN_positive), " & ",
      fmt_prob(r$PET_any), " & ",
      fmt_prob(r$PET_low_biomarker_count), " & ",
      fmt_prob(r$Enter_enrichment), " & ",
      fmt_num(r$Mean_total_sample_size, 1), " & ",
      fmt_num(r$Mean_positive_sample_size, 1), " \\\\"
    )
  }, character(1L))

  oc_filename <- paste0(
    "table_proposed_oc_exact_theta_pos_", theta_tag, ".tex"
  )
  oc_tex <- make_table_environment(
    caption = paste0(
      "Exact operating characteristics of the proposed design for ",
      "$\\theta_+=", fmt_param(theta_pos_value, 1), "$."
    ),
    label = paste0("tab:proposed_oc_exact_", theta_tag),
    column_spec = "ccrrrrrrrr",
    header_lines = paste0(
      "$\\theta_-$ & $\\pi$ & PRN-any & PRN-all & PRN-positive & ",
      "PET-any & Stop: $n_+<10$ & Enrichment & Mean $N$ & Mean $N_+$ \\\\"
    ),
    body_lines = oc_body,
    note_lines = paste0(
      "PRN, PET, low-positive-count stopping, and enrichment are percentages. ",
      "Under the no-bridging rule, an all-comer futility decision with fewer ",
      "than 10 biomarker-positive patients terminates the trial. Mean $N$ and ",
      "mean $N_+$ are the mean total and biomarker-positive sample sizes."
    ),
    resize = TRUE
  )
  write_tex(oc_tex, oc_filename)
  power_files <- c(power_files, oc_filename)
}

## -----------------------------------------------------------------------------
## 7. No-bridging enrichment-path tables
## -----------------------------------------------------------------------------

alt_proposed <- exact[
  exact$scenario_type == "Alternative" & exact$Method == proposed_method,
  , drop = FALSE
]
enrichment_files <- character()

for (theta_pos_value in sort(unique(alt_proposed$theta_pos))) {
  sub <- alt_proposed[
    abs(alt_proposed$theta_pos - theta_pos_value) < 1e-12,
    , drop = FALSE
  ]
  sub <- sub[order(-sub$theta_neg, sub$pi), , drop = FALSE]

  path_body <- vapply(seq_len(nrow(sub)), function(i) {
    r <- sub[i, , drop = FALSE]
    paste0(
      fmt_param(r$theta_neg, 1), " & ",
      fmt_param(r$pi, 1), " & ",
      fmt_prob(r$PET_all), " & ",
      fmt_prob(r$PET_low_biomarker_count), " & ",
      fmt_prob(r$Reassess_positive), " & ",
      fmt_prob(r$PET_positive_initial), " & ",
      fmt_prob(r$Enter_enrichment), " & ",
      fmt_prob(r$Reach_positive_look20), " & ",
      fmt_prob(r$PET_positive_look20), " & ",
      fmt_prob(r$Reach_positive_look30), " & ",
      fmt_prob(r$PET_positive_look30), " \\\\"
    )
  }, character(1L))

  theta_tag <- gsub("\\.", "p", format(theta_pos_value, nsmall = 1))
  filename <- paste0(
    "table_proposed_enrichment_path_exact_theta_pos_", theta_tag, ".tex"
  )

  path_tex <- make_table_environment(
    caption = paste0(
      "No-bridging enrichment-path operating characteristics for ",
      "$\\theta_+=", fmt_param(theta_pos_value, 1), "$."
    ),
    label = paste0("tab:proposed_enrichment_path_exact_", theta_tag),
    column_spec = "ccrrrrrrrrr",
    header_lines = paste0(
      "$\\theta_-$ & $\\pi$ & All-comer futility & Stop: $n_+<10$ & ",
      "Positive reassessment & Initial positive futility & Entered & ",
      "Reach 20 & Stop at 20 & Reach 30 & Stop at 30 \\\\"
    ),
    body_lines = path_body,
    note_lines = paste0(
      "All quantities are unconditional percentages among initiated trials. ",
      "All-comer futility partitions into stopping because $n_+<10$ and ",
      "immediate positive reassessment. Positive reassessment partitions into ",
      "initial positive futility and entry into enrichment. Analyses already ",
      "reached when enrichment begins are not repeated."
    ),
    resize = TRUE
  )
  write_tex(path_tex, filename)
  enrichment_files <- c(enrichment_files, filename)
}

## -----------------------------------------------------------------------------
## 8. Proposed boundary tables
## -----------------------------------------------------------------------------

extract_boundary <- function(df) {
  sample_columns <- names(df)[grepl("^[0-9]+$", names(df))]
  sample_sizes <- as.integer(sample_columns)
  ord <- order(sample_sizes)
  list(
    sample_sizes = sample_sizes[ord],
    boundaries = as.integer(df[1L, sample_columns[ord], drop = TRUE])
  )
}

compress_boundaries <- function(sample_sizes, boundaries) {
  if (length(sample_sizes) == 0L) {
    return(data.frame())
  }
  groups <- integer(length(sample_sizes))
  groups[1L] <- 1L
  if (length(sample_sizes) > 1L) {
    for (i in 2:length(sample_sizes)) {
      same_boundary <- boundaries[i] == boundaries[i - 1L]
      consecutive <- sample_sizes[i] == sample_sizes[i - 1L] + 1L
      groups[i] <- groups[i - 1L] + as.integer(!(same_boundary && consecutive))
    }
  }
  pieces <- split(seq_along(sample_sizes), groups)
  do.call(rbind, lapply(pieces, function(idx) {
    lo <- min(sample_sizes[idx])
    hi <- max(sample_sizes[idx])
    data.frame(
      Sample_size = if (lo == hi) as.character(lo) else paste0(lo, "--", hi),
      Boundary = boundaries[idx[1L]],
      stringsAsFactors = FALSE
    )
  }))
}

all_b <- extract_boundary(boundary_all)
pos_b <- extract_boundary(boundary_pos)
compact_all <- compress_boundaries(all_b$sample_sizes, all_b$boundaries)
compact_all$Population <- "All-comer"
compact_pos <- compress_boundaries(pos_b$sample_sizes, pos_b$boundaries)
compact_pos$Population <- "Biomarker-positive"
compact <- rbind(compact_all, compact_pos)

compact_body <- vapply(seq_len(nrow(compact)), function(i) {
  rule <- if (compact$Boundary[i] < 0L) {
    "No futility stopping"
  } else {
    paste0("Responses $\\le ", compact$Boundary[i], "$")
  }
  paste0(
    compact$Population[i], " & ", compact$Sample_size[i], " & ",
    rule, " \\\\"
  )
}, character(1L))

compact_tex <- make_table_environment(
  caption = "Compact decision boundaries for the proposed design.",
  label = "tab:proposed_boundaries_compact",
  column_spec = "lcl",
  header_lines = paste0(
    "Population & Number evaluated & Futility or no-claim rule \\\\"
  ),
  body_lines = compact_body,
  note_lines = paste0(
    "At an interim analysis, satisfying the listed rule results in futility. ",
    "At the final analysis, satisfying the listed rule results in no efficacy ",
    "claim. Biomarker-positive boundaries apply at the immediate assessment ",
    "and at subsequent scheduled analyses that have not already been passed."
  )
)
write_tex(compact_tex, "table_proposed_boundaries_compact.tex")

## -----------------------------------------------------------------------------
## 9. Standalone comparator calibration tables
## -----------------------------------------------------------------------------

comparator_all <- comparator_all[order(comparator_all$analysis_n), , drop = FALSE]
all_body <- vapply(seq_len(nrow(comparator_all)), function(i) {
  paste0(
    comparator_all$analysis_n[i], " & ",
    comparator_all$futility_or_no_claim_boundary[i], " & ",
    ifelse(comparator_all$is_final[i], "Final", "Interim"), " & ",
    fmt_prob(comparator_all$standalone_type1[i], 2), " & ",
    fmt_prob(comparator_all$standalone_power[i], 2), " \\\\"
  )
}, character(1L))

all_tex <- make_table_environment(
  caption = paste0(
    "Standalone calibration of the conventional all-comer BOP2 component."
  ),
  label = "tab:comparator_allcomer_bop2",
  column_spec = "cccrr",
  header_lines = paste0(
    "Analysis sample size & Boundary & Analysis type & Type I error (\\%) & ",
    "Power (\\%) \\\\"
  ),
  body_lines = all_body,
  note_lines = paste0(
    "The conventional all-comer BOP2 was calibrated independently at ",
    "$p_0=0.20$ and $p_1=0.40$. Among boundary rules with type I error not ",
    "exceeding 0.10, the rule closest to 0.10 was selected, with power used ",
    "to resolve ties."
  )
)
write_tex(all_tex, "table_comparator_allcomer_bop2.tex")

comparator_pos <- comparator_pos[order(comparator_pos$entry_n), , drop = FALSE]
pos_body <- vapply(seq_len(nrow(comparator_pos)), function(i) {
  paste0(
    comparator_pos$entry_n[i], " & ",
    comparator_pos$analysis_schedule[i], " & ",
    fmt_param(comparator_pos$lambda[i], 3), " & ",
    fmt_param(comparator_pos$gamma[i], 3), " & ",
    fmt_prob(comparator_pos$standalone_type1[i], 2), " & ",
    fmt_prob(comparator_pos$standalone_power[i], 2), " & ",
    fmt_integer_or_dash(comparator_pos$boundary_at_entry[i]), " & ",
    fmt_integer_or_dash(comparator_pos$boundary_20[i]), " & ",
    fmt_integer_or_dash(comparator_pos$boundary_30[i]), " & ",
    fmt_integer_or_dash(comparator_pos$boundary_40[i]), " \\\\"
  )
}, character(1L))

pos_tex <- make_table_environment(
  caption = paste0(
    "Entry-specific standalone calibration of the conventional ",
    "biomarker-positive BOP2 components."
  ),
  label = "tab:comparator_positive_bop2_library",
  column_spec = "clrrrrrrrr",
  header_lines = paste0(
    "Entry $m$ & Analysis schedule & $\\lambda$ & $\\gamma$ & ",
    "Type I error (\\%) & Power (\\%) & $b_m$ & $b_{20}$ & ",
    "$b_{30}$ & $b_{40}$ \\\\"
  ),
  body_lines = pos_body,
  note_lines = paste0(
    "For each attainable first positive-population assessment size ",
    "$m=10,\\ldots,30$, a conventional BOP2 was calibrated independently ",
    "using the schedule consisting of $m$, future scheduled analyses among ",
    "20 and 30, and the final analysis at 40. A dash denotes an analysis not ",
    "included in that entry-specific schedule."
  ),
  resize = TRUE,
  font_command = "\\scriptsize"
)
write_tex(pos_tex, "table_comparator_positive_bop2_library.tex")

## -----------------------------------------------------------------------------
## 10. Input lists and manifest
## -----------------------------------------------------------------------------

main_table_files <- c(
  "table_proposed_boundaries_compact.tex",
  "table_type1_exact.tex",
  "table_power_comparison_exact_theta_pos_0p4.tex",
  "table_power_comparison_exact_theta_pos_0p5.tex",
  "table_power_comparison_exact_theta_pos_0p6.tex"
)
write_tex(
  c(
    "% Recommended tables for the main manuscript.",
    paste0("\\input{", main_table_files, "}")
  ),
  "main_tables_input.tex"
)

supplemental_table_files <- c(
  "table_selected_designs.tex",
  "table_comparator_allcomer_bop2.tex",
  "table_comparator_positive_bop2_library.tex",
  enrichment_files
)
write_tex(
  c(
    "% Recommended tables for the Supplementary Material.",
    paste0("\\input{", supplemental_table_files, "}")
  ),
  "supplemental_tables_input.tex"
)

tex_files <- sort(list.files(table_dir, pattern = "\\.tex$", full.names = FALSE))
tex_files <- setdiff(tex_files, "all_tables_input.tex")
write_tex(
  c(
    "% Automatically generated input list. Remove tables not required.",
    paste0("\\input{", tex_files, "}")
  ),
  "all_tables_input.tex"
)

manifest <- c(
  "Corrected enrichment BOP2 TeX table files",
  "=========================================",
  "",
  "Required binary-result files:",
  "  - analysis_settings.csv",
  "  - selected_designs.csv",
  "  - exact_operating_characteristics.csv",
  "  - proposed_boundary_table_allcomer.csv",
  "  - proposed_boundary_table_positive.csv",
  "  - comparator_allcomer_bop2_boundary.csv",
  "  - comparator_positive_bop2_library.csv",
  "",
  "Generated TeX files:",
  paste0("  - ", sort(list.files(table_dir, pattern = "\\.tex$"))),
  "",
  "The wide tables require graphicx; the booktabs commands require booktabs."
)
writeLines(manifest, file.path(table_dir, "tables_manifest.txt"))

message("Completed. Corrected exact-result TeX tables are in: ", table_dir)
