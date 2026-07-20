# Create LaTeX tables from exact enrichment BOP2 operating characteristics
#
# Run through run_main.R from the repository root. The analysis results are
# read from results/main and LaTeX tables are written to tables/main.

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## -----------------------------------------------------------------------------
## 1. Paths
## -----------------------------------------------------------------------------

project_dir <- path.expand(Sys.getenv(
  "ENRICHMENT_BOP2_ROOT",
  unset = normalizePath(getwd(), mustWork = TRUE)
))
quick_test <- identical(Sys.getenv("BOP2_QUICK"), "1")
result_dir <- file.path(project_dir, "results", "main")
table_dir <- file.path(project_dir, "tables", "main")
if (quick_test) {
  result_dir <- file.path(result_dir, "quick_test")
  table_dir <- file.path(table_dir, "quick_test")
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Remove legacy Monte Carlo and exact-versus-Monte-Carlo tables so stale files
# are not included in the exact-only output manifest.
legacy_tex_patterns <- c(
  "_mc_",
  "monte_carlo",
  "exact_vs_mc"
)
legacy_tex <- list.files(table_dir, pattern = "\\.tex$", full.names = TRUE)
if (length(legacy_tex) > 0L) {
  legacy_name <- basename(legacy_tex)
  remove_legacy <- Reduce(
    `|`,
    lapply(legacy_tex_patterns, function(pattern) grepl(pattern, legacy_name))
  )
  unlink(legacy_tex[remove_legacy], force = TRUE)
}

if (!dir.exists(result_dir)) {
  stop("Result directory does not exist: ", result_dir)
}

message("Reading results from: ", result_dir)
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

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

fmt_prob <- function(x, digits = 1L) {
  ifelse(is.na(x), "--", formatC(100 * x, format = "f", digits = digits))
}

fmt_prob_bold <- function(x, digits = 1L) {
  value <- fmt_prob(x, digits = digits)
  ifelse(value == "--", value, paste0("\\textbf{", value, "}"))
}

fmt_num <- function(x, digits = 1L) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

fmt_param <- function(x, digits = 3L) {
  ifelse(is.na(x), "--", sub("0+$", "", sub("\\.$", "", formatC(
    x, format = "f", digits = digits
  ))))
}

method_display <- function(x) {
  out <- x
  out[out == "Globally calibrated"] <- "Proposed"
  out[out == "Componentwise calibrated"] <- "Componentwise"
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
    stop("Expected one row but found ", nrow(out),
         " for method = ", method)
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
    "\\hline",
    header_lines,
    "\\hline",
    body_lines,
    "\\hline",
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
## 3. Read exact result files
## -----------------------------------------------------------------------------

settings <- read_required_csv("analysis_settings.csv")
selected <- read_required_csv("selected_designs.csv")
exact <- read_required_csv("exact_operating_characteristics.csv")
boundary_all <- read_required_csv("proposed_boundary_table_allcomer.csv")
boundary_pos <- read_required_csv("proposed_boundary_table_positive.csv")

required_methods <- c("Globally calibrated", "Componentwise calibrated")
if (!all(required_methods %in% unique(exact$Method))) {
  stop("The expected methods were not found in the exact results.")
}

required_post_enrichment_columns <- c(
  "First_positive_assessment",
  "Bridge_to_minimum",
  "PET_positive_initial",
  "Enter_enrichment",
  "Analyzed_positive_look20",
  "PET_positive_look20",
  "Analyzed_positive_look30",
  "PET_positive_look30"
)
if (!all(required_post_enrichment_columns %in% names(exact))) {
  stop(
    "The exact results do not contain the post-enrichment monitoring columns. ",
    "Run the exact analysis program first."
  )
}

## 4. Selected design table
## -----------------------------------------------------------------------------

selected$Method_display <- method_display(selected$Method)
selected <- selected[match(required_methods, selected$Method), , drop = FALSE]

selected_body <- vapply(seq_len(nrow(selected)), function(i) {
  paste(
    selected$Method_display[i],
    fmt_param(selected$lambda_A[i], 3),
    fmt_param(selected$gamma_A[i], 3),
    fmt_param(selected$lambda_positive[i], 3),
    fmt_param(selected$gamma_positive[i], 3),
    fmt_prob(selected$max_global_type1[i]),
    fmt_prob(selected$average_power[i]),
    fmt_prob(selected$minimum_power[i]),
    sep = " & "
  ) |> paste0(" \\\\")
}, character(1L))

selected_tex <- make_table_environment(
  caption = paste0(
    "Selected tuning parameters and exact operating characteristics used for ",
    "design calibration."
  ),
  label = "tab:selected_designs",
  column_spec = "lrrrrrrr",
  header_lines = paste0(
    "Method & $\\lambda_A$ & $\\gamma_A$ & $\\lambda_+$ & ",
    "$\\gamma_+$ & Max. global type I error (\\%) & ",
    "Average calibration power (\\%) & Minimum calibration power (\\%) \\\\"
  ),
  body_lines = selected_body,
  note_lines = paste0(
    "Average and minimum calibration power are calculated over the ",
    "prespecified biomarker-prevalence values under the single alternative ",
    "configuration $(\\theta_+,\\theta_-)=(p_1,p_0)$. The componentwise ",
    "comparator was calibrated separately for the all-comer and ",
    "biomarker-positive efficacy-claim components and was not constrained to ",
    "control their union."
  ),
  resize = TRUE
)
write_tex(selected_tex, "table_selected_designs.tex")

## -----------------------------------------------------------------------------
## 5. Global type I error tables
## -----------------------------------------------------------------------------

make_type1_table <- function(df, source_name, filename, label_suffix) {
  null_df <- df[df$scenario_type == "Null", , drop = FALSE]
  pi_values <- sort(unique(null_df$pi))

  body <- vapply(pi_values, function(pi_value) {
    prop <- get_method_row(null_df, "Globally calibrated", pi_value)
    comp <- get_method_row(null_df, "Componentwise calibrated", pi_value)
    paste(
      fmt_param(pi_value, 1),
      fmt_prob_bold(prop$PRN_any),
      fmt_prob(prop$PRN_all),
      fmt_prob(prop$PRN_positive),
      fmt_prob_bold(comp$PRN_any),
      fmt_prob(comp$PRN_all),
      fmt_prob(comp$PRN_positive),
      sep = " & "
    ) |> paste0(" \\\\")
  }, character(1L))

  header <- c(
    paste0(
      "$\\pi$ & \\multicolumn{3}{c}{Proposed} & ",
      "\\multicolumn{3}{c}{Componentwise} \\\\"
    ),
    paste0(
      " & \\textbf{PRN-any} & PRN-all & PRN-positive & ",
      "\\textbf{PRN-any} & PRN-all & PRN-positive \\\\"
    )
  )

  tex <- make_table_environment(
    caption = paste0(
      "Global and component-specific type I error rates under the global null ",
      "hypothesis (", source_name, ")."
    ),
    label = paste0("tab:type1_", label_suffix),
    column_spec = "crrrrrr",
    header_lines = header,
    body_lines = body,
    note_lines = c(
      "Values are percentages. Boldface indicates PRN-any, the primary global type I error measure, defined as the probability of making an efficacy claim in either population.",
      "PRN-all and PRN-positive are the corresponding path-specific false-positive probabilities."
    ),
    resize = TRUE
  )
  write_tex(tex, filename)
}

make_type1_table(
  exact,
  source_name = "exact enumeration",
  filename = "table_type1_exact.tex",
  label_suffix = "exact"
)

## -----------------------------------------------------------------------------
## 6. Power comparison and proposed-design operating-characteristic tables
## -----------------------------------------------------------------------------

make_power_comparison_tables <- function(df, source_name, source_tag) {
  alt <- df[df$scenario_type == "Alternative", , drop = FALSE]
  theta_pos_values <- sort(unique(alt$theta_pos))

  generated <- character()

  for (theta_pos_value in theta_pos_values) {
    sub <- alt[abs(alt$theta_pos - theta_pos_value) < 1e-12, , drop = FALSE]
    key <- unique(sub[, c("pi", "theta_neg")])
    key <- key[order(-key$theta_neg, key$pi), , drop = FALSE]

    body_comparison <- vapply(seq_len(nrow(key)), function(i) {
      pi_value <- key$pi[i]
      theta_neg_value <- key$theta_neg[i]
      prop <- get_method_row(
        sub, "Globally calibrated", pi_value,
        theta_pos_value, theta_neg_value
      )
      comp <- get_method_row(
        sub, "Componentwise calibrated", pi_value,
        theta_pos_value, theta_neg_value
      )
      paste(
        fmt_param(theta_neg_value, 1),
        fmt_param(pi_value, 1),
        fmt_prob_bold(prop$PRN_any),
        fmt_prob(prop$PRN_all),
        fmt_prob(prop$PRN_positive),
        fmt_prob_bold(comp$PRN_any),
        fmt_prob(comp$PRN_all),
        fmt_prob(comp$PRN_positive),
        sep = " & "
      ) |> paste0(" \\\\")
    }, character(1L))

    theta_tag <- gsub("\\.", "p", format(theta_pos_value, nsmall = 1))
    comparison_filename <- paste0(
      "table_power_comparison_", source_tag,
      "_theta_pos_", theta_tag, ".tex"
    )

    comparison_tex <- make_table_environment(
      caption = paste0(
        "Exact efficacy-claim probabilities for $\\theta_+=",
        fmt_param(theta_pos_value, 1), "$."
      ),
      label = paste0(
        "tab:power_comparison_", source_tag, "_", theta_tag
      ),
      column_spec = "ccrrrrrr",
      header_lines = c(
        paste0(
          "$\\theta_-$ & $\\pi$ & \\multicolumn{3}{c}{Proposed} & ",
          "\\multicolumn{3}{c}{Componentwise} \\\\"
        ),
        paste0(
          " & & \\textbf{PRN-any} & PRN-all & PRN-positive & ",
          "\\textbf{PRN-any} & PRN-all & PRN-positive \\\\"
        )
      ),
      body_lines = body_comparison,
      note_lines = c(
        "Values are percentages. Boldface indicates PRN-any, the primary overall power measure, defined as the probability of making an efficacy claim in either population.",
        "PRN-all and PRN-positive are the corresponding all-comer and biomarker-positive efficacy-claim probabilities. The componentwise comparator does not control the global type I error rate."
      ),
      resize = TRUE
    )
    write_tex(comparison_tex, comparison_filename)
    generated <- c(generated, comparison_filename)

    prop_sub <- sub[sub$Method == "Globally calibrated", , drop = FALSE]
    prop_sub <- prop_sub[order(-prop_sub$theta_neg, prop_sub$pi), , drop = FALSE]

    body_oc <- vapply(seq_len(nrow(prop_sub)), function(i) {
      r <- prop_sub[i, ]
      paste(
        fmt_param(r$theta_neg, 1),
        fmt_param(r$pi, 1),
        fmt_prob(r$PRN_any),
        fmt_prob(r$PRN_all),
        fmt_prob(r$PRN_positive),
        fmt_prob(r$PET_any),
        fmt_prob(r$PET_positive),
        fmt_prob(r$Enter_enrichment),
        fmt_num(r$Mean_total_sample_size, 1),
        fmt_num(r$Mean_positive_sample_size, 1),
        sep = " & "
      ) |> paste0(" \\\\")
    }, character(1L))

    oc_filename <- paste0(
      "table_proposed_oc_", source_tag,
      "_theta_pos_", theta_tag, ".tex"
    )

    oc_tex <- make_table_environment(
      caption = paste0(
        "Operating characteristics of the proposed design for $\\theta_+=",
        fmt_param(theta_pos_value, 1), "$ (", source_name, ")."
      ),
      label = paste0("tab:proposed_oc_", source_tag, "_", theta_tag),
      column_spec = "ccrrrrrrrr",
      header_lines = paste0(
        "$\\theta_-$ & $\\pi$ & PRN-any & PRN-all & PRN-positive & ",
        "PET-any & PET-positive & Enrichment & Mean $N$ & Mean $N_+$ \\\\"
      ),
      body_lines = body_oc,
      note_lines = paste0(
        "PRN and PET quantities and the enrichment probability are reported as ",
        "percentages. Because an all-comer futility decision transitions to the ",
        "biomarker-positive-only branch rather than terminating the trial, PET-any ",
        "equals PET-positive under this design. Mean $N$ is the mean total sample ",
        "size, and mean $N_+$ is the mean number of biomarker-positive patients ",
        "evaluated."
      ),
      resize = TRUE
    )
    write_tex(oc_tex, oc_filename)
    generated <- c(generated, oc_filename)
  }

  generated
}

power_files <- make_power_comparison_tables(
  exact,
  source_name = "exact enumeration",
  source_tag = "exact"
)

## -----------------------------------------------------------------------------
## 7. Biomarker-positive enrichment-path tables
## -----------------------------------------------------------------------------

make_enrichment_path_tables <- function(df, source_name, source_tag) {
  alt <- df[df$scenario_type == "Alternative" &
              df$Method == "Globally calibrated", , drop = FALSE]
  theta_pos_values <- sort(unique(alt$theta_pos))
  generated <- character()

  for (theta_pos_value in theta_pos_values) {
    sub <- alt[abs(alt$theta_pos - theta_pos_value) < 1e-12, , drop = FALSE]
    sub <- sub[order(-sub$theta_neg, sub$pi), , drop = FALSE]

    tolerance <- 1e-9
    if (any(abs(
      sub$First_positive_assessment -
        sub$PET_positive_initial - sub$Enter_enrichment
    ) > tolerance)) {
      stop("First assessment does not equal initial futility plus entry.")
    }
    if (any(sub$Bridge_to_minimum - sub$First_positive_assessment > tolerance) ||
        any(sub$PET_positive_look20 - sub$Analyzed_positive_look20 > tolerance) ||
        any(sub$PET_positive_look30 - sub$Analyzed_positive_look30 > tolerance)) {
      stop("Inconsistent detailed biomarker-positive path probabilities.")
    }

    body <- vapply(seq_len(nrow(sub)), function(i) {
      r <- sub[i, ]
      paste(
        fmt_param(r$theta_neg, 1),
        fmt_param(r$pi, 1),
        fmt_prob(r$First_positive_assessment),
        fmt_prob(r$Bridge_to_minimum),
        fmt_prob(r$PET_positive_initial),
        fmt_prob(r$Enter_enrichment),
        fmt_prob(r$Analyzed_positive_look20),
        fmt_prob(r$PET_positive_look20),
        fmt_prob(r$Analyzed_positive_look30),
        fmt_prob(r$PET_positive_look30),
        sep = " & "
      ) |> paste0(" \\\\")
    }, character(1L))

    theta_tag <- gsub("\\.", "p", format(theta_pos_value, nsmall = 1))
    filename <- paste0(
      "table_proposed_enrichment_path_", source_tag,
      "_theta_pos_", theta_tag, ".tex"
    )

    tex <- make_table_environment(
      caption = paste0(
        "Biomarker-positive enrichment-path operating characteristics for ",
        "$\\theta_+=", fmt_param(theta_pos_value, 1), "$ (", source_name, ")."
      ),
      label = paste0(
        "tab:proposed_enrichment_path_", source_tag, "_", theta_tag
      ),
      column_spec = "ccrrrrrrrr",
      header_lines = paste0(
        "$\\theta_-$ & $\\pi$ & First assessment & Bridged & ",
        "Initial futility & Entered & Analyzed at 20 & Stop at 20 & ",
        "Analyzed at 30 & Stop at 30 \\\\"
      ),
      body_lines = body,
      note_lines = paste0(
        "All quantities are unconditional percentages among all initiated trials. ",
        "First assessment is the probability that the trial reaches the first ",
        "biomarker-positive assessment after crossing an all-comer futility ",
        "boundary. Bridged is the probability that additional biomarker-positive ",
        "enrollment is required to reach $n_{+,\\min}$ before that assessment. ",
        "Initial futility is the probability of stopping at the first ",
        "biomarker-positive assessment, whether performed immediately or after ",
        "bridging. Entered is the probability of passing that assessment and ",
        "continuing on the enrichment path. Analyzed at 20 and Analyzed at 30 ",
        "indicate that the corresponding biomarker-positive analysis is conducted, ",
        "including when the first assessment itself occurs at that sample size. ",
        "Stop at 20 and Stop at 30 refer to futility at a subsequent scheduled ",
        "analysis and exclude stopping at the first assessment. An analysis ",
        "already passed at the first assessment is skipped. Values are rounded ",
        "to one decimal place; therefore, a displayed value of 0.0 does not ",
        "necessarily indicate an exact probability of zero."
      ),
      resize = TRUE
    )
    write_tex(tex, filename)
    generated <- c(generated, filename)
  }

  generated
}

enrichment_files <- make_enrichment_path_tables(
  exact,
  source_name = "exact enumeration",
  source_tag = "exact"
)

## -----------------------------------------------------------------------------
## 8. Proposed-design boundary tables
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

boundary_display <- function(x) {
  ifelse(x < 0L, "--", as.character(x))
}

all_b <- extract_boundary(boundary_all)
pos_b <- extract_boundary(boundary_pos)

make_boundary_block <- function(title, sample_sizes, boundaries) {
  c(
    paste0(
      "\\multicolumn{", length(sample_sizes) + 1L,
      "}{l}{\\textit{", title, "}} \\\\"
    ),
    paste0(
      "Number of patients evaluated & ",
      paste(sample_sizes, collapse = " & "), " \\\\"
    ),
    paste0(
      "Futile/no efficacy claim if responses $\\le$ & ",
      paste(boundary_display(boundaries), collapse = " & "), " \\\\"
    )
  )
}

# Split the biomarker-positive table so it remains readable on a portrait page.
pos_split_1 <- which(pos_b$sample_sizes <= 20L)
pos_split_2 <- which(pos_b$sample_sizes > 20L)

boundary_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Decision boundaries of the proposed globally calibrated enrichment BOP2 design.}",
  "\\label{tab:proposed_boundaries}",
  "\\small",
  paste0("\\begin{tabular}{l", paste(rep("c", length(all_b$sample_sizes)), collapse = ""), "}"),
  "\\hline",
  make_boundary_block(
    "All-comer population",
    all_b$sample_sizes,
    all_b$boundaries
  ),
  "\\hline",
  "\\end{tabular}",
  "\\vspace{0.7em}",
  "",
  paste0("\\resizebox{\\textwidth}{!}{%"),
  paste0("\\begin{tabular}{l", paste(rep("c", length(pos_split_1)), collapse = ""), "}"),
  "\\hline",
  make_boundary_block(
    "Biomarker-positive population",
    pos_b$sample_sizes[pos_split_1],
    pos_b$boundaries[pos_split_1]
  ),
  "\\hline",
  "\\end{tabular}",
  "}",
  "\\vspace{0.7em}",
  "",
  paste0("\\begin{tabular}{l", paste(rep("c", length(pos_split_2)), collapse = ""), "}"),
  "\\hline",
  make_boundary_block(
    "Biomarker-positive population (continued)",
    pos_b$sample_sizes[pos_split_2],
    pos_b$boundaries[pos_split_2]
  ),
  "\\hline",
  "\\end{tabular}",
  "\\begin{minipage}{0.98\\textwidth}",
  "\\footnotesize",
  paste0(
    "A dash indicates that no response count can satisfy the futility rule at ",
    "that sample size. Biomarker-positive boundaries apply at the first ",
    "biomarker-positive assessment and, when applicable, at subsequent analyses ",
    "of 20 and 30 cumulative biomarker-positive patients and at the final ",
    "analysis of 40 biomarker-positive patients. At an interim analysis, ",
    "satisfying the listed rule results in stopping for futility. At the final ",
    "analysis, satisfying the listed rule results in no efficacy claim."
  ),
  "\\end{minipage}",
  "\\end{table}"
)
write_tex(boundary_lines, "table_proposed_boundaries.tex")

compress_boundaries <- function(sample_sizes, boundaries) {
  stopifnot(length(sample_sizes) == length(boundaries))
  if (length(sample_sizes) == 0L) {
    return(data.frame())
  }

  groups <- integer(length(sample_sizes))
  groups[1L] <- 1L
  group_id <- 1L
  if (length(sample_sizes) > 1L) {
    for (i in 2:length(sample_sizes)) {
      same_boundary <- boundaries[i] == boundaries[i - 1L]
      consecutive <- sample_sizes[i] == sample_sizes[i - 1L] + 1L
      if (!(same_boundary && consecutive)) {
        group_id <- group_id + 1L
      }
      groups[i] <- group_id
    }
  }

  out <- lapply(split(seq_along(sample_sizes), groups), function(idx) {
    lo <- min(sample_sizes[idx])
    hi <- max(sample_sizes[idx])
    sample_label <- if (lo == hi) as.character(lo) else paste0(lo, "--", hi)
    boundary <- boundaries[idx[1L]]
    rule <- if (boundary < 0L) {
      "No futility stopping"
    } else {
      paste0("\\# responses $\\le ", boundary, "$")
    }
    data.frame(
      Sample_size = sample_label,
      Boundary = boundary,
      Rule = rule,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

compact_all <- compress_boundaries(all_b$sample_sizes, all_b$boundaries)
compact_all$Population <- "All-comer"
compact_pos <- compress_boundaries(pos_b$sample_sizes, pos_b$boundaries)
compact_pos$Population <- "Biomarker-positive"
compact <- rbind(compact_all, compact_pos)

compact_body <- vapply(seq_len(nrow(compact)), function(i) {
  paste(
    compact$Population[i],
    compact$Sample_size[i],
    compact$Rule[i],
    sep = " & "
  ) |> paste0(" \\\\")
}, character(1L))

compact_tex <- make_table_environment(
  caption = paste0(
    "Compact presentation of the decision boundaries for the proposed design."
  ),
  label = "tab:proposed_boundaries_compact",
  column_spec = "lcl",
  header_lines = paste0(
    "Population & Number of patients evaluated & Futility or no-claim rule \\\\"
  ),
  body_lines = compact_body,
  note_lines = paste0(
    paste0(
      "For the biomarker-positive population, the boundaries apply at the first ",
      "biomarker-positive assessment and at subsequent analyses of 20 and 30 ",
      "cumulative biomarker-positive patients when applicable. At an interim ",
      "analysis, satisfying the listed ",
      "rule results in stopping for futility. At the final analysis, satisfying ",
      "the listed rule results in no efficacy claim."
    )
  ),
  resize = FALSE
)
write_tex(compact_tex, "table_proposed_boundaries_compact.tex")

## -----------------------------------------------------------------------------
## 9. Manifest and optional input file
## -----------------------------------------------------------------------------

tex_files <- sort(list.files(table_dir, pattern = "\\.tex$", full.names = FALSE))
tex_files <- setdiff(tex_files, "all_tables_input.tex")

manifest_lines <- c(
  "Enrichment BOP2 exact TeX table files",
  "=====================================",
  "",
  "Primary result files read from results/main:",
  "  - exact_operating_characteristics.csv: exact operating characteristics",
  "  - exact_type1_error.csv: exact type I error results",
  "  - exact_power.csv: exact power results",
  "  - selected_designs.csv: selected tuning parameters and calibration summary",
  "  - proposed_boundary_table_allcomer.csv: proposed all-comer boundaries",
  "  - proposed_boundary_table_positive.csv: proposed biomarker-positive boundaries",
  "",
  "Generated TeX files:",
  paste0("  - ", tex_files),
  "",
  "All reported operating characteristics are obtained by exact recursive",
  "enumeration; no Monte Carlo tables are generated. The wide tables use",
  "\\\\resizebox, so include graphicx in the manuscript preamble."
)
writeLines(manifest_lines, file.path(table_dir, "tables_manifest.txt"))

input_lines <- c(
  "% Automatically generated input list. Remove tables not needed in the manuscript.",
  paste0("\\input{", tex_files, "}")
)
write_tex(input_lines, "all_tables_input.tex")

main_table_files <- c(
  "table_proposed_boundaries_compact.tex",
  "table_type1_exact.tex",
  "table_power_comparison_exact_theta_pos_0p4.tex",
  "table_power_comparison_exact_theta_pos_0p5.tex",
  "table_power_comparison_exact_theta_pos_0p6.tex"
)
main_input_lines <- c(
  "% Recommended tables for the main manuscript.",
  paste0("\\input{", main_table_files, "}")
)
write_tex(main_input_lines, "main_tables_input.tex")

supplemental_table_files <- c(
  "table_proposed_enrichment_path_exact_theta_pos_0p4.tex",
  "table_proposed_enrichment_path_exact_theta_pos_0p5.tex",
  "table_proposed_enrichment_path_exact_theta_pos_0p6.tex"
)
supplemental_input_lines <- c(
  "% Recommended enrichment-path tables for the Supplementary Material.",
  paste0("\\input{", supplemental_table_files, "}")
)
write_tex(supplemental_input_lines, "supplemental_tables_input.tex")

message("Completed. Exact-result TeX tables are available in: ", table_dir)
