# Create the five LaTeX tables reported in the main manuscript.

options(stringsAsFactors = FALSE, scipen = 999)

repo_root <- normalizePath(
  Sys.getenv("ENRICHMENT_BOP2_ROOT", unset = getwd()),
  mustWork = FALSE
)
result_dir <- file.path(repo_root, "results", "main")
table_dir <- file.path(repo_root, "tables", "main")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) stop("Missing result file: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tex <- function(lines, name) {
  writeLines(lines, file.path(table_dir, name), useBytes = TRUE)
}

fmt_pct <- function(x) formatC(100 * x, format = "f", digits = 1L)
fmt_bold_pct <- function(x) paste0("\\textbf{", fmt_pct(x), "}")
fmt_value <- function(x, digits = 1L) {
  y <- formatC(x, format = "f", digits = digits)
  if (digits > 0L) y <- sub("0+$", "", sub("\\.$", "", y))
  y
}

get_row <- function(df, method, pi, theta_pos = NULL, theta_neg = NULL) {
  keep <- df$Method == method & abs(df$pi - pi) < 1e-12
  if (!is.null(theta_pos)) keep <- keep & abs(df$theta_pos - theta_pos) < 1e-12
  if (!is.null(theta_neg)) keep <- keep & abs(df$theta_neg - theta_neg) < 1e-12
  out <- df[keep, , drop = FALSE]
  if (nrow(out) != 1L) stop("Expected one matching row; found ", nrow(out))
  out
}

make_table <- function(caption, label, spec, header, body, note, resize = FALSE) {
  tabular <- c(
    paste0("\\begin{tabular}{", spec, "}"),
    "\\hline", header, "\\hline", body, "\\hline", "\\end{tabular}"
  )
  if (resize) tabular <- c("\\resizebox{\\textwidth}{!}{%", tabular, "}")
  c(
    "\\begin{table}[htbp]", "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), "\\small", tabular,
    "\\begin{minipage}{0.98\\textwidth}", "\\footnotesize", note,
    "\\end{minipage}", "\\end{table}"
  )
}

exact <- read_csv("exact_operating_characteristics.csv")
boundary_all <- read_csv("proposed_boundary_table_allcomer.csv")
boundary_pos <- read_csv("proposed_boundary_table_positive.csv")

extract_boundary <- function(df) {
  cols <- names(df)[grepl("^[0-9]+$", names(df))]
  sizes <- as.integer(cols)
  ord <- order(sizes)
  list(sizes = sizes[ord], boundaries = as.integer(df[1L, cols[ord], drop = TRUE]))
}

compress_boundaries <- function(sizes, boundaries) {
  group <- integer(length(sizes)); group[1L] <- 1L
  if (length(sizes) > 1L) {
    for (i in 2:length(sizes)) {
      same <- boundaries[i] == boundaries[i - 1L]
      consecutive <- sizes[i] == sizes[i - 1L] + 1L
      group[i] <- group[i - 1L] + as.integer(!(same && consecutive))
    }
  }
  do.call(rbind, lapply(split(seq_along(sizes), group), function(idx) {
    lo <- min(sizes[idx]); hi <- max(sizes[idx]); b <- boundaries[idx[1L]]
    data.frame(
      Sample = if (lo == hi) as.character(lo) else paste0(lo, "--", hi),
      Rule = if (b < 0L) "No futility stopping" else paste0("Responses $\\le ", b, "$"),
      stringsAsFactors = FALSE
    )
  }))
}

all_b <- extract_boundary(boundary_all)
pos_b <- extract_boundary(boundary_pos)
compact_all <- compress_boundaries(all_b$sizes, all_b$boundaries)
compact_all$Population <- "All-comer"
compact_pos <- compress_boundaries(pos_b$sizes, pos_b$boundaries)
compact_pos$Population <- "Biomarker-positive"
compact <- rbind(compact_all, compact_pos)

boundary_body <- vapply(seq_len(nrow(compact)), function(i) {
  paste0(compact$Population[i], " & ", compact$Sample[i], " & ", compact$Rule[i], " \\\\")
}, character(1L))

write_tex(
  make_table(
    "Compact presentation of the futility boundaries for the proposed design.",
    "tab:proposed_boundaries_compact", "lcl",
    "Population & Number of patients evaluated & Decision boundary \\\\",
    boundary_body,
    paste0(
      "For the biomarker-positive population, the boundaries at 20 and 30 patients ",
      "also apply to post-enrichment interim looks when applicable. At the final ",
      "look, crossing the listed boundary results in no efficacy claim."
    )
  ),
  "table_proposed_boundaries_compact.tex"
)

null_df <- exact[exact$scenario_type == "Null", , drop = FALSE]
type1_body <- vapply(sort(unique(null_df$pi)), function(pi) {
  proposed <- get_row(null_df, "Globally calibrated", pi)
  comparator <- get_row(null_df, "Componentwise calibrated", pi)
  paste(
    fmt_value(pi),
    fmt_bold_pct(proposed$PRN_any), fmt_pct(proposed$PRN_all), fmt_pct(proposed$PRN_positive),
    fmt_bold_pct(comparator$PRN_any), fmt_pct(comparator$PRN_all), fmt_pct(comparator$PRN_positive),
    sep = " & "
  ) |> paste0(" \\\\")
}, character(1L))

write_tex(
  make_table(
    "Global and component-specific type I error rates under the global null hypothesis (exact enumeration).",
    "tab:type1_exact", "crrrrrr",
    c(
      "$\\pi$ & \\multicolumn{3}{c}{Proposed} & \\multicolumn{3}{c}{Componentwise} \\\\",
      " & \\textbf{PRN-any} & PRN-all & PRN-positive & \\textbf{PRN-any} & PRN-all & PRN-positive \\\\") ,
    type1_body,
    c(
      "Values are percentages. Boldface indicates PRN-any, the primary global type I error measure, defined as the probability of making an efficacy claim in either population.",
      "PRN-all and PRN-positive are the corresponding path-specific false-positive probabilities."
    ),
    resize = TRUE
  ),
  "table_type1_exact.tex"
)

alt <- exact[exact$scenario_type == "Alternative", , drop = FALSE]
for (theta_pos in sort(unique(alt$theta_pos))) {
  sub <- alt[abs(alt$theta_pos - theta_pos) < 1e-12, , drop = FALSE]
  key <- unique(sub[, c("theta_neg", "pi")])
  key <- key[order(-key$theta_neg, key$pi), , drop = FALSE]
  body <- vapply(seq_len(nrow(key)), function(i) {
    proposed <- get_row(sub, "Globally calibrated", key$pi[i], theta_pos, key$theta_neg[i])
    comparator <- get_row(sub, "Componentwise calibrated", key$pi[i], theta_pos, key$theta_neg[i])
    paste(
      fmt_value(key$theta_neg[i]), fmt_value(key$pi[i]),
      fmt_bold_pct(proposed$PRN_any), fmt_pct(proposed$PRN_all), fmt_pct(proposed$PRN_positive),
      fmt_bold_pct(comparator$PRN_any), fmt_pct(comparator$PRN_all), fmt_pct(comparator$PRN_positive),
      sep = " & "
    ) |> paste0(" \\\\")
  }, character(1L))
  tag <- gsub("\\.", "p", format(theta_pos, nsmall = 1L))
  write_tex(
    make_table(
      paste0("Exact efficacy-claim probabilities for $\\theta_+=", fmt_value(theta_pos), "$.") ,
      paste0("tab:power_comparison_exact_", tag), "ccrrrrrr",
      c(
        "$\\theta_-$ & $\\pi$ & \\multicolumn{3}{c}{Proposed} & \\multicolumn{3}{c}{Componentwise} \\\\",
        " & & \\textbf{PRN-any} & PRN-all & PRN-positive & \\textbf{PRN-any} & PRN-all & PRN-positive \\\\") ,
      body,
      c(
        "Values are percentages. Boldface indicates PRN-any, the primary overall power measure, defined as the probability of making an efficacy claim in either population.",
        "PRN-all and PRN-positive are the corresponding all-comer and biomarker-positive efficacy-claim probabilities. The componentwise comparator does not control the global type I error rate."
      ),
      resize = TRUE
    ),
    paste0("table_power_comparison_exact_theta_pos_", tag, ".tex")
  )
}

write_tex(
  c(
    "\\input{table_proposed_boundaries_compact.tex}",
    "\\input{table_type1_exact.tex}",
    "\\input{table_power_comparison_exact_theta_pos_0p4.tex}",
    "\\input{table_power_comparison_exact_theta_pos_0p5.tex}",
    "\\input{table_power_comparison_exact_theta_pos_0p6.tex}"
  ),
  "main_tables_input.tex"
)

message("Main-manuscript tables written to: ", table_dir)
