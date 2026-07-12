root <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(ENRICHMENT_BOP2_ROOT = root)
for (script in c(
  "01_main_binary_exact.R",
  "02_main_tables.R",
  "03_supp_complex_type1.R",
  "04_supp_tables.R"
)) {
  sys.source(file.path(root, "R", script), envir = new.env(parent = globalenv()))
}
