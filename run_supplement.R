root <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(ENRICHMENT_BOP2_ROOT = root)
sys.source(file.path(root, "R", "03_supp_complex_type1.R"), envir = new.env(parent = globalenv()))
sys.source(file.path(root, "R", "04_supp_tables.R"), envir = new.env(parent = globalenv()))
