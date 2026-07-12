root <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(ENRICHMENT_BOP2_ROOT = root)
sys.source(file.path(root, "R", "01_main_binary_exact.R"), envir = new.env(parent = globalenv()))
sys.source(file.path(root, "R", "02_main_tables.R"), envir = new.env(parent = globalenv()))
