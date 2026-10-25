# Repository-relative entry point for reviewer-requested supplementary analyses.
# Usage: Rscript run_revision.R <r1-benchmarks|r1-reverse|r2-claims|r2-adjusted|r2-all|main-timing|case-study>
# These are separate runs. They DO NOT modify manuscript files or overwrite
# the immutable reference_results/ snapshots.
args <- commandArgs(trailingOnly=TRUE)
valid <- c("r1-benchmarks", "r1-reverse", "r2-claims", "r2-adjusted",
           "r2-all", "main-timing", "case-study")
if (length(args) != 1L || !args[1L] %in% valid) {
  cat("Usage: Rscript run_revision.R <",paste(valid, collapse="|"),">\n",sep="")
  quit(status=if (length(args)==0L) 0L else 2L)
}
all_args <- commandArgs(trailingOnly=FALSE)
script_arg <- grep("^--file=",all_args,value=TRUE)
if (length(script_arg) == 0L) {
  stop("Run as 'Rscript run_revision.R ...' from the repository directory.")
}
root <- normalizePath(dirname(sub("^--file=", "", script_arg[1L])), mustWork=TRUE)
rev <- file.path(root, "R", "revision")
res <- file.path(root, "results", "revision")
dir.create(res,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(ENRICHMENT_BOP2_ROOT=root, BOP2_REVISION_RESULTS_ROOT=res)
mode <- args[1L]
if (mode=="r1-benchmarks") {
  # Default gamma maximum 8. Set BOP2_PROPOSED_GAMMA_MAX=12 for the
  # independently explored extended-grid sensitivity analysis.
  Sys.setenv(BOP2_R1M1_OUTPUT_DIR=file.path(res,paste0("R1M1_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())))
  setwd(root)
  source(file.path(rev,"R1M1_additional_settings_exact.R"))
} else if (mode=="r1-reverse") {
  Sys.setenv(BOP2_R1M2_OUTPUT_DIR=file.path(res,paste0("R1M2_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())))
  setwd(root)
  source(file.path(rev,"R1M2_reverse_order_exact.R"))
} else if (mode %in% c("r2-claims", "r2-adjusted", "r2-all")) {
  setwd(rev)
  filename <- switch(mode,
    "r2-claims"="R2_Comments2_3_subgroup_only_exact.R",
    "r2-adjusted"="R2_Comment4_error_matched_BOP2.R",
    "r2-all"="R2_run_all.R")
  source(filename)
} else if (mode=="main-timing") {
  Sys.setenv(BOP2_MAIN_TIMING_OUTPUT_DIR=file.path(res,"main_timing"))
  setwd(root)
  source(file.path(rev,"main_binary_boundary_calibration_timing.R"))
} else if (mode=="case-study") {
  Sys.setenv(BOP2_CASE_STUDY_OUTPUT_DIR=file.path(res,"case_study"))
  setwd(root)
  source(file.path(rev,"case_study_boundary_calibration_timing.R"))
}
