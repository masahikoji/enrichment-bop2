# Reviewer 2 Comments 2 and 3: additional subgroup-only scenarios and claim targeting.
# Place R2_exact_common.R in the SAME directory, then run:
#   setwd("/path/to/unzipped/folder")
#   source("R2_Comments2_3_subgroup_only_exact.R")
# Base R. Exact finite-state enumeration, NOT Monte Carlo simulation.

source("R2_exact_common.R")
R2_verify_original()
outdir <- R2_new_output_dir("R2_Comments2_3_exact")
scenarios <- R2_scenarios()
designs <- list(R2_proposed_design, R2_independent_design)

# The same positive kernel is reused for different prevalence values.
kernels <- new.env(parent = emptyenv())
get_kernel <- function(d, tp) {
  key <- paste0(d$name, "_", format(tp, nsmall=5, digits=9))
  if (!exists(key, envir=kernels, inherits=FALSE)) {
    assign(key, R2_positive_kernel(d, tp), envir=kernels)
  }
  get(key, envir=kernels, inherits=FALSE)
}

results <- vector("list", nrow(scenarios)*length(designs))
idx <- 0L
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  st20 <- R2_joint_states(20L, s$pi, s$theta_positive, s$theta_negative)
  st30 <- R2_joint_states(30L, s$pi, s$theta_positive, s$theta_negative)
  for (d in designs) {
    idx <- idx+1L
    a <- R2_evaluate_scenario(s$pi, s$theta_positive, s$theta_negative, d,
                              kernel=get_kernel(d,s$theta_positive), st20=st20, st30=st30)
    a$scenario <- s$scenario
    results[[idx]] <- a
  }
}
res <- do.call(rbind,results)
rownames(res) <- NULL
stopifnot(nrow(res) == 2L*nrow(scenarios))

subgroup_only <- res[res$scenario == "Subgroup_only", ,drop=FALSE]
stopifnot(all(subgroup_only$null_allcomer_true), all(!subgroup_only$null_positive_true))

percent <- res
prob_columns <- c("PRN_any", "PRN_all", "PRN_positive", "P_false_allcomer", "P_false_positive",
                  "pointwise_FWER", "P_valid_allcomer", "P_valid_positive", "P_at_least_one_valid_claim",
                  "P_correct_positive_target", "P_potentially_overbroad_allcomer")
for (nm in prob_columns) percent[[paste0(nm,"_pct")]] <- round(100*percent[[nm]],4)
write.csv(res, file.path(outdir,"R2_Comments2_3_all_exact.csv"),row.names=FALSE)
write.csv(subgroup_only, file.path(outdir,"R2_C2_subgroup_only_exact.csv"),row.names=FALSE)
write.csv(percent[percent$scenario=="Subgroup_only", ,drop=FALSE],
          file.path(outdir,"R2_C2_subgroup_only_percent.csv"),row.names=FALSE)
write.csv(percent[percent$scenario=="Original_30_alternatives", ,drop=FALSE],
          file.path(outdir,"R2_C3_original_30_claim_targeting_percent.csv"),row.names=FALSE)
write.csv(percent[percent$scenario=="Reverse_order_R1M2", ,drop=FALSE],
          file.path(outdir,"R2_C3_reverse_order_claim_targeting_percent.csv"),row.names=FALSE)

# Terminology:
#  - P_false_allcomer = P(R_A) only when theta_A<=p0 (statistical false claim).
#  - P_potentially_overbroad_allcomer = P(R_A) if theta_pos>p0 and theta_neg<=p0.
#    This is a CLINICAL overbreadth descriptor, not necessarily a type I error:
#    under the 30 original alternatives theta_A>p0, even though theta_neg<=p0.
message("=== Subgroup-only: exact PRN outcomes (%) ===")
print(percent[percent$scenario=="Subgroup_only",
              c("pi","theta_positive","theta_negative","theta_all","method",
                "PRN_all_pct","PRN_positive_pct","PRN_any_pct","pointwise_FWER_pct")],row.names=FALSE)
message("=== Commentary and interpretation ===")
message("P_correct_positive_target is relevant when positive efficacy is expected but negative efficacy is not.")
message("P_potentially_overbroad_allcomer is NOT the same as false all-comer claim if theta_A>p0.")
message("Tables S9 and/or a response table should be built from the CSVs after checking the output.")
message("Results directory: ",outdir)
