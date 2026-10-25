# Reviewer 2 Comment 4: globally error-matched independent BOP2 comparator.
# Place R2_exact_common.R in the SAME directory, then run:
#   setwd("/path/to/unzipped/folder")
#   source("R2_Comment4_error_matched_BOP2.R")
#
# Base R, exact finite-state enumeration. No Monte Carlo simulation.
#
# This comparator keeps the original independent standalone BOP2 structure:
#   - one standalone all-comer BOP2, and
#   - a separate standalone positive BOP2 for each entry m=10,...,30.
# For each common individual-component nominal cap tau, the conventional BOP2
# selection criterion is applied to each standalone component separately.
# After embedding all 22 rules in the SAME no-bridging enrichment procedure,
# the exact maximum point-global-null FWER over the five prespecified prevalences
# is evaluated. Among the investigated component caps with maximum
# point-global-null error <=0.10, select the cap achieving the largest
# maximum global error. The proposed method does not enter tau selection.
# No assumption of monotonicity in tau is made (discrete boundaries may jump).
# The result is NOT a claim of strong FWER under composite nulls.

source("R2_exact_common.R")
R2_verify_original()
outdir <- R2_new_output_dir("R2_Comment4_error_matched")

# Fine grid from the original independently calibrated BOP2 comparator.
lambda_grid <- seq(0.001, 1.000, by=0.0025)
gamma_grid <- seq(0.010, 10.000, by=0.025)

R2_make_fine_candidates <- function(looks) {
  grid <- expand.grid(lambda=lambda_grid, gamma=gamma_grid, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  b <- matrix(0L, nrow=nrow(grid), ncol=length(looks), dimnames=list(NULL,paste0("b_",looks)))
  for (j in seq_along(looks)) {
    n <- as.integer(looks[j])
    post <- pbeta(R2_p0, shape1=0.5+(0:n), shape2=0.5+n-(0:n))
    cutoffs <- 1-grid$lambda * (n/40)^grid$gamma
    # Counts post_prob > cutoff, minus one. Strict inequality matches submitted R.
    b[,j] <- as.integer(n - findInterval(cutoffs, rev(post)))
  }
  first <- !duplicated(b)
  list(B=b[first, ,drop=FALSE], lambda=grid$lambda[first], gamma=grid$gamma[first])
}

R2_precompute_standalone <- function(B, looks, full_indices=NULL) {
  keep <- !duplicated(B)
  if (is.null(full_indices)) full_indices <- seq_len(nrow(B))
  B <- B[keep, ,drop=FALSE]
  ix <- full_indices[keep]
  list(B=B, looks=looks, full_index=ix,
       type1=R2_standalone_prob_matrix(B, looks, R2_p0),
       power=R2_standalone_prob_matrix(B, looks, R2_p1))
}

R2_choose_standalone <- function(pre, nominal_cap, tolerance=1e-12) {
  eligible <- which(pre$type1 <= nominal_cap+tolerance)
  if (!length(eligible)) stop("No standalone BOP2 candidate at component nominal cap ",nominal_cap)
  best_type1 <- max(pre$type1[eligible])
  choice <- eligible[abs(pre$type1[eligible]-best_type1) <= tolerance]
  best_power <- max(pre$power[choice])
  choice <- choice[abs(pre$power[choice]-best_power) <= tolerance]
  as.integer(choice[1L])
}

message("Generating conventional BOP2 boundary tables on the original fine grid...")
Afull <- R2_make_fine_candidates(c(20L,30L,40L))
Pfull <- R2_make_fine_candidates(c(10:30,40L))
Apre <- R2_precompute_standalone(Afull$B,c(20L,30L,40L))
Ppre <- vector("list",21L)
names(Ppre) <- as.character(10:30)
for (m in 10:30) {
  looks <- c(m,c(20L,30L,40L)[c(20L,30L,40L)>m])
  cols <- match(paste0("b_",looks),colnames(Pfull$B))
  Ppre[[as.character(m)]] <- R2_precompute_standalone(Pfull$B[,cols,drop=FALSE],looks,
                                                        full_indices=seq_len(nrow(Pfull$B)))
}
message("All-comer unique full tables: ", nrow(Afull$B), "; positive unique full tables: ",nrow(Pfull$B))

# Verify exact reproduction of the comparator already reported in the manuscript.
original_A_idx <- R2_choose_standalone(Apre,R2_alpha)
original_A <- setNames(as.integer(Apre$B[original_A_idx,]),c("20","30","40"))
if (!identical(original_A,R2_independent_A)) stop("Original all-comer BOP2 boundaries were not reproduced")
for (m in 10:30) {
  pre <- Ppre[[as.character(m)]]
  ix <- R2_choose_standalone(pre,R2_alpha)
  row <- as.integer(pre$B[ix,])
  orig <- R2_independent_library[m-9L,,drop=FALSE]
  if (row[1L]!=orig$b_entry || row[length(row)]!=orig$b_40) stop("Original positive entry library mismatch at m=",m)
  if (m<20L && row[which(pre$looks==20L)]!=orig$b_20) stop("Original positive b20 mismatch at m=",m)
  if (m<30L && row[which(pre$looks==30L)]!=orig$b_30) stop("Original positive b30 mismatch at m=",m)
}
message("Original fine-grid comparator boundary reproduction: ALL 22 COMPONENTS PASSED")

# Cache every candidate's exact positive-claim probability conditional on each
# entry state (m,x) under the point global null. The library may vary with m.
message("Caching standalone positive-entry success probabilities under p0...")
q_table <- lapply(10:30, function(m) {
  pre <- Ppre[[as.character(m)]]
  t(vapply(seq_len(nrow(pre$B)),function(i) {
    R2_positive_success_at_entry(m,as.integer(pre$B[i,]),R2_p0)
  },numeric(m+1L)))
})
names(q_table) <- as.character(10:30)

R2_selected_design_from_cap <- function(tau) {
  ai <- R2_choose_standalone(Apre,tau)
  ii <- vapply(Ppre, R2_choose_standalone, integer(1L), nominal_cap=tau)
  list(A_index=ai, P_indices=ii)
}

R2_kernel_from_selection <- function(ii) {
  kernel <- numeric(31L*32L/2L)
  for (m in 10:30) {
    j <- m-9L
    ids <- as.integer(m*(m+1L)/2L+(0:m)+1L)
    kernel[ids] <- q_table[[j]][ii[j],]
  }
  kernel
}

null_st20 <- lapply(R2_prevalences, function(pi) R2_joint_states(20L,pi,R2_p0,R2_p0))
null_st30 <- lapply(R2_prevalences, function(pi) R2_joint_states(30L,pi,R2_p0,R2_p0))
all_cache <- new.env(hash=TRUE,parent=emptyenv())
R2_get_cached_all <- function(A_index, pi_index) {
  key <- paste(A_index,pi_index,sep="_")
  if (!exists(key,envir=all_cache,inherits=FALSE)) {
    A <- setNames(as.integer(Apre$B[A_index,]),c("20","30","40"))
    result <- R2_all_components(A,null_st20[[pi_index]],null_st30[[pi_index]],R2_p0)
    assign(key,result,envir=all_cache)
  }
  get(key,envir=all_cache,inherits=FALSE)
}

R2_global_null_by_cap <- function(tau) {
  selected <- R2_selected_design_from_cap(tau)
  kern <- R2_kernel_from_selection(selected$P_indices)
  p <- vapply(seq_along(R2_prevalences),function(j) {
    z <- R2_get_cached_all(selected$A_index,j)
    z$PRN_all + sum(z$weights*kern)
  }, numeric(1L))
  list(tau=tau, selected=selected, point_null=p,
       max_global=max(p), min_global=min(p), mean_global=mean(p))
}

proposed_null <- vapply(R2_prevalences,function(pi) {
  R2_evaluate_scenario(pi,R2_p0,R2_p0,R2_proposed_design)$PRN_any
},numeric(1L))
proposed_max <- max(proposed_null)
if (abs(proposed_max - 0.09894542782786911)>3e-9) stop("Proposed maximum error regression check failed")

# Each standalone component changes only when its null rejection probability
# crosses a candidate's attainable null error. Include all these breakpoints,
# not just a coarse tau grid, and do not use bisection/assume monotonicity.
crit <- unlist(c(list(Apre$type1),lapply(Ppre,function(x)x$type1)),use.names=FALSE)
crit <- crit[crit>=0.015 & crit<=R2_alpha+1e-12]
tau_grid <- sort(unique(c(round(crit,13),seq(0.02,0.10,by=0.001),0.05,0.10)))
tau_grid <- tau_grid[tau_grid>=0.015 & tau_grid<=R2_alpha+1e-12]
message("Checking ",length(tau_grid)," distinct common nominal component caps under all five point-global-null prevalences...")
search <- vector("list",length(tau_grid))
for (i in seq_along(tau_grid)) {
  z <- R2_global_null_by_cap(tau_grid[i])
  search[[i]] <- data.frame(component_cap=tau_grid[i], max_point_null=z$max_global,
                            min_point_null=z$min_global, mean_point_null=z$mean_global,
                            absolute_gap_vs_proposed=abs(z$max_global-proposed_max),
                            valid_nominal_10pct=z$max_global<=R2_alpha+1e-12,
                            bA20=Apre$B[z$selected$A_index,1],
                            bA30=Apre$B[z$selected$A_index,2],
                            bA40=Apre$B[z$selected$A_index,3])
}
search_df <- do.call(rbind,search)
write.csv(search_df,file.path(outdir,"R2_C4_component_cap_search.csv"),row.names=FALSE)
feasible <- which(search_df$valid_nominal_10pct)
if (!length(feasible)) stop("No independently calibrated BOP2 comparator satisfies the global point-null criterion")
best_global_error <- max(search_df$max_point_null[feasible])
closest <- feasible[abs(search_df$max_point_null[feasible] - best_global_error) <= 1e-13]
# If different component caps yield the same attained maximum global error,
# use the largest cap as a prespecified deterministic tie-breaker.
best_id <- closest[which.max(search_df$component_cap[closest])]
tau_best <- search_df$component_cap[best_id]
matched <- R2_global_null_by_cap(tau_best)
bonferroni_style <- R2_global_null_by_cap(0.05) # illustrative component cap; verify global rate explicitly

R2_design_from_selection <- function(sel,label) {
  A <- setNames(as.integer(Apre$B[sel$A_index,]),c("20","30","40"))
  rows <- lapply(10:30,function(m) {
    pre <- Ppre[[as.character(m)]]
    row <- as.integer(pre$B[sel$P_indices[m-9L],])
    looks <- pre$looks
    data.frame(entry_n=m, b_entry=row[1L],
               b_20=if (20L %in% looks) row[match(20L,looks)] else NA_integer_,
               b_30=if (30L %in% looks) row[match(30L,looks)] else NA_integer_,
               b_40=row[length(row)],
               standalone_type1=pre$type1[sel$P_indices[m-9L]],
               standalone_power=pre$power[sel$P_indices[m-9L]],
               stringsAsFactors=FALSE)
  })
  lib <- do.call(rbind,rows)
  R2_make_design(A,library=lib,name=label)
}

matched_design <- R2_design_from_selection(matched$selected,"Independent_BOP2_global_error_matched")
component05_design <- R2_design_from_selection(bonferroni_style$selected,"Independent_BOP2_component_cap_0p05")
write.csv(matched_design$library,file.path(outdir,"R2_C4_error_matched_positive_library.csv"),row.names=FALSE)
write.csv(component05_design$library,file.path(outdir,"R2_C4_cap_0p05_positive_library.csv"),row.names=FALSE)
write.csv(data.frame(method=c(matched_design$name,component05_design$name),
                     common_component_cap=c(tau_best,.05),
                     all_b20=c(matched_design$A["20"],component05_design$A["20"]),
                     all_b30=c(matched_design$A["30"],component05_design$A["30"]),
                     all_b40=c(matched_design$A["40"],component05_design$A["40"]),
                     max_global_point_null=c(matched$max_global,bonferroni_style$max_global)),
          file.path(outdir,"R2_C4_selected_allcomer_and_error_rate.csv"),row.names=FALSE)

# Exact operating characteristics, including previous 30 alternative cases,
# subgroup-only efficacy cases, reverse-order cases, and five point nulls.
scenarios <- R2_scenarios()
designs <- list(R2_proposed_design,R2_independent_design,matched_design,component05_design)
all_results <- vector("list",nrow(scenarios)*length(designs))
kernels <- new.env(hash=TRUE,parent=emptyenv())
get_kernel <- function(d,tp) {
  key <- paste(d$name,format(tp,nsmall=5,digits=9),sep="_")
  if (!exists(key,envir=kernels,inherits=FALSE)) assign(key,R2_positive_kernel(d,tp),envir=kernels)
  get(key,envir=kernels,inherits=FALSE)
}
count <- 0L
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i,]
  st20 <- R2_joint_states(20L,s$pi,s$theta_positive,s$theta_negative)
  st30 <- R2_joint_states(30L,s$pi,s$theta_positive,s$theta_negative)
  for (d in designs) {
    count <- count+1L
    z <- R2_evaluate_scenario(s$pi,s$theta_positive,s$theta_negative,d,
                              kernel=get_kernel(d,s$theta_positive),st20=st20,st30=st30)
    z$scenario <- s$scenario
    all_results[[count]] <- z
  }
}
res <- do.call(rbind,all_results)
rownames(res) <- NULL
write.csv(res,file.path(outdir,"R2_C4_all_four_methods_exact.csv"),row.names=FALSE)

percent <- res
prob_columns <- c("PRN_any","PRN_all","PRN_positive","P_false_allcomer","P_false_positive",
                  "pointwise_FWER","P_valid_allcomer","P_valid_positive","P_at_least_one_valid_claim",
                  "P_correct_positive_target","P_potentially_overbroad_allcomer")
for (nm in prob_columns) percent[[paste0(nm,"_pct")]] <- round(100*percent[[nm]],4)
for (group in unique(percent$scenario)) {
  write.csv(percent[percent$scenario==group,,drop=FALSE],
            file.path(outdir,paste0("R2_C4_",group,"_percent.csv")),row.names=FALSE)
}

# Double-check the calculated adjusted global point-null maximum.
adjusted_null <- res[res$scenario=="Point_global_null" & res$method==matched_design$name,,drop=FALSE]
stopifnot(nrow(adjusted_null)==length(R2_prevalences),
          abs(max(adjusted_null$PRN_any)-matched$max_global)<1e-10,
          matched$max_global<=R2_alpha+1e-12)

message("=== Exact error-matched comparator results ===")
message("Proposed actual maximum point-null global error: ",format(proposed_max,digits=12))
message("Selected common BOP2 component nominal cap: ",format(tau_best,digits=12))
message("Adjusted comparator maximum point-null global error: ",format(matched$max_global,digits=12))
message("Absolute difference in maximum error rate: ",format(abs(matched$max_global-proposed_max),digits=10))
message("Adjusted all-comer boundaries: ",paste(matched_design$A,collapse=","))
message("Illustrative component cap=0.05, actual global point-null max: ",format(bonferroni_style$max_global,digits=10))
message("This calibration is across the FIVE prespecified point-global-null prevalences only; it does not establish strong FWER control.")
message("=== Subgroup-only scenarios: comparison across all four methods (%) ===")
print(percent[percent$scenario=="Subgroup_only",
              c("pi","theta_positive","theta_negative","method","PRN_all_pct","PRN_positive_pct","PRN_any_pct")],
      row.names=FALSE)
message("=== Point-global-null scenarios: comparison across four methods (%) ===")
print(percent[percent$scenario=="Point_global_null",
              c("pi","method","PRN_any_pct","PRN_all_pct","PRN_positive_pct")],row.names=FALSE)
message("Files saved in: ",outdir)
