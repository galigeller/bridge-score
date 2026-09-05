#!/usr/bin/env Rscript
# =============================================================================
# s1_sensitivity.R — sensitivity of the read-out to the flow parameters
#
# Context: Tian et al. (2023) give no rule for the number of Ricci-flow
# iterations T — they hand-tune it and leave systematic selection to future
# work. The main S1 grid uses the thesis settings (T = 8, delta_prune = 0.1).
# This script re-runs the BASE CONDITION ONLY (pi = (.5,.5), 10+10 split,
# p_within = .8, p_between = .05 — condition 2 of s1_run.R) on the SAME 200
# networks (same seeds), varying one parameter at a time:
#
#   T_iter      in {4, 16, 30}   (delta_prune held at 0.1)
#   delta_prune in {0.05, 0.2}   (T_iter held at 8)
#
# The T = 8 / delta_prune = 0.1 reference row is read from the main run's
# results/cond_2.csv — run s1_run.R first (at least condition 2).
#
# Usage:    Rscript s1_sensitivity.R           (S1_REPS / S1_CORES as in s1_run.R)
# Output:   results/sens_<param>_<value>.csv (resumable, one per setting),
#           results/sens_summary.csv, and a printed summary table.
# Runtime:  5 settings x 200 networks; ~45 min with S1_CORES=6 on an M-series Mac.
# =============================================================================

source("curvature_CD.R")
source("s1_functions.R")

REPS      <- as.integer(Sys.getenv("S1_REPS", "200"))
CORES     <- as.integer(Sys.getenv("S1_CORES", "1"))
BASE_COND <- 2L          # seed block of the base condition in s1_run.R

settings <- rbind(
  data.frame(param = "T_iter",      value = c(4, 16, 30),
             T_iter = c(4, 16, 30), delta_prune = 0.1),
  data.frame(param = "delta_prune", value = c(0.05, 0.2),
             T_iter = 8,            delta_prune = c(0.05, 0.2))
)

dir.create("results", showWarnings = FALSE)

run_one <- function(setting, rep_idx) {
  seed <- 10000L * BASE_COND + rep_idx        # identical networks to cond 2
  gen <- generate_network(10, 10, c(0.5, 0.5),
                          p_within = 0.8, p_between = 0.05, seed = seed)
  truth <- true_bridge_score(gen$A, gen$Pi)
  g <- igraph::graph_from_adjacency_matrix(gen$A, mode = "undirected")
  det <- run_curvature_cd(g, T_iter = setting$T_iter,
                          delta_prune = setting$delta_prune)
  sc <- bridge_scores_detected(det$membership, igraph::degree(g))
  info <- node_types(gen$A, gen$Pi)
  data.frame(param = setting$param, value = setting$value, rep = rep_idx,
             seed = seed, node = info$node, type = info$type,
             b_true = truth$b_true, b_corrected = sc$b_corrected,
             b_raw = sc$b_raw, k_detected = sc$k, delta_cut = det$delta_cut,
             stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  out_file <- file.path("results", sprintf("sens_%s_%s.csv", s$param, s$value))
  if (file.exists(out_file)) { message(out_file, " exists, skipping"); next }
  message(sprintf("setting %d/%d: %s = %s, %d reps",
                  i, nrow(settings), s$param, s$value, REPS))
  rows <- if (CORES > 1) {
    parallel::mclapply(seq_len(REPS), function(r) run_one(s, r), mc.cores = CORES)
  } else {
    lapply(seq_len(REPS), function(r) {
      if (r %% 25 == 0) message("  rep ", r, "/", REPS)
      run_one(s, r)
    })
  }
  failed <- vapply(rows, function(x) !is.data.frame(x), logical(1))
  if (any(failed)) stop(s$param, "=", s$value, ": ", sum(failed),
                        " replications failed")
  write.csv(do.call(rbind, rows), out_file, row.names = FALSE)
}

## ---- summary table ---------------------------------------------------------
ref_file <- file.path("results", "cond_2.csv")
if (!file.exists(ref_file))
  stop("results/cond_2.csv not found — run s1_run.R (condition 2) first; ",
       "it is the T = 8 / delta_prune = 0.1 reference")
ref <- read.csv(ref_file, stringsAsFactors = FALSE)

summarise_setting <- function(d, param, value) {
  per_net <- split(d, d$rep)
  k2   <- mean(vapply(per_net, function(x) x$k_detected[1] == 2, logical(1)))
  mk   <- mean(vapply(per_net, function(x) x$k_detected[1], numeric(1)))
  top  <- mean(vapply(per_net, function(x)
                x$node[which.max(x$b_corrected)] == "B", logical(1)))
  br   <- d[d$type == "bridge", ]
  data.frame(param = param, value = value, n_networks = length(per_net),
             share_k2 = round(k2, 3), mean_k = round(mk, 2),
             bridge_mae = round(mean(abs(br$b_corrected - br$b_true)), 3),
             bridge_top_ranked = round(top, 3))
}

summ <- rbind(
  summarise_setting(ref, "T_iter (reference)", 8),
  do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
    s <- settings[i, ]
    d <- read.csv(file.path("results",
                            sprintf("sens_%s_%s.csv", s$param, s$value)),
                  stringsAsFactors = FALSE)
    summarise_setting(d, s$param, s$value)
  }))
)
write.csv(summ, file.path("results", "sens_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
message("summary written to results/sens_summary.csv")
