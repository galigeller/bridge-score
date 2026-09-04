#!/usr/bin/env Rscript
# =============================================================================
# s1_run.R — Supplement S1 simulation: the read-out on known networks
#
#   Grid mirrors the thesis's simulation conditions:
#     Part A  5 separability levels x 3 bridge memberships, balanced 10+10
#             (p_within, p_between) in {(.8,0), (.8,.05), (.75,.10),
#             (.7,.15), (.6,.30)} — the thesis's levels
#     Part B  2 unequal splits (13+7, 15+5) x 3 bridge memberships at the
#             baseline separability (.8, .05)
#   = 21 conditions x REPS networks. Thesis curvature detection (T_iter = 8,
#   delta_prune = 0.1, automatic delta_cut) on every network, corrected + raw
#   scores vs b_true. The whole-count floor is computed at the paper's base
#   condition (10+10, pi = (.5,.5), baseline separability) — the condition
#   Appendix B refers to.
#
# Usage:    Rscript s1_run.R                       (defaults: 200 reps, 1 core)
#           S1_REPS=20 Rscript s1_run.R            (quick pass)
#           S1_CORES=8 Rscript s1_run.R            (parallel over replications;
#                                                   results are identical —
#                                                   every rep has its own seed)
# Output:   results/cond_<i>.csv (one per condition, skipped if present —
#           delete a file to re-run its condition), then results/s1_results.csv
# Runtime:  ~2-6 s per network with CRAN 'transport'; ~5-7 h single-core for
#           the full 21 x 200, ~1 h with S1_CORES=8.
# =============================================================================

source("curvature_CD.R")
source("s1_functions.R")

REPS    <- as.integer(Sys.getenv("S1_REPS", "200"))
CORES   <- as.integer(Sys.getenv("S1_CORES", "1"))
T_ITER  <- 8      # thesis settings (final_thesis, Methods)
D_PRUNE <- 0.1

# The thesis's five separability levels (final_thesis, Simulation Conditions):
# p_within is reduced alongside p_between to keep the within/between ratio realistic.
separability <- data.frame(
  p_within  = c(0.80, 0.80, 0.75, 0.70, 0.60),
  p_between = c(0.00, 0.05, 0.10, 0.15, 0.30)
)
BASE_PW <- 0.80
BASE_PB <- 0.05   # the paper's / thesis's baseline

part_a <- merge(separability,
                data.frame(pi1 = c(0.5, 0.7, 0.9)))   # 15 cells, balanced
part_a$split <- "10+10"
part_b <- expand.grid(pi1 = c(0.5, 0.7, 0.9),
                      split = c("13+7", "15+5"),
                      stringsAsFactors = FALSE)       # 6 cells, baseline separability
part_b$p_within <- BASE_PW
part_b$p_between <- BASE_PB
conditions <- rbind(part_a[, c("pi1", "split", "p_within", "p_between")],
                    part_b[, c("pi1", "split", "p_within", "p_between")])
rownames(conditions) <- NULL
splits <- list("10+10" = c(10, 10), "13+7" = c(13, 7), "15+5" = c(15, 5))

dir.create("results", showWarnings = FALSE)

run_one <- function(cond_idx, rep_idx) {
  cond <- conditions[cond_idx, ]
  ns <- splits[[cond$split]]
  pi_bridge <- c(cond$pi1, 1 - cond$pi1)
  seed <- 10000L * cond_idx + rep_idx

  gen <- generate_network(ns[1], ns[2], pi_bridge,
                          p_within = cond$p_within, p_between = cond$p_between,
                          seed = seed)
  truth <- true_bridge_score(gen$A, gen$Pi)

  g <- igraph::graph_from_adjacency_matrix(gen$A, mode = "undirected")
  t0 <- Sys.time()
  det <- run_curvature_cd(g, T_iter = T_ITER, delta_prune = D_PRUNE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sc <- bridge_scores_detected(det$membership, igraph::degree(g))

  base_condition <- cond$pi1 == 0.5 && cond$split == "10+10" &&
    cond$p_within == BASE_PW && cond$p_between == BASE_PB
  fl <- if (base_condition) {
    whole_count_floor(gen$A, gen$Pi, truth$b_true)
  } else {
    list(floor_gap = rep(NA_real_, nrow(gen$Pi)),
         closest_b = rep(NA_real_, nrow(gen$Pi)),
         n_free_edges = NA_integer_, exact = NA)
  }

  info <- node_types(gen$A, gen$Pi)
  data.frame(
    cond        = cond_idx,
    pi1         = cond$pi1,
    split       = cond$split,
    p_within    = cond$p_within,
    p_between   = cond$p_between,
    rep         = rep_idx,
    seed        = seed,
    node        = info$node,
    type        = info$type,
    community   = info$community,
    b_true      = truth$b_true,
    b_corrected = sc$b_corrected,
    b_raw       = sc$b_raw,
    k_detected  = sc$k,
    delta_cut   = det$delta_cut,
    floor_gap   = fl$floor_gap,
    floor_b     = fl$closest_b,
    curv_secs   = secs,
    stringsAsFactors = FALSE
  )
}

for (cond_idx in seq_len(nrow(conditions))) {
  out_file <- file.path("results", sprintf("cond_%d.csv", cond_idx))
  if (file.exists(out_file)) {
    message("condition ", cond_idx, " exists, skipping")
    next
  }
  cond <- conditions[cond_idx, ]
  message(sprintf(
    "condition %d/%d: pi = (%.1f, %.1f), split %s, (pw, pb) = (%.2f, %.2f), %d reps",
    cond_idx, nrow(conditions), cond$pi1, 1 - cond$pi1,
    cond$split, cond$p_within, cond$p_between, REPS))
  rows <- if (CORES > 1) {
    parallel::mclapply(seq_len(REPS), function(r) run_one(cond_idx, r),
                       mc.cores = CORES)
  } else {
    lapply(seq_len(REPS), function(r) {
      if (r %% 25 == 0) message("  rep ", r, "/", REPS)
      run_one(cond_idx, r)
    })
  }
  failed <- vapply(rows, function(x) !is.data.frame(x), logical(1))
  if (any(failed)) stop("condition ", cond_idx, ": ", sum(failed),
                        " replications failed; first error: ",
                        conditionMessage(attr(rows[[which(failed)[1]]], "condition")))
  write.csv(do.call(rbind, rows), out_file, row.names = FALSE)
}

all_files <- file.path("results", sprintf("cond_%d.csv", seq_len(nrow(conditions))))
results <- do.call(rbind, lapply(all_files, read.csv, stringsAsFactors = FALSE))
write.csv(results, file.path("results", "s1_results.csv"), row.names = FALSE)
message("done: ", nrow(results), " node-level rows -> results/s1_results.csv")
