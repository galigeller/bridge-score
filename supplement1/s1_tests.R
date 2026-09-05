#!/usr/bin/env Rscript
# =============================================================================
# s1_tests.R — correctness anchors for Supplement S1. Run before s1_run.R.
#
#   1. b_true reproduces manuscript Table B1 exactly (0.96 / 0.65 / 0.36).
#   2. The corrected score reproduces the worked example of the
#      cluster-size-correction note (bridge ~1.00, small-cluster pure ~0.30).
#   3. Raw score sanity: y = (0.75, 0.25) -> b = 0.75.
#   4. Edge conventions: k = 1 -> all scores 0.
#   5. Generator determinism: same seed -> identical network.
#   6. End-to-end smoke: thesis curvature runs on one 10+10 network;
#      membership rows sum to 1; all scores in [0, 1].
#   7. Whole-count floor on the Table B1 network matches Appendix B's
#      enumeration for node 1 (achievable scores 0 and 0.56).
#
# All tests use stopifnot: silence up to the final message means pass.
# =============================================================================

source("curvature_CD.R")
source("s1_functions.R")

near <- function(x, y, tol = 1e-10) all(abs(x - y) < tol)

# --- the Appendix B network -------------------------------------------------
# Community 1 = {1,2,3,4} complete, Community 2 = {5,6,7,8} complete,
# bridge B ~ {1,2,3,5,7}, between-community edge {4,6}.
appendix_b_network <- function() {
  p <- 9
  A <- matrix(0L, p, p)
  edge <- function(i, j) { A[i, j] <<- 1L; A[j, i] <<- 1L }
  for (i in 1:3) for (j in (i + 1):4) edge(i, j)
  for (i in 5:7) for (j in (i + 1):8) edge(i, j)
  for (j in c(1, 2, 3, 5, 7)) edge(9, j)
  edge(4, 6)
  Pi <- rbind(matrix(c(1, 0), 4, 2, byrow = TRUE),
              matrix(c(0, 1), 4, 2, byrow = TRUE),
              c(0.5, 0.5))
  rownames(Pi) <- rownames(A) <- colnames(A) <-
    c(paste0("P1_", 1:4), paste0("P2_", 1:4), "B")
  list(A = A, Pi = Pi)
}

## Test 1: Table B1 ----------------------------------------------------------
net <- appendix_b_network()
tr <- true_bridge_score(net$A, net$Pi)
b <- setNames(tr$b_true, rownames(net$Pi))
stopifnot(near(b["B"], 0.96))                       # bridge: 4 * .6 * .4
stopifnot(near(round(b["P1_4"], 2), 0.65))          # node 4: 756/1156 = 0.6540
stopifnot(near(b["P1_1"], 0.36))                    # node 1: 4 * .9 * .1
stopifnot(near(tr$counts["B", ], c(3, 2)))          # whole counts (pure neighbors)
stopifnot(near(tr$counts["P1_1", ], c(3.5, 0.5)))   # fractional (bridge neighbor)
# raw score on the true structure: y_raw = counts/deg, participation of that.
b_rt <- setNames(tr$b_raw_true, rownames(net$Pi))
stopifnot(near(b_rt["B"], 0.96))                    # (3,2)/5 -> 4*.6*.4 (masses equal)
stopifnot(near(b_rt["P1_1"], 4 * (3.5/4) * (0.5/4)))# = 0.4375
stopifnot(near(b_rt["P1_4"], 4 * 0.75 * 0.25))      # (3,1)/4 -> 0.75
message("test 1 (Table B1 + raw-on-truth) .......... ok")

## Test 2: cluster-size-correction worked example -----------------------------
# 4 pure small nodes y=(1,0); one small node y=(.75,.25) with d=(3,1);
# 15 pure large nodes y=(0,1); bridge y=(.25,.75) with d=(2,6).
memb <- rbind(matrix(c(1, 0), 4, 2, byrow = TRUE),
              c(0.75, 0.25),
              matrix(c(0, 1), 15, 2, byrow = TRUE),
              c(0.25, 0.75))
degs <- c(rep(4, 4), 4, rep(4, 15), 8)              # d = y * deg stays whole
sc <- bridge_scores_detected(memb, degs)
stopifnot(near(colSums(memb), c(5, 16)))            # soft masses from the note
stopifnot(near(round(sc$b_corrected[21], 2), 1.00)) # bridge ~ 1.00  (0.99884)
stopifnot(near(round(sc$b_corrected[5], 2), 0.30))  # small pure ~ 0.30 (0.30286)
stopifnot(near(sc$b_raw[5], 0.75), near(sc$b_raw[21], 0.75))  # raw ties them
message("test 2 (correction worked example) ........ ok")

## Test 3: raw score ----------------------------------------------------------
stopifnot(near(participation(matrix(c(0.75, 0.25), 1), 2), 0.75))
message("test 3 (raw participation) ................ ok")

## Test 4: k = 1 convention ---------------------------------------------------
one <- bridge_scores_detected(matrix(1, 5, 1), rep(3, 5))
stopifnot(all(one$b_corrected == 0), all(one$b_raw == 0))
message("test 4 (k = 1 -> 0) ....................... ok")

## Test 5: seed determinism ---------------------------------------------------
g1 <- generate_network(10, 10, c(0.5, 0.5), seed = 42)
g2 <- generate_network(10, 10, c(0.5, 0.5), seed = 42)
stopifnot(identical(g1$A, g2$A))
message("test 5 (determinism) ...................... ok")

## Test 6: end-to-end smoke ---------------------------------------------------
g <- igraph::graph_from_adjacency_matrix(g1$A, mode = "undirected")
det <- run_curvature_cd(g, T_iter = 8, delta_prune = 0.1)
stopifnot(ncol(det$membership) >= 1)
stopifnot(near(rowSums(det$membership), rep(1, 21), tol = 1e-9))
sc <- bridge_scores_detected(det$membership, igraph::degree(g))
stopifnot(all(sc$b_corrected >= -1e-12 & sc$b_corrected <= 1 + 1e-12))
stopifnot(all(sc$b_raw >= -1e-12 & sc$b_raw <= 1 + 1e-12))
tr <- true_bridge_score(g1$A, g1$Pi)
stopifnot(all(tr$b_true >= 0 & tr$b_true <= 1))
message("test 6 (end-to-end smoke, k = ", sc$k, ") .......... ok")

## Test 7: whole-count floor vs Appendix B enumeration ------------------------
fl <- whole_count_floor(net$A, net$Pi, true_bridge_score(net$A, net$Pi)$b_true)
# 6 free edges (5 bridge-incident + the 4-6 edge) -> 64 assignments, exact.
# Appendix B enumerates node 1's two options: B's edge counted with community 1
# gives counts (4,0) and a score of 0 (gap 0.36); counted with community 2 the
# score moves toward the target. So the best gap is at most 0.36 and — the
# manuscript's claim — strictly positive: whole counts cannot hit 0.36 exactly.
stopifnot(fl$exact, fl$n_free_edges == 6)
stopifnot(all(is.finite(fl$floor_gap)), all(fl$floor_gap >= 0))
stopifnot(fl$floor_gap[1] <= 0.36 + 1e-9)
stopifnot(fl$floor_gap[1] > 1e-6)                   # cannot equal the target exactly
message("test 7 (whole-count floor, node-1 gap = ",
        round(fl$floor_gap[1], 3), ") ... ok")

message("all tests passed")
