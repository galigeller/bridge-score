# =============================================================================
# ricci_flow_tests.R — hand-verified ORC and Ricci flow tests (thesis code)
#
# Copied VERBATIM from the thesis repository
#   galigeller/Thesis-sensitivity-of-SBM-prior-to-MMSBM-data
#   Tests/testthat/ricci_flow_tests.R
# with one adaptation: the source() path points to the local curvature_CD.R
# (the thesis file resolves it relative to the repo root).
#
# Run with:  Rscript ricci_flow_tests.R
# =============================================================================

library(testthat)
library(igraph)

# Hand-verified ORC and Ricci flow tests.
#
# All ORC values are for the paper's m_v = uniform over neighbors (no idleness),
# kappa = 1 - W1(m_u, m_v) / d(u,v).
# compute_orc() is called on the original graph directly (not the line graph)
# to test the mathematical formula; this matches how curvature_tests.R works.
#
# Edge ordering in all test graphs matches igraph's insertion order,
# which is also what as_edgelist() and compute_orc() iterate over.

source("curvature_CD.R")   # adapted path (thesis: algorithm_implementations/curvature_CD.R)

# ---------------------------------------------------------------------------
# Graph constructors
# ---------------------------------------------------------------------------

# Two triangles {1,2,3} and {4,5,6} joined by bridge {3,4}.
# Edge insertion order (1-indexed):
#   e1={1,2}  e2={1,3}  e3={2,3}  e4={3,4}  e5={4,5}  e6={4,6}  e7={5,6}
# Node degrees: deg(1)=deg(2)=deg(5)=deg(6)=2,  deg(3)=deg(4)=3.
make_two_triangle_bridge <- function() {
  g <- make_empty_graph(n = 6, directed = FALSE)
  add_edges(g, c(1,2, 1,3, 2,3, 3,4, 4,5, 4,6, 5,6))
}

# Bowtie: two triangles {1,2,3} and {3,4,5} sharing node 3.
# Edge insertion order:
#   e1={1,2}  e2={1,3}  e3={2,3}  e4={3,4}  e5={3,5}  e6={4,5}
# Node degrees: deg(1)=deg(2)=deg(4)=deg(5)=2,  deg(3)=4.
make_bowtie <- function() {
  g <- make_empty_graph(n = 5, directed = FALSE)
  add_edges(g, c(1,2, 1,3, 2,3, 3,4, 3,5, 4,5))
}

# Set all edge weights to 1 and run compute_orc on the original graph.
orc_on_graph <- function(g) {
  g <- set_edge_attr(g, "weight", value = rep(1.0, ecount(g)))
  compute_orc(g)
}

# ===========================================================================
# Section A: ORC on hand-verified graphs
# ===========================================================================

# ---------------------------------------------------------------------------
# Example 1 / Example 3: Triangle K_3 and complete graphs K_n
#
# For K_n (n >= 3):
#   N(u) and N(v) share n-2 nodes (matched at cost 0).
#   The only mismatch is 1/(n-1) of mass moved distance 1.
#   W1 = 1/(n-1),  kappa = 1 - 1/(n-1) = (n-2)/(n-1).
# ---------------------------------------------------------------------------

test_that("K_3: all 3 edges have kappa = 1/2", {
  g      <- make_full_graph(3, directed = FALSE)
  kappas <- orc_on_graph(g)
  expect_length(kappas, 3)
  expect_equal(kappas, rep(0.5, 3), tolerance = 1e-9)
})

test_that("K_4: all 6 edges have kappa = 2/3", {
  g      <- make_full_graph(4, directed = FALSE)
  kappas <- orc_on_graph(g)
  expect_length(kappas, 6)
  expect_equal(kappas, rep(2/3, 6), tolerance = 1e-9)
  expect_true(sd(kappas) < 1e-12)
})

test_that("K_5: all 10 edges have kappa = 3/4", {
  g      <- make_full_graph(5, directed = FALSE)
  kappas <- orc_on_graph(g)
  expect_length(kappas, 10)
  expect_equal(kappas, rep(3/4, 10), tolerance = 1e-9)
  expect_true(sd(kappas) < 1e-12)
})

test_that("K_n: kappa = (n-2)/(n-1) for n = 3..6, all edges equal", {
  for (n in 3:6) {
    g        <- make_full_graph(n, directed = FALSE)
    kappas   <- orc_on_graph(g)
    expected <- (n - 2) / (n - 1)
    expect_true(
      all(abs(kappas - expected) < 1e-9),
      info = sprintf("K_%d: expected %.6f, got range [%.6f, %.6f]",
                     n, expected, min(kappas), max(kappas))
    )
    expect_true(sd(kappas) < 1e-12,
                info = sprintf("K_%d curvature is not uniform", n))
  }
})

# ---------------------------------------------------------------------------
# Example 4: Even cycle C_n — zero curvature everywhere
#
# For edge {u,v} in C_4:  N(u)={a,w},  N(v)={u's-other, v's-other}.
# Every pair (x in N(u), y in N(v)) has d(x,y)=1 in C_4,
# so the cost matrix is all-ones and W1 = 1 = d(u,v),  kappa = 0.
# The same holds for all even cycles C_n (n >= 4).
# ---------------------------------------------------------------------------

test_that("C_4: all 4 edges have kappa = 0", {
  g      <- make_ring(4, directed = FALSE)
  kappas <- orc_on_graph(g)
  expect_length(kappas, 4)
  expect_equal(kappas, rep(0.0, 4), tolerance = 1e-9)
})

test_that("Even cycles C_6 and C_8: all edges have kappa = 0", {
  for (n in c(6, 8)) {
    g      <- make_ring(n, directed = FALSE)
    kappas <- orc_on_graph(g)
    expect_equal(kappas, rep(0.0, n), tolerance = 1e-9,
                 info = sprintf("C_%d: some edges deviate from zero", n))
  }
})

# ---------------------------------------------------------------------------
# Example 2: Two triangles connected by a bridge
#
# Graph: {1,2,3} triangle + {3,4} bridge + {4,5,6} triangle.
# Degrees: deg(1)=deg(2)=deg(5)=deg(6)=2,  deg(3)=deg(4)=3.
#
# Bridge edge e4={3,4}  (hand-computed below):
#   m_3 = (1/3, 1/3, 1/3) on {1, 2, 4}
#   m_4 = (1/3, 1/3, 1/3) on {3, 5, 6}
#   Cost matrix (rows=N(3), cols=N(4)):
#         3    5    6
#   1  [  1    3    3 ]
#   2  [  1    3    3 ]
#   4  [  1    1    1 ]
#   Optimal plan: 1/3 from 1->3 (cost 1), 1/3 from 4->5 (cost 1), 1/3 from 2->6 (cost 3).
#   W1 = (1 + 1 + 3)/3 = 5/3,  kappa = 1 - 5/3 = -2/3.
#
# Pure-internal edges e1={1,2} and e7={5,6}  (both endpoints have degree 2):
#   Same calculation as an isolated K_3:  kappa = 1/2.
#
# Semi-internal edges e2={1,3}, e3={2,3}, e5={4,5}, e6={4,6}
#   (one endpoint has degree 3 due to the bridge):
#   Worked example for e2={1,3}:
#     m_1 = (1/2, 1/2) on {2, 3}
#     m_3 = (1/3, 1/3, 1/3) on {1, 2, 4}
#     Set x_{2,2}=1/3 (cost 0), then minimise over remaining mass ->
#     W1 = 2/3,  kappa = 1/3.
#   e3, e5, e6 follow by symmetry.
# ---------------------------------------------------------------------------

test_that("Two-triangle-bridge: bridge edge e4={3,4} has kappa = -2/3", {
  g      <- make_two_triangle_bridge()
  kappas <- orc_on_graph(g)
  expect_length(kappas, 7)
  expect_equal(kappas[4], -2/3, tolerance = 1e-9)
})

test_that("Two-triangle-bridge: pure-internal edges e1={1,2} and e7={5,6} have kappa = 1/2", {
  g      <- make_two_triangle_bridge()
  kappas <- orc_on_graph(g)
  expect_equal(kappas[1], 0.5, tolerance = 1e-9,
               info = "edge e1={1,2}: both endpoints degree 2")
  expect_equal(kappas[7], 0.5, tolerance = 1e-9,
               info = "edge e7={5,6}: both endpoints degree 2")
})

test_that("Two-triangle-bridge: semi-internal edges e2,e3,e5,e6 have kappa = 1/3", {
  g      <- make_two_triangle_bridge()
  kappas <- orc_on_graph(g)
  # e2={1,3}, e3={2,3}: one endpoint is node 3 (degree 3)
  # e5={4,5}, e6={4,6}: one endpoint is node 4 (degree 3)
  for (i in c(2, 3, 5, 6)) {
    expect_equal(kappas[i], 1/3, tolerance = 1e-9,
                 info = sprintf("edge index %d should have kappa = 1/3, got %.6f",
                                i, kappas[i]))
  }
})

test_that("Two-triangle-bridge: bridge curvature is strictly below every internal curvature", {
  g      <- make_two_triangle_bridge()
  kappas <- orc_on_graph(g)
  expect_true(
    kappas[4] < min(kappas[-4]),
    info = sprintf("bridge kappa (%.4f) should be below min internal (%.4f)",
                   kappas[4], min(kappas[-4]))
  )
})

test_that("Two-triangle-bridge: all internal curvatures are positive", {
  g      <- make_two_triangle_bridge()
  kappas <- orc_on_graph(g)
  expect_true(all(kappas[-4] > 0),
              info = "all non-bridge edges should have positive curvature")
})

# ---------------------------------------------------------------------------
# Example 5: Bowtie (two triangles sharing node 3, no bridge edge)
#
# Graph: triangle {1,2,3} and triangle {3,4,5}, shared node 3.
# Degrees: deg(1)=deg(2)=deg(4)=deg(5)=2,  deg(3)=4.
#
# Pure-internal edges e1={1,2} and e6={4,5}  (both endpoints degree 2):
#   Same as K_3:  kappa = 1/2.
#
# Edges incident on node 3: e2={1,3}, e3={2,3}, e4={3,4}, e5={3,5}.
#   Worked example for e2={1,3}:
#     m_1 = (1/2, 1/2) on {2, 3}
#     m_3 = (1/4, 1/4, 1/4, 1/4) on {1, 2, 4, 5}
#     Optimal: x_{2,2}=1/4 (cost 0), x_{2,1}=1/4 (cost 1),
#              x_{3,4}=1/4 (cost 1), x_{3,5}=1/4 (cost 1).
#     W1 = 3/4,  kappa = 1 - 3/4 = 1/4.
#   e3, e4, e5 follow by symmetry (deg-4 node is always node 3).
# ---------------------------------------------------------------------------

test_that("Bowtie: pure-internal edges e1={1,2} and e6={4,5} have kappa = 1/2", {
  g      <- make_bowtie()
  kappas <- orc_on_graph(g)
  expect_length(kappas, 6)
  expect_equal(kappas[1], 0.5, tolerance = 1e-9,
               info = "edge e1={1,2}: both endpoints degree 2")
  expect_equal(kappas[6], 0.5, tolerance = 1e-9,
               info = "edge e6={4,5}: both endpoints degree 2")
})

test_that("Bowtie: edges e2,e3,e4,e5 incident on shared node 3 have kappa = 1/4", {
  g      <- make_bowtie()
  kappas <- orc_on_graph(g)
  for (i in c(2, 3, 4, 5)) {
    expect_equal(kappas[i], 0.25, tolerance = 1e-9,
                 info = sprintf("edge %d (incident on node 3, deg 4): got %.6f",
                                i, kappas[i]))
  }
})

test_that("Bowtie: exactly two distinct curvature values (1/4 and 1/2), all positive", {
  g             <- make_bowtie()
  kappas        <- orc_on_graph(g)
  unique_kappas <- sort(unique(round(kappas, 9)))
  expect_equal(length(unique_kappas), 2,
               info = "Bowtie should have exactly 2 curvature levels")
  expect_equal(unique_kappas[1], 0.25, tolerance = 1e-9)
  expect_equal(unique_kappas[2], 0.50, tolerance = 1e-9)
  expect_true(all(kappas > 0),
              info = "Bowtie has no bridge edge — all curvatures should be positive")
})

test_that("Bowtie: (1 - kappa) factor is larger for node-3 edges than for pure-internal edges", {
  # This verifies the Ricci flow update direction:
  # edges with lower kappa (incident on node 3) receive a larger (1-kappa) multiplier
  # and therefore grow relative to the pure-internal edges.
  g               <- make_bowtie()
  kappas          <- orc_on_graph(g)
  factor_internal <- mean(1 - kappas[c(1, 6)])    # edges {1,2} and {4,5}
  factor_incident <- mean(1 - kappas[c(2, 3, 4, 5)])  # edges on node 3
  expect_true(
    factor_incident > factor_internal,
    info = sprintf("node-3 edges factor (%.3f) should exceed internal factor (%.3f)",
                   factor_incident, factor_internal)
  )
})

# ===========================================================================
# Section B: Ricci flow qualitative and community-detection tests
# ===========================================================================

test_that("ricci_flow: evolved line graph has no more edges than original line graph", {
  g       <- make_two_triangle_bridge()
  evolved <- ricci_flow(g, T_iter = 1, delta_prune = 0.0)
  lg      <- make_line_graph(g)
  expect_lte(ecount(evolved), ecount(lg))
})

test_that("ricci_flow: all evolved edge weights are positive and finite", {
  g       <- make_two_triangle_bridge()
  evolved <- ricci_flow(g, T_iter = 10, delta_prune = 0.0)
  ws      <- E(evolved)$weight
  expect_true(all(is.finite(ws)),   info = "all evolved weights must be finite")
  expect_true(all(ws > 0),          info = "all evolved weights must be positive")
})

test_that("ricci_flow: normalisation keeps mean edge weight = 1 after each step", {
  # Eq. 4.3: sum(w_new) = M (number of edges), so mean = 1.
  g       <- make_two_triangle_bridge()
  evolved <- ricci_flow(g, T_iter = 5, delta_prune = 0.0)
  ws      <- E(evolved)$weight
  # Mean should be 1 (after normalisation the sum equals the current edge count).
  expect_equal(mean(ws), 1.0, tolerance = 1e-9)
})

test_that("Two-triangle-bridge: run_curvature_cd finds >= 2 communities", {
  g      <- make_two_triangle_bridge()
  set.seed(42)
  result <- run_curvature_cd(g, T_iter = 50, delta_prune = 0.0)
  n_comm <- ncol(result$membership)
  expect_true(n_comm >= 2,
              info = sprintf("Expected >= 2 communities, got %d", n_comm))
})

test_that("Two-triangle-bridge: nodes 1&2 co-cluster; nodes 5&6 co-cluster; the two groups differ", {
  g      <- make_two_triangle_bridge()
  set.seed(42)
  result <- run_curvature_cd(g, T_iter = 50, delta_prune = 0.0)
  skip_if(ncol(result$membership) < 2, "Algorithm did not produce >= 2 communities")

  mem <- result$membership
  dom <- function(v) which.max(mem[v, ])

  # Nodes 1 and 2 should have their dominant community in common.
  expect_equal(dom(1), dom(2),
               info = "Nodes 1 and 2 should share a dominant community")

  # Nodes 5 and 6 should have their dominant community in common.
  expect_equal(dom(5), dom(6),
               info = "Nodes 5 and 6 should share a dominant community")

  # The two groups should be in different communities.
  expect_false(dom(1) == dom(5),
               info = "Nodes {1,2} and nodes {5,6} should be in different communities")
})

test_that("Two-triangle-bridge: bridge-endpoint nodes 3&4 each have non-zero membership in >= 1 community", {
  g      <- make_two_triangle_bridge()
  set.seed(42)
  result <- run_curvature_cd(g, T_iter = 50, delta_prune = 0.0)
  skip_if(ncol(result$membership) < 1, "No communities detected")

  mem <- result$membership
  expect_true(sum(mem[3, ] > 0) >= 1,
              info = "Node 3 (bridge endpoint) should belong to at least one community")
  expect_true(sum(mem[4, ] > 0) >= 1,
              info = "Node 4 (bridge endpoint) should belong to at least one community")
})
