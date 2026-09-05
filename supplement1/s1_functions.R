# =============================================================================
# s1_functions.R — everything Supplement S1 needs besides the thesis algorithm
#
# Contents:
#   planted_memberships / generate_network  the paper's data-generating networks
#   true_bridge_score                       the target b_true (manuscript Eq. B1–B4)
#   bridge_scores_detected                  corrected + raw score from a detected
#                                           membership matrix (manuscript Eq. 1–5)
#   whole_count_floor                       best score any whole-edge-count
#                                           read-out can achieve vs b_true
#
# Conventions (stated in the supplement text):
#   * K = 2 planted communities; community 1 is the LARGER one in unequal
#     splits, and pi_bridge[1] is the bridge's membership in community 1
#     (so an asymmetric bridge leans toward the larger community; swap the
#     entries of pi_bridge to flip this).
#   * k = 1 detected community  => all scores 0 (no boundary, no bridge evidence).
#   * d_l(v) = 0                => r_l(v) = 0 directly (never 0/0).
#   * zero-degree nodes are avoided by regeneration in generate_network.
# =============================================================================

# ---------------------------------------------------------------------------
# Generator: planted mixed-membership SBM (Tian et al. 2023, Definition 1)
# ---------------------------------------------------------------------------

planted_memberships <- function(n_comm1, n_comm2, pi_bridge) {
  stopifnot(length(pi_bridge) == 2, abs(sum(pi_bridge) - 1) < 1e-12)
  Pi <- rbind(
    matrix(c(1, 0), nrow = n_comm1, ncol = 2, byrow = TRUE),
    matrix(c(0, 1), nrow = n_comm2, ncol = 2, byrow = TRUE),
    matrix(pi_bridge, nrow = 1)
  )
  rownames(Pi) <- c(
    paste0("P1_", seq_len(n_comm1)),
    paste0("P2_", seq_len(n_comm2)),
    "B"
  )
  colnames(Pi) <- c("C1", "C2")
  Pi
}

#' Draw one network: edge {i,j} present with prob pi(i)' B pi(j).
#' Regenerates (up to max_tries) if any node has degree zero.
generate_network <- function(n_comm1, n_comm2, pi_bridge,
                             p_within = 0.8, p_between = 0.05,
                             seed = NULL, max_tries = 50) {
  if (!is.null(seed)) set.seed(seed)
  Pi <- planted_memberships(n_comm1, n_comm2, pi_bridge)
  B <- matrix(c(p_within, p_between, p_between, p_within), 2, 2)
  P <- Pi %*% B %*% t(Pi)
  p <- nrow(Pi)
  upper <- which(upper.tri(P))
  for (attempt in seq_len(max_tries)) {
    A <- matrix(0L, p, p)
    A[upper] <- rbinom(length(upper), 1, P[upper])
    A <- A + t(A)
    if (min(rowSums(A)) > 0) {
      dimnames(A) <- list(rownames(Pi), rownames(Pi))
      return(list(A = A, Pi = Pi, attempts = attempt))
    }
  }
  stop("generate_network: zero-degree node persisted after ", max_tries, " tries")
}

# ---------------------------------------------------------------------------
# Ground truth: manuscript Eq. (B1)-(B4) / Eq. (7)
# ---------------------------------------------------------------------------

#' b_true for every node: true memberships in place of detected ones.
#' A: p x p adjacency (0/1). Pi: p x K true membership matrix.
true_bridge_score <- function(A, Pi) {
  K <- ncol(Pi)
  p <- nrow(Pi)
  m <- colSums(Pi)                                   # community size = total membership
  Cnt <- A %*% Pi                                    # (B1) c_l(v), fractional
  denom <- matrix(m, p, K, byrow = TRUE) - Pi        # (B2) size excluding v's own share
  stopifnot(all(denom > 0 | Cnt == 0))
  R <- Cnt / denom                                   # (B2) rates
  R[Cnt == 0] <- 0
  s <- rowSums(R)
  Y <- R / ifelse(s > 0, s, 1)                       # (B3) corrected memberships
  b <- (K / (K - 1)) * (1 - rowSums(Y^2))            # (B4)
  b[s == 0] <- 0
  # The RAW score evaluated on the same true structure (participation of the
  # uncorrected memberships y_l = c_l/deg, manuscript Eq. 2 on Eq. 1). The pair
  # (b_raw_true, b_true) isolates the size correction from detection entirely:
  # both use the true communities, they differ only in the normalization.
  deg <- rowSums(A)
  Yraw <- Cnt / ifelse(deg > 0, deg, 1)
  b_raw <- (K / (K - 1)) * (1 - rowSums(Yraw^2))
  b_raw[deg == 0] <- 0
  list(b_true = as.numeric(b), b_raw_true = as.numeric(b_raw),
       y_true = Y, counts = Cnt)
}

# ---------------------------------------------------------------------------
# Scores from a detected membership matrix (manuscript Eq. 1-5 and Eq. 2 raw)
# ---------------------------------------------------------------------------

participation <- function(U, k) {
  if (k <= 1) return(rep(0, nrow(U)))
  (k / (k - 1)) * (1 - rowSums(U^2))
}

#' Corrected and raw bridge scores from the detection output.
#' membership: p x k matrix with rows y(v) = d(v)/deg(v)  (thesis Eq. 4.4,
#'             i.e. run_curvature_cd()$membership). degrees: degree in G.
bridge_scores_detected <- function(membership, degrees) {
  p <- nrow(membership)
  k <- ncol(membership)
  if (k <= 1) {
    return(list(b_corrected = rep(0, p), b_raw = rep(0, p), k = k))
  }
  D <- membership * degrees                          # d_l(v), whole by construction
  mass <- colSums(membership)                        # Eq. (3): soft community mass
  denom <- matrix(mass, p, k, byrow = TRUE) - membership
  stopifnot(all(denom > 0 | D == 0))
  R <- D / denom                                     # Eq. (4): rates
  R[D == 0] <- 0
  s <- rowSums(R)
  Ytilde <- R / ifelse(s > 0, s, 1)                  # Eq. (5)
  b_corrected <- participation(Ytilde, k)            # Eq. (2) on corrected memberships
  b_corrected[s == 0] <- 0
  b_raw <- participation(membership, k)              # Eq. (2) on raw memberships
  list(b_corrected = as.numeric(b_corrected), b_raw = as.numeric(b_raw),
       k = k, y_corrected = Ytilde)
}

# ---------------------------------------------------------------------------
# Discretization floor: whole counts vs the fractional target
# ---------------------------------------------------------------------------

#' The target's counts are fractional wherever a neighbor's membership is
#' fractional (Appendix B), but any edge-partition read-out assigns each edge
#' wholly to one community. This computes, per node, the closest score a
#' whole-count read-out can achieve under the correct two-community structure:
#' edges between same-community pure nodes are fixed to that community; every
#' other edge (bridge-incident and pure-pure between-community) is free, and
#' all joint assignments of the free edges are enumerated (exhaustively up to
#' 2^max_exact free edges, Monte Carlo above that — MC makes the reported
#' floor conservative, i.e. never smaller than the true one).
whole_count_floor <- function(A, Pi, b_true, max_exact = 14, n_mc = 4096) {
  p <- nrow(A)
  k <- 2L
  el <- which(upper.tri(A) & A == 1, arr.ind = TRUE)  # edge list (u < v)
  hard <- apply(Pi, 1, which.max)
  pure <- apply(Pi, 1, max) == 1
  u <- el[, 1]; v <- el[, 2]
  fixed <- pure[u] & pure[v] & hard[u] == hard[v]
  free <- which(!fixed)
  n_free <- length(free)

  # base counts from the fixed edges
  D0 <- matrix(0, p, k)
  for (e in which(fixed)) {
    cl <- hard[u[e]]
    D0[u[e], cl] <- D0[u[e], cl] + 1
    D0[v[e], cl] <- D0[v[e], cl] + 1
  }

  assignments <- if (n_free <= max_exact) {
    as.matrix(expand.grid(rep(list(1:2), n_free)))
  } else {
    matrix(sample(1:2, n_mc * n_free, replace = TRUE), n_mc, n_free)
  }
  if (n_free == 0) assignments <- matrix(integer(0), 1, 0)

  degrees <- rowSums(A)
  best_gap <- rep(Inf, p)
  best_b <- rep(NA_real_, p)
  for (a in seq_len(nrow(assignments))) {
    D <- D0
    for (idx in seq_len(n_free)) {
      e <- free[idx]
      cl <- assignments[a, idx]
      D[u[e], cl] <- D[u[e], cl] + 1
      D[v[e], cl] <- D[v[e], cl] + 1
    }
    memb <- D / degrees
    sc <- bridge_scores_detected(memb, degrees)$b_corrected
    gap <- abs(sc - b_true)
    improved <- gap < best_gap
    best_gap[improved] <- gap[improved]
    best_b[improved] <- sc[improved]
  }
  list(floor_gap = best_gap, closest_b = best_b, n_free_edges = n_free,
       exact = n_free <= max_exact)
}

# ---------------------------------------------------------------------------
# Node bookkeeping for results tables
# ---------------------------------------------------------------------------

node_types <- function(A, Pi) {
  p <- nrow(Pi)
  bridge <- which(rownames(Pi) == "B")
  type <- rep("pure", p)
  type[bridge] <- "bridge"
  type[A[, bridge] == 1 & seq_len(p) != bridge] <- "neighbor"
  community <- c("C1", "C2", "bridge")[apply(cbind(Pi[, 1] == 1, Pi[, 2] == 1, TRUE), 1,
                                             which.max)]
  data.frame(node = rownames(Pi), type = type, community = community,
             stringsAsFactors = FALSE)
}
