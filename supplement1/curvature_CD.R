# =============================================================================
# curvature_CD.R — line-graph curvature community detection (thesis code)
#
# Copied VERBATIM from the thesis repository
#   galigeller/Thesis-sensitivity-of-SBM-prior-to-MMSBM-data
#   algorithm_implementations/curvature_CD.R
# Line-graph pipeline only (compute_orc, select_cut_threshold, ricci_flow,
# infer_membership, run_curvature_cd). The original-graph (G) variant in that
# file is not used in the paper and is omitted here. DO NOT EDIT — any change
# to the algorithm belongs in the thesis repo, not in the supplement.
# =============================================================================

library(igraph)
library(lpSolve)
library(transport)

#' Step 1 / Sections 2.2–2.3 (Curvature definition)
compute_orc <- function(line_graph) {
  # Section 2.2 introduces nodes' neighborhoods; we grab all shortest-path
  # distances because Eq. 2.3 defines the transport cost via d(u', v') between
  # neighbors u' ∈ N(u) and v' ∈ N(v) of the line-graph edge {u,v}.
  distance_matrix <- distances(line_graph, weights = E(line_graph)$weight)
  adjacency_lists <- as_adj_list(line_graph)
  line_edges <- as_edgelist(line_graph, names = FALSE)
  curvatures <- numeric(ecount(line_graph))

  for (edge_idx in seq_len(ecount(line_graph))) {
    node_u <- line_edges[edge_idx, 1]
    node_v <- line_edges[edge_idx, 2]
    neighborhood_u <- adjacency_lists[[node_u]]
    neighborhood_v <- adjacency_lists[[node_v]]

    if (length(neighborhood_u) == 0 || length(neighborhood_v) == 0) {
      curvatures[edge_idx] <- 0
    } else {
      # match Eq. 2.2 notation: deg_u := degree(u)
      deg_u <- degree(line_graph, node_u)
      deg_v <- degree(line_graph, node_v)
      # uniform probability m_u assigns weight 1/deg(u) to each neighbor
      measure_u <- rep(1 / deg_u, deg_u)
      measure_v <- rep(1 / deg_v, deg_v)

    transport_cost <- distance_matrix[neighborhood_u, neighborhood_v, drop = FALSE]
    if (!all(is.finite(transport_cost))) {
      finite_entries <- transport_cost[is.finite(transport_cost)]
      penalty <- if (length(finite_entries) > 0) 10 * max(finite_entries) else 1e6
      transport_cost[!is.finite(transport_cost)] <- penalty
    }

    # Eq. 2.3: compute Wasserstein-1 distance using transport::wasserstein.
    # It solves the same linear transport problem with uniform marginals.
    w1 <- transport::wasserstein(
      a = measure_u,
      b = measure_v,
      cost = transport_cost
    )

    if (!is.na(w1)) {
      edge_length <- E(line_graph)$weight[edge_idx]
      if (is.finite(edge_length) && edge_length > 0) {
        # Eq. 2.2 final ratio W1 / d(u,v)
        curvatures[edge_idx] <- 1 - (w1 / edge_length)
      } else {
        curvatures[edge_idx] <- 0
      }
    } else {
      curvatures[edge_idx] <- 0
    }
    }
  }

  return(curvatures)
}

#' Step 2b / Section 4.2 (Delta selection via modularity heuristic)
select_cut_threshold <- function(evolved_lg, probs = seq(0.70, 0.995, by = 0.005)) {
  # Pulls the current Ricci-flow-evolved edge weights of the line graph L^T.
  ws <- E(evolved_lg)$weight
  if (length(ws) == 0) {
    return(0)
  }

  # Builds candidate cut values from quantiles of the weights (e.g. 70% up to 99.5%);
  #idea: bridges are high-weight, so try cutting at high quantiles.
  candidates <- sort(unique(as.numeric(quantile(ws, probs = probs))))

  # Also consider the actual observed weight levels. With the inclusive cut rule (>=),
  # this is sufficient to reproduce "just below a weight" behavior without eps hacks.
  candidates <- sort(unique(c(candidates, ws)))

  # Initializes a table to store, for each candidate cut: the cut value, the modularity score, and number of components.
  scores <- data.frame(cut = numeric(), modularity = numeric(), n_comp = integer())

  # Iterates all candidate thresholds.
  for (cut in candidates) {
    cut_lg <- delete_edges(evolved_lg, which(E(evolved_lg)$weight >= cut))
    memb <- components(cut_lg)$membership

    # Number of components (since membership labels are 1..n_comp).
    n_comp <- max(memb)
    if (n_comp >= 2) {
      # Section 4.2 motivates selecting Delta to maximize modularity of the line-graph partition.
      q <- modularity(evolved_lg, membership = memb, weights = ws)
      scores <- rbind(scores, data.frame(cut = cut, modularity = q, n_comp = n_comp))
    }
  }

  if (nrow(scores) == 0) {
    return(as.numeric(quantile(ws, 0.99)))
  }

  scores <- scores[order(-scores$modularity, scores$n_comp), ]
  return(scores$cut[1])
}


###################

# MUST RUN SENSITIVITY ANALYSIS FOR DELTA_PRUNE

##################

#' Step 2a / Algorithm 4.1 & Eq. 4.2–4.3 (Ricci flow evolution)
ricci_flow <- function(graph, T_iter = 50, delta_prune = 0.1) {
  line_graph <- make_line_graph(graph)
  V(line_graph)$name <- as.character(seq_len(vcount(line_graph)))

  # because we work with unweighed graphs
  E(line_graph)$weight <- 1.0

  for (iteration in seq_len(T_iter)) {
    if (ecount(line_graph) == 0) break
    # Eq. 4.2: compute κ for each edge, then shrink high-curvature edges slower.
    curvature <- compute_orc(line_graph)
    relaxed_weights <- (1 - curvature) * E(line_graph)$weight
    if (!all(is.finite(relaxed_weights)) || sum(relaxed_weights) <= 0) break

    # Normalize according to Eq. 4.3 and update weights.
    E(line_graph)$weight <- (ecount(line_graph) * relaxed_weights) / sum(relaxed_weights)

    # Prune edges whose weight drops below the coarsening threshold.
    to_prune <- which(E(line_graph)$weight < delta_prune)
    if (length(to_prune) > 0 && length(to_prune) < ecount(line_graph)) {
      line_graph <- delete_edges(line_graph, to_prune)
    }
  }

  return(line_graph)
}


###################

# MUST RUN SENSITIVITY ANALYSIS FOR DELTA_CUT

##################

#
#' Step 3 / Section 4.1 Eq. (4.4) (normalize cluster counts by degree)
#' @param delta_cut The threshold to cut bridges (high weight edges) per Section 4.1.
#' If NULL, Section 4.2 modularity heuristic is reused.
infer_membership <- function(graph, evolved_line_graph, delta_cut = NULL) {
  # Ensures original graph vertices have names; later used as rownames.
  if (is.null(V(graph)$name)) {
    V(graph)$name <- as.character(seq_len(vcount(graph)))
  }

  # If user didn’t provide delta_cut, pick it via the modularity heuristic above.
  if (is.null(delta_cut)) {
    delta_cut <- select_cut_threshold(evolved_line_graph)
  }

  # Builds the cut line graph: removes “bridge” edges (high weight) as in the pipeline logic.
  partitioned_line_graph <- delete_edges(
    evolved_line_graph,
    which(E(evolved_line_graph)$weight >= delta_cut)
  )

  # Community label for each vertex of the line graph, i.e. each original edge of graph.
  edge_clusters <- components(partitioned_line_graph)$membership

  num_vertices <- vcount(graph)
  num_clusters <- if (length(edge_clusters) == 0) 0 else max(edge_clusters)
  if (num_clusters == 0) {
    membership_matrix <- matrix(0, num_vertices, 0)
    rownames(membership_matrix) <- V(graph)$name
    return(membership_matrix)
  }

  # Initializes the matrix Y where Y[v, c] will count how many incident edges of v belong to edge-community c.
  membership_matrix <- matrix(0, num_vertices, num_clusters)
  graph_edges <- as_edgelist(graph, names = FALSE)

  # Tries to map each line-graph vertex back to an original edge ID using the vertex name field.
  # (This assumes V(line_graph)$name was set to 1..|E(G)| earlier.)
  edge_to_node <- suppressWarnings(as.integer(V(partitioned_line_graph)$name))

  # If that mapping is inconsistent/missing, fall back to a simple 1..min(...) mapping.
  if (length(edge_to_node) != nrow(graph_edges) || any(is.na(edge_to_node))) {
    edge_to_node <- rep(NA_integer_, vcount(partitioned_line_graph))
    usable <- seq_len(min(vcount(partitioned_line_graph), nrow(graph_edges)))
    edge_to_node[usable] <- usable
  }

  # Iterates each line-graph vertex (each original edge).
  for (line_vertex_idx in seq_len(vcount(partitioned_line_graph))) {
    original_edge_id <- edge_to_node[line_vertex_idx]
    if (!is.na(original_edge_id)) {
      endpoints <- graph_edges[original_edge_id, , drop = FALSE]
      cluster <- edge_clusters[line_vertex_idx]

      membership_matrix[endpoints[1], cluster] <- membership_matrix[endpoints[1], cluster] + 1
      membership_matrix[endpoints[2], cluster] <- membership_matrix[endpoints[2], cluster] + 1
    }
  }

  # Eq. 4.4: normalize the indicator sums by the node degree to obtain mixed-membership scores.
  degrees <- degree(graph)
  degrees[degrees == 0] <- 1
  membership_matrix <- membership_matrix / degrees
  rownames(membership_matrix) <- V(graph)$name
  colnames(membership_matrix) <- paste0("C", seq_len(num_clusters))
  membership_matrix
}

###################

# MUST RUN SENSITIVITY ANALYSIS FOR DELTA_CUT

##################

#' Step 4 / Appendix workflow: run Ricci flow → select threshold → infer membership.
run_curvature_cd <- function(g, T_iter = 50, delta_prune = 0.1, delta_cut = NULL) {
  # Use the paper’s suggested flow: (a) evolve line graph, (b) cut high-weight edges, (c) infer membership.
  evolved_lg <- ricci_flow(g, T_iter = T_iter, delta_prune = delta_prune)

  if (is.null(delta_cut)) {
    delta_cut <- select_cut_threshold(evolved_lg)
  }

  membership <- infer_membership(g, evolved_lg, delta_cut = delta_cut)
  partitioned_lg <- delete_edges(evolved_lg, which(E(evolved_lg)$weight >= delta_cut))
  edge_partition <- components(partitioned_lg)$membership

  list(
    membership = membership,
    edge_partition = edge_partition,
    evolved_lg = evolved_lg,
    delta_cut = delta_cut
  )
}
