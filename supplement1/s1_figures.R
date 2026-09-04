#!/usr/bin/env Rscript
# =============================================================================
# s1_figures.R — Supplement S1 figures and summary tables from results/s1_results.csv
#
#   fig_s1_tracking.pdf      estimated (corrected) vs true score,
#                            separability x bridge strength (balanced 10+10)
#   fig_s1_separability.pdf  accuracy vs p_between by node type (thesis Fig. 6 analog)
#   fig_s1_invariance.pdf    bridge and pure-node scores across splits,
#                            raw vs corrected, at baseline separability
#   fig_s1_floor.pdf         whole-count floor at the base condition
#   s1_summary.csv           MSE / bias by condition and node type
#   s1_correlations.csv      Spearman(est, true) and mean detected k by condition
# =============================================================================

library(ggplot2)

res <- read.csv(file.path("results", "s1_results.csv"), stringsAsFactors = FALSE)
res$split <- factor(res$split, levels = c("10+10", "13+7", "15+5"))
res$pi_lab <- sprintf("pi == '(%.1f, %.1f)'", res$pi1, 1 - res$pi1)
res$sep_lab <- sprintf("p[b] == %.2f", res$p_between)
res$type_lab <- factor(res$type, levels = c("bridge", "neighbor", "pure"),
                       labels = c("bridge", "neighbor of bridge", "other nodes"))
BASE_PB <- 0.05

dir.create("results", showWarnings = FALSE)

## Figure 1: tracking across separability x strength (balanced networks) -------
bal <- res[res$split == "10+10", ]
p1 <- ggplot(bal, aes(b_true, b_corrected, colour = type_lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(alpha = 0.2, size = 0.6) +
  facet_grid(sep_lab ~ pi_lab, labeller = label_parsed) +
  scale_colour_manual(values = c("#D5A021", "#4477AA", "grey40")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = expression("true score" ~ b[true]),
       y = "estimated score (corrected), known network", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path("results", "fig_s1_tracking.pdf"), p1, width = 6.5, height = 10)

## Figure 2: accuracy vs separability (thesis Fig. 6 analog, known networks) ---
acc <- do.call(rbind, lapply(
  split(bal, list(bal$p_between, bal$pi1, bal$type_lab), drop = TRUE),
  function(d) data.frame(p_between = d$p_between[1], pi1 = d$pi1[1],
                         type_lab = d$type_lab[1],
                         mse_corrected = mean((d$b_corrected - d$b_true)^2),
                         mse_raw = mean((d$b_raw - d$b_true)^2))))
p2 <- ggplot(acc, aes(p_between, mse_corrected, colour = type_lab)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
  facet_wrap(~sprintf("pi == '(%.1f, %.1f)'", pi1, 1 - pi1),
             labeller = label_parsed, nrow = 1) +
  scale_colour_manual(values = c("#D5A021", "#4477AA", "grey40")) +
  labs(x = expression("between-community edge probability" ~ p[between]),
       y = expression("MSE of corrected score vs" ~ b[true]), colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path("results", "fig_s1_separability.pdf"), p2, width = 7.5, height = 3.6)

## Figure 3: size invariance at baseline separability, raw vs corrected --------
inv_dat <- res[res$p_between == BASE_PB, ]
inv <- do.call(rbind, lapply(
  split(inv_dat, list(inv_dat$pi1, inv_dat$split, inv_dat$type_lab,
                      inv_dat$community), drop = TRUE),
  function(d) data.frame(pi1 = d$pi1[1], split = d$split[1],
                         type_lab = d$type_lab[1], community = d$community[1],
                         raw = mean(d$b_raw), corrected = mean(d$b_corrected),
                         truth = mean(d$b_true),
                         raw_sd = sd(d$b_raw), corrected_sd = sd(d$b_corrected))))
inv_long <- rbind(
  transform(inv, score = raw, sd = raw_sd, estimator = "raw"),
  transform(inv, score = corrected, sd = corrected_sd, estimator = "corrected")
)
inv_long$grp <- ifelse(inv_long$type_lab == "bridge", "bridge",
                       paste0("pure, ", ifelse(inv_long$community == "C1",
                                               "larger community", "smaller community")))
inv_long <- inv_long[inv_long$type_lab != "neighbor of bridge", ]
p3 <- ggplot(inv_long, aes(split, score, colour = estimator, group = estimator)) +
  geom_line(linewidth = 0.5) +
  geom_pointrange(aes(ymin = score - sd, ymax = score + sd),
                  size = 0.3, position = position_dodge(width = 0.15)) +
  geom_point(aes(y = truth), colour = "black", shape = 4, size = 2) +
  facet_grid(grp ~ sprintf("pi == '(%.1f, %.1f)'", pi1, 1 - pi1),
             labeller = labeller(.cols = label_parsed)) +
  scale_colour_manual(values = c(corrected = "#4477AA", raw = "#CC6677")) +
  labs(x = "community split", y = "mean score (± SD); x = mean true score",
       colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path("results", "fig_s1_invariance.pdf"), p3, width = 7, height = 6.5)

## Figure 4: whole-count floor at the base condition ---------------------------
fl <- res[!is.na(res$floor_gap), ]
p4 <- ggplot(fl, aes(floor_gap, fill = type_lab)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2) +
  facet_wrap(~type_lab, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#D5A021", "#4477AA", "grey40"), guide = "none") +
  labs(x = expression("smallest achievable" ~ "|" * b - b[true] * "|" ~
                        "with whole counts"),
       y = "nodes") +
  theme_minimal(base_size = 10)
ggsave(file.path("results", "fig_s1_floor.pdf"), p4, width = 5.5, height = 5)

## Summary tables --------------------------------------------------------------
cell_key <- list(res$pi1, res$split, res$p_between, res$type_lab)
summ <- do.call(rbind, lapply(split(res, cell_key, drop = TRUE), function(d) {
  data.frame(pi1 = d$pi1[1], split = as.character(d$split[1]),
             p_within = d$p_within[1], p_between = d$p_between[1],
             type = as.character(d$type_lab[1]), n = nrow(d),
             mse_corrected = mean((d$b_corrected - d$b_true)^2),
             mse_raw = mean((d$b_raw - d$b_true)^2),
             bias_corrected = mean(d$b_corrected - d$b_true),
             bias_raw = mean(d$b_raw - d$b_true),
             mean_true = mean(d$b_true), stringsAsFactors = FALSE)
}))
spear <- do.call(rbind, lapply(
  split(res, list(res$pi1, res$split, res$p_between), drop = TRUE),
  function(d) data.frame(pi1 = d$pi1[1], split = as.character(d$split[1]),
                         p_between = d$p_between[1],
                         spearman_corrected = suppressWarnings(
                           cor(d$b_corrected, d$b_true, method = "spearman")),
                         spearman_raw = suppressWarnings(
                           cor(d$b_raw, d$b_true, method = "spearman")),
                         mean_k_detected = mean(d$k_detected),
                         share_k2 = mean(d$k_detected == 2),
                         stringsAsFactors = FALSE)))
write.csv(summ, file.path("results", "s1_summary.csv"), row.names = FALSE)
write.csv(spear, file.path("results", "s1_correlations.csv"), row.names = FALSE)

fl_summ <- c(median = median(fl$floor_gap), q95 = quantile(fl$floor_gap, 0.95),
             max = max(fl$floor_gap))
message("floor gap: median ", signif(fl_summ[1], 3),
        ", 95th pct ", signif(fl_summ[2], 3), ", max ", signif(fl_summ[3], 3))
message("figures + tables written to results/")
