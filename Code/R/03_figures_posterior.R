# =============================================================================
# 03_figures_posterior.R
# Stage 1: which summary of the posterior should be the point estimate of
# the mode? Compares the posterior mean, median and mode at one starting age
# (starting ages are compared only in stage 2, against the other methods).
# Run 01_simulation.R first.
#
# Figures/posterior/
#   post_scenarios_N<N>_start<age>.pdf   panels A, B, C (C: posterior mean,
#                                        median and mode across replications)
#   post_bias.pdf, post_mae.pdf, post_sd.pdf
#                                        by population size and scenario
#   post_coverage.pdf                    coverage of the 95% credible interval
# Results/simulation/
#   posterior_estimator_ranking.csv      errors by estimator and N
#   best_posterior_estimator.txt         estimator with the lowest MAE,
#                                        read by 04_figures_comparison.R
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
source("sim_plot_utils.R")

stage1_start = 30      # starting age used in stage 1
stage1_N = 1e4         # population size for the A/B/C figure

fig_dir = file.path("Figures", "posterior")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

methods = c("posterior_mean", "posterior_median", "posterior_mode")
cols = setNames(c("#0072B2", "#E69F00", "#009E73"), method_labels[methods])
shapes = setNames(c(16, 17, 15), method_labels[methods])

# ---- Scenario figure --------------------------------------------------------

save_scenario_figures(methods, cols, shapes, fig_dir, "post_scenarios",
                      Ns = stage1_N, starts = stage1_start)

# ---- Bias, MAE, SD by population size ---------------------------------------

ggsave(file.path(fig_dir, "post_bias.pdf"),
       error_plot_by_N("bias", "Bias (years)", methods, cols, shapes, stage1_start),
       width = 11, height = 3.8)
ggsave(file.path(fig_dir, "post_mae.pdf"),
       error_plot_by_N("mae", "Mean absolute error (years)", methods, cols, shapes,
                       stage1_start),
       width = 11, height = 3.8)
ggsave(file.path(fig_dir, "post_sd.pdf"),
       error_plot_by_N("sd", "SD of the estimates (years)", methods, cols, shapes,
                       stage1_start),
       width = 11, height = 3.8)

# ---- Coverage of the 95% credible interval ----------------------------------

cv = summ[summ$method == "posterior_mean" & summ$start_age == stage1_start, ]
cv$scen_lab = factor(short_labels[cv$scenario], levels = short_labels)
p_cov = ggplot(cv, aes(N, coverage, colour = scen_lab, shape = scen_lab)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_log10(breaks = sort(unique(cv$N)),
                labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73", "#D55E00"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 4), name = NULL) +
  coord_cartesian(ylim = c(min(0.8, min(cv$coverage)), 1)) +
  labs(x = "Population size N", y = "Coverage of the 95% credible interval",
       caption = sprintf("Starting age %d.", stage1_start))
ggsave(file.path(fig_dir, "post_coverage.pdf"), p_cov, width = 7, height = 4)

# ---- Ranking ----------------------------------------------------------------

r = summ[summ$method %in% methods & summ$start_age == stage1_start, ]
rank_N = aggregate(cbind(mae, rmse, abs_bias = abs(bias), sd) ~ method + N, data = r, FUN = mean)
rank_all = aggregate(cbind(mae, rmse, abs_bias = abs(bias), sd) ~ method, data = r, FUN = mean)
rank_all$N = "all"
ranking = rbind(rank_N, rank_all[, names(rank_N)])
write.csv(ranking, file.path(in_dir, "posterior_estimator_ranking.csv"), row.names = FALSE)

best = rank_all$method[which.min(rank_all$mae)]
writeLines(best, file.path(in_dir, "best_posterior_estimator.txt"))
print(ranking, row.names = FALSE)
cat("Lowest mean absolute error overall:", method_labels[best], "\n")
cat("Figures saved to", fig_dir, "\n")
