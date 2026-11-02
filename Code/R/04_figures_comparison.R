# =============================================================================
# 04_figures_comparison.R
# Stage 2: the posterior mean, median and mode against Gompertz, P-spline
# and DDT, at every population size and starting age. Run 01_simulation.R
# first.
#
# Figures/comparison/
#   comp_scenarios_N<N>_start<age>.pdf   panels A, B, C (C: the three
#                                        posterior estimates and the three
#                                        benchmarks)
#   comp_bias.pdf, comp_mae.pdf, comp_sd.pdf
#                                        by starting age and scenario
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
source("sim_plot_utils.R")

fig_dir = file.path("Figures", "comparison")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

methods = c("posterior_mean", "posterior_median", "posterior_mode",
            "gompertz", "pspline", "ddt")
cols = setNames(c("#0072B2", "#56B4E9", "#CC79A7", "#E69F00", "#009E73", "#D55E00"),
                method_labels[methods])
shapes = setNames(c(16, 1, 18, 17, 15, 4), method_labels[methods])

save_scenario_figures(methods, cols, shapes, fig_dir, "comp_scenarios")

ggsave(file.path(fig_dir, "comp_bias.pdf"),
       error_plot("bias", "Bias (years)", methods, cols, shapes), width = 11, height = 9)
ggsave(file.path(fig_dir, "comp_mae.pdf"),
       error_plot("mae", "Mean absolute error (years)", methods, cols, shapes),
       width = 11, height = 9)
ggsave(file.path(fig_dir, "comp_sd.pdf"),
       error_plot("sd", "SD of the estimates (years)", methods, cols, shapes),
       width = 11, height = 9)
cat("Figures saved to", fig_dir, "\n")
