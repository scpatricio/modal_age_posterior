# =============================================================================
# 03_figures_simulation.R
# Figures for the Monte Carlo experiment (run 01_simulation.R first).
#
#   Figures/sim_scenarios_N<N>_start<age>.pdf
#     One row per scenario.
#     A: true density of age at death (from the formula), true mode(s).
#     B: posterior distribution of the mode, one replication (bars) and the
#        average over replications (line and points).
#     C: point estimates of every method across replications.
#   Figures/sim_bias.pdf, Figures/sim_mae.pdf
#     Bias and mean absolute error (years) by starting age, method and
#     population size. In scenario (c) the error is measured to the nearest
#     of the two true modes.
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
library(ggplot2)
library(patchwork)

in_dir = file.path("Results", "simulation")
fig_dir = "Figures"
dir.create(fig_dir, showWarnings = FALSE)

tr = readRDS(file.path(in_dir, "sim_truth.rds"))
est = readRDS(file.path(in_dir, "sim_estimates.rds"))
post = readRDS(file.path(in_dir, "sim_posteriors.rds"))

# which cell to show in the scenario figure
show_N = 1e5
show_start = 40
xlim = c(40, 105)

method_labels = c(posterior_mean = "Posterior (mean)",
                  posterior_map = "Posterior (most probable age)",
                  gompertz = "Gompertz", pspline = "P-spline", ddt = "DDT")
method_cols = c("Posterior (mean)" = "#0072B2",
                "Posterior (most probable age)" = "#56B4E9",
                "Gompertz" = "#E69F00", "P-spline" = "#009E73", "DDT" = "#D55E00")
method_shapes = c("Posterior (mean)" = 16, "Posterior (most probable age)" = 1,
                  "Gompertz" = 17, "P-spline" = 15, "DDT" = 4)
scen_labels = setNames(tr$scenarios$label, tr$scenarios$scenario)
est$method_lab = factor(method_labels[est$method], levels = method_labels)

theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  strip.text = element_text(face = "bold", hjust = 0),
                  legend.position = "bottom"))

true_mode_df = function(start_age) {
  do.call(rbind, lapply(tr$truth, function(t) {
    p = t$peaks[t$peaks$x >= start_age, ]
    data.frame(scenario = t$scenario, x = p$x[p$f >= max(p$f) * (1 - 1e-6)])
  }))
}

# ---- Scenario figure --------------------------------------------------------

make_scenario_figure = function(N, start_age) {
  tm = true_mode_df(start_age)
  rows = list()
  for (k in seq_along(tr$truth)) {
    sc = tr$truth[[k]]$scenario
    dens = tr$truth[[k]]$density
    dens = dens[dens$age >= xlim[1] & dens$age <= xlim[2], ]
    tmk = tm[tm$scenario == sc, ]

    pA = ggplot(dens, aes(age, f)) +
      geom_line(linewidth = 0.6) +
      geom_vline(data = tmk, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
      coord_cartesian(xlim = xlim) +
      labs(x = NULL, y = "Density", subtitle = scen_labels[sc])

    pm = post$mean[post$mean$scenario == sc & post$mean$N == N &
                     post$mean$start_age == start_age, ]
    p1 = post$individual[post$individual$scenario == sc & post$individual$N == N &
                           post$individual$start_age == start_age &
                           post$individual$rep == 1, ]
    pB = ggplot() +
      geom_col(data = p1, aes(age + 0.5, prob), fill = "#0072B2", alpha = 0.35,
               width = 0.9) +
      geom_line(data = pm, aes(age + 0.5, prob), linewidth = 0.4) +
      geom_point(data = pm, aes(age + 0.5, prob), size = 0.8) +
      geom_vline(data = tmk, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
      coord_cartesian(xlim = xlim) +
      labs(x = NULL, y = "P(mode = x | data)")

    ek = est[est$scenario == sc & est$N == N & est$start_age == start_age, ]
    pC = ggplot(ek, aes(estimate, method_lab, colour = method_lab, shape = method_lab)) +
      geom_vline(data = tmk, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
      geom_jitter(height = 0.2, width = 0, size = 0.9, alpha = 0.5) +
      scale_colour_manual(values = method_cols, guide = "none") +
      scale_shape_manual(values = method_shapes, guide = "none") +
      scale_y_discrete(limits = rev(levels(est$method_lab))) +
      coord_cartesian(xlim = xlim) +
      labs(x = NULL, y = NULL)

    rows[[k]] = pA | pB | pC
  }
  fig = wrap_plots(rows, ncol = 1) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "", "", "", "", "", "")),
      caption = sprintf("Population size N = %s, starting age %d. Dashed lines: true mode(s). Age at death on the horizontal axis.",
                        format(N, big.mark = ",", scientific = FALSE), start_age))
  fig
}

for (N in unique(est$N)) {
  for (st in c(40, 70)) {
    if (!any(est$N == N & est$start_age == st)) next
    ggsave(file.path(fig_dir, sprintf("sim_scenarios_N%s_start%d.pdf",
                                      format(N, scientific = FALSE), st)),
           make_scenario_figure(N, st), width = 11, height = 8)
  }
}

# ---- Bias and MAE by starting age -------------------------------------------

summ = read.csv(file.path(in_dir, "sim_summary.csv"))
summ$method_lab = factor(method_labels[summ$method], levels = method_labels)
short_labels = c(a = "(a) Unimodal", b = "(b) Lower secondary mode", c = "(c) Equal modes")
summ$scen_lab = factor(short_labels[summ$scenario], levels = short_labels)
summ$N_lab = factor(paste("N =", format(summ$N, big.mark = ",", scientific = FALSE)),
                    levels = paste("N =", format(sort(unique(summ$N)), big.mark = ",",
                                                 scientific = FALSE)))

err_plot = function(var, ylab) {
  ggplot(summ, aes(start_age, .data[[var]], colour = method_lab, shape = method_lab)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(scen_lab ~ N_lab, scales = "free_y") +
    scale_colour_manual(values = method_cols, name = NULL) +
    scale_shape_manual(values = method_shapes, name = NULL) +
    scale_x_continuous(breaks = sort(unique(summ$start_age))) +
    labs(x = "Starting age", y = ylab)
}

ggsave(file.path(fig_dir, "sim_bias.pdf"), err_plot("bias", "Bias (years)"),
       width = 10, height = 8)
ggsave(file.path(fig_dir, "sim_mae.pdf"), err_plot("mae", "Mean absolute error (years)"),
       width = 10, height = 8)
cat("Figures saved to", fig_dir, "\n")
