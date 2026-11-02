# =============================================================================
# plot_utils.R
# Plotting helpers for the Monte Carlo figures (02_figures_simulation.R).
# Reads the simulation output in results/simulation/ when sourced.
# =============================================================================

library(ggplot2)
library(patchwork)

in_dir = file.path("results", "simulation")
tr = readRDS(file.path(in_dir, "sim_truth.rds"))
est = readRDS(file.path(in_dir, "sim_estimates.rds"))
summ = read.csv(file.path(in_dir, "sim_summary.csv"))

scen_labels = setNames(tr$scenarios$label, tr$scenarios$scenario)
short_labels = setNames(ifelse(tr$scenarios$w == 1,
                              sprintf("(%s) Unimodal", tr$scenarios$scenario),
                              sprintf("(%s) %.0f%% earlier", tr$scenarios$scenario,
                                      100 * (1 - tr$scenarios$w))),
                       tr$scenarios$scenario)

fmt_N = function(N) paste("N =", format(N, big.mark = ",", scientific = FALSE))

theme_set(theme_mode())

# true mode of one scenario on [start_age, 110], and the other peak(s)
true_modes_sc = function(sc, start_age) {
  p = tr$truth[[sc]]$peaks
  p = p[p$x >= start_age, ]
  p$x[which.max(p$f)]
}
other_peaks_sc = function(sc, start_age) {
  p = tr$truth[[sc]]$peaks
  p = p[p$x >= start_age, ]
  p$x[-which.max(p$f)]
}

# Bias, error or SD by population size, at one starting age. Uses the
# method_labels vector defined by the calling script.
error_plot_by_N = function(var, ylab, methods, cols, shapes, start_age,
                           scenarios = tr$scenarios$scenario) {
  d = summ[summ$method %in% methods & summ$scenario %in% scenarios &
             summ$start_age == start_age, ]
  d$method_lab = factor(method_labels[d$method], levels = method_labels[methods])
  d$scen_lab = factor(short_labels[d$scenario], levels = short_labels)
  ggplot(d, aes(N, .data[[var]], colour = method_lab, shape = method_lab)) +
    geom_hline(yintercept = if (var == "bias") 0 else NA_real_, colour = "grey60",
               linewidth = 0.3, na.rm = TRUE) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_wrap(~ scen_lab, nrow = 1, scales = "free_y") +
    scale_x_log10(breaks = sort(unique(d$N)),
                  labels = function(x) parse(text = paste0("10^", log10(x)))) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    labs(x = "Population size N", y = ylab)
}
