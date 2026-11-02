# =============================================================================
# sim_plot_utils.R
# Plotting helpers shared by 03_figures_posterior.R and
# 04_figures_comparison.R. Reads the simulation output on source.
# =============================================================================

library(ggplot2)
library(patchwork)

in_dir = file.path("Results", "simulation")
tr = readRDS(file.path(in_dir, "sim_truth.rds"))
est = readRDS(file.path(in_dir, "sim_estimates.rds"))
post = readRDS(file.path(in_dir, "sim_posteriors.rds"))
summ = read.csv(file.path(in_dir, "sim_summary.csv"))

method_labels = c(posterior_mean = "Posterior mean",
                  posterior_median = "Posterior median",
                  posterior_mode = "Posterior mode",
                  gompertz = "Gompertz", pspline = "P-spline", ddt = "DDT")
scen_labels = setNames(tr$scenarios$label, tr$scenarios$scenario)
short_labels = setNames(ifelse(tr$scenarios$w == 1,
                              sprintf("(%s) Unimodal", tr$scenarios$scenario),
                              sprintf("(%s) %.0f%% in earlier component", tr$scenarios$scenario,
                                      100 * (1 - tr$scenarios$w))),
                       tr$scenarios$scenario)

fmt_N = function(N) paste("N =", format(N, big.mark = ",", scientific = FALSE))

theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  strip.text = element_text(face = "bold", hjust = 0),
                  legend.position = "bottom"))

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

# age window holding the posterior: ages with mean posterior >= min_prob
post_window = function(sc, N, start_age, min_prob = 0.001) {
  pm = post$mean[post$mean$scenario == sc & post$mean$N == N &
                   post$mean$start_age == start_age, ]
  a = pm$age[pm$prob >= min_prob]
  c(min(a) - 1, max(a) + 2)
}

# One row of the scenario figure.
# A: true density, with the window of B and C shaded
# B: posterior of the mode in that window, averaged over replications
#    (a discrete distribution over single-year ages, shown as bars)
# C: point estimates across replications for the chosen methods
scenario_row = function(sc, N, start_age, methods, cols, shapes) {
  tm = data.frame(x = true_modes_sc(sc, start_age))
  op = data.frame(x = other_peaks_sc(sc, start_age))
  win = post_window(sc, N, start_age)
  # widen the window so every point estimate in panel C is visible
  ek = est[est$scenario == sc & est$N == N & est$start_age == start_age &
             est$method %in% methods, ]
  ek = ek[is.finite(ek$estimate), ]
  win = c(min(win[1], floor(min(ek$estimate))), max(win[2], ceiling(max(ek$estimate))))

  dens = tr$truth[[sc]]$density
  dens = dens[dens$age >= start_age & dens$age <= 105, ]
  pA = ggplot(dens, aes(age, f)) +
    annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_line(linewidth = 0.6) +
    geom_vline(data = tm, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
    geom_vline(data = op, aes(xintercept = x), linetype = 3, linewidth = 0.3,
               colour = "grey40") +
    labs(x = NULL, y = "Density", subtitle = scen_labels[sc])

  pm = post$mean[post$mean$scenario == sc & post$mean$N == N &
                   post$mean$start_age == start_age, ]
  pm = pm[pm$age + 0.5 >= win[1] & pm$age + 0.5 <= win[2], ]
  # bar for age x covers the interval [x, x + 1)
  pB = ggplot() +
    geom_col(data = pm, aes(age + 0.5, prob), fill = "#0072B2", width = 0.9) +
    geom_vline(data = tm, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
    coord_cartesian(xlim = win) +
    labs(x = NULL, y = "P(mode = x | data)")

  ek$method_lab = factor(method_labels[ek$method], levels = method_labels[methods])
  pC = ggplot(ek, aes(estimate, method_lab, colour = method_lab, shape = method_lab)) +
    geom_vline(data = tm, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
    geom_jitter(height = 0.2, width = 0, size = 1.2, alpha = 0.6) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_shape_manual(values = shapes, guide = "none") +
    scale_y_discrete(limits = rev(levels(ek$method_lab))) +
    coord_cartesian(xlim = win) +
    labs(x = NULL, y = NULL)
  pC = pC + geom_boxplot(aes(group = method_lab), fill = NA, colour = "grey20",
                         outlier.shape = NA, width = 0.5, linewidth = 0.3)

  list(A = pA, B = pB, C = pC)
}

# Whole figure: one row per scenario run at this starting age
scenario_figure = function(N, start_age, methods, cols, shapes) {
  scs = tr$scenarios$scenario[sapply(tr$scenarios$scenario, function(sc) {
    any(est$scenario == sc & est$N == N & est$start_age == start_age)
  })]
  rows = lapply(scs, function(sc) scenario_row(sc, N, start_age, methods, cols, shapes))
  plots = unlist(lapply(rows, function(r) list(r$A, r$B, r$C)), recursive = FALSE)
  tags = c("A", "B", "C", rep("", 3 * (length(scs) - 1)))
  wrap_plots(plots, ncol = 3, widths = c(1, 1, 1.2)) +
    plot_annotation(
      tag_levels = list(tags),
      caption = sprintf(paste0("%s, starting age %d. A: true density (shaded: ages shown in B and C). ",
                               "B: posterior of the mode, averaged over replications.\n",
                               "C: point estimates across replications. ",
                               "Dashed line: true (global) mode; dotted line in A: secondary peak. ",
                               "Age at death on the horizontal axes."),
                        fmt_N(N), start_age),
      theme = theme(plot.caption = element_text(hjust = 0)))
}

# Ns, starts: population sizes and starting ages to draw (NULL = all)
save_scenario_figures = function(methods, cols, shapes, out_dir, prefix,
                                 Ns = NULL, starts = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(Ns)) Ns = sort(unique(est$N))
  if (is.null(starts)) starts = sort(unique(est$start_age))
  for (N in Ns) {
    for (st in starts) {
      n_rows = sum(sapply(tr$scenarios$scenario, function(sc) {
        any(est$scenario == sc & est$N == N & est$start_age == st)
      }))
      ggsave(file.path(out_dir, sprintf("%s_N%s_start%d.pdf", prefix,
                                        format(N, scientific = FALSE), st)),
             scenario_figure(N, st, methods, cols, shapes),
             width = 11, height = 2.7 * n_rows + 0.6)
    }
  }
}

# Bias, error or SD by starting age
error_plot = function(var, ylab, methods, cols, shapes,
                      scenarios = tr$scenarios$scenario) {
  d = summ[summ$method %in% methods & summ$scenario %in% scenarios, ]
  d$method_lab = factor(method_labels[d$method], levels = method_labels[methods])
  d$scen_lab = factor(short_labels[d$scenario], levels = short_labels)
  d$N_lab = factor(fmt_N(d$N), levels = fmt_N(sort(unique(d$N))))
  p = ggplot(d, aes(start_age, .data[[var]], colour = method_lab, shape = method_lab)) +
    geom_hline(yintercept = if (var == "bias") 0 else NA_real_, colour = "grey60",
               linewidth = 0.3, na.rm = TRUE) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(scen_lab ~ N_lab, scales = "free_y") +
    scale_colour_manual(values = cols, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    scale_x_continuous(breaks = sort(unique(d$start_age))) +
    labs(x = "Starting age", y = ylab)
  p
}

# Bias, error or SD by population size, at one starting age
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
    labs(x = "Population size N", y = ylab,
         caption = sprintf("Starting age %d.", start_age))
}
