# =============================================================================
# 03_figures_rw.R
# Manuscript figures for the Monte Carlo experiment (run 01_simulation_rw.R
# first). Starting age 40. Every figure is saved to random_walk/Figures/
# simulation/ and to Manuscript/figures/.
#
#   mc_scenarios.pdf   one row per scenario (N = 10^4):
#                      A true density, B posterior of the mode with
#                      independent gamma priors (raw), C posterior of the
#                      mode with the RW2 prior, D point estimates across
#                      replications (box plots)
#   mc_mae.pdf         mean absolute error by population size (log scale)
#   mc_bias.pdf        bias by population size
#   mc_coverage.pdf    coverage and width of the 95% credible interval,
#                      raw and RW2 posteriors
# =============================================================================

code_dir = Sys.getenv("MODE_RW_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk")
setwd(code_dir)
source("../sim_plot_utils.R")      # theme, truth, estimates, summaries, helpers

post_all = readRDS(file.path(in_dir, "sim_posteriors_all.rds"))$mean

fig_dir = file.path("Figures", "simulation")
ms_dir = Sys.getenv("MODE_MS_FIG_DIR", file.path("..", "..", "..", "Manuscript", "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)
save_both = function(p, name, width, height) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height)
  ggsave(file.path(ms_dir, name), p, width = width, height = height)
}

start = 40
show_N = 1e4

method_labels = c(rw2_mean = "RW2 posterior mean", rw2_median = "RW2 posterior median",
                  rw2_mode = "RW2 posterior mode", raw_mean = "Raw posterior mean",
                  raw_median = "Raw posterior median", raw_mode = "Raw posterior mode",
                  gompertz = "Gompertz", pspline = "P-spline", ddt = "DDT")
methods = c("rw2_mean", "rw2_median", "rw2_mode", "raw_mean", "gompertz", "pspline", "ddt")
cols = setNames(c("#0072B2", "#56B4E9", "#CC79A7", "#999999", "#E69F00", "#009E73",
                  "#D55E00"), method_labels[methods])
shapes = setNames(c(16, 1, 18, 3, 17, 15, 4), method_labels[methods])

# ---- Scenario figure --------------------------------------------------------

scenario_row_rw = function(sc, first) {
  tm = data.frame(x = true_modes_sc(sc, start))
  op = data.frame(x = other_peaks_sc(sc, start))
  pm = post_all[post_all$scenario == sc & post_all$N == show_N &
                  post_all$start_age == start, ]
  ek = est[est$scenario == sc & est$N == show_N & est$start_age == start &
             est$method %in% methods & is.finite(est$estimate), ]

  # common age window for B, C, D
  a = pm$age[pm$prob >= 0.001]
  win = c(min(a, floor(min(ek$estimate))) - 1, max(a + 1, ceiling(max(ek$estimate))) + 1)

  dens = tr$truth[[sc]]$density
  dens = dens[dens$age >= start & dens$age <= 105, ]
  pA = ggplot(dens, aes(age, f)) +
    annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5) +
    geom_line(linewidth = 0.6) +
    geom_vline(data = tm, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
    geom_vline(data = op, aes(xintercept = x), linetype = 3, linewidth = 0.3,
               colour = "grey40") +
    labs(x = NULL, y = "Density", subtitle = scen_labels[sc])

  post_panel = function(fam, title) {
    d = pm[pm$family == fam & pm$age + 0.5 >= win[1] & pm$age + 0.5 <= win[2], ]
    ggplot(d, aes(age + 0.5, prob)) +
      geom_col(fill = if (fam == "rw2") "#0072B2" else "#999999", width = 0.9) +
      geom_vline(data = tm, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
      coord_cartesian(xlim = win) +
      labs(x = NULL, y = "P(mode = x | data)", subtitle = if (first) title else NULL)
  }
  pB = post_panel("raw", "Raw posterior")
  pC = post_panel("rw2", "RW2 posterior")

  ek$method_lab = factor(method_labels[ek$method], levels = method_labels[methods])
  pD = ggplot(ek, aes(estimate, method_lab, colour = method_lab, shape = method_lab)) +
    geom_vline(data = tm, aes(xintercept = x), linetype = 2, linewidth = 0.3) +
    geom_jitter(height = 0.2, width = 0, size = 0.8, alpha = 0.4) +
    geom_boxplot(aes(group = method_lab), fill = NA, colour = "grey20",
                 outlier.shape = NA, width = 0.5, linewidth = 0.3) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_shape_manual(values = shapes, guide = "none") +
    scale_y_discrete(limits = rev(levels(ek$method_lab))) +
    coord_cartesian(xlim = win) +
    labs(x = NULL, y = NULL, subtitle = if (first) "Point estimates" else NULL)
  list(pA, pB, pC, pD)
}

scs = tr$scenarios$scenario
plots = unlist(lapply(seq_along(scs), function(k) scenario_row_rw(scs[k], k == 1)),
               recursive = FALSE)
fig = wrap_plots(plots, ncol = 4, widths = c(1, 1, 1, 1.3)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", rep("", 4 * (length(scs) - 1)))),
    caption = sprintf(paste0("%s, starting age %d, %d replications. A: true density ",
                             "(shaded: ages shown in B-D). B, C: posterior of the mode averaged ",
                             "over replications.\nD: point estimates across replications. ",
                             "Dashed line: true mode; dotted line in A: secondary peak. ",
                             "Age at death on the horizontal axes."),
                      fmt_N(show_N), start, tr$settings$n_rep),
    theme = theme(plot.caption = element_text(hjust = 0)))
save_both(fig, "mc_scenarios.pdf", 13, 2.7 * length(scs) + 0.8)

# ---- Errors by population size ----------------------------------------------

save_both(error_plot_by_N("mae", "Mean absolute error (years)", methods, cols, shapes,
                          start) + scale_y_log10(),
          "mc_mae.pdf", 11, 4)
save_both(error_plot_by_N("bias", "Bias (years)", methods[methods != "gompertz"],
                          cols, shapes, start) +
            labs(caption = sprintf("Starting age %d. Gompertz not shown (bias of 3 to 7 years with a secondary peak).", start)),
          "mc_bias.pdf", 11, 4)

# ---- Coverage and width of the 95% credible interval ------------------------

fam_lab = c(raw_mean = "Raw posterior", rw2_mean = "RW2 posterior")
cv = summ[summ$method %in% names(fam_lab) & summ$start_age == start, ]
cv$fam = factor(fam_lab[cv$method], levels = fam_lab)
cv$scen_lab = factor(short_labels[cv$scenario], levels = short_labels)
base = function(var) {
  ggplot(cv, aes(N, .data[[var]], colour = fam, shape = fam)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_wrap(~ scen_lab, nrow = 1) +
    scale_x_log10(breaks = sort(unique(cv$N)),
                  labels = function(x) parse(text = paste0("10^", log10(x)))) +
    scale_colour_manual(values = setNames(c("#999999", "#0072B2"), fam_lab), name = NULL) +
    scale_shape_manual(values = setNames(c(3, 16), fam_lab), name = NULL)
}
p1 = base("coverage") +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50", linewidth = 0.3) +
  coord_cartesian(ylim = c(min(0.8, min(cv$coverage)), 1)) +
  labs(x = NULL, y = "Coverage")
p2 = base("mean_width") +
  scale_y_log10() +
  labs(x = "Population size N", y = "Width (years)")
save_both(p1 / p2 + plot_layout(guides = "collect") &
            theme(legend.position = "bottom"),
          "mc_coverage.pdf", 11, 6)
cat("Figures saved to", fig_dir, "and", ms_dir, "\n")
