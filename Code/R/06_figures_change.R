# =============================================================================
# 06_figures_change.R
# Figures for 05_simulation_change.R.
#
# Figures/change/
#   change_detection.pdf   share of replications with P(M2 > M1 | data) > 0.95,
#                          by true change delta, N and overdispersion cv
#                          (delta = 0: false-positive rate, nominal 0.05)
#   change_calibration.pdf distribution of P(M2 > M1 | data) when delta = 0
#                          (uniform when the posterior is calibrated)
#   change_sign.pdf        share of replications in which each method
#                          (posterior mean, median, mode and the three
#                          benchmarks) gets the sign of the change right
#                          (cv = 0), with the posterior decision rule for
#                          reference
#   overdisp_coverage.pdf  coverage of the 95% credible interval by cv and N
#   overdisp_mae.pdf       mean absolute error by cv, method and N
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
library(ggplot2)

in_dir = file.path("Results", "change")
fig_dir = file.path("Figures", "change")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cs = read.csv(file.path(in_dir, "change_summary.csv"))
od = read.csv(file.path(in_dir, "overdisp_summary.csv"))
reps = readRDS(file.path(in_dir, "change_reps.rds"))

method_labels = c(posterior_mean = "Posterior mean", posterior_median = "Posterior median",
                  posterior_mode = "Posterior mode", gompertz = "Gompertz",
                  pspline = "P-spline", ddt = "DDT")
methods = c("posterior_mean", "posterior_median", "posterior_mode",
            "gompertz", "pspline", "ddt")
cols = setNames(c("#0072B2", "#56B4E9", "#CC79A7", "#E69F00", "#009E73", "#D55E00"),
                method_labels[methods])
shapes = setNames(c(16, 1, 18, 17, 15, 4), method_labels[methods])
n_cols = c("#0072B2", "#E69F00", "#009E73", "#D55E00")

fmt_N = function(N) paste("N =", format(N, big.mark = ",", scientific = FALSE))
fmt_cv = function(cv) sprintf("Extra noise CV = %g%%", 100 * cv)
N_levels = fmt_N(sort(unique(cs$N)))
cv_levels = fmt_cv(sort(unique(cs$cv)))

theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  strip.text = element_text(face = "bold", hjust = 0),
                  legend.position = "bottom"))

# ---- Detection rate ---------------------------------------------------------

d = cs[cs$method == "posterior_mean", ]       # detection does not depend on method
d$N_lab = factor(fmt_N(d$N), levels = N_levels)
d$cv_lab = factor(fmt_cv(d$cv), levels = cv_levels)
p = ggplot(d, aes(delta, detect, colour = N_lab, shape = N_lab)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ cv_lab, nrow = 1) +
  scale_colour_manual(values = n_cols, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 4), name = NULL) +
  scale_x_continuous(breaks = sort(unique(d$delta))) +
  labs(x = "True change in the modal age (years)",
       y = "Share with P(M2 > M1 | data) > 0.95")
ggsave(file.path(fig_dir, "change_detection.pdf"), p, width = 11, height = 3.8)

# ---- Calibration under no change --------------------------------------------

r0 = reps[reps$delta == 0 & reps$method == "posterior_mean", ]
r0$N_lab = factor(fmt_N(r0$N), levels = N_levels)
r0$cv_lab = factor(fmt_cv(r0$cv), levels = cv_levels)
p = ggplot(r0, aes(p_increase)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "#0072B2", colour = "white",
                 linewidth = 0.3) +
  facet_grid(N_lab ~ cv_lab) +
  labs(x = "P(M2 > M1 | data) when the true change is 0", y = "Replications")
ggsave(file.path(fig_dir, "change_calibration.pdf"), p, width = 11, height = 8)

# ---- Sign accuracy ----------------------------------------------------------

sg = cs[cs$cv == 0 & cs$delta > 0 & cs$method %in% methods, ]
sg$method_lab = factor(method_labels[sg$method], levels = method_labels[methods])
dec = cs[cs$cv == 0 & cs$delta > 0 & cs$method == "posterior_mean", ]
dec$delta_lab = sprintf("True change = %g years", dec$delta)
sg$delta_lab = sprintf("True change = %g years", sg$delta)
p = ggplot(sg, aes(N, sign_correct, colour = method_lab, shape = method_lab)) +
  geom_line(data = dec, aes(N, detect), inherit.aes = FALSE, linetype = 2,
            colour = "grey30", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ delta_lab, nrow = 1) +
  scale_x_log10(breaks = sort(unique(sg$N)),
                labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Population size N", y = "Share of replications",
       caption = paste("Coloured lines: estimated change has the correct sign.",
                       "Dashed line: posterior decision P(M2 > M1 | data) > 0.95."))
ggsave(file.path(fig_dir, "change_sign.pdf"), p, width = 11, height = 4)

# ---- Overdispersion ---------------------------------------------------------

cv = od[od$method == "posterior_mean", ]
cv$N_lab = factor(fmt_N(cv$N), levels = N_levels)
p = ggplot(cv, aes(100 * cv, coverage, colour = N_lab, shape = N_lab)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(values = n_cols, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 4), name = NULL) +
  scale_x_continuous(breaks = 100 * sort(unique(cv$cv))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Coefficient of variation of the extra noise (%)",
       y = "Coverage of the 95% credible interval")
ggsave(file.path(fig_dir, "overdisp_coverage.pdf"), p, width = 6, height = 4)

om = od[od$method %in% methods, ]
om$method_lab = factor(method_labels[om$method], levels = method_labels[methods])
om$N_lab = factor(fmt_N(om$N), levels = N_levels)
p = ggplot(om, aes(100 * cv, mae, colour = method_lab, shape = method_lab)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ N_lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_x_continuous(breaks = 100 * sort(unique(om$cv))) +
  labs(x = "Coefficient of variation of the extra noise (%)",
       y = "Mean absolute error (years)")
ggsave(file.path(fig_dir, "overdisp_mae.pdf"), p, width = 11, height = 3.8)
cat("Figures saved to", fig_dir, "\n")
