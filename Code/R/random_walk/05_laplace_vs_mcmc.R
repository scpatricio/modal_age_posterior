# =============================================================================
# 05_laplace_vs_mcmc.R
# Checks the Laplace approximation used by mode_posterior_rw() against exact
# posterior sampling of the same RW2 model, for a few HMD populations at
# starting age 50: a small (Denmark) and a large (Japan) population, women
# and men, an early and a recent year.
#
# Exact sampling: independence Metropolis-Hastings (rw_imh_cpp() in
# rw_core.cpp) with the Laplace mixture over the grid of lambda as proposal.
# The accept/reject step corrects any error of the approximation, so the
# chain targets the exact posterior of the model (Poisson likelihood, RW2
# prior, flat prior on the grid of log10 lambda). Its acceptance rate is a
# direct measure of the quality of the approximation (1 if it were exact).
#
# Outputs:
#   Results/mcmc_check/laplace_vs_mcmc.csv     posterior summaries, both methods
#   Manuscript/tables/app_laplace_mcmc.tex      the same as a LaTeX table
#   Figures/mcmc_check/app_laplace_mcmc.pdf     P(mode = x | data), both methods
#   (the figure is also saved to Manuscript/figures/)
#
# Run from Code/R/random_walk.
# =============================================================================

code_dir = Sys.getenv("MODE_RW_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk")
setwd(code_dir)
for (f in c("../utils.R", "../mode_posterior.R", "mode_posterior_rw.R")) source(f)
library(ggplot2)
if (!rw_use_cpp) stop("This check needs the C++ core (rw_core.cpp); see mode_posterior_rw.R.")

data_dir = file.path("..", "..", "data")
start_age = 50
max_age = 110
S = 20000                       # Laplace draws
n_iter = 100000                 # Metropolis-Hastings iterations
log10_lambdas = seq(-1, 8, by = 0.25)
drop = 12
cases = data.frame(country = c("Denmark", "Denmark", "Japan", "Japan",
                               "Denmark", "Denmark", "Japan", "Japan"),
                   sex = rep(c("Female", "Male"), each = 4),
                   year = rep(c(1960, 2019), 4))
seed = 20261108

out_dir = file.path("Results", "mcmc_check")
fig_dir = file.path("Figures", "mcmc_check")
ms_fig = Sys.getenv("MODE_MS_FIG_DIR", file.path("..", "..", "..", "Manuscript", "figures"))
ms_tab = Sys.getenv("MODE_MS_TAB_DIR", file.path("..", "..", "..", "Manuscript", "tables"))
for (d in c(out_dir, fig_dir, ms_fig, ms_tab)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

hmd = do.call(rbind, lapply(unique(cases$country), load_country, data_dir = data_dir))

# life-table deaths from draws of the log rates (draws in rows); last column
# the open interval, as in mode_posterior_rw()
d_from_eta = function(eta) {
  d = cbind(-expm1(-exp(eta)), 1)
  l = rep(1, nrow(d))
  for (i in 1:ncol(d)) {
    qi = d[, i]
    d[, i] = l * qi
    l = l * (1 - qi)
  }
  d
}

# effective sample size of a chain (initial positive sequence estimator)
ess = function(x) {
  n = length(x)
  r = acf(x, lag.max = min(n - 1, 2000), plot = FALSE)$acf[-1]
  k = which(r[seq(1, length(r) - 1, by = 2)] + r[seq(2, length(r), by = 2)] < 0)[1]
  m = if (is.na(k)) length(r) else 2 * (k - 1)
  n / (1 + 2 * sum(r[seq_len(m)]))
}

summ = list()
probs = list()
for (k in seq_len(nrow(cases))) {
  cs = cases[k, ]
  y = hmd[hmd$country == cs$country & hmd$sex == cs$sex & hmd$year == cs$year, ]
  s = prep_schedule(y$age, y$deaths, y$exposure, start_age, max_age)
  nc = length(s$age) - 1
  D = s$deaths[1:nc]
  E = s$exposure[1:nc]
  set.seed(seed + k)

  # Laplace approximation (the method used in the paper)
  lap = mode_posterior_rw(y$age, y$deaths, y$exposure, start_age, max_age, S = S,
                          log10_lambdas = log10_lambdas)

  # exact posterior: independence Metropolis-Hastings with the same proposal
  g = rw_grid_cpp(D, E, 2L, log10_lambdas, drop)
  w = exp(g$logml - max(g$logml))
  w = w / sum(w)
  ch = rw_imh_cpp(D, E, g$eta, g$L, g$log10_lambda, w, n_iter, 2L)
  mc = posterior_summary(d_from_eta(ch$eta), s$age, start_age, n_iter, 0.95, 0.1, 0.5,
                         return_draws = TRUE)
  ll_chain = g$log10_lambda[ch$cell + 1]

  lab = sprintf("%s, %s, %d", cs$country, tolower(cs$sex), cs$year)
  for (meth in c("Laplace", "MCMC")) {
    p = if (meth == "Laplace") lap else mc
    summ[[paste(k, meth)]] = data.frame(cs, case = lab, method = meth,
      mean = p$mean, median = p$median, q025 = p$q025, q975 = p$q975, map = p$map,
      p_map = p$p_map,
      log10_lambda = if (meth == "Laplace") lap$post_log10_lambda else mean(ll_chain),
      acceptance = if (meth == "MCMC") ch$acceptance else NA,
      ess_mode = if (meth == "MCMC") ess(mc$draws_cont) else NA)
    probs[[paste(k, meth)]] = data.frame(case = lab, method = meth, age = p$post$age,
                                         prob = p$post$prob)
  }
  cat(lab, ": acceptance", round(ch$acceptance, 3), "\n")
}
summ = do.call(rbind, summ)
probs = do.call(rbind, probs)
write.csv(summ, file.path(out_dir, "laplace_vs_mcmc.csv"), row.names = FALSE)
print(summ[, c("case", "method", "mean", "q025", "q975", "map", "p_map", "log10_lambda",
               "acceptance", "ess_mode")], row.names = FALSE)

L = summ[summ$method == "Laplace", ]
M = summ[summ$method == "MCMC", ]
cat("\nLargest differences, Laplace vs MCMC:",
    "\n  posterior mean of the mode:", max(abs(L$mean - M$mean)),
    "\n  2.5% and 97.5% limits:", max(abs(L$q025 - M$q025)), max(abs(L$q975 - M$q975)),
    "\n  probability of the most probable age:", max(abs(L$p_map - M$p_map)),
    "\n  same most probable age in all cases:", all(L$map == M$map),
    "\n  posterior mean of log10 lambda:", max(abs(L$log10_lambda - M$log10_lambda)),
    "\n  acceptance rate, range:", range(M$acceptance),
    "\n  effective sample size of the mode, min:", min(M$ess_mode), "\n")

# ---- Table for the appendix -------------------------------------------------

f2 = function(x) formatC(x, format = "f", digits = 2)
rows = sprintf("%s & %s & %s (%s, %s) & %s (%s, %s) & %s \\\\",
               L$country, paste0(tolower(L$sex), ", ", L$year),
               f2(L$mean), f2(L$q025), f2(L$q975), f2(M$mean), f2(M$q025), f2(M$q975),
               f2(M$acceptance))
tab = c("\\begin{table}[htb!]",
        "\\centering",
        "\\caption{RW2 posterior of the modal age at death computed with the Laplace approximation and with exact sampling (independence Metropolis--Hastings): posterior mean of $\\tilde{M}$ with 95\\% credible interval, and acceptance rate of the sampler. Starting age 50.}",
        "\\label{tab:app_mcmc}",
        "\\small",
        "\\begin{tabular}{llccc}",
        "\\hline",
        " & & \\multicolumn{2}{c}{Modal age at death} & Acceptance \\\\",
        "Country & Population & Laplace & Exact & rate \\\\",
        "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(tab, file.path(out_dir, "app_laplace_mcmc.tex"))
writeLines(tab, file.path(ms_tab, "app_laplace_mcmc.tex"))

# ---- Figure -----------------------------------------------------------------

probs = probs[probs$prob >= 0.001, ]
probs$case = factor(probs$case, levels = unique(summ$case))
probs$method = factor(ifelse(probs$method == "MCMC", "Exact (MCMC)", "Laplace"),
                      levels = c("Laplace", "Exact (MCMC)"))
p = ggplot(probs, aes(age + 0.5, prob, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  facet_wrap(~ case, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("Laplace" = "#0072B2", "Exact (MCMC)" = "#E69F00"),
                    name = NULL) +
  labs(x = "Age", y = "P(mode = x | data)",
       caption = sprintf("RW2 posterior, starting age %d.", start_age)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", hjust = 0))
ggsave(file.path(fig_dir, "app_laplace_mcmc.pdf"), p, width = 11, height = 5.5)
ggsave(file.path(ms_fig, "app_laplace_mcmc.pdf"), p, width = 11, height = 5.5)
cat("Saved to", out_dir, ",", fig_dir, "and the Manuscript folders\n")
