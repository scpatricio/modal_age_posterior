# =============================================================================
# 06_laplace_vs_mcmc.R
# Checks the Laplace approximation of the RW2 posterior against posterior
# sampling of the same model (Appendix D): women in Denmark, France, Japan
# and Sweden in 1960 and 2019 (starting age 50), and four small simulated
# populations from the simulation study (scenarios a and d, N = 1000,
# replications 1 and 2, starting age 40).
#
# Sampling: independence Metropolis-Hastings (rw_imh_cpp() in rw_core.cpp)
# with the Laplace mixture over the grid of lambda as proposal. Every grid
# value is evaluated and has a positive proposal probability (Laplace weight
# mixed with 1% uniform), so the chain targets the posterior of the model on
# the full grid (Poisson likelihood, RW2 prior, flat prior on the grid of
# log10 lambda). The acceptance rate measures the quality of the
# approximation (1 if it were exact).
#
# Outputs:
#   results/mcmc_check/laplace_vs_mcmc.csv   posterior summaries, both methods
#   results/mcmc_check/app_laplace_mcmc.tex  the same as a LaTeX table
#   figures/mcmc_check/app_laplace_mcmc.pdf  P(mode = x | data), both methods
# (table and figure also go to the manuscript, see utils.R). The numbers
# quoted in Appendix D are printed at the end.
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")
for (f in c("utils.R", "mode_posterior.R", "mode_posterior_rw.R", "simulation_truth.R")) source(f)
library(ggplot2)
if (!rw_use_cpp) stop("This check needs the C++ core (rw_core.cpp); see mode_posterior_rw.R.")

data_dir = "data"
start_age = 50
max_age = 110
S = 20000                       # Laplace draws
n_iter = 100000                 # Metropolis-Hastings iterations
log10_lambdas = seq(-1, 8, by = 0.25)
drop = 12
cases = data.frame(code = rep(names(hmd_countries), each = 2),
                   country = rep(unname(hmd_countries), each = 2),
                   sex = "Female",
                   year = rep(c(1960, 2019), length(hmd_countries)),
                   start = start_age, scenario = NA, N = NA, rep = NA)
# small simulated populations, same data sets as in 01_simulation.R
sim_cases = expand.grid(rep = 1:2, scenario = c("a", "d"), stringsAsFactors = FALSE)
sim_cases = data.frame(code = NA, country = sprintf("Simulated (%s)", sim_cases$scenario),
                       sex = NA, year = NA, start = 40, scenario = sim_cases$scenario,
                       N = 1e3, rep = sim_cases$rep)
cases = rbind(cases, sim_cases)
sim_weights = c(a = 1, b = 0.80, c = 0.70, d = 0.60)
sim_data = function(cs) {
  # task id as in 01_simulation.R (500 replications, N = 10^3 ... 10^6)
  id = cs$rep + 500 * 0 + 500 * 4 * (match(cs$scenario, names(sim_weights)) - 1)
  set.seed(20261104 + id)
  sch = expected_schedule(sim_weights[[cs$scenario]], cs$N)
  data.frame(age = sch$age, deaths = rpois(nrow(sch), sch$m * sch$exposure),
             exposure = sch$exposure)
}
seed = 20261108

out_dir = file.path("results", "mcmc_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hmd = do.call(rbind, lapply(unique(na.omit(cases$code)), load_country, data_dir = data_dir))

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
  y = if (is.na(cs$scenario)) {
    hmd[hmd$country == cs$country & hmd$sex == cs$sex & hmd$year == cs$year, ]
  } else {
    sim_data(cs)
  }
  st = cs$start
  s = prep_schedule(y$age, y$deaths, y$exposure, st, max_age)
  nc = length(s$age) - 1
  D = s$deaths[1:nc]
  E = s$exposure[1:nc]
  set.seed(seed + k)

  # Laplace approximation (the method used in the paper)
  lap = mode_posterior_rw(y$age, y$deaths, y$exposure, st, max_age, S = S,
                          log10_lambdas = log10_lambdas)

  # posterior sampling: independence Metropolis-Hastings; all grid values
  # evaluated (drop = Inf) and proposed with positive probability
  g = rw_grid_cpp(D, E, 2L, log10_lambdas, Inf)
  w = exp(g$logml - max(g$logml))
  w = 0.99 * w / sum(w) + 0.01 / length(w)
  ch = rw_imh_cpp(D, E, g$eta, g$L, g$log10_lambda, w, n_iter, 2L)
  mc = posterior_summary(d_from_eta(ch$eta), s$age, st, n_iter, 0.95, 0.1, 0.5,
                         return_draws = TRUE)
  ll_chain = g$log10_lambda[ch$cell + 1]

  lab = if (is.na(cs$scenario)) sprintf("%s, %d", cs$country, cs$year) else
    sprintf("%s, N = %g, replication %d", cs$country, cs$N, cs$rep)
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

f2 = function(x) formatC(x, format = "f", digits = 4)
rows = sprintf("%s & %s & %s (%s, %s) & %s (%s, %s) & %s \\\\",
               L$country, ifelse(is.na(L$year), sprintf("rep. %d", L$rep), L$year),
               f2(L$mean), f2(L$q025), f2(L$q975), f2(M$mean), f2(M$q025), f2(M$q975),
               f2(M$acceptance))
tab = c("\\begin{table}[htb!]",
        "\\centering",
        "\\caption{RW2 posterior of the modal age at death computed with the Laplace approximation and by posterior sampling (independence Metropolis--Hastings): posterior mean of $\\tilde{M}$ with 95\\% credible interval, and acceptance rate of the sampler. Women, starting age 50; simulated populations (scenarios a and d of Section~\\ref{sec:simulation}, $N = 10^3$), starting age 40.}",
        "\\label{tab:app_mcmc}",
        "\\small",
        "\\begin{tabular}{llccc}",
        "\\hline",
        " & & \\multicolumn{2}{c}{Modal age at death} & Acceptance \\\\",
        "Population & Year & Laplace & Sampling & rate \\\\",
        "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{table}")
save_table(tab, "app_laplace_mcmc.tex", "mcmc_check")

# ---- Figure -----------------------------------------------------------------

probs = probs[probs$prob >= 0.001 & !grepl("^Simulated", probs$case), ]
probs$case = factor(probs$case, levels = unique(summ$case))
probs$method = factor(ifelse(probs$method == "MCMC", "Sampling (MCMC)", "Laplace"),
                      levels = c("Laplace", "Sampling (MCMC)"))
p = ggplot(probs, aes(age + 0.5, prob, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  facet_wrap(~ case, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("Laplace" = "#0072B2", "Sampling (MCMC)" = "#E69F00"),
                    name = NULL) +
  labs(x = "Age", y = "P(mode = x | data)") +
  theme_mode()
save_figure(p, "app_laplace_mcmc.pdf", "mcmc_check", 6)
cat("Saved to results/mcmc_check and figures/mcmc_check\n")
