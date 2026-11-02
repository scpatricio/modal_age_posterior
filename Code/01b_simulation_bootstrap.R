# =============================================================================
# 01b_simulation_bootstrap.R
# Bootstrap confidence intervals for the benchmark estimators in the Monte
# Carlo experiment (Section 3 of the manuscript). Run after 01_simulation.R.
#
# The script regenerates exactly the data sets of 01_simulation.R (same seed
# and task order). For each data set and each benchmark (Gompertz, P-spline,
# DDT) it draws B bootstrap samples D*_x ~ Poisson(D_x) with the exposures
# held fixed, re-estimates the mode on each, and takes the 2.5% and 97.5%
# quantiles as the 95% percentile interval. For the P-spline the smoothing
# parameter is fixed at the value chosen by BIC on the original data, which
# keeps the computation feasible.
#
# The intervals are scored in the same way as the posterior credible
# intervals: coverage of the true mode and mean width.
#
# Outputs (results/simulation/):
#   sim_bootstrap.rds        one row per replication x method
#   sim_bootstrap_summary.csv  coverage and mean width by scenario, N, method
#
# Uses forked parallelism (macOS / Linux). The number of bootstrap samples
# can be changed with the environment variable MODE_NBOOT (default 200).
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")
for (f in c("utils.R", "mode_gompertz.R", "mode_pspline.R", "mode_ddt.R",
            "simulation_truth.R")) source(f)

# ---- Settings (must match 01_simulation.R) ----------------------------------

n_rep = as.integer(Sys.getenv("MODE_NREP", 500))
n_boot = as.integer(Sys.getenv("MODE_NBOOT", 200))
sizes = c(1e3, 1e4, 1e5, 1e6)
start_age = 40
weights = c(a = 1, b = 0.80, c = 0.70, d = 0.60)
seed0 = 20261104
n_cores = get_cores()

out_dir = file.path("results", "simulation")
tr = readRDS(file.path(out_dir, "sim_truth.rds"))
true_mode = sapply(tr$truth, function(t) {
  p = t$peaks[t$peaks$x >= start_age, ]
  p$x[which.max(p$f)]
})

schedules = list()
for (sc in names(weights)) {
  for (N in sizes) schedules[[paste(sc, N)]] = expected_schedule(weights[[sc]], N)
}

methods = c("gompertz", "pspline", "ddt")

# ---- One replication --------------------------------------------------------

run_one = function(task) {
  set.seed(seed0 + task$id)
  sch = schedules[[paste(task$scenario, task$N)]]
  D = rpois(nrow(sch), sch$m * sch$exposure)
  E = sch$exposure
  age = sch$age
  safe = function(expr) tryCatch(expr, error = function(e) NA_real_)

  ps = tryCatch(mode_pspline(age, D, E, start_age, max_age), error = function(e) NULL)
  ll = if (is.null(ps)) NA_real_ else ps$log10_lambda
  fit = list(
    gompertz = function(d) safe(mode_gompertz(age, d, E, start_age, max_age)$mode),
    pspline = function(d) {
      if (!is.finite(ll)) return(NA_real_)
      safe(mode_pspline(age, d, E, start_age, max_age, log10_lambdas = ll)$mode)
    },
    ddt = function(d) safe(mode_ddt(age, d, E, start_age, max_age)$mode))

  point = sapply(methods, function(m) fit[[m]](D))

  # bootstrap samples, shared by the three methods
  set.seed(seed0 + 1e6 + task$id)
  boot = matrix(NA_real_, n_boot, length(methods), dimnames = list(NULL, methods))
  for (b in seq_len(n_boot)) {
    Db = rpois(length(D), D)
    for (m in methods) boot[b, m] = fit[[m]](Db)
  }

  tm = true_mode[[task$scenario]]
  do.call(rbind, lapply(methods, function(m) {
    v = boot[is.finite(boot[, m]), m]
    q = if (length(v) >= 0.9 * n_boot) quantile(v, c(0.025, 0.975), names = FALSE) else c(NA, NA)
    data.frame(scenario = task$scenario, N = task$N, rep = task$rep, method = m,
               estimate = point[[m]], lower = q[1], upper = q[2],
               n_boot_ok = length(v), true_mode = tm,
               cover = q[1] <= tm & tm <= q[2], width = q[2] - q[1])
  }))
}

# ---- Run --------------------------------------------------------------------

tasks = expand.grid(rep = seq_len(n_rep), N = sizes, scenario = names(weights),
                    stringsAsFactors = FALSE)
tasks$id = seq_len(nrow(tasks))
task_list = split(tasks, tasks$id)

cat("Bootstrap:", nrow(tasks), "replications x", n_boot, "samples on", n_cores, "cores\n")
t0 = Sys.time()
res = par_lapply(task_list, run_one, n_cores, "Bootstrap")
cat("Done in", format(Sys.time() - t0), "\n")

failed = sapply(res, inherits, "try-error")
if (any(failed)) warning(sum(failed), " replications failed; first error: ", res[[which(failed)[1]]])
boot = do.call(rbind, res[!failed])
saveRDS(boot, file.path(out_dir, "sim_bootstrap.rds"))

# check: point estimates should equal those of 01_simulation.R
est = readRDS(file.path(out_dir, "sim_estimates.rds"))
key = function(z) paste(z$scenario, z$N, z$rep, z$method)
m = match(key(boot), key(est[est$start_age == start_age, ]))
dif = abs(boot$estimate - est$estimate[m])
cat("Largest difference from the point estimates of 01_simulation.R:",
    signif(max(dif, na.rm = TRUE), 3), "\n")

g = split(boot, list(boot$scenario, boot$N, boot$method), drop = TRUE)
summ = do.call(rbind, lapply(g, function(z) {
  data.frame(scenario = z$scenario[1], N = z$N[1], start_age = start_age,
             method = z$method[1], n = nrow(z), n_no_interval = sum(!is.finite(z$width)),
             coverage = mean(z$cover, na.rm = TRUE),
             mean_width = mean(z$width, na.rm = TRUE))
}))
summ = summ[order(summ$scenario, summ$N, summ$method), ]
write.csv(summ, file.path(out_dir, "sim_bootstrap_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
cat("Saved results to", out_dir, "\n")
