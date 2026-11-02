# =============================================================================
# 05_simulation_change.R
# Two experiments on the unimodal Gompertz distribution (mode 90, b = 0.14):
#
# 1. Detecting change. Year 1 has mode 90, year 2 has mode 90 + delta.
#    From the two posteriors (independent populations), P(M2 > M1 | data)
#    is the share of paired draws with M2 > M1. We record how often it
#    exceeds 0.95 (delta = 0 gives the false-positive rate) and, for every
#    method, how often the estimated change has the right sign.
#
# 2. Overdispersion. True rates are multiplied by independent noise
#    eps_x ~ Gamma(1/cv^2, 1/cv^2) (mean 1, coefficient of variation cv)
#    before drawing D_x ~ Poisson(m_x eps_x E_x). The posterior assumes no
#    such noise, so this checks how coverage and accuracy degrade when the
#    counts vary more than Poisson. cv = 0 is the pure Poisson case.
#
# All methods use starting age 40. The year-1 data set is shared by the
# four values of delta within a replication.
#
# Outputs (Results/change/):
#   change_reps.rds          one row per replication x delta
#   year1_reps.rds           year-1 estimates (overdispersion experiment)
#   change_summary.csv       detection rate, mean P(M2 > M1), sign accuracy
#   overdisp_summary.csv     bias, MAE, SD, coverage by cv, N and method
#
# Run from Code/R. Uses forked parallelism (macOS / Linux).
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
for (f in c("utils.R", "mode_posterior.R", "mode_gompertz.R",
            "mode_pspline.R", "mode_ddt.R", "sim_truth_utils.R")) source(f)

# ---- Settings ---------------------------------------------------------------

n_rep = as.integer(Sys.getenv("MODE_NREP", 500))
S = 5000
sizes = c(1e3, 1e4, 1e5, 1e6)
deltas = c(0, 0.25, 0.5, 1)
cvs = c(0, 0.02, 0.05, 0.10)
start_age = 40
M0 = M1                 # baseline mode (90)
b0 = b1                 # slope (0.14)
threshold = 0.95
seed0 = 20261103
n_cores = get_cores()

out_dir = file.path("Results", "change")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# expected schedules: baseline and shifted modes, for every N
schedules = list()
for (N in sizes) {
  for (M in M0 + c(0, deltas[deltas > 0])) {
    schedules[[paste(M, N)]] = expected_schedule_S(function(t) gomp_S(t, M, b0), N)
  }
}

methods = c("posterior_mean", "posterior_median", "posterior_mode",
            "gompertz", "pspline", "ddt")

# ---- One replication --------------------------------------------------------

simulate_D = function(sch, cv) {
  eps = if (cv > 0) rgamma(nrow(sch), shape = 1 / cv^2, rate = 1 / cv^2) else 1
  rpois(nrow(sch), sch$m * eps * sch$exposure)
}

fit_all = function(sch, D) {
  safe = function(expr) tryCatch(expr, error = function(e) NA_real_)
  p = mode_posterior(sch$age, D, sch$exposure, start_age, max_age, S = S,
                     return_draws = TRUE)
  list(post = p,
       est = c(p$mean, p$median, p$mode,
               safe(mode_gompertz(sch$age, D, sch$exposure, start_age, max_age)$mode),
               safe(mode_pspline(sch$age, D, sch$exposure, start_age, max_age)$mode),
               safe(mode_ddt(sch$age, D, sch$exposure, start_age, max_age)$mode)))
}

run_one = function(task) {
  set.seed(seed0 + task$id)
  sch1 = schedules[[paste(M0, task$N)]]
  f1 = fit_all(sch1, simulate_D(sch1, task$cv))
  p1 = f1$post

  year1 = data.frame(N = task$N, cv = task$cv, rep = task$rep, method = methods,
                     estimate = f1$est,
                     cover = p1$q025 <= M0 && M0 <= p1$q975,
                     width = p1$q975 - p1$q025)

  change = lapply(deltas, function(dl) {
    sch2 = schedules[[paste(M0 + dl, task$N)]]
    f2 = fit_all(sch2, simulate_D(sch2, task$cv))
    data.frame(N = task$N, cv = task$cv, rep = task$rep, delta = dl,
               p_increase = mean(f2$post$draws_cont > p1$draws_cont),
               method = methods, change = f2$est - f1$est)
  })
  list(year1 = year1, change = do.call(rbind, change))
}

# ---- Run --------------------------------------------------------------------

tasks = expand.grid(rep = seq_len(n_rep), cv = cvs, N = sizes)
tasks$id = seq_len(nrow(tasks))
task_list = split(tasks, tasks$id)

cat("Running", nrow(tasks), "replications (", length(deltas) + 1,
    "data sets each) on", n_cores, "cores\n")
t0 = Sys.time()
res = par_lapply(task_list, run_one, n_cores, "Change experiment")
cat("Done in", format(Sys.time() - t0), "\n")

failed = sapply(res, inherits, "try-error")
if (any(failed)) {
  warning(sum(failed), " replications failed; first error: ", res[[which(failed)[1]]])
  res = res[!failed]
}

year1 = do.call(rbind, lapply(res, `[[`, "year1"))
year1$error = year1$estimate - M0
change = do.call(rbind, lapply(res, `[[`, "change"))
saveRDS(year1, file.path(out_dir, "year1_reps.rds"))
saveRDS(change, file.path(out_dir, "change_reps.rds"))

# ---- Summaries --------------------------------------------------------------

g = split(change, list(change$N, change$cv, change$delta, change$method), drop = TRUE)
change_summ = do.call(rbind, lapply(g, function(z) {
  data.frame(N = z$N[1], cv = z$cv[1], delta = z$delta[1], method = z$method[1],
             n = nrow(z),
             detect = mean(z$p_increase > threshold),      # posterior decision
             mean_p_increase = mean(z$p_increase),
             sign_correct = if (z$delta[1] > 0) mean(z$change > 0, na.rm = TRUE) else NA,
             mean_change = mean(z$change, na.rm = TRUE),
             sd_change = sd(z$change, na.rm = TRUE))
}))
change_summ = change_summ[order(change_summ$N, change_summ$cv, change_summ$delta,
                                change_summ$method), ]
write.csv(change_summ, file.path(out_dir, "change_summary.csv"), row.names = FALSE)

g = split(year1, list(year1$N, year1$cv, year1$method), drop = TRUE)
od_summ = do.call(rbind, lapply(g, function(z) {
  data.frame(N = z$N[1], cv = z$cv[1], method = z$method[1], n = nrow(z),
             bias = mean(z$error, na.rm = TRUE),
             mae = mean(abs(z$error), na.rm = TRUE),
             sd = sd(z$estimate, na.rm = TRUE),
             coverage = mean(z$cover),
             mean_width = mean(z$width))
}))
od_summ = od_summ[order(od_summ$N, od_summ$cv, od_summ$method), ]
write.csv(od_summ, file.path(out_dir, "overdisp_summary.csv"), row.names = FALSE)
cat("Saved results to", out_dir, "\n")
