# =============================================================================
# 01_simulation.R
# Monte Carlo experiment: can each method recover the modal age at death
# when the age-at-death distribution has a secondary peak?
#
# Age at death follows a mixture of two Gompertz distributions with the same
# slope b,
#   f(x) = w g1(x) + (1 - w) g2(x),
# g1: main component (mode 90), g2: earlier component (mode 65). The mixture
# describes a population with two subgroups of different mortality levels.
# Scenarios differ in the share w of deaths in the main component:
#   (a) w = 1      unimodal
#   (b) w = 0.80   20% of deaths in the earlier component
#   (c) w = 0.70   30%
#   (d) w = 0.60   40%
# The global mode stays at 90 in every scenario (the height of the lower
# peak relative to the main one is saved in sim_truth.rds, peak_ratio).
#
# Data: for population size N (radix at age 0), exposures are the expected
# person-years E_x = N L_x, and deaths are D_x ~ Poisson(m_x E_x) with m_x
# the true central rates.
#
# Outputs (Results/simulation/):
#   sim_truth.rds       scenario parameters, true density, peaks, valleys
#   sim_estimates.rds   one row per replication x start age x method
#   sim_posteriors.rds  posterior of the mode: mean over replications and
#                       the first few individual replications
#   sim_summary.csv     bias, mean absolute error, RMSE, SD, variance,
#                       coverage of the 95% credible interval of the
#                       continuous mode (2.5% and 97.5% quantiles)
#
# Run from Code/R. Uses forked parallelism (macOS / Linux).
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
for (f in c("utils.R", "mode_posterior.R", "mode_gompertz.R",
            "mode_pspline.R", "mode_ddt.R", "sim_truth_utils.R")) source(f)

# ---- Settings ---------------------------------------------------------------

n_rep = as.integer(Sys.getenv("MODE_NREP", 500))   # replications per cell
S = 5000                                            # posterior draws
sizes = c(1e3, 1e4, 1e5, 1e6)                       # radix N
start_ages = c(30, 40, 50, 60)
weights = c(a = 1, b = 0.80, c = 0.70, d = 0.60)
n_keep_post = 5                                     # individual posteriors kept
seed0 = 20261102
n_cores = get_cores()

out_dir = file.path("Results", "simulation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- True distributions -----------------------------------------------------

scenarios = data.frame(scenario = names(weights), w = unname(weights))
scenarios$label = sprintf("(%s) %s", scenarios$scenario,
                          ifelse(scenarios$w == 1, "Unimodal",
                                 sprintf("%.0f%% of deaths in the earlier component",
                                         100 * (1 - scenarios$w))))

truth = lapply(seq_len(nrow(scenarios)), function(k) {
  w = scenarios$w[k]
  pv = find_peaks(w)
  grid = data.frame(scenario = scenarios$scenario[k], age = seq(0, max_age, by = 0.05))
  grid$f = mix_f(grid$age, w)
  list(scenario = scenarios$scenario[k], w = w, peaks = pv$peaks,
       valleys = pv$valleys, density = grid)
})
names(truth) = scenarios$scenario
scenarios$peak_ratio = sapply(truth, function(t) {
  if (nrow(t$peaks) < 2) NA else min(t$peaks$f) / max(t$peaks$f)
})

saveRDS(list(scenarios = scenarios, truth = truth,
             settings = list(M1 = M1, b1 = b1, M2 = M2, b2 = b2, sizes = sizes,
                             start_ages = start_ages, n_rep = n_rep, S = S)),
        file.path(out_dir, "sim_truth.rds"))

schedules = list()
for (k in seq_len(nrow(scenarios))) {
  for (N in sizes) {
    schedules[[paste(scenarios$scenario[k], N)]] = expected_schedule(scenarios$w[k], N)
  }
}

# true mode on [start_age, max_age] for each scenario and start age
true_mode = sapply(truth, function(t) sapply(start_ages, function(st) {
  p = t$peaks[t$peaks$x >= st, ]
  p$x[which.max(p$f)]
}))
rownames(true_mode) = start_ages

# ---- One replication --------------------------------------------------------

methods = c("posterior_mean", "posterior_median", "posterior_mode",
            "gompertz", "pspline", "ddt")

run_one = function(task) {
  set.seed(seed0 + task$id)
  sch = schedules[[paste(task$scenario, task$N)]]
  D = rpois(nrow(sch), sch$m * sch$exposure)
  tr = truth[[task$scenario]]
  safe = function(expr) tryCatch(expr, error = function(e) NA_real_)

  pp = mode_posterior(sch$age, D, sch$exposure, start_ages, max_age, S = S)
  n_st = length(start_ages)
  est = matrix(NA_real_, n_st, length(methods), dimnames = list(NULL, methods))
  extra = matrix(NA_real_, n_st, 4,
                 dimnames = list(NULL, c("cover", "cred_size", "p_low", "exp_major_modes")))
  post_probs = vector("list", n_st)
  for (i in seq_len(n_st)) {
    st = start_ages[i]
    p = pp[[i]]
    est[i, ] = c(p$mean, p$median, p$mode,
                 safe(mode_gompertz(sch$age, D, sch$exposure, st, max_age)$mode),
                 safe(mode_pspline(sch$age, D, sch$exposure, st, max_age)$mode),
                 safe(mode_ddt(sch$age, D, sch$exposure, st, max_age)$mode))
    valley = tr$valleys[tr$valleys >= st]
    tmi = true_mode[i, task$scenario]
    extra[i, ] = c(p$q025 <= tmi && tmi <= p$q975,
                   p$cred_size,
                   if (length(valley) > 0) sum(p$post$prob[p$post$age < valley[1]]) else NA,
                   p$exp_major_modes)
    post_probs[[i]] = p$post
  }
  list(est = est, extra = extra, post = post_probs)
}

# ---- Run --------------------------------------------------------------------

tasks = expand.grid(rep = seq_len(n_rep), N = sizes, scenario = scenarios$scenario,
                    stringsAsFactors = FALSE)
tasks$id = seq_len(nrow(tasks))
task_list = split(tasks, tasks$id)

cat("Running", nrow(tasks), "replications on", n_cores, "cores\n")
t0 = Sys.time()
res = par_lapply(task_list, run_one, n_cores, "Simulation")
cat("Done in", format(Sys.time() - t0), "\n")

failed = sapply(res, inherits, "try-error")
if (any(failed)) {
  warning(sum(failed), " replications failed; first error: ", res[[which(failed)[1]]])
}
ok = which(!failed)

# ---- Assemble ---------------------------------------------------------------

n_st = length(start_ages)
est = do.call(rbind, lapply(ok, function(j) {
  r = res[[j]]
  tk = tasks[j, ]
  data.frame(scenario = tk$scenario, N = tk$N, rep = tk$rep,
             start_age = rep(start_ages, times = length(methods)),
             method = rep(methods, each = n_st),
             estimate = as.vector(r$est),
             cover = rep(r$extra[, "cover"], times = length(methods)),
             cred_size = rep(r$extra[, "cred_size"], times = length(methods)),
             p_low = rep(r$extra[, "p_low"], times = length(methods)),
             exp_major_modes = rep(r$extra[, "exp_major_modes"], times = length(methods)))
}))
is_post = startsWith(est$method, "posterior")
est[!is_post, c("cover", "cred_size", "p_low", "exp_major_modes")] = NA
est$true_mode = true_mode[cbind(match(est$start_age, start_ages),
                                match(est$scenario, scenarios$scenario))]
est$error = est$estimate - est$true_mode
saveRDS(est, file.path(out_dir, "sim_estimates.rds"))

# posteriors: mean over replications (fast with rowsum) and a few individual ones
post_long = do.call(rbind, lapply(ok, function(j) {
  tk = tasks[j, ]
  do.call(rbind, lapply(seq_len(n_st), function(i) {
    p = res[[j]]$post[[i]]
    data.frame(scenario = tk$scenario, N = tk$N, rep = tk$rep,
               start_age = start_ages[i], age = p$age, prob = p$prob)
  }))
}))
key = paste(post_long$scenario, post_long$N, post_long$start_age, post_long$age)
n_ok = table(paste(tasks$scenario[ok], tasks$N[ok]))
sums = rowsum(post_long$prob, key)
parts = do.call(rbind, strsplit(rownames(sums), " "))
post_mean = data.frame(scenario = parts[, 1], N = as.numeric(parts[, 2]),
                       start_age = as.numeric(parts[, 3]), age = as.numeric(parts[, 4]),
                       prob = sums[, 1])
post_mean$prob = post_mean$prob / as.numeric(n_ok[paste(post_mean$scenario, post_mean$N)])
post_mean = post_mean[order(post_mean$scenario, post_mean$N, post_mean$start_age,
                            post_mean$age), ]
saveRDS(list(mean = post_mean, individual = post_long[post_long$rep <= n_keep_post, ]),
        file.path(out_dir, "sim_posteriors.rds"))

# ---- Summary table ----------------------------------------------------------

g = split(est, list(est$scenario, est$N, est$start_age, est$method), drop = TRUE)
summ = do.call(rbind, lapply(g, function(z) {
  data.frame(scenario = z$scenario[1], N = z$N[1], start_age = z$start_age[1],
             method = z$method[1], true_mode = z$true_mode[1], n = nrow(z),
             n_failed = sum(!is.finite(z$estimate)),
             bias = mean(z$error, na.rm = TRUE),
             mae = mean(abs(z$error), na.rm = TRUE),
             rmse = sqrt(mean(z$error^2, na.rm = TRUE)),
             sd = sd(z$estimate, na.rm = TRUE),
             variance = var(z$estimate, na.rm = TRUE),
             coverage = mean(z$cover),
             mean_cred_size = mean(z$cred_size),
             mean_p_low = mean(z$p_low))
}))
summ = summ[order(summ$scenario, summ$N, summ$start_age, summ$method), ]
write.csv(summ, file.path(out_dir, "sim_summary.csv"), row.names = FALSE)
cat("Saved results to", out_dir, "\n")
