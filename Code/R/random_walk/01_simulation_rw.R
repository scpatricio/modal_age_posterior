# =============================================================================
# 01_simulation_rw.R
# Monte Carlo experiment with the random-walk posteriors.
#
# Same design as ../01_simulation.R: age at death is a mixture of two
# Gompertz distributions (modes 90 and 65, b = 0.14) with
#   (a) w = 1 (unimodal), (b) w = 0.80, (c) w = 0.70, (d) w = 0.60
# of deaths in the main component; D_x ~ Poisson(m_x E_x) with expected
# exposures for radix N.
#
# Methods (starting age 40):
#   raw_*   posterior with independent gamma priors (../mode_posterior.R)
#   rw2_*   posterior with a second-order random-walk prior
#           (* = mean, median, mode of the posterior)
#   gompertz, pspline, ddt
#
# Outputs (random_walk/Results/simulation/):
#   sim_truth.rds          scenario parameters, true density, peaks
#   sim_estimates.rds      one row per replication x start age x method
#   sim_posteriors.rds     RW2 posterior of the mode (mean over replications,
#                          and a few individual replications); used in panel B
#   sim_posteriors_all.rds the same for raw and RW2
#   sim_summary.csv        bias, MAE, RMSE, SD, variance, coverage of the
#                          95% credible interval, credible-interval width
#
# Run from Code/R/random_walk. Uses forked parallelism (macOS / Linux).
# =============================================================================

code_dir = Sys.getenv("MODE_RW_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk")
setwd(code_dir)
for (f in c("../utils.R", "../mode_posterior.R", "../mode_gompertz.R",
            "../mode_pspline.R", "../mode_ddt.R", "../sim_truth_utils.R",
            "mode_posterior_rw.R")) source(f)

# ---- Settings ---------------------------------------------------------------

n_rep = as.integer(Sys.getenv("MODE_NREP", 500))
S = 5000
sizes = c(1e3, 1e4, 1e5, 1e6)
start_ages = 40
weights = c(a = 1, b = 0.80, c = 0.70, d = 0.60)
n_keep_post = 5
seed0 = 20261104
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
true_mode = matrix(sapply(truth, function(t) sapply(start_ages, function(st) {
  p = t$peaks[t$peaks$x >= st, ]
  p$x[which.max(p$f)]
})), nrow = length(start_ages), dimnames = list(start_ages, names(truth)))

# ---- One replication --------------------------------------------------------

families = c("raw", "rw2")
methods = c(paste0(rep(families, each = 3), c("_mean", "_median", "_mode")),
            "gompertz", "pspline", "ddt")

run_one = function(task) {
  set.seed(seed0 + task$id)
  sch = schedules[[paste(task$scenario, task$N)]]
  D = rpois(nrow(sch), sch$m * sch$exposure)
  safe = function(expr) tryCatch(expr, error = function(e) NA_real_)

  # one list element per starting age (the functions return the result
  # itself when there is only one)
  by_start = function(p) if (length(start_ages) == 1) list(p) else p
  post = list(raw = by_start(mode_posterior(sch$age, D, sch$exposure, start_ages,
                                            max_age, S = S)),
              rw2 = by_start(mode_posterior_rw(sch$age, D, sch$exposure, start_ages,
                                               max_age, S = S)))
  n_st = length(start_ages)
  est = matrix(NA_real_, n_st, length(methods), dimnames = list(NULL, methods))
  fam_stats = array(NA_real_, c(n_st, length(families), 4),
                    dimnames = list(NULL, families,
                                    c("cover", "width", "cred_size", "log10_lambda")))
  post_probs = list()
  for (i in seq_len(n_st)) {
    st = start_ages[i]
    tm = true_mode[i, task$scenario]
    pe = c()
    for (fm in families) {
      p = post[[fm]][[i]]
      pe = c(pe, p$mean, p$median, p$mode)
      fam_stats[i, fm, ] = c(p$q025 <= tm && tm <= p$q975, p$q975 - p$q025, p$cred_size,
                             if (is.null(p$post_log10_lambda)) NA else p$post_log10_lambda)
      post_probs[[paste(fm, st)]] = data.frame(family = fm, start_age = st,
                                               age = p$post$age, prob = p$post$prob)
    }
    est[i, ] = c(pe,
                 safe(mode_gompertz(sch$age, D, sch$exposure, st, max_age)$mode),
                 safe(mode_pspline(sch$age, D, sch$exposure, st, max_age)$mode),
                 safe(mode_ddt(sch$age, D, sch$exposure, st, max_age)$mode))
  }
  list(est = est, fam_stats = fam_stats, post = do.call(rbind, post_probs))
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
  fam = sub("_.*", "", methods)
  fam_idx = match(fam, families)
  stat = function(name) {
    v = matrix(NA_real_, n_st, length(methods))
    for (m in which(!is.na(fam_idx))) v[, m] = r$fam_stats[, fam_idx[m], name]
    as.vector(v)
  }
  data.frame(scenario = tk$scenario, N = tk$N, rep = tk$rep,
             start_age = rep(start_ages, times = length(methods)),
             method = rep(methods, each = n_st),
             estimate = as.vector(r$est),
             cover = stat("cover"), width = stat("width"), cred_size = stat("cred_size"),
             log10_lambda = stat("log10_lambda"))
}))
est$true_mode = true_mode[cbind(match(est$start_age, start_ages),
                                match(est$scenario, scenarios$scenario))]
est$error = est$estimate - est$true_mode
saveRDS(est, file.path(out_dir, "sim_estimates.rds"))

post_long = do.call(rbind, lapply(ok, function(j) {
  tk = tasks[j, ]
  data.frame(scenario = tk$scenario, N = tk$N, rep = tk$rep, res[[j]]$post)
}))
key = paste(post_long$family, post_long$scenario, post_long$N, post_long$start_age,
            post_long$age)
n_ok = table(paste(tasks$scenario[ok], tasks$N[ok]))
sums = rowsum(post_long$prob, key)
parts = do.call(rbind, strsplit(rownames(sums), " "))
post_mean = data.frame(family = parts[, 1], scenario = parts[, 2],
                       N = as.numeric(parts[, 3]), start_age = as.numeric(parts[, 4]),
                       age = as.numeric(parts[, 5]), prob = sums[, 1])
post_mean$prob = post_mean$prob / as.numeric(n_ok[paste(post_mean$scenario, post_mean$N)])
post_mean = post_mean[order(post_mean$family, post_mean$scenario, post_mean$N,
                            post_mean$start_age, post_mean$age), ]
post_ind = post_long[post_long$rep <= n_keep_post, ]
saveRDS(list(mean = post_mean, individual = post_ind),
        file.path(out_dir, "sim_posteriors_all.rds"))
# RW2 only, in the layout read by ../sim_plot_utils.R (panel B)
saveRDS(list(mean = post_mean[post_mean$family == "rw2", -1],
             individual = post_ind[post_ind$family == "rw2", names(post_ind) != "family"]),
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
             mean_width = mean(z$width),
             mean_cred_size = mean(z$cred_size),
             mean_log10_lambda = mean(z$log10_lambda))
}))
summ = summ[order(summ$scenario, summ$N, summ$start_age, summ$method), ]
write.csv(summ, file.path(out_dir, "sim_summary.csv"), row.names = FALSE)
cat("Saved results to", out_dir, "\n")
