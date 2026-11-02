# =============================================================================
# 08_review_checks.R
# Additional checks requested in review. Run after 01_simulation.R,
# 01b_simulation_bootstrap.R and 04_empirical.R. Each part can be switched
# off with the environment variable named in its header.
#
#   Part 1  Change between 2019 and 2020 (Section 4): posterior of the
#           difference in the modal age, and the posterior distribution of
#           the mode in each year. Writes Table 2 and Figure 6 of the paper.
#   Part 2  Interpolation without noise (Section 2.3): the interpolated mode
#           of the true schedule against the true mode of the density, for
#           true modes off the integer grid.
#   Part 3  Grid of the smoothing parameter: posterior weight at the ends of
#           the grid of log10 lambda, in the application.
#   Part 4  P-spline bootstrap with the smoothing parameter re-selected in
#           every bootstrap sample, against the bootstrap with the smoothing
#           parameter fixed (subset of the simulation). MODE_SKIP_PSPLINE=1
#           skips it (about 30 minutes on 8 cores).
#
# Outputs: results/review_checks/ (and the table and figure in the
# manuscript folder). The numbers are also written to
# results/review_checks/review_numbers.txt.
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")
for (f in c("utils.R", "mode_posterior.R", "mode_posterior_rw.R", "mode_pspline.R",
            "simulation_truth.R")) source(f)
library(ggplot2)
theme_set(theme_mode())

out_dir = file.path("results", "review_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
n_cores = get_cores()
max_age = 110
log_env = new.env()
log_env$lines = c()
say = function(...) {
  txt = paste0(...)
  log_env$lines = c(log_env$lines, txt)
  cat(txt, "\n")
}
r2 = function(x) round(x, 2)
r3 = function(x) round(x, 3)

# ---- Part 1: change between 2019 and 2020 ----------------------------------

data_dir = "data"
start_age = 40
S = as.numeric(Sys.getenv("MODE_CHECK_S", 1e5))
hmd = do.call(rbind, lapply(names(hmd_countries), load_country, data_dir = data_dir))

say("== Part 1: change between 2019 and 2020 (RW2, starting age ", start_age, ", S = ", S, ") ==")
change = list()
post = list()
for (k in seq_along(hmd_countries)) {
  cc = unname(hmd_countries[k])
  fits = lapply(c(2019, 2020), function(yr) {
    set.seed(20261110 + 10 * k + yr - 2019)
    y = hmd[hmd$country == cc & hmd$year == yr, ]
    mode_posterior_rw(y$age, y$deaths, y$exposure, start_age, max_age, S = S,
                      return_draws = TRUE)
  })
  # the two years are fitted separately, so independent draws are paired
  delta = fits[[2]]$draws_cont - fits[[1]]$draws_cont
  change[[cc]] = data.frame(country = cc,
    map_2019 = fits[[1]]$map, p_map_2019 = fits[[1]]$p_map,
    mean_2019 = fits[[1]]$mean, lo_2019 = fits[[1]]$q025, hi_2019 = fits[[1]]$q975,
    map_2020 = fits[[2]]$map, p_map_2020 = fits[[2]]$p_map,
    mean_2020 = fits[[2]]$mean, lo_2020 = fits[[2]]$q025, hi_2020 = fits[[2]]$q975,
    delta_mean = mean(delta), delta_lo = quantile(delta, 0.025, names = FALSE),
    delta_hi = quantile(delta, 0.975, names = FALSE), p_decline = mean(delta < 0))
  for (i in 1:2) {
    post[[paste(cc, i)]] = data.frame(country = cc, year = c(2019, 2020)[i],
                                      age = fits[[i]]$post$age, prob = fits[[i]]$post$prob)
  }
  z = change[[cc]]
  say(sprintf("%s: 2019 most probable age %d (prob %s), mean %s (%s-%s); 2020 most probable age %d (prob %s), mean %s (%s-%s); change %s (%s to %s); P(decline) %s",
              cc, z$map_2019, r2(z$p_map_2019), r2(z$mean_2019), r2(z$lo_2019), r2(z$hi_2019),
              z$map_2020, r2(z$p_map_2020), r2(z$mean_2020), r2(z$lo_2020), r2(z$hi_2020),
              r2(z$delta_mean), r2(z$delta_lo), r2(z$delta_hi), r3(z$p_decline)))
}
change = do.call(rbind, change)
post = do.call(rbind, post)
write.csv(change, file.path(out_dir, "change_2019_2020.csv"), row.names = FALSE)

# Table 2
f2 = function(x) {
  out = formatC(x, format = "f", digits = 4)
  ifelse(substr(out, 1, 1) == "-", paste0("$-$", substring(out, 2)), out)
}
rows = sprintf("%s & %d (%s) & %d (%s) & %s (%s, %s) & %s \\\\", change$country,
               change$map_2019, f2(change$p_map_2019), change$map_2020, f2(change$p_map_2020),
               f2(change$delta_mean), f2(change$delta_lo), f2(change$delta_hi),
               f2(change$p_decline))
tab = c("\\begin{table}[htb!]",
        "\\centering",
        sprintf(paste0("\\caption{Modal age at death of women in 2019 and 2020 (RW2 posterior, starting age %d): ",
                       "most probable single-year age with its posterior probability, posterior mean of the ",
                       "change in $\\tilde{M}$ from 2019 to 2020 with 95\\%% credible interval, and posterior ",
                       "probability that the modal age declined.}"), start_age),
        "\\label{tab:change}",
        "\\small",
        "\\begin{tabular}{lcccc}",
        "\\hline",
        " & \\multicolumn{2}{c}{Most probable age (probability)} & Change, 2019 to 2020 & \\\\",
        "Country & 2019 & 2020 & mean (95\\% interval) & P(decline) \\\\",
        "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{table}")
save_table(tab, "emp_change.tex", "review_checks")

# Figure: posterior distribution of the mode, 2019 and 2020
write.csv(post, file.path(out_dir, "posterior_2019_2020.csv"), row.names = FALSE)
post = post[post$prob >= 0.001, ]
post$year = factor(post$year)
p = ggplot(post, aes(age + 0.5, prob, fill = year)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  facet_wrap(~ country, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c("2019" = "#0072B2", "2020" = "#E69F00"), name = NULL) +
  scale_x_continuous(breaks = function(l) seq(ceiling(l[1]), floor(l[2]), by = 1)) +
  labs(x = "Age (bars: one-year age intervals)", y = "P(mode = x | data)")
save_figure(p, "emp_posterior_2019_2020.pdf", "review_checks", 4)

# ---- Part 2: interpolation without noise -----------------------------------

say("\n== Part 2: interpolated mode of the true schedule (Gompertz, b = ", b1, ", starting age 40) ==")
noise_free = do.call(rbind, lapply(c(90, 90.25, 90.5, 90.75), function(M) {
  Sfun = function(t) gomp_S(t, M, b1)
  sch = expected_schedule_S(Sfun, 1)
  keep = sch$age >= 40 & sch$age < max_age
  ag = sch$age[keep]
  # constant-force life table from the true central rates (as in the method)
  m = sch$m[keep]
  q = -expm1(-m)
  l = exp(-c(0, cumsum(m)[-length(m)]))
  d_cf = l * q
  # exact interval probabilities of the density
  d_ex = Sfun(ag) - Sfun(ag + 1)
  interp = function(d) {
    i = which.max(d)
    parabola_mode(ag[i], d[i - 1], d[i], d[i + 1])
  }
  data.frame(true_mode = M, modal_interval = ag[which.max(d_cf)],
             interp_constant_force = interp(d_cf), interp_exact = interp(d_ex))
}))
noise_free$error = noise_free$interp_constant_force - noise_free$true_mode
write.csv(noise_free, file.path(out_dir, "interpolation_noise_free.csv"), row.names = FALSE)
for (i in seq_len(nrow(noise_free))) {
  z = noise_free[i, ]
  say(sprintf("true mode %s: modal interval [%d, %d); interpolated %s (exact interval probabilities %s); error %s",
              z$true_mode, z$modal_interval, z$modal_interval + 1, round(z$interp_constant_force, 4),
              round(z$interp_exact, 4), round(z$error, 4)))
}
say("largest absolute error: ", round(max(abs(noise_free$error)), 4))

# ---- Part 3: ends of the grid of the smoothing parameter -------------------

say("\n== Part 3: posterior weight at the ends of the grid of log10 lambda (application, starting age 40) ==")
grid = seq(-1, 8, by = 0.25)
tasks = unique(hmd[hmd$year >= 1960, c("country", "year")])
ends = par_lapply(split(tasks, seq_len(nrow(tasks))), function(tk) {
  y = hmd[hmd$country == tk$country & hmd$year == tk$year, ]
  set.seed(1)
  h = mode_posterior_rw(y$age, y$deaths, y$exposure, start_age, max_age, S = 100,
                        log10_lambdas = grid)$hyper
  data.frame(tk, w_lower = h$weight[h$log10_lambda == min(grid)],
             w_upper = if (max(h$log10_lambda) == max(grid)) h$weight[h$log10_lambda == max(grid)] else 0,
             post_mean = sum(h$weight * h$log10_lambda),
             range_lo = min(h$log10_lambda[h$weight > 1e-3]),
             range_hi = max(h$log10_lambda[h$weight > 1e-3]))
}, n_cores, "Grid check")
ends = do.call(rbind, ends)
write.csv(ends, file.path(out_dir, "lambda_grid_ends.csv"), row.names = FALSE)
say("population-years: ", nrow(ends))
say("largest weight at log10 lambda = ", min(grid), ": ", signif(max(ends$w_lower), 3),
    "; at log10 lambda = ", max(grid), ": ", signif(max(ends$w_upper), 3))
say("posterior mean of log10 lambda: range ", paste(r2(range(ends$post_mean)), collapse = " to "),
    "; grid values with weight > 0.001 lie within ", min(ends$range_lo), " to ", max(ends$range_hi))
sim = readRDS(file.path("results", "simulation", "sim_estimates.rds"))
sl = sim$log10_lambda[sim$method == "rw2_mean"]
say("simulation, posterior mean of log10 lambda per replication: range ",
    paste(r2(range(sl, na.rm = TRUE)), collapse = " to "))

# ---- Part 4: P-spline bootstrap with the smoothing parameter re-selected ---

if (Sys.getenv("MODE_SKIP_PSPLINE") != "1") {
  n_rep = as.integer(Sys.getenv("MODE_CHECK_NREP", 200))
  n_boot = as.integer(Sys.getenv("MODE_CHECK_NBOOT", 500))
  seed0 = 20261104                       # as in 01_simulation.R
  sizes = c(1e3, 1e4, 1e5, 1e6)          # as in 01_simulation.R (task ids)
  weights = c(a = 1, b = 0.80, c = 0.70, d = 0.60)
  check = expand.grid(rep = seq_len(n_rep), N = c(1e3, 1e6), scenario = c("a", "d"),
                      stringsAsFactors = FALSE)
  n_rep_sim = 500
  check$id = check$rep + n_rep_sim * (match(check$N, sizes) - 1) +
    n_rep_sim * length(sizes) * (match(check$scenario, names(weights)) - 1)
  tr = readRDS(file.path("results", "simulation", "sim_truth.rds"))
  true_mode = sapply(tr$truth, function(t) {
    p = t$peaks[t$peaks$x >= start_age, ]
    p$x[which.max(p$f)]
  })
  schedules = list()
  for (sc in c("a", "d")) for (N in c(1e3, 1e6)) {
    schedules[[paste(sc, N)]] = expected_schedule(weights[[sc]], N)
  }

  say("\n== Part 4: P-spline bootstrap, smoothing parameter fixed vs re-selected (",
      n_rep, " replications, ", n_boot, " bootstrap samples) ==")
  res = par_lapply(split(check, seq_len(nrow(check))), function(tk) {
    set.seed(seed0 + tk$id)                       # same data as 01_simulation.R
    sch = schedules[[paste(tk$scenario, tk$N)]]
    D = rpois(nrow(sch), sch$m * sch$exposure)
    E = sch$exposure
    fit0 = mode_pspline(sch$age, D, E, start_age, max_age)
    set.seed(seed0 + 2e6 + tk$id)
    b_fixed = b_resel = rep(NA_real_, n_boot)
    for (b in seq_len(n_boot)) {
      Db = rpois(length(D), D)
      b_fixed[b] = tryCatch(mode_pspline(sch$age, Db, E, start_age, max_age,
                                         log10_lambdas = fit0$log10_lambda)$mode,
                            error = function(e) NA_real_)
      b_resel[b] = tryCatch(mode_pspline(sch$age, Db, E, start_age, max_age)$mode,
                            error = function(e) NA_real_)
    }
    tm = true_mode[[tk$scenario]]
    do.call(rbind, lapply(c("fixed", "reselected"), function(v) {
      x = if (v == "fixed") b_fixed else b_resel
      q = quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      data.frame(tk, smoothing = v, cover = q[1] <= tm & tm <= q[2], width = q[2] - q[1])
    }))
  }, n_cores, "P-spline bootstrap")
  res = do.call(rbind, res)
  write.csv(res, file.path(out_dir, "pspline_bootstrap_reselect.csv"), row.names = FALSE)
  g = aggregate(cbind(cover, width) ~ scenario + N + smoothing, data = res, FUN = mean)
  g$mc_se = sqrt(g$cover * (1 - g$cover) / n_rep)
  print(g)
  for (i in seq_len(nrow(g))) {
    say(sprintf("scenario %s, N = %g, smoothing %s: coverage %s (MC s.e. %s), mean width %s",
                g$scenario[i], g$N[i], g$smoothing[i], r3(g$cover[i]), r3(g$mc_se[i]),
                r2(g$width[i])))
  }
}

writeLines(log_env$lines, file.path(out_dir, "review_numbers.txt"))
cat("\nSaved to", out_dir, "\n")
