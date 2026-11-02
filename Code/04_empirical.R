# =============================================================================
# 04_empirical.R
# Empirical application (Section 4): modal age at death of women in Denmark,
# France, Japan and Sweden by year, from 1960 to the last year in the data
# files, with the RW2 posterior (mean, median, mode, 95% credible interval)
# and the Gompertz, P-spline and DDT benchmarks, at starting ages 40, 50, 60
# and 70.
#
# Data: data/<code>/STATS/Deaths_1x1.txt and Exposures_1x1.txt from the
# Human Mortality Database (www.mortality.org); codes DNK, FRATNP (France,
# total population), JPN and SWE.
#
# Outputs (results/empirical/):
#   point_estimates.rds / .csv    one row per country x sex x year x start
#                                 age x method (estimate; 95% interval for
#                                 the RW2 rows)
#   posterior_summary.rds / .csv  RW2 posterior summaries
#   posterior_int.rds             RW2 P(mode = x | data) on integer ages
#
# Uses forked parallelism (macOS / Linux); a few minutes on 8 cores.
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")
for (f in c("utils.R", "mode_posterior.R", "mode_posterior_rw.R", "mode_gompertz.R",
            "mode_pspline.R", "mode_ddt.R")) source(f)

# ---- Settings ---------------------------------------------------------------

data_dir = "data"
countries = names(hmd_countries)             # DNK, FRATNP, JPN, SWE (utils.R)
first_year = 1960
last_year = as.integer(Sys.getenv("MODE_LAST_YEAR", 9999))
start_ages = c(40, 50, 60, 70)
max_age = 110
S = 1e4
seed0 = 20261107
n_cores = get_cores()

out_dir = file.path("results", "empirical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hmd = do.call(rbind, lapply(countries, load_country, data_dir = data_dir))
hmd = hmd[hmd$year >= first_year & hmd$year <= last_year, ]
tasks = unique(hmd[, c("country", "sex", "year")])
tasks = tasks[order(tasks$country, tasks$sex, tasks$year), ]
tasks$id = seq_len(nrow(tasks))
rownames(tasks) = NULL

# ---- One population-year ----------------------------------------------------

run_one = function(k) {
  tk = tasks[k, ]
  set.seed(seed0 + tk$id)
  y = hmd[hmd$country == tk$country & hmd$sex == tk$sex & hmd$year == tk$year, ]
  key = data.frame(country = tk$country, sex = tk$sex, year = tk$year)
  safe = function(expr) tryCatch(expr, error = function(e) NA_real_)
  pp = mode_posterior_rw(y$age, y$deaths, y$exposure, start_ages, max_age, S = S)
  pe = list()
  ps = list()
  pi = list()
  for (st in start_ages) {
    p = pp[[as.character(st)]]
    i = as.character(st)
    pe[[i]] = data.frame(key, start_age = st,
      method = c("rw2_mean", "rw2_median", "rw2_mode", "gompertz", "pspline", "ddt"),
      estimate = c(p$mean, p$median, p$mode,
                   safe(mode_gompertz(y$age, y$deaths, y$exposure, st, max_age)$mode),
                   safe(mode_pspline(y$age, y$deaths, y$exposure, st, max_age)$mode),
                   safe(mode_ddt(y$age, y$deaths, y$exposure, st, max_age)$mode)),
      lower = c(rep(p$q025, 3), NA, NA, NA),
      upper = c(rep(p$q975, 3), NA, NA, NA))
    ps[[i]] = data.frame(key, start_age = st, map = p$map, p_map = p$p_map,
      mean = p$mean, median = p$median, mode = p$mode, q025 = p$q025, q975 = p$q975,
      cred_lo = min(p$cred_set), cred_hi = max(p$cred_set), cred_size = p$cred_size,
      cred_contiguous = p$cred_contiguous, exp_local_modes = p$exp_local_modes,
      exp_major_modes = p$exp_major_modes, p_multimodal = p$p_multimodal,
      post_log10_lambda = p$post_log10_lambda, last_age = p$last_age)
    post = p$post[p$post$prob > 0, ]
    pi[[i]] = data.frame(key, start_age = st, age = post$age, prob = post$prob)
  }
  list(pe = do.call(rbind, pe), ps = do.call(rbind, ps), pi = do.call(rbind, pi))
}

# ---- Run --------------------------------------------------------------------

cat("Running", nrow(tasks), "population-years on", n_cores, "cores\n")
t0 = Sys.time()
res = par_lapply(seq_len(nrow(tasks)), run_one, n_cores, "HMD")
cat("Done in", format(Sys.time() - t0), "\n")

failed = sapply(res, inherits, "try-error")
if (any(failed)) {
  warning(sum(failed), " population-years failed: ",
          paste(apply(tasks[failed, 1:3], 1, paste, collapse = " "), collapse = "; "))
  res = res[!failed]
}
grab = function(name) do.call(rbind, lapply(res, `[[`, name))
point = grab("pe")
summ = grab("ps")
saveRDS(point, file.path(out_dir, "point_estimates.rds"))
write.csv(point, file.path(out_dir, "point_estimates.csv"), row.names = FALSE)
saveRDS(summ, file.path(out_dir, "posterior_summary.rds"))
write.csv(summ, file.path(out_dir, "posterior_summary.csv"), row.names = FALSE)
saveRDS(grab("pi"), file.path(out_dir, "posterior_int.rds"))
cat("Saved results to", out_dir, "\n")
