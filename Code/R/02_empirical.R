# =============================================================================
# 02_empirical.R
# Modal age at death from HMD data with the four methods, at each starting
# age, for every country, sex and year from first_year to the last year
# available in each country's files.
#
# Outputs (Results/empirical/):
#   point_estimates.rds / .csv   one row per country x sex x year x start
#                                age x method, with method diagnostics
#   posterior_summary.rds / .csv summaries of the posterior of the mode
#   posterior_int.rds            P(mode = x | data) on integer ages
#                                (rows with probability 0 dropped)
#   posterior_cont.rds           posterior of the continuous (parabola)
#                                mode on a 0.1-year grid
#
# The posterior files keep the full distributions, so comparisons such as
# P(M_female > M_male) or P(M_t > M_{t-1}) can be computed afterwards
# (years and sexes are independent populations).
#
# Run from Code/R. Uses forked parallelism (macOS / Linux).
# =============================================================================

code_dir = Sys.getenv("MODE_CODE_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R")
setwd(code_dir)
for (f in c("utils.R", "mode_posterior.R", "mode_gompertz.R",
            "mode_pspline.R", "mode_ddt.R")) source(f)

# ---- Settings ---------------------------------------------------------------

data_dir = file.path("..", "data")
countries = c("Denmark", "France", "Italy", "Japan", "Netherlands", "US")
first_year = 1960
last_year = as.integer(Sys.getenv("MODE_LAST_YEAR", 9999))   # for testing
start_ages = c(30, 40, 50, 60)
max_age = 110
S = 1e4
seed0 = 20261102
n_cores = get_cores()

out_dir = file.path("Results", "empirical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data -------------------------------------------------------------------

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
  pe = list()
  ps = list()
  pi = list()
  pc = list()
  # one set of posterior draws serves every starting age
  pp = mode_posterior(y$age, y$deaths, y$exposure, start_ages, max_age, S = S)
  for (st in start_ages) {
    p = pp[[as.character(st)]]
    g = tryCatch(mode_gompertz(y$age, y$deaths, y$exposure, st, max_age),
                 error = function(e) list(mode = NA_real_, a = NA_real_, b = NA_real_))
    sp = tryCatch(mode_pspline(y$age, y$deaths, y$exposure, st, max_age),
                  error = function(e) list(mode = NA_real_, log10_lambda = NA_real_,
                                           boundary = NA))
    dd = tryCatch(mode_ddt(y$age, y$deaths, y$exposure, st, max_age),
                  error = function(e) list(mode = NA_real_, n_candidates = NA_integer_,
                                           n_conditions = NA_integer_))
    i = as.character(st)
    pe[[i]] = data.frame(key, start_age = st,
      method = c("posterior_mean", "posterior_median", "posterior_mode",
                 "gompertz", "pspline", "ddt"),
      estimate = c(p$mean, p$median, p$mode, g$mode, sp$mode, dd$mode),
      gompertz_a = c(NA, NA, NA, g$a, NA, NA),
      gompertz_b = c(NA, NA, NA, g$b, NA, NA),
      pspline_log10_lambda = c(NA, NA, NA, NA, sp$log10_lambda, NA),
      pspline_boundary = c(NA, NA, NA, NA, sp$boundary, NA),
      ddt_n_candidates = c(NA, NA, NA, NA, NA, dd$n_candidates),
      ddt_n_conditions = c(NA, NA, NA, NA, NA, dd$n_conditions))
    ps[[i]] = data.frame(key, start_age = st, map = p$map, p_map = p$p_map,
      mean = p$mean, median = p$median, mode = p$mode, q025 = p$q025, q975 = p$q975,
      cred_lo = min(p$cred_set), cred_hi = max(p$cred_set),
      cred_size = p$cred_size, cred_contiguous = p$cred_contiguous,
      exp_local_modes = p$exp_local_modes, exp_major_modes = p$exp_major_modes,
      p_multimodal = p$p_multimodal, last_age = p$last_age)
    post = p$post[p$post$prob > 0, ]
    pi[[i]] = data.frame(key, start_age = st, age = post$age, prob = post$prob)
    pc[[i]] = data.frame(key, start_age = st, mode = p$post_cont$mode,
                         prob = p$post_cont$prob)
  }
  list(pe = do.call(rbind, pe), ps = do.call(rbind, ps),
       pi = do.call(rbind, pi), pc = do.call(rbind, pc))
}

# ---- Run --------------------------------------------------------------------

cat("Running", nrow(tasks), "population-years x", length(start_ages),
    "starting ages on", n_cores, "cores\n")
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
saveRDS(grab("pc"), file.path(out_dir, "posterior_cont.rds"))
cat("Saved results to", out_dir, "\n")
