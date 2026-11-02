# =============================================================================
# 07_text_numbers.R
# Recomputes every number quoted in the text of the manuscript (abstract,
# Sections 3 and 4, conclusions) from the saved results. Run after
# 01_simulation.R, 01b_simulation_bootstrap.R, 04_empirical.R and
# 05_figures_empirical.R. The numbers of
# Appendix D are printed by 06_laplace_vs_mcmc.R.
#
# Output: results/text_numbers.txt (also printed).
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")

out_file = file.path("results", "text_numbers.txt")
log_env = new.env()
log_env$lines = c()
say = function(...) {
  txt = paste0(...)
  log_env$lines = c(log_env$lines, txt)
  cat(txt, "\n")
}
r2 = function(x) round(x, 2)
r3 = function(x) round(x, 3)

# ---- Section 3: Monte Carlo experiment --------------------------------------

tr = readRDS(file.path("results", "simulation", "sim_truth.rds"))
s = read.csv(file.path("results", "simulation", "sim_summary.csv"))
v = function(var, m, sc, N) s[[var]][s$method == m & s$scenario == sc & s$N == N]

say("== Section 3, design ==")
for (k in seq_len(nrow(tr$scenarios))) {
  t = tr$truth[[k]]
  say(sprintf("scenario %s (w = %.2f): peaks at %s; lower/main peak height = %s",
              tr$scenarios$scenario[k], tr$scenarios$w[k],
              paste(r2(t$peaks$x), collapse = ", "), r2(tr$scenarios$peak_ratio[k])))
}
say("replications: ", tr$settings$n_rep, "; posterior draws: ", tr$settings$S)

say("\n== Section 3, choice of posterior summary (MAE, years) ==")
for (sc in c("a", "b", "c", "d")) {
  for (N in c(1e3, 1e4)) {
    say(sprintf("scenario %s, N = %g: RW2 mean %s, median %s, mode %s", sc, N,
                r2(v("mae", "rw2_mean", sc, N)), r2(v("mae", "rw2_median", sc, N)),
                r2(v("mae", "rw2_mode", sc, N))))
  }
}
say(sprintf("scenario d, N = 1000, bias: RW2 mean %s, median %s",
            r2(v("bias", "rw2_mean", "d", 1e3)), r2(v("bias", "rw2_median", "d", 1e3))))

say("\n== Section 3, RW2 vs raw posterior ==")
ratio = with(s[s$method %in% c("raw_mean", "rw2_mean"), ],
             tapply(mae, list(paste(scenario, N), method), identity))
ratio = ratio[, "raw_mean"] / ratio[, "rw2_mean"]
say("MAE ratio raw/RW2 (posterior means), all cells: ", r2(min(ratio)), " to ", r2(max(ratio)))
say("  unimodal (a): ", paste(r2(ratio[grep("^a ", names(ratio))]), collapse = ", "))
say("  secondary peak (b-d), excluding d at N = 1000: ",
    paste(r2(range(ratio[grep("^[bcd] ", names(ratio))][names(ratio)[grep("^[bcd] ", names(ratio))] != "d 1000"])), collapse = " to "))
say(sprintf("example a, N = 1e4: RW2 %s vs raw %s", r2(v("mae", "rw2_mean", "a", 1e4)),
            r2(v("mae", "raw_mean", "a", 1e4))))
say(sprintf("d, N = 1000: raw mean MAE %s, RW2 mean MAE %s", r2(v("mae", "raw_mean", "d", 1e3)),
            r2(v("mae", "rw2_mean", "d", 1e3))))
say(sprintf("raw mean bias, a, N = 1000: %s", r2(v("bias", "raw_mean", "a", 1e3))))
b = s[s$method == "rw2_mean", ]
say("RW2 mean bias, all cells except d at N = 1000: max |bias| = ",
    r2(max(abs(b$bias[!(b$scenario == "d" & b$N == 1e3)]))))
cov = s[s$method %in% c("raw_mean", "rw2_mean"), ]
say("coverage raw: ", paste(r3(range(cov$coverage[cov$method == "raw_mean"])), collapse = " to "),
    "; RW2: ", paste(r3(range(cov$coverage[cov$method == "rw2_mean"])), collapse = " to "))
say(sprintf("width at N = 1e4, unimodal: RW2 %s vs raw %s", r2(v("mean_width", "rw2_mean", "a", 1e4)),
            r2(v("mean_width", "raw_mean", "a", 1e4))))
say("width at N = 1e4, b-d: RW2 ", paste(r2(sapply(c("b", "c", "d"), function(sc) v("mean_width", "rw2_mean", sc, 1e4))), collapse = ", "),
    "; raw ", paste(r2(sapply(c("b", "c", "d"), function(sc) v("mean_width", "raw_mean", sc, 1e4))), collapse = ", "))

say("\n== Section 3, comparison with other estimators ==")
for (N in c(1e3, 1e6)) {
  say(sprintf("unimodal, N = %g: RW2 mean %s, P-spline %s, Gompertz %s", N,
              r2(v("mae", "rw2_mean", "a", N)), r2(v("mae", "pspline", "a", N)),
              r2(v("mae", "gompertz", "a", N))))
}
say("Gompertz bias b, c, d (range over N): ",
    paste(sapply(c("b", "c", "d"), function(sc) paste(r2(range(s$bias[s$method == "gompertz" & s$scenario == sc])), collapse = "/")), collapse = "; "))
say("DDT bias, unimodal: ", paste(r2(range(s$bias[s$method == "ddt" & s$scenario == "a"])), collapse = " to "))
say("N = 1e4, b-d: P-spline ", paste(r2(sapply(c("b", "c", "d"), function(sc) v("mae", "pspline", sc, 1e4))), collapse = ", "),
    "; RW2 mean ", paste(r2(sapply(c("b", "c", "d"), function(sc) v("mae", "rw2_mean", sc, 1e4))), collapse = ", "))
say(sprintf("N = 1000, b: P-spline %s vs RW2 mean %s", r2(v("mae", "pspline", "b", 1e3)),
            r2(v("mae", "rw2_mean", "b", 1e3))))

say("DDT bias, all scenarios: ", paste(r2(range(s$bias[s$method == "ddt"])), collapse = " to "),
    "; DDT MAE by scenario at N = 1e4: ",
    paste(r2(sapply(c("a", "b", "c", "d"), function(sc) v("mae", "ddt", sc, 1e4))), collapse = ", "))
say("MAE at N = 1e4, all methods (a, b, c, d):")
for (m in c("rw2_mean", "raw_mean", "pspline", "gompertz", "ddt")) {
  say(sprintf("  %s: %s", m, paste(r2(sapply(c("a", "b", "c", "d"), function(sc) v("mae", m, sc, 1e4))), collapse = ", ")))
}

boot_file = file.path("results", "simulation", "sim_bootstrap_summary.csv")
if (file.exists(boot_file)) {
  bs = read.csv(boot_file)
  w = function(var, m, sc, N) bs[[var]][bs$method == m & bs$scenario == sc & bs$N == N]
  say("\n== Section 3, bootstrap intervals of the benchmarks ==")
  for (m in c("gompertz", "pspline", "ddt")) {
    for (sc in c("a", "b", "c", "d")) {
      say(sprintf("%s, scenario %s: coverage %s; width %s (N = 1e3 ... 1e6); no interval %s", m, sc,
                  paste(r3(sapply(sort(unique(bs$N)), function(N) w("coverage", m, sc, N))), collapse = ", "),
                  paste(r2(sapply(sort(unique(bs$N)), function(N) w("mean_width", m, sc, N))), collapse = ", "),
                  sum(bs$n_no_interval[bs$method == m & bs$scenario == sc])))
    }
  }
  say("RW2 for reference, coverage by scenario at N = 1e3 ... 1e6:")
  for (sc in c("a", "b", "c", "d")) {
    say(sprintf("  rw2 %s: coverage %s; width %s", sc,
                paste(r3(sapply(sort(unique(s$N)), function(N) v("coverage", "rw2_mean", sc, N))), collapse = ", "),
                paste(r2(sapply(sort(unique(s$N)), function(N) v("mean_width", "rw2_mean", sc, N))), collapse = ", ")))
  }
} else {
  say("\n(sim_bootstrap_summary.csv not found: run 01b_simulation_bootstrap.R)")
}

# ---- Section 4: empirical application (women) -------------------------------

pe = readRDS(file.path("results", "empirical", "point_estimates.rds"))
ps = readRDS(file.path("results", "empirical", "posterior_summary.rds"))
st = read.csv(file.path("results", "empirical", "start_age_consistency.csv"))
base = min(ps$start_age)
s0 = ps[ps$start_age == base, ]
s0$width = s0$q975 - s0$q025

say("\n== Section 4, data ==")
say("countries: ", paste(unique(pe$country), collapse = ", "), "; sexes: ",
    paste(unique(pe$sex), collapse = ", "))
say("years: ", paste(range(pe$year), collapse = " to "), "; population-years: ",
    nrow(unique(pe[, c("country", "sex", "year")])))
say("last year by country: ", paste(names(tapply(pe$year, pe$country, max)),
                                     tapply(pe$year, pe$country, max), collapse = ", "))
say("starting ages: ", paste(sort(unique(pe$start_age)), collapse = ", "),
    "; reference: ", base)
lowest = s0[which.min(s0$mean), ]
say("lowest RW2 posterior mean at the reference starting age: ", r2(lowest$mean),
    " (", lowest$country, ", ", lowest$year, ")")

say(sprintf("\n== Section 4, trends (RW2 posterior mean, starting age %d) ==", base))
for (cc in unique(s0$country)) {
  z = s0[s0$country == cc, ]
  first = z[z$year == min(z$year), ]
  last = z[z$year == max(z$year), ]
  say(sprintf("%s: %d %s (%s-%s) -> %d %s (%s-%s); change %s", cc, first$year, r2(first$mean),
              r2(first$q025), r2(first$q975), last$year, r2(last$mean), r2(last$q025),
              r2(last$q975), r2(last$mean - first$mean)))
}
say("width of 95% interval: median ", r2(median(s0$width)), ", range ",
    paste(r2(range(s0$width)), collapse = " to "))
wd = tapply(s0$width, s0$country, median)
say("median width by country: ", paste(names(wd), r2(wd), collapse = ", "))
say("median P(more than one major local mode): ", r3(median(s0$p_multimodal)))

w = reshape(pe[pe$start_age == base, c("country", "sex", "year", "method", "estimate")],
            idvar = c("country", "sex", "year"), timevar = "method", direction = "wide")
w = merge(w, s0[, c("country", "sex", "year", "q025", "q975")])
for (m in c("pspline", "ddt", "gompertz")) {
  x = w[[paste0("estimate.", m)]]
  say(sprintf("%s: mean difference from RW2 mean %s, mean absolute difference %s, inside RW2 interval %s%%",
              m, r2(mean(x - w$estimate.rw2_mean)), r2(mean(abs(x - w$estimate.rw2_mean))),
              round(100 * mean(x >= w$q025 & x <= w$q975), 1)))
}

say(sprintf("\n== Section 4, 2019 vs 2020 (RW2, starting age %d) ==", base))
for (cc in unique(s0$country)) {
  z = s0[s0$country == cc & s0$year %in% c(2019, 2020), ]
  if (nrow(z) < 2) next
  z = z[order(z$year), ]
  say(sprintf("%s: 2019 %s (%s-%s); 2020 %s (%s-%s); intervals overlap: %s", cc,
              r2(z$mean[1]), r2(z$q025[1]), r2(z$q975[1]), r2(z$mean[2]), r2(z$q025[2]),
              r2(z$q975[2]), z$q025[2] <= z$q975[1] & z$q975[2] >= z$q025[1]))
}

say("\n== Section 4, starting age (also Table 1) ==")
a = st[st$country == "All", ]
for (i in seq_len(nrow(a))) {
  say(sprintf("%s, %s: mean %s, mean abs %s, max abs %s", a$method[i], a$comparison[i],
              r3(a$mean_change[i]), r3(a$mean_abs_change[i]), r3(a$max_abs_change[i])))
}
last_cmp = tail(sort(unique(st$comparison)), 1)
b = st[st$country != "All" & st$comparison == last_cmp, ]
for (m in c("rw2_mean", "pspline", "gompertz")) {
  z = b[b$method == m, ]
  say(sprintf("%s, %s, mean abs change by country: %s", m, last_cmp,
              paste(z$country, r2(z$mean_abs_change), collapse = ", ")))
}

writeLines(log_env$lines, out_file)
cat("\nSaved to", out_file, "\n")
