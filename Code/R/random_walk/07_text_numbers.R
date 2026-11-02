# =============================================================================
# 07_text_numbers.R
# Recomputes every number quoted in the text of the manuscript (abstract,
# Sections 3 and 4, conclusions) from the saved results. Run after
# 01_simulation_rw.R, 02_empirical_rw.R and 04_figures_empirical.R. The
# numbers of Appendix D (Laplace vs exact sampling) are printed by
# 05_laplace_vs_mcmc.R.
#
# Output: Results/text_numbers.txt (also printed).
# Run from Code/R/random_walk.
# =============================================================================

code_dir = Sys.getenv("MODE_RW_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk")
setwd(code_dir)

out_file = file.path("Results", "text_numbers.txt")
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

tr = readRDS(file.path("Results", "simulation", "sim_truth.rds"))
s = read.csv(file.path("Results", "simulation", "sim_summary.csv"))
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

# ---- Section 4: empirical application ---------------------------------------

pe = readRDS(file.path("Results", "empirical", "point_estimates.rds"))
ps = readRDS(file.path("Results", "empirical", "posterior_summary.rds"))
st = read.csv(file.path("Results", "empirical", "start_age_consistency.csv"))
base = 50
s50 = ps[ps$start_age == base, ]
s50$width = s50$q975 - s50$q025

say("\n== Section 4, data ==")
say("years: ", paste(range(pe$year), collapse = " to "), "; population-years: ",
    nrow(unique(pe[, c("country", "sex", "year")])))
say("last year by country: ", paste(names(tapply(pe$year, pe$country, max)),
                                     tapply(pe$year, pe$country, max), collapse = ", "))
say("lowest RW2 posterior mean at starting age 50: ", r2(min(s50$mean)))

say("\n== Section 4, trends (RW2 posterior mean, starting age 50) ==")
gain = do.call(rbind, lapply(split(s50, list(s50$country, s50$sex)), function(z) {
  data.frame(country = z$country[1], sex = z$sex[1], first = z$mean[z$year == 1960],
             last = z$mean[z$year == max(z$year)], last_year = max(z$year))
}))
gain$gain = gain$last - gain$first
for (sx in c("Female", "Male")) {
  g = gain[gain$sex == sx, ]
  say(sprintf("%s: increase %s (%s) to %s (%s)", sx, r2(min(g$gain)), g$country[which.min(g$gain)],
              r2(max(g$gain)), g$country[which.max(g$gain)]))
}
say("width of 95% interval: median ", r2(median(s50$width)), ", range ",
    paste(r2(range(s50$width)), collapse = " to "))
wd = tapply(s50$width, list(s50$country, s50$sex), median)
say("median width by country and sex:")
for (i in seq_len(nrow(wd))) say(sprintf("  %s: women %s, men %s", rownames(wd)[i], r2(wd[i, 1]), r2(wd[i, 2])))
say("median P(more than one major local mode): ", r3(median(s50$p_multimodal)))

w = reshape(pe[pe$start_age == base, c("country", "sex", "year", "method", "estimate")],
            idvar = c("country", "sex", "year"), timevar = "method", direction = "wide")
w = merge(w, s50[, c("country", "sex", "year", "q025", "q975")])
for (m in c("pspline", "ddt", "gompertz")) {
  x = w[[paste0("estimate.", m)]]
  say(sprintf("%s: mean difference from RW2 mean %s, mean absolute difference %s, inside RW2 interval %s%%",
              m, r2(mean(x - w$estimate.rw2_mean)), r2(mean(abs(x - w$estimate.rw2_mean))),
              round(100 * mean(x >= w$q025 & x <= w$q975), 1)))
}

say("\n== Section 4, changes around 2020 (RW2, starting age 50) ==")
for (cs in list(c("US", "Male"), c("US", "Female"), c("Italy", "Male"), c("France", "Male"),
                c("Denmark", "Female"), c("Denmark", "Male"))) {
  z = s50[s50$country == cs[1] & s50$sex == cs[2] & s50$year %in% c(2019, 2020), ]
  z = z[order(z$year), ]
  say(sprintf("%s %s: 2019 %s (%s-%s); 2020 %s (%s-%s); intervals overlap: %s", cs[1], cs[2],
              r2(z$mean[1]), r2(z$q025[1]), r2(z$q975[1]), r2(z$mean[2]), r2(z$q025[2]),
              r2(z$q975[2]), z$q025[2] <= z$q975[1] & z$q975[2] >= z$q025[1]))
}
z = s50[s50$country == "Italy" & s50$sex == "Male" & s50$year %in% c(1967, 1969), ]
say(sprintf("Italy men 1969: %s (%s-%s); 1967: %s", r2(z$mean[z$year == 1969]),
            r2(z$q025[z$year == 1969]), r2(z$q975[z$year == 1969]), r2(z$mean[z$year == 1967])))

say("\n== Section 4, starting age (also Table 1) ==")
a = st[st$country == "All", ]
for (i in seq_len(nrow(a))) {
  say(sprintf("%s, %s: mean %s, mean abs %s, max abs %s", a$method[i], a$comparison[i],
              r3(a$mean_change[i]), r3(a$mean_abs_change[i]), r3(a$max_abs_change[i])))
}
b = st[st$country != "All" & st$comparison == "70 vs 50", ]
for (m in c("rw2_mean", "pspline", "gompertz")) {
  z = b[b$method == m, ]
  say(sprintf("%s, 70 vs 50, mean abs change by country: %s (%s) to %s (%s)", m,
              r2(min(z$mean_abs_change)), z$country[which.min(z$mean_abs_change)],
              r2(max(z$mean_abs_change)), z$country[which.max(z$mean_abs_change)]))
}

writeLines(log_env$lines, out_file)
cat("\nSaved to", out_file, "\n")
