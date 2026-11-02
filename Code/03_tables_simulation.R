# =============================================================================
# 03_tables_simulation.R
# Appendix tables with the full Monte Carlo results (Appendix E). Run
# 01_simulation.R first. Tables are saved to results/simulation/tables/ and
# to the manuscript (see utils.R).
#
#   app_mc_mae.tex        mean absolute error, all methods
#   app_mc_bias.tex       bias, all methods
#   app_mc_sd.tex         standard deviation of the estimates, all methods
#   app_mc_coverage.tex   coverage and width of the 95% intervals: credible
#                         (raw, RW2) and bootstrap (Gompertz, P-spline, DDT)
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")
source("utils.R")

in_dir = file.path("results", "simulation")
summ = read.csv(file.path(in_dir, "sim_summary.csv"))
tr = readRDS(file.path(in_dir, "sim_truth.rds"))
start = 40
summ = summ[summ$start_age == start, ]

methods = c(rw2_mean = "Mean", rw2_median = "Median", rw2_mode = "Mode",
            raw_mean = "Mean", raw_median = "Median", raw_mode = "Mode",
            gompertz = "Gompertz", pspline = "P-spline", ddt = "DDT")
scen = tr$scenarios
scen_lab = ifelse(scen$w == 1, "(a) unimodal",
                  sprintf("(%s) %.0f\\%% earlier", scen$scenario, 100 * (1 - scen$w)))
names(scen_lab) = scen$scenario
N_lab = function(N) sprintf("$10^%d$", round(log10(N)))

fmt = function(x, digits) {
  out = formatC(x, format = "f", digits = digits)
  out = ifelse(substr(out, 1, 1) == "-", paste0("$-$", substring(out, 2)), out)
  ifelse(is.finite(x), out, "--")
}

save_tab = function(lines, name) save_table(lines, name, file.path("simulation", "tables"))

method_table = function(var, caption, label, digits = 4) {
  head = c("\\begin{table}[htb!]",
           "\\centering",
           sprintf("\\caption{%s}", caption),
           sprintf("\\label{%s}", label),
           "\\scriptsize",
           "\\setlength{\\tabcolsep}{4pt}",
           "\\begin{tabular}{llrrrrrrrrr}",
           "\\hline",
           " & & \\multicolumn{3}{c}{RW2 posterior} & \\multicolumn{3}{c}{Raw posterior} & & & \\\\",
           "Scenario & $N$ & Mean & Median & Mode & Mean & Median & Mode & Gompertz & P-spline & DDT \\\\",
           "\\hline")
  body = c()
  for (sc in scen$scenario) {
    for (N in sort(unique(summ$N))) {
      z = summ[summ$scenario == sc & summ$N == N, ]
      v = sapply(names(methods), function(m) z[[var]][z$method == m])
      first = N == min(summ$N)
      body = c(body, paste0(if (first) scen_lab[sc] else "", " & ", N_lab(N), " & ",
                            paste(fmt(v, digits), collapse = " & "), " \\\\"))
    }
    if (sc != tail(scen$scenario, 1)) body = c(body, "")
  }
  c(head, body, "\\hline", "\\end{tabular}", "\\end{table}")
}

# turn the empty separator lines into extra space after the preceding row
tidy = function(lines) {
  out = c()
  for (ln in lines) {
    if (ln == "" && length(out) > 0) {
      out[length(out)] = sub("\\\\\\\\$", "\\\\\\\\[2pt]", out[length(out)])
    } else {
      out = c(out, ln)
    }
  }
  out
}

cap_note = sprintf("Starting age %d, %d replications, true mode 90.", start, tr$settings$n_rep)
save_tab(tidy(method_table("mae", paste("Simulation study: mean absolute error (years) of the estimated modal age.", cap_note), "tab:app_mae")),
         "app_mc_mae.tex")
save_tab(tidy(method_table("bias", paste("Simulation study: bias (years) of the estimated modal age.", cap_note), "tab:app_bias")),
         "app_mc_bias.tex")
save_tab(tidy(method_table("sd", paste("Simulation study: standard deviation (years) of the estimated modal age across replications.", cap_note), "tab:app_sd")),
         "app_mc_sd.tex")

# coverage and width: credible intervals (raw, RW2) and bootstrap intervals
boot_file = file.path(in_dir, "sim_bootstrap_summary.csv")
bs = if (file.exists(boot_file)) read.csv(boot_file) else summ[0, ]
cmeth = c("raw_mean", "rw2_mean", "gompertz", "pspline", "ddt")
clab = c("Raw", "RW2", "Gomp.", "P-spl.", "DDT")
head = c("\\begin{table}[htb!]",
         "\\centering",
         sprintf("\\caption{Simulation study: coverage and mean width (years) of the 95\\%% intervals for the modal age: credible intervals of the raw and RW2 posteriors (posterior of $\\tilde{M}$) and bootstrap percentile intervals of the Gompertz, P-spline and DDT estimators. %s}", cap_note),
         "\\label{tab:app_coverage}",
         "\\footnotesize",
         "\\setlength{\\tabcolsep}{4pt}",
         "\\begin{tabular}{ll*{10}{r}}",
         "\\hline",
         " & & \\multicolumn{5}{c}{Coverage} & \\multicolumn{5}{c}{Width} \\\\",
         paste0("Scenario & $N$ & ", paste(c(clab, clab), collapse = " & "), " \\\\"),
         "\\hline")
body = c()
for (sc in scen$scenario) {
  for (N in sort(unique(summ$N))) {
    z = rbind(summ[summ$scenario == sc & summ$N == N & summ$method %in% c("raw_mean", "rw2_mean"),
                   c("method", "coverage", "mean_width")],
              bs[bs$scenario == sc & bs$N == N, c("method", "coverage", "mean_width")])
    g = function(m, v) {
      x = as.numeric(z[[v]][z$method == m])
      if (length(x) == 0) NA_real_ else x[1]
    }
    first = N == min(summ$N)
    body = c(body, paste0(if (first) scen_lab[sc] else "", " & ", N_lab(N), " & ",
                          paste(fmt(vapply(cmeth, g, numeric(1), v = "coverage"), 4), collapse = " & "), " & ",
                          paste(fmt(vapply(cmeth, g, numeric(1), v = "mean_width"), 4), collapse = " & "),
                          " \\\\"))
  }
  if (sc != tail(scen$scenario, 1)) body = c(body, "")
}
save_tab(tidy(c(head, body, "\\hline", "\\end{tabular}", "\\end{table}")), "app_mc_coverage.tex")
cat("Tables saved to results/simulation/tables\n")
