# =============================================================================
# 06_tables_appendix.R
# Appendix tables with the full Monte Carlo results (run 01_simulation_rw.R
# first). Writes LaTeX tables to Manuscript/tables/ (input in main.tex) and
# copies to random_walk/Results/simulation/tables/.
#
#   app_mc_mae.tex        mean absolute error, all methods
#   app_mc_bias.tex       bias, all methods
#   app_mc_sd.tex         standard deviation of the estimates, all methods
#   app_mc_coverage.tex   coverage and width of the 95% credible interval,
#                         raw and RW2 posteriors
# =============================================================================

code_dir = Sys.getenv("MODE_RW_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk")
setwd(code_dir)

in_dir = file.path("Results", "simulation")
tab_dir = file.path(in_dir, "tables")
ms_dir = Sys.getenv("MODE_MS_TAB_DIR", file.path("..", "..", "..", "Manuscript", "tables"))
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)

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

save_tab = function(lines, name) {
  writeLines(lines, file.path(tab_dir, name))
  writeLines(lines, file.path(ms_dir, name))
}

method_table = function(var, caption, label, digits = 2) {
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
save_tab(tidy(method_table("mae", paste("Monte Carlo experiment: mean absolute error (years) of the estimated modal age.", cap_note), "tab:app_mae")),
         "app_mc_mae.tex")
save_tab(tidy(method_table("bias", paste("Monte Carlo experiment: bias (years) of the estimated modal age.", cap_note), "tab:app_bias")),
         "app_mc_bias.tex")
save_tab(tidy(method_table("sd", paste("Monte Carlo experiment: standard deviation (years) of the estimated modal age across replications.", cap_note), "tab:app_sd")),
         "app_mc_sd.tex")

# coverage and width
head = c("\\begin{table}[htb!]",
         "\\centering",
         sprintf("\\caption{Monte Carlo experiment: coverage and mean width (years) of the 95\\%% credible interval of $\\tilde{M}$. %s}", cap_note),
         "\\label{tab:app_coverage}",
         "\\small",
         "\\begin{tabular}{llrrrr}",
         "\\hline",
         " & & \\multicolumn{2}{c}{Coverage} & \\multicolumn{2}{c}{Width} \\\\",
         "Scenario & $N$ & Raw & RW2 & Raw & RW2 \\\\",
         "\\hline")
body = c()
for (sc in scen$scenario) {
  for (N in sort(unique(summ$N))) {
    z = summ[summ$scenario == sc & summ$N == N, ]
    g = function(m, v) z[[v]][z$method == m]
    first = N == min(summ$N)
    body = c(body, paste0(if (first) scen_lab[sc] else "", " & ", N_lab(N), " & ",
                          fmt(g("raw_mean", "coverage"), 3), " & ", fmt(g("rw2_mean", "coverage"), 3),
                          " & ", fmt(g("raw_mean", "mean_width"), 2), " & ",
                          fmt(g("rw2_mean", "mean_width"), 2), " \\\\"))
  }
  if (sc != tail(scen$scenario, 1)) body = c(body, "")
}
save_tab(tidy(c(head, body, "\\hline", "\\end{tabular}", "\\end{table}")), "app_mc_coverage.tex")
cat("Tables saved to", tab_dir, "and", ms_dir, "\n")
