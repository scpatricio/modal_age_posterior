# =============================================================================
# 04_figures_empirical.R
# Manuscript figures and tables for the empirical application (run
# 02_empirical_rw.R first). Every figure is saved to random_walk/Figures/
# empirical/ and to Manuscript/figures/.
#
#   emp_trends.pdf             modal age at death by year, country and sex,
#                              starting age 50: RW2 posterior mean with 95%
#                              interval, Gompertz, P-spline, DDT
#   emp_start_age.pdf          change in each method's estimate when the
#                              starting age moves from 50 to 60 and to 70
#   emp_start_age_female.pdf,  estimates by year at the three starting ages,
#   emp_start_age_male.pdf     one row per method
# Results/empirical/
#   start_age_consistency.csv  mean and mean absolute change by method,
#                              comparison and country
# Manuscript/tables/
#   emp_start_age.tex          Table 1 of the manuscript (all countries)
# =============================================================================

code_dir = Sys.getenv("MODE_RW_DIR", "~/Desktop/fertility/PDF_of_the_mode/Code/R/random_walk")
setwd(code_dir)
library(ggplot2)

in_dir = file.path("Results", "empirical")
fig_dir = file.path("Figures", "empirical")
ms_dir = Sys.getenv("MODE_MS_FIG_DIR", file.path("..", "..", "..", "Manuscript", "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)
save_both = function(p, name, width, height) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height)
  ggsave(file.path(ms_dir, name), p, width = width, height = height)
}

pe = readRDS(file.path(in_dir, "point_estimates.rds"))

method_labels = c(rw2_mean = "RW2 posterior mean", rw2_median = "RW2 posterior median",
                  rw2_mode = "RW2 posterior mode", gompertz = "Gompertz",
                  pspline = "P-spline", ddt = "DDT")
cols = setNames(c("#0072B2", "#56B4E9", "#CC79A7", "#E69F00", "#009E73", "#D55E00"),
                method_labels)
shapes = setNames(c(16, 1, 18, 17, 15, 4), method_labels)
pe$method_lab = factor(method_labels[pe$method], levels = method_labels)

theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  strip.text = element_text(face = "bold", hjust = 0),
                  legend.position = "bottom"))

# ---- Trends at starting age 50 ----------------------------------------------

base_start = 50
tm = c("rw2_mean", "gompertz", "pspline", "ddt")
d = pe[pe$start_age == base_start & pe$method %in% tm, ]
band = d[d$method == "rw2_mean", ]
p = ggplot(d, aes(year, estimate, colour = method_lab)) +
  geom_ribbon(data = band, aes(ymin = lower, ymax = upper), fill = "#0072B2",
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.4) +
  facet_grid(sex ~ country) +
  scale_colour_manual(values = cols[method_labels[tm]], name = NULL) +
  labs(x = "Year", y = "Modal age at death",
       caption = sprintf(paste("Starting age %d. Shaded: 95%% credible interval of the",
                               "RW2 posterior."), base_start))
save_both(p, "emp_trends.pdf", 12, 6)

# ---- Consistency across starting ages ---------------------------------------

key = paste(pe$country, pe$sex, pe$year, pe$method)
ref = pe[pe$start_age == base_start, ]
pe$ref = ref$estimate[match(key, paste(ref$country, ref$sex, ref$year, ref$method))]
dd = pe[pe$start_age != base_start, ]
dd$change = dd$estimate - dd$ref
dd$comparison = factor(sprintf("%d vs %d", dd$start_age, base_start),
                       levels = sprintf("%d vs %d", sort(unique(dd$start_age)), base_start))

p = ggplot(dd, aes(method_lab, change, fill = comparison)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3,
               position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ country, nrow = 2) +
  scale_fill_manual(values = c("#BBBBBB", "#555555"), name = "Starting age") +
  labs(x = NULL, y = "Change in the estimated mode (years)",
       caption = sprintf(paste("Estimate at the later starting age minus the estimate",
                               "at %d; all years and both sexes."), base_start)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p, "emp_start_age.pdf", 12, 7)

for (sx in c("Female", "Male")) {
  d = pe[pe$sex == sx, ]
  d$start_lab = factor(paste("Starting age", d$start_age))
  p = ggplot(d, aes(year, estimate, colour = start_lab)) +
    geom_line(linewidth = 0.35) +
    facet_grid(method_lab ~ country) +
    scale_colour_manual(values = c("#0072B2", "#E69F00", "#D55E00"), name = NULL) +
    labs(x = "Year", y = "Modal age at death", caption = paste0(sx, "s.")) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
  save_both(p, sprintf("emp_start_age_%s.pdf", tolower(sx)), 12, 10)
}

# ---- Table ------------------------------------------------------------------

agg = function(z) {
  data.frame(n = sum(is.finite(z$change)),
             mean_change = mean(z$change, na.rm = TRUE),
             mean_abs_change = mean(abs(z$change), na.rm = TRUE),
             max_abs_change = max(abs(z$change), na.rm = TRUE))
}
g = split(dd, list(dd$method, dd$comparison), drop = TRUE)
tab_all = do.call(rbind, lapply(g, function(z) {
  data.frame(method = z$method[1], comparison = z$comparison[1], country = "All", agg(z))
}))
g = split(dd, list(dd$method, dd$comparison, dd$country), drop = TRUE)
tab_c = do.call(rbind, lapply(g, function(z) {
  data.frame(method = z$method[1], comparison = z$comparison[1], country = z$country[1],
             agg(z))
}))
tab = rbind(tab_all, tab_c)
tab$method = factor(tab$method, levels = names(method_labels))
tab = tab[order(tab$country != "All", tab$country, tab$comparison, tab$method), ]
write.csv(tab, file.path(in_dir, "start_age_consistency.csv"), row.names = FALSE)
print(tab_all[order(tab_all$comparison, tab_all$method), ], row.names = FALSE)
cat("Figures saved to", fig_dir, "and", ms_dir, "\n")

# ---- LaTeX table (Table 1 of the manuscript) --------------------------------

ms_tab = Sys.getenv("MODE_MS_TAB_DIR", file.path("..", "..", "..", "Manuscript", "tables"))
dir.create(ms_tab, recursive = TRUE, showWarnings = FALSE)
f3 = function(x) {
  out = formatC(x, format = "f", digits = 3)
  out = ifelse(abs(x) < 0.0005, formatC(0, format = "f", digits = 3), out)
  ifelse(substr(out, 1, 1) == "-", paste0("$-$", substring(out, 2)), out)
}
row_methods = c("rw2_mean", "rw2_median", "rw2_mode", "pspline", "gompertz", "ddt")
comps = levels(dd$comparison)
rows = sapply(row_methods, function(m) {
  v = unlist(lapply(comps, function(cp) {
    z = tab_all[tab_all$method == m & tab_all$comparison == cp, ]
    c(z$mean_change, z$mean_abs_change, z$max_abs_change)
  }))
  paste0(method_labels[m], " & ", paste(f3(v), collapse = " & "), " \\\\")
})
n_pop = length(unique(paste(pe$country, pe$sex, pe$year)))
tab = c("\\begin{table}[htb!]",
        "\\centering",
        sprintf(paste0("\\caption{Change in the estimated modal age at death when the starting age moves ",
                       "from %d to %s: mean change, mean absolute change and largest absolute change ",
                       "(years), over %d countries, both sexes and all years (%d population-years).}"),
                base_start, paste(sort(unique(dd$start_age)), collapse = " and to "),
                length(unique(pe$country)), n_pop),
        "\\label{tab:start_age}",
        "\\small",
        "\\begin{tabular}{lrrrrrr}",
        "\\hline",
        paste0(" & ", paste(sprintf("\\multicolumn{3}{c}{%s}", comps), collapse = " & "), " \\\\"),
        "Method & Mean & Mean abs. & Max abs. & Mean & Mean abs. & Max abs. \\\\",
        "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(tab, file.path(in_dir, "emp_start_age.tex"))
writeLines(tab, file.path(ms_tab, "emp_start_age.tex"))
cat("Table saved to", ms_tab, "\n")
