# =============================================================================
# 05_figures_empirical.R
# Figures and tables of the empirical application (Section 4, Table 1 and
# Appendix F). Women in Denmark, France, Japan and Sweden. Run 04_empirical.R
# first. Figures are saved to figures/empirical/, tables to
# results/empirical/, and both to the manuscript (see utils.R).
#
#   emp_trends.pdf             modal age at death by year and country,
#                              starting age 40: RW2 post. mean with 95%
#                              interval, Gompertz, P-spline, DDT
#   emp_start_age.pdf          Denmark: estimates by year for each method
#                              (one panel per method), one line per starting
#                              age (40, 50, 60, 70)
#   emp_start_age_all.pdf      the same for the four countries (Appendix F)
#   emp_start_age.tex          Table 1: change in each method's estimate when
#                              the starting age moves from 40 to 50, 60, 70
#   start_age_consistency.csv  mean, mean absolute and largest absolute
#                              change by method, comparison and country
# =============================================================================

if (!file.exists("utils.R")) stop("Run this script with the Code folder as working directory (open Code.Rproj).")
source("utils.R")
library(ggplot2)

in_dir = file.path("results", "empirical")
save_both = function(p, name, height) save_figure(p, name, "empirical", height)

pe = readRDS(file.path(in_dir, "point_estimates.rds"))
base_start = min(pe$start_age)           # reference starting age (40)
show_country = "Denmark"                 # country of the main starting-age figure

method_labels = c(rw2_mean = "RW2 post. mean", rw2_median = "RW2 post. median",
                  rw2_mode = "RW2 post. mode", gompertz = "Gompertz",
                  pspline = "P-spline", ddt = "DDT")
cols = setNames(c("#0072B2", "#56B4E9", "#CC79A7", "#E69F00", "#009E73", "#D55E00"),
                method_labels)
pe$method_lab = factor(method_labels[pe$method], levels = method_labels)
start_levels = paste("Starting age", sort(unique(pe$start_age)))
pe$start_lab = factor(paste("Starting age", pe$start_age), levels = start_levels)
start_cols = setNames(c("#0072B2", "#E69F00", "#009E73", "#D55E00")[seq_along(start_levels)],
                      start_levels)

theme_set(theme_mode())

# ---- Trends at the reference starting age -----------------------------------

tm = c("rw2_mean", "gompertz", "pspline", "ddt")
d = pe[pe$start_age == base_start & pe$method %in% tm, ]
band = d[d$method == "rw2_mean", ]
p = ggplot(d, aes(year, estimate, colour = method_lab)) +
  geom_ribbon(data = band, aes(ymin = lower, ymax = upper), fill = "#0072B2",
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ country, nrow = 1) +
  scale_colour_manual(values = cols[method_labels[tm]], name = NULL) +
  labs(x = "Year", y = "Modal age at death")
save_both(p, "emp_trends.pdf", 4.5)

# ---- Estimates at the four starting ages ------------------------------------

start_plot = function(d, facet) {
  ggplot(d, aes(year, estimate, colour = start_lab)) +
    geom_line(linewidth = 0.35) +
    facet +
    scale_colour_manual(values = start_cols, name = NULL) +
    labs(x = "Year", y = "Modal age at death")
}

d = pe[pe$country == show_country, ]
p = start_plot(d, facet_wrap(~ method_lab, nrow = 2))
save_both(p, "emp_start_age.pdf", 7)

p = start_plot(pe, facet_grid(method_lab ~ country)) +
  theme(strip.text.y = element_text(angle = -90, hjust = 0.5))
save_both(p, "emp_start_age_all.pdf", 12)
cat("Figures saved to figures/empirical\n")

# ---- Changes relative to the reference starting age -------------------------

key = paste(pe$country, pe$sex, pe$year, pe$method)
ref = pe[pe$start_age == base_start, ]
pe$ref = ref$estimate[match(key, paste(ref$country, ref$sex, ref$year, ref$method))]
dd = pe[pe$start_age != base_start, ]
dd$change = dd$estimate - dd$ref
dd$comparison = factor(sprintf("%d vs %d", dd$start_age, base_start),
                       levels = sprintf("%d vs %d", sort(unique(dd$start_age)), base_start))

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

# ---- LaTeX table (Table 1 of the manuscript) --------------------------------

f2 = function(x, digits = 4) {
  out = formatC(x, format = "f", digits = digits)
  out = ifelse(abs(x) < 0.5 * 10^-digits, formatC(0, format = "f", digits = digits), out)
  ifelse(substr(out, 1, 1) == "-", paste0("$-$", substring(out, 2)), out)
}
row_methods = c("rw2_mean", "rw2_median", "rw2_mode", "pspline", "gompertz", "ddt")
comps = levels(dd$comparison)
rows = sapply(row_methods, function(m) {
  v = unlist(lapply(comps, function(cp) {
    z = tab_all[tab_all$method == m & tab_all$comparison == cp, ]
    c(z$mean_abs_change, z$max_abs_change)
  }))
  paste0(method_labels[m], " & ", paste(f2(v), collapse = " & "), " \\\\")
})
n_pop = length(unique(paste(pe$country, pe$sex, pe$year)))
later = sort(unique(dd$start_age))
tab = c("\\begin{table}[htb!]",
        "\\centering",
        sprintf(paste0("\\caption{Change in the estimated modal age at death when the starting age moves ",
                       "from %d to %s: mean absolute change and largest absolute change (years), ",
                       "women in %d countries, all years (%d population-years). ",
                       "DDT uses the rates at each age separately and does not depend on the starting age.}"),
                base_start, paste0(paste(head(later, -1), collapse = ", "), " and ", tail(later, 1)),
                length(unique(pe$country)), n_pop),
        "\\label{tab:start_age}",
        "\\small",
        sprintf("\\begin{tabular}{l%s}", strrep("rr", length(comps))),
        "\\hline",
        paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{%s}", comps), collapse = " & "), " \\\\"),
        paste0("Method & ", paste(rep("Mean abs. & Max abs.", length(comps)), collapse = " & "),
               " \\\\"),
        "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{table}")
save_table(tab, "emp_start_age.tex", "empirical")
cat("Table saved to results/empirical\n")
