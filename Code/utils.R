# =============================================================================
# utils.R
# Shared helpers: HMD input, schedule preparation, life table, parallel
# computation and saving of outputs. Sourced by every script.
#
# All scripts are run with the Code folder as working directory (open
# Code.Rproj). Outputs go to figures/ and results/. Figures and tables that
# appear in the manuscript are also copied to ../Manuscript/figures and
# ../Manuscript/tables, from where main.tex includes them (skipped if the
# Manuscript folder is not there).
# =============================================================================

# ---- Outputs ----------------------------------------------------------------

manuscript_dir = Sys.getenv("MODE_MANUSCRIPT_DIR", file.path("..", "Manuscript"))

# ---- Figure settings: change these two numbers to restyle every figure ------

fig_width = 11       # width of every figure, in inches
fig_text_size = 14   # base font size of every figure, in points

# common theme of all figures
theme_mode = function(base_size = fig_text_size) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "bottom")
}

# save a ggplot to figures/<sub>/ (and to the manuscript, if present); every
# figure has width fig_width
save_figure = function(p, name, sub, height, width = fig_width) {
  d = file.path("figures", sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(d, name), p, width = width, height = height)
  if (dir.exists(manuscript_dir)) {
    md = file.path(manuscript_dir, "figures")
    dir.create(md, showWarnings = FALSE)
    ggplot2::ggsave(file.path(md, name), p, width = width, height = height)
  }
}

# save the lines of a LaTeX table to results/<sub>/ (and to the manuscript)
save_table = function(lines, name, sub) {
  d = file.path("results", sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(d, name))
  if (dir.exists(manuscript_dir)) {
    md = file.path(manuscript_dir, "tables")
    dir.create(md, showWarnings = FALSE)
    writeLines(lines, file.path(md, name))
  }
}

# ---- HMD input --------------------------------------------------------------

# Read an HMD 1x1 file (Deaths_1x1 or Exposures_1x1). The open interval
# "110+" is coded as age 110.
read_hmd = function(file) {
  x = read.table(file, skip = 2, header = TRUE, na.strings = ".",
                 stringsAsFactors = FALSE)
  x$Age = as.integer(sub("+", "", x$Age, fixed = TRUE))
  x
}

# HMD country codes (folder names in data/) and the names used in the paper
hmd_countries = c(DNK = "Denmark", FRATNP = "France", JPN = "Japan", SWE = "Sweden")

# Deaths and exposures for one country, read from the HMD country folder
# data/<code>/STATS/ (Deaths_1x1.txt, Exposures_1x1.txt), in long format
# (country, sex, year, age, deaths, exposure). Women only by default.
load_country = function(code, data_dir = "data", sexes = "Female") {
  D = read_hmd(file.path(data_dir, code, "STATS", "Deaths_1x1.txt"))
  E = read_hmd(file.path(data_dir, code, "STATS", "Exposures_1x1.txt"))
  E = E[match(paste(D$Year, D$Age), paste(E$Year, E$Age)), ]
  name = if (code %in% names(hmd_countries)) hmd_countries[[code]] else code
  out = list()
  for (sex in sexes) {
    out[[sex]] = data.frame(country = name, sex = sex,
                            year = D$Year, age = D$Age,
                            deaths = D[[sex]], exposure = E[[sex]])
  }
  do.call(rbind, out)
}

# ---- Schedule preparation ---------------------------------------------------

# Keep ages start_age..max_age. The last age kept is treated as an open
# interval (everyone alive dies there). If an age has zero or missing
# exposure, the schedule is truncated just before it, so the previous age
# becomes the open interval.
prep_schedule = function(age, deaths, exposure, start_age, max_age = 110) {
  o = order(age)
  age = age[o]
  deaths = deaths[o]
  exposure = exposure[o]
  keep = age >= start_age & age <= max_age
  age = age[keep]
  deaths = deaths[keep]
  exposure = exposure[keep]
  bad = which(!is.finite(exposure) | exposure <= 0 | !is.finite(deaths))
  if (length(bad) > 0) {
    last = bad[1] - 1
    if (last < 3) stop("Too few ages with positive exposure above start_age.")
    age = age[1:last]
    deaths = deaths[1:last]
    exposure = exposure[1:last]
  }
  list(age = age, deaths = deaths, exposure = exposure)
}

# ---- Life table -------------------------------------------------------------

# Life-table deaths from death rates on single-year ages, assuming a constant
# force of mortality within each interval (q_x = 1 - exp(-m_x), Eq. 1 of the
# manuscript), l at the first age = 1, and the last age as open interval.
# m can be a vector (one schedule) or a matrix (ages in rows, draws in
# columns). Returns d with the same shape.
lifetable_d = function(m) {
  m = as.matrix(m)
  n = nrow(m)
  q = -expm1(-m)
  q[n, ] = 1
  l = matrix(1, nrow = n, ncol = ncol(m))
  for (i in 2:n) {
    l[i, ] = l[i - 1, ] * (1 - q[i - 1, ])
  }
  l * q
}

# Continuous mode from single-year life-table deaths: vertex of the parabola
# through the deaths at ages x-1, x, x+1, where x is the age with most deaths
# and d_x refers to the interval [x, x+1). Algebraically this is
# Kannisto's (2001) formula.
parabola_mode = function(x, d_prev, d_x, d_next) {
  A = d_x - d_prev
  B = d_x - d_next
  x + A / (A + B)
}

# ---- Parallel setup ---------------------------------------------------------

# Number of cores to use. Forked parallelism (mclapply) works on macOS and
# Linux. On Windows this falls back to one core.
get_cores = function(reserve = 1) {
  if (.Platform$OS.type == "windows") return(1)
  max(1, parallel::detectCores() - reserve)
}

# Parallel lapply with a progress bar. Tasks run in chunks of a few tasks per
# core (forked with mclapply); after each chunk the bar, elapsed time and an
# estimate of the remaining time are updated. Failed tasks come back as
# "try-error" objects, as with mclapply.
par_lapply = function(X, FUN, n_cores, label = "Progress") {
  n = length(X)
  chunk = n_cores * max(1, ceiling(n / (100 * n_cores)))
  starts = seq(1, n, by = chunk)
  out = vector("list", n)
  t0 = Sys.time()
  fmt = function(sec) {
    sec = round(sec)
    sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
  }
  for (s in starts) {
    i = s:min(n, s + chunk - 1)
    out[i] = parallel::mclapply(X[i], FUN, mc.cores = n_cores, mc.preschedule = TRUE)
    done = max(i)
    el = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    width = 30
    filled = round(width * done / n)
    cat(sprintf("\r%s [%s%s] %3.0f%% (%d/%d) elapsed %s, remaining ~%s   ",
                label, strrep("=", filled), strrep(" ", width - filled),
                100 * done / n, done, n, fmt(el), fmt(el / done * (n - done))))
    flush.console()
  }
  cat("\n")
  names(out) = names(X)
  out
}
