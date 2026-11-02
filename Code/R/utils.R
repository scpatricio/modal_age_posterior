# =============================================================================
# utils.R
# Shared helpers: HMD input, schedule preparation, life table, parallel setup.
# Sourced by every method file and by the main scripts.
# =============================================================================

# ---- HMD input --------------------------------------------------------------

# Read an HMD 1x1 file (Deaths_1x1 or Exposures_1x1). The open interval
# "110+" is coded as age 110.
read_hmd = function(file) {
  x = read.table(file, skip = 2, header = TRUE, na.strings = ".",
                 stringsAsFactors = FALSE)
  x$Age = as.integer(sub("+", "", x$Age, fixed = TRUE))
  x
}

# Deaths and exposures for one country, in long format
# (country, sex, year, age, deaths, exposure).
load_country = function(country, data_dir) {
  D = read_hmd(file.path(data_dir, paste0(country, "_Deaths_1x1.txt")))
  E = read_hmd(file.path(data_dir, paste0(country, "_Exposures_1x1.txt")))
  out = list()
  for (sex in c("Female", "Male")) {
    out[[sex]] = data.frame(country = country, sex = sex,
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
