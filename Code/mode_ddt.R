# =============================================================================
# mode_ddt.R
# Discretized derivative tests (DDT), Vazquez-Castillo, Bergeron-Boucher &
# Missov (2024, Demographic Research 50(11), Sec. 3.3, Eq. 9).
#
# With d(x) the life-table deaths and M(x) = D_x / E_x the observed rates:
#   (9.1) d(x+1) - d(x-1)                    = 0   (first derivative)
#   (9.2) d(x+1) - 2 d(x) + d(x-1)           < 0   (second derivative)
#   (9.3) ln M(x) - ln M(x-1) - M(x)         = 0   (LAR = mu at the mode)
#   (9.4) M(x+1) - M(x-1) - 2 M(x)^2         = 0   (mu' = mu^2 at the mode)
# For 9.1, 9.3 and 9.4 we find the integer ages where the left-hand side
# changes sign from positive to negative and interpolate linearly to the
# zero. Each zero of 9.1 that passes 9.2 is a candidate local maximum; the
# zeros of 9.3 and 9.4 closest to it (within `window` years) are averaged
# with it. Among candidates, the mode is the one with the largest d(x).
#
# Age convention. The published method treats d(x) as the value of the
# death density at exact age x. Life-table deaths d(x) are deaths between
# exact ages x and x+1, so they are centred at x + 1/2, and the zeros found
# above sit half a year below the mode of the density. We therefore add
# `shift` = 0.5 to the estimate (set shift = 0 for the published version).
# The other estimators use the same convention: the Gompertz and P-spline
# fits evaluate the hazard at x + 1/2, and the posterior uses Kannisto's
# formula for grouped data.
#
# Implementation choice (not fixed by the paper): how zeros of 9.3 and 9.4
# are matched to a candidate when the conditions cross zero more than once.
# Here: nearest zero within `window` years; if none, that condition is
# dropped for that candidate and the number of conditions used is reported.
# =============================================================================

zero_crossings = function(x, f) {
  n = length(f)
  i = which(is.finite(f[-n]) & is.finite(f[-1]) & f[-n] > 0 & f[-1] <= 0)
  x[i] + f[i] / (f[i] - f[i + 1])
}

mode_ddt = function(age, deaths, exposure, start_age = 40, max_age = 110,
                    window = 2, shift = 0.5) {
  s = prep_schedule(age, deaths, exposure, start_age, max_age)
  n = length(s$age)
  M = s$deaths / s$exposure
  d = as.vector(lifetable_d(M))

  # closed ages only; conditions defined for interior ages
  ne = n - 1
  i = 2:(ne - 1)
  xs = s$age[i]
  f1 = d[i + 1] - d[i - 1]
  f2 = d[i + 1] - 2 * d[i] + d[i - 1]
  lM = suppressWarnings(log(M))
  f3 = lM[i] - lM[i - 1] - M[i]
  f4 = M[i + 1] - M[i - 1] - 2 * M[i]^2

  r1 = zero_crossings(xs, f1)
  r3 = zero_crossings(xs, f3)
  r4 = zero_crossings(xs, f4)

  empty = list(mode = NA_real_, n_candidates = 0L, n_conditions = NA_integer_,
               r1 = NA_real_, r3 = NA_real_, r4 = NA_real_,
               start_age = start_age)
  if (length(r1) == 0) return(empty)

  cand = lapply(r1, function(r) {
    # integer age with more deaths among the two around the zero
    ages_r = c(floor(r), ceiling(r))
    ages_r = ages_r[ages_r %in% xs]
    if (length(ages_r) == 0) return(NULL)
    a_r = ages_r[which.max(d[match(ages_r, s$age)])]
    if (!(f2[match(a_r, xs)] < 0)) return(NULL)
    near = function(z) {
      if (length(z) == 0) return(NA_real_)
      dz = abs(z - r)
      if (min(dz) > window) return(NA_real_)
      z[which.min(dz)]
    }
    z3 = near(r3)
    z4 = near(r4)
    est = c(r, z3, z4)
    data.frame(r1 = r, r3 = z3, r4 = z4, mode = mean(est, na.rm = TRUE),
               n_conditions = sum(is.finite(est)),
               height = d[match(a_r, s$age)])
  })
  cand = do.call(rbind, cand)
  if (is.null(cand) || nrow(cand) == 0) return(empty)

  best = cand[which.max(cand$height), ]
  list(mode = best$mode + shift, n_candidates = nrow(cand),
       n_conditions = best$n_conditions,
       r1 = best$r1, r3 = best$r3, r4 = best$r4, start_age = start_age)
}
