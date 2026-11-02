# =============================================================================
# mode_gompertz.R
# Gompertz benchmark: mu(x) = a exp(b x), fitted by Poisson maximum
# likelihood to ages start_age..(last closed age), with the rate in [x, x+1)
# taken as mu(x + 1/2). This is a log-linear Poisson GLM with offset
# log(E_x). Mode = ln(b / a) / b (Canudas-Romo 2008; Missov et al. 2015).
# =============================================================================

mode_gompertz = function(age, deaths, exposure, start_age = 40, max_age = 110) {
  s = prep_schedule(age, deaths, exposure, start_age, max_age)
  n = length(s$age)
  x = s$age[1:(n - 1)] + 0.5
  D = s$deaths[1:(n - 1)]
  E = s$exposure[1:(n - 1)]

  X = cbind(1, x)
  # quasipoisson gives the Poisson ML estimates and accepts the
  # non-integer death counts found in HMD
  fit = glm.fit(X, D, offset = log(E), family = quasipoisson())
  a = exp(fit$coefficients[1])
  b = fit$coefficients[2]
  mode = if (b > a) log(b / a) / b else NA_real_

  list(mode = unname(mode), a = unname(a), b = unname(b),
       converged = fit$converged, start_age = start_age)
}
