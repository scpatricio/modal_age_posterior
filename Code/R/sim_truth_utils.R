# =============================================================================
# sim_truth_utils.R
# True distributions and expected data for the simulation experiments.
# Sourced by 01_simulation.R and 05_simulation_change.R.
# =============================================================================

max_age = 110

# mixture components (both Gompertz with the same slope)
M1 = 90
b1 = 0.14
M2 = 65
b2 = 0.14

# Gompertz with modal age M and slope b: a = b exp(-b M)
gomp_a = function(M, b) b * exp(-b * M)
gomp_S = function(x, M, b) exp(-(gomp_a(M, b) / b) * (exp(b * x) - 1))
gomp_f = function(x, M, b) gomp_a(M, b) * exp(b * x) * gomp_S(x, M, b)

mix_f = function(x, w) w * gomp_f(x, M1, b1) + (1 - w) * gomp_f(x, M2, b2)
mix_S = function(x, w) w * gomp_S(x, M1, b1) + (1 - w) * gomp_S(x, M2, b2)

# local maxima and minima of the mixture density on [20, 110]
find_peaks = function(w) {
  x = seq(20, 110, by = 0.01)
  fx = mix_f(x, w)
  s = diff(sign(diff(fx)))
  refine = function(i, max) {
    optimize(function(t) mix_f(t, w), c(x[i - 1], x[i + 1]), maximum = max)
  }
  pk = which(s == -2) + 1
  vl = which(s == 2) + 1
  peaks = data.frame(x = sapply(pk, function(i) refine(i, TRUE)$maximum))
  peaks$f = mix_f(peaks$x, w)
  valleys = unlist(lapply(vl, function(i) refine(i, FALSE)$minimum))
  list(peaks = peaks, valleys = valleys)
}

# Expected exposures and true central rates for ages 0..max_age (last age
# open), for a survival function Sfun and radix N.
expected_schedule_S = function(Sfun, N) {
  age = 0:max_age
  L = sapply(age, function(x) {
    upper = if (x == max_age) x + 30 else x + 1
    integrate(Sfun, x, upper)$value
  })
  dx = Sfun(age) - c(Sfun(age[-1]), 0)
  data.frame(age = age, exposure = N * L, m = dx / L)
}

expected_schedule = function(w, N) {
  expected_schedule_S(function(t) mix_S(t, w), N)
}
