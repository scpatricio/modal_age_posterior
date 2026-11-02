# =============================================================================
# mode_posterior.R
# Proposed method: posterior distribution of the modal age at death.
#
# D_x | m_x ~ Poisson(m_x E_x), m_x ~ Gamma(alpha, beta) independently,
# so m_x | D ~ Gamma(alpha + D_x, beta + E_x). Default alpha = 1/2,
# beta = 0 (Jeffreys). Each posterior draw of the rates is passed through
# the life table and the age with most life-table deaths is recorded
# (Algorithm 1 in the manuscript). For each draw we also record a
# continuous mode, the vertex of the parabola through d at x-1, x, x+1.
# Point estimates returned: posterior mean, median and mode of the
# continuous mode, plus the most probable integer age (map).
#
# The mode is searched over start_age..(last closed age); the open
# interval is never eligible.
#
# start_age can be a vector. The rates have independent priors, so the
# posterior of the rates above any starting age does not depend on the
# starting age, and the same draws serve every starting age (the life-table
# deaths above start_age differ only by a factor per draw, which does not
# move the mode). With several starting ages the result is a list named by
# starting age; with one it is the result itself.
# =============================================================================

mode_posterior = function(age, deaths, exposure, start_age = 40, max_age = 110,
                          S = 1e4, alpha = 0.5, beta = 0, cred = 0.95,
                          grid = 0.1, major = 0.5, return_draws = FALSE) {
  s = prep_schedule(age, deaths, exposure, min(start_age), max_age)
  n = length(s$age)

  # posterior draws of the rates: draws in rows, ages in columns
  m = matrix(rgamma(S * n, shape = rep(alpha + s$deaths, each = S),
                    rate = rep(beta + s$exposure, each = S)), nrow = S, ncol = n)

  # life-table deaths, column by column (constant force: q = 1 - exp(-m))
  d = -expm1(-m)
  d[, n] = 1
  l = rep(1, S)
  for (i in 1:n) {
    qi = d[, i]
    d[, i] = l * qi
    l = l * (1 - qi)
  }

  out = lapply(start_age, function(st) {
    posterior_summary(d, s$age, st, S, cred, grid, major, return_draws)
  })
  names(out) = as.character(start_age)
  if (length(start_age) == 1) out[[1]] else out
}

# Summaries for one starting age, from the matrix of life-table deaths
posterior_summary = function(d, ages, start_age, S, cred, grid, major, return_draws) {
  n = length(ages)
  cols = which(ages >= start_age & seq_len(n) < n)     # eligible ages
  ne = length(cols)
  ag = ages[cols]
  if (!exists("mode_stats_cpp", mode = "function")) de = d[, cols, drop = FALSE]

  if (exists("mode_stats_cpp", mode = "function")) {
    # C++ version (rw_core.cpp), same definitions as below
    ms = mode_stats_cpp(d, as.integer(ages), min(cols) - 1L, max(cols) - 1L, major)
    idx = ms$idx
    mode_int = ag[idx]
    mode_cont = ms$mode_cont
    n_local = ms$n_local
    n_major = ms$n_major
  } else {
    idx = max.col(de, ties.method = "first")
    rows = seq_len(S)
    dmax = de[cbind(rows, idx)]
    mode_int = ag[idx]

    # continuous (parabola) mode per draw; boundary draws get x + 1/2
    mode_cont = mode_int + 0.5
    inner = which(idx > 1 & idx < ne)
    if (length(inner) > 0) {
      mode_cont[inner] = parabola_mode(mode_int[inner],
                                       de[cbind(inner, idx[inner] - 1)],
                                       dmax[inner],
                                       de[cbind(inner, idx[inner] + 1)])
    }

    # interior local maxima per draw, all and "major" (>= major * max)
    mid = de[, 2:(ne - 1), drop = FALSE]
    is_max = mid > de[, 1:(ne - 2), drop = FALSE] & mid >= de[, 3:ne, drop = FALSE]
    n_local = rowSums(is_max)
    n_major = rowSums(is_max & mid >= major * dmax)
  }

  # posterior over integer ages
  probs = tabulate(idx, nbins = ne) / S
  post = data.frame(age = ag, prob = probs)

  # smallest set of ages with probability >= cred
  o = order(probs, decreasing = TRUE)
  k = which(cumsum(probs[o]) >= cred - 1e-12)[1]
  cred_set = sort(ag[o[1:k]])

  # posterior of the continuous mode on a grid (bin lower edges)
  b0 = floor(min(mode_cont) / grid + 1e-9)
  bins = floor(mode_cont / grid + 1e-9) - b0 + 1
  cnt = tabulate(bins)
  keep = cnt > 0
  post_cont = data.frame(mode = ((which(keep) - 1) + b0) * grid, prob = cnt[keep] / S)

  # posterior mode: maximum of a kernel density estimate of the draws
  kd = density(mode_cont, n = 2048)

  out = list(
    map = ag[which.max(probs)],          # most probable integer age
    p_map = max(probs),                  # its posterior probability
    mean = mean(mode_cont),              # posterior mean, continuous mode
    median = median(mode_cont),          # posterior median
    mode = kd$x[which.max(kd$y)],        # posterior mode (kernel density)
    q025 = unname(quantile(mode_cont, 0.025)),
    q975 = unname(quantile(mode_cont, 0.975)),
    cred_set = cred_set,
    cred_size = length(cred_set),
    cred_contiguous = all(diff(cred_set) == 1),
    exp_local_modes = mean(n_local),     # all interior local maxima
    exp_major_modes = mean(n_major),     # local maxima >= major * max
    p_multimodal = mean(n_major > 1),
    post = post,
    post_cont = post_cont,
    start_age = start_age,
    last_age = ages[n]
  )
  if (return_draws) {
    out$draws_int = mode_int
    out$draws_cont = mode_cont
  }
  out
}
