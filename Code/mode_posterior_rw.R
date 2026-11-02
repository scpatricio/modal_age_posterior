# =============================================================================
# mode_posterior_rw.R
# Posterior distribution of the modal age at death with a second-order
# random-walk (RW2) prior on the log death rates.
#
# Model, for the closed ages x = start_age, ..., last closed age:
#   D_x | m_x ~ Poisson(m_x E_x)
#   eta_x = log m_x
#   eta_x - 2 eta_{x-1} + eta_{x-2} ~ N(0, 1 / lambda)
# The rate of increase at age x is close to the one at x - 1, so the prior
# shrinks towards a straight line in log m_x (a Gompertz curve). It is the
# Bayesian counterpart of the second-order difference penalty of
# Whittaker-Henderson graduation and P-splines.
#
# Computation (Laplace approximation). For each lambda on a grid, the
# posterior mode of eta is found by penalized Fisher scoring and the
# posterior is approximated by N(mode, H^-1), with H = W + lambda P the
# curvature of the log posterior at the mode (W = diag(mu), P = D'D). Each
# lambda is weighted by its Laplace approximation of the marginal likelihood
# (flat prior on log10 lambda over the grid), so the uncertainty about the
# smoothing strength is carried into the draws. Draws of eta are passed
# through the life table and summarised exactly as in mode_posterior()
# (posterior_summary() in mode_posterior.R).
#
# The fit borrows strength across ages, so the posterior depends on the age
# range: with several starting ages the model is refitted for each one.
# =============================================================================

# R version of one Laplace fit (used when the C++ core cannot be compiled)
rw_fit = function(D, E, P, rank_P, lambda, eta, max_iter = 50, tol = 1e-7) {
  for (it in 1:max_iter) {
    mu = E * exp(eta)
    H = lambda * P
    diag(H) = diag(H) + mu
    R = chol(H)
    eta_new = backsolve(R, forwardsolve(t(R), mu * eta + D - mu))
    done = max(abs(eta_new - eta)) < tol
    eta = eta_new
    if (done) break
  }
  mu = E * exp(eta)
  H = lambda * P
  diag(H) = diag(H) + mu
  R = chol(H)
  loglik = sum(D * log(mu) - mu - lgamma(D + 1))
  pen = as.numeric(crossprod(eta, P %*% eta))
  logml = loglik - lambda * pen / 2 + rank_P / 2 * log(lambda) - sum(log(diag(R)))
  list(eta = eta, R = R, logml = logml)
}

# R version of the grid search and draws. The search stops once the log
# marginal likelihood is more than `drop` below the best value (past the
# peak), where the weights are negligible.
rw_draws_r = function(D, E, log10_lambdas, drop, S, n) {
  nc = length(D)
  P = crossprod(diff(diag(nc), differences = 2))
  rank_P = nc - 2
  cells = list()
  eta = log((D + 0.5) / E)
  best = -Inf
  for (ll in log10_lambdas) {
    f = rw_fit(D, E, P, rank_P, 10^ll, eta)
    eta = f$eta
    cells[[length(cells) + 1]] = list(log10_lambda = ll, eta = f$eta, R = f$R,
                                      logml = f$logml)
    best = max(best, f$logml)
    if (f$logml < best - drop) break
  }
  logml = sapply(cells, `[[`, "logml")
  w = exp(logml - max(logml))
  w = w / sum(w)
  counts = as.vector(rmultinom(1, S, w))
  cell = sample(rep(seq_along(w), counts))
  eta_draws = matrix(NA_real_, S, nc)
  for (j in unique(cell)) {
    rows = which(cell == j)
    z = matrix(rnorm(nc * length(rows)), nc, length(rows))
    eta_draws[rows, ] = t(cells[[j]]$eta + backsolve(cells[[j]]$R, z))
  }
  m = exp(eta_draws)
  d = cbind(-expm1(-m), 1)                    # constant force: q = 1 - exp(-m)
  l = rep(1, S)
  for (i in 1:n) {
    qi = d[, i]
    d[, i] = l * qi
    l = l * (1 - qi)
  }
  list(d = d, hyper = data.frame(log10_lambda = sapply(cells, `[[`, "log10_lambda"),
                                 weight = w))
}

# Compile the C++ core once when this file is sourced (needs Rcpp and a
# compiler: Xcode command line tools on macOS). Falls back to R otherwise.
# Set the environment variable MODE_RW_NO_CPP=1 to force the R version.
rw_use_cpp = FALSE
if (Sys.getenv("MODE_RW_NO_CPP") != "1" && requireNamespace("Rcpp", quietly = TRUE)) {
  rw_use_cpp = tryCatch({
    Rcpp::sourceCpp("rw_core.cpp")
    TRUE
  }, error = function(e) {
    message("Could not compile rw_core.cpp; using the R version. ", conditionMessage(e))
    FALSE
  })
}
message("mode_posterior_rw: using the ", if (rw_use_cpp) "C++" else "R", " core")

mode_posterior_rw = function(age, deaths, exposure, start_age = 40, max_age = 110,
                             S = 5000, log10_lambdas = seq(-1, 8, by = 0.25),
                             cred = 0.95, grid = 0.1, major = 0.5,
                             return_draws = FALSE) {
  # several starting ages: refit for each one
  if (length(start_age) > 1) {
    out = lapply(start_age, function(st) {
      mode_posterior_rw(age, deaths, exposure, st, max_age, S, log10_lambdas,
                        cred, grid, major, return_draws)
    })
    names(out) = as.character(start_age)
    return(out)
  }

  s = prep_schedule(age, deaths, exposure, start_age, max_age)
  n = length(s$age)
  nc = n - 1                                  # closed ages (rates modelled)
  D = s$deaths[1:nc]
  E = s$exposure[1:nc]
  drop = 12

  if (rw_use_cpp) {
    # C++ core (rw_core.cpp): banded Cholesky, grid over lambda, draws,
    # life table
    g = rw_grid_cpp(D, E, 2L, log10_lambdas, drop)
    w = exp(g$logml - max(g$logml))
    w = w / sum(w)
    counts = as.vector(rmultinom(1, S, w))
    cell = sample(rep(seq_along(w) - 1L, counts))      # random order of draws
    d = rw_draw_d_cpp(g$eta, g$L, cell, 2L)
    hyper = data.frame(log10_lambda = g$log10_lambda, weight = w)
  } else {
    r = rw_draws_r(D, E, log10_lambdas, drop, S, n)
    d = r$d
    hyper = r$hyper
  }

  out = posterior_summary(d, s$age, start_age, S, cred, grid, major, return_draws)
  out$post_log10_lambda = sum(hyper$weight * hyper$log10_lambda)
  out$hyper = hyper
  out
}
