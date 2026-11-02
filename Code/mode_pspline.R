# =============================================================================
# mode_pspline.R
# P-spline benchmark (Ouellette & Bourbeau 2011; Vazquez-Castillo et al.
# 2024). Poisson P-spline for the death counts with offset log(E_x), as in
# MortalitySmooth::Mort1Dsmooth (Camarda 2012): cubic B-splines on equally
# spaced knots, second-order difference penalty, smoothing parameter chosen
# by BIC. MortalitySmooth is archived on CRAN, so the fit is coded here.
#
# The fitted log-hazard is a smooth function of continuous age. The mode is
# the maximum of f(t) = mu(t) exp(-int_{start_age}^t mu(u) du), found on a
# fine grid of ages (0.01 years) and refined by a parabola through the
# three grid points around the maximum.
# =============================================================================

bspline_basis = function(x, xl, xr, ndx, deg = 3) {
  dx = (xr - xl) / ndx
  knots = seq(xl - deg * dx, xr + deg * dx, by = dx)
  splines::splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
}

# penalized IRLS for one lambda
pspline_irls = function(B, D, E, P, lambda, a, max_iter = 100, tol = 1e-8) {
  off = log(E)
  for (it in 1:max_iter) {
    eta = as.vector(B %*% a) + off
    mu = exp(eta)
    z = (D - mu) / mu + eta - off
    BWB = crossprod(B, B * mu)
    a_new = solve(BWB + lambda * P, crossprod(B, mu * z))
    if (max(abs(a_new - a)) < tol) {
      a = a_new
      break
    }
    a = a_new
  }
  mu = exp(as.vector(B %*% a) + off)
  BWB = crossprod(B, B * mu)
  ed = sum(diag(solve(BWB + lambda * P, BWB)))
  dev = 2 * sum(ifelse(D > 0, D * log(D / mu), 0) - (D - mu))
  list(a = as.vector(a), mu = mu, ed = ed, dev = dev,
       bic = dev + log(length(D)) * ed, iter = it)
}

mode_pspline = function(age, deaths, exposure, start_age = 40, max_age = 110,
                        ndx = NULL, deg = 3, pord = 2,
                        log10_lambdas = seq(-2, 7, by = 0.5), grid = 0.01) {
  s = prep_schedule(age, deaths, exposure, start_age, max_age)
  n = length(s$age)
  x = s$age[1:(n - 1)] + 0.5
  D = s$deaths[1:(n - 1)]
  E = s$exposure[1:(n - 1)]
  xl = s$age[1]
  xr = s$age[n]
  if (is.null(ndx)) ndx = max(4, floor((xr - xl) / 5))

  B = bspline_basis(x, xl, xr, ndx, deg)
  K = ncol(B)
  Dm = diff(diag(K), differences = pord)
  P = crossprod(Dm)

  # starting values from a ridge regression of the log rates
  y0 = log((D + 0.5) / E)
  a = as.vector(solve(crossprod(B) + 1e-4 * diag(K), crossprod(B, y0)))

  # BIC over a grid of lambdas (warm starts)
  best = NULL
  for (ll in log10_lambdas) {
    f = pspline_irls(B, D, E, P, 10^ll, a)
    a = f$a
    if (is.null(best) || f$bic < best$bic) {
      best = f
      best$log10_lambda = ll
    }
  }

  # density of age at death on a fine grid, conditional on survival to xl
  tg = seq(xl, xr, by = grid)
  Bg = bspline_basis(tg, xl, xr, ndx, deg)
  mug = exp(as.vector(Bg %*% best$a))
  H = c(0, cumsum((mug[-1] + mug[-length(mug)]) / 2 * grid))
  fg = mug * exp(-H)
  j = which.max(fg)

  # refine with the vertex of the parabola through the three grid points
  mode = tg[j]
  if (j > 1 && j < length(tg)) {
    A = fg[j] - fg[j - 1]
    B = fg[j] - fg[j + 1]
    mode = tg[j] + grid * (A / (A + B) - 0.5)
  }

  list(mode = mode, boundary = j == 1 || j == length(tg),
       log10_lambda = best$log10_lambda, ed = best$ed,
       start_age = start_age)
}
