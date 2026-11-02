// =============================================================================
// rw_core.cpp
// C++ core of mode_posterior_rw(): grid of Laplace fits over the smoothing
// parameter lambda and posterior draws. The posterior precision
// H = W + lambda P is banded (bandwidth 2 for RW2), so every Cholesky
// factorisation and triangular solve is done in band storage: O(n) work per
// fit instead of O(n^3).
//
// Band storage: element (i, k) of a lower-banded matrix with i - k <= p is
// stored at position i * (p + 1) + (i - k).
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

static inline int bidx(int i, int k, int p) { return i * (p + 1) + (i - k); }

// Cholesky H = L L' in band storage; false if H is not positive definite
static bool chol_band(const std::vector<double>& A, std::vector<double>& L, int n, int p) {
  for (int i = 0; i < n; i++) {
    int j0 = std::max(0, i - p);
    for (int j = j0; j <= i; j++) {
      double s = A[bidx(i, j, p)];
      int k0 = std::max(j0, j - p);
      for (int k = k0; k < j; k++) s -= L[bidx(i, k, p)] * L[bidx(j, k, p)];
      if (i == j) {
        if (s <= 0) return false;
        L[bidx(i, i, p)] = std::sqrt(s);
      } else {
        L[bidx(i, j, p)] = s / L[bidx(j, j, p)];
      }
    }
  }
  return true;
}

// solve L y = b, then L' x = y (x = H^-1 b)
static void solve_band(const std::vector<double>& L, const std::vector<double>& b,
                       std::vector<double>& x, int n, int p) {
  std::vector<double> y(n);
  for (int i = 0; i < n; i++) {
    double s = b[i];
    for (int k = std::max(0, i - p); k < i; k++) s -= L[bidx(i, k, p)] * y[k];
    y[i] = s / L[bidx(i, i, p)];
  }
  for (int i = n - 1; i >= 0; i--) {
    double s = y[i];
    for (int k = i + 1; k <= std::min(n - 1, i + p); k++) s -= L[bidx(k, i, p)] * x[k];
    x[i] = s / L[bidx(i, i, p)];
  }
}

// penalty matrix P = D'D (D: differences of order p) in band storage
static std::vector<double> penalty_band(int n, int p) {
  std::vector<double> c(p + 1);
  for (int j = 0; j <= p; j++) {
    double v = 1;
    for (int t = 0; t < j; t++) v = v * (p - t) / (t + 1);
    c[j] = ((p - j) % 2 == 0) ? v : -v;     // coefficients of the p-th difference
  }
  std::vector<double> P(n * (p + 1), 0.0);
  for (int r = 0; r + p < n; r++) {         // row r of D touches columns r..r+p
    for (int a = 0; a <= p; a++) {
      for (int b = 0; b <= a; b++) {
        int i = r + a;
        int k = r + b;
        P[bidx(i, k, p)] += c[a] * c[b];
      }
    }
  }
  return P;
}

static double penalty_value(const std::vector<double>& eta, int n, int p) {
  // ||D eta||^2 by repeated differencing
  std::vector<double> d(eta);
  int m = n;
  for (int o = 0; o < p; o++) {
    for (int i = 0; i < m - 1; i++) d[i] = d[i + 1] - d[i];
    m--;
  }
  double s = 0;
  for (int i = 0; i < m; i++) s += d[i] * d[i];
  return s;
}

// [[Rcpp::export]]
List rw_grid_cpp(NumericVector D, NumericVector E, int order,
                 NumericVector log10_lambdas, double drop,
                 int max_iter = 50, double tol = 1e-7) {
  int n = D.size();
  int p = order;
  int rank_P = n - p;
  std::vector<double> P = penalty_band(n, p);
  std::vector<double> eta(n);
  for (int i = 0; i < n; i++) eta[i] = std::log((D[i] + 0.5) / E[i]);

  std::vector<double> out_ll, out_logml;
  std::vector<std::vector<double> > out_eta, out_L;
  std::vector<double> H(n * (p + 1)), L(n * (p + 1)), mu(n), rhs(n), eta_new(n);

  double best = R_NegInf;
  for (int li = 0; li < log10_lambdas.size(); li++) {
    double lambda = std::pow(10.0, log10_lambdas[li]);
    bool ok = true;
    // penalized Fisher scoring (Poisson, log link), warm start from the
    // previous lambda
    for (int it = 0; it < max_iter; it++) {
      for (int i = 0; i < n; i++) {
        mu[i] = E[i] * std::exp(eta[i]);
        rhs[i] = mu[i] * eta[i] + D[i] - mu[i];
      }
      for (int t = 0; t < n * (p + 1); t++) H[t] = lambda * P[t];
      for (int i = 0; i < n; i++) H[bidx(i, i, p)] += mu[i];
      if (!chol_band(H, L, n, p)) { ok = false; break; }
      solve_band(L, rhs, eta_new, n, p);
      double dmax = 0;
      for (int i = 0; i < n; i++) dmax = std::max(dmax, std::fabs(eta_new[i] - eta[i]));
      eta = eta_new;
      if (dmax < tol) break;
    }
    if (!ok) continue;
    // curvature H = W + lambda P and Laplace log marginal likelihood
    double loglik = 0;
    for (int i = 0; i < n; i++) {
      mu[i] = E[i] * std::exp(eta[i]);
      loglik += D[i] * std::log(mu[i]) - mu[i] - R::lgammafn(D[i] + 1);
    }
    for (int t = 0; t < n * (p + 1); t++) H[t] = lambda * P[t];
    for (int i = 0; i < n; i++) H[bidx(i, i, p)] += mu[i];
    if (!chol_band(H, L, n, p)) continue;
    double logdet_half = 0;
    for (int i = 0; i < n; i++) logdet_half += std::log(L[bidx(i, i, p)]);
    double logml = loglik - lambda * penalty_value(eta, n, p) / 2 +
      rank_P / 2.0 * std::log(lambda) - logdet_half;

    out_ll.push_back(log10_lambdas[li]);
    out_logml.push_back(logml);
    out_eta.push_back(eta);
    out_L.push_back(L);
    best = std::max(best, logml);
    if (logml < best - drop) break;          // past the peak, negligible weight
  }

  int nc = out_ll.size();
  NumericMatrix eta_m(n, nc), L_m(n * (p + 1), nc);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < n; i++) eta_m(i, j) = out_eta[j][i];
    for (int t = 0; t < n * (p + 1); t++) L_m(t, j) = out_L[j][t];
  }
  return List::create(_["log10_lambda"] = wrap(out_ll), _["logml"] = wrap(out_logml),
                      _["eta"] = eta_m, _["L"] = L_m);
}

// Draws of eta ~ N(eta_hat, H^-1) for the grid cell of each draw, passed
// through the life table. Returns life-table deaths: draws in rows, ages in
// columns; the last column is the open interval (q = 1).
// [[Rcpp::export]]
NumericMatrix rw_draw_d_cpp(NumericMatrix eta, NumericMatrix L, IntegerVector cell,
                            int order) {
  int n = eta.nrow();
  int p = order;
  int S = cell.size();
  NumericMatrix d(S, n + 1);
  std::vector<double> z(n), x(n), Lj(n * (p + 1));
  int last = -1;
  for (int s = 0; s < S; s++) {
    int j = cell[s];
    if (j != last) {
      for (int t = 0; t < n * (p + 1); t++) Lj[t] = L(t, j);
      last = j;
    }
    for (int i = 0; i < n; i++) z[i] = R::norm_rand();
    // x = L'^-1 z, so that cov(x) = H^-1
    for (int i = n - 1; i >= 0; i--) {
      double v = z[i];
      for (int k = i + 1; k <= std::min(n - 1, i + p); k++) v -= Lj[bidx(k, i, p)] * x[k];
      x[i] = v / Lj[bidx(i, i, p)];
    }
    double l = 1;
    for (int i = 0; i < n; i++) {
      double m = std::exp(eta(i, j) + x[i]);
      double q = -std::expm1(-m);           // constant force: q = 1 - exp(-m)
      d(s, i) = l * q;
      l *= 1 - q;
    }
    d(s, n) = l;
  }
  return d;
}

// Per-draw summaries of the life-table deaths over the eligible columns
// c0..c1 (0-based): index of the largest value (first in case of ties),
// parabola-vertex mode, number of interior local maxima and of those with
// height >= major * max. Ages are ages[c0..c1]. Same definitions as
// posterior_summary() in mode_posterior.R.
// [[Rcpp::export]]
List mode_stats_cpp(NumericMatrix d, IntegerVector ages, int c0, int c1, double major) {
  int S = d.nrow();
  IntegerVector idx(S), n_local(S), n_major(S);
  NumericVector mode_cont(S);
  int ne = c1 - c0 + 1;
  for (int s = 0; s < S; s++) {
    int best = c0;
    double dmax = d(s, c0);
    for (int j = c0 + 1; j <= c1; j++) {
      if (d(s, j) > dmax) { dmax = d(s, j); best = j; }
    }
    int b = best - c0;                     // position within eligible ages
    idx[s] = b + 1;                        // 1-based, as max.col
    double x = ages[best];
    if (b > 0 && b < ne - 1) {
      double A = dmax - d(s, best - 1);
      double B = dmax - d(s, best + 1);
      mode_cont[s] = x + A / (A + B);
    } else {
      mode_cont[s] = x + 0.5;
    }
    int nl = 0, nm = 0;
    for (int j = c0 + 1; j < c1; j++) {
      if (d(s, j) > d(s, j - 1) && d(s, j) >= d(s, j + 1)) {
        nl++;
        if (d(s, j) >= major * dmax) nm++;
      }
    }
    n_local[s] = nl;
    n_major[s] = nm;
  }
  return List::create(_["idx"] = idx, _["mode_cont"] = mode_cont,
                      _["n_local"] = n_local, _["n_major"] = n_major);
}

// Exact posterior sampling for checking the Laplace approximation:
// independence Metropolis-Hastings with the Laplace mixture as proposal.
// Target: p(j, eta | D) proportional to (1/J) * lambda_j^((K-2)/2)
//   * exp(-lambda_j/2 eta'P eta) * exp(sum D eta - E exp(eta)),
// i.e. the same model as mode_posterior_rw() (flat prior on the grid of
// log10 lambda). Proposal: q(j, eta) = w_j N(eta; eta_hat_j, H_j^-1).
// The accept/reject step corrects any error of the approximation, so the
// chain targets the exact posterior; its acceptance rate measures how close
// the approximation is (1 when it is exact).
// [[Rcpp::export]]
List rw_imh_cpp(NumericVector D, NumericVector E, NumericMatrix eta_hat, NumericMatrix L,
                NumericVector log10_lambda, NumericVector w, int n_iter, int order) {
  int n = eta_hat.nrow();
  int J = eta_hat.ncol();
  int p = order;
  std::vector<double> logdet_half(J, 0.0), cumw(J);
  double cs = 0;
  for (int j = 0; j < J; j++) {
    for (int i = 0; i < n; i++) logdet_half[j] += std::log(L(bidx(i, i, p), j));
    cs += w[j];
    cumw[j] = cs;
  }
  // log target and log proposal density (up to constants common to all states)
  auto log_target = [&](int j, const std::vector<double>& eta) {
    double lambda = std::pow(10.0, log10_lambda[j]);
    double ll = 0;
    for (int i = 0; i < n; i++) ll += D[i] * eta[i] - E[i] * std::exp(eta[i]);
    return (n - p) / 2.0 * std::log(lambda) - lambda / 2.0 * penalty_value(eta, n, p) + ll;
  };
  std::vector<double> z(n), x(n);
  auto propose = [&](int& j, std::vector<double>& eta, double& logq) {
    double u = R::unif_rand() * cs;
    j = 0;
    while (j < J - 1 && cumw[j] < u) j++;
    double zsq = 0;
    for (int i = 0; i < n; i++) { z[i] = R::norm_rand(); zsq += z[i] * z[i]; }
    for (int i = n - 1; i >= 0; i--) {          // x = L'^-1 z
      double v = z[i];
      for (int k = i + 1; k <= std::min(n - 1, i + p); k++) v -= L(bidx(k, i, p), j) * x[k];
      x[i] = v / L(bidx(i, i, p), j);
    }
    for (int i = 0; i < n; i++) eta[i] = eta_hat(i, j) + x[i];
    logq = std::log(w[j]) + logdet_half[j] - 0.5 * zsq;
  };

  NumericMatrix draws(n_iter, n);
  IntegerVector cell(n_iter);
  std::vector<double> eta_cur(n), eta_new(n);
  int j_cur, j_new;
  double logq_cur, logq_new;
  propose(j_cur, eta_cur, logq_cur);
  double logp_cur = log_target(j_cur, eta_cur);
  int accepted = 0;
  for (int t = 0; t < n_iter; t++) {
    propose(j_new, eta_new, logq_new);
    double logp_new = log_target(j_new, eta_new);
    double log_alpha = (logp_new - logq_new) - (logp_cur - logq_cur);
    if (std::log(R::unif_rand()) < log_alpha) {
      eta_cur = eta_new;
      j_cur = j_new;
      logq_cur = logq_new;
      logp_cur = logp_new;
      accepted++;
    }
    for (int i = 0; i < n; i++) draws(t, i) = eta_cur[i];
    cell[t] = j_cur;
  }
  return List::create(_["eta"] = draws, _["cell"] = cell,
                      _["acceptance"] = (double) accepted / n_iter);
}
