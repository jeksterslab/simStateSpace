// -----------------------------------------------------------------------------
// edit .setup/cpp/001-transformation-helpers.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat AccumulatedCovariance(const arma::mat& a, const arma::mat& q,
                                const double delta_t) {
  const arma::uword n = a.n_rows;
  // Van Loan's block contains both A and -A. Even a stable drift can
  // overflow on a long interval, so integrate a short step and compose
  // covariances using Q(2h) = Q(h) + exp(Ah) Q(h) exp(Ah)'.
  const double scaled_norm = arma::norm(a, "inf") * delta_t;
  if (!std::isfinite(scaled_norm)) {
    // # nocov start
    Rcpp::stop("Scaled drift must be finite for covariance integration.");
    // # nocov end
  }
  const int steps = scaled_norm > 1.0
                        ? static_cast<int>(std::ceil(std::log2(scaled_norm)))
                        : 0;
  const double h = std::ldexp(delta_t, -steps);
  arma::mat block = arma::zeros(2 * n, 2 * n);

  block.submat(0, 0, n - 1, n - 1) = -a;
  block.submat(0, n, n - 1, 2 * n - 1) = q;
  block.submat(n, n, 2 * n - 1, 2 * n - 1) = a.t();

  arma::mat exp_block;
  if (!arma::expmat(exp_block, block * h)) {
    // # nocov start
    Rcpp::stop("Matrix exponential for accumulated covariance failed.");
    // # nocov end
  }

  arma::mat exp_a;
  if (!arma::expmat(exp_a, a * h)) {
    Rcpp::stop("Matrix exponential for companion drift failed.");  // # nocov
  }

  arma::mat out = exp_a * exp_block.submat(0, n, n - 1, 2 * n - 1);
  out = 0.5 * (out + out.t());
  for (int j = 0; j < steps; ++j) {
    out = out + exp_a * out * exp_a.t();
    out = 0.5 * (out + out.t());
    if (j + 1 < steps) exp_a = exp_a * exp_a;
  }

  return out;
}

// A semidefinite lower factor, without adding artificial process noise.
arma::mat LowerCholPSD(const arma::mat& x, const std::string& name) {
  const arma::mat sym = 0.5 * (x + x.t());
  if (!sym.is_finite()) Rcpp::stop(name + " must be finite.");
  const arma::uword n = sym.n_rows;
  arma::mat out;
  if (arma::chol(out, sym, "lower")) return out;
  const double scale = arma::norm(sym, "inf");
  out.zeros(n, n);
  if (scale == 0.0) return out;
  const double tol = 100.0 * std::numeric_limits<double>::epsilon() * n * scale;
  arma::vec values;
  arma::mat vectors;
  if (!arma::eig_sym(values, vectors, sym) || values.min() < -tol) {
    Rcpp::stop(name + " is not positive semidefinite.");
  }
  // Unpivoted Cholesky amplifies rounding on rank-deficient input. Instead,
  // factor sym = root * root' with clamped eigenvalues; the QR step
  // root' = Q R gives sym = R' R, so R' is a lower triangular factor.
  const arma::mat root =
      vectors * arma::diagmat(arma::sqrt(arma::clamp(values, 0.0, scale)));
  arma::mat q, r;
  if (!arma::qr(q, r, root.t())) {
    Rcpp::stop(name + " factorization failed.");  // # nocov
  }
  out = arma::trimatl(r.t());
  for (arma::uword j = 0; j < n; ++j) {
    if (out(j, j) < 0.0) out.col(j) *= -1.0;
  }
  if (arma::norm(out * out.t() - sym, "inf") > 10.0 * tol) {
    Rcpp::stop(name + " factorization residual is too large.");  // # nocov
  }
  return out;
}

// Integral from 0 to h of exp(A s), valid even when A is singular.
arma::mat IntegratedExponential(const arma::mat& a, const double h) {
  const arma::uword n = a.n_rows;
  arma::mat block(2 * n, 2 * n, arma::fill::zeros);
  block.submat(0, 0, n - 1, n - 1) = a;
  block.submat(0, n, n - 1, 2 * n - 1) = arma::eye(n, n);
  arma::mat result;
  if (!arma::expmat(result, block * h) || !result.is_finite()) {
    Rcpp::stop("Integrated matrix exponential failed.");  // # nocov
  }
  return result.submat(0, n, n - 1, 2 * n - 1);
}

// Exact discretization of d eta = (iota + phi eta) dt + sigma_l dW over
// delta_t. The only discretization path: valid for a singular drift, a
// semidefinite diffusion, and eigenvalue pairs with lambda_i + lambda_j = 0.
void DiscretizeLinSDE(const arma::vec& iota, const arma::mat& phi,
                      const arma::mat& sigma_l, const double delta_t,
                      arma::vec& alpha, arma::mat& beta, arma::mat& psi_l) {
  if (!arma::expmat(beta, phi * delta_t) || !beta.is_finite()) {
    Rcpp::stop("Matrix exponential of the drift failed.");  // # nocov
  }
  alpha = iota.is_zero() ? iota : IntegratedExponential(phi, delta_t) * iota;
  psi_l = sigma_l.is_zero()
              ? arma::mat(sigma_l)
              : LowerCholPSD(
                    AccumulatedCovariance(phi, sigma_l * sigma_l.t(), delta_t),
                    "psi");
}

// Derivative state is (y, y', ..., y^(p-1)), unlike the VAR lag companion.
arma::mat DerivativeCompanion(const Rcpp::List& phi) {
  const arma::mat first = Rcpp::as<arma::mat>(phi[0]);
  const arma::uword q = first.n_rows, p = phi.size(), n = p * q;
  arma::mat a(n, n, arma::fill::zeros);
  for (arma::uword j = 0; j < p; ++j) {
    a.submat(n - q, j * q, n - 1, (j + 1) * q - 1) =
        Rcpp::as<arma::mat>(phi[j]);
    if (j + 1 < p)
      a.submat(j * q, (j + 1) * q, (j + 1) * q - 1, (j + 2) * q - 1) =
          arma::eye(q, q);
  }
  return a;
}
