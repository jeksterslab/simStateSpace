// -----------------------------------------------------------------------------
// edit .setup/cpp/000-forward-declarations.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

bool TestPhi(const arma::mat& phi, const double margin,
             const double auto_ubound);

bool TestStability(const arma::mat& x, const double margin);

bool TestStationarity(const arma::mat& x, const double margin);

double SpectralRadius(const arma::mat& x);

arma::mat ProjectToStability(const arma::mat& x, const double margin,
                             const double tol);

double SpectralAbscissa(const arma::mat& x);

arma::mat ProjectToHurwitz(const arma::mat& x, const double margin);

Rcpp::List LinSDE2SSM(const arma::vec&, const arma::mat&, const arma::mat&,
                      double);
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
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-companion-mat-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.CompanionMat)]]
arma::mat CompanionMat(const Rcpp::List& x) {
  const R_xlen_t order_r = x.size();
  Rcpp::NumericMatrix x_1_r(x[0]);
  const arma::uword q = static_cast<arma::uword>(x_1_r.nrow());
  const arma::uword order = static_cast<arma::uword>(order_r);
  const arma::uword n = q * order;
  arma::mat out = arma::zeros(n, n);
  const arma::mat identity = arma::eye(q, q);
  for (arma::uword i = 0; i < order; i++) {
    Rcpp::NumericMatrix x_i_r(x[i]);
    const arma::mat x_i(x_i_r.begin(), x_i_r.nrow(), x_i_r.ncol(), false);
    out.submat(0, i * q, q - 1, ((i + 1) * q) - 1) = x_i;
  }
  for (arma::uword i = 1; i < order; i++) {
    out.submat(i * q, (i - 1) * q, ((i + 1) * q) - 1, (i * q) - 1) = identity;
  }
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-exact-ctvar-2-ssm-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.CTVAR2SSM)]]
Rcpp::List CTVAR2SSM(const arma::vec& iota, const Rcpp::List& phi,
                     const arma::mat& sigma_l, const double delta_t) {
  const arma::mat a = DerivativeCompanion(phi);
  const arma::uword n = a.n_rows, q = iota.n_elem;
  arma::vec b(n, arma::fill::zeros);
  b.tail(q) = iota;
  arma::mat g(n, n, arma::fill::zeros);
  g.submat(n - q, n - q, n - 1, n - 1) = sigma_l;
  Rcpp::List out = LinSDE2SSM(b, a, g, delta_t);
  arma::mat lambda(q, n, arma::fill::zeros);
  lambda.submat(0, 0, q - 1, q - 1) = arma::eye(q, q);
  out["lambda"] = lambda;
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-2-ssm-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDE2SSM)]]
Rcpp::List LinSDE2SSM(const arma::vec& iota, const arma::mat& phi,
                      const arma::mat& sigma_l, const double delta_t) {
  arma::vec alpha;
  arma::mat beta, psi_l;
  DiscretizeLinSDE(iota, phi, sigma_l, delta_t, alpha, beta, psi_l);
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi_l") = psi_l);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-exact-ssm-2-lin-sde-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSM2LinSDE)]]
Rcpp::List SSM2LinSDE(const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_l, const double delta_t,
                      const double tol) {
  const arma::uword p = beta.n_rows;
  arma::mat I = arma::eye(p, p);
  arma::cx_mat log_beta;
  if (!arma::logmat(log_beta, beta)) {
    Rcpp::stop("Matrix logarithm of beta failed.");
  }
  if (!log_beta.is_finite()) {
    Rcpp::stop("Matrix logarithm of beta is not finite.");
  }
  const double max_imag = arma::abs(arma::imag(log_beta)).max();
  if (!(max_imag <= tol)) {
    Rcpp::stop("Matrix logarithm of beta has a non-negligible imaginary part.");
  }
  arma::mat phi = arma::real(log_beta) / delta_t;
  arma::vec iota;
  if (!arma::solve(iota, IntegratedExponential(phi, delta_t), alpha,
                   arma::solve_opts::no_approx)) {
    // # nocov start
    Rcpp::stop("Continuous-time intercept is not uniquely recoverable.");
    // # nocov end
  }
  const arma::mat psi = psi_l * psi_l.t();
  const arma::mat generator = arma::kron(I, phi) + arma::kron(phi, I);
  arma::vec sigma_vec;
  if (!arma::solve(sigma_vec, IntegratedExponential(generator, delta_t),
                   arma::vectorise(psi), arma::solve_opts::no_approx)) {
    // # nocov start
    Rcpp::stop("Continuous-time diffusion is not uniquely recoverable.");
    // # nocov end
  }
  const arma::mat sigma_l =
      LowerCholPSD(arma::reshape(sigma_vec, p, p), "sigma");
  return Rcpp::List::create(Rcpp::Named("iota") = iota,
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("sigma_l") = sigma_l);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-exp-mat-log-mat.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix Exponential
//'
//' Computes the matrix exponential of a square numeric matrix using the
//' Armadillo linear algebra library.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x A square numeric matrix.
//'
//' @return
//' A numeric matrix containing the matrix exponential of `x`.
//'
//' @examples
//' a <- matrix(
//'   data = c(
//'     -0.40,  0.10,
//'     -0.05, -0.30
//'   ),
//'   nrow = 2,
//'   ncol = 2
//' )
//'
//' b <- ExpMat(x = a)
//' b
//'
//' a <- LogMat(x = b)
//' a
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
arma::mat ExpMat(const arma::mat& x) {
  if (x.n_rows != x.n_cols) {
    Rcpp::stop("`x` must be a square matrix.");
  }
  arma::mat out;
  const bool success = arma::expmat(out, x);
  if (!success) {
    // # nocov start
    Rcpp::stop("Armadillo failed to compute the matrix exponential.");
    // # nocov end
  }
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-2-ssm-sets-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]

// Discretize the individual-varying linear SDE once per parameter set.
// Person i uses iota, phi, and sigma_l at i modulo each list's length, so
// the sets repeat with the least common multiple m of the lengths (at most
// n), and person i uses set i modulo m. With `ou = true`, `iota` holds mu
// and the intercept is -phi mu.
// [[Rcpp::export(.LinSDE2SSMSets)]]
Rcpp::List LinSDE2SSMSets(const arma::uword n, const Rcpp::List& iota,
                          const Rcpp::List& phi, const Rcpp::List& sigma_l,
                          const double delta_t, const bool ou) {
  arma::uword m = 1;
  for (const arma::uword length : {static_cast<arma::uword>(iota.size()),
                                   static_cast<arma::uword>(phi.size()),
                                   static_cast<arma::uword>(sigma_l.size())}) {
    m = std::lcm(m, length);
    if (m >= n) {
      m = n;
      break;
    }
  }
  Rcpp::List alpha(m), beta(m), psi_l(m);
  for (arma::uword s = 0; s < m; ++s) {
    const arma::mat phi_s = phi[s % phi.size()];
    arma::vec iota_s = Rcpp::as<arma::vec>(iota[s % iota.size()]);
    if (ou && !iota_s.is_zero()) iota_s = -phi_s * iota_s;
    arma::vec alpha_s;
    arma::mat beta_s, psi_l_s;
    DiscretizeLinSDE(iota_s, phi_s,
                     Rcpp::as<arma::mat>(sigma_l[s % sigma_l.size()]), delta_t,
                     alpha_s, beta_s, psi_l_s);
    alpha[s] = alpha_s;
    beta[s] = beta_s;
    psi_l[s] = psi_l_s;
  }
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi_l") = psi_l);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-cov-delta-t-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDECovDeltaT)]]
arma::mat LinSDECovDeltaT(const arma::mat& phi, const arma::mat& sigma,
                          const double delta_t) {
  return AccumulatedCovariance(phi, sigma, delta_t);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-cov-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDECovEta)]]
arma::mat LinSDECovEta(const arma::mat phi, const arma::mat sigma) {
  arma::mat X;
  arma::sylvester(X, phi, phi.t(), sigma);
  return ((X + X.t()) / 2);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-cov-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDECovY)]]
arma::mat LinSDECovY(const arma::mat& lambda, const arma::mat& theta,
                     const arma::mat& cov_eta) {
  arma::mat X = lambda * cov_eta * lambda.t() + theta;
  return ((X + X.t()) / 2);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-intercept-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDEInterceptEta)]]
Rcpp::NumericVector LinSDEInterceptEta(const arma::mat& phi,
                                       const arma::vec& mean_eta) {
  arma::vec output = -phi * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-mean-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDEMeanEta)]]
Rcpp::NumericVector LinSDEMeanEta(const arma::mat phi, const arma::vec iota) {
  arma::vec output = arma::solve(-phi, iota);
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-mean-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDEMeanY)]]
Rcpp::NumericVector LinSDEMeanY(const arma::vec& nu, const arma::mat& lambda,
                                const arma::vec& mean_eta) {
  arma::vec output = nu + lambda * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-log-mat.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix Logarithm
//'
//' Computes the real-valued principal matrix logarithm of a square numeric
//' matrix using the Armadillo linear algebra library.
//'
//' @param x A square numeric matrix.
//' @param tol A positive numeric scalar giving the tolerance used to determine
//'   whether the imaginary component of the matrix logarithm is negligible.
//'   The default is `1e-10`.
//'
//' @return
//' A numeric matrix containing the real part of the matrix logarithm of `x`.
//'
//' @examples
//' a <- matrix(
//'   data = c(
//'     -0.40,  0.10,
//'     -0.05, -0.30
//'   ),
//'   nrow = 2,
//'   ncol = 2
//' )
//'
//' b <- ExpMat(x = a)
//' b
//'
//' a <- LogMat(x = b)
//' a
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
arma::mat LogMat(const arma::mat& x, const double tol = 1e-10) {
  if (x.n_rows != x.n_cols) {
    Rcpp::stop("`x` must be a square matrix.");
  }
  if (!x.is_finite()) {
    Rcpp::stop("`x` must contain only finite values.");
  }
  if (!std::isfinite(tol)) {
    Rcpp::stop("`tol` must be finite.");
  }
  if (tol <= 0.0) {
    Rcpp::stop("`tol` must be positive.");
  }
  arma::cx_mat out;
  const bool success = arma::logmat(out, x);
  if (!success) {
    Rcpp::stop("Armadillo failed to compute the matrix logarithm.");  // # nocov
  }
  const double imag_norm = arma::norm(arma::imag(out), "inf");
  if (!std::isfinite(imag_norm) || imag_norm > tol) {
    Rcpp::stop("Matrix logarithm has a non-negligible imaginary component.");
  }
  return arma::real(out);
}
// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-project-to-hurwitz.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Project Matrix to Hurwitz Stability
//'
//' Shifts a square matrix left on the real axis
//' so that its spectral abscissa (maximum real part of the eigenvalues)
//' is strictly less than `-margin`.
//' This is useful for ensuring that continuous-time drift matrices
//' (e.g. in linear SDEs/state-space models) are Hurwitz-stable.
//' If the matrix already satisfies the margin,
//' it is returned unchanged.
//'
//' The projection is performed by subtracting a multiple of the identity:
//' \deqn{x^\star = x - (\alpha + \text{margin}) I,}
//' where \eqn{\alpha = \max \Re\{\lambda_i(x)\}} is the spectral abscissa.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//' @param margin Positive numeric.
//'   Target buffer inside the Hurwitz region;
//'   the result satisfies
//'   \eqn{\max \Re\{\lambda_i(x^\star)\} \le -\text{margin}}
//'   (default `1e-3`).
//'
//' @return A numeric matrix of the same dimensions as `x`,
//'   shifted if necessary to satisfy the Hurwitz stability constraint.
//'
//' @examples
//' # Unstable (spectral abscissa >= 0):
//' x <- matrix(
//'   data = c(
//'     0.10, -0.40,
//'     0.50, 0.20
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # >= 0
//' SpectralAbscissa(x = ProjectToHurwitz(x = x)) # <= -1e-3 (default margin)
//'
//' # Already Hurwitz-stable is returned unchanged up to numerics:
//' x <- matrix(
//'   data = c(
//'     -0.50, -0.20,
//'      1.00, -0.30
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # < 0
//' identical(ProjectToHurwitz(x = x), x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace stability linsde
//' @export
// [[Rcpp::export]]
arma::mat ProjectToHurwitz(const arma::mat& x, const double margin = 1e-3) {
  arma::mat output = x;
  const double alpha = SpectralAbscissa(output);
  if (alpha >= -margin) {
    output -= (alpha + margin) * arma::eye(output.n_rows, output.n_cols);
  }
  return output;
}
// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-project-to-stability.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Project Matrix to Stability
//'
//' Scales a square matrix
//' so that its spectral radius
//' is strictly less than 1 by a specified stability margin.
//' This is useful for ensuring that transition matrices
//' in state space or vector autoregressive (VAR) models are stationary.
//' If the matrix is already within the margin,
//' it is returned unchanged.
//'
//' The projection is performed
//' by multiplying the matrix by a constant factor
//' \eqn{c = \frac{\text{margin}}{\rho + \text{tol}}},
//' where \eqn{\rho} is the spectral radius and `tol`
//' is a small positive number to prevent division by zero.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//' @param margin Double in \eqn{(0, 1)}.
//'   Target upper bound for the spectral
//'   radius (default = 0.98).
//' @param tol Small positive double
//'   added to the denominator in the scaling
//'   factor to avoid division by zero
//'   (default `1e-12`).
//'
//' @return A numeric matrix of the same dimensions as `x`,
//'   scaled if necessary to satisfy the stability constraint.
//'
//' @examples
//' # Matrix with eigenvalues greater than 1
//' x <- matrix(
//'   data = c(
//'     1.2, 0.3,
//'     0.4, 0.9
//'   ),
//'   nrow = 2
//' )
//' SpectralRadius(x = x) # > 1
//' SpectralRadius(x = ProjectToStability(x = x))  # < 1
//'
//' # Matrix already stable is returned unchanged
//' x <- matrix(
//'   data = c(
//'     0.5, 0.3,
//'     0.2, 0.4
//'   ),
//'   nrow = 2
//' )
//' identical(ProjectToStability(x = x), x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace stability ssm
//' @export
// [[Rcpp::export]]
arma::mat ProjectToStability(const arma::mat& x, const double margin = 0.98,
                             const double tol = 1e-12) {
  double rho = SpectralRadius(x);
  if (rho < margin) return x;
  double c = margin / (rho + tol);
  return c * x;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n-0-inline.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

inline bool in_bounds_element(double x, double lb, double ub, bool has_lb,
                              bool has_ub) {
  // Treat non-finite bounds (NA/NaN/Inf) as "no bound" on that side
  bool lb_ok = true, ub_ok = true;
  if (has_lb && std::isfinite(lb)) lb_ok = (x >= lb);
  if (has_ub && std::isfinite(ub)) ub_ok = (x <= ub);
  return lb_ok && ub_ok;
}

inline bool matrix_in_bounds(const arma::mat& x, const arma::mat* lb_ptr,
                             const arma::mat* ub_ptr, bool has_lb,
                             bool has_ub) {
  if (!has_lb && !has_ub) return true;  // nothing to check
  const arma::uword nr = x.n_rows, nc = x.n_cols;
  for (arma::uword i = 0; i < nr; ++i) {
    for (arma::uword j = 0; j < nc; ++j) {
      double lb = has_lb ? (*lb_ptr)(i, j) : NA_REAL;
      double ub = has_ub ? (*ub_ptr)(i, j) : NA_REAL;
      if (!in_bounds_element(x(i, j), lb, ub, has_lb, has_ub)) return false;
    }
  }
  return true;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n-2.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Transition Matrices
//' from the Multivariate Normal Distribution
//' and Project to Stability
//'
//' This function simulates random transition matrices
//' from the multivariate normal distribution
//' then projects each draw to the stability region
//' using [ProjectToStability()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param beta Numeric matrix.
//'   The transition matrix (\eqn{\boldsymbol{\beta}}).
//' @param vcov_beta_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_beta_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\beta} \right)}.
//' @param margin Double in \eqn{(0, 1)}.
//'   Target upper bound for the spectral radius
//'   (default = 0.98).
//' @param tol Small positive double added to the denominator
//'   in the scaling factor to avoid division by zero
//'   (default = 1e-12).
//' @return Returns a list of random transition matrices.
//'
//' @examples
//' n <- 10
//' beta <- matrix(
//'   data = c(
//'     0.7, 0.5, -0.1,
//'     0.0, 0.6, 0.4,
//'     0, 0, 0.5
//'   ),
//'   nrow = 3
//' )
//' vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
//' SimBetaN2(n = n, beta = beta, vcov_beta_vec_l = vcov_beta_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimBetaN2(const arma::uword& n, const arma::mat& beta,
                     const arma::mat& vcov_beta_vec_l,
                     const double margin = 0.98, const double tol = 1e-12) {
  Rcpp::List output(n);
  arma::vec beta_vec = arma::vectorise(beta);
  arma::vec beta_vec_i(beta.n_rows * beta.n_cols, arma::fill::none);
  arma::mat beta_i(beta.n_rows, beta.n_cols, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    beta_vec_i =
        beta_vec + (vcov_beta_vec_l * arma::randn(beta.n_rows * beta.n_cols));
    beta_i = arma::reshape(beta_vec_i, beta.n_rows, beta.n_cols);
    beta_i = ProjectToStability(beta_i, margin, tol);
    output[i] = beta_i;
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n-covariate.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

//' Simulate Transition Matrices with a Covariate
//' from the Multivariate Normal Distribution
//'
//' This function simulates random transition matrices from a multivariate
//' normal distribution, allowing the mean transition matrix to vary as a
//' linear function of a covariate.
//' The function ensures that the generated transition matrices are stationary
//' using [TestStationarity()] with a rejection sampling approach.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param beta0 Numeric matrix.
//'   Baseline transition matrix \eqn{\boldsymbol{\beta}_0}
//'   corresponding to \eqn{\mathbf{x} = \mathbf{0}}.
//' @param vcov_beta_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_beta_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\beta} \right)}.
//' @param beta1 Numeric matrix.
//'   Matrix of covariate effects mapping \eqn{\mathbf{x}} to
//'   \eqn{\mathrm{vec}(\boldsymbol{\beta})}.
//' @param x List of numeric vectors.
//'   Covariate values.
//' @param margin Numeric scalar specifying the stationarity threshold.
//'   Values less than 1 indicate stricter stationarity criteria.
//' @param beta_lbound Optional numeric matrix of same dim as `beta`.
//'   Use NA for no lower bound.
//' @param beta_ubound Optional numeric matrix of same dim as `beta`.
//'   Use NA for no upper bound.
//' @param bound Logical;
//'   if TRUE, resample until all elements respect bounds (NA bounds ignored).
//' @param max_iter Safety cap on resampling attempts per draw.
//' @return Returns a list of random transition matrices.
//'
//' @examples
//' n <- 5
//' beta0 <- matrix(
//'   data = c(
//'     0.7, 0.5, -0.1,
//'     0.0, 0.6, 0.4,
//'     0, 0, 0.5
//'   ),
//'   nrow = 3
//' )
//' vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
//' # One scalar covariate per replication
//' beta1 <- matrix(data = 0, nrow = 9, ncol = 1)
//' beta1[1, 1] <- 0.10  # x shifts beta[1,1]
//' x <- list(c(0), c(1), c(-1), c(0.5), c(2))
//'
//' SimBetaNCovariate(
//'   n = n,
//'   beta0 = beta0,
//'   vcov_beta_vec_l = vcov_beta_vec_l,
//'   beta1 = beta1,
//'   x = x
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimBetaNCovariate(
    const arma::uword& n, const arma::mat& beta0,
    const arma::mat& vcov_beta_vec_l, const arma::mat& beta1,
    const Rcpp::List& x, const double margin = 1.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> beta_lbound = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> beta_ubound = R_NilValue,
    const bool bound = false, const arma::uword max_iter = 100000) {
  const arma::uword nr = beta0.n_rows, nc = beta0.n_cols;
  const arma::uword p = nr * nc;

  arma::mat lb, ub;
  arma::umat has_lb_el(nr, nc, arma::fill::zeros);
  arma::umat has_ub_el(nr, nc, arma::fill::zeros);

  const bool has_lb = beta_lbound.isNotNull();
  const bool has_ub = beta_ubound.isNotNull();

  if (has_lb) {
    lb = Rcpp::as<arma::mat>(beta_lbound);
    if (lb.n_rows != nr || lb.n_cols != nc)
      Rcpp::stop("beta_lbound dims must match beta0.");
    // Build mask: 1 where finite lower bound exists
    for (arma::uword i = 0; i < nr; ++i) {
      for (arma::uword j = 0; j < nc; ++j) {
        has_lb_el(i, j) = std::isfinite(lb(i, j)) ? 1u : 0u;
      }
    }
  } else {
    lb.set_size(nr, nc);
    lb.fill(0.0);
  }

  if (has_ub) {
    ub = Rcpp::as<arma::mat>(beta_ubound);
    if (ub.n_rows != nr || ub.n_cols != nc)
      Rcpp::stop("beta_ubound dims must match beta0.");
    // Build mask: 1 where finite upper bound exists
    for (arma::uword i = 0; i < nr; ++i) {
      for (arma::uword j = 0; j < nc; ++j) {
        has_ub_el(i, j) = std::isfinite(ub(i, j)) ? 1u : 0u;
      }
    }
  } else {
    ub.set_size(nr, nc);
    ub.fill(0.0);
  }

  // Vectorized bounds check with masks; quick-reject
  auto bounds_ok = [&](const arma::mat& x) -> bool {
    if (!bound) return true;
    arma::umat low_violate = (x < lb) % has_lb_el;
    arma::umat high_violate = (x > ub) % has_ub_el;
    return !(arma::any(arma::vectorise(low_violate)) ||
             arma::any(arma::vectorise(high_violate)));
  };

  Rcpp::List out(n);
  const arma::vec beta_vec = arma::vectorise(beta0);
  arma::vec z(p, arma::fill::none), beta_vec_i(p, arma::fill::none);
  arma::mat beta_i(nr, nc, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword iter = 0;
    arma::vec x_i = x[i];
    for (;;) {
      if (iter++ >= max_iter) {
        Rcpp::stop(
            "SimBetaNCovariate: exceeded max_iter while drawing beta_i (i=%u). "
            "Relax bounds or stationarity.",
            i + 1);
      }
      z.randn();
      beta_vec_i = beta_vec + (beta1 * x_i) + (vcov_beta_vec_l * z);
      beta_i = arma::reshape(beta_vec_i, nr, nc);

      if (!bounds_ok(beta_i)) continue;

      if (!TestStationarity(beta_i, margin)) continue;

      out[i] = beta_i;
      break;
    }
  }

  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

//' Simulate Transition Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random transition matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated transition matrices are stationary
//' using [TestStationarity()] with a rejection sampling approach.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param beta Numeric matrix.
//'   The transition matrix (\eqn{\boldsymbol{\beta}}).
//' @param vcov_beta_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_beta_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\beta} \right)}.
//' @param margin Numeric scalar specifying the stationarity threshold.
//'   Values less than 1 indicate stricter stationarity criteria.
//' @param beta_lbound Optional numeric matrix of same dim as `beta`.
//'   Use NA for no lower bound.
//' @param beta_ubound Optional numeric matrix of same dim as `beta`.
//'   Use NA for no upper bound.
//' @param bound Logical;
//'   if TRUE, resample until all elements respect bounds (NA bounds ignored).
//' @param max_iter Safety cap on resampling attempts per draw.
//' @return Returns a list of random transition matrices.
//'
//' @examples
//' n <- 10
//' beta <- matrix(
//'   data = c(
//'     0.7, 0.5, -0.1,
//'     0.0, 0.6, 0.4,
//'     0, 0, 0.5
//'   ),
//'   nrow = 3
//' )
//' vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
//' SimBetaN(n = n, beta = beta, vcov_beta_vec_l = vcov_beta_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimBetaN(
    const arma::uword& n, const arma::mat& beta,
    const arma::mat& vcov_beta_vec_l, const double margin = 1.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> beta_lbound = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> beta_ubound = R_NilValue,
    const bool bound = false, const arma::uword max_iter = 100000) {
  const arma::uword nr = beta.n_rows, nc = beta.n_cols;
  const arma::uword p = nr * nc;

  arma::mat lb, ub;
  arma::umat has_lb_el(nr, nc, arma::fill::zeros);
  arma::umat has_ub_el(nr, nc, arma::fill::zeros);

  const bool has_lb = beta_lbound.isNotNull();
  const bool has_ub = beta_ubound.isNotNull();

  if (has_lb) {
    lb = Rcpp::as<arma::mat>(beta_lbound);
    if (lb.n_rows != nr || lb.n_cols != nc)
      Rcpp::stop("beta_lbound dims must match beta.");
    // Build mask: 1 where finite lower bound exists
    for (arma::uword i = 0; i < nr; ++i) {
      for (arma::uword j = 0; j < nc; ++j) {
        has_lb_el(i, j) = std::isfinite(lb(i, j)) ? 1u : 0u;
      }
    }
  } else {
    lb.set_size(nr, nc);
    lb.fill(0.0);
  }

  if (has_ub) {
    ub = Rcpp::as<arma::mat>(beta_ubound);
    if (ub.n_rows != nr || ub.n_cols != nc)
      Rcpp::stop("beta_ubound dims must match beta.");
    // Build mask: 1 where finite upper bound exists
    for (arma::uword i = 0; i < nr; ++i) {
      for (arma::uword j = 0; j < nc; ++j) {
        has_ub_el(i, j) = std::isfinite(ub(i, j)) ? 1u : 0u;
      }
    }
  } else {
    ub.set_size(nr, nc);
    ub.fill(0.0);
  }

  // Vectorized bounds check with masks; quick-reject
  auto bounds_ok = [&](const arma::mat& x) -> bool {
    if (!bound) return true;
    arma::umat low_violate = (x < lb) % has_lb_el;
    arma::umat high_violate = (x > ub) % has_ub_el;
    return !(arma::any(arma::vectorise(low_violate)) ||
             arma::any(arma::vectorise(high_violate)));
  };

  Rcpp::List out(n);
  const arma::vec beta_vec = arma::vectorise(beta);
  arma::vec z(p, arma::fill::none), beta_vec_i(p, arma::fill::none);
  arma::mat beta_i(nr, nc, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword iter = 0;
    for (;;) {
      if (iter++ >= max_iter) {
        Rcpp::stop(
            "SimBetaN: exceeded max_iter while drawing beta_i (i=%u). Relax "
            "bounds or stationarity.",
            i + 1);
      }
      z.randn();
      beta_vec_i = beta_vec + (vcov_beta_vec_l * z);
      beta_i = arma::reshape(beta_vec_i, nr, nc);

      if (!bounds_ok(beta_i)) continue;

      if (!TestStationarity(beta_i, margin)) continue;

      out[i] = beta_i;
      break;
    }
  }

  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-cov-diag-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Diagonal Covariance Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random diagonal covariance matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated covariance matrices
//' are positive semi-definite.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param sigma_diag Numeric matrix.
//'   The covariance matrix
//'   (\eqn{\boldsymbol{\Sigma}}).
//' @param vcov_sigma_diag_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_sigma_vech))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vech} \left( \boldsymbol{\Sigma} \right)}.
//' @return Returns a list of random diagonal covariance matrices.
//'
//' @examples
//' n <- 10
//' sigma_diag <- c(1, 1, 1)
//' vcov_sigma_diag_l <- t(chol(0.001 * diag(3)))
//' SimCovDiagN(
//'   n = n,
//'   sigma_diag = sigma_diag,
//'   vcov_sigma_diag_l = vcov_sigma_diag_l
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimCovDiagN(const arma::uword& n, const arma::vec& sigma_diag,
                       const arma::mat& vcov_sigma_diag_l) {
  Rcpp::List output(n);
  arma::uword p = sigma_diag.n_rows;
  arma::vec sigma_diag_i(p, arma::fill::none);
  arma::mat sigma_i(p, p, arma::fill::zeros);  // diagonal only
  for (arma::uword i = 0; i < n; i++) {
    // generate data
    sigma_diag_i = sigma_diag + (vcov_sigma_diag_l * arma::randn(p));
    // make positive semi definite
    sigma_diag_i.transform([](double val) { return std::max(val, 1e-8); });
    sigma_i.diag() = sigma_diag_i;
    output[i] = sigma_i;
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-cov-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Covariance Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random covariance matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated covariance matrices
//' are positive semi-definite.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param sigma Numeric matrix.
//'   The covariance matrix
//'   (\eqn{\boldsymbol{\Sigma}}).
//' @param vcov_sigma_vech_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_sigma_vech))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vech} \left( \boldsymbol{\Sigma} \right)}.
//' @return Returns a list of random covariance matrices.
//'
//' @examples
//' n <- 10
//' sigma <- matrix(
//'   data = c(
//'     1.0, 0.5, 0.3,
//'     0.5, 1.0, 0.4,
//'     0.3, 0.4, 1.0
//'   ),
//'   nrow = 3
//' )
//' vcov_sigma_vech_l <- t(
//'   chol(
//'     0.001 * diag(3 * (3 + 1) / 2)
//'   )
//' )
//' SimCovN(
//'   n = n,
//'   sigma = sigma,
//'   vcov_sigma_vech_l = vcov_sigma_vech_l
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimCovN(const arma::uword& n, const arma::mat& sigma,
                   const arma::mat& vcov_sigma_vech_l) {
  Rcpp::List output(n);
  arma::uword p = sigma.n_rows;
  arma::uword q = p * (p + 1) / 2;
  arma::vec vech(q, arma::fill::none);
  arma::uword idx = 0;
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = j; i < p; ++i) {
      vech(idx++) = sigma(i, j);
    }
  }
  arma::mat sigma_i(p, p, arma::fill::none);
  arma::vec vech_i(q, arma::fill::none);
  arma::vec eigval;
  arma::mat eigvec;
  for (arma::uword i = 0; i < n; i++) {
    // generate data
    vech_i = vech + (vcov_sigma_vech_l * arma::randn(q));
    // make positive semi definite
    idx = 0;
    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword j = i; j < p; ++j) {
        sigma_i(i, j) = vech_i(idx);
        sigma_i(j, i) = vech_i(idx);
        idx++;
      }
    }
    arma::eig_sym(eigval, eigvec, sigma_i);
    eigval.transform([](double val) { return std::max(val, 1e-8); });
    sigma_i = eigvec * arma::diagmat(eigval) * eigvec.t();
    output[i] = sigma_i;
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-mvn.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimMVN)]]
Rcpp::List SimMVN(const arma::uword& n, const arma::vec& mu,
                  const arma::mat& sigma_l) {
  Rcpp::List output(n);
  arma::vec mu_i(mu.n_rows, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    mu_i = mu + (sigma_l * arma::randn(mu.n_rows));
    output[i] = Rcpp::NumericVector(mu_i.begin(), mu_i.end());
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi-n-2.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices
//' from the Multivariate Normal Distribution
//' and Project to Hurwitz
//'
//' This function simulates random dirft matrices
//' from the multivariate normal distribution
//' then projects each draw to the Hurwitz-stable region
//' using [ProjectToHurwitz()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//' @param margin Positive numeric.
//'   Target buffer so that the spectral abscissa
//'   is \eqn{\le -\text{margin}} (default `1e-3`).
//' @return Returns a list of random drift matrices.
//'
//' @examples
//' n <- 10
//' phi <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
//' SimPhiN2(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhiN2(const arma::uword& n, const arma::mat& phi,
                    const arma::mat& vcov_phi_vec_l,
                    const double margin = 1e-3) {
  Rcpp::List output(n);
  arma::vec phi_vec = arma::vectorise(phi);
  arma::vec phi_vec_i(phi.n_rows * phi.n_cols, arma::fill::none);
  arma::mat phi_i(phi.n_rows, phi.n_cols, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    phi_vec_i =
        phi_vec + (vcov_phi_vec_l * arma::randn(phi.n_rows * phi.n_cols));
    phi_i = arma::reshape(phi_vec_i, phi.n_rows, phi.n_cols);
    phi_i = ProjectToHurwitz(phi_i, margin);
    output[i] = phi_i;
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi-n-covariate.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices with a Covariate
//' from the Multivariate Normal Distribution
//'
//' This function simulates random drift matrices
//' from the multivariate normal distribution,
//' allowing the mean drift matrix to vary
//' as a linear function of a covariate
//' The function ensures that the generated drift matrices are stable
//' using [TestPhi()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi0 Numeric matrix.
//'   Baseline drift matrix (\eqn{\boldsymbol{\Phi}_0}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//' @param phi1 Numeric matrix.
//'   Matrix of covariate effects mapping \eqn{\mathbf{x}} to
//'   \eqn{\mathrm{vec}(\boldsymbol{\Phi})}.
//' @param x List of numeric vectors.
//'   Covariate values.
//' @param margin Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//' @param auto_ubound Numeric scalar specifying the upper bound
//'   for the diagonal elements of \eqn{\boldsymbol{\Phi}}.
//'   Default is `0.0`, requiring all diagonal values to be \eqn{\leq 0}.
//' @param phi_lbound Optional numeric matrix of same dim as `phi`.
//'   Use NA for no lower bound.
//' @param phi_ubound Optional numeric matrix of same dim as `phi`.
//'   Use NA for no upper bound.
//' @param bound Logical;
//'   if TRUE, resample until all elements respect bounds (NA bounds ignored).
//' @param max_iter Safety cap on resampling attempts per draw.
//' @return Returns a list of random drift matrices.
//'
//' @examples
//' n <- 5
//' phi0 <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
//' # One scalar covariate per replication
//' phi1 <- matrix(data = 0, nrow = 9, ncol = 1)
//' phi1[1, 1] <- 0.10  # x shifts phi[1,1]
//' x <- list(c(0), c(1), c(-1), c(0.5), c(2))
//' SimPhiNCovariate(
//'   n = n,
//'   phi0 = phi0,
//'   vcov_phi_vec_l = vcov_phi_vec_l,
//'   phi1 = phi1,
//'   x = x
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhiNCovariate(
    const arma::uword& n, const arma::mat& phi0,
    const arma::mat& vcov_phi_vec_l, const arma::mat& phi1, const Rcpp::List& x,
    const double margin = 0.0, const double auto_ubound = 0.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> phi_lbound = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> phi_ubound = R_NilValue,
    const bool bound = false, const arma::uword max_iter = 100000) {
  const arma::uword nr = phi0.n_rows, nc = phi0.n_cols;
  const arma::uword p = nr * nc;

  // Bounds & masks
  arma::mat lb, ub;
  arma::umat has_lb_el(nr, nc, arma::fill::zeros);
  arma::umat has_ub_el(nr, nc, arma::fill::zeros);

  const bool has_lb = phi_lbound.isNotNull();
  const bool has_ub = phi_ubound.isNotNull();

  if (has_lb) {
    lb = Rcpp::as<arma::mat>(phi_lbound);
    if (lb.n_rows != nr || lb.n_cols != nc)
      Rcpp::stop("phi_lbound dims must match phi0.");
    for (arma::uword i = 0; i < nr; ++i)
      for (arma::uword j = 0; j < nc; ++j)
        has_lb_el(i, j) = std::isfinite(lb(i, j)) ? 1u : 0u;
  } else {
    lb.set_size(nr, nc);
    lb.fill(0.0);
  }

  if (has_ub) {
    ub = Rcpp::as<arma::mat>(phi_ubound);
    if (ub.n_rows != nr || ub.n_cols != nc)
      Rcpp::stop("phi_ubound dims must match phi0.");
    for (arma::uword i = 0; i < nr; ++i)
      for (arma::uword j = 0; j < nc; ++j)
        has_ub_el(i, j) = std::isfinite(ub(i, j)) ? 1u : 0u;
  } else {
    ub.set_size(nr, nc);
    ub.fill(0.0);
  }

  auto bounds_ok = [&](const arma::mat& x) -> bool {
    if (!bound) return true;
    arma::umat low_violate = (x < lb) % has_lb_el;
    arma::umat high_violate = (x > ub) % has_ub_el;
    return !(arma::any(arma::vectorise(low_violate)) ||
             arma::any(arma::vectorise(high_violate)));
  };

  Rcpp::List output(n);
  const arma::vec phi_vec = arma::vectorise(phi0);
  arma::vec z(p, arma::fill::none), phi_vec_i(p, arma::fill::none);
  arma::mat phi_i(nr, nc, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword iter = 0;
    arma::vec x_i = x[i];
    for (;;) {
      if (iter++ >= max_iter) {
        Rcpp::stop(
            "SimPhiN: exceeded max_iter while drawing phi_i (i=%u). Relax "
            "bounds or TestPhi().",
            i + 1);
      }
      z.randn();
      phi_vec_i = phi_vec + (phi1 * x_i) + (vcov_phi_vec_l * z);
      phi_i = arma::reshape(phi_vec_i, nr, nc);

      if (!bounds_ok(phi_i)) continue;

      if (!TestPhi(phi_i, margin, auto_ubound)) continue;

      output[i] = phi_i;
      break;
    }
  }

  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random drift matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated drift matrices are stable
//' using [TestPhi()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//' @param margin Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//' @param auto_ubound Numeric scalar specifying the upper bound
//'   for the diagonal elements of \eqn{\boldsymbol{\Phi}}.
//'   Default is `0.0`, requiring all diagonal values to be \eqn{\leq 0}.
//' @param phi_lbound Optional numeric matrix of same dim as `phi`.
//'   Use NA for no lower bound.
//' @param phi_ubound Optional numeric matrix of same dim as `phi`.
//'   Use NA for no upper bound.
//' @param bound Logical;
//'   if TRUE, resample until all elements respect bounds (NA bounds ignored).
//' @param max_iter Safety cap on resampling attempts per draw.
//' @return Returns a list of random drift matrices.
//'
//' @examples
//' n <- 10
//' phi <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
//' SimPhiN(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhiN(const arma::uword& n, const arma::mat& phi,
                   const arma::mat& vcov_phi_vec_l, const double margin = 0.0,
                   const double auto_ubound = 0.0,
                   Rcpp::Nullable<Rcpp::NumericMatrix> phi_lbound = R_NilValue,
                   Rcpp::Nullable<Rcpp::NumericMatrix> phi_ubound = R_NilValue,
                   const bool bound = false,
                   const arma::uword max_iter = 100000) {
  const arma::uword nr = phi.n_rows, nc = phi.n_cols;
  const arma::uword p = nr * nc;

  // Bounds & masks
  arma::mat lb, ub;
  arma::umat has_lb_el(nr, nc, arma::fill::zeros);
  arma::umat has_ub_el(nr, nc, arma::fill::zeros);

  const bool has_lb = phi_lbound.isNotNull();
  const bool has_ub = phi_ubound.isNotNull();

  if (has_lb) {
    lb = Rcpp::as<arma::mat>(phi_lbound);
    if (lb.n_rows != nr || lb.n_cols != nc)
      Rcpp::stop("phi_lbound dims must match phi.");
    for (arma::uword i = 0; i < nr; ++i)
      for (arma::uword j = 0; j < nc; ++j)
        has_lb_el(i, j) = std::isfinite(lb(i, j)) ? 1u : 0u;
  } else {
    lb.set_size(nr, nc);
    lb.fill(0.0);
  }

  if (has_ub) {
    ub = Rcpp::as<arma::mat>(phi_ubound);
    if (ub.n_rows != nr || ub.n_cols != nc)
      Rcpp::stop("phi_ubound dims must match phi.");
    for (arma::uword i = 0; i < nr; ++i)
      for (arma::uword j = 0; j < nc; ++j)
        has_ub_el(i, j) = std::isfinite(ub(i, j)) ? 1u : 0u;
  } else {
    ub.set_size(nr, nc);
    ub.fill(0.0);
  }

  auto bounds_ok = [&](const arma::mat& x) -> bool {
    if (!bound) return true;
    arma::umat low_violate = (x < lb) % has_lb_el;
    arma::umat high_violate = (x > ub) % has_ub_el;
    return !(arma::any(arma::vectorise(low_violate)) ||
             arma::any(arma::vectorise(high_violate)));
  };

  Rcpp::List output(n);
  const arma::vec phi_vec = arma::vectorise(phi);
  arma::vec z(p, arma::fill::none), phi_vec_i(p, arma::fill::none);
  arma::mat phi_i(nr, nc, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword iter = 0;
    for (;;) {
      if (iter++ >= max_iter) {
        Rcpp::stop(
            "SimPhiN: exceeded max_iter while drawing phi_i (i=%u). Relax "
            "bounds or TestPhi().",
            i + 1);
      }
      z.randn();
      phi_vec_i = phi_vec + (vcov_phi_vec_l * z);
      phi_i = arma::reshape(phi_vec_i, nr, nc);

      if (!bounds_ok(phi_i)) continue;

      if (!TestPhi(phi_i, margin, auto_ubound)) continue;

      output[i] = phi_i;
      break;
    }
  }

  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-kernel-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// Parameter sets as Armadillo matrices, converted once;
// vectors become column vectors.
std::vector<arma::mat> SimSSMSets(const Rcpp::List& x) {
  std::vector<arma::mat> out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (Rf_isMatrix(x[i])) {
      out[i] = Rcpp::as<arma::mat>(x[i]);
    } else {
      out[i] = Rcpp::as<arma::vec>(x[i]);
    }
  }
  return out;
}

// One simulation kernel for every state space model.
// Each parameter is a list of parameter sets;
// person i uses set i modulo the length of the list.
// Latent-only models pass empty `nu`, `lambda`, and `theta_l` (y = eta).
// Models without covariates pass empty `x` and `gamma`,
// and models without covariates in the measurement model pass empty `kappa`.
// [[Rcpp::export(.SimSSMKernel)]]
Rcpp::List SimSSMKernel(const arma::uword n, const arma::uword time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& x,
                        const Rcpp::List& gamma, const Rcpp::List& kappa) {
  const bool latent = nu.size() == 0;
  const bool covariates = gamma.size() > 0;
  const bool covariates_y = kappa.size() > 0;
  const std::vector<arma::mat> mu0_s = SimSSMSets(mu0);
  const std::vector<arma::mat> sigma0_l_s = SimSSMSets(sigma0_l);
  const std::vector<arma::mat> alpha_s = SimSSMSets(alpha);
  const std::vector<arma::mat> beta_s = SimSSMSets(beta);
  const std::vector<arma::mat> psi_l_s = SimSSMSets(psi_l);
  const std::vector<arma::mat> nu_s = SimSSMSets(nu);
  const std::vector<arma::mat> lambda_s = SimSSMSets(lambda);
  const std::vector<arma::mat> theta_l_s = SimSSMSets(theta_l);
  const std::vector<arma::mat> x_s = SimSSMSets(x);
  const std::vector<arma::mat> gamma_s = SimSSMSets(gamma);
  const std::vector<arma::mat> kappa_s = SimSSMSets(kappa);
  const arma::uword p = mu0_s[0].n_elem;
  const arma::uword k = latent ? 0 : nu_s[0].n_elem;
  // arma::randn draws normal deviates in pairs, so a draw of odd length
  // discards its last deviate. Padding the state and the measurement block
  // to even length makes one draw per person reproduce the draws of earlier
  // versions, one call per block and time point (same seed, same data).
  const arma::uword p2 = p + p % 2;
  const arma::uword k2 = k + k % 2;
  const arma::span state(0, p - 1);
  const arma::vec time_vec = arma::linspace(0, (time - 1) * delta_t, time);
  Rcpp::List output(n);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::mat& beta_i = beta_s[i % beta_s.size()];
    const arma::mat z = arma::randn(p2 + k2, time);
    // everything but beta eta_{t - 1}, then the recursion
    arma::mat eta(p, time);
    eta.col(0) = mu0_s[i % mu0_s.size()] +
                 sigma0_l_s[i % sigma0_l_s.size()] * z(state, 0);
    if (time > 1) {
      const arma::span later(1, time - 1);
      eta.cols(later) = psi_l_s[i % psi_l_s.size()] * z(state, later);
      eta.cols(later).each_col() += arma::vec(alpha_s[i % alpha_s.size()]);
    }
    arma::mat x_i;
    if (covariates) {
      x_i = x_s[i % x_s.size()];
      eta += gamma_s[i % gamma_s.size()] * x_i;
    }
    // plain loops: small matrix-vector products without temporaries
    const double* b = beta_i.memptr();
    double* e = eta.memptr();
    for (arma::uword t = 1; t < time; ++t) {
      const double* prev = e + (t - 1) * p;
      double* curr = e + t * p;
      for (arma::uword c = 0; c < p; ++c) {
        const double prev_c = prev[c];
        const double* b_c = b + c * p;
        for (arma::uword r = 0; r < p; ++r) curr[r] += b_c[r] * prev_c;
      }
    }
    const arma::mat eta_t = eta.t();
    arma::mat y_t;
    if (latent) {
      y_t = eta_t;
    } else {
      y_t = eta_t * lambda_s[i % lambda_s.size()].t() +
            z.rows(p2, p2 + k - 1).t() * theta_l_s[i % theta_l_s.size()].t();
      y_t.each_row() += arma::rowvec(nu_s[i % nu_s.size()].t());
      if (covariates_y) {
        y_t += x_i.t() * kappa_s[i % kappa_s.size()].t();
      }
    }
    arma::vec id(time);
    id.fill(i + 1);
    if (covariates) {
      output[i] = Rcpp::List::create(
          Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
          Rcpp::Named("y") = y_t, Rcpp::Named("eta") = eta_t,
          Rcpp::Named("x") = x_i.t());
    } else {
      output[i] = Rcpp::List::create(
          Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
          Rcpp::Named("y") = y_t, Rcpp::Named("eta") = eta_t);
    }
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-solve-lya-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SolveLya)]]
arma::mat SolveLya(const arma::mat A, const arma::mat Q) {
  arma::mat X;
  arma::sylvester(X, A, A.t(), Q);
  return ((X + X.t()) / 2);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-solve-syl-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SolveSyl)]]
arma::mat SolveSyl(const arma::mat A, const arma::mat B, const arma::mat C) {
  arma::mat X;
  arma::sylvester(X, A, B, C);
  return X;
}
// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-spectral-abscissa.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Spectral Abscissa
//'
//' Returns the maximum real part of the eigenvalues of a square matrix.
//' For continuous-time stability (Hurwitz),
//' a matrix is stable if the spectral abscissa is strictly less than 0.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//'
//' @return Numeric value \eqn{\alpha(x) = \max \Re(\lambda_i(x))}.
//'
//' @examples
//' # Hurwitz-stable (spectral abscissa < 0):
//' x <- matrix(
//'   data = c(
//'     -0.5, -0.2,
//'      1.0, -0.3
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # < 0
//'
//' # Unstable (spectral abscissa > 0):
//' x <- matrix(
//'   data = c(
//'      0.10, 0.50,
//'     -0.40, 0.20
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # > 0
//'
//' @keywords simStateSpace stability linsde
//' @export
// [[Rcpp::export]]
double SpectralAbscissa(const arma::mat& x) {
  arma::cx_vec ev = arma::eig_gen(x);
  return arma::max(arma::real(ev));
}
// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-spectral-radius.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Spectral Radius
//'
//' Computes the spectral radius of a square matrix,
//' defined as the maximum modulus (absolute value) of its eigenvalues.
//' The spectral radius is often used to assess the stability of systems
//' such as vector autoregressive (VAR) models:
//' a system is considered stationary
//' if the spectral radius of its transition matrix is strictly less than 1.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//'
//' @return Numeric value representing the spectral radius of `x`.
//'
//' @examples
//' # Matrix with eigenvalues less than 1
//' x <- matrix(
//'   data = c(
//'     0.5, 0.3,
//'     0.2, 0.4
//'   ),
//'   nrow = 2
//' )
//' SpectralRadius(x)
//'
//' # Matrix with eigenvalues greater than 1
//' y <- matrix(
//'   data = c(
//'     1.2, 0.3,
//'     0.4, 0.9
//'   ),
//'   nrow = 2
//' )
//' SpectralRadius(y)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace stability ssm
//' @export
// [[Rcpp::export]]
double SpectralRadius(const arma::mat& x) {
  arma::cx_vec ev = arma::eig_gen(x);
  return arma::abs(ev).max();
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-cov-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMCovEta)]]
arma::mat SSMCovEta(const arma::mat& beta, const arma::mat& psi) {
  arma::vec vec_sigma = arma::solve(
      arma::eye(beta.n_rows * beta.n_rows, beta.n_rows * beta.n_rows) -
          arma::kron(beta, beta),
      arma::vectorise(psi));
  arma::mat X = arma::reshape(vec_sigma, beta.n_rows, beta.n_rows);
  return ((X + X.t()) / 2);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-cov-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMCovY)]]
arma::mat SSMCovY(const arma::mat& lambda, const arma::mat& theta,
                  const arma::mat& cov_eta) {
  arma::mat X = lambda * cov_eta * lambda.t() + theta;
  return ((X + X.t()) / 2);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-intercept-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMInterceptEta)]]
Rcpp::NumericVector SSMInterceptEta(const arma::mat& beta,
                                    const arma::vec& mean_eta) {
  arma::vec output = mean_eta - beta * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-intercept-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMInterceptY)]]
Rcpp::NumericVector SSMInterceptY(const arma::vec& mean_y,
                                  const arma::vec& mean_eta,
                                  const arma::mat& lambda) {
  arma::vec output = mean_y - lambda * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-mean-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMMeanEta)]]
Rcpp::NumericVector SSMMeanEta(const arma::mat& beta, const arma::vec& alpha) {
  arma::vec output =
      arma::solve(arma::eye(beta.n_rows, beta.n_rows) - beta, alpha);
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-mean-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMMeanY)]]
Rcpp::NumericVector SSMMeanY(const arma::vec& nu, const arma::mat& lambda,
                             const arma::vec& mean_eta) {
  arma::vec output = nu + lambda * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-phi-hurwitz.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Hurwitz Stability of a Drift Matrix
//'
//' Returns `TRUE` iff the drift matrix \eqn{\boldsymbol{\Phi}}
//' is Hurwitz-stable,
//' i.e., all eigenvalues have real parts strictly less than `-eps`.
//' Setting `eps = 0` enforces the usual strict condition
//' \eqn{\max \Re\{\lambda_i(\boldsymbol{\Phi})\} < 0}.
//' A small positive `eps` (e.g., `1e-12`) can be used
//' to guard against floating-point round-off.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param eps Nonnegative numeric tolerance (default `0.0`).
//'   The test checks \eqn{\Re(\lambda_i) < -\text{eps}} for all eigenvalues.
//'
//' @examples
//' # Unstable example (spectral abscissa >= 0):
//' phi <- matrix(
//'   data = c(
//'     0.10, -0.40,
//'     0.50, 0.20
//'   ),
//'   nrow = 2
//' )
//' TestPhiHurwitz(phi = phi) # FALSE
//'
//' # Stable example (all real parts < 0):
//' phi <- matrix(
//'   data = c(
//'     -0.50, -0.20,
//'      1.00, -0.30
//'   ),
//'   nrow = 2
//' )
//' TestPhiHurwitz(phi = phi) # TRUE
//' TestPhiHurwitz(phi = phi, eps = 1e-12) # also TRUE with tolerance
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestPhiHurwitz(const arma::mat& phi, const double eps = 0.0) {
  arma::cx_vec ev = arma::eig_gen(phi);
  return arma::all(arma::real(ev) < -eps);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-phi.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test the Drift Matrix
//'
//' Both have to be true for the function to return `TRUE`.
//'   - Test that the real part of all eigenvalues of \eqn{\boldsymbol{\Phi}}
//'     are less than zero.
//'   - Test that the diagonal values of \eqn{\boldsymbol{\Phi}}
//'     are between 0 to negative inifinity.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param margin Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//' @param auto_ubound Numeric scalar specifying the upper bound
//'   for the diagonal elements of \eqn{\boldsymbol{\Phi}}.
//'   Default is `0.0`, requiring all diagonal values to be \eqn{\leq 0}.
//'
//' @examples
//' phi <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' TestPhi(phi = phi)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestPhi(const arma::mat& phi, const double margin = 0.0,
             const double auto_ubound = 0.0) {
  arma::vec phi_diag = phi.diag(0);
  arma::cx_vec eigenvalues_phi = arma::eig_gen(phi);
  return arma::all(arma::real(eigenvalues_phi) < margin) &&
         arma::all(phi_diag <= auto_ubound);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-stability.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Stability
//'
//' The function computes the eigenvalues of the input matrix `x`.
//' It checks if the real part of all eigenvalues is negative.
//' If all eigenvalues have negative real parts,
//' the system is considered stable.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric matrix.
//' @param margin Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//'
//' @examples
//' x <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' TestStability(x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestStability(const arma::mat& x, const double margin = 0.0) {
  arma::cx_vec eigenvalues = arma::eig_gen(x);
  return arma::all(arma::real(eigenvalues) < margin);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-stationarity.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Stationarity
//'
//' The function computes the eigenvalues of the input matrix `x`.
//' It checks if all eigenvalues have moduli less than 1.
//' If all eigenvalues have moduli less than 1,
//' the system is considered stationary.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric matrix.
//' @param margin Numeric scalar specifying the stationarity threshold.
//'   Values less than 1 indicate stricter stationarity criteria.
//'
//' @examples
//' x <- matrix(
//'   data = c(0.5, 0.3, 0.2, 0.4),
//'   nrow = 2
//' )
//' TestStationarity(x)
//'
//' x <- matrix(
//'   data = c(0.9, -0.5, 0.8, 0.7),
//'   nrow = 2
//' )
//' TestStationarity(x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test ssm
//' @export
// [[Rcpp::export]]
bool TestStationarity(const arma::mat& x, const double margin = 1.0) {
  arma::cx_vec eigenvalues = arma::eig_gen(x);
  return arma::all(arma::abs(eigenvalues) < margin);
}
