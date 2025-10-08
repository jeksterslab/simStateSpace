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
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-2-ssm.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Convert Parameters from the Linear Stochastic Differential Equation Model
//' to State Space Model Parameterization
//'
//' This function converts parameters from
//' the linear stochastic differential equation model
//' to state space model parameterization.
//'
//' @details Let the linear stochastic equation model be given by
//'   \deqn{
//'     \mathrm{d}
//'     \boldsymbol{\eta}_{i, t}
//'     =
//'     \left(
//'       \boldsymbol{\iota}
//'       +
//'       \boldsymbol{\Phi}
//'       \boldsymbol{\eta}_{i, t}
//'     \right)
//'     \mathrm{d} t
//'     +
//'     \boldsymbol{\Sigma}^{\frac{1}{2}}
//'     \mathrm{d}
//'     \mathbf{W}_{i, t}
//'   }
//'   for individual \eqn{i} and time \eqn{t}.
//'   The discrete-time state space model
//'   given below
//'   represents the discrete-time solution
//'   for the linear stochastic differential equation.
//'   \deqn{
//'     \boldsymbol{\eta}_{i, t_{{l_{i}}}}
//'     =
//'     \boldsymbol{\alpha}_{\Delta t_{{l_{i}}}}
//'     +
//'     \boldsymbol{\beta}_{\Delta t_{{l_{i}}}}
//'     \boldsymbol{\eta}_{i, t_{l_{i} - 1}}
//'     +
//'     \boldsymbol{\zeta}_{i, t_{{l_{i}}}},
//'     \quad
//'     \mathrm{with}
//'     \quad
//'     \boldsymbol{\zeta}_{i, t_{{l_{i}}}}
//'     \sim
//'     \mathcal{N}
//'     \left(
//'     \mathbf{0},
//'     \boldsymbol{\Psi}_{\Delta t_{{l_{i}}}}
//'     \right)
//'   }
//'   with
//'   \deqn{
//'       \boldsymbol{\beta}_{\Delta t_{{l_{i}}}}
//'       =
//'       \exp{
//'         \left(
//'           \Delta t
//'           \boldsymbol{\Phi}
//'         \right)
//'       },
//'   }
//'
//'   \deqn{
//'       \boldsymbol{\alpha}_{\Delta t_{{l_{i}}}}
//'       =
//'       \boldsymbol{\Phi}^{-1}
//'       \left(
//'         \boldsymbol{\beta} - \mathbf{I}_{p}
//'       \right)
//'       \boldsymbol{\iota}, \quad \mathrm{and}
//'   }
//'
//'   \deqn{
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Psi}_{\Delta t_{{l_{i}}}}
//'       \right)
//'       =
//'       \left[
//'         \left(
//'           \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'         \right)
//'         +
//'         \left(
//'           \mathbf{I}_{p} \otimes \boldsymbol{\Phi}
//'         \right)
//'       \right]
//'       \left[
//'         \exp
//'         \left(
//'           \left[
//'             \left(
//'               \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'             \right)
//'             +
//'             \left(
//'               \mathbf{I}_{p} \otimes \boldsymbol{\Phi}
//'             \right)
//'           \right]
//'           \Delta t
//'         \right)
//'         -
//'         \mathbf{I}_{p \times p}
//'       \right]
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Sigma}
//'       \right)
//'   }
//'   where \eqn{t} denotes continuous-time processes
//'   that can be defined by any arbitrary time point,
//'   \eqn{t_{l_{i}}} the \eqn{l^\mathrm{th}}
//'   observed measurement occassion for individual \eqn{i},
//'   \eqn{p} the number of latent variables and
//'   \eqn{\Delta t} the time interval.
//'
//' @references
//'   Harvey, A. C. (1990).
//'   Forecasting, structural time series models and the Kalman filter.
//'   Cambridge University Press.
//'   \doi{10.1017/cbo9781107049994}
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param iota Numeric vector.
//'   An unobserved term that is constant over time
//'   (\eqn{\boldsymbol{\iota}}).
//' @param phi Numeric matrix.
//'   The drift matrix
//'   which represents the rate of change of the solution
//'   in the absence of any random fluctuations
//'   (\eqn{\boldsymbol{\Phi}}).
//' @param sigma_l Numeric matrix.
//'   Cholesky factorization (`t(chol(sigma))`)
//'   of the covariance matrix of volatility
//'   or randomness in the process
//'   (\eqn{\boldsymbol{\Sigma}}).
//' @param delta_t Numeric.
//'   Time interval
//'   (\eqn{\Delta_t}).
//'
//' @return Returns a list of state space parameters:
//'   - `alpha`: Numeric vector.
//'     Vector of constant values for the dynamic model
//'     (\eqn{\boldsymbol{\alpha}}).
//'   - `beta`: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     from the previous time point to the current time point.
//'     (\eqn{\boldsymbol{\beta}}).
//'   - `psi_l`: Numeric matrix.
//'     Cholesky factorization (`t(chol(psi))`)
//'     of the process noise covariance matrix
//'     \eqn{\boldsymbol{\Psi}}.
//'
//' @examples
//' p <- 2
//' iota <- c(0.317, 0.230)
//' phi <- matrix(
//'   data = c(
//'    -0.10,
//'    0.05,
//'    0.05,
//'    -0.10
//'  ),
//'  nrow = p
//' )
//' sigma <- matrix(
//'   data = c(
//'     2.79,
//'     0.06,
//'     0.06,
//'     3.27
//'   ),
//'   nrow = p
//' )
//' sigma_l <- t(chol(sigma))
//' delta_t <- 0.10
//'
//' LinSDE2SSM(
//'   iota = iota,
//'   phi = phi,
//'   sigma_l = sigma_l,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace transformation linsde
//' @export
// [[Rcpp::export]]
Rcpp::List LinSDE2SSM(const arma::vec& iota, const arma::mat& phi,
                      const arma::mat& sigma_l, const double delta_t) {
  arma::mat I = arma::eye(iota.n_elem, iota.n_elem);
  arma::mat beta = arma::expmat(phi * delta_t);
  arma::vec alpha(iota.n_elem, arma::fill::none);
  arma::mat psi_l = arma::mat(iota.n_elem, iota.n_elem, arma::fill::none);
  // alpha = iota.is_zero() ? iota : arma::inv(phi) * (beta - I) * iota;
  alpha = iota.is_zero() ? iota : arma::solve(phi, (beta - I) * iota);
  if (sigma_l.is_zero()) {
    psi_l = sigma_l;
  } else {
    arma::mat J =
        arma::eye(iota.n_elem * iota.n_elem, iota.n_elem * iota.n_elem);
    arma::mat phi_hashtag = arma::kron(phi, I) + arma::kron(I, phi);
    arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
    // arma::vec psi_vec = arma::inv(phi_hashtag) * (arma::expmat(phi_hashtag *
    // delta_t) - J) * sigma_vec;
    arma::vec psi_vec = arma::solve(
        phi_hashtag, (arma::expmat(phi_hashtag * delta_t) - J) * sigma_vec);
    psi_l =
        arma::chol(arma::reshape(psi_vec, iota.n_elem, iota.n_elem), "lower");
  }
  // output
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi_l") = psi_l);
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
// edit .setup/cpp/simStateSpace-sim-alpha-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Intercept Vectors
//' in a Discrete-Time Vector Autoregressive Model
//' from the Multivariate Normal Distribution
//'
//' This function simulates random intercept vectors
//' in a discrete-time vector autoregressive model
//' from the multivariate normal distribution.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param alpha Numeric vector.
//'   Intercept (\eqn{\boldsymbol{\alpha}}).
//' @param vcov_alpha_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_alpha))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\boldsymbol{\alpha}}.
//' @return Returns a list of random intercept vectors.
//'
//' @examples
//' n <- 10
//' alpha <- c(0, 0, 0)
//' vcov_alpha_l <- t(chol(0.001 * diag(3)))
//' SimAlphaN(n = n, alpha = alpha, vcov_alpha_l = vcov_alpha_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimAlphaN(const arma::uword& n, const arma::vec& alpha,
                     const arma::mat& vcov_alpha_l) {
  Rcpp::List output(n);
  arma::vec alpha_i(alpha.n_rows, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    alpha_i = alpha + (vcov_alpha_l * arma::randn(alpha.n_rows));
    output[i] = Rcpp::NumericVector(alpha_i.begin(), alpha_i.end());
  }
  return output;
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
// edit .setup/cpp/simStateSpace-sim-beta-n.cpp
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

  // if (vcov_beta_vec_l.n_rows != p || vcov_beta_vec_l.n_cols != p) {
  //   Rcpp::stop("vcov_beta_vec_l must be p x p with p = nrow(beta) *
  //   ncol(beta).");
  // }

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
// edit .setup/cpp/simStateSpace-sim-iota-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Intercept Vectors
//' in a Continuous-Time Vector Autoregressive Model
//' from the Multivariate Normal Distribution
//'
//' This function simulates random intercept vectors
//' in a continuous-time vector autoregressive model
//' from the multivariate normal distribution.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param iota Numeric vector.
//'   Intercept (\eqn{\boldsymbol{\iota}}).
//' @param vcov_iota_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_iota))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\boldsymbol{\iota}}.
//' @return Returns a list of random intercept vectors.
//'
//' @examples
//' n <- 10
//' iota <- c(0, 0, 0)
//' vcov_iota_l <- t(chol(0.001 * diag(3)))
//' SimIotaN(n = n, iota = iota, vcov_iota_l = vcov_iota_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimIotaN(const arma::uword& n, const arma::vec& iota,
                    const arma::mat& vcov_iota_l) {
  Rcpp::List output(n);
  arma::vec iota_i(iota.n_rows, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    iota_i = iota + (vcov_iota_l * arma::randn(iota.n_rows));
    output[i] = Rcpp::NumericVector(iota_i.begin(), iota_i.end());
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-nu-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Intercept Vectors
//' in a Discrete-Time Vector Autoregressive Model
//' from the Multivariate Normal Distribution
//'
//' This function simulates random intercept vectors
//' in a discrete-time vector autoregressive model
//' from the multivariate normal distribution.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param nu Numeric vector.
//'   Intercept (\eqn{\boldsymbol{\nu}}).
//' @param vcov_nu_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_nu))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\boldsymbol{\nu}}.
//' @return Returns a list of random intercept vectors.
//'
//' @examples
//' n <- 10
//' nu <- c(0, 0, 0)
//' vcov_nu_l <- t(chol(0.001 * diag(3)))
//' SimNuN(n = n, nu = nu, vcov_nu_l = vcov_nu_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimNuN(const arma::uword& n, const arma::vec& nu,
                  const arma::mat& vcov_nu_l) {
  Rcpp::List output(n);
  arma::vec nu_i(nu.n_rows, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    nu_i = nu + (vcov_nu_l * arma::randn(nu.n_rows));
    output[i] = Rcpp::NumericVector(nu_i.begin(), nu_i.end());
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

  // if (vcov_phi_vec_l.n_rows != p || vcov_phi_vec_l.n_cols != p) {
  //   Rcpp::stop("vcov_phi_vec_l must be p x p with p = nrow(phi) *
  //   ncol(phi).");
  // }

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
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed0)]]
Rcpp::List SimSSMFixed0(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(mu0.n_elem));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) =
          alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(mu0.n_elem));
      y.col(t) =
          nu + (lambda * eta.col(t)) + (theta_l * arma::randn(nu.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed1)]]
Rcpp::List SimSSMFixed1(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const Rcpp::List& x,
                        const arma::mat& gamma) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) =
        mu0 + (sigma0_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(t));
      y.col(t) =
          nu + (lambda * eta.col(t)) + (theta_l * arma::randn(nu.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed2)]]
Rcpp::List SimSSMFixed2(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const Rcpp::List& x,
                        const arma::mat& gamma, const arma::mat& kappa) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) =
        mu0 + (sigma0_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem)) +
               (kappa * x_i.col(0));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(nu.n_elem)) + (kappa * x_i.col(t));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-i-vary-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMIVary0)]]
Rcpp::List SimSSMIVary0(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = Rcpp::as<arma::vec>(alpha[i]);
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-i-vary-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMIVary1)]]
Rcpp::List SimSSMIVary1(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& x,
                        const Rcpp::List& gamma) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = Rcpp::as<arma::vec>(alpha[i]);
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-i-vary-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMIVary2)]]
Rcpp::List SimSSMIVary2(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& x,
                        const Rcpp::List& gamma, const Rcpp::List& kappa) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = Rcpp::as<arma::vec>(alpha[i]);
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];
    arma::mat kappa_i = kappa[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) +
               (theta_l_i * arma::randn(nu_i.n_elem)) + (kappa_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem)) +
                 (kappa_i * x_i.col(t));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-fixed-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatFixed0)]]
Rcpp::List SimSSMLatFixed0(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const arma::vec& mu0,
                           const arma::mat& sigma0_l, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_l) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(mu0.n_elem));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) =
          alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(mu0.n_elem));
      y.col(t) = eta.col(t);
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-fixed-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatFixed1)]]
Rcpp::List SimSSMLatFixed1(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const arma::vec& mu0,
                           const arma::mat& sigma0_l, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_l,
                           const Rcpp::List& x, const arma::mat& gamma) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Generate initial condition
    eta.col(0) =
        mu0 + (sigma0_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(t));
      y.col(t) = eta.col(t);
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-i-vary-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatIVary0)]]
Rcpp::List SimSSMLatIVary0(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = Rcpp::as<arma::vec>(alpha[i]);
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = eta.col(t);
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-i-vary-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatIVary1)]]
Rcpp::List SimSSMLatIVary1(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l,
                           const Rcpp::List& x, const Rcpp::List& gamma) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = Rcpp::as<arma::vec>(alpha[i]);
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = eta.col(t);
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary0)]]
Rcpp::List SimSSMLinSDEIVary0(const arma::uword& n, const arma::uword& time,
                              const double delta_t, const Rcpp::List& mu0,
                              const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l,
                              const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = Rcpp::as<arma::vec>(iota[i]);
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lin-sde-i-vary-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary1)]]
Rcpp::List SimSSMLinSDEIVary1(const arma::uword& n, const arma::uword& time,
                              const double delta_t, const Rcpp::List& mu0,
                              const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l, const Rcpp::List& x,
                              const Rcpp::List& gamma, const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = Rcpp::as<arma::vec>(iota[i]);
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lin-sde-i-vary-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary2)]]
Rcpp::List SimSSMLinSDEIVary2(
    const arma::uword& n, const arma::uword& time, const double delta_t,
    const Rcpp::List& mu0, const Rcpp::List& sigma0_l, const Rcpp::List& iota,
    const Rcpp::List& phi, const Rcpp::List& sigma_l, const Rcpp::List& nu,
    const Rcpp::List& lambda, const Rcpp::List& theta_l, const Rcpp::List& x,
    const Rcpp::List& gamma, const Rcpp::List& kappa, const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = Rcpp::as<arma::vec>(iota[i]);
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];
    arma::mat kappa_i = kappa[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) +
               (theta_l_i * arma::randn(nu_i.n_elem)) + (kappa_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem)) +
                 (kappa_i * x_i.col(t));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
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
