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
//'       \boldsymbol{\gamma}
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
//'   for individual \eqn{i} and time {t}.
//'   The state space parameters
//'   as a function of the linear stochastic differential equation model
//'   parameters
//'   are given by
//'   \deqn{
//'       \boldsymbol{\beta}
//'       =
//'       \exp{
//'         \left(
//'           \boldsymbol{\Phi}
//'           \Delta_{t}
//'         \right)
//'       }
//'   }
//'
//'   \deqn{
//'       \boldsymbol{\alpha}
//'       =
//'       \boldsymbol{\Phi}^{-1}
//'       \left(
//'         \boldsymbol{\beta} - \mathbf{I}_{p}
//'       \right)
//'       \boldsymbol{\gamma}
//'   }
//'
//'   \deqn{
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Psi}
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
//'           \Delta_{t}
//'         \right)
//'         -
//'         \mathbf{I}_{p \times p}
//'       \right]
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Sigma}
//'       \right)
//'   }
//'   where \eqn{p} is the number of latent variables and
//'   \eqn{\Delta_{t}} is the time interval.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param gamma Numeric vector.
//'   An unobserved term that is constant over time
//'   (\eqn{\boldsymbol{\gamma}}).
//' @param phi Numeric matrix.
//'   The drift matrix
//'   which represents the rate of change of the solution
//'   in the absence of any random fluctuations
//'   (\eqn{\boldsymbol{\Phi}}).
//' @param sigma_l Numeric matrix.
//'   Cholesky factorization (`t(chol(sigma))`)
//'   of the covariance matrix of volatility
//'   or randomness in the process
//'   \eqn{\boldsymbol{\Sigma}}.
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
//' gamma <- c(0.317, 0.230)
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
//'   gamma = gamma,
//'   phi = phi,
//'   sigma_l = sigma_l,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim linsde
//' @export
// [[Rcpp::export]]
Rcpp::List LinSDE2SSM(const arma::vec& gamma, const arma::mat& phi, const arma::mat& sigma_l, const double delta_t) {
  int p = gamma.n_elem;
  arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat J = arma::eye<arma::mat>(p * p, p * p);
  arma::mat phi_hashtag = arma::kron(phi, I) + arma::kron(I, phi);
  arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
  arma::vec psi_vec = arma::inv(phi_hashtag) * (arma::expmat(phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi_l = arma::chol(arma::reshape(psi_vec, p, p), "lower");
  arma::mat beta = arma::expmat(phi * delta_t);
  arma::vec alpha = arma::inv(phi) * (beta - I) * gamma;
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha, Rcpp::Named("beta") = beta, Rcpp::Named("psi_l") = psi_l);
}
