// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-linear-sde-2-ssm.cpp
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
//' @details The state space parameters
//'   as a function of the linear stochastic differential equation model
//'   parameters
//'   are given by
//'   \deqn{
//'       \boldsymbol{\beta}
//'       =
//'       \exp{
//'         \left(
//'           \mathbf{A}
//'           \Delta_{t}
//'         \right)
//'       }
//'   }
//'
//'   \deqn{
//'       \boldsymbol{\alpha}
//'       =
//'       \mathbf{A}^{-1}
//'       \left(
//'         \boldsymbol{\beta} - \mathbf{I}_{p}
//'       \right)
//'       \mathbf{b}
//'   }
//'
//'   \deqn{
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Psi}
//'       \right)
//'       =
//'       \left\{
//'         \left[
//'           \left(
//'             \mathbf{A} \otimes \mathbf{I}_{p}
//'           \right)
//'           +
//'           \left(
//'             \mathbf{I}_{p} \otimes \mathbf{A}
//'           \right)
//'         \right]
//'         \left[
//'           \exp
//'           \left(
//'             \left[
//'               \left(
//'                 \mathbf{A} \otimes \mathbf{I}_{p}
//'               \right)
//'               +
//'               \left(
//'                 \mathbf{I}_{p} \otimes \mathbf{A}
//'               \right)
//'             \right]
//'             \Delta_{t}
//'         \right)
//'         -
//'         \mathbf{I}_{p \times p}
//'       \right]
//'       \mathrm{vec}
//'       \left(
//'         \mathbf{Q}
//'       \right)
//'     \right\}
//'   }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @inheritParams SimSSMLinSDE
//'
//' @return Returns a list of state space parameters:
//'   - `alpha`: Numeric vector.
//'     Vector of intercepts for the dynamic model
//'     (\eqn{\boldsymbol{\alpha}}).
//'   - `beta`: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     at time `t_k - 1` to those at time `t_k`
//'     (\eqn{\boldsymbol{\beta}}).
//'   - `psi`: Numeric matrix.
//'     The process noise covariance matrix
//'     (\eqn{\boldsymbol{\Psi}}).
//'
//' @examples
//' p <- k <- 2
//' a <- matrix(
//'   data = c(
//'    -0.10,
//'    0.05,
//'    0.05,
//'    -0.10
//'  ),
//'  nrow = p
//' )
//' b <- c(0.317, 0.230)
//' q <- matrix(
//'   data = c(2.79, 0.06, 0.06, 3.27),
//'   nrow = p
//' )
//' delta_t <- 0.10
//'
//' LinSDE2SSM(
//'   b = b,
//'   a = a,
//'   q = q,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim linsde
//' @export
// [[Rcpp::export]]
Rcpp::List LinSDE2SSM(const arma::vec& b, const arma::mat& a,
                      const arma::mat& q, const double delta_t) {
  // Step 1: Determine indices
  int num_latent_vars = b.n_elem;

  // Step 2: Get state space parameters
  arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
  arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                     num_latent_vars * num_latent_vars);
  // 2.1 beta
  arma::mat beta = arma::expmat(a * delta_t);
  // 2.2 alpha
  arma::vec alpha = arma::inv(a) * (beta - I) * b;
  // 2.3 psi
  arma::mat a_hashtag = arma::kron(a, I) + arma::kron(I, a);
  arma::vec q_vec = arma::vectorise(q);
  arma::vec psi_vec =
      arma::inv(a_hashtag) * (arma::expmat(a_hashtag * delta_t) - J) * q_vec;
  arma::mat psi =
      arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

  // Step 3: Return state space parameters in a list
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi") = psi);
}
