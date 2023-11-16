// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ou-2-ssm.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Convert Parameters from the Ornstein–Uhlenbeck Model
//' to State Space Model Parameterization
//'
//' This function converts parameters from the Ornstein–Uhlenbeck model
//' to state space model parameterization.
//' See details for more information.
//'
//' @details The state space parameters
//'   as a function of the  Ornstein–Uhlenbeck model parameters
//'   are given by
//'   \deqn{
//'       \boldsymbol{\beta}
//'       =
//'       \exp{
//'         \left(
//'           - \boldsymbol{\Phi}
//'           \Delta_{t}
//'         \right)
//'       }
//'   }
//'
//'   \deqn{
//'       \boldsymbol{\alpha}
//'       =
//'       - \boldsymbol{\Phi}^{-1}
//'       \left(
//'         \boldsymbol{\beta} - \mathbf{I}_{p}
//'       \right)
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
//'             - \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'           \right)
//'           +
//'           \left(
//'             \mathbf{I}_{p} \otimes - \boldsymbol{\Phi}
//'           \right)
//'         \right]
//'         \left[
//'           \exp
//'           \left(
//'             \left[
//'               \left(
//'                 - \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'               \right)
//'               +
//'               \left(
//'                 \mathbf{I}_{p} \otimes - \boldsymbol{\Phi}
//'               \right)
//'             \right]
//'             \Delta_{t}
//'         \right)
//'         -
//'         \mathbf{I}_{p \times p}
//'       \right]
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Sigma}
//'       \right)
//'     \right\}
//'   }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @inheritParams SimSSMOU
//'
//' @return Returns a list of state space parameters:
//'   - alpha: Numeric vector.
//'     Vector of intercepts for the dynamic model
//'     (\eqn{\boldsymbol{\alpha}}).
//'   - beta: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     at time `t - 1` to those at time `t`
//'     (\eqn{\boldsymbol{\beta}}).
//'   - psi: Numeric matrix.
//'     The process noise covariance matrix
//'     (\eqn{\boldsymbol{\Psi}}).
//'
//' @examples
//' p <- k <- 2
//' mu <- c(5.76, 5.18)
//' phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
//' sigma_sqrt <- chol(
//'   matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
//' )
//' delta_t <- 0.10
//'
//' OU2SSM(
//'   mu = mu,
//'   phi = phi,
//'   sigma_sqrt = sigma_sqrt,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim
//' @export
// [[Rcpp::export]]
Rcpp::List OU2SSM(const arma::vec& mu, const arma::mat& phi,
                  const arma::mat& sigma_sqrt, const double delta_t) {
  // Step 1: Determine indices
  int num_latent_vars = mu.n_elem;

  // Step 2: Get state space parameters
  arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
  arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                     num_latent_vars * num_latent_vars);
  arma::mat neg_phi = -1 * phi;
  // 2.1 beta
  arma::mat beta = arma::expmat(neg_phi * delta_t);  // A(Delta t)
  // 2.2 alpha
  arma::vec alpha = arma::inv(neg_phi) * (beta - I) * phi * mu;  // b(Delta t)
  // 2.3 psi
  arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
  arma::vec sigma_vec = arma::vectorise(sigma_sqrt * sigma_sqrt.t());
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi = arma::reshape(psi_vec, num_latent_vars, num_latent_vars);

  // Step 3: Return state space parameters in a list
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi") = psi);
}
