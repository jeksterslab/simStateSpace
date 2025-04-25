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
    output[i] = alpha_i;
  }
  return output;
}
