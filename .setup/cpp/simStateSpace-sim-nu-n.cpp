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
