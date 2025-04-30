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
