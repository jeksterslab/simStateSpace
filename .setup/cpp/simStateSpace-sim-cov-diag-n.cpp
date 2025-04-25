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
