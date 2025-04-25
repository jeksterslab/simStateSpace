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
