// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Transition Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random transition matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated transition matrices are stationary
//' using [TestStationarity()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param beta Numeric matrix.
//'   The transition matrix (\eqn{\boldsymbol{\beta}}).
//' @param vcov_beta_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_beta_vec))`)
//'   of the sampling variance-covariance matrix
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\beta} \right)}.
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
Rcpp::List SimBetaN(const arma::uword& n, const arma::mat& beta,
                    const arma::mat& vcov_beta_vec_l) {
  Rcpp::List output(n);
  int p = beta.n_rows;
  int q = p * p;
  arma::vec beta_vec = arma::vectorise(beta);
  for (arma::uword i = 0; i < n; i++) {
    bool run = true;
    while (run) {
      arma::vec beta_vec_i = beta_vec + (vcov_beta_vec_l * arma::randn(q));
      arma::mat beta_i = arma::reshape(beta_vec_i, p, p);
      if (TestStationarity(beta_i)) {
        run = false;
      }
      if (!run) {
        output[i] = beta_i;
      }
    }
  }
  return output;
}
