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
//' @param margin Double in \eqn{(0, 1)}. Target upper bound for the spectral
//'   radius (default = 0.98).
//' @param tol Small positive double added to the denominator in the scaling
//'   factor to avoid division by zero (default = 1e-12).
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
