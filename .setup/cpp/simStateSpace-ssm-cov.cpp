// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-cov.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' State Covariance Matrix for the
//' State Space Model
//'
//' This function calculates the state covariance matrix
//' for the state space model
//' given by
//' \deqn{
//'   \mathrm{vec}
//'   \left(
//'     \mathrm{Cov} \left( \boldsymbol{\eta} \right)
//'   \right)
//'   =
//'   \left(
//'     \mathbf{I} - \boldsymbol{\beta} \otimes \boldsymbol{\beta}
//'   \right)^{-1}
//'   \mathrm{vec} \left( \boldsymbol{\Psi} \right) .
//' }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param beta Numeric matrix.
//'   The transition matrix (\eqn{\boldsymbol{\beta}}).
//' @param psi Numeric matrix.
//'   The covariance matrix
//'   of the process noise
//'   (\eqn{\boldsymbol{\Psi}}).
//'
//' @examples
//' beta <- 0.50 * diag(3)
//' psi <- 0.001 * diag(3)
//' SSMCov(beta = beta, psi = psi)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
arma::mat SSMCov(const arma::mat& beta, const arma::mat& psi) {
  arma::vec vec_sigma = arma::solve(
      arma::eye(beta.n_rows * beta.n_rows, beta.n_rows * beta.n_rows) -
          arma::kron(beta, beta),
      arma::vectorise(psi));
  return arma::reshape(vec_sigma, beta.n_rows, beta.n_rows);
}
