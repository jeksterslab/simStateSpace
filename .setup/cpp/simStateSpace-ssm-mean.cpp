// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-mean.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' State Mean Vector for the
//' State Space Model
//'
//' This function calculates the state mean vector
//' for the state space model
//' given by
//' \deqn{
//'   \mathrm{Mean} \left( \boldsymbol{\eta} \right)
//'   =
//'   \left(
//'     \mathbf{I} - \boldsymbol{\beta}
//'   \right)^{-1}
//'   \boldsymbol{\alpha} .
//' }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param beta Numeric matrix.
//'   The transition matrix (\eqn{\boldsymbol{\beta}}).
//' @param alpha Numeric vector.
//'   Vector of constant values for the dynamic model
//'   (\eqn{\boldsymbol{\alpha}}).
//' @examples
//' beta <- 0.50 * diag(3)
//' alpha <- rep(x = 0.001, times = 3)
//' SSMMean(beta = beta, alpha = alpha)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector SSMMean(const arma::mat& beta, const arma::vec& alpha) {
  arma::vec output =
      arma::solve(arma::eye(beta.n_rows, beta.n_rows) - beta, alpha);
  return Rcpp::NumericVector(output.begin(), output.end());
}
