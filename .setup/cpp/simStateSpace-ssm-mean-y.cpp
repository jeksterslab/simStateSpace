// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-mean-y.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Observed Variable Mean Vector for the
//' State Space Model
//'
//' This function calculates the observed variable mean vector
//' for the state space model
//' given by
//' \deqn{
//'   \mathrm{Mean} \left( \mathbf{y} \right)
//'   =
//'   \boldsymbol{\nu}
//'   +
//'   \boldsymbol{\Lambda}
//'   \mathrm{Mean} \left( \boldsymbol{\eta} \right) .
//' }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param nu Numeric vector.
//'   Vector of constant values for the measurement model
//'   (\eqn{\boldsymbol{\nu}}).
//' @param lambda Numeric matrix.
//'   The factor loadings matrix (\eqn{\boldsymbol{\Lambda}}).
//' @param mean_eta Numeric vector.
//'   State mean vector
//'   \eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}.
//'
//' @examples
//' beta <- 0.50 * diag(3)
//' alpha <- rep(x = 0.001, times = 3)
//' nu <- rep(x = 0.03, times = 3)
//' lambda <- diag(3)
//' mean_eta <- SSMMeanEta(beta = beta, alpha = alpha)
//' SSMMeanY(
//'   nu = nu,
//'   lambda = lambda,
//'   mean_eta = mean_eta
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector SSMMeanY(const arma::vec& nu, const arma::mat& lambda,
                             const arma::vec& mean_eta) {
  arma::vec output = nu + lambda * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
