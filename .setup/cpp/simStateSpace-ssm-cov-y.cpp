// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-cov-y.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Observed Variable Covariance Matrix for the
//' State Space Model
//'
//' This function calculates the observed variable covariance matrix
//' for the state space model
//' given by
//' \deqn{
//'   \mathrm{Cov} \left( \mathbf{y} \right)
//'   =
//'   \boldsymbol{\Lambda}
//'   \mathrm{Cov} \left( \boldsymbol{\eta} \right)
//'   \boldsymbol{\Lambda}^{\prime}
//'   +
//'   \boldsymbol{\Theta} .
//' }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param lambda Numeric matrix.
//'   The factor loadings matrix (\eqn{\boldsymbol{\Lambda}}).
//' @param theta Numeric matrix.
//'   The covariance matrix
//'   of the measurement error
//'   (\eqn{\boldsymbol{\Theta}}).
//' @param cov_eta Numeric matrix.
//'   State covariance matrix
//'   \eqn{\mathrm{Cov} \left( \boldsymbol{\eta} \right)}.
//'
//' @examples
//' beta <- 0.50 * diag(3)
//' psi <- 0.001 * diag(3)
//' lambda <- diag(3)
//' theta <- 0.02 * diag(3)
//' cov_eta <- SSMCovEta(beta = beta, psi = psi)
//' SSMCovY(
//'   lambda = lambda,
//'   theta = theta,
//'   cov_eta = cov_eta
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
arma::mat SSMCovY(const arma::mat& lambda, const arma::mat& theta,
                  const arma::mat& cov_eta) {
  return lambda * cov_eta * lambda.t() + theta;
}
