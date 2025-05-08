// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-cov-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMCovY)]]
arma::mat SSMCovY(const arma::mat& lambda, const arma::mat& theta,
                  const arma::mat& cov_eta) {
  arma::mat X = lambda * cov_eta * lambda.t() + theta;
  return ((X + X.t()) / 2);
}
