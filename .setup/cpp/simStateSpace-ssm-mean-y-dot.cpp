// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-mean-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMMeanY)]]
Rcpp::NumericVector SSMMeanY(const arma::vec& nu, const arma::mat& lambda,
                             const arma::vec& mean_eta) {
  arma::vec output = nu + lambda * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
