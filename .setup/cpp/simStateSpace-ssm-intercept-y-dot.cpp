// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-intercept-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMInterceptY)]]
Rcpp::NumericVector SSMInterceptY(const arma::vec& mean_y, const arma::vec& mean_eta, const arma::mat& lambda) {
  arma::vec output = mean_y + lambda * mean_eta;
  return Rcpp::NumericVector(
    output.begin(),
    output.end()
  );
}
