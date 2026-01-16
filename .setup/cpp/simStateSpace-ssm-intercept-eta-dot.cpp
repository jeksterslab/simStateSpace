// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-intercept-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMInterceptEta)]]
Rcpp::NumericVector SSMInterceptEta(const arma::mat& beta,
                                    const arma::vec& mean_eta) {
  arma::vec output = mean_eta - beta * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
