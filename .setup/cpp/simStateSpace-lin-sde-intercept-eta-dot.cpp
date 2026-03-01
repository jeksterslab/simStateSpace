// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-intercept-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDEInterceptEta)]]
Rcpp::NumericVector LinSDEInterceptEta(const arma::mat& phi,
                                       const arma::vec& mean_eta) {
  arma::vec output = -phi * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
