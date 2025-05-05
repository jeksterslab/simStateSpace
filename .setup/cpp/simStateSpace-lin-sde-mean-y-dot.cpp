// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-mean-y-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDEMeanY)]]
Rcpp::NumericVector LinSDEMeanY(const arma::vec& nu, const arma::mat& lambda,
                                const arma::vec& mean_eta) {
  arma::vec output = nu + lambda * mean_eta;
  return Rcpp::NumericVector(output.begin(), output.end());
}
