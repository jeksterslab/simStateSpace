// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-mean-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDEMeanEta)]]
Rcpp::NumericVector LinSDEMeanEta(const arma::mat phi, const arma::vec iota) {
  arma::vec output = arma::solve(-phi, iota);
  return Rcpp::NumericVector(output.begin(), output.end());
}
