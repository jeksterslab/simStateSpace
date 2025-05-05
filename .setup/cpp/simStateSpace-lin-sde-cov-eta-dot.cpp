// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-cov-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDECovEta)]]
arma::mat LinSDECovEta(const arma::mat phi, const arma::mat sigma) {
  arma::mat X;
  arma::syl(X, phi, phi.t(), sigma);
  return ((X + X.t()) / 2);
}
