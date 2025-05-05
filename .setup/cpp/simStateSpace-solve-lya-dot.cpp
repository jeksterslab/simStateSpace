// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-solve-lya-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SolveLya)]]
arma::mat SolveLya(const arma::mat A, const arma::mat Q) {
  arma::mat X;
  arma::syl(X, A, A.t(), Q);
  return ((X + X.t()) / 2);
}
