// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-solve-syl-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SolveSyl)]]
arma::mat SolveSyl(const arma::mat A, const arma::mat B, const arma::mat C) {
  arma::mat X;
  arma::syl(X, A, B, C);
  return X;
}
