// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-solve-syl-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SolveSyl)]]
arma::mat SolveSyl(arma::mat A, arma::mat B, arma::mat C) {
  arma::mat X;
  arma::syl(X, A, B, C);
  return X;
}
