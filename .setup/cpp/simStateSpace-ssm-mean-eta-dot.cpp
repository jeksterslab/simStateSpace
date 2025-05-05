// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-mean-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMMeanEta)]]
Rcpp::NumericVector SSMMeanEta(const arma::mat& beta, const arma::vec& alpha) {
  arma::vec output =
      arma::solve(arma::eye(beta.n_rows, beta.n_rows) - beta, alpha);
  return Rcpp::NumericVector(output.begin(), output.end());
}
