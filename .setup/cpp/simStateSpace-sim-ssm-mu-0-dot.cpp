// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-mu-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.Mu0)]]
Rcpp::List Mu0(const arma::vec& alpha, const arma::mat& beta, const arma::vec& nu) {
  int p = beta.n_rows;
  arma::vec mu_eta = arma::inv(arma::eye(p, p) - beta) * alpha;
  arma::vec mu_y = nu;
  return Rcpp::List::create(Rcpp::Named("mu_y") = mu_y, Rcpp::Named("mu_eta") = mu_eta);
}
