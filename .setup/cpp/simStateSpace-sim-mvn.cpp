// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-mvn.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimMVN)]]
Rcpp::List SimMVN(const arma::uword& n, const arma::vec& mu,
                  const arma::mat& sigma_l) {
  Rcpp::List output(n);
  arma::vec mu_i(mu.n_rows, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    mu_i = mu + (sigma_l * arma::randn(mu.n_rows));
    output[i] = Rcpp::NumericVector(mu_i.begin(), mu_i.end());
  }
  return output;
}
