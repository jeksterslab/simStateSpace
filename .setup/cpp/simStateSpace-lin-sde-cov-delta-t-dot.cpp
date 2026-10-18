// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-cov-delta-t-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.LinSDECovDeltaT)]]
arma::mat LinSDECovDeltaT(const arma::mat& phi, const arma::mat& sigma,
                          const double delta_t) {
  return AccumulatedCovariance(phi, sigma, delta_t);
}
