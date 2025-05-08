// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ssm-cov-eta-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSMCovEta)]]
arma::mat SSMCovEta(const arma::mat& beta, const arma::mat& psi) {
  arma::vec vec_sigma = arma::solve(
      arma::eye(beta.n_rows * beta.n_rows, beta.n_rows * beta.n_rows) -
          arma::kron(beta, beta),
      arma::vectorise(psi));
  arma::mat X = arma::reshape(vec_sigma, beta.n_rows, beta.n_rows);
  return ((X + X.t()) / 2);
}
