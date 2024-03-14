// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-sigma-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.Sigma0)]]
Rcpp::List Sigma0(const arma::mat& beta, const arma::mat& psi_l,
                  const arma::mat& lambda, const arma::mat& theta_l) {
  int p = beta.n_rows;
  int q = p * p;
  arma::vec psi_vec = arma::vectorise(psi_l * psi_l.t());
  arma::vec sigma_vec(p * p);
  sigma_vec = arma::inv(arma::eye(q, q) - arma::kron(beta, beta)) * psi_vec;
  arma::mat sigma_eta = arma::reshape(sigma_vec, p, p);
  arma::mat sigma_y = lambda * (theta_l * theta_l.t()) * lambda.t();
  return Rcpp::List::create(Rcpp::Named("sigma_y") = sigma_y,
                            Rcpp::Named("sigma_eta") = sigma_eta);
}
