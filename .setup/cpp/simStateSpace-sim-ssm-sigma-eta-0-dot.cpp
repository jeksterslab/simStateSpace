// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-sigma-eta-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Covariance Matrix of the Latent Variables
//' for a Stationary Process
//'
//' @details The column vector of the elements
//'   of the covariance matrix of the latent variables
//'   \eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta}, \boldsymbol{\eta}}}
//'   is given by
//'   \deqn{
//'     \mathrm{Vec}
//'     \left(
//'       \boldsymbol{\Sigma}_{\boldsymbol{\eta}, \boldsymbol{\eta}}
//'     \right)
//'     =
//'     \left(
//'       \mathbf{I} - \boldsymbol{\beta} \otimes \boldsymbol{\beta}
//'     \right)^{-1}
//'     \mathrm{Vec}
//'     \left(
//'       \boldsymbol{\Psi}
//'     \right)
//'   }
//' @inheritParams SimSSMFixed
//' @examples
//' p <- 3
//' beta <- 0.50 * diag(p)
//' psi <- 0.001 * diag(p)
//' psi_l <- t(chol(psi))
//' SigmaEta0(beta = beta, psi_l = psi_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim ssm
//' @export
// [[Rcpp::export]]
arma::mat SigmaEta0(const arma::mat& beta, const arma::mat& psi_l) {
  int p = beta.n_rows;
  int q = p * p;
  arma::vec psi_vec = arma::vectorise(psi_l * psi_l.t());
  arma::vec sigma_vec(p * p);
  sigma_vec = arma::inv(arma::eye(q, q) - arma::kron(beta, beta)) * psi_vec;
  return arma::reshape(sigma_vec, p, p);
}
