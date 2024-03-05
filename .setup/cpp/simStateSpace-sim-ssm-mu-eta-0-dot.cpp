// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-mu-eta-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Mean Vector of the Latent Variables
//' for a Stationary Process
//'
//' @details The mean vector of the latent variables
//'   \eqn{\boldsymbol{\mu}_{\boldsymbol{\eta}, \boldsymbol{\eta}}}
//'   is given by
//'   \deqn{
//'     \boldsymbol{\mu}_{\boldsymbol{\eta}, \boldsymbol{\eta}}
//'     =
//'     \left(
//'       \mathbf{I} - \boldsymbol{\beta}
//'     \right)^{-1}
//'     \boldsymbol{\alpha}
//'   }
//' @inheritParams SimSSMFixed
//' @examples
//' p <- 3
//' beta <- 0.50 * diag(p)
//' alpha <- rep(x = 0.50, times = p)
//' MuEta0(beta = beta, alpha = alpha)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim ssm
//' @export
// [[Rcpp::export]]
arma::vec MuEta0(const arma::mat& beta, const arma::vec& alpha) {
  int p = beta.n_rows;
  return arma::inv(arma::eye(p, p) - beta) * alpha;
}
