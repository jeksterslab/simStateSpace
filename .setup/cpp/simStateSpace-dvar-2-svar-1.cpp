// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-dvar-2-svar-1.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Convert Parameters from the First Order
//' Discrete-Time Vector Autoregressive Model
//' to Structural Vector Autoregressive Model
//'
//' This function converts parameters from
//' the first order discrete-time vector autoregressive model
//' to structural vector autoregressive model.
//'
//' @inheritParams SimSSMVARFixed
//' @inheritParams LinSDE2SSM
//'
//' @return Returns a list of structural vector autoregressive model parameters:
//'   - `betastar_0`: Numeric matrix.
//'     Matrix of contemporaneous effects.
//'     (\eqn{\boldsymbol{\beta}_{0}}).
//'   - `betastar_1`: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     from the previous time point to the current time point.
//'     (\eqn{\boldsymbol{\beta}_{1}}).
//'   - `psistar_l`: Numeric matrix.
//'     Cholesky factorization (`t(chol(psistar))`)
//'     of the process noise covariance matrix
//'     \eqn{\boldsymbol{\Psi}^{\ast}}.
//'
//' @examples
//' beta <- matrix(
//'   data = c(
//'     0.7,
//'     0.5,
//'     -0.1,
//'     0.0,
//'     0.6,
//'     0.4,
//'     0,
//'     0,
//'     0.5
//'   ),
//'   nrow = 3
//' )
//' psi <- 0.1 * diag(3)
//' psi_l <- t(chol(psi))
//' delta_t <- 0.10
//' DVAR2SVAR1(
//'   beta = beta,
//'   psi_l = psi_l,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace transformation var svar
//' @export
// [[Rcpp::export]]
Rcpp::List DVAR2SVAR1(const arma::mat& beta, const arma::mat& psi_l, const double delta_t) {
  int p = beta.n_rows;
  int order = 1;
  arma::mat I = arma::eye(p, p);
  // betastar
  arma::mat betastar_0 = arma::mat(p, p);
  arma::mat betastar_1 = arma::mat(p, p);
  arma::mat diff = arma::mat(p, p);
  betastar_0 = I - (
    pow(-1, order + 1) * pow(delta_t, -1 * order) * arma::inv(beta)
  );
  diff = I - betastar_0;
  betastar_1 = diff * beta;
  // psistar_l
  arma::mat psistar_l = arma::mat(p, p);
  if (arma::all(arma::vectorise(psi_l) == 0)) {
    psistar_l = psi_l;
  } else{
    psistar_l = arma::chol(diff * (psi_l * psi_l.t()) * diff.t());
  }
  // output
  return Rcpp::List::create(Rcpp::Named("betastar_0") = betastar_0, Rcpp::Named("betastar_1") = betastar_1, Rcpp::Named("psistar_l") = psistar_l);
}
