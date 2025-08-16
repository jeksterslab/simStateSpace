// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi-n-2.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices
//' from the Multivariate Normal Distribution
//' and Project to Hurwitz
//'
//' This function simulates random dirft matrices
//' from the multivariate normal distribution
//' then projects each draw to the Hurwitz-stable region
//' using [ProjectToHurwitz()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//' @param margin Positive numeric.
//'   Target buffer so that the spectral abscissa
//'   is \eqn{\le -\text{margin}} (default `1e-3`).
//' @return Returns a list of random drift matrices.
//'
//' @examples
//' n <- 10
//' phi <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
//' SimPhiN2(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhiN2(const arma::uword& n, const arma::mat& phi,
                    const arma::mat& vcov_phi_vec_l,
                    const double margin = 1e-3) {
  Rcpp::List output(n);
  arma::vec phi_vec = arma::vectorise(phi);
  arma::vec phi_vec_i(phi.n_rows * phi.n_cols, arma::fill::none);
  arma::mat phi_i(phi.n_rows, phi.n_cols, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    phi_vec_i =
        phi_vec + (vcov_phi_vec_l * arma::randn(phi.n_rows * phi.n_cols));
    phi_i = arma::reshape(phi_vec_i, phi.n_rows, phi.n_cols);
    phi_i = ProjectToHurwitz(phi_i, margin);
    output[i] = phi_i;
  }
  return output;
}
