// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random drift matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated drift matrices are stable
//' using [TestPhi()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//'
//' @examples
//' phi <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.10 * diag(9)))
//' SimPhi(n = 10, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhi(const arma::uword& n, const arma::mat& phi,
                  const arma::mat& vcov_phi_vec_l) {
  Rcpp::List output(n);
  int p = phi.n_rows;
  int q = p * p;
  arma::vec phi_vec = arma::vectorise(phi);
  for (arma::uword i = 0; i < n; i++) {
    bool run = true;
    while (run) {
      arma::vec phi_vec_i = phi_vec + (vcov_phi_vec_l * arma::randn(q));
      arma::mat phi_i = arma::reshape(phi_vec_i, p, p);
      if (TestPhi(phi_i)) {
        run = false;
      }
      if (!run) {
        output[i] = phi_i;
      }
    }
  }
  return output;
}
