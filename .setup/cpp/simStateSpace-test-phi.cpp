// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-phi.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test the Drift Matrix
//'
//' Both have to be true for the function to return `TRUE`.
//'   - Test that the real part of all eigenvalues of \eqn{\boldsymbol{\Phi}}
//'     are less than zero.
//'   - Test that the diagonal values of \eqn{\boldsymbol{\Phi}}
//'     are between 0 to negative inifinity.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param a_target Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//' @param auto_ubound Numeric scalar specifying the upper bound
//'   for the diagonal elements of \eqn{\boldsymbol{\Phi}}.
//'   Default is `0.0`, requiring all diagonal values to be \eqn{\leq 0}.
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
//' TestPhi(phi = phi)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestPhi(const arma::mat& phi, const double a_target = 0.0,
             const double auto_ubound = 0.0) {
  arma::vec phi_diag = phi.diag(0);
  arma::cx_vec eigenvalues_phi = arma::eig_gen(phi);
  return arma::all(arma::real(eigenvalues_phi) < a_target) &&
         arma::all(phi_diag <= auto_ubound);
}
