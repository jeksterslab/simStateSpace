// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-phi-hurwitz.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Hurwitz Stability of a Drift Matrix
//'
//' Returns `TRUE` iff the drift matrix \eqn{\boldsymbol{\Phi}}
//' is Hurwitz-stable,
//' i.e., all eigenvalues have real parts strictly less than `-eps`.
//' Setting `eps = 0` enforces the usual strict condition
//' \eqn{\max \Re\{\lambda_i(\boldsymbol{\Phi})\} < 0}.
//' A small positive `eps` (e.g., `1e-12`) can be used
//' to guard against floating-point round-off.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param eps Nonnegative numeric tolerance (default `0.0`).
//'   The test checks \eqn{\Re(\lambda_i) < -\text{eps}} for all eigenvalues.
//'
//' @examples
//' # Unstable example (spectral abscissa >= 0):
//' phi <- matrix(
//'   data = c(
//'     0.10, -0.40,
//'     0.50, 0.20
//'   ),
//'   nrow = 2
//' )
//' TestPhiHurwitz(phi = phi) # FALSE
//'
//' # Stable example (all real parts < 0):
//' phi <- matrix(
//'   data = c(
//'     -0.50, -0.20,
//'      1.00, -0.30
//'   ),
//'   nrow = 2
//' )
//' TestPhiHurwitz(phi = phi) # TRUE
//' TestPhiHurwitz(phi = phi, eps = 1e-12) # also TRUE with tolerance
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestPhiHurwitz(const arma::mat& phi, const double eps = 0.0) {
  arma::cx_vec ev = arma::eig_gen(phi);
  return arma::all(arma::real(ev) < -eps);
}
