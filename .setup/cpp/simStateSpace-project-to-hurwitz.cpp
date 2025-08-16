// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-project-to-hurwitz.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Project Matrix to Hurwitz Stability
//'
//' Shifts a square matrix left on the real axis
//' so that its spectral abscissa (maximum real part of the eigenvalues)
//' is strictly less than `-margin`.
//' This is useful for ensuring that continuous-time drift matrices
//' (e.g. in linear SDEs/state-space models) are Hurwitz-stable.
//' If the matrix already satisfies the margin,
//' it is returned unchanged.
//'
//' The projection is performed by subtracting a multiple of the identity:
//' \deqn{x^\star = x - (\alpha + \text{margin}) I,}
//' where \eqn{\alpha = \max \Re\{\lambda_i(x)\}} is the spectral abscissa.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//' @param margin Positive numeric.
//'   Target buffer inside the Hurwitz region;
//'   the result satisfies
//'   \eqn{\max \Re\{\lambda_i(x^\star)\} \le -\text{margin}}
//'   (default `1e-3`).
//'
//' @return A numeric matrix of the same dimensions as `x`,
//'   shifted if necessary to satisfy the Hurwitz stability constraint.
//'
//' @examples
//' # Unstable (spectral abscissa >= 0):
//' x <- matrix(
//'   data = c(
//'     0.10, -0.40,
//'     0.50, 0.20
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # >= 0
//' SpectralAbscissa(x = ProjectToHurwitz(x = x)) # <= -1e-3 (default margin)
//'
//' # Already Hurwitz-stable is returned unchanged up to numerics:
//' x <- matrix(
//'   data = c(
//'     -0.50, -0.20,
//'      1.00, -0.30
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # < 0
//' identical(ProjectToHurwitz(x = x), x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace stability linsde
//' @export
// [[Rcpp::export]]
arma::mat ProjectToHurwitz(const arma::mat& x, const double margin = 1e-3) {
  arma::mat output = x;
  const double alpha = SpectralAbscissa(output);
  if (alpha >= -margin) {
    output -= (alpha + margin) * arma::eye(output.n_rows, output.n_cols);
  }
  return output;
}
