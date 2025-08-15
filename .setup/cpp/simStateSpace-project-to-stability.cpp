// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-project-to-stability.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Project Matrix to Stability
//'
//' Scales a square matrix so that its spectral radius is strictly less than
//' 1 by a specified stability margin. This is useful for ensuring that
//' transition matrices in state space or vector autoregressive (VAR) models
//' are stationary. If the matrix is already within the margin, it is returned
//' unchanged.
//'
//' The projection is performed by multiplying the matrix by a constant factor
//' \eqn{c = \frac{\text{margin}}{\rho + \text{tol}}}, where \eqn{\rho} is the
//' spectral radius and \code{tol} is a small positive number to prevent
//' division by zero.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//' @param margin Double in \eqn{(0, 1)}. Target upper bound for the spectral
//'   radius (default = 0.98).
//' @param tol Small positive double added to the denominator in the scaling
//'   factor to avoid division by zero (default = 1e-12).
//'
//' @return A numeric matrix of the same dimensions as \code{x}, scaled if
//'   necessary to satisfy the stability constraint.
//'
//' @examples
//' # Matrix with eigenvalues greater than 1
//' A <- matrix(
//'   data = c(
//'     1.2, 0.3,
//'     0.4, 0.9
//'   ),
//'   nrow = 2
//' )
//' SpectralRadius(A)  # > 1
//' A_stable <- ProjectToStability(A)
//' SpectralRadius(A_stable)  # < 1
//'
//' # Matrix already stable is returned unchanged
//' B <- matrix(
//'   data = c(
//'     0.5, 0.3,
//'     0.2, 0.4
//'   ),
//'   nrow = 2
//' )
//' identical(ProjectToStability(B), B)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace stability ssm
//' @export
// [[Rcpp::export]]
arma::mat ProjectToStability(const arma::mat& x, const double margin = 0.98,
                             const double tol = 1e-12) {
  double rho = SpectralRadius(x);
  if (rho < margin) return x;
  double c = margin / (rho + tol);
  return c * x;
}
