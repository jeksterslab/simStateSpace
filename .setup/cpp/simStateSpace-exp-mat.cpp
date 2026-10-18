// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-exp-mat-log-mat.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix Exponential
//'
//' Computes the matrix exponential of a square numeric matrix using the
//' Armadillo linear algebra library.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x A square numeric matrix.
//'
//' @return
//' A numeric matrix containing the matrix exponential of `x`.
//'
//' @examples
//' a <- matrix(
//'   data = c(
//'     -0.40,  0.10,
//'     -0.05, -0.30
//'   ),
//'   nrow = 2,
//'   ncol = 2
//' )
//'
//' b <- ExpMat(x = a)
//' b
//'
//' a <- LogMat(x = b)
//' a
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
arma::mat ExpMat(const arma::mat& x) {
  if (x.n_rows != x.n_cols) {
    Rcpp::stop("`x` must be a square matrix.");
  }
  arma::mat out;
  const bool success = arma::expmat(out, x);
  if (!success) {
    // # nocov start
    Rcpp::stop("Armadillo failed to compute the matrix exponential.");
    // # nocov end
  }
  return out;
}
