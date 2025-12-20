// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi-n-covariate.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices with a Covariate
//' from the Multivariate Normal Distribution
//'
//' This function simulates random drift matrices
//' from the multivariate normal distribution,
//' allowing the mean drift matrix to vary
//' as a linear function of a covariate
//' The function ensures that the generated drift matrices are stable
//' using [TestPhi()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi0 Numeric matrix.
//'   Baseline drift matrix (\eqn{\boldsymbol{\Phi}_0}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//' @param phi1 Numeric matrix.
//'   Matrix of covariate effects mapping \eqn{\mathbf{x}} to
//'   \eqn{\mathrm{vec}(\boldsymbol{\Phi})}.
//' @param x List of numeric vectors.
//'   Covariate values.
//' @param margin Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//' @param auto_ubound Numeric scalar specifying the upper bound
//'   for the diagonal elements of \eqn{\boldsymbol{\Phi}}.
//'   Default is `0.0`, requiring all diagonal values to be \eqn{\leq 0}.
//' @param phi_lbound Optional numeric matrix of same dim as `phi`.
//'   Use NA for no lower bound.
//' @param phi_ubound Optional numeric matrix of same dim as `phi`.
//'   Use NA for no upper bound.
//' @param bound Logical;
//'   if TRUE, resample until all elements respect bounds (NA bounds ignored).
//' @param max_iter Safety cap on resampling attempts per draw.
//' @return Returns a list of random drift matrices.
//'
//' @examples
//' n <- 5
//' phi0 <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
//' # One scalar covariate per replication
//' phi1 <- matrix(data = 0, nrow = 9, ncol = 1)
//' phi1[1, 1] <- 0.10  # x shifts phi[1,1]
//' x <- list(c(0), c(1), c(-1), c(0.5), c(2))
//' SimPhiNCovariate(
//'   n = n,
//'   phi0 = phi0,
//'   vcov_phi_vec_l = vcov_phi_vec_l,
//'   phi1 = phi1,
//'   x = x
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhiNCovariate(
    const arma::uword& n, const arma::mat& phi0,
    const arma::mat& vcov_phi_vec_l, const arma::mat& phi1, const Rcpp::List& x,
    const double margin = 0.0, const double auto_ubound = 0.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> phi_lbound = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> phi_ubound = R_NilValue,
    const bool bound = false, const arma::uword max_iter = 100000) {
  const arma::uword nr = phi0.n_rows, nc = phi0.n_cols;
  const arma::uword p = nr * nc;

  // if (vcov_phi_vec_l.n_rows != p || vcov_phi_vec_l.n_cols != p) {
  //   Rcpp::stop("vcov_phi_vec_l must be p x p with p = nrow(phi0) *
  //   ncol(phi0).");
  // }

  // Bounds & masks
  arma::mat lb, ub;
  arma::umat has_lb_el(nr, nc, arma::fill::zeros);
  arma::umat has_ub_el(nr, nc, arma::fill::zeros);

  const bool has_lb = phi_lbound.isNotNull();
  const bool has_ub = phi_ubound.isNotNull();

  if (has_lb) {
    lb = Rcpp::as<arma::mat>(phi_lbound);
    if (lb.n_rows != nr || lb.n_cols != nc)
      Rcpp::stop("phi_lbound dims must match phi0.");
    for (arma::uword i = 0; i < nr; ++i)
      for (arma::uword j = 0; j < nc; ++j)
        has_lb_el(i, j) = std::isfinite(lb(i, j)) ? 1u : 0u;
  } else {
    lb.set_size(nr, nc);
    lb.fill(0.0);
  }

  if (has_ub) {
    ub = Rcpp::as<arma::mat>(phi_ubound);
    if (ub.n_rows != nr || ub.n_cols != nc)
      Rcpp::stop("phi_ubound dims must match phi0.");
    for (arma::uword i = 0; i < nr; ++i)
      for (arma::uword j = 0; j < nc; ++j)
        has_ub_el(i, j) = std::isfinite(ub(i, j)) ? 1u : 0u;
  } else {
    ub.set_size(nr, nc);
    ub.fill(0.0);
  }

  auto bounds_ok = [&](const arma::mat& x) -> bool {
    if (!bound) return true;
    arma::umat low_violate = (x < lb) % has_lb_el;
    arma::umat high_violate = (x > ub) % has_ub_el;
    return !(arma::any(arma::vectorise(low_violate)) ||
             arma::any(arma::vectorise(high_violate)));
  };

  Rcpp::List output(n);
  const arma::vec phi_vec = arma::vectorise(phi0);
  arma::vec z(p, arma::fill::none), phi_vec_i(p, arma::fill::none);
  arma::mat phi_i(nr, nc, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword iter = 0;
    arma::vec x_i = x[i];
    for (;;) {
      if (iter++ >= max_iter) {
        Rcpp::stop(
            "SimPhiN: exceeded max_iter while drawing phi_i (i=%u). Relax "
            "bounds or TestPhi().",
            i + 1);
      }
      z.randn();
      phi_vec_i = phi_vec + (phi1 * x_i) + (vcov_phi_vec_l * z);
      phi_i = arma::reshape(phi_vec_i, nr, nc);

      if (!bounds_ok(phi_i)) continue;

      if (!TestPhi(phi_i, margin, auto_ubound)) continue;

      output[i] = phi_i;
      break;
    }
  }

  return output;
}
