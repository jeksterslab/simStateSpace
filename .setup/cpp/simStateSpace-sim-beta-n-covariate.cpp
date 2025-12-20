// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n-covariate.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

//' Simulate Transition Matrices with a Covariate
//' from the Multivariate Normal Distribution
//'
//' This function simulates random transition matrices from a multivariate
//' normal distribution, allowing the mean transition matrix to vary as a
//' linear function of a covariate.
//' The function ensures that the generated transition matrices are stationary
//' using [TestStationarity()] with a rejection sampling approach.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param beta0 Numeric matrix.
//'   Baseline transition matrix \eqn{\boldsymbol{\beta}_0}
//'   corresponding to \eqn{\mathbf{x} = \mathbf{0}}.
//' @param vcov_beta_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_beta_vec))`)
//'   of the sampling variance-covariance matrix of
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\beta} \right)}.
//' @param beta1 Numeric matrix.
//'   Matrix of covariate effects mapping \eqn{\mathbf{x}} to
//'   \eqn{\mathrm{vec}(\boldsymbol{\beta})}.
//' @param x List of numeric vectors.
//'   Covariate values.
//' @param margin Numeric scalar specifying the stationarity threshold.
//'   Values less than 1 indicate stricter stationarity criteria.
//' @param beta_lbound Optional numeric matrix of same dim as `beta`.
//'   Use NA for no lower bound.
//' @param beta_ubound Optional numeric matrix of same dim as `beta`.
//'   Use NA for no upper bound.
//' @param bound Logical;
//'   if TRUE, resample until all elements respect bounds (NA bounds ignored).
//' @param max_iter Safety cap on resampling attempts per draw.
//' @return Returns a list of random transition matrices.
//'
//' @examples
//' n <- 5
//' beta0 <- matrix(
//'   data = c(
//'     0.7, 0.5, -0.1,
//'     0.0, 0.6, 0.4,
//'     0, 0, 0.5
//'   ),
//'   nrow = 3
//' )
//' vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
//' # One scalar covariate per replication
//' beta1 <- matrix(data = 0, nrow = 9, ncol = 1)
//' beta1[1, 1] <- 0.10  # x shifts beta[1,1]
//' x <- list(c(0), c(1), c(-1), c(0.5), c(2))
//'
//' SimBetaNCovariate(
//'   n = n,
//'   beta0 = beta0,
//'   vcov_beta_vec_l = vcov_beta_vec_l,
//'   beta1 = beta1,
//'   x = x
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimBetaNCovariate(
    const arma::uword& n, const arma::mat& beta0,
    const arma::mat& vcov_beta_vec_l, const arma::mat& beta1,
    const Rcpp::List& x, const double margin = 1.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> beta_lbound = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> beta_ubound = R_NilValue,
    const bool bound = false, const arma::uword max_iter = 100000) {
  const arma::uword nr = beta0.n_rows, nc = beta0.n_cols;
  const arma::uword p = nr * nc;

  // if (vcov_beta_vec_l.n_rows != p || vcov_beta_vec_l.n_cols != p) {
  //   Rcpp::stop("vcov_beta_vec_l must be p x p with p = nrow(beta0) *
  //   ncol(beta0).");
  // }

  arma::mat lb, ub;
  arma::umat has_lb_el(nr, nc, arma::fill::zeros);
  arma::umat has_ub_el(nr, nc, arma::fill::zeros);

  const bool has_lb = beta_lbound.isNotNull();
  const bool has_ub = beta_ubound.isNotNull();

  if (has_lb) {
    lb = Rcpp::as<arma::mat>(beta_lbound);
    if (lb.n_rows != nr || lb.n_cols != nc)
      Rcpp::stop("beta_lbound dims must match beta0.");
    // Build mask: 1 where finite lower bound exists
    for (arma::uword i = 0; i < nr; ++i) {
      for (arma::uword j = 0; j < nc; ++j) {
        has_lb_el(i, j) = std::isfinite(lb(i, j)) ? 1u : 0u;
      }
    }
  } else {
    lb.set_size(nr, nc);
    lb.fill(0.0);
  }

  if (has_ub) {
    ub = Rcpp::as<arma::mat>(beta_ubound);
    if (ub.n_rows != nr || ub.n_cols != nc)
      Rcpp::stop("beta_ubound dims must match beta0.");
    // Build mask: 1 where finite upper bound exists
    for (arma::uword i = 0; i < nr; ++i) {
      for (arma::uword j = 0; j < nc; ++j) {
        has_ub_el(i, j) = std::isfinite(ub(i, j)) ? 1u : 0u;
      }
    }
  } else {
    ub.set_size(nr, nc);
    ub.fill(0.0);
  }

  // Vectorized bounds check with masks; quick-reject
  auto bounds_ok = [&](const arma::mat& x) -> bool {
    if (!bound) return true;
    arma::umat low_violate = (x < lb) % has_lb_el;
    arma::umat high_violate = (x > ub) % has_ub_el;
    return !(arma::any(arma::vectorise(low_violate)) ||
             arma::any(arma::vectorise(high_violate)));
  };

  Rcpp::List out(n);
  const arma::vec beta_vec = arma::vectorise(beta0);
  arma::vec z(p, arma::fill::none), beta_vec_i(p, arma::fill::none);
  arma::mat beta_i(nr, nc, arma::fill::none);

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword iter = 0;
    arma::vec x_i = x[i];
    for (;;) {
      if (iter++ >= max_iter) {
        Rcpp::stop(
            "SimBetaNCovariate: exceeded max_iter while drawing beta_i (i=%u). "
            "Relax bounds or stationarity.",
            i + 1);
      }
      z.randn();
      beta_vec_i = beta_vec + (beta1 * x_i) + (vcov_beta_vec_l * z);
      beta_i = arma::reshape(beta_vec_i, nr, nc);

      if (!bounds_ok(beta_i)) continue;

      if (!TestStationarity(beta_i, margin)) continue;

      out[i] = beta_i;
      break;
    }
  }

  return out;
}
