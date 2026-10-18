// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-exact-ssm-2-lin-sde-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SSM2LinSDE)]]
Rcpp::List SSM2LinSDE(const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_l, const double delta_t,
                      const double tol) {
  const arma::uword p = beta.n_rows;
  arma::mat I = arma::eye(p, p);
  arma::cx_mat log_beta;
  if (!arma::logmat(log_beta, beta)) {
    Rcpp::stop("Matrix logarithm of beta failed.");
  }
  if (!log_beta.is_finite()) {
    Rcpp::stop("Matrix logarithm of beta is not finite.");
  }
  const double max_imag = arma::abs(arma::imag(log_beta)).max();
  if (!(max_imag <= tol)) {
    Rcpp::stop("Matrix logarithm of beta has a non-negligible imaginary part.");
  }
  arma::mat phi = arma::real(log_beta) / delta_t;
  arma::vec iota;
  if (!arma::solve(iota, IntegratedExponential(phi, delta_t), alpha,
                   arma::solve_opts::no_approx)) {
    // # nocov start
    Rcpp::stop("Continuous-time intercept is not uniquely recoverable.");
    // # nocov end
  }
  const arma::mat psi = psi_l * psi_l.t();
  const arma::mat generator = arma::kron(I, phi) + arma::kron(phi, I);
  arma::vec sigma_vec;
  if (!arma::solve(sigma_vec, IntegratedExponential(generator, delta_t),
                   arma::vectorise(psi), arma::solve_opts::no_approx)) {
    // # nocov start
    Rcpp::stop("Continuous-time diffusion is not uniquely recoverable.");
    // # nocov end
  }
  const arma::mat sigma_l =
      LowerCholPSD(arma::reshape(sigma_vec, p, p), "sigma");
  return Rcpp::List::create(Rcpp::Named("iota") = iota,
                            Rcpp::Named("phi") = phi,
                            Rcpp::Named("sigma_l") = sigma_l);
}
