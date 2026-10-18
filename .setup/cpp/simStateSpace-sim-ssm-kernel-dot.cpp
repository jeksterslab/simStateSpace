// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-kernel-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// Parameter sets as Armadillo matrices, converted once;
// vectors become column vectors.
std::vector<arma::mat> SimSSMSets(const Rcpp::List& x) {
  std::vector<arma::mat> out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (Rf_isMatrix(x[i])) {
      out[i] = Rcpp::as<arma::mat>(x[i]);
    } else {
      out[i] = Rcpp::as<arma::vec>(x[i]);
    }
  }
  return out;
}

// One simulation kernel for every state space model.
// Each parameter is a list of parameter sets;
// person i uses set i modulo the length of the list.
// Latent-only models pass empty `nu`, `lambda`, and `theta_l` (y = eta).
// Models without covariates pass empty `x` and `gamma`,
// and models without covariates in the measurement model pass empty `kappa`.
// [[Rcpp::export(.SimSSMKernel)]]
Rcpp::List SimSSMKernel(const arma::uword n, const arma::uword time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& x,
                        const Rcpp::List& gamma, const Rcpp::List& kappa) {
  const bool latent = nu.size() == 0;
  const bool covariates = gamma.size() > 0;
  const bool covariates_y = kappa.size() > 0;
  const std::vector<arma::mat> mu0_s = SimSSMSets(mu0);
  const std::vector<arma::mat> sigma0_l_s = SimSSMSets(sigma0_l);
  const std::vector<arma::mat> alpha_s = SimSSMSets(alpha);
  const std::vector<arma::mat> beta_s = SimSSMSets(beta);
  const std::vector<arma::mat> psi_l_s = SimSSMSets(psi_l);
  const std::vector<arma::mat> nu_s = SimSSMSets(nu);
  const std::vector<arma::mat> lambda_s = SimSSMSets(lambda);
  const std::vector<arma::mat> theta_l_s = SimSSMSets(theta_l);
  const std::vector<arma::mat> x_s = SimSSMSets(x);
  const std::vector<arma::mat> gamma_s = SimSSMSets(gamma);
  const std::vector<arma::mat> kappa_s = SimSSMSets(kappa);
  const arma::uword p = mu0_s[0].n_elem;
  const arma::uword k = latent ? 0 : nu_s[0].n_elem;
  // arma::randn draws normal deviates in pairs, so a draw of odd length
  // discards its last deviate. Padding the state and the measurement block
  // to even length makes one draw per person reproduce the draws of earlier
  // versions, one call per block and time point (same seed, same data).
  const arma::uword p2 = p + p % 2;
  const arma::uword k2 = k + k % 2;
  const arma::span state(0, p - 1);
  const arma::vec time_vec = arma::linspace(0, (time - 1) * delta_t, time);
  Rcpp::List output(n);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::mat& beta_i = beta_s[i % beta_s.size()];
    const arma::mat z = arma::randn(p2 + k2, time);
    // everything but beta eta_{t - 1}, then the recursion
    arma::mat eta(p, time);
    eta.col(0) = mu0_s[i % mu0_s.size()] +
                 sigma0_l_s[i % sigma0_l_s.size()] * z(state, 0);
    if (time > 1) {
      const arma::span later(1, time - 1);
      eta.cols(later) = psi_l_s[i % psi_l_s.size()] * z(state, later);
      eta.cols(later).each_col() += arma::vec(alpha_s[i % alpha_s.size()]);
    }
    arma::mat x_i;
    if (covariates) {
      x_i = x_s[i % x_s.size()];
      eta += gamma_s[i % gamma_s.size()] * x_i;
    }
    // plain loops: small matrix-vector products without temporaries
    const double* b = beta_i.memptr();
    double* e = eta.memptr();
    for (arma::uword t = 1; t < time; ++t) {
      const double* prev = e + (t - 1) * p;
      double* curr = e + t * p;
      for (arma::uword c = 0; c < p; ++c) {
        const double prev_c = prev[c];
        const double* b_c = b + c * p;
        for (arma::uword r = 0; r < p; ++r) curr[r] += b_c[r] * prev_c;
      }
    }
    const arma::mat eta_t = eta.t();
    arma::mat y_t;
    if (latent) {
      y_t = eta_t;
    } else {
      y_t = eta_t * lambda_s[i % lambda_s.size()].t() +
            z.rows(p2, p2 + k - 1).t() * theta_l_s[i % theta_l_s.size()].t();
      y_t.each_row() += arma::rowvec(nu_s[i % nu_s.size()].t());
      if (covariates_y) {
        y_t += x_i.t() * kappa_s[i % kappa_s.size()].t();
      }
    }
    arma::vec id(time);
    id.fill(i + 1);
    if (covariates) {
      output[i] = Rcpp::List::create(
          Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
          Rcpp::Named("y") = y_t, Rcpp::Named("eta") = eta_t,
          Rcpp::Named("x") = x_i.t());
    } else {
      output[i] = Rcpp::List::create(
          Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
          Rcpp::Named("y") = y_t, Rcpp::Named("eta") = eta_t);
    }
  }
  return output;
}
