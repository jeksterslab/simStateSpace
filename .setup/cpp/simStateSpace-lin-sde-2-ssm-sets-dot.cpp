// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-2-ssm-sets-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]

// Discretize the individual-varying linear SDE once per parameter set.
// Person i uses iota, phi, and sigma_l at i modulo each list's length, so
// the sets repeat with the least common multiple m of the lengths (at most
// n), and person i uses set i modulo m. With `ou = true`, `iota` holds mu
// and the intercept is -phi mu.
// [[Rcpp::export(.LinSDE2SSMSets)]]
Rcpp::List LinSDE2SSMSets(const arma::uword n, const Rcpp::List& iota,
                          const Rcpp::List& phi, const Rcpp::List& sigma_l,
                          const double delta_t, const bool ou) {
  arma::uword m = 1;
  for (const arma::uword length : {static_cast<arma::uword>(iota.size()),
                                   static_cast<arma::uword>(phi.size()),
                                   static_cast<arma::uword>(sigma_l.size())}) {
    m = std::lcm(m, length);
    if (m >= n) {
      m = n;
      break;
    }
  }
  Rcpp::List alpha(m), beta(m), psi_l(m);
  for (arma::uword s = 0; s < m; ++s) {
    const arma::mat phi_s = phi[s % phi.size()];
    arma::vec iota_s = Rcpp::as<arma::vec>(iota[s % iota.size()]);
    if (ou && !iota_s.is_zero()) iota_s = -phi_s * iota_s;
    arma::vec alpha_s;
    arma::mat beta_s, psi_l_s;
    DiscretizeLinSDE(iota_s, phi_s,
                     Rcpp::as<arma::mat>(sigma_l[s % sigma_l.size()]), delta_t,
                     alpha_s, beta_s, psi_l_s);
    alpha[s] = alpha_s;
    beta[s] = beta_s;
    psi_l[s] = psi_l_s;
  }
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi_l") = psi_l);
}
