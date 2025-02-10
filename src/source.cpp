// -----------------------------------------------------------------------------
// edit .setup/cpp/000-forward-declarations.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

bool TestPhi(const arma::mat& phi);

bool TestStability(const arma::mat& x);

bool TestStationarity(const arma::mat& x);
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-lin-sde-2-ssm.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Convert Parameters from the Linear Stochastic Differential Equation Model
//' to State Space Model Parameterization
//'
//' This function converts parameters from
//' the linear stochastic differential equation model
//' to state space model parameterization.
//'
//' @details Let the linear stochastic equation model be given by
//'   \deqn{
//'     \mathrm{d}
//'     \boldsymbol{\eta}_{i, t}
//'     =
//'     \left(
//'       \boldsymbol{\iota}
//'       +
//'       \boldsymbol{\Phi}
//'       \boldsymbol{\eta}_{i, t}
//'     \right)
//'     \mathrm{d} t
//'     +
//'     \boldsymbol{\Sigma}^{\frac{1}{2}}
//'     \mathrm{d}
//'     \mathbf{W}_{i, t}
//'   }
//'   for individual \eqn{i} and time \eqn{t}.
//'   The discrete-time state space model
//'   given below
//'   represents the discrete-time solution
//'   for the linear stochastic differential equation.
//'   \deqn{
//'     \boldsymbol{\eta}_{i, t_{{l_{i}}}}
//'     =
//'     \boldsymbol{\alpha}_{\Delta t_{{l_{i}}}}
//'     +
//'     \boldsymbol{\beta}_{\Delta t_{{l_{i}}}}
//'     \boldsymbol{\eta}_{i, t_{l_{i} - 1}}
//'     +
//'     \boldsymbol{\zeta}_{i, t_{{l_{i}}}},
//'     \quad
//'     \mathrm{with}
//'     \quad
//'     \boldsymbol{\zeta}_{i, t_{{l_{i}}}}
//'     \sim
//'     \mathcal{N}
//'     \left(
//'     \mathbf{0},
//'     \boldsymbol{\Psi}_{\Delta t_{{l_{i}}}}
//'     \right)
//'   }
//'   with
//'   \deqn{
//'       \boldsymbol{\beta}_{\Delta t_{{l_{i}}}}
//'       =
//'       \exp{
//'         \left(
//'           \Delta t
//'           \boldsymbol{\Phi}
//'         \right)
//'       },
//'   }
//'
//'   \deqn{
//'       \boldsymbol{\alpha}_{\Delta t_{{l_{i}}}}
//'       =
//'       \boldsymbol{\Phi}^{-1}
//'       \left(
//'         \boldsymbol{\beta} - \mathbf{I}_{p}
//'       \right)
//'       \boldsymbol{\iota}, \quad \mathrm{and}
//'   }
//'
//'   \deqn{
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Psi}_{\Delta t_{{l_{i}}}}
//'       \right)
//'       =
//'       \left[
//'         \left(
//'           \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'         \right)
//'         +
//'         \left(
//'           \mathbf{I}_{p} \otimes \boldsymbol{\Phi}
//'         \right)
//'       \right]
//'       \left[
//'         \exp
//'         \left(
//'           \left[
//'             \left(
//'               \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'             \right)
//'             +
//'             \left(
//'               \mathbf{I}_{p} \otimes \boldsymbol{\Phi}
//'             \right)
//'           \right]
//'           \Delta t
//'         \right)
//'         -
//'         \mathbf{I}_{p \times p}
//'       \right]
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Sigma}
//'       \right)
//'   }
//'   where \eqn{t} denotes continuous-time processes
//'   that can be defined by any arbitrary time point,
//'   \eqn{t_{l_{i}}} the \eqn{l^\mathrm{th}}
//'   observed measurement occassion for individual \eqn{i},
//'   \eqn{p} the number of latent variables and
//'   \eqn{\Delta t} the time interval.
//'
//' @references
//'   Harvey, A. C. (1990).
//'   Forecasting, structural time series models and the Kalman filter.
//'   Cambridge University Press.
//'   \doi{10.1017/cbo9781107049994}
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param iota Numeric vector.
//'   An unobserved term that is constant over time
//'   (\eqn{\boldsymbol{\iota}}).
//' @param phi Numeric matrix.
//'   The drift matrix
//'   which represents the rate of change of the solution
//'   in the absence of any random fluctuations
//'   (\eqn{\boldsymbol{\Phi}}).
//' @param sigma_l Numeric matrix.
//'   Cholesky factorization (`t(chol(sigma))`)
//'   of the covariance matrix of volatility
//'   or randomness in the process
//'   (\eqn{\boldsymbol{\Sigma}}).
//' @param delta_t Numeric.
//'   Time interval
//'   (\eqn{\Delta_t}).
//'
//' @return Returns a list of state space parameters:
//'   - `alpha`: Numeric vector.
//'     Vector of constant values for the dynamic model
//'     (\eqn{\boldsymbol{\alpha}}).
//'   - `beta`: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     from the previous time point to the current time point.
//'     (\eqn{\boldsymbol{\beta}}).
//'   - `psi_l`: Numeric matrix.
//'     Cholesky factorization (`t(chol(psi))`)
//'     of the process noise covariance matrix
//'     \eqn{\boldsymbol{\Psi}}.
//'
//' @examples
//' p <- 2
//' iota <- c(0.317, 0.230)
//' phi <- matrix(
//'   data = c(
//'    -0.10,
//'    0.05,
//'    0.05,
//'    -0.10
//'  ),
//'  nrow = p
//' )
//' sigma <- matrix(
//'   data = c(
//'     2.79,
//'     0.06,
//'     0.06,
//'     3.27
//'   ),
//'   nrow = p
//' )
//' sigma_l <- t(chol(sigma))
//' delta_t <- 0.10
//'
//' LinSDE2SSM(
//'   iota = iota,
//'   phi = phi,
//'   sigma_l = sigma_l,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace transformation linsde
//' @export
// [[Rcpp::export]]
Rcpp::List LinSDE2SSM(const arma::vec& iota, const arma::mat& phi,
                      const arma::mat& sigma_l, const double delta_t) {
  arma::mat I = arma::eye(iota.n_elem, iota.n_elem);
  arma::mat beta = arma::expmat(phi * delta_t);
  arma::vec alpha(iota.n_elem, arma::fill::none);
  arma::mat psi_l = arma::mat(iota.n_elem, iota.n_elem, arma::fill::none);
  // alpha = iota.is_zero() ? iota : arma::inv(phi) * (beta - I) * iota;
  alpha = iota.is_zero() ? iota : arma::solve(phi, (beta - I) * iota);
  if (sigma_l.is_zero()) {
    psi_l = sigma_l;
  } else {
    arma::mat J =
        arma::eye(iota.n_elem * iota.n_elem, iota.n_elem * iota.n_elem);
    arma::mat phi_hashtag = arma::kron(phi, I) + arma::kron(I, phi);
    arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
    // arma::vec psi_vec = arma::inv(phi_hashtag) * (arma::expmat(phi_hashtag *
    // delta_t) - J) * sigma_vec;
    arma::vec psi_vec = arma::solve(
        phi_hashtag, (arma::expmat(phi_hashtag * delta_t) - J) * sigma_vec);
    psi_l =
        arma::chol(arma::reshape(psi_vec, iota.n_elem, iota.n_elem), "lower");
  }
  // output
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi_l") = psi_l);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Transition Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random transition matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated transition matrices are stationary
//' using [TestStationarity()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param beta Numeric matrix.
//'   The transition matrix (\eqn{\boldsymbol{\beta}}).
//' @param vcov_beta_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_beta_vec))`)
//'   of the sampling variance-covariance matrix
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\beta} \right)}.
//'
//' @examples
//' n <- 10
//' beta <- matrix(
//'   data = c(
//'     0.7, 0.5, -0.1,
//'     0.0, 0.6, 0.4,
//'     0, 0, 0.5
//'   ),
//'   nrow = 3
//' )
//' vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
//' SimBetaN(n = n, beta = beta, vcov_beta_vec_l = vcov_beta_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace ssm
//' @export
// [[Rcpp::export]]
Rcpp::List SimBetaN(const arma::uword& n, const arma::mat& beta,
                    const arma::mat& vcov_beta_vec_l) {
  Rcpp::List output(n);
  arma::vec beta_vec = arma::vectorise(beta);
  arma::vec beta_vec_i(beta.n_rows * beta.n_cols, arma::fill::none);
  arma::mat beta_i(beta.n_rows, beta.n_cols, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    bool run = true;
    while (run) {
      beta_vec_i =
          beta_vec + (vcov_beta_vec_l * arma::randn(beta.n_rows * beta.n_rows));
      beta_i = arma::reshape(beta_vec_i, beta.n_rows, beta.n_cols);
      if (TestStationarity(beta_i)) {
        run = false;
      }
      if (!run) {
        output[i] = beta_i;
      }
    }
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-phi-n.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Random Drift Matrices
//' from the Multivariate Normal Distribution
//'
//' This function simulates random drift matrices
//' from the multivariate normal distribution.
//' The function ensures that the generated drift matrices are stable
//' using [TestPhi()].
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of replications.
//' @param phi Numeric matrix.
//'   The drift matrix (\eqn{\boldsymbol{\Phi}}).
//' @param vcov_phi_vec_l Numeric matrix.
//'   Cholesky factorization (`t(chol(vcov_phi_vec))`)
//'   of the sampling variance-covariance matrix
//'   \eqn{\mathrm{vec} \left( \boldsymbol{\Phi} \right)}.
//'
//' @examples
//' n <- 10
//' phi <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
//' SimPhiN(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
Rcpp::List SimPhiN(const arma::uword& n, const arma::mat& phi,
                   const arma::mat& vcov_phi_vec_l) {
  Rcpp::List output(n);
  arma::vec phi_vec = arma::vectorise(phi);
  arma::vec phi_vec_i(phi.n_rows * phi.n_cols, arma::fill::none);
  arma::mat phi_i(phi.n_rows, phi.n_cols, arma::fill::none);
  for (arma::uword i = 0; i < n; i++) {
    bool run = true;
    while (run) {
      phi_vec_i =
          phi_vec + (vcov_phi_vec_l * arma::randn(phi.n_rows * phi.n_cols));
      phi_i = arma::reshape(phi_vec_i, phi.n_rows, phi.n_cols);
      if (TestPhi(phi_i)) {
        run = false;
      }
      if (!run) {
        output[i] = phi_i;
      }
    }
  }
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed0)]]
Rcpp::List SimSSMFixed0(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(mu0.n_elem));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) =
          alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(mu0.n_elem));
      y.col(t) =
          nu + (lambda * eta.col(t)) + (theta_l * arma::randn(nu.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed1)]]
Rcpp::List SimSSMFixed1(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const Rcpp::List& x,
                        const arma::mat& gamma) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) =
        mu0 + (sigma0_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(t));
      y.col(t) =
          nu + (lambda * eta.col(t)) + (theta_l * arma::randn(nu.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed2)]]
Rcpp::List SimSSMFixed2(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const Rcpp::List& x,
                        const arma::mat& gamma, const arma::mat& kappa) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) =
        mu0 + (sigma0_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem)) +
               (kappa * x_i.col(0));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(nu.n_elem)) + (kappa * x_i.col(t));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-i-vary-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMIVary0)]]
Rcpp::List SimSSMIVary0(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-i-vary-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMIVary1)]]
Rcpp::List SimSSMIVary1(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& x,
                        const Rcpp::List& gamma) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-i-vary-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMIVary2)]]
Rcpp::List SimSSMIVary2(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& x,
                        const Rcpp::List& gamma, const Rcpp::List& kappa) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];
    arma::mat kappa_i = kappa[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) +
               (theta_l_i * arma::randn(nu_i.n_elem)) + (kappa_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem)) +
                 (kappa_i * x_i.col(t));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-fixed-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatFixed0)]]
Rcpp::List SimSSMLatFixed0(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const arma::vec& mu0,
                           const arma::mat& sigma0_l, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_l) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(mu0.n_elem));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) =
          alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(mu0.n_elem));
      y.col(t) = eta.col(t);
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-fixed-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatFixed1)]]
Rcpp::List SimSSMLatFixed1(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const arma::vec& mu0,
                           const arma::mat& sigma0_l, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_l,
                           const Rcpp::List& x, const arma::mat& gamma) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Generate initial condition
    eta.col(0) =
        mu0 + (sigma0_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(mu0.n_elem)) + (gamma * x_i.col(t));
      y.col(t) = eta.col(t);
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-i-vary-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatIVary0)]]
Rcpp::List SimSSMLatIVary0(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = eta.col(t);
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-i-vary-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatIVary1)]]
Rcpp::List SimSSMLatIVary1(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l,
                           const Rcpp::List& x, const Rcpp::List& gamma) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = eta.col(t);
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary0)]]
Rcpp::List SimSSMLinSDEIVary0(const arma::uword& n, const arma::uword& time,
                              const double delta_t, const Rcpp::List& mu0,
                              const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l,
                              const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = iota[i];
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lin-sde-i-vary-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary1)]]
Rcpp::List SimSSMLinSDEIVary1(const arma::uword& n, const arma::uword& time,
                              const double delta_t, const Rcpp::List& mu0,
                              const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l, const Rcpp::List& x,
                              const Rcpp::List& gamma, const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = iota[i];
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lin-sde-i-vary-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary2)]]
Rcpp::List SimSSMLinSDEIVary2(
    const arma::uword& n, const arma::uword& time, const double delta_t,
    const Rcpp::List& mu0, const Rcpp::List& sigma0_l, const Rcpp::List& iota,
    const Rcpp::List& phi, const Rcpp::List& sigma_l, const Rcpp::List& nu,
    const Rcpp::List& lambda, const Rcpp::List& theta_l, const Rcpp::List& x,
    const Rcpp::List& gamma, const Rcpp::List& kappa, const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::mat x_i = x[i];
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = iota[i];
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_i = gamma[i];
    arma::mat kappa_i = kappa[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem)) +
                 (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) +
               (theta_l_i * arma::randn(nu_i.n_elem)) + (kappa_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem)) +
                   (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem)) +
                 (kappa_i * x_i.col(t));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-solve-syl-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SolveSyl)]]
arma::mat SolveSyl(arma::mat A, arma::mat B, arma::mat C) {
  arma::mat X;
  arma::syl(X, A, B, C);
  return X;
}
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
bool TestPhi(const arma::mat& phi) {
  arma::vec phi_diag = phi.diag(0);
  arma::cx_vec eigenvalues_phi = arma::eig_gen(phi);
  return arma::all(arma::real(eigenvalues_phi) < 0) && arma::all(phi_diag <= 0);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-stability.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Stability
//'
//' The function computes the eigenvalues of the input matrix `x`.
//' It checks if the real part of all eigenvalues is negative.
//' If all eigenvalues have negative real parts,
//' the system is considered stable.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric matrix.
//'
//' @examples
//' x <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' TestStability(x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestStability(const arma::mat& x) {
  arma::cx_vec eigenvalues = arma::eig_gen(x);
  return arma::all(arma::real(eigenvalues) < 0);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-stationarity.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Stationarity
//'
//' The function computes the eigenvalues of the input matrix `x`.
//' It checks if all eigenvalues have moduli less than 1.
//' If all eigenvalues have moduli less than 1,
//' the system is considered stationary.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric matrix.
//'
//' @examples
//' x <- matrix(
//'   data = c(0.5, 0.3, 0.2, 0.4),
//'   nrow = 2
//' )
//' TestStationarity(x)
//'
//' x <- matrix(
//'   data = c(0.9, -0.5, 0.8, 0.7),
//'   nrow = 2
//' )
//' TestStationarity(x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test ssm
//' @export
// [[Rcpp::export]]
bool TestStationarity(const arma::mat& x) {
  arma::cx_vec eigenvalues = arma::eig_gen(x);
  return arma::all(arma::abs(eigenvalues) < 1.0);
}
