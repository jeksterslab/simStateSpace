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
//'   \eqn{\boldsymbol{\Sigma}}.
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
//' @keywords simStateSpace sim linsde
//' @export
// [[Rcpp::export]]
Rcpp::List LinSDE2SSM(const arma::vec& iota, const arma::mat& phi,
                      const arma::mat& sigma_l, const double delta_t) {
  int p = iota.n_elem;
  arma::mat I = arma::eye(p, p);
  // beta
  arma::mat beta = arma::expmat(phi * delta_t);
  // alpha
  arma::vec alpha = arma::inv(phi) * (beta - I) * iota;
  // psi_l
  arma::mat psi_l = arma::mat(p, p);
  if (arma::all(arma::vectorise(sigma_l) == 0)) {
    psi_l = sigma_l;
  } else {
    arma::mat J = arma::eye(p * p, p * p);
    arma::mat phi_hashtag = arma::kron(phi, I) + arma::kron(I, phi);
    arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
    arma::vec psi_vec = arma::inv(phi_hashtag) *
                        (arma::expmat(phi_hashtag * delta_t) - J) * sigma_vec;
    psi_l = arma::chol(arma::reshape(psi_vec, p, p), "lower");
  }
  // output
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi_l") = psi_l);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed0)]]
Rcpp::List SimSSMFixed0(const int n, const int time, const double delta_t,
                        const arma::vec& mu0, const arma::mat& sigma0_l,
                        const arma::vec& alpha, const arma::mat& beta,
                        const arma::mat& psi_l, const arma::vec& nu,
                        const arma::mat& lambda, const arma::mat& theta_l) {
  // Step 1: Determine dimensions
  int p = mu0.n_elem;  // number of latent variables
  int k = nu.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(p));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(k));
    // Step 3.3: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(p));
      y.col(t) = nu + (lambda * eta.col(t)) + (theta_l * arma::randn(k));
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
Rcpp::List SimSSMFixed1(const int n, const int time, const double delta_t,
                        const arma::vec& mu0, const arma::mat& sigma0_l,
                        const arma::vec& alpha, const arma::mat& beta,
                        const arma::mat& psi_l, const arma::vec& nu,
                        const arma::mat& lambda, const arma::mat& theta_l,
                        const Rcpp::List& x, const arma::mat& gamma) {
  // Step 1: Determine dimensions
  int p = mu0.n_elem;  // number of latent variables
  int k = nu.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(p)) + (gamma * x_i.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(k));
    // Step 3.3: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(p)) +
                   (gamma * x_i.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) + (theta_l * arma::randn(k));
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
Rcpp::List SimSSMFixed2(const int n, const int time, const double delta_t,
                        const arma::vec& mu0, const arma::mat& sigma0_l,
                        const arma::vec& alpha, const arma::mat& beta,
                        const arma::mat& psi_l, const arma::vec& nu,
                        const arma::mat& lambda, const arma::mat& theta_l,
                        const Rcpp::List& x, const arma::mat& gamma,
                        const arma::mat& kappa) {
  // Step 1: Determine dimensions
  int p = mu0.n_elem;  // number of latent variables
  int k = nu.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(p)) + (gamma * x_i.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(k)) +
               (kappa * x_i.col(0));
    // Step 3.3: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(p)) +
                   (gamma * x_i.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) + (theta_l * arma::randn(k)) +
                 (kappa * x_i.col(t));
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
Rcpp::List SimSSMIVary0(const int n, const int time, const double delta_t,
                        const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                        const Rcpp::List& alpha, const Rcpp::List& beta,
                        const Rcpp::List& psi_l, const Rcpp::List& nu,
                        const Rcpp::List& lambda, const Rcpp::List& theta_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::vec id(time, arma::fill::zeros);
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
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) =
          alpha_i + (beta_i * eta.col(t - 1)) + (psi_l_i * arma::randn(p));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k));
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
Rcpp::List SimSSMIVary1(const int n, const int time, const double delta_t,
                        const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                        const Rcpp::List& alpha, const Rcpp::List& beta,
                        const Rcpp::List& psi_l, const Rcpp::List& nu,
                        const Rcpp::List& lambda, const Rcpp::List& theta_l,
                        const Rcpp::List& x, const Rcpp::List& gamma) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
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
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p)) + (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(p)) + (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k));
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
Rcpp::List SimSSMIVary2(const int n, const int time, const double delta_t,
                        const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                        const Rcpp::List& alpha, const Rcpp::List& beta,
                        const Rcpp::List& psi_l, const Rcpp::List& nu,
                        const Rcpp::List& lambda, const Rcpp::List& theta_l,
                        const Rcpp::List& x, const Rcpp::List& gamma,
                        const Rcpp::List& kappa) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
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
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p)) + (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k)) +
               (kappa_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(p)) + (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k)) +
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
Rcpp::List SimSSMLatFixed0(const int n, const int time, const double delta_t,
                           const arma::vec& mu0, const arma::mat& sigma0_l,
                           const arma::vec& alpha, const arma::mat& beta,
                           const arma::mat& psi_l) {
  // Step 1: Determine dimensions
  int p = mu0.n_elem;  // number of latent variables
  int k = p;           // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(p));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(p));
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
Rcpp::List SimSSMLatFixed1(const int n, const int time, const double delta_t,
                           const arma::vec& mu0, const arma::mat& sigma0_l,
                           const arma::vec& alpha, const arma::mat& beta,
                           const arma::mat& psi_l, const Rcpp::List& x,
                           const arma::mat& gamma) {
  // Step 1: Determine dimensions
  int p = mu0.n_elem;  // number of latent variables
  int k = p;           // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(p)) + (gamma * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(p)) +
                   (gamma * x_i.col(t));
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
Rcpp::List SimSSMLatIVary0(const int n, const int time, const double delta_t,
                           const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                           const Rcpp::List& alpha, const Rcpp::List& beta,
                           const Rcpp::List& psi_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = p;             // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) =
          alpha_i + (beta_i * eta.col(t - 1)) + (psi_l_i * arma::randn(p));
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
Rcpp::List SimSSMLatIVary1(const int n, const int time, const double delta_t,
                           const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                           const Rcpp::List& alpha, const Rcpp::List& beta,
                           const Rcpp::List& psi_l, const Rcpp::List& x,
                           const Rcpp::List& gamma) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = p;             // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = alpha[i];
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];
    arma::mat gamma_i = gamma[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p)) + (gamma_i * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(p)) + (gamma_i * x_i.col(t));
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
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lin-sde-i-vary-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary0)]]
Rcpp::List SimSSMLinSDEIVary0(const int n, const int time, const double delta_t,
                              const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l,
                              const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat J = arma::eye<arma::mat>(p * p, p * p);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::vec id(time, arma::fill::zeros);
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
    if (ou) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
                          (arma::expmat(phi_hashtag_i * delta_t) - J) *
                          sigma_vec_i;
    arma::mat psi_l_i = arma::chol(arma::reshape(psi_vec_i, p, p), "lower");
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) =
          alpha_i + (beta_i * eta.col(t - 1)) + (psi_l_i * arma::randn(p));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k));
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
Rcpp::List SimSSMLinSDEIVary1(const int n, const int time, const double delta_t,
                              const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l, const Rcpp::List& x,
                              const Rcpp::List& gamma, const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat J = arma::eye<arma::mat>(p * p, p * p);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
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
    if (ou) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
                          (arma::expmat(phi_hashtag_i * delta_t) - J) *
                          sigma_vec_i;
    arma::mat psi_l_i = arma::chol(arma::reshape(psi_vec_i, p, p), "lower");
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p)) + (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(p)) + (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k));
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
Rcpp::List SimSSMLinSDEIVary2(const int n, const int time, const double delta_t,
                              const Rcpp::List& mu0, const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l, const Rcpp::List& x,
                              const Rcpp::List& gamma, const Rcpp::List& kappa,
                              const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat J = arma::eye<arma::mat>(p * p, p * p);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
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
    if (ou) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
                          (arma::expmat(phi_hashtag_i * delta_t) - J) *
                          sigma_vec_i;
    arma::mat psi_l_i = arma::chol(arma::reshape(psi_vec_i, p, p), "lower");
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(p)) + (gamma_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k)) +
               (kappa_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(p)) + (gamma_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k)) +
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
// edit .setup/cpp/simStateSpace-sim-ssm-mu-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.Mu0)]]
Rcpp::List Mu0(const arma::vec& alpha, const arma::mat& beta,
               const arma::vec& nu) {
  int p = beta.n_rows;
  arma::vec mu_eta = arma::inv(arma::eye(p, p) - beta) * alpha;
  arma::vec mu_y = nu;
  return Rcpp::List::create(Rcpp::Named("mu_y") = mu_y,
                            Rcpp::Named("mu_eta") = mu_eta);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-sigma-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.Sigma0)]]
Rcpp::List Sigma0(const arma::mat& beta, const arma::mat& psi_l,
                  const arma::mat& lambda, const arma::mat& theta_l) {
  int p = beta.n_rows;
  int q = p * p;
  arma::vec psi_vec = arma::vectorise(psi_l * psi_l.t());
  arma::vec sigma_vec(p * p);
  sigma_vec = arma::inv(arma::eye(q, q) - arma::kron(beta, beta)) * psi_vec;
  arma::mat sigma_eta = arma::reshape(sigma_vec, p, p);
  arma::mat sigma_y = lambda * (theta_l * theta_l.t()) * lambda.t();
  return Rcpp::List::create(Rcpp::Named("sigma_y") = sigma_y,
                            Rcpp::Named("sigma_eta") = sigma_eta);
}
