# Simulate Data from the Ornstein-Uhlenbeck Model using a State Space Model Parameterization (Individual-Varying Parameters)

This function simulates data from the Ornstein-Uhlenbeck model using a
state space model parameterization. It assumes that the parameters can
vary across individuals.

## Usage

``` r
SimSSMOUIVary(
  n,
  time,
  delta_t = 1,
  mu0,
  sigma0_l,
  mu,
  phi,
  sigma_l,
  nu,
  lambda,
  theta_l,
  type = 0,
  x = NULL,
  gamma = NULL,
  kappa = NULL
)
```

## Arguments

- n:

  Positive integer. Number of individuals.

- time:

  Positive integer. Number of time points.

- delta_t:

  Numeric. Time interval. The default value is `1.0` with an option to
  use a numeric value for the discretized state space model
  parameterization of the linear stochastic differential equation model.

- mu0:

  List of numeric vectors. Each element of the list is the mean of
  initial latent variable values (\\\boldsymbol{\mu}\_{\boldsymbol{\eta}
  \mid 0}\\).

- sigma0_l:

  List of numeric matrices. Each element of the list is the Cholesky
  factorization (`t(chol(sigma0))`) of the covariance matrix of initial
  latent variable values (\\\boldsymbol{\Sigma}\_{\boldsymbol{\eta} \mid
  0}\\).

- mu:

  List of numeric vectors. Each element of the list is the long-term
  mean or equilibrium level (\\\boldsymbol{\mu}\\).

- phi:

  List of numeric matrix. Each element of the list is the drift matrix
  which represents the rate of change of the solution in the absence of
  any random fluctuations (\\\boldsymbol{\Phi}\\). It also represents
  the rate of mean reversion, determining how quickly the variable
  returns to its mean.

- sigma_l:

  List of numeric matrix. Each element of the list is the Cholesky
  factorization (`t(chol(sigma))`) of the covariance matrix of
  volatility or randomness in the process \\\boldsymbol{\Sigma}\\.

- nu:

  List of numeric vectors. Each element of the list is the vector of
  intercept values for the measurement model (\\\boldsymbol{\nu}\\).

- lambda:

  List of numeric matrices. Each element of the list is the factor
  loading matrix linking the latent variables to the observed variables
  (\\\boldsymbol{\Lambda}\\).

- theta_l:

  List of numeric matrices. Each element of the list is the Cholesky
  factorization (`t(chol(theta))`) of the covariance matrix of the
  measurement error (\\\boldsymbol{\Theta}\\).

- type:

  Integer. State space model type. See Details in
  [`SimSSMOUFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md)
  for more information.

- x:

  List. Each element of the list is a matrix of covariates for each
  individual `i` in `n`. The number of columns in each matrix should be
  equal to `time`.

- gamma:

  List of numeric matrices. Each element of the list is the matrix
  linking the covariates to the latent variables at current time point
  (\\\boldsymbol{\Gamma}\\).

- kappa:

  List of numeric matrices. Each element of the list is the matrix
  linking the covariates to the observed variables at current time point
  (\\\boldsymbol{\kappa}\\).

## Value

Returns an object of class `simstatespace` which is a list with the
following elements:

- `call`: Function call.

- `args`: Function arguments.

- `data`: Generated data which is a list of length `n`. Each element of
  `data` is a list with the following elements:

  - `id`: A vector of ID numbers with length `l`, where `l` is the value
    of the function argument `time`.

  - `time`: A vector time points of length `l`.

  - `y`: A `l` by `k` matrix of values for the manifest variables.

  - `eta`: A `l` by `p` matrix of values for the latent variables.

  - `x`: A `l` by `j` matrix of values for the covariates (when
    covariates are included).

- `fun`: Function used.

## Details

Parameters can vary across individuals by providing a list of parameter
values. If the length of any of the parameters (`mu0`, `sigma0_l`, `mu`,
`phi`, `sigma_l`, `nu`, `lambda`, `theta_l`, `gamma`, or `kappa`) is
less the `n`, the function will cycle through the available values.

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, 17(2), 303-332.
[doi:10.1080/10705511003661553](https://doi.org/10.1080/10705511003661553)

Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
Continuous-time dynamic models: Connections to structural equation
models and other discrete-time models. In R. H. Hoyle (Ed.), Handbook of
structural equation modeling (2nd ed.). The Guilford Press.

Harvey, A. C. (1990). Forecasting, structural time series models and the
Kalman filter. Cambridge University Press.
[doi:10.1017/cbo9781107049994](https://doi.org/10.1017/cbo9781107049994)

Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011). A hierarchical
latent stochastic differential equation model for affective dynamics.
Psychological Methods, 16 (4), 468-490.
[doi:10.1037/a0024375](https://doi.org/10.1037/a0024375)

Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the
brownian motion. Physical Review, 36 (5), 823-841.
[doi:10.1103/physrev.36.823](https://doi.org/10.1103/physrev.36.823)

## See also

Other Simulation of State Space Models Data Functions:
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md),
[`LinSDECovEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovEta.md),
[`LinSDECovY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovY.md),
[`LinSDEMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanEta.md),
[`LinSDEMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanY.md),
[`ProjectToHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToHurwitz.md),
[`ProjectToStability()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToStability.md),
[`SSMCovEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovEta.md),
[`SSMCovY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovY.md),
[`SSMInterceptEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMInterceptEta.md),
[`SSMInterceptY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMInterceptY.md),
[`SSMMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanEta.md),
[`SSMMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanY.md),
[`SimAlphaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimAlphaN.md),
[`SimBetaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN.md),
[`SimBetaN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN2.md),
[`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md),
[`SimCovDiagN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovDiagN.md),
[`SimCovN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovN.md),
[`SimIotaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimIotaN.md),
[`SimNuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimNuN.md),
[`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md),
[`SimPhiN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN2.md),
[`SimPhiNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiNCovariate.md),
[`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md),
[`SimSSMIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMIVary.md),
[`SimSSMLinGrowth()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowth.md),
[`SimSSMLinGrowthIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowthIVary.md),
[`SimSSMLinSDEFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEFixed.md),
[`SimSSMLinSDEIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEIVary.md),
[`SimSSMOUFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md),
[`SimSSMVARFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARFixed.md),
[`SimSSMVARIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARIVary.md),
[`SpectralRadius()`](https://github.com/jeksterslab/simStateSpace/reference/SpectralRadius.md),
[`TestPhi()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhi.md),
[`TestPhiHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhiHurwitz.md),
[`TestStability()`](https://github.com/jeksterslab/simStateSpace/reference/TestStability.md),
[`TestStationarity()`](https://github.com/jeksterslab/simStateSpace/reference/TestStationarity.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# prepare parameters
# In this example, phi varies across individuals.
set.seed(42)
## number of individuals
n <- 5
## time points
time <- 50
delta_t <- 0.10
## dynamic structure
p <- 2
mu0 <- list(
  c(-3.0, 1.5)
)
sigma0 <- 0.001 * diag(p)
sigma0_l <- list(
  t(chol(sigma0))
)
mu <- list(
  c(5.76, 5.18)
)
phi <- list(
  -0.1 * diag(p),
  -0.2 * diag(p),
  -0.3 * diag(p),
  -0.4 * diag(p),
  -0.5 * diag(p)
)
sigma <- matrix(
  data = c(
    2.79,
    0.06,
    0.06,
    3.27
  ),
  nrow = p
)
sigma_l <- list(
  t(chol(sigma))
)
## measurement model
k <- 2
nu <- list(
  rep(x = 0, times = k)
)
lambda <- list(
  diag(k)
)
theta <- 0.001 * diag(k)
theta_l <- list(
  t(chol(theta))
)
## covariates
j <- 2
x <- lapply(
  X = seq_len(n),
  FUN = function(i) {
    matrix(
      data = stats::rnorm(n = time * j),
      nrow = j,
      ncol = time
    )
  }
)
gamma <- list(
  diag(x = 0.10, nrow = p, ncol = j)
)
kappa <- list(
  diag(x = 0.10, nrow = k, ncol = j)
)

# Type 0
ssm <- SimSSMOUIVary(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  mu = mu,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 0
)

plot(ssm)



# Type 1
ssm <- SimSSMOUIVary(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  mu = mu,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 1,
  x = x,
  gamma = gamma
)

plot(ssm)



# Type 2
ssm <- SimSSMOUIVary(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  mu = mu,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 2,
  x = x,
  gamma = gamma,
  kappa = kappa
)

plot(ssm)


```
