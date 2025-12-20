# Simulate Data from a State Space Model (Individual-Varying Parameters)

This function simulates data using a state space model. It assumes that
the parameters can vary across individuals.

## Usage

``` r
SimSSMIVary(
  n,
  time,
  delta_t = 1,
  mu0,
  sigma0_l,
  alpha,
  beta,
  psi_l,
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

- alpha:

  List of numeric vectors. Each element of the list is the vector of
  constant values for the dynamic model (\\\boldsymbol{\alpha}\\).

- beta:

  List of numeric matrices. Each element of the list is the transition
  matrix relating the values of the latent variables at the previous to
  the current time point (\\\boldsymbol{\beta}\\).

- psi_l:

  List of numeric matrices. Each element of the list is the Cholesky
  factorization (`t(chol(psi))`) of the covariance matrix of the process
  noise (\\\boldsymbol{\Psi}\\).

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
  [`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md)
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
values. If the length of any of the parameters (`mu0`, `sigma0_l`,
`alpha`, `beta`, `psi_l`, `nu`, `lambda`, `theta_l`, `gamma`, or
`kappa`) is less the `n`, the function will cycle through the available
values.

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, 17(2), 303–332.
[doi:10.1080/10705511003661553](https://doi.org/10.1080/10705511003661553)

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
[`SimSSMLinGrowth()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowth.md),
[`SimSSMLinGrowthIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowthIVary.md),
[`SimSSMLinSDEFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEFixed.md),
[`SimSSMLinSDEIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEIVary.md),
[`SimSSMOUFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md),
[`SimSSMOUIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUIVary.md),
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
# In this example, beta varies across individuals.
set.seed(42)
## number of individuals
n <- 5
## time points
time <- 50
## dynamic structure
p <- 3
mu0 <- list(
  rep(x = 0, times = p)
)
sigma0 <- 0.001 * diag(p)
sigma0_l <- list(
  t(chol(sigma0))
)
alpha <- list(
  rep(x = 0, times = p)
)
beta <- list(
  0.1 * diag(p),
  0.2 * diag(p),
  0.3 * diag(p),
  0.4 * diag(p),
  0.5 * diag(p)
)
psi <- 0.001 * diag(p)
psi_l <- list(
  t(chol(psi))
)
## measurement model
k <- 3
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
ssm <- SimSSMIVary(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 0
)

plot(ssm)




# Type 1
ssm <- SimSSMIVary(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 1,
  x = x,
  gamma = gamma
)

plot(ssm)




# Type 2
ssm <- SimSSMIVary(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
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
