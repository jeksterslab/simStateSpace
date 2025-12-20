# Simulate Data from the Ornstein–Uhlenbeck Model using a State Space Model Parameterization (Fixed Parameters)

This function simulates data from the Ornstein–Uhlenbeck (OU) model
using a state space model parameterization. It assumes that the
parameters remain constant across individuals and over time.

## Usage

``` r
SimSSMOUFixed(
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

  Numeric. Time interval (\\\Delta_t\\).

- mu0:

  Numeric vector. Mean of initial latent variable values
  (\\\boldsymbol{\mu}\_{\boldsymbol{\eta} \mid 0}\\).

- sigma0_l:

  Numeric matrix. Cholesky factorization (`t(chol(sigma0))`) of the
  covariance matrix of initial latent variable values
  (\\\boldsymbol{\Sigma}\_{\boldsymbol{\eta} \mid 0}\\).

- mu:

  Numeric vector. The long-term mean or equilibrium level
  (\\\boldsymbol{\mu}\\).

- phi:

  Numeric matrix. The drift matrix which represents the rate of change
  of the solution in the absence of any random fluctuations
  (\\\boldsymbol{\Phi}\\). It also represents the rate of mean
  reversion, determining how quickly the variable returns to its mean.

- sigma_l:

  Numeric matrix. Cholesky factorization (`t(chol(sigma))`) of the
  covariance matrix of volatility or randomness in the process
  (\\\boldsymbol{\Sigma}\\).

- nu:

  Numeric vector. Vector of intercept values for the measurement model
  (\\\boldsymbol{\nu}\\).

- lambda:

  Numeric matrix. Factor loading matrix linking the latent variables to
  the observed variables (\\\boldsymbol{\Lambda}\\).

- theta_l:

  Numeric matrix. Cholesky factorization (`t(chol(theta))`) of the
  covariance matrix of the measurement error (\\\boldsymbol{\Theta}\\).

- type:

  Integer. State space model type. See Details for more information.

- x:

  List. Each element of the list is a matrix of covariates for each
  individual `i` in `n`. The number of columns in each matrix should be
  equal to `time`.

- gamma:

  Numeric matrix. Matrix linking the covariates to the latent variables
  at current time point (\\\boldsymbol{\Gamma}\\).

- kappa:

  Numeric matrix. Matrix linking the covariates to the observed
  variables at current time point (\\\boldsymbol{\kappa}\\).

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

### Type 0

The measurement model is given by \$\$ \mathbf{y}\_{i, t} =
\boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\eta}\_{i, t} +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( \mathbf{0},
\boldsymbol{\Theta} \right) \$\$ where \\\mathbf{y}\_{i, t}\\,
\\\boldsymbol{\eta}\_{i, t}\\, and \\\boldsymbol{\varepsilon}\_{i, t}\\
are random variables and \\\boldsymbol{\nu}\\, \\\boldsymbol{\Lambda}\\,
and \\\boldsymbol{\Theta}\\ are model parameters. \\\mathbf{y}\_{i, t}\\
represents a vector of observed random variables,
\\\boldsymbol{\eta}\_{i, t}\\ a vector of latent random variables, and
\\\boldsymbol{\varepsilon}\_{i, t}\\ a vector of random measurement
errors, at time \\t\\ and individual \\i\\. \\\boldsymbol{\nu}\\ denotes
a vector of intercepts, \\\boldsymbol{\Lambda}\\ a matrix of factor
loadings, and \\\boldsymbol{\Theta}\\ the covariance matrix of
\\\boldsymbol{\varepsilon}\\.

An alternative representation of the measurement error is given by \$\$
\boldsymbol{\varepsilon}\_{i, t} = \boldsymbol{\Theta}^{\frac{1}{2}}
\mathbf{z}\_{i, t}, \quad \mathrm{with} \quad \mathbf{z}\_{i, t} \sim
\mathcal{N} \left( \mathbf{0}, \mathbf{I} \right) \$\$ where
\\\mathbf{z}\_{i, t}\\ is a vector of independent standard normal random
variables and \\ \left( \boldsymbol{\Theta}^{\frac{1}{2}} \right) \left(
\boldsymbol{\Theta}^{\frac{1}{2}} \right)^{\prime} = \boldsymbol{\Theta}
. \\

The dynamic structure is given by \$\$ \mathrm{d} \boldsymbol{\eta}\_{i,
t} = \boldsymbol{\Phi} \left( \boldsymbol{\eta}\_{i, t} -
\boldsymbol{\mu} \right) \mathrm{d}t + \boldsymbol{\Sigma}^{\frac{1}{2}}
\mathrm{d} \mathbf{W}\_{i, t} \$\$ where \\\boldsymbol{\mu}\\ is the
long-term mean or equilibrium level, \\\boldsymbol{\Phi}\\ is the rate
of mean reversion, determining how quickly the variable returns to its
mean, \\\boldsymbol{\Sigma}\\ is the matrix of volatility or randomness
in the process, and \\\mathrm{d}\boldsymbol{W}\\ is a Wiener process or
Brownian motion, which represents random fluctuations.

### Type 1

The measurement model is given by \$\$ \mathbf{y}\_{i, t} =
\boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\eta}\_{i, t} +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( \mathbf{0},
\boldsymbol{\Theta} \right) . \$\$

The dynamic structure is given by \$\$ \mathrm{d} \boldsymbol{\eta}\_{i,
t} = \boldsymbol{\Phi} \left( \boldsymbol{\eta}\_{i, t} -
\boldsymbol{\mu} \right) \mathrm{d}t + \boldsymbol{\Gamma}
\mathbf{x}\_{i, t} + \boldsymbol{\Sigma}^{\frac{1}{2}} \mathrm{d}
\mathbf{W}\_{i, t} \$\$ where \\\mathbf{x}\_{i, t}\\ represents a vector
of covariates at time \\t\\ and individual \\i\\, and
\\\boldsymbol{\Gamma}\\ the coefficient matrix linking the covariates to
the latent variables.

### Type 2

The measurement model is given by \$\$ \mathbf{y}\_{i, t} =
\boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\eta}\_{i, t} +
\boldsymbol{\kappa} \mathbf{x}\_{i, t} + \boldsymbol{\varepsilon}\_{i,
t}, \quad \mathrm{with} \quad \boldsymbol{\varepsilon}\_{i, t} \sim
\mathcal{N} \left( \mathbf{0}, \boldsymbol{\Theta} \right) \$\$ where
\\\boldsymbol{\kappa}\\ represents the coefficient matrix linking the
covariates to the observed variables.

The dynamic structure is given by \$\$ \mathrm{d} \boldsymbol{\eta}\_{i,
t} = \boldsymbol{\Phi} \left( \boldsymbol{\eta}\_{i, t} -
\boldsymbol{\mu} \right) \mathrm{d}t + \boldsymbol{\Gamma}
\mathbf{x}\_{i, t} + \boldsymbol{\Sigma}^{\frac{1}{2}} \mathrm{d}
\mathbf{W}\_{i, t} . \$\$

### The OU model as a linear stochastic differential equation model

The OU model is a first-order linear stochastic differential equation
model in the form of

\$\$ \mathrm{d} \boldsymbol{\eta}\_{i, t} = \left( \boldsymbol{\iota} +
\boldsymbol{\Phi} \boldsymbol{\eta}\_{i, t} \right) \mathrm{d}t +
\boldsymbol{\Sigma}^{\frac{1}{2}} \mathrm{d} \mathbf{W}\_{i, t} \$\$
where \\\boldsymbol{\mu} = - \boldsymbol{\Phi}^{-1} \boldsymbol{\iota}\\
and, equivalently \\\boldsymbol{\iota} = - \boldsymbol{\Phi}
\boldsymbol{\mu}\\.

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, 17(2), 303–332.
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
Psychological Methods, 16 (4), 468–490.
[doi:10.1037/a0024375](https://doi.org/10.1037/a0024375)

Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the
brownian motion. Physical Review, 36 (5), 823–841.
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
set.seed(42)
## number of individuals
n <- 5
## time points
time <- 50
delta_t <- 0.10
## dynamic structure
p <- 2
mu0 <- c(-3.0, 1.5)
sigma0 <- 0.001 * diag(p)
sigma0_l <- t(chol(sigma0))
mu <- c(5.76, 5.18)
phi <- matrix(
  data = c(
    -0.10,
    0.05,
    0.05,
    -0.10
  ),
  nrow = p
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
sigma_l <- t(chol(sigma))
## measurement model
k <- 2
nu <- rep(x = 0, times = k)
lambda <- diag(k)
theta <- 0.001 * diag(k)
theta_l <- t(chol(theta))
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
gamma <- diag(x = 0.10, nrow = p, ncol = j)
kappa <- diag(x = 0.10, nrow = k, ncol = j)

# Type 0
ssm <- SimSSMOUFixed(
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
ssm <- SimSSMOUFixed(
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
ssm <- SimSSMOUFixed(
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
