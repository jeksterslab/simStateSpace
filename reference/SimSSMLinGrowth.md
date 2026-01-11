# Simulate Data from the Linear Growth Curve Model

This function simulates data from the linear growth curve model.

## Usage

``` r
SimSSMLinGrowth(
  n,
  time,
  mu0,
  sigma0_l,
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

- mu0:

  Numeric vector. A vector of length two. The first element is the mean
  of the intercept, and the second element is the mean of the slope.

- sigma0_l:

  Numeric matrix. Cholesky factorization (`t(chol(sigma0))`) of the
  covariance matrix of the intercept and the slope.

- theta_l:

  Numeric. Square root of the common measurement error variance.

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

The measurement model is given by \$\$ Y\_{i, t} = \left(
\begin{array}{cc} 1 & 0 \\ \end{array} \right) \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\ \eta\_{1\_{i, t}} \\ \end{array} \right) +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( 0, \theta
\right) \$\$ where \\Y\_{i, t}\\, \\\eta\_{0\_{i, t}}\\, \\\eta\_{1\_{i,
t}}\\, and \\\boldsymbol{\varepsilon}\_{i, t}\\ are random variables and
\\\theta\\ is a model parameter. \\Y\_{i, t}\\ is the observed random
variable at time \\t\\ and individual \\i\\, \\\eta\_{0\_{i, t}}\\
(intercept) and \\\eta\_{1\_{i, t}}\\ (slope) form a vector of latent
random variables at time \\t\\ and individual \\i\\, and
\\\boldsymbol{\varepsilon}\_{i, t}\\ a vector of random measurement
errors at time \\t\\ and individual \\i\\. \\\theta\\ is the variance of
\\\boldsymbol{\varepsilon}\\.

The dynamic structure is given by \$\$ \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\ \eta\_{1\_{i, t}} \\ \end{array} \right) = \left(
\begin{array}{cc} 1 & 1 \\ 0 & 1 \\ \end{array} \right) \left(
\begin{array}{c} \eta\_{0\_{i, t - 1}} \\ \eta\_{1\_{i, t - 1}} \\
\end{array} \right) . \$\$

The mean vector and covariance matrix of the intercept and slope are
captured in the mean vector and covariance matrix of the initial
condition given by \$\$ \boldsymbol{\mu}\_{\boldsymbol{\eta} \mid 0} =
\left( \begin{array}{c} \mu\_{\eta\_{0}} \\ \mu\_{\eta\_{1}} \\
\end{array} \right) \quad \mathrm{and,} \$\$

\$\$ \boldsymbol{\Sigma}\_{\boldsymbol{\eta} \mid 0} = \left(
\begin{array}{cc} \sigma^{2}\_{\eta\_{0}} & \sigma\_{\eta\_{0},
\eta\_{1}} \\ \sigma\_{\eta\_{1}, \eta\_{0}} & \sigma^{2}\_{\eta\_{1}}
\\ \end{array} \right) . \$\$

### Type 1

The measurement model is given by \$\$ Y\_{i, t} = \left(
\begin{array}{cc} 1 & 0 \\ \end{array} \right) \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\ \eta\_{1\_{i, t}} \\ \end{array} \right) +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( 0, \theta
\right) . \$\$

The dynamic structure is given by \$\$ \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\ \eta\_{1\_{i, t}} \\ \end{array} \right) = \left(
\begin{array}{cc} 1 & 1 \\ 0 & 1 \\ \end{array} \right) \left(
\begin{array}{c} \eta\_{0\_{i, t - 1}} \\ \eta\_{1\_{i, t - 1}} \\
\end{array} \right) + \boldsymbol{\Gamma} \mathbf{x}\_{i, t} \$\$ where
\\\mathbf{x}\_{i, t}\\ represents a vector of covariates at time \\t\\
and individual \\i\\, and \\\boldsymbol{\Gamma}\\ the coefficient matrix
linking the covariates to the latent variables.

### Type 2

The measurement model is given by \$\$ Y\_{i, t} = \left(
\begin{array}{cc} 1 & 0 \\ \end{array} \right) \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\ \eta\_{1\_{i, t}} \\ \end{array} \right) +
\boldsymbol{\kappa} \mathbf{x}\_{i, t} + \boldsymbol{\varepsilon}\_{i,
t}, \quad \mathrm{with} \quad \boldsymbol{\varepsilon}\_{i, t} \sim
\mathcal{N} \left( 0, \theta \right) \$\$ where \\\boldsymbol{\kappa}\\
represents the coefficient matrix linking the covariates to the observed
variables.

The dynamic structure is given by \$\$ \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\ \eta\_{1\_{i, t}} \\ \end{array} \right) = \left(
\begin{array}{cc} 1 & 1 \\ 0 & 1 \\ \end{array} \right) \left(
\begin{array}{c} \eta\_{0\_{i, t - 1}} \\ \eta\_{1\_{i, t - 1}} \\
\end{array} \right) + \boldsymbol{\Gamma} \mathbf{x}\_{i, t} . \$\$

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, 17(2), 303-332.
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
[`SimSSMIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMIVary.md),
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
set.seed(42)
## number of individuals
n <- 5
## time points
time <- 5
## dynamic structure
p <- 2
mu0 <- c(0.615, 1.006)
sigma0 <- matrix(
  data = c(
    1.932,
    0.618,
    0.618,
    0.587
  ),
  nrow = p
)
sigma0_l <- t(chol(sigma0))
## measurement model
k <- 1
theta <- 0.50
theta_l <- sqrt(theta)
## covariates
j <- 2
x <- lapply(
  X = seq_len(n),
  FUN = function(i) {
    matrix(
      data = rnorm(n = j * time),
      nrow = j
    )
  }
)
gamma <- diag(x = 0.10, nrow = p, ncol = j)
kappa <- diag(x = 0.10, nrow = k, ncol = j)

# Type 0
ssm <- SimSSMLinGrowth(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  theta_l = theta_l,
  type = 0
)

plot(ssm)


# Type 1
ssm <- SimSSMLinGrowth(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  theta_l = theta_l,
  type = 1,
  x = x,
  gamma = gamma
)

plot(ssm)


# Type 2
ssm <- SimSSMLinGrowth(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  theta_l = theta_l,
  type = 2,
  x = x,
  gamma = gamma,
  kappa = kappa
)

plot(ssm)

```
