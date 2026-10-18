# Convert State Space Parameters to Linear SDE Parameters

This function converts parameters from the discrete-time state space
model parameterization to the linear stochastic differential equation
model using the exact matrix-logarithm mapping.

## Usage

``` r
SSM2LinSDE(
  alpha,
  beta,
  psi_l,
  delta_t,
  tol = sqrt(.Machine$double.eps),
  tol_psi_l = sqrt(.Machine$double.eps)
)
```

## Arguments

- alpha:

  Numeric vector of length \\(q\\). Intercept of the discrete-time state
  equation.

- beta:

  Numeric \\(q \times q\\) matrix. Transition matrix of the
  discrete-time state equation.

- psi_l:

  Numeric \\(q \times q\\) matrix. Lower triangular factor of the
  interval process noise covariance: `psi_l %*% t(psi_l)`. For a
  positive-definite covariance matrix `psi`, use `t(chol(psi))`. Zero
  and positive-semidefinite covariances are allowed.

- delta_t:

  Positive finite numeric scalar. Sampling interval \\(\Delta t\\).

- tol:

  Non-negative finite numeric scalar. Tolerance for the largest absolute
  imaginary component of the matrix logarithm.

- tol_psi_l:

  Non-negative finite numeric scalar. Tolerance for checking that
  `psi_l` is lower triangular with a non-negative diagonal.

## Value

Returns a list of linear stochastic differential equation model
parameters:

- `iota`: Numeric vector. An unobserved term that is constant over time
  \\(\boldsymbol{\iota}\\).

- `phi`: Numeric matrix. The drift matrix which represents the rate of
  change of the solution in the absence of any random fluctuations
  \\(\boldsymbol{\Phi}\\).

- `sigma_l`: Numeric matrix. Lower triangular factor of the covariance
  matrix of volatility or randomness in the process
  \\(\boldsymbol{\Sigma}\\).

## Details

The state equation and its continuous-time counterpart are described in
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md).
The drift matrix is computed using the principal matrix logarithm: \$\$
\boldsymbol{\Phi} = \frac{1}{\Delta t}\log(\boldsymbol{\beta}). \$\$ The
intercept and diffusion covariance are recovered by solving \$\$
\boldsymbol{\alpha} = \left\[\int_0^{\Delta t}
\exp(s\boldsymbol{\Phi})\\,\mathrm{d}s\right\]\boldsymbol{\iota} \$\$
and \$\$ \boldsymbol{\Psi} = \int_0^{\Delta t}
\exp(s\boldsymbol{\Phi})\boldsymbol{\Sigma}
\exp(s\boldsymbol{\Phi}^{\prime})\\,\mathrm{d}s, \$\$ where
`psi_l %*% t(psi_l)` gives \\(\boldsymbol{\Psi}\\). The integrated
exponential operators support zero and singular drifts.

The principal matrix-logarithm branch does not search all possible
continuous-time embeddings or resolve temporal aliasing. Inputs are
rejected if the logarithm is not finite (singular `beta`), has an
imaginary component exceeding `tol`, the integral inversions are not
unique, or the recovered diffusion covariance is not positive
semidefinite. All inputs must contain finite values and have compatible
dimensions.

## References

Harvey, A. C. (1990). Forecasting, structural time series models and the
Kalman filter. Cambridge University Press.
[doi:10.1017/cbo9781107049994](https://doi.org/10.1017/cbo9781107049994)

## See also

[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md),
[`CTVAR2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/CTVAR2SSM.md)

Other Parameter Transformation Functions:
[`CTVAR2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/CTVAR2SSM.md),
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# Example 1: Recover the parameters of a bivariate linear SDE
iota <- c(0.32, 0.23)

phi <- matrix(
  data = c(
    -0.10, 0.05,
    0.05, -0.10
  ),
  nrow = 2
)

sigma <- matrix(
  data = c(
    2.79, 0.06,
    0.06, 3.27
  ),
  nrow = 2
)

sigma_l <- t(chol(sigma))
delta_t <- 1.0

ssm <- LinSDE2SSM(
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  delta_t = delta_t
)

SSM2LinSDE(
  alpha = ssm$alpha,
  beta = ssm$beta,
  psi_l = ssm$psi_l,
  delta_t = delta_t
)
#> $iota
#>      [,1]
#> [1,] 0.32
#> [2,] 0.23
#> 
#> $phi
#>       [,1]  [,2]
#> [1,] -0.10  0.05
#> [2,]  0.05 -0.10
#> 
#> $sigma_l
#>            [,1]     [,2]
#> [1,] 1.67032931 0.000000
#> [2,] 0.03592106 1.807957
#> 
```
