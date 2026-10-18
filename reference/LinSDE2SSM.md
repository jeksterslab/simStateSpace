# Convert Linear SDE Parameters to State Space Parameters

This function converts linear stochastic differential equation (SDE)
parameters to discrete-time state space model parameters by exact
discretization over a specified sampling interval.

## Usage

``` r
LinSDE2SSM(iota, phi, sigma_l, delta_t = 1)
```

## Arguments

- iota:

  Numeric vector of length \\(q\\). An unobserved term that is constant
  over time (\\(\boldsymbol{\iota}\\)).

- phi:

  Numeric \\(q \times q\\) matrix. The drift matrix which represents the
  rate of change of the solution in the absence of any random
  fluctuations (\\(\boldsymbol{\Phi}\\)).

- sigma_l:

  Numeric \\(q \times q\\) matrix. Lower triangular factor of the
  continuous-time diffusion covariance: `sigma_l %*% t(sigma_l)`. For a
  positive-definite covariance matrix `sigma`, use `t(chol(sigma))`. The
  factor must be lower triangular with a non-negative diagonal. Zero and
  positive-semidefinite covariances are allowed.

- delta_t:

  Positive finite numeric scalar. Sampling interval \\(\Delta t\\).

## Value

Returns a list of state space parameters:

- `alpha`: Numeric vector. Vector of constant values for the dynamic
  model (\\(\boldsymbol{\alpha}\\)).

- `beta`: Numeric matrix. Transition matrix relating the values of the
  latent variables from the previous time point to the current time
  point (\\(\boldsymbol{\beta}\\)).

- `psi_l`: Numeric matrix. Lower triangular factor of the process noise
  covariance matrix \\(\boldsymbol{\Psi}\\).

## Details

The linear stochastic differential equation is \$\$
\mathrm{d}\boldsymbol{\eta}(t) = (\boldsymbol{\iota} +
\boldsymbol{\Phi}\boldsymbol{\eta}(t)) \mathrm{d}t +
\mathbf{L}\\,\mathrm{d}\mathbf{W}(t), \$\$ where \\(\mathbf{L}\\) is
`sigma_l` and \\(\boldsymbol{\Sigma} = \mathbf{L}\mathbf{L}^{\prime}\\).
Over a sampling interval \\(\Delta t\\), the state equation is \$\$
\boldsymbol{\eta}(t + \Delta t) = \boldsymbol{\alpha} +
\boldsymbol{\beta}\boldsymbol{\eta}(t) + \boldsymbol{\zeta}(t), \qquad
\boldsymbol{\zeta}(t) \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Psi}).
\$\$ The exact discrete-time parameters are \$\$ \boldsymbol{\beta} =
\exp(\Delta t\boldsymbol{\Phi}), \qquad \boldsymbol{\alpha} =
\left\[\int_0^{\Delta t}
\exp(s\boldsymbol{\Phi})\\,\mathrm{d}s\right\]\boldsymbol{\iota}, \$\$
and \$\$ \boldsymbol{\Psi} = \int_0^{\Delta t}
\exp(s\boldsymbol{\Phi})\boldsymbol{\Sigma}
\exp(s\boldsymbol{\Phi}^{\prime})\\,\mathrm{d}s. \$\$ The implementation
evaluates these integrals using block matrix exponentials, so zero and
singular drift matrices are supported. All inputs must contain finite
values and have compatible dimensions.

## References

Harvey, A. C. (1990). Forecasting, structural time series models and the
Kalman filter. Cambridge University Press.
[doi:10.1017/cbo9781107049994](https://doi.org/10.1017/cbo9781107049994)

## See also

[`SSM2LinSDE()`](https://github.com/jeksterslab/simStateSpace/reference/SSM2LinSDE.md),
[`CTVAR2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/CTVAR2SSM.md)

Other Parameter Transformation Functions:
[`CTVAR2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/CTVAR2SSM.md),
[`SSM2LinSDE()`](https://github.com/jeksterslab/simStateSpace/reference/SSM2LinSDE.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# Example 1: Exact sampling of a bivariate linear SDE
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

LinSDE2SSM(
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  delta_t = delta_t
)
#> $alpha
#>           [,1]
#> [1,] 0.3100258
#> [2,] 0.2264506
#> 
#> $beta
#>            [,1]       [,2]
#> [1,] 0.90596870 0.04526072
#> [2,] 0.04526072 0.90596870
#> 
#> $psi_l
#>           [,1]     [,2]
#> [1,] 1.5923846 0.000000
#> [2,] 0.1176299 1.719563
#> 

# Example 2: Zero drift (Brownian motion with a constant drift term)
LinSDE2SSM(
  iota = 0.2,
  phi = matrix(0),
  sigma_l = matrix(0.5),
  delta_t = 0.5
)
#> $alpha
#>      [,1]
#> [1,]  0.1
#> 
#> $beta
#>      [,1]
#> [1,]    1
#> 
#> $psi_l
#>           [,1]
#> [1,] 0.3535534
#> 
```
