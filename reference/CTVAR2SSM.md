# Convert CTVAR Parameters to State Space Parameters

This function converts continuous-time vector autoregressive (CTVAR)
parameters to discrete-time state space model parameters by exact
sampling of the derivative companion state.

## Usage

``` r
CTVAR2SSM(iota, phi, sigma_l, delta_t = 1)
```

## Arguments

- iota:

  Numeric vector of length \\(q\\). Constant term entering the highest
  derivative of the CTVAR model.

- phi:

  List of length \\(p\\). Numeric \\(q \times q\\) matrices containing
  the CTVAR derivative coefficient matrices \\(\mathbf{S}\_{0}, \ldots,
  \mathbf{S}\_{p - 1}\\), in increasing derivative order.

- sigma_l:

  Numeric \\(q \times q\\) matrix. Lower triangular factor of the
  continuous-time diffusion covariance: `sigma_l %*% t(sigma_l)`. Noise
  enters the highest derivative. For a positive-definite covariance
  matrix `sigma`, use `t(chol(sigma))`. Zero and positive-semidefinite
  covariances are allowed.

- delta_t:

  Positive finite numeric scalar. Sampling interval \\(\Delta t\\).

## Value

Returns a list of state space model parameters:

- `alpha`: Numeric vector of length \\(pq\\). Intercept of the sampled
  derivative state.

- `beta`: Numeric \\(pq \times pq\\) matrix. Transition matrix of the
  sampled derivative state.

- `psi_l`: Numeric \\(pq \times pq\\) matrix. Lower triangular factor of
  the interval process noise covariance, reconstructed as
  `psi_l %*% t(psi_l)`.

- `lambda`: Numeric \\(q \times pq\\) matrix. Observation matrix
  selecting the first \\(q\\) state components.

## Details

The transformation uses Equations (4)–(7) and Table 1, block A, of Chow
et al. (2023). For \\(q\\) observed variables and derivative order
\\(p\\), the state has length \\(pq\\) and is ordered as \\((\mathbf{y},
\mathbf{y}^{(1)}, \ldots, \mathbf{y}^{(p-1)})\\). The drift coefficients
are \\(\mathbf{S}\_0, \ldots, \mathbf{S}\_{p-1}\\); the intercept and
diffusion enter the last block of the state.

For orders above one, the sampled observed process generally requires a
state space representation. For order one, this transformation
corresponds to
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md)
with `phi = phi[[1]]` and an identity observation matrix.

## References

Chow, S.-M., Losardo, D., Park, J. J., & Molenaar, P. C. M. (2023).
Continuous-time dynamic models: Connections to structural equation
models and other discrete-time models. In R. H. Hoyle (Ed.), *Handbook
of structural equation modeling* (2nd ed.). The Guilford Press.

## See also

[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md),
[`SSM2LinSDE()`](https://github.com/jeksterslab/simStateSpace/reference/SSM2LinSDE.md)

Other Parameter Transformation Functions:
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md),
[`SSM2LinSDE()`](https://github.com/jeksterslab/simStateSpace/reference/SSM2LinSDE.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# Example 1: CTVAR(1) to a state space model
iota <- c(0.32, 0.23)
phi <- list(
  matrix(
    data = c(
      -0.10, 0.05,
      0.05, -0.10
    ),
    nrow = 2
  )
)
sigma <- matrix(
  data = c(
    2.79, 0.06,
    0.06, 3.27
  ),
  nrow = 2
)
sigma_l <- t(chol(sigma))
CTVAR2SSM(
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  delta_t = 1.0
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
#> $lambda
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 

# Example 2: CTVAR(2) to a state space model
iota <- c(0.32, 0.23)
phi <- list(
  matrix(
    data = c(
      -0.10, 0.05,
      0.05, -0.10
    ),
    nrow = 2
  ),
  matrix(
    data = c(
      -0.30, 0.02,
      0.02, -0.30
    ),
    nrow = 2
  )
)
sigma <- matrix(
  data = c(
    2.79, 0.06,
    0.06, 3.27
  ),
  nrow = 2
)
sigma_l <- t(chol(sigma))
CTVAR2SSM(
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  delta_t = 1.0
)
#> $alpha
#>           [,1]
#> [1,] 0.1450384
#> [2,] 0.1049712
#> [3,] 0.2754204
#> [4,] 0.2003077
#> 
#> $beta
#>             [,1]        [,2]       [,3]       [,4]
#> [1,]  0.95524659  0.02202526 0.84977288 0.01518730
#> [2,]  0.02202526  0.95524659 0.01518730 0.84977288
#> [3,] -0.08421792  0.04096991 0.70061848 0.03446453
#> [4,]  0.04096991 -0.08421792 0.03446453 0.70061848
#> 
#> $psi_l
#>            [,1]       [,2]       [,3]     [,4]
#> [1,] 0.85637465 0.00000000 0.00000000 0.000000
#> [2,] 0.04062557 0.92617013 0.00000000 0.000000
#> [3,] 1.17763926 0.01543529 0.80688762 0.000000
#> [4,] 0.06937570 1.27291298 0.01905466 0.873332
#> 
#> $lambda
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    1    0    0
#> 
```
