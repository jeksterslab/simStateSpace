# Print Method for an Object of Class `simstatespace`

Print Method for an Object of Class `simstatespace`

## Usage

``` r
# S3 method for class 'simstatespace'
print(x, ...)
```

## Arguments

- x:

  Object of Class `simstatespace`.

- ...:

  Additional arguments.

## Value

Prints simulated data in long format.

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
## dynamic structure
p <- 3
mu0 <- rep(x = 0, times = p)
sigma0 <- diag(p)
sigma0_l <- t(chol(sigma0))
alpha <- rep(x = 0, times = p)
beta <- 0.50 * diag(p)
psi <- diag(p)
psi_l <- t(chol(psi))
## measurement model
k <- 3
nu <- rep(x = 0, times = k)
lambda <- diag(k)
theta <- 0.50 * diag(k)
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
ssm <- SimSSMFixed(
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

print(ssm)
#> Call:
#> SimSSMFixed(n = n, time = time, mu0 = mu0, sigma0_l = sigma0_l, 
#>     alpha = alpha, beta = beta, psi_l = psi_l, nu = nu, lambda = lambda, 
#>     theta_l = theta_l, type = 0)
#> Use `as.data.frame` or `as.matrix` on the output of `SimSSMFixed`
#> to convert it to a data frame or a matrix.
#> 

# Type 1
ssm <- SimSSMFixed(
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

print(ssm)
#> Call:
#> SimSSMFixed(n = n, time = time, mu0 = mu0, sigma0_l = sigma0_l, 
#>     alpha = alpha, beta = beta, psi_l = psi_l, nu = nu, lambda = lambda, 
#>     theta_l = theta_l, type = 1, x = x, gamma = gamma)
#> Use `as.data.frame` or `as.matrix` on the output of `SimSSMFixed`
#> to convert it to a data frame or a matrix.
#> 

# Type 2
ssm <- SimSSMFixed(
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

print(ssm)
#> Call:
#> SimSSMFixed(n = n, time = time, mu0 = mu0, sigma0_l = sigma0_l, 
#>     alpha = alpha, beta = beta, psi_l = psi_l, nu = nu, lambda = lambda, 
#>     theta_l = theta_l, type = 2, x = x, gamma = gamma, kappa = kappa)
#> Use `as.data.frame` or `as.matrix` on the output of `SimSSMFixed`
#> to convert it to a data frame or a matrix.
#> 
```
