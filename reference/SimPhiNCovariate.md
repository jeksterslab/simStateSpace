# Simulate Random Drift Matrices with a Covariate from the Multivariate Normal Distribution

This function simulates random drift matrices from the multivariate
normal distribution, allowing the mean drift matrix to vary as a linear
function of a covariate The function ensures that the generated drift
matrices are stable using
[`TestPhi()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhi.md).

## Usage

``` r
SimPhiNCovariate(
  n,
  phi0,
  vcov_phi_vec_l,
  phi1,
  x,
  margin = 0,
  auto_ubound = 0,
  phi_lbound = NULL,
  phi_ubound = NULL,
  bound = FALSE,
  max_iter = 100000L
)
```

## Arguments

- n:

  Positive integer. Number of replications.

- phi0:

  Numeric matrix. Baseline drift matrix (\\\boldsymbol{\Phi}\_0\\).

- vcov_phi_vec_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_phi_vec))`) of
  the sampling variance-covariance matrix of \\\mathrm{vec} \left(
  \boldsymbol{\Phi} \right)\\.

- phi1:

  Numeric matrix. Matrix of covariate effects mapping \\\mathbf{x}\\ to
  \\\mathrm{vec}(\boldsymbol{\Phi})\\.

- x:

  List of numeric vectors. Covariate values.

- margin:

  Numeric scalar specifying the stability threshold for the real part of
  the eigenvalues. The default `0.0` corresponds to the imaginary axis;
  values less than `0.0` enforce a stricter stability margin.

- auto_ubound:

  Numeric scalar specifying the upper bound for the diagonal elements of
  \\\boldsymbol{\Phi}\\. Default is `0.0`, requiring all diagonal values
  to be \\\leq 0\\.

- phi_lbound:

  Optional numeric matrix of same dim as `phi`. Use NA for no lower
  bound.

- phi_ubound:

  Optional numeric matrix of same dim as `phi`. Use NA for no upper
  bound.

- bound:

  Logical; if TRUE, resample until all elements respect bounds (NA
  bounds ignored).

- max_iter:

  Safety cap on resampling attempts per draw.

## Value

Returns a list of random drift matrices.

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
[`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md),
[`SimSSMIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMIVary.md),
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
n <- 5
phi0 <- matrix(
  data = c(
    -0.357, 0.771, -0.450,
    0.0, -0.511, 0.729,
    0, 0, -0.693
  ),
  nrow = 3
)
vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
# One scalar covariate per replication
phi1 <- matrix(data = 0, nrow = 9, ncol = 1)
phi1[1, 1] <- 0.10  # x shifts phi[1,1]
x <- list(c(0), c(1), c(-1), c(0.5), c(2))
SimPhiNCovariate(
  n = n,
  phi0 = phi0,
  vcov_phi_vec_l = vcov_phi_vec_l,
  phi1 = phi1,
  x = x
)
#> [[1]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3263810  0.02748248 -0.01996047
#> [2,]  0.7990942 -0.54293567  0.01845080
#> [3,] -0.4872297  0.74994066 -0.71458031
#> 
#> [[2]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.2733589  0.04119382 -0.012222221
#> [2,]  0.7460655 -0.48859568  0.007253938
#> [3,] -0.4587105  0.76231040 -0.718280621
#> 
#> [[3]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.4041169 -0.01825355 -0.053012342
#> [2,]  0.7637093 -0.48724678 -0.004845984
#> [3,] -0.4206867  0.73122412 -0.665425689
#> 
#> [[4]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3141659  0.07068859 -0.04942143
#> [2,]  0.7726766 -0.52295359  0.03284134
#> [3,] -0.4364785  0.69182719 -0.71915302
#> 
#> [[5]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.1954552  0.04566694 -0.03711423
#> [2,]  0.8050238 -0.54217002 -0.04592696
#> [3,] -0.4695914  0.71841152 -0.73917385
#> 
```
