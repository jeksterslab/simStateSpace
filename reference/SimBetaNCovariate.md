# Simulate Transition Matrices with a Covariate from the Multivariate Normal Distribution

This function simulates random transition matrices from a multivariate
normal distribution, allowing the mean transition matrix to vary as a
linear function of a covariate. The function ensures that the generated
transition matrices are stationary using
[`TestStationarity()`](https://github.com/jeksterslab/simStateSpace/reference/TestStationarity.md)
with a rejection sampling approach.

## Usage

``` r
SimBetaNCovariate(
  n,
  beta0,
  vcov_beta_vec_l,
  beta1,
  x,
  margin = 1,
  beta_lbound = NULL,
  beta_ubound = NULL,
  bound = FALSE,
  max_iter = 100000L
)
```

## Arguments

- n:

  Positive integer. Number of replications.

- beta0:

  Numeric matrix. Baseline transition matrix \\\boldsymbol{\beta}\_0\\
  corresponding to \\\mathbf{x} = \mathbf{0}\\.

- vcov_beta_vec_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_beta_vec))`) of
  the sampling variance-covariance matrix of \\\mathrm{vec} \left(
  \boldsymbol{\beta} \right)\\.

- beta1:

  Numeric matrix. Matrix of covariate effects mapping \\\mathbf{x}\\ to
  \\\mathrm{vec}(\boldsymbol{\beta})\\.

- x:

  List of numeric vectors. Covariate values.

- margin:

  Numeric scalar specifying the stationarity threshold. Values less than
  1 indicate stricter stationarity criteria.

- beta_lbound:

  Optional numeric matrix of same dim as `beta`. Use NA for no lower
  bound.

- beta_ubound:

  Optional numeric matrix of same dim as `beta`. Use NA for no upper
  bound.

- bound:

  Logical; if TRUE, resample until all elements respect bounds (NA
  bounds ignored).

- max_iter:

  Safety cap on resampling attempts per draw.

## Value

Returns a list of random transition matrices.

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
beta0 <- matrix(
  data = c(
    0.7, 0.5, -0.1,
    0.0, 0.6, 0.4,
    0, 0, 0.5
  ),
  nrow = 3
)
vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
# One scalar covariate per replication
beta1 <- matrix(data = 0, nrow = 9, ncol = 1)
beta1[1, 1] <- 0.10  # x shifts beta[1,1]
x <- list(c(0), c(1), c(-1), c(0.5), c(2))

SimBetaNCovariate(
  n = n,
  beta0 = beta0,
  vcov_beta_vec_l = vcov_beta_vec_l,
  beta1 = beta1,
  x = x
)
#> [[1]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.74284297 0.01446261 -0.002043949
#> [2,]  0.48430026 0.58948165  0.001537672
#> [3,] -0.07907666 0.37758361  0.495462488
#> 
#> [[2]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.8043460 0.004087614 -0.02657330
#> [2,]  0.5009306 0.610956032  0.06268561
#> [3,] -0.1370391 0.420064618  0.53629846
#> 
#> [[3]]
#>            [,1]       [,2]       [,3]
#> [1,]  0.6340902 0.01974464 0.02682031
#> [2,]  0.5014228 0.60876051 0.04162174
#> [3,] -0.1113340 0.35350421 0.52992448
#> 
#> [[4]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7248926 -0.07374598 0.002912312
#> [2,]  0.5376623  0.58398763 0.038392832
#> [3,] -0.1181411  0.35794732 0.469290227
#> 
#> [[5]]
#>            [,1]        [,2]         [,3]
#> [1,]  0.8453100 -0.01530127 0.0595277334
#> [2,]  0.4995039  0.64748942 0.0006538527
#> [3,] -0.1704366  0.39839450 0.5065986906
#> 
```
