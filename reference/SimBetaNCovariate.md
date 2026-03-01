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
[`LinSDEInterceptEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEInterceptEta.md),
[`LinSDEInterceptY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEInterceptY.md),
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
[`SimMVN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMVN.md),
[`SimMuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMuN.md),
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
#>            [,1]        [,2]         [,3]
#> [1,]  0.7209233 -0.02241639 -0.004537512
#> [2,]  0.5144626  0.59795605 -0.016906034
#> [3,] -0.1105184  0.40153767  0.504346010
#> 
#> [[2]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.76296094 0.02006462 0.036298455
#> [2,]  0.50408761 0.57342670 0.004114273
#> [3,] -0.08904397 0.46268561 0.534090175
#> 
#> [[3]]
#>             [,1]        [,2]         [,3]
#> [1,]  0.58866605 -0.04649579  0.029924476
#> [2,]  0.51974464  0.62682031 -0.009478689
#> [3,] -0.09123949  0.44162174  0.474892572
#> 
#> [[4]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7318589 -0.04205268 -0.03070977
#> [2,]  0.4262540  0.60291231 -0.01042189
#> [3,] -0.1160124  0.43839283  0.44531002
#> 
#> [[5]]
#>             [,1]         [,2]         [,3]
#> [1,]  0.82956340 -0.001605503  0.006598691
#> [2,]  0.48469873  0.659527733 -0.053251434
#> [3,] -0.05251058  0.400653853  0.518960943
#> 
```
