# Simulate Transition Matrices from the Multivariate Normal Distribution and Project to Stability

This function simulates random transition matrices from the multivariate
normal distribution then projects each draw to the stability region
using
[`ProjectToStability()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToStability.md).

## Usage

``` r
SimBetaN2(n, beta, vcov_beta_vec_l, margin = 0.98, tol = 1e-12)
```

## Arguments

- n:

  Positive integer. Number of replications.

- beta:

  Numeric matrix. The transition matrix (\\\boldsymbol{\beta}\\).

- vcov_beta_vec_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_beta_vec))`) of
  the sampling variance-covariance matrix of \\\mathrm{vec} \left(
  \boldsymbol{\beta} \right)\\.

- margin:

  Double in \\(0, 1)\\. Target upper bound for the spectral radius
  (default = 0.98).

- tol:

  Small positive double added to the denominator in the scaling factor
  to avoid division by zero (default = 1e-12).

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
[`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md),
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
n <- 10
beta <- matrix(
  data = c(
    0.7, 0.5, -0.1,
    0.0, 0.6, 0.4,
    0, 0, 0.5
  ),
  nrow = 3
)
vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
SimBetaN2(n = n, beta = beta, vcov_beta_vec_l = vcov_beta_vec_l)
#> [[1]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.7220518 0.03249946 -0.01232484
#> [2,]  0.4702043 0.61665058  0.04992734
#> [3,] -0.1305185 0.38812133  0.55434223
#> 
#> [[2]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.67454086 -0.01985415  0.01570994
#> [2,]  0.54074997  0.59686953 -0.02832438
#> [3,] -0.09213611  0.34792962  0.50969799
#> 
#> [[3]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.69062398 0.005225181 -0.02010547
#> [2,]  0.52625414 0.595154569 -0.02030443
#> [3,] -0.09222055 0.414382456  0.47769312
#> 
#> [[4]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.68802190 0.04630413 -0.00768968
#> [2,]  0.52123165 0.61984513 -0.03364059
#> [3,] -0.08394862 0.41551940  0.48020576
#> 
#> [[5]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.73415159 0.02122083 -0.015559204
#> [2,]  0.55256328 0.65052020  0.004007384
#> [3,] -0.08774542 0.36608981  0.552494039
#> 
#> [[6]]
#>             [,1]        [,2]       [,3]
#> [1,]  0.73639433 -0.03453232 0.02412758
#> [2,]  0.50179583  0.61801095 0.03740576
#> [3,] -0.05347199  0.40766860 0.50818978
#> 
#> [[7]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.70905180 0.02610376 -0.003110754
#> [2,]  0.41579241 0.59660088 -0.006949438
#> [3,] -0.05697283 0.41507708  0.512086060
#> 
#> [[8]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.72310965 0.07649833 -0.02875536
#> [2,]  0.51979430 0.57440904  0.01673866
#> [3,] -0.04677507 0.39292892  0.50473418
#> 
#> [[9]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.66243996 0.05015854 0.008272799
#> [2,]  0.55001538 0.53426108 0.011466010
#> [3,] -0.09489223 0.45338640 0.535675170
#> 
#> [[10]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.75195525 -0.0343205 -0.01862077
#> [2,]  0.47508497  0.6059596 -0.01193983
#> [3,] -0.08074849  0.4401671  0.53462435
#> 
```
