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
#>            [,1]       [,2]       [,3]
#> [1,]  0.6776931 0.02123165 0.01984513
#> [2,]  0.5245851 0.61605138 0.01551940
#> [3,] -0.1119781 0.44630413 0.49231032
#> 
#> [[2]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.68020576 0.05256328  0.05052020
#> [2,]  0.50736615 0.61225458 -0.03391019
#> [3,] -0.06584841 0.42122083  0.48444080
#> 
#> [[3]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.75249404 0.001795825 0.018010949
#> [2,]  0.45034649 0.646528005 0.007668595
#> [3,] -0.06360567 0.365467682 0.524127584
#> 
#> [[4]]
#>            [,1]        [,2]         [,3]
#> [1,]  0.7081898 -0.08420759 -0.003399117
#> [2,]  0.4501561  0.64302717  0.015077080
#> [3,] -0.0909482  0.42610376  0.496889246
#> 
#> [[5]]
#>             [,1]      [,2]         [,3]
#> [1,]  0.71208606 0.0197943 -0.025590957
#> [2,]  0.45524973 0.6532249 -0.007071082
#> [3,] -0.07689035 0.4764983  0.471244635
#> 
#> [[6]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.7047342 0.05001538 -0.06573892
#> [2,]  0.5120590 0.60510777  0.05338640
#> [3,] -0.1375600 0.45015854  0.50827280
#> 
#> [[7]]
#>             [,1]        [,2]      [,3]
#> [1,]  0.73567517 -0.02491503 0.0059596
#> [2,]  0.51117277  0.61925151 0.0401671
#> [3,] -0.04804475  0.36567950 0.4813792
#> 
#> [[8]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.73462435 0.02649681 -0.01423570
#> [2,]  0.49158655 0.60046216 -0.03370255
#> [3,] -0.06963177 0.41429155  0.51462578
#> 
#> [[9]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.68667642 0.03964171  0.03603814
#> [2,]  0.46425196 0.57565413 -0.04383784
#> [3,] -0.08653574 0.41219682  0.47388074
#> 
#> [[10]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.7067041 0.07280703 0.027558781
#> [2,]  0.5144374 0.56172021 0.006540289
#> [3,] -0.1310961 0.40764590 0.537142750
#> 
```
