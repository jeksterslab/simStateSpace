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
#>            [,1]        [,2]        [,3]
#> [1,]  0.6951546 -0.02030443 -0.01197810
#> [2,]  0.5143825  0.57769312  0.02123165
#> [3,] -0.1201055  0.42458508  0.51605138
#> 
#> [[2]]
#>            [,1]        [,2]       [,3]
#> [1,]  0.7198451 -0.03364059 0.03415159
#> [2,]  0.5155194  0.58020576 0.05256328
#> [3,] -0.1076897  0.40736615 0.51225458
#> 
#> [[3]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7505202 0.004007384 0.036394327
#> [2,]  0.4660898 0.652494039 0.001795825
#> [3,] -0.1155592 0.350346489 0.546528005
#> 
#> [[4]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.71801095 0.03740576  0.009051799
#> [2,]  0.50766860 0.60818978 -0.084207593
#> [3,] -0.07587242 0.35015607  0.543027168
#> 
#> [[5]]
#>            [,1]         [,2]       [,3]
#> [1,]  0.6966009 -0.006949438 0.02310965
#> [2,]  0.5150771  0.612086060 0.01979430
#> [3,] -0.1031108  0.355249733 0.55322493
#> 
#> [[6]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.6744090 0.01673866 -0.03756004
#> [2,]  0.4929289 0.60473418  0.05001538
#> [3,] -0.1287554 0.41205901  0.50510777
#> 
#> [[7]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.6342611 0.01146601  0.05195525
#> [2,]  0.5533864 0.63567517 -0.02491503
#> [3,] -0.0917272 0.41117277  0.51925151
#> 
#> [[8]]
#>            [,1]        [,2]       [,3]
#> [1,]  0.7059596 -0.01193983 0.03036823
#> [2,]  0.5401671  0.63462435 0.02649681
#> [3,] -0.1186208  0.39158655 0.50046216
#> 
#> [[9]]
#>             [,1]        [,2]       [,3]
#> [1,]  0.68576430 -0.01684809 0.01346426
#> [2,]  0.46629745  0.58667642 0.03964171
#> [3,] -0.08537422  0.36425196 0.47565413
#> 
#> [[10]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.7360381 0.04816789 -0.03109606
#> [2,]  0.4561622 0.60670407  0.07280703
#> [3,] -0.1261193 0.41443740  0.46172021
#> 
```
