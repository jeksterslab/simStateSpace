# Simulate Random Drift Matrices from the Multivariate Normal Distribution and Project to Hurwitz

This function simulates random dirft matrices from the multivariate
normal distribution then projects each draw to the Hurwitz-stable region
using
[`ProjectToHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToHurwitz.md).

## Usage

``` r
SimPhiN2(n, phi, vcov_phi_vec_l, margin = 0.001)
```

## Arguments

- n:

  Positive integer. Number of replications.

- phi:

  Numeric matrix. The drift matrix (\\\boldsymbol{\Phi}\\).

- vcov_phi_vec_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_phi_vec))`) of
  the sampling variance-covariance matrix of \\\mathrm{vec} \left(
  \boldsymbol{\Phi} \right)\\.

- margin:

  Positive numeric. Target buffer so that the spectral abscissa is \\\le
  -\text{margin}\\ (default `1e-3`).

## Value

Returns a list of random drift matrices.

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
[`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md),
[`SimCovDiagN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovDiagN.md),
[`SimCovN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovN.md),
[`SimIotaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimIotaN.md),
[`SimMVN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMVN.md),
[`SimMuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMuN.md),
[`SimNuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimNuN.md),
[`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md),
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
phi <- matrix(
  data = c(
    -0.357, 0.771, -0.450,
    0.0, -0.511, 0.729,
    0, 0, -0.693
  ),
  nrow = 3
)
vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
SimPhiN2(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
#> [[1]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3898905 -0.00553865 -0.01778982
#> [2,]  0.7938146 -0.53407038  0.03052566
#> [3,] -0.5122352  0.74059043 -0.69089644
#> 
#> [[2]]
#>            [,1]       [,2]        [,3]
#> [1,] -0.3676836 -0.0039967 -0.06034504
#> [2,]  0.7744181 -0.4793175  0.03005414
#> [3,] -0.4750508  0.7039543 -0.65819457
#> 
#> [[3]]
#>            [,1]         [,2]         [,3]
#> [1,] -0.3648512  0.005999648  0.004565715
#> [2,]  0.7460663 -0.545095278  0.002514716
#> [3,] -0.4587547  0.718689032 -0.655266780
#> 
#> [[4]]
#>            [,1]         [,2]          [,3]
#> [1,] -0.3646226  0.002226638  0.0001713893
#> [2,]  0.7623997 -0.513546670 -0.0401962759
#> [3,] -0.4712947  0.759137005 -0.7331935022
#> 
#> [[5]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3374115 -0.03631547 -0.03522833
#> [2,]  0.7829520 -0.51424701 -0.02122889
#> [3,] -0.4593276  0.71458045 -0.72769393
#> 
#> [[6]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3639378 -0.02178019  0.03261244
#> [2,]  0.7968141 -0.51485211  0.00289716
#> [3,] -0.4635005  0.69815458 -0.66238104
#> 
#> [[7]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3942297  0.02094066 -0.02158031
#> [2,]  0.7984825 -0.53096047  0.03250162
#> [3,] -0.4819357  0.74745080 -0.70935891
#> 
#> [[8]]
#>            [,1]       [,2]        [,3]
#> [1,] -0.3657105  0.0333104 -0.02528062
#> [2,]  0.8121938 -0.5232222 -0.04990546
#> [3,] -0.4275957  0.7362539 -0.64011685
#> 
#> [[9]]
#>            [,1]         [,2]        [,3]
#> [1,] -0.3276867  0.002224122  0.02757431
#> [2,]  0.7527464 -0.564012342 -0.02370732
#> [3,] -0.4262468  0.724154016 -0.70016594
#> 
#> [[10]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3434785 -0.03717281 -0.026153021
#> [2,]  0.8416886 -0.56042143  0.001184558
#> [3,] -0.4619536  0.76184134 -0.731455246
#> 
```
