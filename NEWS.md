# simStateSpace 1.2.11

## Patch

* Added the `SimAlphaN()` and `SimIotaN()` functions.
* Added the `LinSDECovEta()`, `LinSDECovY()`, `LinSDEMeanEta()` and `LinSDEMeanY()` functions. Removed the `LinSDECov()` and `LinSDEMean()` functions.
* Added the `SSMCovEta()`, `SSMCovY()`, `SSMMeanEta()` and `SSMMeanY()` functions.
* Added the `SimCovN()` and `SimCovDiagN()` functions.
* Added the `SimBetaN2()` and `SimPhiN2()` functions.
* 

# simStateSpace 1.2.10

## Patch

* Make sure the result of `LinSDECov()` is symmetric.

# simStateSpace 1.2.9

## Patch

* Added the `LinSDECov()` and `LinSDEMean()` functions.

# simStateSpace 1.2.8

## Patch

* Moved bootstrap components to a dedicated package, `bootStateSpace`, for better modularity.

# simStateSpace 1.2.7

## Patch

* Added the `PBSSMFixed()`, `PBSSMVARFixed()`, `PBSSMLinSDEFixed()`, and `PBSSMOUFixed()` functions.
* Added `SystemRequirements: GSL (>= 2.6)` in the DESCRIPTION file for the `dynr` package.

# simStateSpace 1.2.3

## Patch

* Replaced `arma::inv` with `arma::solve`.

# simStateSpace 1.2.2

## Patch

* Added the `TestStationarity()`, `TestStability()`, and `TestPhi()` functions.
* Added the `SimBetaN()` and `SimPhiN()` functions.
* The `SimSSMLinSDEIVary()` and `SimSSMOUIVary()` functions can now accept a matrix of zeros for the argument `sigma_l`.
* The `SimSSMLinSDEIVary()` function can now accept a vector of zeros for the argument `iota`.
* The `SimSSMOUIVary()` function can now accept a vector of zeros for the argument `mu`.

# simStateSpace 1.2.1

## Patch

* The `LinSDE2SSM` function can now accept a matrix of zeros for the argument `sigma_l`.
* The `LinSDE2SSM` function can now accept a vector of zeros for the argument `iota`.

# simStateSpace 1.2.0

## Patch

* Added functions to generate linear stochastic differential equation models.

# simStateSpace 1.1.0

## Minor

* Added functions to generate data for models with covariates.
* Added functions to generate data for linear growth curve models.
* Added `print`, `as.data.frame`, `as.matrix`, and `plot` methods.

# simStateSpace 1.0.1

## Patch

* Updates to package documentation.

# simStateSpace 1.0.0

* Initial release.
