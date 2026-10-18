# simStateSpace 1.2.17

## Patch

* Added `ExpMat()` and `LogMat()` for matrix exponentials and real principal matrix logarithms, and `CompanionMat()` for VAR block companion matrices.
* Added `CTVAR2SSM()` for exact sampling of higher-order derivative states and `SSM2LinSDE()` for continuous-time recovery using the principal matrix logarithm, with checks for uniquely recoverable intercepts and positive-semidefinite diffusion.
* Updated `LinSDE2SSM()` to integrate intercepts and diffusion with block matrix exponentials, supporting singular drifts and semidefinite noise. Covariance integration uses interval subdivision to avoid overflow for strongly decaying drifts.
* Added `LinSDECovDeltaT()` for the process noise covariance accumulated over a sampling interval in the linear SDE model.
* `SSM2LinSDE()` rejects a singular `beta` (non-finite matrix logarithm) and has a separate tolerance, `tol_psi_l`, for checking that `psi_l` is lower triangular; `tol` now applies only to the imaginary part of the matrix logarithm.
* Semidefinite noise covariances are factored with an eigendecomposition, so `LinSDE2SSM()` and `SSM2LinSDE()` no longer reject rank-deficient covariances because of rounding error.
* Added input validation to the linear SDE mean, covariance, and intercept functions and expanded transformation tests with analytic reference cases.
* Fixed `SSMInterceptY()` and `LinSDEInterceptY()`, which returned `mean_y + lambda %*% mean_eta` instead of `mean_y - lambda %*% mean_eta` (since 1.2.15).
* Replaced the Monte Carlo and smoke tests with exact references: the simulators are checked against a plain R state space recursion, and the moment functions against their defining equations.
* Fixed `SimSSMLinSDEIVary()` and `SimSSMOUIVary()` failing for singular drifts, semidefinite diffusion, and undamped oscillators: they now discretize with `LinSDE2SSM()`, once per distinct parameter set.
* All `Sim*` functions now share one simulation kernel and one set of input checks, with informative errors for invalid dimensions, non-finite values, and missing covariates. The same seed gives the same data as before.

# simStateSpace 1.2.16

## Patch

* Added the `SimMVN()` function.
* Consolidated the `SimAlphaN()`, `SimIotaN()`, `SimNuN()`, and  `SimMuN()` functions.
* Added the `LinSDEInterceptEta()` and `LinSDEInterceptY()` functions.

# simStateSpace 1.2.15

## Patch

* Added the `SSMInterceptEta()` and `SSMInterceptY()` functions.

# simStateSpace 1.2.14

## Patch

* `CXX_STD = CXX17` in Makevars

# simStateSpace 1.2.13

## Patch

* Added the `SimBetaNCovariate()` and `SimPhiNCovariate()` functions.

# simStateSpace 1.2.12

## Patch

* Added citation to Pesigan, I. J. A., Russell, M. A., & Chow, S.-M. (2025). Inferences and effect sizes for direct, indirect, and total effects in continuous-time mediation models. *Psychological Methods*. https://doi.org/10.1037/met0000779.

# simStateSpace 1.2.11

## Patch

* Added the `SimAlphaN()`, `SimIotaN()`, and  `SimNuN()` functions.
* Added the `LinSDECovEta()`, `LinSDECovY()`, `LinSDEMeanEta()` and `LinSDEMeanY()` functions. Removed the `LinSDECov()` and `LinSDEMean()` functions.
* Added the `SSMCovEta()`, `SSMCovY()`, `SSMMeanEta()` and `SSMMeanY()` functions.
* Added the `SimCovN()` and `SimCovDiagN()` functions.
* Added the `SimBetaN2()` and `SimPhiN2()` functions.
* Transition matrices from `SimBetaN()` can now be generated with optional bounds.
* Drift matrices from `SimPhiN()` can now be generated with optional bounds.

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
