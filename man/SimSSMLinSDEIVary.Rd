% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simStateSpace-sim-ssm-lin-sde-i-vary.R
\name{SimSSMLinSDEIVary}
\alias{SimSSMLinSDEIVary}
\title{Simulate Data from the
Linear Stochastic Differential Equation Model
using a State Space Model Parameterization
(Individual-Varying Parameters)}
\usage{
SimSSMLinSDEIVary(
  n,
  time,
  delta_t = 1,
  mu0,
  sigma0_l,
  iota,
  phi,
  sigma_l,
  nu,
  lambda,
  theta_l,
  type = 0,
  x = NULL,
  gamma = NULL,
  kappa = NULL
)
}
\arguments{
\item{n}{Positive integer.
Number of individuals.}

\item{time}{Positive integer.
Number of time points.}

\item{delta_t}{Numeric.
Time interval.
The default value is \code{1.0}
with an option to use a numeric value
for the discretized state space model
parameterization of the
linear stochastic differential equation model.}

\item{mu0}{List of numeric vectors.
Each element of the list
is the mean of initial latent variable values
(\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).}

\item{sigma0_l}{List of numeric matrices.
Each element of the list
is the Cholesky factorization (\code{t(chol(sigma0))})
of the covariance matrix
of initial latent variable values
(\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).}

\item{iota}{List of numeric vectors.
Each element of the list
is an unobserved term that is constant over time
(\eqn{\boldsymbol{\iota}}).}

\item{phi}{List of numeric matrix.
Each element of the list
is the drift matrix
which represents the rate of change of the solution
in the absence of any random fluctuations
(\eqn{\boldsymbol{\Phi}}).}

\item{sigma_l}{List of numeric matrix.
Each element of the list
is the Cholesky factorization (\code{t(chol(sigma))})
of the covariance matrix of volatility
or randomness in the process
\eqn{\boldsymbol{\Sigma}}.}

\item{nu}{List of numeric vectors.
Each element of the list
is the vector of intercept values for the measurement model
(\eqn{\boldsymbol{\nu}}).}

\item{lambda}{List of numeric matrices.
Each element of the list
is the factor loading matrix linking the latent variables
to the observed variables
(\eqn{\boldsymbol{\Lambda}}).}

\item{theta_l}{List of numeric matrices.
Each element of the list
is the Cholesky factorization (\code{t(chol(theta))})
of the covariance matrix
of the measurement error
(\eqn{\boldsymbol{\Theta}}).}

\item{type}{Integer.
State space model type.
See Details in \code{\link[=SimSSMLinSDEFixed]{SimSSMLinSDEFixed()}} for more information.}

\item{x}{List.
Each element of the list is a matrix of covariates
for each individual \code{i} in \code{n}.
The number of columns in each matrix
should be equal to \code{time}.}

\item{gamma}{List of numeric matrices.
Each element of the list
is the matrix linking the covariates to the latent variables
at current time point
(\eqn{\boldsymbol{\Gamma}}).}

\item{kappa}{List of numeric matrices.
Each element of the list
is the matrix linking the covariates to the observed variables
at current time point
(\eqn{\boldsymbol{\kappa}}).}
}
\value{
Returns an object of class \code{simstatespace}
which is a list with the following elements:
\itemize{
\item \code{call}: Function call.
\item \code{args}: Function arguments.
\item \code{data}: Generated data which is a list of length \code{n}.
Each element of \code{data} is a list with the following elements:
\itemize{
\item \code{id}: A vector of ID numbers with length \code{l},
where \code{l} is the value of the function argument \code{time}.
\item \code{time}: A vector time points of length \code{l}.
\item \code{y}: A \code{l} by \code{k} matrix of values for the manifest variables.
\item \code{eta}: A \code{l} by \code{p} matrix of values for the latent variables.
\item \code{x}: A \code{l} by \code{j} matrix of values for the covariates
(when covariates are included).
}
\item \code{fun}: Function used.
}
}
\description{
This function simulates data from the
linear stochastic differential equation model
using a state space model parameterization.
It assumes that the parameters can vary
across individuals.
}
\details{
Parameters can vary across individuals
by providing a list of parameter values.
If the length of any of the parameters
(\code{mu0},
\code{sigma0_l},
\code{iota},
\code{phi},
\code{sigma_l},
\code{nu},
\code{lambda},
\code{theta_l},
\code{gamma}, or
\code{kappa})
is less the \code{n},
the function will cycle through the available values.
}
\examples{
# prepare parameters
# In this example, phi varies across individuals.
set.seed(42)
## number of individuals
n <- 5
## time points
time <- 50
delta_t <- 0.10
## dynamic structure
p <- 2
mu0 <- list(
  c(-3.0, 1.5)
)
sigma0 <- 0.001 * diag(p)
sigma0_l <- list(
  t(chol(sigma0))
)
iota <- list(
  c(0.317, 0.230)
)
phi <- list(
  -0.1 * diag(p),
  -0.2 * diag(p),
  -0.3 * diag(p),
  -0.4 * diag(p),
  -0.5 * diag(p)
)
sigma <- matrix(
  data = c(
    2.79,
    0.06,
    0.06,
    3.27
  ),
  nrow = p
)
sigma_l <- list(
  t(chol(sigma))
)
## measurement model
k <- 2
nu <- list(
  rep(x = 0, times = k)
)
lambda <- list(
  diag(k)
)
theta <- 0.001 * diag(k)
theta_l <- list(
  t(chol(theta))
)
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
gamma <- list(
  diag(x = 0.10, nrow = p, ncol = j)
)
kappa <- list(
  diag(x = 0.10, nrow = k, ncol = j)
)

# Type 0
ssm <- SimSSMLinSDEIVary(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 0
)

plot(ssm)

# Type 1
ssm <- SimSSMLinSDEIVary(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 1,
  x = x,
  gamma = gamma
)

plot(ssm)

# Type 2
ssm <- SimSSMLinSDEIVary(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 2,
  x = x,
  gamma = gamma,
  kappa = kappa
)

plot(ssm)

}
\references{
Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling
and state-space modeling techniques.
\emph{Structural Equation Modeling: A Multidisciplinary Journal},
17(2), 303–332.
\doi{10.1080/10705511003661553}

Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
Continuous-time dynamic models:
Connections to structural equation models and other discrete-time models.
In R. H. Hoyle (Ed.),
Handbook of structural equation modeling (2nd ed.).
The Guilford Press.

Harvey, A. C. (1990).
Forecasting, structural time series models and the Kalman filter.
Cambridge University Press.
\doi{10.1017/cbo9781107049994}
}
\seealso{
Other Simulation of State Space Models Data Functions: 
\code{\link{LinSDE2SSM}()},
\code{\link{LinSDECovEta}()},
\code{\link{LinSDECovY}()},
\code{\link{LinSDEMeanEta}()},
\code{\link{LinSDEMeanY}()},
\code{\link{SSMCovEta}()},
\code{\link{SSMCovY}()},
\code{\link{SSMMeanEta}()},
\code{\link{SSMMeanY}()},
\code{\link{SimAlphaN}()},
\code{\link{SimBetaN}()},
\code{\link{SimCovDiagN}()},
\code{\link{SimCovN}()},
\code{\link{SimIotaN}()},
\code{\link{SimPhiN}()},
\code{\link{SimSSMFixed}()},
\code{\link{SimSSMIVary}()},
\code{\link{SimSSMLinGrowth}()},
\code{\link{SimSSMLinGrowthIVary}()},
\code{\link{SimSSMLinSDEFixed}()},
\code{\link{SimSSMOUFixed}()},
\code{\link{SimSSMOUIVary}()},
\code{\link{SimSSMVARFixed}()},
\code{\link{SimSSMVARIVary}()},
\code{\link{TestPhi}()},
\code{\link{TestStability}()},
\code{\link{TestStationarity}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{Simulation of State Space Models Data Functions}
\keyword{linsde}
\keyword{sim}
\keyword{simStateSpace}
