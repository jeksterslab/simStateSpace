% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SimSSMOUFixed}
\alias{SimSSMOUFixed}
\title{Simulate Data from an Ornstein–Uhlenbeck Model
using a State Space Model Parameterization
for n > 1 Individuals (Fixed Parameters)}
\usage{
SimSSMOUFixed(
  n,
  mu0,
  sigma0_sqrt,
  mu,
  phi,
  sigma_sqrt,
  nu,
  lambda,
  theta_sqrt,
  delta_t,
  time,
  burn_in
)
}
\arguments{
\item{n}{Positive integer.
Number of individuals.}

\item{mu0}{Numeric vector.
Mean of initial latent variable values
(\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).}

\item{sigma0_sqrt}{Numeric matrix.
Cholesky decomposition of the covariance matrix
of initial latent variable values
(\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).}

\item{mu}{Numeric vector.
The long-term mean or equilibrium level
(\eqn{\boldsymbol{\mu}}).}

\item{phi}{Numeric matrix.
The rate of mean reversion,
determining how quickly the variable returns to its mean
(\eqn{\boldsymbol{\Phi}}).}

\item{sigma_sqrt}{Numeric matrix.
Cholesky decomposition of the matrix of volatility
or randomness in the process
(\eqn{\boldsymbol{\Sigma}}).}

\item{nu}{Numeric vector.
Vector of intercepts for the measurement model
(\eqn{\boldsymbol{\nu}}).}

\item{lambda}{Numeric matrix.
Factor loading matrix linking the latent variables
to the observed variables
(\eqn{\boldsymbol{\Lambda}}).}

\item{theta_sqrt}{Numeric matrix.
Cholesky decomposition of the measurement error covariance matrix
(\eqn{\boldsymbol{\Theta}}).}

\item{delta_t}{Numeric.
Time interval (\eqn{\delta_t}).}

\item{time}{Positive integer.
Number of time points to simulate.}

\item{burn_in}{Positive integer.
Number of burn-in points to exclude before returning the results.}
}
\value{
Returns a list of length \code{n}.
Each element is a list with the following elements:
\itemize{
\item \code{y}: A \code{t} by \code{k} matrix of values for the manifest variables.
\item \code{eta}: A \code{t} by \code{p} matrix of values for the latent variables.
\item \code{time}: A vector of continuous time points of length \code{t}
starting from 0 with \code{delta_t} increments.
\item \code{id}: A vector of ID numbers of length \code{t}.
\item \code{n}: Number of individuals.
}
}
\description{
This function simulates data from an Ornstein–Uhlenbeck model
using a state space model parameterization
for \code{n > 1} individuals.
In this model,
the parameters are invariant across individuals.
See details for more information.
}
\details{
The measurement model is given by
\deqn{
    \mathbf{y}_{i, t}
    =
    \boldsymbol{\nu}
    +
    \boldsymbol{\Lambda}
    \boldsymbol{\eta}_{i, t}
    +
    \boldsymbol{\varepsilon}_{i, t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\varepsilon}_{i, t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Theta}
    \right)
  }
where \eqn{\mathbf{y}_{i, t}}, \eqn{\boldsymbol{\eta}_{i, t}},
and \eqn{\boldsymbol{\varepsilon}_{i, t}}
are random variables and \eqn{\boldsymbol{\nu}},
\eqn{\boldsymbol{\Lambda}},
and \eqn{\boldsymbol{\Theta}} are model parameters.
\eqn{\mathbf{y}_{i, t}} is a vector of observed random variables
at time \eqn{t} and individual \eqn{i},
\eqn{\boldsymbol{\eta}_{i, t}} is a vector of latent random variables
at time \eqn{t} and individual \eqn{i},
and \eqn{\boldsymbol{\varepsilon}_{i, t}}
is a vector of random measurement errors
at time \eqn{t} and individual \eqn{i},
while \eqn{\boldsymbol{\nu}} is a vector of intercept,
\eqn{\boldsymbol{\Lambda}} is a matrix of factor loadings,
and \eqn{\boldsymbol{\Theta}} is the covariance matrix of
\eqn{\boldsymbol{\varepsilon}}.

The dynamic structure is given by
\deqn{
    \mathrm{d} \boldsymbol{\eta}_{i, t}
    =
    \boldsymbol{\Phi}
    \left(
    \boldsymbol{\mu}
    -
    \boldsymbol{\eta}_{i, t}
    \right)
    \mathrm{d}t
    +
    \boldsymbol{\Sigma}^{\frac{1}{2}}
    \mathrm{d}
    \mathbf{W}_{i, t}
  }
where \eqn{\boldsymbol{\mu}} is the long-term mean or equilibrium level,
\eqn{\boldsymbol{\Phi}} is the rate of mean reversion,
determining how quickly the variable returns to its mean,
\eqn{\boldsymbol{\Sigma}} is the matrix of volatility
or randomness in the process, and \eqn{\mathrm{d}\boldsymbol{W}}
is a Wiener process or Brownian motion,
which represents random fluctuations.
}
\examples{
# prepare parameters
set.seed(42)
p <- k <- 2
I <- diag(p)
I_sqrt <- chol(I)
n <- 5
mu0 <- c(-3.0, 1.5)
sigma0_sqrt <- I_sqrt
mu <- c(5.76, 5.18)
phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
sigma_sqrt <- chol(
  matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
)
nu <- rep(x = 0, times = k)
lambda <- diag(k)
theta_sqrt <- chol(diag(x = 0.50, nrow = k))
delta_t <- 0.10
time <- 50
burn_in <- 0

# generate data
ssm <- SimSSMOUFixed(
  n = n,
  mu0 = mu0,
  sigma0_sqrt = sigma0_sqrt,
  mu = mu,
  phi = phi,
  sigma_sqrt = sigma_sqrt,
  nu = nu,
  lambda = lambda,
  theta_sqrt = theta_sqrt,
  delta_t = delta_t,
  time = time,
  burn_in = burn_in
)

str(ssm)

}
\references{
Uhlenbeck, G. E., & Ornstein, L. S. (1930).
On the theory of the brownian motion.
\emph{Physical Review}, \emph{36}(5), 823–841.
\doi{10.1103/physrev.36.823}
}
\seealso{
Other Simulation of State Space Models Data Functions: 
\code{\link{OU2SSM}()},
\code{\link{Sim2Matrix}()},
\code{\link{SimSSM0Fixed}()},
\code{\link{SimSSM0Vary}()},
\code{\link{SimSSM0}()},
\code{\link{SimSSMOUVary}()},
\code{\link{SimSSMOU}()},
\code{\link{SimSSMVARFixed}()},
\code{\link{SimSSMVARVary}()},
\code{\link{SimSSMVAR}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{Simulation of State Space Models Data Functions}
\keyword{sim}
\keyword{simStateSpace}