% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simStateSpace-sim-ssm-i-vary.R
\name{SimSSMIVary}
\alias{SimSSMIVary}
\title{Simulate Data using a State Space Model Parameterization
for n > 1 Individuals (Individual-Varying Parameters)}
\usage{
SimSSMIVary(
  n,
  mu0,
  sigma0_sqrt,
  alpha,
  beta,
  psi_sqrt,
  nu,
  lambda,
  theta_sqrt,
  gamma_y = NULL,
  gamma_eta = NULL,
  x = NULL,
  type,
  time = 0,
  burn_in = 0
)
}
\arguments{
\item{n}{Positive integer.
Number of individuals.}

\item{mu0}{List of numeric vectors.
Each element of the list
is the mean of initial latent variable values
(\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).}

\item{sigma0_sqrt}{List of numeric matrices.
Each element of the list
is the Cholesky decomposition of the covariance matrix
of initial latent variable values
(\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).}

\item{alpha}{List of numeric vectors.
Each element of the list
is the vector of intercepts for the dynamic model
(\eqn{\boldsymbol{\alpha}}).}

\item{beta}{List of numeric matrices.
Each element of the list
is the transition matrix relating the values of the latent variables
at time \code{t - 1} to those at time \code{t}
(\eqn{\boldsymbol{\beta}}).}

\item{psi_sqrt}{List of numeric matrices.
Each element of the list
is the Cholesky decomposition of the process noise covariance matrix
(\eqn{\boldsymbol{\Psi}}).}

\item{nu}{List of numeric vectors.
Each element of the list
is the vector of intercepts for the measurement model
(\eqn{\boldsymbol{\nu}}).}

\item{lambda}{List of numeric matrices.
Each element of the list
is the factor loading matrix linking the latent variables
to the observed variables
(\eqn{\boldsymbol{\Lambda}}).}

\item{theta_sqrt}{List of numeric matrices.
Each element of the list
is the Cholesky decomposition of the measurement error covariance matrix
(\eqn{\boldsymbol{\Theta}}).}

\item{gamma_y}{Numeric matrix.
Matrix relating the values of the covariate matrix
at time \code{t} to the observed variables at time \code{t}
(\eqn{\boldsymbol{\Gamma}_{\mathbf{y}}}).}

\item{gamma_eta}{Numeric matrix.
Matrix relating the values of the covariate matrix
at time \code{t} to the latent variables at time \code{t}
(\eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}}).}

\item{x}{A list of length \code{n} of numeric matrices.
Each element of the list
is a matrix of observed covariates in \code{type = 1} or \code{type = 2}.
The number of rows in each matrix should be equal to \code{time + burn_in}.}

\item{type}{Integer.
State space model type.}

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
\item \code{x}: A \code{t} by \code{j} matrix of values for the covariates.
\item \code{time}: A vector of discrete time points from 1 to \code{t}.
\item \code{id}: A vector of ID numbers of length \code{t}.
}
}
\description{
This function simulates data
using a state space model parameterization
for \code{n > 1} individuals.
In this model,
the parameters can vary across individuals.
}
\details{
Parameters can vary across individuals
by providing a list of parameter values.
If the length of any of the parameters
(\code{mu0},
\code{sigma0_sqrt},
\code{alpha},
\code{beta},
\code{psi_sqrt},
\code{nu},
\code{lambda},
\code{theta_sqrt},
\code{gamma_y}, or
\code{gamma_eta})
is less the \code{n},
the function will cycle through the available values.
}
\examples{
# prepare parameters
# In this example, beta varies across individuals
set.seed(42)
k <- p <- 3
iden <- diag(k)
iden_sqrt <- chol(iden)
null_vec <- rep(x = 0, times = k)
n <- 5
mu0 <- list(null_vec)
sigma0_sqrt <- list(iden_sqrt)
alpha <- list(null_vec)
beta <- list(
  diag(x = 0.1, nrow = k),
  diag(x = 0.2, nrow = k),
  diag(x = 0.3, nrow = k),
  diag(x = 0.4, nrow = k),
  diag(x = 0.5, nrow = k)
)
psi_sqrt <- list(iden_sqrt)
nu <- list(null_vec)
lambda <- list(iden)
theta_sqrt <- list(chol(diag(x = 0.50, nrow = k)))
time <- 50
burn_in <- 0
gamma_y <- gamma_eta <- list(0.10 * diag(k))
x <- lapply(
  X = seq_len(n),
  FUN = function(i) {
    return(
      matrix(
        data = rnorm(n = k * (time + burn_in)),
        ncol = k
      )
    )
  }
)

# Type 0
ssm <- SimSSMIVary(
  n = n,
  mu0 = mu0,
  sigma0_sqrt = sigma0_sqrt,
  alpha = alpha,
  beta = beta,
  psi_sqrt = psi_sqrt,
  nu = nu,
  lambda = lambda,
  theta_sqrt = theta_sqrt,
  type = 0,
  time = time,
  burn_in = burn_in
)

str(ssm)

# Type 1
ssm <- SimSSMIVary(
  n = n,
  mu0 = mu0,
  sigma0_sqrt = sigma0_sqrt,
  alpha = alpha,
  beta = beta,
  psi_sqrt = psi_sqrt,
  nu = nu,
  lambda = lambda,
  theta_sqrt = theta_sqrt,
  gamma_eta = gamma_eta,
  x = x,
  type = 1,
  time = time,
  burn_in = burn_in
)

str(ssm)

# Type 2
ssm <- SimSSMIVary(
  n = n,
  mu0 = mu0,
  sigma0_sqrt = sigma0_sqrt,
  alpha = alpha,
  beta = beta,
  psi_sqrt = psi_sqrt,
  nu = nu,
  lambda = lambda,
  theta_sqrt = theta_sqrt,
  gamma_y = gamma_y,
  gamma_eta = gamma_eta,
  x = x,
  type = 2,
  time = time,
  burn_in = burn_in
)

str(ssm)

}
\references{
Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling
and state-space modeling techniques.
\emph{Structural Equation Modeling: A Multidisciplinary Journal},
17(2), 303–332.
\doi{10.1080/10705511003661553}
}
\seealso{
Other Simulation of State Space Models Data Functions: 
\code{\link{OU2SSM}()},
\code{\link{Sim2Matrix}()},
\code{\link{SimSSMFixed}()},
\code{\link{SimSSMLinGrowthIVary}()},
\code{\link{SimSSMLinGrowth}()},
\code{\link{SimSSMOUFixed}()},
\code{\link{SimSSMOUIVary}()},
\code{\link{SimSSMOU}()},
\code{\link{SimSSMVARFixed}()},
\code{\link{SimSSMVARIVary}()},
\code{\link{SimSSMVAR}()},
\code{\link{SimSSM}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{Simulation of State Space Models Data Functions}
\keyword{sim}
\keyword{simStateSpace}
\keyword{ssm}
