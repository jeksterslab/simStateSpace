% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simStateSpace-sim-ssm.R
\name{SimSSM}
\alias{SimSSM}
\title{Simulate Data from a State Space Model (n = 1)}
\usage{
SimSSM(
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
  type = 0,
  time,
  burn_in = 0
)
}
\arguments{
\item{mu0}{Numeric vector.
Mean of initial latent variable values
(\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).}

\item{sigma0_sqrt}{Numeric matrix.
Cholesky decomposition of the covariance matrix
of initial latent variable values
(\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).}

\item{alpha}{Numeric vector.
Vector of intercepts for the dynamic model
(\eqn{\boldsymbol{\alpha}}).}

\item{beta}{Numeric matrix.
Transition matrix relating the values of the latent variables
at time \code{t - 1} to those at time \code{t}
(\eqn{\boldsymbol{\beta}}).}

\item{psi_sqrt}{Numeric matrix.
Cholesky decomposition of the process noise covariance matrix
(\eqn{\boldsymbol{\Psi}}).}

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

\item{gamma_y}{Numeric matrix.
Matrix relating the values of the covariate matrix
at time \code{t} to the observed variables at time \code{t}
(\eqn{\boldsymbol{\Gamma}_{\mathbf{y}}}).}

\item{gamma_eta}{Numeric matrix.
Matrix relating the values of the covariate matrix
at time \code{t} to the latent variables at time \code{t}
(\eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}}).}

\item{x}{Numeric matrix.
The matrix of observed covariates in \code{type = 1} or \code{type = 2}.
The number of rows should be equal to \code{time + burn_in}.}

\item{type}{Integer.
State space model type.}

\item{time}{Positive integer.
Number of time points to simulate.}

\item{burn_in}{Positive integer.
Number of burn-in points to exclude before returning the results.}
}
\value{
Returns a list with the following elements:
\itemize{
\item \code{y}: A \code{t} by \code{k} matrix of values for the manifest variables.
\item \code{eta}: A \code{t} by \code{p} matrix of values for the latent variables.
\item \code{x}: A \code{t} by \code{j} matrix of values for the covariates.
\item \code{time}: A vector of discrete time points from 0 to \code{t - 1}.
\item \code{id}: A vector of ones.
}
}
\description{
This function simulates data from a state space model.
See details for more information.
}
\details{
\subsection{Type 0}{

The measurement model is given by
\deqn{
    \mathbf{y}_{t}
    =
    \boldsymbol{\nu}
    +
    \boldsymbol{\Lambda}
    \boldsymbol{\eta}_{t}
    +
    \boldsymbol{\varepsilon}_{t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\varepsilon}_{t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Theta}
    \right)
  }
where
\eqn{\mathbf{y}_{t}},
\eqn{\boldsymbol{\eta}_{t}},
and
\eqn{\boldsymbol{\varepsilon}_{t}}
are random variables
and
\eqn{\boldsymbol{\nu}},
\eqn{\boldsymbol{\Lambda}},
and
\eqn{\boldsymbol{\Theta}}
are model parameters.
\eqn{\mathbf{y}_{t}}
is a vector of observed random variables,
\eqn{\boldsymbol{\eta}_{t}}
is a vector of latent random variables,
and
\eqn{\boldsymbol{\varepsilon}_{t}}
is a vector of random measurement errors,
at time \eqn{t}.
\eqn{\boldsymbol{\nu}}
is a vector of intercepts,
\eqn{\boldsymbol{\Lambda}}
is a matrix of factor loadings,
and
\eqn{\boldsymbol{\Theta}}
is the covariance matrix of
\eqn{\boldsymbol{\varepsilon}}.

The dynamic structure is given by
\deqn{
    \boldsymbol{\eta}_{t}
    =
    \boldsymbol{\alpha}
    +
    \boldsymbol{\beta}
    \boldsymbol{\eta}_{t - 1}
    +
    \boldsymbol{\zeta}_{t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\zeta}_{t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Psi}
    \right)
  }
where
\eqn{\boldsymbol{\eta}_{t}},
\eqn{\boldsymbol{\eta}_{t - 1}},
and
\eqn{\boldsymbol{\zeta}_{t}}
are random variables,
and
\eqn{\boldsymbol{\alpha}},
\eqn{\boldsymbol{\beta}},
and
\eqn{\boldsymbol{\Psi}}
are model parameters.
\eqn{\boldsymbol{\eta}_{t}}
is a vector of latent variables
at time \eqn{t},
\eqn{\boldsymbol{\eta}_{t - 1}}
is a vector of latent variables
at time \eqn{t - 1},
and
\eqn{\boldsymbol{\zeta}_{t}}
is a vector of dynamic noise
at time \eqn{t}.
\eqn{\boldsymbol{\alpha}}
is a vector of intercepts,
\eqn{\boldsymbol{\beta}}
is a matrix of autoregression
and cross regression coefficients,
and
\eqn{\boldsymbol{\Psi}}
is the covariance matrix of
\eqn{\boldsymbol{\zeta}_{t}}.
}

\subsection{Type 1}{

The measurement model is given by
\deqn{
    \mathbf{y}_{t}
    =
    \boldsymbol{\nu}
    +
    \boldsymbol{\Lambda}
    \boldsymbol{\eta}_{t}
    +
    \boldsymbol{\varepsilon}_{t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\varepsilon}_{t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Theta}
    \right) .
  }

The dynamic structure is given by
\deqn{
    \boldsymbol{\eta}_{t}
    =
    \boldsymbol{\alpha}
    +
    \boldsymbol{\beta}
    \boldsymbol{\eta}_{t - 1}
    +
    \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
    \mathbf{x}_{t}
    +
    \boldsymbol{\zeta}_{t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\zeta}_{t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Psi}
    \right)
  }
where
\eqn{\mathbf{x}_{t}}
is a vector of covariates
at time \eqn{t},
and
\eqn{\boldsymbol{\Gamma}_{\boldsymbol{\eta}}}
is the coefficient matrix
linking the covariates to the latent variables.
}

\subsection{Type 2}{

The measurement model is given by
\deqn{
    \mathbf{y}_{t}
    =
    \boldsymbol{\nu}
    +
    \boldsymbol{\Lambda}
    \boldsymbol{\eta}_{t}
    +
    \boldsymbol{\Gamma}_{\mathbf{y}}
    \mathbf{x}_{t}
    +
    \boldsymbol{\varepsilon}_{t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\varepsilon}_{t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Theta}
    \right)
  }
where
\eqn{\boldsymbol{\Gamma}_{\mathbf{y}}}
is the coefficient matrix
linking the covariates to the observed variables.

The dynamic structure is given by
\deqn{
    \boldsymbol{\eta}_{t}
    =
    \boldsymbol{\alpha}
    +
    \boldsymbol{\beta}
    \boldsymbol{\eta}_{t - 1}
    +
    \boldsymbol{\Gamma}_{\boldsymbol{\eta}}
    \mathbf{x}_{t}
    +
    \boldsymbol{\zeta}_{t}
    \quad
    \mathrm{with}
    \quad
    \boldsymbol{\zeta}_{t}
    \sim
    \mathcal{N}
    \left(
    \mathbf{0},
    \boldsymbol{\Psi}
    \right) .
  }
}
}
\examples{
# prepare parameters
set.seed(42)
k <- p <- 3
iden <- diag(k)
iden_sqrt <- chol(iden)
null_vec <- rep(x = 0, times = k)
mu0 <- null_vec
sigma0_sqrt <- iden_sqrt
alpha <- null_vec
beta <- diag(x = 0.50, nrow = k)
psi_sqrt <- iden_sqrt
nu <- null_vec
lambda <- iden
theta_sqrt <- chol(diag(x = 0.50, nrow = k))
time <- 50
burn_in <- 0
gamma_y <- gamma_eta <- 0.10 * diag(k)
x <- matrix(
  data = rnorm(n = k * (time + burn_in)),
  ncol = k
)

# Type 0
ssm <- SimSSM(
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
ssm <- SimSSM(
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
ssm <- SimSSM(
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
\code{\link{SimSSMIVary}()},
\code{\link{SimSSMLinGrowthIVary}()},
\code{\link{SimSSMLinGrowth}()},
\code{\link{SimSSMOUFixed}()},
\code{\link{SimSSMOUIVary}()},
\code{\link{SimSSMOU}()},
\code{\link{SimSSMVARFixed}()},
\code{\link{SimSSMVARIVary}()},
\code{\link{SimSSMVAR}()}
}
\concept{Simulation of State Space Models Data Functions}
\keyword{sim}
\keyword{simStateSpace}
\keyword{ssm}
