% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simStateSpace-lin-sde-mean-eta.R
\name{LinSDEMeanEta}
\alias{LinSDEMeanEta}
\title{Steady-State Mean Vector for the
Latent Variables in the
Linear Stochastic Differential Equation Model}
\usage{
LinSDEMeanEta(phi, iota)
}
\arguments{
\item{phi}{Numeric matrix.
The drift matrix
which represents the rate of change of the solution
in the absence of any random fluctuations
(\eqn{\boldsymbol{\Phi}}).}

\item{iota}{Numeric vector.
An unobserved term that is constant over time
(\eqn{\boldsymbol{\iota}}).}
}
\description{
The steady-state mean vector
for the latent variables
in the linear stochastic differential equation model
\eqn{\mathrm{Mean} \left( \boldsymbol{\eta} \right)}
is given by
\deqn{
  \mathrm{Mean} \left( \boldsymbol{\eta} \right)
  =
  -\boldsymbol{\Phi}^{-1} \boldsymbol{\iota}
}
where
\eqn{\boldsymbol{\Phi}}
is the drift matrix,
and
\eqn{\boldsymbol{\iota}}
is an unobserved term that is constant over time.
}
\examples{
phi <- matrix(
  data = c(
    -0.357, 0.771, -0.450,
    0.0, -0.511, 0.729,
    0.0, 0.0, -0.693
  ),
  nrow = 3
)
iota <- rep(x = 1, times = 3)
LinSDEMeanEta(
  phi = phi,
  iota = iota
)

}
\seealso{
Other Simulation of State Space Models Data Functions: 
\code{\link{LinSDE2SSM}()},
\code{\link{LinSDECovEta}()},
\code{\link{LinSDECovY}()},
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
\code{\link{SimSSMLinSDEIVary}()},
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
\keyword{simStateSpace}
