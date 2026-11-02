# The Latent Curve Model with Structured Residuals

## Overview

The latent curve model with structured residuals (LCM-SR) combines two
ideas (Curran et al., 2014). Take alcohol use and depressive symptoms
measured once a year in adolescents. Each adolescent follows a
straight-line trajectory on each variable, with their own intercept and
slope, as in the linear growth curve model. On each occasion, the
adolescent is also above or below the trajectory, and these deviations,
called the structured residuals, carry over to the next occasion and
affect each other. The growth part describes how individuals differ in
their long-term change; the residual part describes the within-person
dynamics around it, with autoregressive and cross-lagged effects.

The LCM-SR is usually fitted as a structural equation model (SEM). It
can also be written as a state space model by combining the state space
forms of the linear growth curve model and of the vector autoregressive
model: the levels and slopes of the two trajectories and the two
residuals are latent states.

This tutorial does three things:

1.  It writes the LCM-SR as a state space model.
2.  It simulates data from the model with
    [`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md).
3.  It fits the model to the simulated data as a state space model with
    `dynr` and `OpenMx` and as a structural equation model (SEM) with
    `lavaan` and `OpenMx`, and compares every set of estimates with the
    values used to simulate the data.

No prior knowledge of state space models is needed. The vignette on the
linear growth curve model explains the growth part in more detail.

## Model

The measurement model is given by \\\[\left( \begin{array}{c} Y\_{1\_{i,
t}} \\\\ Y\_{2\_{i, t}} \\\\ \end{array} \right) = \left(
\begin{array}{cccccc} 1 & 0 & 0 & 0 & 1 & 0 \\\\ 0 & 1 & 0 & 0 & 0 & 1
\\\\ \end{array} \right) \left( \begin{array}{c} \eta\_{1\_{i, t}} \\\\
\eta\_{2\_{i, t}} \\\\ \eta\_{3\_{i, t}} \\\\ \eta\_{4\_{i, t}} \\\\
\eta\_{5\_{i, t}} \\\\ \eta\_{6\_{i, t}} \\\\ \end{array} \right)\\\]
where \\(Y\_{1\_{i, t}}\\) and \\(Y\_{2\_{i, t}}\\) are the two observed
variables at time \\(t\\) and individual \\(i\\), \\(\eta\_{1\_{i,
t}}\\) and \\(\eta\_{2\_{i, t}}\\) are the current levels of the two
trajectories, \\(\eta\_{3\_{i, t}}\\) and \\(\eta\_{4\_{i, t}}\\) are
their slopes, and \\(\eta\_{5\_{i, t}}\\) and \\(\eta\_{6\_{i, t}}\\)
are the structured residuals.

In words, each observed score is the point on the individual’s
trajectory plus the individual’s deviation from it on that occasion. As
in most applications of the LCM-SR, there is no separate measurement
error.

The dynamic structure is given by \\\[\left( \begin{array}{c}
\eta\_{1\_{i, t}} \\\\ \eta\_{2\_{i, t}} \\\\ \eta\_{3\_{i, t}} \\\\
\eta\_{4\_{i, t}} \\\\ \eta\_{5\_{i, t}} \\\\ \eta\_{6\_{i, t}} \\\\
\end{array} \right) = \left( \begin{array}{cccccc} 1 & 0 & 1 & 0 & 0 & 0
\\\\ 0 & 1 & 0 & 1 & 0 & 0 \\\\ 0 & 0 & 1 & 0 & 0 & 0 \\\\ 0 & 0 & 0 & 1
& 0 & 0 \\\\ 0 & 0 & 0 & 0 & \beta\_{55} & \beta\_{56} \\\\ 0 & 0 & 0 &
0 & \beta\_{65} & \beta\_{66} \\\\ \end{array} \right) \left(
\begin{array}{c} \eta\_{1\_{i, t - 1}} \\\\ \eta\_{2\_{i, t - 1}} \\\\
\eta\_{3\_{i, t - 1}} \\\\ \eta\_{4\_{i, t - 1}} \\\\ \eta\_{5\_{i, t -
1}} \\\\ \eta\_{6\_{i, t - 1}} \\\\ \end{array} \right) + \left(
\begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\ 0 \\\\ \zeta\_{5\_{i, t}} \\\\
\zeta\_{6\_{i, t}} \\\\ \end{array} \right), \quad \mathrm{with} \quad
\left( \begin{array}{c} \zeta\_{5\_{i, t}} \\\\ \zeta\_{6\_{i, t}} \\\\
\end{array} \right) \sim \mathcal{N} \left( \mathbf{0}, \left(
\begin{array}{cc} \psi\_{55} & \psi\_{56} \\\\ \psi\_{65} & \psi\_{66}
\\\\ \end{array} \right) \right) .\\\]

In words, the first four rows are the linear growth curve model for each
variable: the level at each occasion is the previous level plus the
slope, and the slope stays the same, with no process noise. The last two
rows are a vector autoregressive model for the residuals:
\\(\beta\_{55}\\) and \\(\beta\_{66}\\) are the autoregressive effects,
how much of a deviation from the trajectory carries over to the next
occasion, and \\(\beta\_{56}\\) and \\(\beta\_{65}\\) are the
cross-lagged effects: \\(\beta\_{56}\\) is the effect of a deviation in
the second variable on the deviation in the first variable at the next
occasion, and \\(\beta\_{65}\\) the reverse. \\(\psi\_{55}\\),
\\(\psi\_{66}\\), and \\(\psi\_{56}\\) are the variances and covariance
of the new part of the deviations.

The initial condition holds the means, variances, and covariances of the
intercepts and slopes (the levels at the first occasion are the
intercepts), and the variances and covariance of the residuals at the
first occasion. The residuals have mean zero, because they are
deviations from the trajectories, and are uncorrelated with the
intercepts and slopes.

### Diagrams

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, *17*(2), 303–332.
<https://doi.org/10.1080/10705511003661553>

Curran, P. J., Howard, A. L., Bainter, S. A., Lane, S. T., & McGinley,
J. S. (2014). The separation of between-person and within-person
components of individual change over time: A latent curve model with
structured residuals. *Journal of Consulting and Clinical Psychology*,
*82*(5), 879–894. <https://doi.org/10.1037/a0035297>

Ou, L., Hunter, M. D., & Chow, S.-M. (2019). What’s for dynr: A package
for linear and nonlinear dynamic modeling in R. *The R Journal*,
*11*(1), 91. <https://doi.org/10.32614/rj-2019-012>

R Core Team. (2026). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software*, *48*(2).
<https://doi.org/10.18637/jss.v048.i02>
