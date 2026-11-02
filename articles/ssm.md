# The State Space Model

## Overview

A state space model describes variables that we cannot observe directly
(latent states) and how they change over time, using repeated,
error-prone measurements of them. It has two parts: a measurement model
that links the observed variables to the latent states, and a dynamic
model that describes how the latent states carry over and influence each
other from one occasion to the next.

This tutorial does three things:

1.  It writes a state space model with a first-order vector
    autoregressive dynamic model.
2.  It simulates data from the model with
    [`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md),
    with every individual sharing the same parameters.
3.  It fits the model to the simulated data with the `dynr` and `OpenMx`
    packages and compares both sets of estimates with the values used to
    simulate the data.

No prior knowledge of state space models is needed. If the model without
measurement error is enough for you, see the vignette on the vector
autoregressive model.

## Model

The measurement model is given by \\\[\mathbf{y}\_{i, t} =
\boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\eta}\_{i, t} +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( \mathbf{0},
\boldsymbol{\Theta} \right)\\\] where \\(\mathbf{y}\_{i, t}\\),
\\(\boldsymbol{\eta}\_{i, t}\\), and \\(\boldsymbol{\varepsilon}\_{i,
t}\\) are random variables and \\(\boldsymbol{\nu}\\),
\\(\boldsymbol{\Lambda}\\), and \\(\boldsymbol{\Theta}\\) are model
parameters. \\(\mathbf{y}\_{i, t}\\) represents a vector of observed
random variables, \\(\boldsymbol{\eta}\_{i, t}\\) a vector of latent
random variables, and \\(\boldsymbol{\varepsilon}\_{i, t}\\) a vector of
random measurement errors, at time \\(t\\) and individual \\(i\\).
\\(\boldsymbol{\nu}\\) denotes a vector of intercepts,
\\(\boldsymbol{\Lambda}\\) a matrix of factor loadings, and
\\(\boldsymbol{\Theta}\\) the covariance matrix of
\\(\boldsymbol{\varepsilon}\\).

In words, each observed variable is an intercept plus a weighted
combination of the latent variables plus measurement error. The latent
variables \\(\boldsymbol{\eta}\_{i, t}\\) are the quantities we care
about (for example, a person’s true stress level), and the observed
variables \\(\mathbf{y}\_{i, t}\\) are noisy measurements of them.
Separating measurement error from true change is the main reason to use
a state space model: measurement error affects only the occasion it
occurs on, while changes in the latent variables carry over in time.

An alternative representation of the measurement error is given by
\\\[\boldsymbol{\varepsilon}\_{i, t} = \boldsymbol{\Theta}^{\frac{1}{2}}
\mathbf{z}\_{i, t}, \quad \mathrm{with} \quad \mathbf{z}\_{i, t} \sim
\mathcal{N} \left( \mathbf{0}, \mathbf{I} \right)\\\] where
\\(\mathbf{z}\_{i, t}\\) is a vector of independent standard normal
random variables and \\(\left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)
\left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)^{\prime} =
\boldsymbol{\Theta}\\) .

The dynamic structure is given by \\\[\boldsymbol{\eta}\_{i, t} =
\boldsymbol{\alpha} + \boldsymbol{\beta} \boldsymbol{\eta}\_{i, t - 1} +
\boldsymbol{\zeta}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\zeta}\_{i, t} \sim \mathcal{N} \left( \mathbf{0},
\boldsymbol{\Psi} \right)\\\] where \\(\boldsymbol{\eta}\_{i, t}\\),
\\(\boldsymbol{\eta}\_{i, t - 1}\\), and \\(\boldsymbol{\zeta}\_{i,
t}\\) are random variables, and \\(\boldsymbol{\alpha}\\),
\\(\boldsymbol{\beta}\\), and \\(\boldsymbol{\Psi}\\) are model
parameters. Here, \\(\boldsymbol{\eta}\_{i, t}\\) is a vector of latent
variables at time \\(t\\) and individual \\(i\\),
\\(\boldsymbol{\eta}\_{i, t - 1}\\) represents a vector of latent
variables at time \\(t - 1\\) and individual \\(i\\), and
\\(\boldsymbol{\zeta}\_{i, t}\\) represents a vector of dynamic noise at
time \\(t\\) and individual \\(i\\). \\(\boldsymbol{\alpha}\\) denotes a
vector of intercepts, \\(\boldsymbol{\beta}\\) a matrix of
autoregression and cross regression coefficients, and
\\(\boldsymbol{\Psi}\\) the covariance matrix of
\\(\boldsymbol{\zeta}\_{i, t}\\).

In words, each latent variable at time \\(t\\) is a weighted sum of all
the latent variables at time \\(t - 1\\), plus an intercept and a random
shock. The weights are the elements of \\(\boldsymbol{\beta}\\):

- A diagonal element \\(\beta\_{jj}\\) is the autoregressive effect, how
  much of variable \\(j\\) carries over to the next occasion. Values
  closer to 1 mean the variable changes slowly (high inertia).
- An off-diagonal element \\(\beta\_{jk}\\) is the cross-lagged effect
  of variable \\(k\\) at time \\(t - 1\\) on variable \\(j\\) at time
  \\(t\\). Note the order: the row is the outcome, the column is the
  predictor.

The shocks \\(\boldsymbol{\zeta}\_{i, t}\\) stand for everything that
moves the latent variables that the model does not explain. Unlike
measurement error, a shock carries over to later occasions through
\\(\boldsymbol{\beta}\\). For the process to settle around a stable
mean, all eigenvalues of \\(\boldsymbol{\beta}\\) must have modulus less
than 1.

An alternative representation of the dynamic noise is given by
\\\[\boldsymbol{\zeta}\_{i, t} = \boldsymbol{\Psi}^{\frac{1}{2}}
\mathbf{z}\_{i, t}, \quad \mathrm{with} \quad \mathbf{z}\_{i, t} \sim
\mathcal{N} \left( \mathbf{0}, \mathbf{I} \right)\\\] where \\(\left(
\boldsymbol{\Psi}^{\frac{1}{2}} \right) \left(
\boldsymbol{\Psi}^{\frac{1}{2}} \right)^{\prime} = \boldsymbol{\Psi}\\)
. This is how both kinds of noise are simulated: independent standard
normal draws are multiplied by a square root of the covariance matrix.
The functions in `simStateSpace` take this square root, the lower
triangular Cholesky factor, as the arguments `psi_l` (for
\\(\boldsymbol{\Psi}\\)) and `theta_l` (for \\(\boldsymbol{\Theta}\\)).

### Diagram

![State space diagram of the state space
model.](fig-vignettes-ssm-ssm.png)

State space diagram of the state space model.

The model in state space form, for the first three occasions (\\(t = 0,
1, 2\\)), with every element of the state drawn. Squares are observed
variables and circles are latent variables. Each latent variable
\\(\eta\_{j, t}\\) has one indicator \\(Y\_{j, t}\\) with a loading of
1, and measurement error \\(\varepsilon\_{j, t}\\) with variance
\\(\theta\_{jj}\\). The horizontal arrows are the autoregressive effects
and the diagonal arrows the cross-lagged effects, the elements of
\\(\boldsymbol{\beta}\\); only the nonzero elements in this example are
drawn (\\(\eta\_{1}\\) affects \\(\eta\_{2}\\) and \\(\eta\_{3}\\), and
\\(\eta\_{2}\\) affects \\(\eta\_{3}\\)). \\(\zeta\_{j, t}\\) is the
process noise, with variance \\(\psi\_{jj}\\), from \\(t = 1\\) on. A
loop on a variable is its variance and a curved two-headed arrow a
covariance. At \\(t = 0\\), the initial condition
\\(\boldsymbol{\eta}\_{0} \sim \mathcal{N} \left( \boldsymbol{\mu}\_{0},
\boldsymbol{\Sigma}\_{0} \right)\\) gives the variances
\\(\sigma\_{jj}\\) and covariances \\(\sigma\_{jk}\\) of the latent
variables. Below the diagram, all the latent variables at one occasion
form the state vector \\(\boldsymbol{\eta}\_{t}\\), and all the observed
variables form \\(\mathbf{y}\_{t}\\): \\(\boldsymbol{\beta}\\) carries
the state from one occasion to the next, \\(\boldsymbol{\Lambda}\\) maps
the state to the observed variables, \\(\boldsymbol{\Psi}\\) is the
covariance matrix of the process noise \\(\boldsymbol{\zeta}\_{t}\\),
and \\(\boldsymbol{\Theta}\\) the covariance matrix of the measurement
error \\(\boldsymbol{\varepsilon}\_{t}\\).
[`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md)
takes `beta`, `lambda`, `mu0`, and the Cholesky factors of
\\(\boldsymbol{\Psi}\\), \\(\boldsymbol{\Theta}\\), and
\\(\boldsymbol{\Sigma}\_{0}\\) (`psi_l`, `theta_l`, and `sigma0_l`).

## Data Generation

### Notation

Let \\(t = 200\\) be the number of time points and \\(n = 100\\) be the
number of individuals.

Let the measurement model intercept vector \\(\boldsymbol{\nu}\\) be
given by

\\\[\boldsymbol{\nu} = \left( \begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\
\end{array} \right) .\\\]

Let the factor loadings matrix \\(\boldsymbol{\Lambda}\\) be given by

\\\[\boldsymbol{\Lambda} = \left( \begin{array}{ccc} 1 & 0 & 0 \\\\ 0 &
1 & 0 \\\\ 0 & 0 & 1 \\\\ \end{array} \right) .\\\]

Let the measurement error covariance matrix \\(\boldsymbol{\Theta}\\) be
given by

\\\[\boldsymbol{\Theta} = \left( \begin{array}{ccc} 0.2 & 0 & 0 \\\\ 0 &
0.2 & 0 \\\\ 0 & 0 & 0.2 \\\\ \end{array} \right) .\\\]

With \\(\boldsymbol{\nu} = \mathbf{0}\\) and \\(\boldsymbol{\Lambda} =
\mathbf{I}\\), each observed variable measures one latent variable, and
\\(\boldsymbol{\Theta}\\) adds measurement error with variance 0.2 to
each.

Let the initial condition \\(\boldsymbol{\eta}\_{0}\\) be given by

\\\[\boldsymbol{\eta}\_{0} \sim \mathcal{N} \left(
\boldsymbol{\mu}\_{0}, \boldsymbol{\Sigma}\_{0} \right)\\\]

\\\[\boldsymbol{\mu}\_{0} = \left( \begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\
\end{array} \right)\\\]

\\\[\boldsymbol{\Sigma}\_{0} = \left( \begin{array}{ccc} 0.1960784 &
0.1183232 & 0.0298539 \\\\ 0.1183232 & 0.3437711 & 0.1381855 \\\\
0.0298539 & 0.1381855 & 0.2663828 \\\\ \end{array} \right) .\\\]

Let the constant vector \\(\boldsymbol{\alpha}\\) be given by

\\\[\boldsymbol{\alpha} = \left( \begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\
\end{array} \right) .\\\]

Let the transition matrix \\(\boldsymbol{\beta}\\) be given by

\\\[\boldsymbol{\beta} = \left( \begin{array}{ccc} 0.7 & 0 & 0 \\\\ 0.5
& 0.6 & 0 \\\\ -0.1 & 0.4 & 0.5 \\\\ \end{array} \right) .\\\]

Read the columns as predictors and the rows as outcomes. Each variable
carries over to the next occasion (autoregressive effects 0.7, 0.6, and
0.5). In addition, \\(\eta\_{1}\\) predicts later \\(\eta\_{2}\\) (0.5),
\\(\eta\_{2}\\) predicts later \\(\eta\_{3}\\) (0.4), and
\\(\eta\_{1}\\) has a small negative effect on later \\(\eta\_{3}\\)
(-0.1). Because \\(\boldsymbol{\beta}\\) is lower triangular, its
eigenvalues are its diagonal elements, all less than 1, so the process
is stable.

Let the dynamic process noise \\(\boldsymbol{\Psi}\\) be given by

\\\[\boldsymbol{\Psi} = \left( \begin{array}{ccc} 0.1 & 0 & 0 \\\\ 0 &
0.1 & 0 \\\\ 0 & 0 & 0.1 \\\\ \end{array} \right) .\\\]

The initial condition describes where each individual starts. Here, the
mean vector and covariance matrix of \\(\boldsymbol{\eta}\_{0}\\) are
the model-implied mean and covariance of the stable process, computed
with
[`SSMMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanEta.md)
and
[`SSMCovEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovEta.md),
so the data are stationary from the first time point.

### R Function Arguments

The simulation functions take covariance matrices as lower triangular
Cholesky factors (the arguments ending in `_l`).

\
`n`\
`#> [1] 100`\
`time`\
`#> [1] 200`\
`mu0`\
`#> [1] 0 0 0`\
`sigma0`\
`#>            [,1]      [,2]       [,3]`\
`#> [1,] 0.19607843 0.1183232 0.02985385`\
`#> [2,] 0.11832319 0.3437711 0.13818551`\
`#> [3,] 0.02985385 0.1381855 0.26638284`\
`sigma0_l`` ``# sigma0_l <- t(chol(sigma0))`\
`#>            [,1]      [,2]     [,3]`\
`#> [1,] 0.44280744 0.0000000 0.000000`\
`#> [2,] 0.26721139 0.5218900 0.000000`\
`#> [3,] 0.06741949 0.2302597 0.456966`\
`alpha`\
`#> [1] 0 0 0`\
`beta`\
`#>      [,1] [,2] [,3]`\
`#> [1,]  0.7  0.0  0.0`\
`#> [2,]  0.5  0.6  0.0`\
`#> [3,] -0.1  0.4  0.5`\
`psi`\
`#>      [,1] [,2] [,3]`\
`#> [1,]  0.1  0.0  0.0`\
`#> [2,]  0.0  0.1  0.0`\
`#> [3,]  0.0  0.0  0.1`\
`psi_l`` ``# psi_l <- t(chol(psi))`\
`#>           [,1]      [,2]      [,3]`\
`#> [1,] 0.3162278 0.0000000 0.0000000`\
`#> [2,] 0.0000000 0.3162278 0.0000000`\
`#> [3,] 0.0000000 0.0000000 0.3162278`\
`nu`\
`#> [1] 0 0 0`\
`lambda`\
`#>      [,1] [,2] [,3]`\
`#> [1,]    1    0    0`\
`#> [2,]    0    1    0`\
`#> [3,]    0    0    1`\
`theta`\
`#>      [,1] [,2] [,3]`\
`#> [1,]  0.2  0.0  0.0`\
`#> [2,]  0.0  0.2  0.0`\
`#> [3,]  0.0  0.0  0.2`\
`theta_l`` ``# theta_l <- t(chol(theta))`\
`#>           [,1]      [,2]      [,3]`\
`#> [1,] 0.4472136 0.0000000 0.0000000`\
`#> [2,] 0.0000000 0.4472136 0.0000000`\
`#> [3,] 0.0000000 0.0000000 0.4472136`

### Visualizing the Dynamics Without Process Noise (n = 5 with Different Initial Condition)

To see what \\(\boldsymbol{\beta}\\) alone does, we first simulate five
individuals with the process noise and the measurement error set to
zero. Each individual starts at a different point, and the trajectories
return to the mean at a speed set by \\(\boldsymbol{\beta}\\).

![](fig-vignettes-ssm-no-error-ssm-1.png)![](fig-vignettes-ssm-no-error-ssm-2.png)![](fig-vignettes-ssm-no-error-ssm-3.png)

### Using the `SimSSMFixed` Function from the `simStateSpace` Package to Simulate Data

Now we add the process noise and the measurement error and simulate the
full data set.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the data in long format, one row per individual and time point, with the
columns `id`, `time`, `y1`, `y2`, and `y3`.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`simStateSpace`](https://github.com/jeksterslab/simStateSpace)`)`\
`sim`` ``<-`` `[`SimSSMFixed`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md)`(`\
`  n ``=`` ``n``,`\
`  time ``=`` ``time``,`\
`  mu0 ``=`` ``mu0``,`\
`  sigma0_l ``=`` ``sigma0_l``,`\
`  alpha ``=`` ``alpha``,`\
`  beta ``=`` ``beta``,`\
`  psi_l ``=`` ``psi_l``,`\
`  nu ``=`` ``nu``,`\
`  lambda ``=`` ``lambda``,`\
`  theta_l ``=`` ``theta_l``,`\
`  type ``=`` ``0`\
`)`\
`data`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data``)`\
`#>   id time          y1         y2         y3`\
`#> 1  1    0 -0.66908807 0.16178434  0.2955693`\
`#> 2  1    1 -0.23021811 0.22873073 -0.2540210`\
`#> 3  1    2  0.93689389 0.09496693 -0.9267506`\
`#> 4  1    3  0.04445794 0.67521289 -0.1792168`\
`#> 5  1    4  0.15413707 0.82591676  0.8976536`\
`#> 6  1    5 -0.09943698 0.67154173  0.2853507`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``data``)`\
`#>        id              time              y1                  y2           `\
`#>  Min.   :  1.00   Min.   :  0.00   Min.   :-2.840534   Min.   :-2.671337  `\
`#>  1st Qu.: 25.75   1st Qu.: 49.75   1st Qu.:-0.419270   1st Qu.:-0.500990  `\
`#>  Median : 50.50   Median : 99.50   Median :-0.001263   Median : 0.007650  `\
`#>  Mean   : 50.50   Mean   : 99.50   Mean   :-0.002412   Mean   : 0.001159  `\
`#>  3rd Qu.: 75.25   3rd Qu.:149.25   3rd Qu.: 0.414442   3rd Qu.: 0.501791  `\
`#>  Max.   :100.00   Max.   :199.00   Max.   : 2.602905   Max.   : 3.043414  `\
`#>        y3           `\
`#>  Min.   :-2.675588  `\
`#>  1st Qu.:-0.454266  `\
`#>  Median :-0.005830  `\
`#>  Mean   :-0.001039  `\
`#>  3rd Qu.: 0.455474  `\
`#>  Max.   : 2.553926`\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-ssm-error-ssm-1.png)![](fig-vignettes-ssm-error-ssm-2.png)![](fig-vignettes-ssm-error-ssm-3.png)

## Model Fitting

The same simulated data are fitted with two programs that estimate state
space models by maximum likelihood with the Kalman filter: `dynr` and
`OpenMx`. The section Population Values and Estimates compares their
estimates with the population values.

### dynr

We fit the model with `dynr`, which estimates state space models by
maximum likelihood using the Kalman filter. The Kalman filter estimates
the latent states occasion by occasion, weighing the model’s prediction
against the new measurement.

A `dynr` model is built from pieces that match the equations above: the
data, the initial condition, the measurement model, the dynamic model,
and the noise. Each parameter gets a name (for example, `beta_2_1` for
\\(\beta\_{21}\\)), and the same name used in two places means one
shared parameter. The starting values are set to the population values
for convenience; with real data, they would be rough guesses.

#### Prepare Data

Tell `dynr` which columns hold the individual, the time, and the
observed variables.

\
`dynr_data`` ``<-`` ``dynr``::`[`dynr.data`](https://rdrr.io/pkg/dynr/man/dynr.data.html)`(`\
`  data ``=`` ``data``,`\
`  id ``=`` ``"id"``,`\
`  time ``=`` ``"time"``,`\
`  observed ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``, ``"y3"``)`\
`)`

#### Prepare Initial Condition

The mean \\(\boldsymbol{\mu}\_{0}\\) and covariance
\\(\boldsymbol{\Sigma}\_{0}\\) of the initial state are estimated. The
covariance matrix is symmetric, so `sigma0_2_1` appears both above and
below the diagonal.

\
`dynr_initial`` ``<-`` ``dynr``::`[`prep.initial`](https://rdrr.io/pkg/dynr/man/prep.initial.html)`(`\
`  values.inistate ``=`` ``mu0``,`\
`  params.inistate ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_1_1"``, ``"mu0_2_1"``, ``"mu0_3_1"``)``,`\
`  values.inicov ``=`` ``sigma0``,`\
`  params.inicov ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"sigma0_1_1"``, ``"sigma0_2_1"``, ``"sigma0_3_1"``,`\
`      ``"sigma0_2_1"``, ``"sigma0_2_2"``, ``"sigma0_3_2"``,`\
`      ``"sigma0_3_1"``, ``"sigma0_3_2"``, ``"sigma0_3_3"`\
`    ``)``,`\
`    nrow ``=`` ``3`\
`  ``)`\
`)`

#### Prepare Measurement Model

The loadings \\(\boldsymbol{\Lambda}\\) are an identity matrix with no
free parameters (`"fixed"`).

\
`dynr_measurement`` ``<-`` ``dynr``::`[`prep.measurement`](https://rdrr.io/pkg/dynr/man/prep.measurement.html)`(`\
`  values.load ``=`` `[`diag`](https://rdrr.io/r/base/diag.html)`(``3``)``,`\
`  params.load ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``3``, ncol ``=`` ``3``)``,`\
`  state.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta_1"``, ``"eta_2"``, ``"eta_3"``)``,`\
`  obs.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``, ``"y3"``)`\
`)`

#### Prepare Dynamic Process

Each formula is one row of the dynamic equation: the intercept
\\(\alpha\_{j}\\) plus the effects \\(\beta\_{jk}\\) of the three latent
variables at the previous occasion. `isContinuousTime = FALSE` makes
this a discrete-time model.

\
`dynr_dynamics`` ``<-`` ``dynr``::`[`prep.formulaDynamics`](https://rdrr.io/pkg/dynr/man/prep.formulaDynamics.html)`(`\
`  formula ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    ``eta_1`` ``~`` ``alpha_1_1`` ``*`` ``1`` ``+`` ``beta_1_1`` ``*`` ``eta_1`` ``+`` ``beta_1_2`` ``*`` ``eta_2`` ``+`` ``beta_1_3`` ``*`` ``eta_3``,`\
`    ``eta_2`` ``~`` ``alpha_2_1`` ``*`` ``1`` ``+`` ``beta_2_1`` ``*`` ``eta_1`` ``+`` ``beta_2_2`` ``*`` ``eta_2`` ``+`` ``beta_2_3`` ``*`` ``eta_3``,`\
`    ``eta_3`` ``~`` ``alpha_3_1`` ``*`` ``1`` ``+`` ``beta_3_1`` ``*`` ``eta_1`` ``+`` ``beta_3_2`` ``*`` ``eta_2`` ``+`` ``beta_3_3`` ``*`` ``eta_3`\
`  ``)``,`\
`  startval ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    alpha_1_1 ``=`` ``alpha``[``1``]``, alpha_2_1 ``=`` ``alpha``[``2``]``, alpha_3_1 ``=`` ``alpha``[``3``]``,`\
`    beta_1_1 ``=`` ``beta``[``1``, ``1``]``, beta_1_2 ``=`` ``beta``[``1``, ``2``]``, beta_1_3 ``=`` ``beta``[``1``, ``3``]``,`\
`    beta_2_1 ``=`` ``beta``[``2``, ``1``]``, beta_2_2 ``=`` ``beta``[``2``, ``2``]``, beta_2_3 ``=`` ``beta``[``2``, ``3``]``,`\
`    beta_3_1 ``=`` ``beta``[``3``, ``1``]``, beta_3_2 ``=`` ``beta``[``3``, ``2``]``, beta_3_3 ``=`` ``beta``[``3``, ``3``]`\
`  ``)``,`\
`  isContinuousTime ``=`` ``FALSE`\
`)`

#### Prepare Process Noise

The process noise covariance \\(\boldsymbol{\Psi}\\) is estimated in
full. The measurement error covariance \\(\boldsymbol{\Theta}\\) is
diagonal: the three error variances are estimated, and the error
covariances are fixed at zero.

\
`dynr_noise`` ``<-`` ``dynr``::`[`prep.noise`](https://rdrr.io/pkg/dynr/man/prep.noise.html)`(`\
`  values.latent ``=`` ``psi``,`\
`  params.latent ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"psi_1_1"``, ``"psi_2_1"``, ``"psi_3_1"``,`\
`      ``"psi_2_1"``, ``"psi_2_2"``, ``"psi_3_2"``,`\
`      ``"psi_3_1"``, ``"psi_3_2"``, ``"psi_3_3"`\
`    ``)``,`\
`    nrow ``=`` ``3`\
`  ``)``,`\
`  values.observed ``=`` ``theta``,`\
`  params.observed ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"theta_1_1"``, ``"fixed"``, ``"fixed"``,`\
`      ``"fixed"``, ``"theta_2_2"``, ``"fixed"``,`\
`      ``"fixed"``, ``"fixed"``, ``"theta_3_3"`\
`    ``)``,`\
`    nrow ``=`` ``3`\
`  ``)`\
`)`

#### Prepare the Model

The pieces are put together into one model. `dynr` writes the model as C
code to `outfile` and compiles it.

\
`model`` ``<-`` ``dynr``::`[`dynr.model`](https://rdrr.io/pkg/dynr/man/dynr.model.html)`(`\
`  data ``=`` ``dynr_data``,`\
`  initial ``=`` ``dynr_initial``,`\
`  measurement ``=`` ``dynr_measurement``,`\
`  dynamics ``=`` ``dynr_dynamics``,`\
`  noise ``=`` ``dynr_noise``,`\
`  outfile ``=`` ``"ssm.c"`\
`)`

The plot below shows the model equations with the parameter names.

![](fig-vignettes-ssm-model-ssm-1.png)

#### Fit the Model

\
`results`` ``<-`` ``dynr``::`[`dynr.cook`](https://rdrr.io/pkg/dynr/man/dynr.cook.html)`(`\
`  ``model``,`\
`  debug_flag ``=`` ``TRUE``,`\
`  verbose ``=`` ``FALSE`\
`)`\
`#> [1] "Get ready!!!!"`\
`#> using C compiler: ‘gcc (Ubuntu 15.2.0-16ubuntu1) 15.2.0’`\
`#> Optimization function called.`\
`#> Starting Hessian calculation ...`\
`#> Finished Hessian calculation.`\
`#> Original exit flag:  3 `\
`#> Modified exit flag:  3 `\
`#> Optimization terminated successfully: ftol_rel or ftol_abs was reached. `\
`#> Original fitted parameters:  -0.0006566288 0.00164596 -0.001510021 0.7052377 `\
`#> -0.006132876 0.0003878288 0.5112233 0.5887748 0.01241123 -0.08998742 0.3907879 `\
`#> 0.4952418 -2.32728 0.0126181 -0.03370301 -2.24728 -0.001445368 -2.352586 `\
`#> -1.611374 -1.660788 -1.586927 -0.0008762037 0.02978035 0.07119585 -1.830427 `\
`#> 0.4567187 0.06185989 -0.7307456 0.254486 -1.438148 `\
`#> `\
`#> Transformed fitted parameters:  -0.0006566288 0.00164596 -0.001510021 0.7052377 `\
`#> -0.006132876 0.0003878288 0.5112233 0.5887748 0.01241123 -0.08998742 0.3907879 `\
`#> 0.4952418 0.0975608 0.001231032 -0.003288092 0.1057019 -0.0001942451 0.0952339 `\
`#> 0.1996131 0.1899893 0.2045533 -0.0008762037 0.02978035 0.07119585 0.1603451 `\
`#> 0.07323263 0.009918933 0.5149965 0.1270778 0.2691672 `\
`#> `\
`#> Doing end processing`\
`#> Successful trial`\
`#> Total Time: 4.779007 `\
`#> Backend Time: 4.771481`

#### Summary

The summary lists each estimate with its standard error, \\(t\\) value,
confidence interval, and \\(p\\) value.

\
[`summary`](https://rdrr.io/r/base/summary.html)`(``results``)`\
`#> Coefficients:`\
`#>              Estimate Std. Error t value   ci.lower   ci.upper`\
`#> alpha_1_1  -0.0006566  0.0024115  -0.272 -0.0053832  0.0040699`\
`#> alpha_2_1   0.0016460  0.0030958   0.532 -0.0044218  0.0077137`\
`#> alpha_3_1  -0.0015100  0.0029970  -0.504 -0.0073839  0.0043639`\
`#> beta_1_1    0.7052377  0.0218228  32.317  0.6624658  0.7480097`\
`#> beta_1_2   -0.0061329  0.0134549  -0.456 -0.0325039  0.0202382`\
`#> beta_1_3    0.0003878  0.0102911   0.038 -0.0197823  0.0205579`\
`#> beta_2_1    0.5112233  0.0234282  21.821  0.4653048  0.5571417`\
`#> beta_2_2    0.5887748  0.0180029  32.704  0.5534897  0.6240599`\
`#> beta_2_3    0.0124112  0.0130672   0.950 -0.0132000  0.0380224`\
`#> beta_3_1   -0.0899874  0.0170005  -5.293 -0.1233079 -0.0566670`\
`#> beta_3_2    0.3907879  0.0168012  23.260  0.3578582  0.4237176`\
`#> beta_3_3    0.4952418  0.0164473  30.111  0.4630056  0.5274780`\
`#> psi_1_1     0.0975608  0.0068530  14.236  0.0841292  0.1109924`\
`#> psi_2_1     0.0012310  0.0024808   0.496 -0.0036313  0.0060933`\
`#> psi_3_1    -0.0032881  0.0023068  -1.425 -0.0078094  0.0012332`\
`#> psi_2_2     0.1057019  0.0060019  17.612  0.0939384  0.1174653`\
`#> psi_3_2    -0.0001942  0.0023068  -0.084 -0.0047156  0.0043271`\
`#> psi_3_3     0.0952339  0.0067865  14.033  0.0819326  0.1085352`\
`#> theta_1_1   0.1996131  0.0057442  34.751  0.1883547  0.2108714`\
`#> theta_2_2   0.1899893  0.0058134  32.681  0.1785952  0.2013833`\
`#> theta_3_3   0.2045533  0.0069223  29.550  0.1909859  0.2181207`\
`#> mu0_1_1    -0.0008762  0.0527703  -0.017 -0.1043041  0.1025517`\
`#> mu0_2_1     0.0297803  0.0812463   0.367 -0.1294595  0.1890202`\
`#> mu0_3_1     0.0711959  0.0666440   1.068 -0.0594241  0.2018158`\
`#> sigma0_1_1  0.1603451  0.0401133   3.997  0.0817244  0.2389659`\
`#> sigma0_2_1  0.0732326  0.0444635   1.647 -0.0139141  0.1603794`\
`#> sigma0_3_1  0.0099189  0.0358029   0.277 -0.0602535  0.0800914`\
`#> sigma0_2_2  0.5149965  0.0931946   5.526  0.3323384  0.6976546`\
`#> sigma0_3_2  0.1270778  0.0552607   2.300  0.0187688  0.2353869`\
`#> sigma0_3_3  0.2691672  0.0632758   4.254  0.1451489  0.3931856`\
`#>                       Pr(>|t|)    `\
`#> alpha_1_1               0.3927    `\
`#> alpha_2_1               0.2975    `\
`#> alpha_3_1               0.3072    `\
`#> beta_1_1   <0.0000000000000002 ***`\
`#> beta_1_2                0.3243    `\
`#> beta_1_3                0.4850    `\
`#> beta_2_1   <0.0000000000000002 ***`\
`#> beta_2_2   <0.0000000000000002 ***`\
`#> beta_2_3                0.1711    `\
`#> beta_3_1   <0.0000000000000002 ***`\
`#> beta_3_2   <0.0000000000000002 ***`\
`#> beta_3_3   <0.0000000000000002 ***`\
`#> psi_1_1    <0.0000000000000002 ***`\
`#> psi_2_1                 0.3099    `\
`#> psi_3_1                 0.0770 .  `\
`#> psi_2_2    <0.0000000000000002 ***`\
`#> psi_3_2                 0.4664    `\
`#> psi_3_3    <0.0000000000000002 ***`\
`#> theta_1_1  <0.0000000000000002 ***`\
`#> theta_2_2  <0.0000000000000002 ***`\
`#> theta_3_3  <0.0000000000000002 ***`\
`#> mu0_1_1                 0.4934    `\
`#> mu0_2_1                 0.3570    `\
`#> mu0_3_1                 0.1427    `\
`#> sigma0_1_1 <0.0000000000000002 ***`\
`#> sigma0_2_1              0.0498 *  `\
`#> sigma0_3_1              0.3909    `\
`#> sigma0_2_2 <0.0000000000000002 ***`\
`#> sigma0_3_2              0.0107 *  `\
`#> sigma0_3_3 <0.0000000000000002 ***`\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> -2 log-likelihood value at convergence = 105691.54`\
`#> AIC = 105751.54`\
`#> BIC = 105988.64`

### OpenMx: State Space Model

`OpenMx` (Neale et al., 2016) fits the same state space model with its
own Kalman filter
([`OpenMx::mxExpectationStateSpace()`](https://rdrr.io/pkg/OpenMx/man/mxExpectationStateSpace.html);
Hunter (2018)). Its matrices have the names used in engineering: `A` is
the transition matrix \\(\boldsymbol{\beta}\\), `B` holds the intercepts
of the dynamic model, `C` is the loading matrix
\\(\boldsymbol{\Lambda}\\), `D` holds the intercepts of the measurement
model, `Q` is the process noise covariance matrix, and `R` is the
measurement error covariance matrix. The parameter labels are the dynr
names; as in dynr, a label used in two places is one shared parameter,
and the starting values are the population values.

Two details differ from dynr. First, OpenMx reads the data as a single
time series. To start the filter again at the first occasion of every
individual, each individual gets a group, and all groups use the
matrices of one parent model, so every individual has the same
parameters. Second, OpenMx starts the filter one occasion before the
first row, so the first row has the state mean \\(\mathbf{A}
\mathbf{x}\_{0} + \mathbf{B} \mathbf{u}\\) and the state covariance
\\(\mathbf{A} \mathbf{P}\_{0} \mathbf{A}^{\prime} + \mathbf{Q}\\), with
the constant input \\(\mathbf{u} = 1\\). Setting \\(\mathbf{x}\_{0} =
\mathbf{A}^{-1} \left( \boldsymbol{\mu}\_{0} - \mathbf{B} \mathbf{u}
\right)\\) and \\(\mathbf{P}\_{0} = \mathbf{A}^{-1} \left(
\boldsymbol{\Sigma}\_{0} - \mathbf{Q} \right) \left( \mathbf{A}^{-1}
\right)^{\prime}\\) gives the first occasion the initial condition of
the model, so the likelihood is the same as in dynr. This needs an
invertible \\(\mathbf{A}\\), which holds here.

The two helper functions below are the same in every vignette of the
package. `MxMat()` makes an OpenMx matrix from values and dynr-style
labels (`"fixed"` marks a fixed element), and `FitMxSSM()` builds the
groups and fits the model.

\
`MxMat`` ``<-`` ``function``(``values``, ``labels``, ``name``)`` ``{`\
`  ``values`` ``<-`` `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``values``)`\
`  ``labels`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` ``labels``,`\
`    nrow ``=`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``values``)``,`\
`    ncol ``=`` `[`ncol`](https://rdrr.io/r/base/nrow.html)`(``values``)`\
`  ``)`\
`  ``free`` ``<-`` ``labels`` ``!=`` ``"fixed"`\
`  ``labels``[``!``free``]`` ``<-`` ``NA`\
`  ``OpenMx``::`[`mxMatrix`](https://rdrr.io/pkg/OpenMx/man/mxMatrix.html)`(`\
`    type ``=`` ``"Full"``,`\
`    values ``=`` ``values``,`\
`    free ``=`` ``free``,`\
`    labels ``=`` ``labels``,`\
`    name ``=`` ``name`\
`  ``)`\
`}`\
`FitMxSSM`` ``<-`` ``function``(``data``, ``observed``, ``matrices``, ``name``)`` ``{`\
`  ``# the groups refer to the matrices of the parent model`\
`  ``shared`` ``<-`` ``function``(``x``)`` ``{`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``name``, ``"."``, ``x``)`\
`  ``}`\
`  ``expectation`` ``<-`` ``OpenMx``::`[`mxExpectationStateSpace`](https://rdrr.io/pkg/OpenMx/man/mxExpectationStateSpace.html)`(`\
`    A ``=`` ``shared``(``"A"``)``,`\
`    B ``=`` ``shared``(``"B"``)``,`\
`    C ``=`` ``shared``(``"C"``)``,`\
`    D ``=`` ``shared``(``"D"``)``,`\
`    Q ``=`` ``shared``(``"Q"``)``,`\
`    R ``=`` ``shared``(``"R"``)``,`\
`    x0 ``=`` ``shared``(``"x0"``)``,`\
`    P0 ``=`` ``shared``(``"P0"``)``,`\
`    u ``=`` ``shared``(``"u"``)`\
`  ``)`\
`  ``# OpenMx needs the names of the observed variables and of the states on C`\
`  ``for`` ``(``i`` ``in`` `[`seq_along`](https://rdrr.io/r/base/seq.html)`(``matrices``)``)`` ``{`\
`    ``if`` ``(``matrices``[[``i``]``]``$``name`` ``==`` ``"C"``)`` ``{`\
`      `[`dimnames`](https://rdrr.io/r/base/dimnames.html)`(``matrices``[[``i``]``]``)`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`        ``observed``,`\
`        `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"eta_"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``matrices``[[``i``]``]``$``values``)``)``)`\
`      ``)`\
`    ``}`\
`  ``}`\
`  ``persons`` ``<-`` `[`split`](https://rdrr.io/r/base/split.html)`(`\
`    x ``=`` ``data``[``, ``observed``, drop ``=`` ``FALSE``]``,`\
`    f ``=`` ``data``[``, ``"id"``]`\
`  ``)`\
`  ``groups`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"id_"``, `[`names`](https://rdrr.io/r/base/names.html)`(``persons``)``)`\
`  ``submodels`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(`\
`    X ``=`` `[`seq_along`](https://rdrr.io/r/base/seq.html)`(``persons``)``,`\
`    FUN ``=`` ``function``(``i``)`` ``{`\
`      ``OpenMx``::`[`mxModel`](https://rdrr.io/pkg/OpenMx/man/mxModel.html)`(`\
`        name ``=`` ``groups``[``i``]``,`\
`        ``expectation``,`\
`        ``OpenMx``::`[`mxFitFunctionML`](https://rdrr.io/pkg/OpenMx/man/mxFitFunctionML.html)`(``)``,`\
`        ``OpenMx``::`[`mxData`](https://rdrr.io/pkg/OpenMx/man/mxData.html)`(`\
`          observed ``=`` ``persons``[[``i``]``]``,`\
`          type ``=`` ``"raw"`\
`        ``)`\
`      ``)`\
`    ``}`\
`  ``)`\
`  ``model`` ``<-`` ``OpenMx``::`[`mxModel`](https://rdrr.io/pkg/OpenMx/man/mxModel.html)`(`\
`    name ``=`` ``name``,`\
`    ``matrices``,`\
`    ``OpenMx``::`[`mxMatrix`](https://rdrr.io/pkg/OpenMx/man/mxMatrix.html)`(`\
`      type ``=`` ``"Unit"``,`\
`      nrow ``=`` ``1``,`\
`      ncol ``=`` ``1``,`\
`      name ``=`` ``"u"`\
`    ``)``,`\
`    ``# the first occasion gets the initial condition of the model`\
`    ``OpenMx``::`[`mxAlgebra`](https://rdrr.io/pkg/OpenMx/man/mxAlgebra.html)`(`\
`      expression ``=`` `[`solve`](https://rdrr.io/r/base/solve.html)`(``A``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``(``mu0`` ``-`` ``B`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``u``)``,`\
`      name ``=`` ``"x0"`\
`    ``)``,`\
`    ``OpenMx``::`[`mxAlgebra`](https://rdrr.io/pkg/OpenMx/man/mxAlgebra.html)`(`\
`      expression ``=`` `[`solve`](https://rdrr.io/r/base/solve.html)`(``A``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``(``sigma0`` ``-`` ``Q``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` `[`t`](https://rdrr.io/r/base/t.html)`(`[`solve`](https://rdrr.io/r/base/solve.html)`(``A``)``)``,`\
`      name ``=`` ``"P0"`\
`    ``)``,`\
`    ``submodels``,`\
`    ``OpenMx``::`[`mxFitFunctionMultigroup`](https://rdrr.io/pkg/OpenMx/man/mxFitFunctionMultigroup.html)`(``groups ``=`` ``groups``)`\
`  ``)`\
`  ``# retries from perturbed starting values only if the first run fails`\
`  ``OpenMx``::`[`mxTryHard`](https://rdrr.io/pkg/OpenMx/man/mxTryHard.html)`(`\
`    model ``=`` ``model``,`\
`    extraTries ``=`` ``5``,`\
`    silent ``=`` ``TRUE`\
`  ``)`\
`}`

The matrices of the model, with the dynr labels.

\
`labels_beta`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"beta_1_1"``, ``"beta_2_1"``, ``"beta_3_1"``,`\
`    ``"beta_1_2"``, ``"beta_2_2"``, ``"beta_3_2"``,`\
`    ``"beta_1_3"``, ``"beta_2_3"``, ``"beta_3_3"`\
`  ``)``,`\
`  nrow ``=`` ``3`\
`)`\
`labels_psi`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"psi_1_1"``, ``"psi_2_1"``, ``"psi_3_1"``,`\
`    ``"psi_2_1"``, ``"psi_2_2"``, ``"psi_3_2"``,`\
`    ``"psi_3_1"``, ``"psi_3_2"``, ``"psi_3_3"`\
`  ``)``,`\
`  nrow ``=`` ``3`\
`)`\
`labels_theta`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"theta_1_1"``, ``"fixed"``, ``"fixed"``,`\
`    ``"fixed"``, ``"theta_2_2"``, ``"fixed"``,`\
`    ``"fixed"``, ``"fixed"``, ``"theta_3_3"`\
`  ``)``,`\
`  nrow ``=`` ``3`\
`)`\
`labels_sigma0`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"sigma0_1_1"``, ``"sigma0_2_1"``, ``"sigma0_3_1"``,`\
`    ``"sigma0_2_1"``, ``"sigma0_2_2"``, ``"sigma0_3_2"``,`\
`    ``"sigma0_3_1"``, ``"sigma0_3_2"``, ``"sigma0_3_3"`\
`  ``)``,`\
`  nrow ``=`` ``3`\
`)`\
`mx_matrices`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`  ``MxMat``(``values ``=`` ``beta``, labels ``=`` ``labels_beta``, name ``=`` ``"A"``)``,`\
`  ``MxMat``(`\
`    values ``=`` ``alpha``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"alpha_1_1"``, ``"alpha_2_1"``, ``"alpha_3_1"``)``,`\
`    name ``=`` ``"B"`\
`  ``)``,`\
`  ``MxMat``(``values ``=`` ``lambda``, labels ``=`` ``"fixed"``, name ``=`` ``"C"``)``,`\
`  ``MxMat``(``values ``=`` ``nu``, labels ``=`` ``"fixed"``, name ``=`` ``"D"``)``,`\
`  ``MxMat``(``values ``=`` ``psi``, labels ``=`` ``labels_psi``, name ``=`` ``"Q"``)``,`\
`  ``MxMat``(``values ``=`` ``theta``, labels ``=`` ``labels_theta``, name ``=`` ``"R"``)``,`\
`  ``MxMat``(`\
`    values ``=`` ``mu0``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_1_1"``, ``"mu0_2_1"``, ``"mu0_3_1"``)``,`\
`    name ``=`` ``"mu0"`\
`  ``)``,`\
`  ``MxMat``(``values ``=`` ``sigma0``, labels ``=`` ``labels_sigma0``, name ``=`` ``"sigma0"``)`\
`)`

``` r
fit_mx_ssm <- FitMxSSM(
  data = data,
  observed = c("y1", "y2", "y3"),
  matrices = mx_matrices,
  name = "ssm"
)
#> 
Beginning initial fit attempt
Fit attempt 0, fit=105691.534423263, new current best! (was 105714.907085852)
                                                                             
#> 
#>  Solution found!  Final fit=105691.53 (started at 105714.91)  (1 attempt(s): 1 valid, 0 errors)
summary(fit_mx_ssm)
#> Summary of ssm 
#>  
#> free parameters:
#>          name matrix row col      Estimate   Std.Error A
#> 1    beta_1_1      A   1   1  0.7053981077 0.021838141  
#> 2    beta_2_1      A   2   1  0.5112620939 0.023493351  
#> 3    beta_3_1      A   3   1 -0.0900410232 0.017085396  
#> 4    beta_1_2      A   1   2 -0.0061800615 0.013468019  
#> 5    beta_2_2      A   2   2  0.5888519508 0.018093205  
#> 6    beta_3_2      A   3   2  0.3908662892 0.016906723  
#> 7    beta_1_3      A   1   3  0.0004417040 0.010298111  
#> 8    beta_2_3      A   2   3  0.0123737293 0.013111506  
#> 9    beta_3_3      A   3   3  0.4952567006 0.016525427  
#> 10  alpha_1_1      B   1   1 -0.0006118823 0.002410796  
#> 11  alpha_2_1      B   2   1  0.0016758350 0.003095680  
#> 12  alpha_3_1      B   3   1 -0.0014686848 0.002996973  
#> 13    psi_1_1      Q   1   1  0.0975097042 0.006853796  
#> 14    psi_2_1      Q   2   1  0.0012306351 0.002479268  
#> 15    psi_3_1      Q   3   1 -0.0032844672 0.002314293  
#> 16    psi_2_2      Q   2   2  0.1056516849 0.006032225  
#> 17    psi_3_2      Q   3   2 -0.0001907060 0.002309151  
#> 18    psi_3_3      Q   3   3  0.0952071001 0.006810719  
#> 19  theta_1_1      R   1   1  0.1996612442 0.005745493  
#> 20  theta_2_2      R   2   2  0.1900372708 0.005845908  
#> 21  theta_3_3      R   3   3  0.2045821252 0.006947999  
#> 22    mu0_1_1    mu0   1   1 -0.0008607865 0.052817445  
#> 23    mu0_2_1    mu0   2   1  0.0298171307 0.080975260  
#> 24    mu0_3_1    mu0   3   1  0.0712427739 0.066588868  
#> 25 sigma0_1_1 sigma0   1   1  0.1603073688 0.039618205  
#> 26 sigma0_2_1 sigma0   2   1  0.0732247803 0.043153542  
#> 27 sigma0_3_1 sigma0   3   1  0.0099103002 0.035256844  
#> 28 sigma0_2_2 sigma0   2   2  0.5149410876 0.093136953  
#> 29 sigma0_3_2 sigma0   3   2  0.1270997132 0.055085499  
#> 30 sigma0_3_3 sigma0   3   3  0.2691550668 0.063077429  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:             30                  59970              105691.5
#>    Saturated:             NA                     NA                    NA
#> Independence:             NA                     NA                    NA
#> Number of observations/statistics: 20000/60000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:      -14248.47               105751.5                 105751.6
#> BIC:     -488220.61               105988.6                 105893.3
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-11-02 06:48:52 
#> Wall clock time: 122.085 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

## Population Values and Estimates

Each table below puts the population values used to simulate the data
next to the estimates of one method. Bias is the estimate minus the
population value. Because each model is fitted to a single simulated
data set, this is the error of one estimate, not a Monte Carlo bias,
which would average over many replications. The 95% confidence interval
is the Wald interval, \\(\hat{\theta} \pm 1.96 \times \mathrm{SE}\\),
unless the caption says otherwise, and `Covers` shows whether it
contains the population value. The parameter names are the dynr names in
every table.

\
`EstimatesTable`` ``<-`` ``function``(``population``, ``estimate``, ``se``, ``ll`` ``=`` ``NULL``, ``ul`` ``=`` ``NULL``)`` ``{`\
`  ``z`` ``<-`` ``stats``::`[`qnorm`](https://rdrr.io/r/stats/Normal.html)`(``p ``=`` ``0.975``)`\
`  ``if`` ``(`[`is.null`](https://rdrr.io/r/base/NULL.html)`(``ll``)``)`` ``{`\
`    ``ll`` ``<-`` ``estimate`` ``-`` ``z`` ``*`` ``se`\
`  ``}`\
`  ``if`` ``(`[`is.null`](https://rdrr.io/r/base/NULL.html)`(``ul``)``)`` ``{`\
`    ``ul`` ``<-`` ``estimate`` ``+`` ``z`` ``*`` ``se`\
`  ``}`\
`  `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`    Parameter ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(``population``)``,`\
`    Population ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``population``)``,`\
`    Estimate ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``estimate``)``,`\
`    SE ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``se``)``,`\
`    Bias ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``estimate`` ``-`` ``population``)``,`\
``     `2.5%`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``ll``)``,`\
``     `97.5%`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``ul``)``,`\
`    Covers ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``population`` ``>=`` ``ll`` ``&`` ``population`` ``<=`` ``ul``)``,`\
`    check.names ``=`` ``FALSE`\
`  ``)`\
`}`

\
`# the population value for each dynr parameter name,`\
`# for example, "beta_1_2" is beta[1, 2] and "mu0_2_1" is mu0[2]`\
`population`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`\
`  X ``=`` `[`strsplit`](https://rdrr.io/r/base/strsplit.html)`(``x ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(`[`coef`](https://rdrr.io/r/stats/coef.html)`(``results``)``)``, split ``=`` ``"_"``)``,`\
`  FUN ``=`` ``function``(``x``)`` ``{`\
`    `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(`[`get`](https://rdrr.io/r/base/get.html)`(``x``[``1``]``)``)``[`[`as.integer`](https://rdrr.io/r/base/integer.html)`(``x``[``2``]``)``, `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``x``[``3``]``)``]`\
`  ``}``,`\
`  FUN.VALUE ``=`` `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``1``)`\
`)`\
[`names`](https://rdrr.io/r/base/names.html)`(``population``)`` ``<-`` `[`names`](https://rdrr.io/r/base/names.html)`(`[`coef`](https://rdrr.io/r/stats/coef.html)`(``results``)``)`

### dynr

\
`dynr_est`` ``<-`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(``results``)`\
`dynr_se`` ``<-`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(`[`diag`](https://rdrr.io/r/base/diag.html)`(`[`vcov`](https://rdrr.io/r/stats/vcov.html)`(``results``)``)``)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` ``EstimatesTable``(`\
`    population ``=`` ``population``,`\
`    estimate ``=`` ``dynr_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    se ``=`` ``dynr_se``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and dynr estimates (one simulated data set)"`\
`)`

| Parameter  | Population | Estimate |    SE |   Bias |   2.5% |  97.5% | Covers |
|:-----------|-----------:|---------:|------:|-------:|-------:|-------:|:-------|
| alpha_1_1  |      0.000 |   -0.001 | 0.002 | -0.001 | -0.005 |  0.004 | TRUE   |
| alpha_2_1  |      0.000 |    0.002 | 0.003 |  0.002 | -0.004 |  0.008 | TRUE   |
| alpha_3_1  |      0.000 |   -0.002 | 0.003 | -0.002 | -0.007 |  0.004 | TRUE   |
| beta_1_1   |      0.700 |    0.705 | 0.022 |  0.005 |  0.662 |  0.748 | TRUE   |
| beta_1_2   |      0.000 |   -0.006 | 0.013 | -0.006 | -0.033 |  0.020 | TRUE   |
| beta_1_3   |      0.000 |    0.000 | 0.010 |  0.000 | -0.020 |  0.021 | TRUE   |
| beta_2_1   |      0.500 |    0.511 | 0.023 |  0.011 |  0.465 |  0.557 | TRUE   |
| beta_2_2   |      0.600 |    0.589 | 0.018 | -0.011 |  0.553 |  0.624 | TRUE   |
| beta_2_3   |      0.000 |    0.012 | 0.013 |  0.012 | -0.013 |  0.038 | TRUE   |
| beta_3_1   |     -0.100 |   -0.090 | 0.017 |  0.010 | -0.123 | -0.057 | TRUE   |
| beta_3_2   |      0.400 |    0.391 | 0.017 | -0.009 |  0.358 |  0.424 | TRUE   |
| beta_3_3   |      0.500 |    0.495 | 0.016 | -0.005 |  0.463 |  0.527 | TRUE   |
| psi_1_1    |      0.100 |    0.098 | 0.007 | -0.002 |  0.084 |  0.111 | TRUE   |
| psi_2_1    |      0.000 |    0.001 | 0.002 |  0.001 | -0.004 |  0.006 | TRUE   |
| psi_3_1    |      0.000 |   -0.003 | 0.002 | -0.003 | -0.008 |  0.001 | TRUE   |
| psi_2_2    |      0.100 |    0.106 | 0.006 |  0.006 |  0.094 |  0.117 | TRUE   |
| psi_3_2    |      0.000 |    0.000 | 0.002 |  0.000 | -0.005 |  0.004 | TRUE   |
| psi_3_3    |      0.100 |    0.095 | 0.007 | -0.005 |  0.082 |  0.109 | TRUE   |
| theta_1_1  |      0.200 |    0.200 | 0.006 |  0.000 |  0.188 |  0.211 | TRUE   |
| theta_2_2  |      0.200 |    0.190 | 0.006 | -0.010 |  0.179 |  0.201 | TRUE   |
| theta_3_3  |      0.200 |    0.205 | 0.007 |  0.005 |  0.191 |  0.218 | TRUE   |
| mu0_1_1    |      0.000 |   -0.001 | 0.053 | -0.001 | -0.104 |  0.103 | TRUE   |
| mu0_2_1    |      0.000 |    0.030 | 0.081 |  0.030 | -0.129 |  0.189 | TRUE   |
| mu0_3_1    |      0.000 |    0.071 | 0.067 |  0.071 | -0.059 |  0.202 | TRUE   |
| sigma0_1_1 |      0.196 |    0.160 | 0.040 | -0.036 |  0.082 |  0.239 | TRUE   |
| sigma0_2_1 |      0.118 |    0.073 | 0.044 | -0.045 | -0.014 |  0.160 | TRUE   |
| sigma0_3_1 |      0.030 |    0.010 | 0.036 | -0.020 | -0.060 |  0.080 | TRUE   |
| sigma0_2_2 |      0.344 |    0.515 | 0.093 |  0.171 |  0.332 |  0.698 | TRUE   |
| sigma0_3_2 |      0.138 |    0.127 | 0.055 | -0.011 |  0.019 |  0.235 | TRUE   |
| sigma0_3_3 |      0.266 |    0.269 | 0.063 |  0.003 |  0.145 |  0.393 | TRUE   |

Population values and dynr estimates (one simulated data set) {.table
style="width:99%;"}

### OpenMx: State Space Model

\
`mx_ssm_est`` ``<-`` ``OpenMx``::`[`omxGetParameters`](https://rdrr.io/pkg/OpenMx/man/omxGetParameters.html)`(``fit_mx_ssm``)`\
`mx_ssm_se`` ``<-`` ``fit_mx_ssm``$``output``$``standardErrors``[``, ``1``]`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` ``EstimatesTable``(`\
`    population ``=`` ``population``,`\
`    estimate ``=`` ``mx_ssm_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    se ``=`` ``mx_ssm_se``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` `[`paste`](https://rdrr.io/r/base/paste.html)`(`\
`    ``"Population values and OpenMx state space estimates"``,`\
`    ``"(one simulated data set)"`\
`  ``)`\
`)`

| Parameter  | Population | Estimate |    SE |   Bias |   2.5% |  97.5% | Covers |
|:-----------|-----------:|---------:|------:|-------:|-------:|-------:|:-------|
| alpha_1_1  |      0.000 |   -0.001 | 0.002 | -0.001 | -0.005 |  0.004 | TRUE   |
| alpha_2_1  |      0.000 |    0.002 | 0.003 |  0.002 | -0.004 |  0.008 | TRUE   |
| alpha_3_1  |      0.000 |   -0.001 | 0.003 | -0.001 | -0.007 |  0.004 | TRUE   |
| beta_1_1   |      0.700 |    0.705 | 0.022 |  0.005 |  0.663 |  0.748 | TRUE   |
| beta_1_2   |      0.000 |   -0.006 | 0.013 | -0.006 | -0.033 |  0.020 | TRUE   |
| beta_1_3   |      0.000 |    0.000 | 0.010 |  0.000 | -0.020 |  0.021 | TRUE   |
| beta_2_1   |      0.500 |    0.511 | 0.023 |  0.011 |  0.465 |  0.557 | TRUE   |
| beta_2_2   |      0.600 |    0.589 | 0.018 | -0.011 |  0.553 |  0.624 | TRUE   |
| beta_2_3   |      0.000 |    0.012 | 0.013 |  0.012 | -0.013 |  0.038 | TRUE   |
| beta_3_1   |     -0.100 |   -0.090 | 0.017 |  0.010 | -0.124 | -0.057 | TRUE   |
| beta_3_2   |      0.400 |    0.391 | 0.017 | -0.009 |  0.358 |  0.424 | TRUE   |
| beta_3_3   |      0.500 |    0.495 | 0.017 | -0.005 |  0.463 |  0.528 | TRUE   |
| psi_1_1    |      0.100 |    0.098 | 0.007 | -0.002 |  0.084 |  0.111 | TRUE   |
| psi_2_1    |      0.000 |    0.001 | 0.002 |  0.001 | -0.004 |  0.006 | TRUE   |
| psi_3_1    |      0.000 |   -0.003 | 0.002 | -0.003 | -0.008 |  0.001 | TRUE   |
| psi_2_2    |      0.100 |    0.106 | 0.006 |  0.006 |  0.094 |  0.117 | TRUE   |
| psi_3_2    |      0.000 |    0.000 | 0.002 |  0.000 | -0.005 |  0.004 | TRUE   |
| psi_3_3    |      0.100 |    0.095 | 0.007 | -0.005 |  0.082 |  0.109 | TRUE   |
| theta_1_1  |      0.200 |    0.200 | 0.006 |  0.000 |  0.188 |  0.211 | TRUE   |
| theta_2_2  |      0.200 |    0.190 | 0.006 | -0.010 |  0.179 |  0.201 | TRUE   |
| theta_3_3  |      0.200 |    0.205 | 0.007 |  0.005 |  0.191 |  0.218 | TRUE   |
| mu0_1_1    |      0.000 |   -0.001 | 0.053 | -0.001 | -0.104 |  0.103 | TRUE   |
| mu0_2_1    |      0.000 |    0.030 | 0.081 |  0.030 | -0.129 |  0.189 | TRUE   |
| mu0_3_1    |      0.000 |    0.071 | 0.067 |  0.071 | -0.059 |  0.202 | TRUE   |
| sigma0_1_1 |      0.196 |    0.160 | 0.040 | -0.036 |  0.083 |  0.238 | TRUE   |
| sigma0_2_1 |      0.118 |    0.073 | 0.043 | -0.045 | -0.011 |  0.158 | TRUE   |
| sigma0_3_1 |      0.030 |    0.010 | 0.035 | -0.020 | -0.059 |  0.079 | TRUE   |
| sigma0_2_2 |      0.344 |    0.515 | 0.093 |  0.171 |  0.332 |  0.697 | TRUE   |
| sigma0_3_2 |      0.138 |    0.127 | 0.055 | -0.011 |  0.019 |  0.235 | TRUE   |
| sigma0_3_3 |      0.266 |    0.269 | 0.063 |  0.003 |  0.146 |  0.393 | TRUE   |

Population values and OpenMx state space estimates (one simulated data
set) {.table style="width:99%;"}

### All Methods

The estimates of every method side by side. The methods maximize the
same likelihood, so their estimates should agree closely.

\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`    Parameter ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(``population``)``,`\
`    Population ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``population``)``,`\
``     `dynr`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``dynr_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``)``,`\
``     `OpenMx (SSM)`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``mx_ssm_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``)``,`\
`    check.names ``=`` ``FALSE`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and the estimates of every method"`\
`)`

| Parameter  | Population |   dynr | OpenMx (SSM) |
|:-----------|-----------:|-------:|-------------:|
| alpha_1_1  |      0.000 | -0.001 |       -0.001 |
| alpha_2_1  |      0.000 |  0.002 |        0.002 |
| alpha_3_1  |      0.000 | -0.002 |       -0.001 |
| beta_1_1   |      0.700 |  0.705 |        0.705 |
| beta_1_2   |      0.000 | -0.006 |       -0.006 |
| beta_1_3   |      0.000 |  0.000 |        0.000 |
| beta_2_1   |      0.500 |  0.511 |        0.511 |
| beta_2_2   |      0.600 |  0.589 |        0.589 |
| beta_2_3   |      0.000 |  0.012 |        0.012 |
| beta_3_1   |     -0.100 | -0.090 |       -0.090 |
| beta_3_2   |      0.400 |  0.391 |        0.391 |
| beta_3_3   |      0.500 |  0.495 |        0.495 |
| psi_1_1    |      0.100 |  0.098 |        0.098 |
| psi_2_1    |      0.000 |  0.001 |        0.001 |
| psi_3_1    |      0.000 | -0.003 |       -0.003 |
| psi_2_2    |      0.100 |  0.106 |        0.106 |
| psi_3_2    |      0.000 |  0.000 |        0.000 |
| psi_3_3    |      0.100 |  0.095 |        0.095 |
| theta_1_1  |      0.200 |  0.200 |        0.200 |
| theta_2_2  |      0.200 |  0.190 |        0.190 |
| theta_3_3  |      0.200 |  0.205 |        0.205 |
| mu0_1_1    |      0.000 | -0.001 |       -0.001 |
| mu0_2_1    |      0.000 |  0.030 |        0.030 |
| mu0_3_1    |      0.000 |  0.071 |        0.071 |
| sigma0_1_1 |      0.196 |  0.160 |        0.160 |
| sigma0_2_1 |      0.118 |  0.073 |        0.073 |
| sigma0_3_1 |      0.030 |  0.010 |        0.010 |
| sigma0_2_2 |      0.344 |  0.515 |        0.515 |
| sigma0_3_2 |      0.138 |  0.127 |        0.127 |
| sigma0_3_3 |      0.266 |  0.269 |        0.269 |

Population values and the estimates of every method {.table}

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, *17*(2), 303–332.
<https://doi.org/10.1080/10705511003661553>

Hunter, M. D. (2018). State space modeling in an open source, modular,
structural equation modeling environment. *Structural Equation Modeling:
A Multidisciplinary Journal*, *25*(2), 307–324.
<https://doi.org/10.1080/10705511.2017.1369354>

Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H., & Boker,
S. M. (2016). OpenMx 2.0: Extended structural equation and statistical
modeling. *Psychometrika*, *81*(2), 535–549.
<https://doi.org/10.1007/s11336-014-9435-8>

Ou, L., Hunter, M. D., & Chow, S.-M. (2019). What’s for dynr: A package
for linear and nonlinear dynamic modeling in R. *The R Journal*,
*11*(1), 91. <https://doi.org/10.32614/rj-2019-012>

R Core Team. (2026). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>
