# The Ornstein-Uhlenbeck Model

## Overview

The Ornstein-Uhlenbeck (OU) model is a continuous-time model for
variables that are pulled back toward an equilibrium level after they
are pushed away from it, for example, mood returning to a person’s usual
level after a stressful event. Continuous time means the model describes
change at every moment, not only at the occasions we measure, so the
same parameters apply whether measurements are one hour or one day
apart, and the spacing can differ between occasions.

The model is a state space model with two parts: a measurement model
that links the observed variables to the latent (unobserved) states, and
a dynamic model, here a stochastic differential equation, that describes
how the latent states change over time.

This tutorial does three things:

1.  It writes the multivariate OU model as a state space model.
2.  It simulates data from the model with
    [`SimSSMOUFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md),
    with every individual sharing the same parameters.
3.  It fits the model to the simulated data with the `dynr` package and
    compares the estimates with the values used to simulate the data.

No prior knowledge of state space models is needed. The vignettes on the
vector autoregressive model and the state space model cover the
discrete-time versions.

## Model

The measurement model is given by \\\[\begin{equation} \mathbf{y}\_{i,
t} = \boldsymbol{\nu} + \boldsymbol{\Lambda} \boldsymbol{\eta}\_{i, t} +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( \mathbf{0},
\boldsymbol{\Theta} \right) \end{equation}\\\] where \\(\mathbf{y}\_{i,
t}\\), \\(\boldsymbol{\eta}\_{i, t}\\), and
\\(\boldsymbol{\varepsilon}\_{i, t}\\) are random variables and
\\(\boldsymbol{\nu}\\), \\(\boldsymbol{\Lambda}\\), and
\\(\boldsymbol{\Theta}\\) are model parameters. \\(\mathbf{y}\_{i, t}\\)
represents a vector of observed random variables,
\\(\boldsymbol{\eta}\_{i, t}\\) a vector of latent random variables, and
\\(\boldsymbol{\varepsilon}\_{i, t}\\) a vector of random measurement
errors, at time \\(t\\) and individual \\(i\\). \\(\boldsymbol{\nu}\\)
denotes a vector of intercepts, \\(\boldsymbol{\Lambda}\\) a matrix of
factor loadings, and \\(\boldsymbol{\Theta}\\) the covariance matrix of
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
\\\[\begin{equation} \boldsymbol{\varepsilon}\_{i, t} =
\boldsymbol{\Theta}^{\frac{1}{2}} \mathbf{z}\_{i, t}, \quad
\mathrm{with} \quad \mathbf{z}\_{i, t} \sim \mathcal{N} \left(
\mathbf{0}, \mathbf{I} \right) \end{equation}\\\] where
\\(\mathbf{z}\_{i, t}\\) is a vector of independent standard normal
random variables and \\(\left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)
\left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)^{\prime} =
\boldsymbol{\Theta}\\) .

The dynamic structure is given by \\\[\begin{equation} \mathrm{d}
\boldsymbol{\eta}\_{i, t} = \boldsymbol{\Phi} \left(
\boldsymbol{\eta}\_{i, t} - \boldsymbol{\mu} \right) \mathrm{d}t +
\boldsymbol{\Sigma}^{\frac{1}{2}} \mathrm{d} \mathbf{W}\_{i, t}
\end{equation}\\\] where \\(\boldsymbol{\mu}\\) is the long-term mean or
equilibrium level, \\(\boldsymbol{\Phi}\\) is the rate of mean
reversion, determining how quickly the variable returns to its mean,
\\(\boldsymbol{\Sigma}\\) is the matrix of volatility or randomness in
the process, and \\(\mathrm{d}\boldsymbol{W}\\) is a Wiener process or
Brownian motion, which represents random fluctuations.

In words, \\(\mathrm{d} \boldsymbol{\eta}\_{i, t}\\) is the change in
the latent variables over a very short time \\(\mathrm{d}t\\). The
change has two parts: a pull toward the equilibrium
\\(\boldsymbol{\mu}\\), proportional to the distance
\\(\boldsymbol{\eta}\_{i, t} - \boldsymbol{\mu}\\), and random
fluctuations with covariance \\(\boldsymbol{\Sigma}\\) per unit time.

- With this sign convention, a diagonal element \\(\phi\_{jj}\\) is
  negative for a variable that returns to its equilibrium. The more
  negative it is, the faster the return.
- An off-diagonal element \\(\phi\_{jk}\\) is the effect of the distance
  of variable \\(k\\) from its equilibrium on the rate of change of
  variable \\(j\\). As in the discrete-time models, the row is the
  outcome and the column is the predictor.

For the process to settle around \\(\boldsymbol{\mu}\\), all eigenvalues
of \\(\boldsymbol{\Phi}\\) must have negative real parts.

Over a time interval \\(\Delta t\\), the OU model implies a
discrete-time vector autoregressive model with transition matrix
\\(\boldsymbol{\beta} = \exp \left( \boldsymbol{\Phi} \Delta t
\right)\\), the matrix exponential.
[`SimSSMOUFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md)
uses this exact discretization to simulate the data, so no approximation
error is introduced.

## Data Generation

### Notation

Let \\(t = 200\\) be the number of time points and \\(n = 10\\) be the
number of individuals.

Let the measurement model intercept vector \\(\boldsymbol{\nu}\\) be
given by

\\\[\begin{equation} \boldsymbol{\nu} = \left( \begin{array}{c} 0 \\\\ 0
\\\\ 0 \\\\ \end{array} \right) . \end{equation}\\\]

Let the factor loadings matrix \\(\boldsymbol{\Lambda}\\) be given by

\\\[\begin{equation} \boldsymbol{\Lambda} = \left( \begin{array}{ccc} 1
& 0 & 0 \\\\ 0 & 1 & 0 \\\\ 0 & 0 & 1 \\\\ \end{array} \right) .
\end{equation}\\\]

Let the measurement error covariance matrix \\(\boldsymbol{\Theta}\\) be
given by

\\\[\begin{equation} \boldsymbol{\Theta} = \left( \begin{array}{ccc} 0.2
& 0 & 0 \\\\ 0 & 0.2 & 0 \\\\ 0 & 0 & 0.2 \\\\ \end{array} \right) .
\end{equation}\\\]

With \\(\boldsymbol{\nu} = \mathbf{0}\\) and \\(\boldsymbol{\Lambda} =
\mathbf{I}\\), each observed variable measures one latent variable, and
\\(\boldsymbol{\Theta}\\) adds measurement error with variance 0.2 to
each.

Let the initial condition \\(\boldsymbol{\eta}\_{0}\\) be given by

\\\[\begin{equation} \boldsymbol{\eta}\_{0} \sim \mathcal{N} \left(
\boldsymbol{\mu}\_{\boldsymbol{\eta} \mid 0},
\boldsymbol{\Sigma}\_{\boldsymbol{\eta} \mid 0} \right)
\end{equation}\\\]

\\\[\begin{equation} \boldsymbol{\mu}\_{\boldsymbol{\eta} \mid 0} =
\left( \begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\ \end{array} \right)
\end{equation}\\\]

\\\[\begin{equation} \boldsymbol{\Sigma}\_{\boldsymbol{\eta} \mid 0} =
\left( \begin{array}{ccc} 0.3425148 & 0.3296023 & 0.0343817 \\\\
0.3296023 & 0.5664625 & 0.2545955 \\\\ 0.0343817 & 0.2545955 & 0.2999909
\\\\ \end{array} \right) . \end{equation}\\\]

Let the long-term mean vector \\(\boldsymbol{\mu}\\) be given by

\\\[\begin{equation} \boldsymbol{\mu} = \left( \begin{array}{c} 0 \\\\ 0
\\\\ 0 \\\\ \end{array} \right) . \end{equation}\\\]

Let the rate of mean reversion matrix \\(\boldsymbol{\Phi}\\) be given
by

\\\[\begin{equation} \boldsymbol{\Phi} = \left( \begin{array}{ccc}
-0.357 & 0 & 0 \\\\ 0.771 & -0.511 & 0 \\\\ -0.45 & 0.729 & -0.693 \\\\
\end{array} \right) . \end{equation}\\\]

Read the columns as predictors and the rows as outcomes. The diagonal
elements are all negative, so each variable returns to its equilibrium.
In addition, \\(\eta\_{1}\\) pushes up the change in \\(\eta\_{2}\\)
(0.771), \\(\eta\_{2}\\) pushes up the change in \\(\eta\_{3}\\)
(0.729), and \\(\eta\_{1}\\) pushes down the change in \\(\eta\_{3}\\)
(-0.45). Because \\(\boldsymbol{\Phi}\\) is lower triangular, its
eigenvalues are its diagonal elements, all negative, so the process is
stable.

Let the dynamic process noise covariance matrix
\\(\boldsymbol{\Sigma}\\) be given by

\\\[\begin{equation} \boldsymbol{\Sigma} = \left( \begin{array}{ccc}
0.2445556 & 0.0220159 & -0.0500476 \\\\ 0.0220159 & 0.070678 & 0.0153946
\\\\ -0.0500476 & 0.0153946 & 0.0755306 \\\\ \end{array} \right) .
\end{equation}\\\]

Let \\(\Delta t = 0.1\\).

The data are measured every 0.1 time units. For example, if the time
unit is one day, the 200 occasions span 19.9 days.

The initial condition describes where each individual starts. Here, the
mean vector and covariance matrix of \\(\boldsymbol{\eta}\_{0}\\) are
the model-implied mean and covariance of the stable process, computed
with
[`LinSDEMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanEta.md)
and
[`LinSDECovEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovEta.md),
so the data are stationary from the first time point.

### R Function Arguments

The simulation functions take covariance matrices as lower triangular
Cholesky factors (the arguments ending in `_l`).

\
`n`\
`#> [1] 10`\
`time`\
`#> [1] 200`\
`delta_t`\
`#> [1] 0.1`\
`mu0`\
`#> [1] 0 0 0`\
`sigma0`\
`#>           [,1]      [,2]      [,3]`\
`#> [1,] 0.3425148 0.3296023 0.0343817`\
`#> [2,] 0.3296023 0.5664625 0.2545955`\
`#> [3,] 0.0343817 0.2545955 0.2999909`\
`sigma0_l`` ``# sigma0_l <- t(chol(sigma0))`\
`#>            [,1]      [,2]      [,3]`\
`#> [1,] 0.58524763 0.0000000 0.0000000`\
`#> [2,] 0.56318429 0.4992855 0.0000000`\
`#> [3,] 0.05874726 0.4436540 0.3157701`\
`mu`\
`#> [1] 0 0 0`\
`phi`\
`#>        [,1]   [,2]   [,3]`\
`#> [1,] -0.357  0.000  0.000`\
`#> [2,]  0.771 -0.511  0.000`\
`#> [3,] -0.450  0.729 -0.693`\
`sigma`\
`#>             [,1]       [,2]        [,3]`\
`#> [1,]  0.24455556 0.02201587 -0.05004762`\
`#> [2,]  0.02201587 0.07067800  0.01539456`\
`#> [3,] -0.05004762 0.01539456  0.07553061`\
`sigma_l`` ``# sigma_l <- t(chol(sigma))`\
`#>             [,1]      [,2]     [,3]`\
`#> [1,]  0.49452559 0.0000000 0.000000`\
`#> [2,]  0.04451917 0.2620993 0.000000`\
`#> [3,] -0.10120330 0.0759256 0.243975`\
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

### Visualizing the Dynamics Without Measurement Error and Process Noise (n = 5 with Different Initial Condition)

To see what \\(\boldsymbol{\Phi}\\) alone does, we first simulate five
individuals with the process noise and the measurement error set to
zero. Each individual starts at a different point, and the trajectories
return to the equilibrium at a speed set by \\(\boldsymbol{\Phi}\\).

![](fig-vignettes-ou-no-error-ou-1.png)![](fig-vignettes-ou-no-error-ou-2.png)![](fig-vignettes-ou-no-error-ou-3.png)

### Using the `SimSSMOUFixed` Function from the `simStateSpace` Package to Simulate Data

Now we add the process noise and the measurement error and simulate the
full data set.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the data in long format, one row per individual and time point, with the
columns `id`, `time`, `y1`, `y2`, and `y3`. The `time` column holds the
actual measurement times (0, 0.1, 0.2, and so on), which the
continuous-time model uses.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`simStateSpace`](https://github.com/jeksterslab/simStateSpace)`)`\
`sim`` ``<-`` `[`SimSSMOUFixed`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md)`(`\
`  n ``=`` ``n``,`\
`  time ``=`` ``time``,`\
`  delta_t ``=`` ``delta_t``,`\
`  mu0 ``=`` ``mu0``,`\
`  sigma0_l ``=`` ``sigma0_l``,`\
`  mu ``=`` ``mu``,`\
`  phi ``=`` ``phi``,`\
`  sigma_l ``=`` ``sigma_l``,`\
`  nu ``=`` ``nu``,`\
`  lambda ``=`` ``lambda``,`\
`  theta_l ``=`` ``theta_l``,`\
`  type ``=`` ``0`\
`)`\
`data`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data``)`\
`#>   id time          y1         y2         y3`\
`#> 1  1  0.0  0.45047365 -1.2185550  0.6708638`\
`#> 2  1  0.1 -0.84190461  0.1242594  0.1813493`\
`#> 3  1  0.2  0.47289651 -0.2398710  0.5995107`\
`#> 4  1  0.3  0.04089817 -0.6547109  0.4674929`\
`#> 5  1  0.4 -1.37222946 -0.2010512 -0.3341885`\
`#> 6  1  0.5 -0.62410106  0.5262251 -0.1892906`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``data``)`\
`#>        id            time              y1                 y2         `\
`#>  Min.   : 1.0   Min.   : 0.000   Min.   :-2.25288   Min.   :-2.7497  `\
`#>  1st Qu.: 3.0   1st Qu.: 4.975   1st Qu.:-0.35833   1st Qu.:-0.3289  `\
`#>  Median : 5.5   Median : 9.950   Median : 0.09071   Median : 0.1856  `\
`#>  Mean   : 5.5   Mean   : 9.950   Mean   : 0.11195   Mean   : 0.1852  `\
`#>  3rd Qu.: 8.0   3rd Qu.:14.925   3rd Qu.: 0.54617   3rd Qu.: 0.6858  `\
`#>  Max.   :10.0   Max.   :19.900   Max.   : 2.83007   Max.   : 3.3403  `\
`#>        y3         `\
`#>  Min.   :-1.9782  `\
`#>  1st Qu.:-0.3253  `\
`#>  Median : 0.1306  `\
`#>  Mean   : 0.1034  `\
`#>  3rd Qu.: 0.5321  `\
`#>  Max.   : 2.4172`\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-ou-error-ou-1.png)![](fig-vignettes-ou-error-ou-2.png)![](fig-vignettes-ou-error-ou-3.png)

## Model Fitting

We fit the model with `dynr`, which estimates state space models by
maximum likelihood. For continuous-time models, `dynr` uses the time
between measurements to get the implied discrete-time model for each
interval, then applies the Kalman filter.

A `dynr` model is built from pieces that match the equations above: the
data, the initial condition, the measurement model, the dynamic model,
and the noise. Each parameter gets a name (for example, `phi_2_1` for
\\(\phi\_{21}\\)\`), and the same name used in two places means one
shared parameter. The starting values are set to the population values
for convenience; with real data, they would be rough guesses.

### Prepare Data

Tell `dynr` which columns hold the individual, the time, and the
observed variables.

\
`dynr_data`` ``<-`` ``dynr``::`[`dynr.data`](https://rdrr.io/pkg/dynr/man/dynr.data.html)`(`\
`  dataframe ``=`` ``data``,`\
`  id ``=`` ``"id"``,`\
`  time ``=`` ``"time"``,`\
`  observed ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``, ``"y3"``)`\
`)`

### Prepare Initial Condition

The mean \\(\boldsymbol{\mu}\_{\boldsymbol{\eta} \mid 0}\\) and
covariance \\(\boldsymbol{\Sigma}\_{\boldsymbol{\eta} \mid 0}\\) of the
initial state are estimated. The covariance matrix is symmetric, so
`sigma0_2_1` appears both above and below the diagonal.

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

### Prepare Measurement Model

The loadings \\(\boldsymbol{\Lambda}\\) are an identity matrix with no
free parameters (`"fixed"`).

\
`dynr_measurement`` ``<-`` ``dynr``::`[`prep.measurement`](https://rdrr.io/pkg/dynr/man/prep.measurement.html)`(`\
`  values.load ``=`` `[`diag`](https://rdrr.io/r/base/diag.html)`(``3``)``,`\
`  params.load ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``3``, ncol ``=`` ``3``)``,`\
`  state.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta_1"``, ``"eta_2"``, ``"eta_3"``)``,`\
`  obs.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``, ``"y3"``)`\
`)`

### Prepare Dynamic Process

Each formula is the rate of change of one latent variable: the effects
\\(\phi\_{jk}\\) of the distances of the three latent variables from
their equilibria \\(\mu\_{k}\\). `isContinuousTime = TRUE` makes this a
continuous-time model.

\
`dynr_dynamics`` ``<-`` ``dynr``::`[`prep.formulaDynamics`](https://rdrr.io/pkg/dynr/man/prep.formulaDynamics.html)`(`\
`  formula ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``  `\
`    ``eta_1`` ``~`` ``(``phi_1_1`` ``*`` ``(``eta_1`` ``-`` ``mu_1_1``)``)`` ``+`` ``(``phi_1_2`` ``*`` ``(``eta_2`` ``-`` ``mu_2_1``)``)`` ``+`` ``(``phi_1_3`` ``*`` ``(``eta_3`` ``-`` ``mu_3_1``)``)``,`\
`    ``eta_2`` ``~`` ``(``phi_2_1`` ``*`` ``(``eta_1`` ``-`` ``mu_1_1``)``)`` ``+`` ``(``phi_2_2`` ``*`` ``(``eta_2`` ``-`` ``mu_2_1``)``)`` ``+`` ``(``phi_2_3`` ``*`` ``(``eta_3`` ``-`` ``mu_3_1``)``)``,`\
`    ``eta_3`` ``~`` ``(``phi_3_1`` ``*`` ``(``eta_1`` ``-`` ``mu_1_1``)``)`` ``+`` ``(``phi_3_2`` ``*`` ``(``eta_2`` ``-`` ``mu_2_1``)``)`` ``+`` ``(``phi_3_3`` ``*`` ``(``eta_3`` ``-`` ``mu_3_1``)``)`\
`  ``)``,`\
`  startval ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    mu_1_1 ``=`` ``mu``[``1``]``, mu_2_1 ``=`` ``mu``[``2``]``, mu_3_1 ``=`` ``mu``[``3``]``,`\
`    phi_1_1 ``=`` ``phi``[``1``, ``1``]``, phi_1_2 ``=`` ``phi``[``1``, ``2``]``, phi_1_3 ``=`` ``phi``[``1``, ``3``]``,`\
`    phi_2_1 ``=`` ``phi``[``2``, ``1``]``, phi_2_2 ``=`` ``phi``[``2``, ``2``]``, phi_2_3 ``=`` ``phi``[``2``, ``3``]``,`\
`    phi_3_1 ``=`` ``phi``[``3``, ``1``]``, phi_3_2 ``=`` ``phi``[``3``, ``2``]``, phi_3_3 ``=`` ``phi``[``3``, ``3``]`\
`  ``)``,`\
`  isContinuousTime ``=`` ``TRUE`\
`)`

### Prepare Process Noise

The diffusion covariance \\(\boldsymbol{\Sigma}\\) is estimated in full.
The measurement error covariance \\(\boldsymbol{\Theta}\\) is diagonal:
the three error variances are estimated, and the error covariances are
fixed at zero.

\
`dynr_noise`` ``<-`` ``dynr``::`[`prep.noise`](https://rdrr.io/pkg/dynr/man/prep.noise.html)`(`\
`  values.latent ``=`` ``sigma``,`\
`  params.latent ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"sigma_1_1"``, ``"sigma_2_1"``, ``"sigma_3_1"``,`\
`      ``"sigma_2_1"``, ``"sigma_2_2"``, ``"sigma_3_2"``,`\
`      ``"sigma_3_1"``, ``"sigma_3_2"``, ``"sigma_3_3"`\
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

### Prepare the Model

The pieces are put together into one model. `dynr` writes the model as C
code to `outfile` and compiles it.

\
`model`` ``<-`` ``dynr``::`[`dynr.model`](https://rdrr.io/pkg/dynr/man/dynr.model.html)`(`\
`  data ``=`` ``dynr_data``,`\
`  initial ``=`` ``dynr_initial``,`\
`  measurement ``=`` ``dynr_measurement``,`\
`  dynamics ``=`` ``dynr_dynamics``,`\
`  noise ``=`` ``dynr_noise``,`\
`  outfile ``=`` ``"ou.c"`\
`)`

Add lower and upper bounds on the elements of \\(\boldsymbol{\Phi}\\) to
keep the optimizer away from implausible values.

\
`model``$``lb``[`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"phi_1_1"``,`\
`    ``"phi_1_2"``,`\
`    ``"phi_1_3"``,`\
`    ``"phi_2_1"``,`\
`    ``"phi_2_2"``,`\
`    ``"phi_2_3"``,`\
`    ``"phi_3_1"``,`\
`    ``"phi_3_2"``,`\
`    ``"phi_3_3"`\
`  ``)`\
`]`` ``<-`` ``-``1.5`\
`model``$``ub``[`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"phi_1_1"``,`\
`    ``"phi_1_2"``,`\
`    ``"phi_1_3"``,`\
`    ``"phi_2_1"``,`\
`    ``"phi_2_2"``,`\
`    ``"phi_2_3"``,`\
`    ``"phi_3_1"``,`\
`    ``"phi_3_2"``,`\
`    ``"phi_3_3"`\
`  ``)`\
`]`` ``<-`` ``+``1.5`

The plot below shows the model equations with the parameter names.

![](fig-vignettes-ou-model-ou-1.png)

### Fit the Model

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
`#> Original fitted parameters:  0.08283144 0.1453763 0.09458133 -0.7229714 `\
`#> 0.3040494 -0.1753704 0.6566508 -0.4967513 0.04628816 -0.1656435 0.4463284 `\
`#> -0.4893126 -1.260946 0.07872375 -0.2825057 -2.53622 0.1628743 -3.298738 `\
`#> -1.633867 -1.639133 -1.59292 0.1528623 0.08530865 0.06489118 -4.928972 7.666848 `\
`#> 4.019886 -3.325513 -1.39442 -7.401412 `\
`#> `\
`#> Transformed fitted parameters:  0.08283144 0.1453763 0.09458133 -0.7229714 `\
`#> 0.3040494 -0.1753704 0.6566508 -0.4967513 0.04628816 -0.1656435 0.4463284 `\
`#> -0.4893126 0.283386 0.0223092 -0.08005815 0.08092138 0.006591483 0.06164672 `\
`#> 0.1951734 0.1941483 0.203331 0.1528623 0.08530865 0.06489118 0.007233938 `\
`#> 0.0554615 0.02907961 0.4611689 0.1728138 0.1874164 `\
`#> `\
`#> Doing end processing`\
`#> Successful trial`\
`#> Total Time: 46.36572 `\
`#> Backend Time: 45.95015`

## Summary

The summary lists each estimate with its standard error, \\(t\\) value,
confidence interval, and \\(p\\) value.

\
[`summary`](https://rdrr.io/r/base/summary.html)`(``results``)`\
`#> Coefficients:`\
`#>             Estimate Std. Error t value  ci.lower  ci.upper            Pr(>|t|)`\
`#> mu_1_1      0.082831   0.096221   0.861 -0.105758  0.271421              0.1947`\
`#> mu_2_1      0.145376   0.152378   0.954 -0.153280  0.444032              0.1701`\
`#> mu_3_1      0.094581   0.096519   0.980 -0.094593  0.283755              0.1636`\
`#> phi_1_1    -0.722971   0.318136  -2.273 -1.346506 -0.099437              0.0116`\
`#> phi_1_2     0.304049   0.281458   1.080 -0.247598  0.855697              0.1401`\
`#> phi_1_3    -0.175370   0.224581  -0.781 -0.615540  0.264800              0.2175`\
`#> phi_2_1     0.656651   0.196747   3.338  0.271034  1.042268              0.0004`\
`#> phi_2_2    -0.496751   0.176656  -2.812 -0.842990 -0.150512              0.0025`\
`#> phi_2_3     0.046288   0.137451   0.337 -0.223112  0.315688              0.3682`\
`#> phi_3_1    -0.165643   0.176727  -0.937 -0.512022  0.180735              0.1744`\
`#> phi_3_2     0.446328   0.157378   2.836  0.137874  0.754783              0.0023`\
`#> phi_3_3    -0.489313   0.123534  -3.961 -0.731434 -0.247191 <0.0000000000000002`\
`#> sigma_1_1   0.283386   0.057191   4.955  0.171293  0.395479 <0.0000000000000002`\
`#> sigma_2_1   0.022309   0.022420   0.995 -0.021634  0.066252              0.1599`\
`#> sigma_3_1  -0.080058   0.022411  -3.572 -0.123982 -0.036134              0.0002`\
`#> sigma_2_2   0.080921   0.016843   4.804  0.047909  0.113934 <0.0000000000000002`\
`#> sigma_3_2   0.006591   0.010786   0.611 -0.014549  0.027732              0.2706`\
`#> sigma_3_3   0.061647   0.014083   4.377  0.034045  0.089249 <0.0000000000000002`\
`#> theta_1_1   0.195173   0.008312  23.480  0.178882  0.211465 <0.0000000000000002`\
`#> theta_2_2   0.194148   0.007035  27.596  0.180359  0.207937 <0.0000000000000002`\
`#> theta_3_3   0.203331   0.007208  28.207  0.189203  0.217459 <0.0000000000000002`\
`#> mu0_1_1     0.152862   0.089210   1.714 -0.021986  0.327710              0.0434`\
`#> mu0_2_1     0.085309   0.225575   0.378 -0.356811  0.527428              0.3527`\
`#> mu0_3_1     0.064891   0.151637   0.428 -0.232311  0.362093              0.3344`\
`#> sigma0_1_1  0.007234   0.015914   0.455 -0.023957  0.038425              0.3247`\
`#> sigma0_2_1  0.055462   0.063327   0.876 -0.068657  0.179580              0.1906`\
`#> sigma0_3_1  0.029080   0.038759   0.750 -0.046886  0.105045              0.2266`\
`#> sigma0_2_2  0.461169   0.228776   2.016  0.012776  0.909562              0.0220`\
`#> sigma0_3_2  0.172814   0.121593   1.421 -0.065503  0.411131              0.0777`\
`#> sigma0_3_3  0.187416   0.102896   1.821 -0.014256  0.389089              0.0343`\
`#>               `\
`#> mu_1_1        `\
`#> mu_2_1        `\
`#> mu_3_1        `\
`#> phi_1_1    *  `\
`#> phi_1_2       `\
`#> phi_1_3       `\
`#> phi_2_1    ***`\
`#> phi_2_2    ** `\
`#> phi_2_3       `\
`#> phi_3_1       `\
`#> phi_3_2    ** `\
`#> phi_3_3    ***`\
`#> sigma_1_1  ***`\
`#> sigma_2_1     `\
`#> sigma_3_1  ***`\
`#> sigma_2_2  ***`\
`#> sigma_3_2     `\
`#> sigma_3_3  ***`\
`#> theta_1_1  ***`\
`#> theta_2_2  ***`\
`#> theta_3_3  ***`\
`#> mu0_1_1    *  `\
`#> mu0_2_1       `\
`#> mu0_3_1       `\
`#> sigma0_1_1    `\
`#> sigma0_2_1    `\
`#> sigma0_3_1    `\
`#> sigma0_2_2 *  `\
`#> sigma0_3_2 .  `\
`#> sigma0_3_3 *  `\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> -2 log-likelihood value at convergence = 8521.43`\
`#> AIC = 8581.43`\
`#> BIC = 8749.46`

### Parameter Estimates

The estimates arranged as the matrices in the model, to compare with the
population values in the Data Generation section. The last matrix,
`beta_var1_hat`, is \\(\exp \left( \hat{\boldsymbol{\Phi}} \right)\\),
the implied discrete-time transition matrix over one time unit.

\
`mu_hat`\
`#> [1] 0.08283144 0.14537631 0.09458133`\
`phi_hat`\
`#>            [,1]       [,2]        [,3]`\
`#> [1,] -0.7229714  0.3040494 -0.17537040`\
`#> [2,]  0.6566508 -0.4967513  0.04628816`\
`#> [3,] -0.1656435  0.4463284 -0.48931259`\
`sigma_hat`\
`#>             [,1]        [,2]         [,3]`\
`#> [1,]  0.28338595 0.022309204 -0.080058147`\
`#> [2,]  0.02230920 0.080921380  0.006591483`\
`#> [3,] -0.08005815 0.006591483  0.061646715`\
`theta_hat`\
`#>           [,1]      [,2]     [,3]`\
`#> [1,] 0.1951734 0.0000000 0.000000`\
`#> [2,] 0.0000000 0.1941483 0.000000`\
`#> [3,] 0.0000000 0.0000000 0.203331`\
`mu0_hat`\
`#> [1] 0.15286229 0.08530865 0.06489118`\
`sigma0_hat`\
`#>             [,1]      [,2]       [,3]`\
`#> [1,] 0.007233938 0.0554615 0.02907961`\
`#> [2,] 0.055461501 0.4611689 0.17281383`\
`#> [3,] 0.029079606 0.1728138 0.18741641`\
`beta_var1_hat`` ``<-`` ``simStateSpace``::`[`ExpMat`](https://github.com/jeksterslab/simStateSpace/reference/ExpMat.md)`(``phi_hat``)`\
`beta_var1_hat`\
`#>              [,1]      [,2]         [,3]`\
`#> [1,]  0.541538515 0.1495829 -0.095629562`\
`#> [2,]  0.369668978 0.6672313 -0.003926669`\
`#> [3,] -0.009502607 0.2683484  0.622663240`

## Table of Parameters and Estimates

The table below puts each population value used to simulate the data
next to its estimate. Bias is the estimate minus the population value.
Because the model is fitted to a single simulated data set, this is the
error of one estimate, not a Monte Carlo bias, which would average over
many replications. The 95% confidence interval is the Wald interval,
\\(\hat{\theta} \pm 1.96 \times \mathrm{SE}\\), and `Covers` shows
whether it contains the population value.

\
`EstimatesTable`` ``<-`` ``function``(``population``, ``estimate``, ``se``)`` ``{`\
`  ``z`` ``<-`` ``stats``::`[`qnorm`](https://rdrr.io/r/stats/Normal.html)`(``p ``=`` ``0.975``)`\
`  ``ll`` ``<-`` ``estimate`` ``-`` ``z`` ``*`` ``se`\
`  ``ul`` ``<-`` ``estimate`` ``+`` ``z`` ``*`` ``se`\
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
`estimates`` ``<-`` `[`coef`](https://rdrr.io/r/stats/coef.html)`(``results``)`\
`se`` ``<-`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(`[`diag`](https://rdrr.io/r/base/diag.html)`(`[`vcov`](https://rdrr.io/r/stats/vcov.html)`(``results``)``)``)`\
`# the population value for each dynr parameter name,`\
`# for example, "beta_1_2" is beta[1, 2] and "mu0_2_1" is mu0[2]`\
`population`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`\
`  X ``=`` `[`strsplit`](https://rdrr.io/r/base/strsplit.html)`(``x ``=`` `[`names`](https://rdrr.io/r/base/names.html)`(``estimates``)``, split ``=`` ``"_"``)``,`\
`  FUN ``=`` ``function``(``x``)`` ``{`\
`    `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(`[`get`](https://rdrr.io/r/base/get.html)`(``x``[``1``]``)``)``[`[`as.integer`](https://rdrr.io/r/base/integer.html)`(``x``[``2``]``)``, `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``x``[``3``]``)``]`\
`  ``}``,`\
`  FUN.VALUE ``=`` `[`numeric`](https://rdrr.io/r/base/numeric.html)`(``1``)`\
`)`\
[`names`](https://rdrr.io/r/base/names.html)`(``population``)`` ``<-`` `[`names`](https://rdrr.io/r/base/names.html)`(``estimates``)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` ``EstimatesTable``(`\
`    population ``=`` ``population``,`\
`    estimate ``=`` ``estimates``,`\
`    se ``=`` ``se`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and dynr estimates (one simulated data set)"`\
`)`

| Parameter  | Population | Estimate |    SE |   Bias |   2.5% |  97.5% | Covers |
|:-----------|-----------:|---------:|------:|-------:|-------:|-------:|:-------|
| mu_1_1     |      0.000 |    0.083 | 0.096 |  0.083 | -0.106 |  0.271 | TRUE   |
| mu_2_1     |      0.000 |    0.145 | 0.152 |  0.145 | -0.153 |  0.444 | TRUE   |
| mu_3_1     |      0.000 |    0.095 | 0.097 |  0.095 | -0.095 |  0.284 | TRUE   |
| phi_1_1    |     -0.357 |   -0.723 | 0.318 | -0.366 | -1.347 | -0.099 | TRUE   |
| phi_1_2    |      0.000 |    0.304 | 0.281 |  0.304 | -0.248 |  0.856 | TRUE   |
| phi_1_3    |      0.000 |   -0.175 | 0.225 | -0.175 | -0.616 |  0.265 | TRUE   |
| phi_2_1    |      0.771 |    0.657 | 0.197 | -0.114 |  0.271 |  1.042 | TRUE   |
| phi_2_2    |     -0.511 |   -0.497 | 0.177 |  0.014 | -0.843 | -0.151 | TRUE   |
| phi_2_3    |      0.000 |    0.046 | 0.137 |  0.046 | -0.223 |  0.316 | TRUE   |
| phi_3_1    |     -0.450 |   -0.166 | 0.177 |  0.284 | -0.512 |  0.181 | TRUE   |
| phi_3_2    |      0.729 |    0.446 | 0.157 | -0.283 |  0.138 |  0.755 | TRUE   |
| phi_3_3    |     -0.693 |   -0.489 | 0.124 |  0.204 | -0.731 | -0.247 | TRUE   |
| sigma_1_1  |      0.245 |    0.283 | 0.057 |  0.039 |  0.171 |  0.395 | TRUE   |
| sigma_2_1  |      0.022 |    0.022 | 0.022 |  0.000 | -0.022 |  0.066 | TRUE   |
| sigma_3_1  |     -0.050 |   -0.080 | 0.022 | -0.030 | -0.124 | -0.036 | TRUE   |
| sigma_2_2  |      0.071 |    0.081 | 0.017 |  0.010 |  0.048 |  0.114 | TRUE   |
| sigma_3_2  |      0.015 |    0.007 | 0.011 | -0.009 | -0.015 |  0.028 | TRUE   |
| sigma_3_3  |      0.076 |    0.062 | 0.014 | -0.014 |  0.034 |  0.089 | TRUE   |
| theta_1_1  |      0.200 |    0.195 | 0.008 | -0.005 |  0.179 |  0.211 | TRUE   |
| theta_2_2  |      0.200 |    0.194 | 0.007 | -0.006 |  0.180 |  0.208 | TRUE   |
| theta_3_3  |      0.200 |    0.203 | 0.007 |  0.003 |  0.189 |  0.217 | TRUE   |
| mu0_1_1    |      0.000 |    0.153 | 0.089 |  0.153 | -0.022 |  0.328 | TRUE   |
| mu0_2_1    |      0.000 |    0.085 | 0.226 |  0.085 | -0.357 |  0.527 | TRUE   |
| mu0_3_1    |      0.000 |    0.065 | 0.152 |  0.065 | -0.232 |  0.362 | TRUE   |
| sigma0_1_1 |      0.343 |    0.007 | 0.016 | -0.335 | -0.024 |  0.038 | FALSE  |
| sigma0_2_1 |      0.330 |    0.055 | 0.063 | -0.274 | -0.069 |  0.180 | FALSE  |
| sigma0_3_1 |      0.034 |    0.029 | 0.039 | -0.005 | -0.047 |  0.105 | TRUE   |
| sigma0_2_2 |      0.566 |    0.461 | 0.229 | -0.105 |  0.013 |  0.910 | TRUE   |
| sigma0_3_2 |      0.255 |    0.173 | 0.122 | -0.082 | -0.066 |  0.411 | TRUE   |
| sigma0_3_3 |      0.300 |    0.187 | 0.103 | -0.113 | -0.014 |  0.389 | TRUE   |

Population values and dynr estimates (one simulated data set) {.table
style="width:99%;"}

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, *17*(2), 303–332.
<https://doi.org/10.1080/10705511003661553>

Chow, S.-M., Losardo, D., Park, J. J., & Molenaar, P. C. M. (2023).
Continuous-time dynamic models: Connections to structural equation
models and other discrete-time models. In R. H. Hoyle (Ed.), *Handbook
of structural equation modeling* (2nd ed.). The Guilford Press.

Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011). A hierarchical
latent stochastic differential equation model for affective dynamics.
*Psychological Methods*, *16*(4), 468–490.
<https://doi.org/10.1037/a0024375>

Ou, L., Hunter, M. D., & Chow, S.-M. (2019). What’s for dynr: A package
for linear and nonlinear dynamic modeling in R. *The R Journal*,
*11*(1), 91. <https://doi.org/10.32614/rj-2019-012>

R Core Team. (2026). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>

Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the
brownian motion. *Physical Review*, *36*(5), 823–841.
<https://doi.org/10.1103/physrev.36.823>
