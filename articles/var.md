# The Vector Autoregressive Model

## Overview

The vector autoregressive (VAR) model describes how a set of variables
measured repeatedly on the same person carry over and influence each
other from one occasion to the next. For example, today’s stress may
predict tomorrow’s stress (carry-over), and today’s stress may also
predict tomorrow’s negative affect (cross-lagged effect).

This tutorial does three things:

1.  It writes the first-order VAR model as a state space model.
2.  It simulates data from the model with
    [`SimSSMVARFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARFixed.md),
    with every individual sharing the same parameters.
3.  It fits the model to the simulated data with the `dynr` and `OpenMx`
    packages and compares both sets of estimates with the values used to
    simulate the data.

No prior knowledge of state space models is needed. A state space model
has two parts: a measurement model that links what we observe to
underlying (latent) states, and a dynamic model that describes how those
states change over time. In the VAR model, the variables are assumed to
be measured without error, so the measurement part is simple, and all
the interest is in the dynamics.

## Model

The measurement model is given by \mathbf{y}\_{i, t} =
\boldsymbol{\eta}\_{i, t} where \mathbf{y}\_{i, t} represents a vector
of observed variables and \boldsymbol{\eta}\_{i, t} a vector of latent
variables for individual i and time t. Since the observed and latent
variables are equal, we only generate data from the dynamic structure.

The dynamic structure is given by \boldsymbol{\eta}\_{i, t} =
\boldsymbol{\alpha} + \boldsymbol{\beta} \boldsymbol{\eta}\_{i, t - 1} +
\boldsymbol{\zeta}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\zeta}\_{i, t} \sim \mathcal{N} \left( \mathbf{0},
\boldsymbol{\Psi} \right) where \boldsymbol{\eta}\_{i, t},
\boldsymbol{\eta}\_{i, t - 1}, and \boldsymbol{\zeta}\_{i, t} are random
variables, and \boldsymbol{\alpha}, \boldsymbol{\beta}, and
\boldsymbol{\Psi} are model parameters. Here, \boldsymbol{\eta}\_{i, t}
is a vector of latent variables at time t and individual i,
\boldsymbol{\eta}\_{i, t - 1} represents a vector of latent variables at
time t - 1 and individual i, and \boldsymbol{\zeta}\_{i, t} represents a
vector of dynamic noise at time t and individual i. \boldsymbol{\alpha}
denotes a vector of intercepts, \boldsymbol{\beta} a matrix of
autoregression and cross regression coefficients, and \boldsymbol{\Psi}
the covariance matrix of \boldsymbol{\zeta}\_{i, t}.

In words, each variable at time t is a weighted sum of all the variables
at time t - 1, plus an intercept and a random shock. The weights are the
elements of \boldsymbol{\beta}:

- A diagonal element \beta\_{jj} is the autoregressive effect, how much
  of variable j carries over to the next occasion. Values closer to 1
  mean the variable changes slowly (high inertia).
- An off-diagonal element \beta\_{jk} is the cross-lagged effect of
  variable k at time t - 1 on variable j at time t. Note the order: the
  row is the outcome, the column is the predictor.

The shocks \boldsymbol{\zeta}\_{i, t} stand for everything that moves
the variables that the model does not explain, for example, events that
happen between occasions. For the process to settle around a stable
mean, all eigenvalues of \boldsymbol{\beta} must have modulus less than
1.

An alternative representation of the dynamic noise is given by
\boldsymbol{\zeta}\_{i, t} = \boldsymbol{\Psi}^{\frac{1}{2}}
\mathbf{z}\_{i, t}, \quad \mathrm{with} \quad \mathbf{z}\_{i, t} \sim
\mathcal{N} \left( \mathbf{0}, \mathbf{I} \right) where \left(
\boldsymbol{\Psi}^{\frac{1}{2}} \right) \left(
\boldsymbol{\Psi}^{\frac{1}{2}} \right)^{\prime} = \boldsymbol{\Psi} .
This is how the noise is simulated: independent standard normal draws
\mathbf{z}\_{i, t} are multiplied by a square root of \boldsymbol{\Psi}.
The functions in `simStateSpace` take this square root, the lower
triangular Cholesky factor of \boldsymbol{\Psi}, as the argument
`psi_l`.

### Diagram

![State space diagram of the VAR model.](fig-vignettes-var-ssm.png)

State space diagram of the VAR model.

The model in state space form, for the first three occasions (t = 0, 1,
2). The observed variables (squares) are the state, \mathbf{y}\_{t} =
\boldsymbol{\eta}\_{t}, so \boldsymbol{\Lambda} is the identity matrix
and there is no measurement error. The horizontal arrows are the
autoregressive effects and the diagonal arrows the cross-lagged effects,
the elements of \boldsymbol{\beta}; only the nonzero elements in this
example are drawn. \zeta\_{j, t} is the process noise, with variance
\psi\_{jj}, from t = 1 on. A loop on a variable is its variance and a
curved two-headed arrow a covariance. At t = 0, the initial condition
\boldsymbol{\eta}\_{0} \sim \mathcal{N} \left( \boldsymbol{\mu}\_{0},
\boldsymbol{\Sigma}\_{0} \right) gives the variances \sigma\_{jj} and
covariances \sigma\_{jk} of the variables.
[`SimSSMVARFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARFixed.md)
takes `beta`, `mu0`, and the Cholesky factors of \boldsymbol{\Psi} and
\boldsymbol{\Sigma}\_{0} (`psi_l` and `sigma0_l`).

## Data Generation

### Notation

Let t = 1000 be the number of time points and n = 100 be the number of
individuals.

Let the initial condition \boldsymbol{\eta}\_{0} be given by

\boldsymbol{\eta}\_{0} \sim \mathcal{N} \left( \boldsymbol{\mu}\_{0},
\boldsymbol{\Sigma}\_{0} \right)

\boldsymbol{\mu}\_{0} = \left( \begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\
\end{array} \right)

\boldsymbol{\Sigma}\_{0} = \left( \begin{array}{ccc} 0.1960784 &
0.1183232 & 0.0298539 \\\\ 0.1183232 & 0.3437711 & 0.1381855 \\\\
0.0298539 & 0.1381855 & 0.2663828 \\\\ \end{array} \right) .

Let the constant vector \boldsymbol{\alpha} be given by

\boldsymbol{\alpha} = \left( \begin{array}{c} 0 \\\\ 0 \\\\ 0 \\\\
\end{array} \right) .

Let the transition matrix \boldsymbol{\beta} be given by

\boldsymbol{\beta} = \left( \begin{array}{ccc} 0.7 & 0 & 0 \\\\ 0.5 &
0.6 & 0 \\\\ -0.1 & 0.4 & 0.5 \\\\ \end{array} \right) .

Read the columns as predictors and the rows as outcomes. Each variable
carries over to the next occasion (autoregressive effects 0.7, 0.6, and
0.5). In addition, \eta\_{1} predicts later \eta\_{2} (0.5), \eta\_{2}
predicts later \eta\_{3} (0.4), and \eta\_{1} has a small negative
effect on later \eta\_{3} (-0.1). Because \boldsymbol{\beta} is lower
triangular, its eigenvalues are its diagonal elements, all less than 1,
so the process is stable.

Let the dynamic process noise \boldsymbol{\Psi} be given by

\boldsymbol{\Psi} = \left( \begin{array}{ccc} 0.1 & 0 & 0 \\\\ 0 & 0.1 &
0 \\\\ 0 & 0 & 0.1 \\\\ \end{array} \right) .

The initial condition describes where each individual starts. Here, the
mean vector and covariance matrix of \boldsymbol{\eta}\_{0} are the
model-implied mean and covariance of the stable process, computed with
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
`#> [1] 1000`\
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
`#> [3,] 0.0000000 0.0000000 0.3162278`

### Visualizing the Dynamics Without Process Noise (n = 5 with Different Initial Condition)

To see what \boldsymbol{\beta} alone does, we first simulate five
individuals with the process noise set to zero. Each individual starts
at a different point, and the trajectories return to the mean at a speed
set by \boldsymbol{\beta}.

![](fig-vignettes-var-no-error-var-1.png)![](fig-vignettes-var-no-error-var-2.png)![](fig-vignettes-var-no-error-var-3.png)

### Using the `SimSSMVARFixed` Function from the `simStateSpace` Package to Simulate Data

Now we add the process noise and simulate the full data set.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the data in long format, one row per individual and time point, with the
columns `id`, `time`, `y1`, `y2`, and `y3`.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`simStateSpace`](https://github.com/jeksterslab/simStateSpace)`)`\
`sim`` ``<-`` `[`SimSSMVARFixed`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARFixed.md)`(`\
`  n ``=`` ``n``,`\
`  time ``=`` ``time``,`\
`  mu0 ``=`` ``mu0``,`\
`  sigma0_l ``=`` ``sigma0_l``,`\
`  alpha ``=`` ``alpha``,`\
`  beta ``=`` ``beta``,`\
`  psi_l ``=`` ``psi_l`\
`)`\
`data`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data``)`\
`#>   id time          y1          y2          y3`\
`#> 1  1    0 -0.81728749  0.01319023  0.57975035`\
`#> 2  1    1 -0.62264147 -0.94877859  0.54422379`\
`#> 3  1    2 -0.06731464 -0.69530903 -0.09626706`\
`#> 4  1    3 -0.09174342 -0.44804477 -0.27983813`\
`#> 5  1    4 -0.04104160 -0.72769549 -0.19238403`\
`#> 6  1    5  0.18907246 -0.81602072 -0.73831208`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``data``)`\
`#>        id              time             y1                   y2            `\
`#>  Min.   :  1.00   Min.   :  0.0   Min.   :-1.9239387   Min.   :-2.6896682  `\
`#>  1st Qu.: 25.75   1st Qu.:249.8   1st Qu.:-0.2991519   1st Qu.:-0.3940385  `\
`#>  Median : 50.50   Median :499.5   Median :-0.0029441   Median : 0.0008316  `\
`#>  Mean   : 50.50   Mean   :499.5   Mean   :-0.0005164   Mean   : 0.0004548  `\
`#>  3rd Qu.: 75.25   3rd Qu.:749.2   3rd Qu.: 0.2957664   3rd Qu.: 0.3976987  `\
`#>  Max.   :100.00   Max.   :999.0   Max.   : 1.9753039   Max.   : 2.3189032  `\
`#>        y3           `\
`#>  Min.   :-2.185100  `\
`#>  1st Qu.:-0.343439  `\
`#>  Median : 0.001899  `\
`#>  Mean   : 0.002713  `\
`#>  3rd Qu.: 0.351819  `\
`#>  Max.   : 2.113439`\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-var-error-var-1.png)![](fig-vignettes-var-error-var-2.png)![](fig-vignettes-var-error-var-3.png)

## Model Fitting

The same simulated data are fitted with two programs that estimate state
space models by maximum likelihood with the Kalman filter: `dynr` and
`OpenMx`. The section Population Values and Estimates compares their
estimates with the population values.

### dynr

We fit the model with `dynr`, which estimates state space models by
maximum likelihood using the Kalman filter. A `dynr` model is built from
pieces that match the equations above: the data, the initial condition,
the measurement model, the dynamic model, and the noise. Each parameter
gets a name (for example, `beta_2_1` for \beta\_{21}), and the same name
used in two places means one shared parameter. The starting values are
set to the population values for convenience; with real data, they would
be rough guesses.

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

The mean \boldsymbol{\mu}\_{0} and covariance \boldsymbol{\Sigma}\_{0}
of the initial state are estimated. The covariance matrix is symmetric,
so `sigma0_2_1` appears both above and below the diagonal.

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

The measurement model \mathbf{y}\_{i, t} = \boldsymbol{\eta}\_{i, t} is
an identity matrix of loadings, with no free parameters (`"fixed"`).

\
`dynr_measurement`` ``<-`` ``dynr``::`[`prep.measurement`](https://rdrr.io/pkg/dynr/man/prep.measurement.html)`(`\
`  values.load ``=`` `[`diag`](https://rdrr.io/r/base/diag.html)`(``3``)``,`\
`  params.load ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``3``, ncol ``=`` ``3``)``,`\
`  state.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta_1"``, ``"eta_2"``, ``"eta_3"``)``,`\
`  obs.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``, ``"y3"``)`\
`)`

#### Prepare Dynamic Process

Each formula is one row of the dynamic equation: the intercept
\alpha\_{j} plus the effects \beta\_{jk} of the three variables at the
previous occasion. `isContinuousTime = FALSE` makes this a discrete-time
model.

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

The latent noise covariance \boldsymbol{\Psi} is estimated. The
measurement error covariance is fixed at zero, because the variables are
measured without error in this model.

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
`  values.observed ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``0``, nrow ``=`` ``3``, ncol ``=`` ``3``)``,`\
`  params.observed ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``3``, ncol ``=`` ``3``)`\
`)`

#### Prepare the Model

The pieces are put together into one model. `dynr` writes the model as C
code to `outfile` and compiles it. The plot below shows the model
equations with the parameter names.

\
`model`` ``<-`` ``dynr``::`[`dynr.model`](https://rdrr.io/pkg/dynr/man/dynr.model.html)`(`\
`  data ``=`` ``dynr_data``,`\
`  initial ``=`` ``dynr_initial``,`\
`  measurement ``=`` ``dynr_measurement``,`\
`  dynamics ``=`` ``dynr_dynamics``,`\
`  noise ``=`` ``dynr_noise``,`\
`  outfile ``=`` ``"var.c"`\
`)`

![](fig-vignettes-var-model-var-1.png)

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
`#> Original fitted parameters:  -0.0002106182 0.000261195 0.0009243715 0.6972056 `\
`#> 0.0006123385 0.002917166 0.4995238 0.6005427 -0.00177082 -0.1015684 0.4031293 `\
`#> 0.4972227 -2.301822 -0.003746233 0.001270917 -2.306638 0.002932055 -2.304454 `\
`#> -0.034054 0.080811 0.1345955 -1.449063 0.6212102 0.0249533 -1.163541 0.3906055 `\
`#> -1.54021 `\
`#> `\
`#> Transformed fitted parameters:  -0.0002106182 0.000261195 0.0009243715 `\
`#> 0.6972056 0.0006123385 0.002917166 0.4995238 0.6005427 -0.00177082 -0.1015684 `\
`#> 0.4031293 0.4972227 0.1000764 -0.0003749095 0.0001271888 0.09959695 `\
`#> 0.0002915431 0.09981433 -0.034054 0.080811 0.1345955 0.2347903 0.1458541 `\
`#> 0.005858792 0.4029842 0.1256562 0.2621428 `\
`#> `\
`#> Doing end processing`\
`#> Successful trial`\
`#> Total Time: 20.71834 `\
`#> Backend Time: 20.70897`

#### Summary

The summary lists each estimate with its standard error, t value,
confidence interval, and p value.

\
[`summary`](https://rdrr.io/r/base/summary.html)`(``results``)`\
`#> Coefficients:`\
`#>              Estimate Std. Error t value   ci.lower   ci.upper`\
`#> alpha_1_1  -0.0002106  0.0010009  -0.210 -0.0021723  0.0017511`\
`#> alpha_2_1   0.0002612  0.0009985   0.262 -0.0016958  0.0022182`\
`#> alpha_3_1   0.0009244  0.0009996   0.925 -0.0010348  0.0028835`\
`#> beta_1_1    0.6972056  0.0025522 273.182  0.6922034  0.7022077`\
`#> beta_1_2    0.0006123  0.0021517   0.285 -0.0036048  0.0048295`\
`#> beta_1_3    0.0029172  0.0021942   1.329 -0.0013834  0.0072177`\
`#> beta_2_1    0.4995238  0.0025463 196.175  0.4945331  0.5045145`\
`#> beta_2_2    0.6005427  0.0021466 279.758  0.5963353  0.6047500`\
`#> beta_2_3   -0.0017708  0.0021888  -0.809 -0.0060609  0.0025192`\
`#> beta_3_1   -0.1015684  0.0025490 -39.846 -0.1065644 -0.0965724`\
`#> beta_3_2    0.4031293  0.0021489 187.595  0.3989175  0.4073411`\
`#> beta_3_3    0.4972227  0.0021913 226.912  0.4929279  0.5015175`\
`#> psi_1_1     0.1000764  0.0004477 223.510  0.0991988  0.1009539`\
`#> psi_2_1    -0.0003749  0.0003159  -1.187 -0.0009940  0.0002441`\
`#> psi_3_1     0.0001272  0.0003162   0.402 -0.0004925  0.0007469`\
`#> psi_2_2     0.0995969  0.0004456 223.504  0.0987236  0.1004703`\
`#> psi_3_2     0.0002915  0.0003154   0.924 -0.0003267  0.0009098`\
`#> psi_3_3     0.0998143  0.0004466 223.503  0.0989390  0.1006896`\
`#> mu0_1_1    -0.0340540  0.0484225  -0.703 -0.1289603  0.0608523`\
`#> mu0_2_1     0.0808110  0.0635000   1.273 -0.0436468  0.2052688`\
`#> mu0_3_1     0.1345955  0.0511762   2.630  0.0342920  0.2348990`\
`#> sigma0_1_1  0.2347903  0.0332779   7.055  0.1695667  0.3000138`\
`#> sigma0_2_1  0.1458541  0.0340614   4.282  0.0790950  0.2126132`\
`#> sigma0_3_1  0.0058588  0.0248510   0.236 -0.0428483  0.0545659`\
`#> sigma0_2_2  0.4029842  0.0571383   7.053  0.2909952  0.5149731`\
`#> sigma0_3_2  0.1256562  0.0349783   3.592  0.0571000  0.1942124`\
`#> sigma0_3_3  0.2621428  0.0371329   7.060  0.1893635  0.3349220`\
`#>                       Pr(>|t|)    `\
`#> alpha_1_1               0.4167    `\
`#> alpha_2_1               0.3968    `\
`#> alpha_3_1               0.1775    `\
`#> beta_1_1   <0.0000000000000002 ***`\
`#> beta_1_2                0.3880    `\
`#> beta_1_3                0.0918 .  `\
`#> beta_2_1   <0.0000000000000002 ***`\
`#> beta_2_2   <0.0000000000000002 ***`\
`#> beta_2_3                0.2093    `\
`#> beta_3_1   <0.0000000000000002 ***`\
`#> beta_3_2   <0.0000000000000002 ***`\
`#> beta_3_3   <0.0000000000000002 ***`\
`#> psi_1_1    <0.0000000000000002 ***`\
`#> psi_2_1                 0.1176    `\
`#> psi_3_1                 0.3437    `\
`#> psi_2_2    <0.0000000000000002 ***`\
`#> psi_3_2                 0.1777    `\
`#> psi_3_3    <0.0000000000000002 ***`\
`#> mu0_1_1                 0.2409    `\
`#> mu0_2_1                 0.1016    `\
`#> mu0_3_1                 0.0043 ** `\
`#> sigma0_1_1 <0.0000000000000002 ***`\
`#> sigma0_2_1 <0.0000000000000002 ***`\
`#> sigma0_3_1              0.4068    `\
`#> sigma0_2_2 <0.0000000000000002 ***`\
`#> sigma0_3_2              0.0002 ***`\
`#> sigma0_3_3 <0.0000000000000002 ***`\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> -2 log-likelihood value at convergence = 160369.89`\
`#> AIC = 160423.89`\
`#> BIC = 160680.73`

### OpenMx: State Space Model

`OpenMx` (Neale et al., 2016) fits the same state space model with its
own Kalman filter
([`OpenMx::mxExpectationStateSpace()`](https://rdrr.io/pkg/OpenMx/man/mxExpectationStateSpace.html);
Hunter (2018)). Its matrices have the names used in engineering: `A` is
the transition matrix \boldsymbol{\beta}, `B` holds the intercepts of
the dynamic model, `C` is the loading matrix \boldsymbol{\Lambda}, `D`
holds the intercepts of the measurement model, `Q` is the process noise
covariance matrix, and `R` is the measurement error covariance matrix.
The parameter labels are the dynr names; as in dynr, a label used in two
places is one shared parameter, and the starting values are the
population values.

Two details differ from dynr. First, OpenMx reads the data as a single
time series. To start the filter again at the first occasion of every
individual, each individual gets a group, and all groups use the
matrices of one parent model, so every individual has the same
parameters. Second, OpenMx starts the filter one occasion before the
first row, so the first row has the state mean \mathbf{A}
\mathbf{x}\_{0} + \mathbf{B} \mathbf{u} and the state covariance
\mathbf{A} \mathbf{P}\_{0} \mathbf{A}^{\prime} + \mathbf{Q}, with the
constant input \mathbf{u} = 1. Setting \mathbf{x}\_{0} = \mathbf{A}^{-1}
\left( \boldsymbol{\mu}\_{0} - \mathbf{B} \mathbf{u} \right) and
\mathbf{P}\_{0} = \mathbf{A}^{-1} \left( \boldsymbol{\Sigma}\_{0} -
\mathbf{Q} \right) \left( \mathbf{A}^{-1} \right)^{\prime} gives the
first occasion the initial condition of the model, so the likelihood is
the same as in dynr. This needs an invertible \mathbf{A}, which holds
here.

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

The matrices of the model, with the dynr labels. The loadings are an
identity matrix, and there is no measurement error (`R` is zero).

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
`  ``MxMat``(``values ``=`` `[`diag`](https://rdrr.io/r/base/diag.html)`(``3``)``, labels ``=`` ``"fixed"``, name ``=`` ``"C"``)``,`\
`  ``MxMat``(``values ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``x ``=`` ``0``, times ``=`` ``3``)``, labels ``=`` ``"fixed"``, name ``=`` ``"D"``)``,`\
`  ``MxMat``(``values ``=`` ``psi``, labels ``=`` ``labels_psi``, name ``=`` ``"Q"``)``,`\
`  ``MxMat``(`\
`    values ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``0``, nrow ``=`` ``3``, ncol ``=`` ``3``)``,`\
`    labels ``=`` ``"fixed"``,`\
`    name ``=`` ``"R"`\
`  ``)``,`\
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
  name = "var"
)
#> 
Beginning initial fit attempt
Fit attempt 0, fit=160369.869850976, new current best! (was 160395.003347781)
                                                                             
#> 
#>  Solution found!  Final fit=160369.87 (started at 160395)  (1 attempt(s): 1 valid, 0 errors)
summary(fit_mx_ssm)
#> Summary of var 
#>  
#> free parameters:
#>          name matrix row col      Estimate    Std.Error A
#> 1    beta_1_1      A   1   1  0.6972237920 0.0025525563  
#> 2    beta_2_1      A   2   1  0.4995472386 0.0025464048  
#> 3    beta_3_1      A   3   1 -0.1015496972 0.0025491922  
#> 4    beta_1_2      A   1   2  0.0006471480 0.0021518767  
#> 5    beta_2_2      A   2   2  0.6005728491 0.0021466896  
#> 6    beta_3_2      A   3   2  0.4031599054 0.0021490392  
#> 7    beta_1_3      A   1   3  0.0029432067 0.0021943277  
#> 8    beta_2_3      A   2   3 -0.0017459428 0.0021890309  
#> 9    beta_3_3      A   3   3  0.4972480365 0.0021914301  
#> 10  alpha_1_1      B   1   1 -0.0001643356 0.0010009442  
#> 11  alpha_2_1      B   2   1  0.0003127692 0.0009985290  
#> 12  alpha_3_1      B   3   1  0.0009726495 0.0009996220  
#> 13    psi_1_1      Q   1   1  0.1000854780 0.0004478198  
#> 14    psi_2_1      Q   2   1 -0.0003684564 0.0003158945  
#> 15    psi_3_1      Q   3   1  0.0001314732 0.0003162384  
#> 16    psi_2_2      Q   2   2  0.0996030177 0.0004456611  
#> 17    psi_3_2      Q   3   2  0.0002954868 0.0003154764  
#> 18    psi_3_3      Q   3   3  0.0998212256 0.0004466375  
#> 19    mu0_1_1    mu0   1   1 -0.0337610679 0.0484615187  
#> 20    mu0_2_1    mu0   2   1  0.0807387686 0.0635130381  
#> 21    mu0_3_1    mu0   3   1  0.1344571751 0.0512344181  
#> 22 sigma0_1_1 sigma0   1   1  0.2346746474 0.0331935664  
#> 23 sigma0_2_1 sigma0   2   1  0.1457558232 0.0340386288  
#> 24 sigma0_3_1 sigma0   3   1  0.0057781197 0.0248193988  
#> 25 sigma0_2_2 sigma0   2   2  0.4028205457 0.0569862803  
#> 26 sigma0_3_2 sigma0   3   2  0.1256301601 0.0348578644  
#> 27 sigma0_3_3 sigma0   3   3  0.2622837822 0.0371009366  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:             27                 299973              160369.9
#>    Saturated:             NA                     NA                    NA
#> Independence:             NA                     NA                    NA
#> Number of observations/statistics: 100000/300000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:      -439576.1               160423.9                 160423.9
#> BIC:     -3293196.9               160680.7                 160594.9
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-11-05 07:52:43 
#> Wall clock time: 503.6521 secs 
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
is the Wald interval, \hat{\theta} \pm 1.96 \times \mathrm{SE}, unless
the caption says otherwise, and `Covers` shows whether it contains the
population value. The parameter names are the dynr names in every table.

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
| alpha_1_1  |      0.000 |    0.000 | 0.001 |  0.000 | -0.002 |  0.002 | TRUE   |
| alpha_2_1  |      0.000 |    0.000 | 0.001 |  0.000 | -0.002 |  0.002 | TRUE   |
| alpha_3_1  |      0.000 |    0.001 | 0.001 |  0.001 | -0.001 |  0.003 | TRUE   |
| beta_1_1   |      0.700 |    0.697 | 0.003 | -0.003 |  0.692 |  0.702 | TRUE   |
| beta_1_2   |      0.000 |    0.001 | 0.002 |  0.001 | -0.004 |  0.005 | TRUE   |
| beta_1_3   |      0.000 |    0.003 | 0.002 |  0.003 | -0.001 |  0.007 | TRUE   |
| beta_2_1   |      0.500 |    0.500 | 0.003 |  0.000 |  0.495 |  0.505 | TRUE   |
| beta_2_2   |      0.600 |    0.601 | 0.002 |  0.001 |  0.596 |  0.605 | TRUE   |
| beta_2_3   |      0.000 |   -0.002 | 0.002 | -0.002 | -0.006 |  0.003 | TRUE   |
| beta_3_1   |     -0.100 |   -0.102 | 0.003 | -0.002 | -0.107 | -0.097 | TRUE   |
| beta_3_2   |      0.400 |    0.403 | 0.002 |  0.003 |  0.399 |  0.407 | TRUE   |
| beta_3_3   |      0.500 |    0.497 | 0.002 | -0.003 |  0.493 |  0.502 | TRUE   |
| psi_1_1    |      0.100 |    0.100 | 0.000 |  0.000 |  0.099 |  0.101 | TRUE   |
| psi_2_1    |      0.000 |    0.000 | 0.000 |  0.000 | -0.001 |  0.000 | TRUE   |
| psi_3_1    |      0.000 |    0.000 | 0.000 |  0.000 |  0.000 |  0.001 | TRUE   |
| psi_2_2    |      0.100 |    0.100 | 0.000 |  0.000 |  0.099 |  0.100 | TRUE   |
| psi_3_2    |      0.000 |    0.000 | 0.000 |  0.000 |  0.000 |  0.001 | TRUE   |
| psi_3_3    |      0.100 |    0.100 | 0.000 |  0.000 |  0.099 |  0.101 | TRUE   |
| mu0_1_1    |      0.000 |   -0.034 | 0.048 | -0.034 | -0.129 |  0.061 | TRUE   |
| mu0_2_1    |      0.000 |    0.081 | 0.064 |  0.081 | -0.044 |  0.205 | TRUE   |
| mu0_3_1    |      0.000 |    0.135 | 0.051 |  0.135 |  0.034 |  0.235 | FALSE  |
| sigma0_1_1 |      0.196 |    0.235 | 0.033 |  0.039 |  0.170 |  0.300 | TRUE   |
| sigma0_2_1 |      0.118 |    0.146 | 0.034 |  0.028 |  0.079 |  0.213 | TRUE   |
| sigma0_3_1 |      0.030 |    0.006 | 0.025 | -0.024 | -0.043 |  0.055 | TRUE   |
| sigma0_2_2 |      0.344 |    0.403 | 0.057 |  0.059 |  0.291 |  0.515 | TRUE   |
| sigma0_3_2 |      0.138 |    0.126 | 0.035 | -0.013 |  0.057 |  0.194 | TRUE   |
| sigma0_3_3 |      0.266 |    0.262 | 0.037 | -0.004 |  0.189 |  0.335 | TRUE   |

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
| alpha_1_1  |      0.000 |    0.000 | 0.001 |  0.000 | -0.002 |  0.002 | TRUE   |
| alpha_2_1  |      0.000 |    0.000 | 0.001 |  0.000 | -0.002 |  0.002 | TRUE   |
| alpha_3_1  |      0.000 |    0.001 | 0.001 |  0.001 | -0.001 |  0.003 | TRUE   |
| beta_1_1   |      0.700 |    0.697 | 0.003 | -0.003 |  0.692 |  0.702 | TRUE   |
| beta_1_2   |      0.000 |    0.001 | 0.002 |  0.001 | -0.004 |  0.005 | TRUE   |
| beta_1_3   |      0.000 |    0.003 | 0.002 |  0.003 | -0.001 |  0.007 | TRUE   |
| beta_2_1   |      0.500 |    0.500 | 0.003 |  0.000 |  0.495 |  0.505 | TRUE   |
| beta_2_2   |      0.600 |    0.601 | 0.002 |  0.001 |  0.596 |  0.605 | TRUE   |
| beta_2_3   |      0.000 |   -0.002 | 0.002 | -0.002 | -0.006 |  0.003 | TRUE   |
| beta_3_1   |     -0.100 |   -0.102 | 0.003 | -0.002 | -0.107 | -0.097 | TRUE   |
| beta_3_2   |      0.400 |    0.403 | 0.002 |  0.003 |  0.399 |  0.407 | TRUE   |
| beta_3_3   |      0.500 |    0.497 | 0.002 | -0.003 |  0.493 |  0.502 | TRUE   |
| psi_1_1    |      0.100 |    0.100 | 0.000 |  0.000 |  0.099 |  0.101 | TRUE   |
| psi_2_1    |      0.000 |    0.000 | 0.000 |  0.000 | -0.001 |  0.000 | TRUE   |
| psi_3_1    |      0.000 |    0.000 | 0.000 |  0.000 |  0.000 |  0.001 | TRUE   |
| psi_2_2    |      0.100 |    0.100 | 0.000 |  0.000 |  0.099 |  0.100 | TRUE   |
| psi_3_2    |      0.000 |    0.000 | 0.000 |  0.000 |  0.000 |  0.001 | TRUE   |
| psi_3_3    |      0.100 |    0.100 | 0.000 |  0.000 |  0.099 |  0.101 | TRUE   |
| mu0_1_1    |      0.000 |   -0.034 | 0.048 | -0.034 | -0.129 |  0.061 | TRUE   |
| mu0_2_1    |      0.000 |    0.081 | 0.064 |  0.081 | -0.044 |  0.205 | TRUE   |
| mu0_3_1    |      0.000 |    0.134 | 0.051 |  0.134 |  0.034 |  0.235 | FALSE  |
| sigma0_1_1 |      0.196 |    0.235 | 0.033 |  0.039 |  0.170 |  0.300 | TRUE   |
| sigma0_2_1 |      0.118 |    0.146 | 0.034 |  0.027 |  0.079 |  0.212 | TRUE   |
| sigma0_3_1 |      0.030 |    0.006 | 0.025 | -0.024 | -0.043 |  0.054 | TRUE   |
| sigma0_2_2 |      0.344 |    0.403 | 0.057 |  0.059 |  0.291 |  0.515 | TRUE   |
| sigma0_3_2 |      0.138 |    0.126 | 0.035 | -0.013 |  0.057 |  0.194 | TRUE   |
| sigma0_3_3 |      0.266 |    0.262 | 0.037 | -0.004 |  0.190 |  0.335 | TRUE   |

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
| alpha_1_1  |      0.000 |  0.000 |        0.000 |
| alpha_2_1  |      0.000 |  0.000 |        0.000 |
| alpha_3_1  |      0.000 |  0.001 |        0.001 |
| beta_1_1   |      0.700 |  0.697 |        0.697 |
| beta_1_2   |      0.000 |  0.001 |        0.001 |
| beta_1_3   |      0.000 |  0.003 |        0.003 |
| beta_2_1   |      0.500 |  0.500 |        0.500 |
| beta_2_2   |      0.600 |  0.601 |        0.601 |
| beta_2_3   |      0.000 | -0.002 |       -0.002 |
| beta_3_1   |     -0.100 | -0.102 |       -0.102 |
| beta_3_2   |      0.400 |  0.403 |        0.403 |
| beta_3_3   |      0.500 |  0.497 |        0.497 |
| psi_1_1    |      0.100 |  0.100 |        0.100 |
| psi_2_1    |      0.000 |  0.000 |        0.000 |
| psi_3_1    |      0.000 |  0.000 |        0.000 |
| psi_2_2    |      0.100 |  0.100 |        0.100 |
| psi_3_2    |      0.000 |  0.000 |        0.000 |
| psi_3_3    |      0.100 |  0.100 |        0.100 |
| mu0_1_1    |      0.000 | -0.034 |       -0.034 |
| mu0_2_1    |      0.000 |  0.081 |        0.081 |
| mu0_3_1    |      0.000 |  0.135 |        0.134 |
| sigma0_1_1 |      0.196 |  0.235 |        0.235 |
| sigma0_2_1 |      0.118 |  0.146 |        0.146 |
| sigma0_3_1 |      0.030 |  0.006 |        0.006 |
| sigma0_2_2 |      0.344 |  0.403 |        0.403 |
| sigma0_3_2 |      0.138 |  0.126 |        0.126 |
| sigma0_3_3 |      0.266 |  0.262 |        0.262 |

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
