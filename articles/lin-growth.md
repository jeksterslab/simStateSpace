# The Linear Growth Curve Model

## Overview

The linear growth curve model describes individuals who each follow
their own straight-line trajectory over time, for example, reading
scores that start at different levels and grow at different rates. Each
individual has a latent intercept (the level at the first occasion) and
a latent slope (the change per occasion). The model estimates the
average intercept and slope, how much individuals differ in them, and
how much the measurements scatter around each trajectory.

The growth curve model is usually fitted as a structural equation model
(SEM). It can also be written as a state space model: the intercept and
slope are latent states that are set at the first occasion and then
follow a fixed rule.

This tutorial does three things:

1.  It writes the linear growth curve model as a state space model.
2.  It simulates data from the model with
    [`SimSSMLinGrowth()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowth.md).
3.  It fits the model to the simulated data as a state space model with
    `dynr` and `OpenMx` and as a structural equation model (SEM) with
    `lavaan` and `OpenMx`, and as a multilevel model with `lme4`, and
    compares every set of estimates with the values used to simulate the
    data.

No prior knowledge of state space models is needed.

## Model

The measurement model is given by \\\[Y\_{i, t} = \left(
\begin{array}{cc} 1 & 0 \\\\ \end{array} \right) \left( \begin{array}{c}
\eta\_{0\_{i, t}} \\\\ \eta\_{1\_{i, t}} \\\\ \end{array} \right) +
\boldsymbol{\varepsilon}\_{i, t}, \quad \mathrm{with} \quad
\boldsymbol{\varepsilon}\_{i, t} \sim \mathcal{N} \left( 0, \theta
\right)\\\] where \\(Y\_{i, t}\\), \\(\eta\_{0\_{i, t}}\\),
\\(\eta\_{1\_{i, t}}\\), and \\(\boldsymbol{\varepsilon}\_{i, t}\\) are
random variables and \\(\theta\\) is a model parameter. \\(Y\_{i, t}\\)
is the observed random variable at time \\(t\\) and individual \\(i\\),
\\(\eta\_{0\_{i, t}}\\) (intercept) and \\(\eta\_{1\_{i, t}}\\) (slope)
form a vector of latent random variables at time \\(t\\) and individual
\\(i\\), and \\(\boldsymbol{\varepsilon}\_{i, t}\\) a vector of random
measurement errors at time \\(t\\) and individual \\(i\\). \\(\theta\\)
is the variance of \\(\boldsymbol{\varepsilon}\\).

In words, the observed score is the current level of the trajectory,
\\(\eta\_{0\_{i, t}}\\), plus measurement error. The loading of the
slope is 0, because the slope affects the score only through its effect
on the level over time.

The dynamic structure is given by \\\[\left( \begin{array}{c}
\eta\_{0\_{i, t}} \\\\ \eta\_{1\_{i, t}} \\\\ \end{array} \right) =
\left( \begin{array}{cc} 1 & 1 \\\\ 0 & 1 \\\\ \end{array} \right)
\left( \begin{array}{c} \eta\_{0\_{i, t - 1}} \\\\ \eta\_{1\_{i, t - 1}}
\\\\ \end{array} \right) .\\\]

In words, the level at each occasion is the previous level plus the
slope, \\(\eta\_{0\_{i, t}} = \eta\_{0\_{i, t - 1}} + \eta\_{1\_{i, t -
1}}\\), and the slope stays the same, \\(\eta\_{1\_{i, t}} =
\eta\_{1\_{i, t - 1}}\\). There is no process noise, so once an
individual’s intercept and slope are set, the trajectory is a straight
line. All the variation in the trajectories between individuals comes
from the initial condition below.

The mean vector and covariance matrix of the intercept and slope are
captured in the mean vector and covariance matrix of the initial
condition given by \\\[\boldsymbol{\mu}\_{0} = \left( \begin{array}{c}
\mu\_{\eta\_{0}} \\\\ \mu\_{\eta\_{1}} \\\\ \end{array} \right) \quad
\mathrm{and,}\\\]

\\\[\boldsymbol{\Sigma}\_{0} = \left( \begin{array}{cc}
\sigma^{2}\_{\eta\_{0}} & \sigma\_{\eta\_{0}, \eta\_{1}} \\\\
\sigma\_{\eta\_{1}, \eta\_{0}} & \sigma^{2}\_{\eta\_{1}} \\\\
\end{array} \right) .\\\]

Here, \\(\mu\_{\eta\_{0}}\\) and \\(\mu\_{\eta\_{1}}\\) are the average
intercept and slope, \\(\sigma^{2}\_{\eta\_{0}}\\) and
\\(\sigma^{2}\_{\eta\_{1}}\\) are how much individuals differ in them,
and \\(\sigma\_{\eta\_{0}, \eta\_{1}}\\) is their covariance, for
example, positive if individuals who start higher also grow faster.

### Diagrams

![Path diagram of the linear growth curve
model.](fig-vignettes-lin-growth-path.png)

Path diagram of the linear growth curve model.

The model as a path diagram, for the occasions \\(t = 0, 1, \ldots, T -
1\\). The intercept \\(\eta\_{0}\\) loads on every occasion with a
loading of 1, and the slope \\(\eta\_{1}\\) with loadings \\(0, 1, 2,
\ldots, T - 1\\), so the slope is the change per occasion.
\\(\varepsilon\_{t}\\) is the measurement error, with variance
\\(\theta\\). A loop on a variable is its variance and a curved
two-headed arrow a covariance. \\(\sigma^{2}\_{\eta\_{0}}\\) and
\\(\sigma^{2}\_{\eta\_{1}}\\) are the variances of the intercept and the
slope, and \\(\sigma\_{\eta\_{0}, \eta\_{1}}\\) their covariance.

![State space form of the linear growth curve
model.](fig-vignettes-lin-growth-ssm.png)

State space form of the linear growth curve model.

The same model in state space form, for the first three occasions (\\(t
= 0, 1, 2\\)), with every element of the state drawn. The state holds
the current level \\(\eta\_{0, t}\\) and the slope \\(\eta\_{1, t}\\).
The transition adds the slope to the level at each occasion and keeps
the slope, so \\(\eta\_{1, t} = \eta\_{1}\\) for all \\(t\\) and
\\(\eta\_{0, t} = \eta\_{0} + t \eta\_{1}\\); the level is measured by
\\(Y\_{t}\\) with a loading of 1 and measurement error
\\(\varepsilon\_{t}\\), with variance \\(\theta\\). There is no process
noise: the loadings of the path diagram are produced by the transition,
and the differences between individuals come from the initial condition.
At \\(t = 0\\), the initial condition \\(\boldsymbol{\eta}\_{0} \sim
\mathcal{N} \left( \boldsymbol{\mu}\_{0}, \boldsymbol{\Sigma}\_{0}
\right)\\) holds the means, variances, and covariance of the intercept
and the slope. Below the diagram, \\(\boldsymbol{\beta}\\),
\\(\boldsymbol{\Lambda}\\), and \\(\boldsymbol{\Theta}\\) are the same
matrices in compact form.
[`SimSSMLinGrowth()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowth.md)
fixes \\(\boldsymbol{\beta}\\) and \\(\boldsymbol{\Lambda}\\) and takes
`mu0` and the Cholesky factors of \\(\boldsymbol{\Sigma}\_{0}\\) and of
\\(\theta\\) (`sigma0_l` and `theta_l`).

## Data Generation

### Notation

Let \\(t = 5\\) be the number of time points and \\(n = 2000\\) be the
number of individuals.

Let the measurement error variance \\(\theta\\) be given by

\\\[\theta = 0.6 .\\\]

Let the initial condition \\(\boldsymbol{\eta}\_{0}\\) be given by

\\\[\boldsymbol{\eta}\_{0} \sim \mathcal{N} \left(
\boldsymbol{\mu}\_{0}, \boldsymbol{\Sigma}\_{0} \right)\\\]

\\\[\boldsymbol{\mu}\_{0} = \left( \begin{array}{c} 0.615 \\\\ 1.006
\\\\ \end{array} \right)\\\]

\\\[\boldsymbol{\Sigma}\_{0} = \left( \begin{array}{cc} 1.932 & 0.618
\\\\ 0.618 & 0.587 \\\\ \end{array} \right) .\\\]

### R Function Arguments

The simulation functions take covariance matrices as lower triangular
Cholesky factors (the arguments ending in `_l`).

\
`n`\
`#> [1] 2000`\
`time`\
`#> [1] 5`\
`mu0`\
`#> [1] 0.615 1.006`\
`sigma0`\
`#>       [,1]  [,2]`\
`#> [1,] 1.932 0.618`\
`#> [2,] 0.618 0.587`\
`sigma0_l`` ``# sigma0_l <- t(chol(sigma0))`\
`#>           [,1]      [,2]`\
`#> [1,] 1.3899640 0.0000000`\
`#> [2,] 0.4446158 0.6239525`\
`theta`\
`#> [1] 0.6`\
`theta_l`` ``# theta_l <- t(chol(theta))`\
`#> [1] 0.7745967`

### Visualizing the Dynamics Without Measurement Error (n = 5 with Different Initial Condition)

Without measurement error, each of the five simulated individuals
follows an exact straight line, with its own intercept and slope.

![](fig-vignettes-lin-growth-no-error-lin-growth-1.png)

### Using the `SimSSMLinGrowth` Function from the `simStateSpace` Package to Simulate Data

Now we add the measurement error and simulate the full data set.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the data in long format, one row per individual and time point, with the
columns `id`, `time`, and `y1`.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`simStateSpace`](https://github.com/jeksterslab/simStateSpace)`)`\
`sim`` ``<-`` `[`SimSSMLinGrowth`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowth.md)`(`\
`  n ``=`` ``n``,`\
`  time ``=`` ``time``,`\
`  mu0 ``=`` ``mu0``,`\
`  sigma0_l ``=`` ``sigma0_l``,`\
`  theta_l ``=`` ``theta_l``,`\
`  type ``=`` ``0`\
`)`\
`data`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data``)`\
`#>   id time       y1`\
`#> 1  1    0 1.723564`\
`#> 2  1    1 2.148702`\
`#> 3  1    2 3.384974`\
`#> 4  1    3 6.225426`\
`#> 5  1    4 5.993965`\
`#> 6  2    0 1.951861`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``data``)`\
`#>        id              time         y1         `\
`#>  Min.   :   1.0   Min.   :0   Min.   :-8.1963  `\
`#>  1st Qu.: 500.8   1st Qu.:1   1st Qu.: 0.4226  `\
`#>  Median :1000.5   Median :2   Median : 2.1950  `\
`#>  Mean   :1000.5   Mean   :2   Mean   : 2.5922  `\
`#>  3rd Qu.:1500.2   3rd Qu.:3   3rd Qu.: 4.4219  `\
`#>  Max.   :2000.0   Max.   :4   Max.   :17.0616`\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-lin-growth-error-lin-growth-1.png)

## Model Fitting

The same simulated data are fitted as a state space model with `dynr`
and `OpenMx` and as a structural equation model (SEM) with `lavaan` and
`OpenMx`, and as a multilevel model with `lme4`. All of them estimate
the parameters by maximum likelihood. The section Population Values and
Estimates compares their estimates with the population values.

### dynr

We fit the model with `dynr`, which estimates state space models by
maximum likelihood using the Kalman filter.

A `dynr` model is built from pieces that match the equations above: the
data, the initial condition, the measurement model, the dynamic model,
and the noise. Each parameter gets a name (for example, `sigma0_01` for
\\(\sigma\_{\eta\_{0}, \eta\_{1}}\\)\`), and the same name used in two
places means one shared parameter. The starting values are set to the
population values for convenience; with real data, they would be rough
guesses. The intercept and slope are indexed 0 and 1 in the parameter
names.

#### Prepare Data

Tell `dynr` which columns hold the individual, the time, and the
observed variables.

\
`dynr_data`` ``<-`` ``dynr``::`[`dynr.data`](https://rdrr.io/pkg/dynr/man/dynr.data.html)`(`\
`  dataframe ``=`` ``data``,`\
`  id ``=`` ``"id"``,`\
`  time ``=`` ``"time"``,`\
`  observed ``=`` ``"y1"`\
`)`

#### Prepare Initial Condition

The mean and covariance of the initial state are the means, variances,
and covariance of the intercept and slope, so these are the main
parameters of the growth model. The covariance matrix is symmetric, so
`sigma0_01` appears both above and below the diagonal.

\
`dynr_initial`` ``<-`` ``dynr``::`[`prep.initial`](https://rdrr.io/pkg/dynr/man/prep.initial.html)`(`\
`  values.inistate ``=`` ``mu0``,`\
`  params.inistate ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_0"``, ``"mu0_1"``)``,`\
`  values.inicov ``=`` ``sigma0``,`\
`  params.inicov ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"sigma0_00"``, ``"sigma0_01"``,`\
`      ``"sigma0_01"``, ``"sigma0_11"`\
`    ``)``,`\
`    nrow ``=`` ``2`\
`  ``)`\
`)`

#### Prepare Measurement Model

The loadings \\(\left( 1, 0 \right)\\) are fixed.

\
`dynr_measurement`` ``<-`` ``dynr``::`[`prep.measurement`](https://rdrr.io/pkg/dynr/man/prep.measurement.html)`(`\
`  values.load ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``0``)``, nrow ``=`` ``1``)``,`\
`  params.load ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``1``, ncol ``=`` ``2``)``,`\
`  state.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta_0"``, ``"eta_1"``)``,`\
`  obs.names ``=`` ``"y1"`\
`)`

#### Prepare Dynamic Process

The two formulas are the two rows of the dynamic equation. They have no
free parameters.

\
`dynr_dynamics`` ``<-`` ``dynr``::`[`prep.formulaDynamics`](https://rdrr.io/pkg/dynr/man/prep.formulaDynamics.html)`(`\
`  formula ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    ``eta_0`` ``~`` ``eta_0`` ``+`` ``eta_1``,`\
`    ``eta_1`` ``~`` ``eta_1`\
`  ``)``,`\
`  isContinuousTime ``=`` ``FALSE`\
`)`\
`#> Warning in dynr::prep.formulaDynamics(formula = list(eta_0 ~ eta_0 + eta_1, :`\
`#> You provided no start values: length(startval)==0. If you have no free`\
`#> parameters, keep calm and carry on.`

#### Prepare Process Noise

The process noise is fixed at zero, and the measurement error variance
\\(\theta\\) is estimated.

\
`dynr_noise`` ``<-`` ``dynr``::`[`prep.noise`](https://rdrr.io/pkg/dynr/man/prep.noise.html)`(`\
`  values.latent ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``0``, nrow ``=`` ``2``, ncol ``=`` ``2``)``,`\
`  params.latent ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``2``, ncol ``=`` ``2``)``,`\
`  values.observed ``=`` ``theta``,`\
`  params.observed ``=`` ``"theta"`\
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
`  outfile ``=`` ``"lin-growth.c"`\
`)`

The plot below shows the model equations with the parameter names.

![](fig-vignettes-lin-growth-model-lin-growth-1.png)

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
`#> Original fitted parameters:  -0.5214641 0.6200474 0.9859686 0.5990437 0.3061706 `\
`#> -0.901686 `\
`#> `\
`#> Transformed fitted parameters:  0.5936507 0.6200474 0.9859686 1.820377 `\
`#> 0.5573459 0.5765277 `\
`#> `\
`#> Doing end processing`\
`#> Successful trial`\
`#> Total Time: 5.789685 `\
`#> Backend Time: 5.362604`

#### Summary

The summary lists each estimate with its standard error, \\(t\\) value,
confidence interval, and \\(p\\) value.

\
[`summary`](https://rdrr.io/r/base/summary.html)`(``results``)`\
`#> Coefficients:`\
`#>           Estimate Std. Error t value ci.lower ci.upper            Pr(>|t|)    `\
`#> theta      0.59365    0.01084   54.78  0.57241  0.61489 <0.0000000000000002 ***`\
`#> mu0_0      0.62005    0.03299   18.79  0.55538  0.68471 <0.0000000000000002 ***`\
`#> mu0_1      0.98597    0.01783   55.30  0.95102  1.02091 <0.0000000000000002 ***`\
`#> sigma0_00  1.82038    0.06912   26.34  1.68491  1.95584 <0.0000000000000002 ***`\
`#> sigma0_01  0.55735    0.02815   19.80  0.50216  0.61253 <0.0000000000000002 ***`\
`#> sigma0_11  0.57653    0.02014   28.63  0.53706  0.61599 <0.0000000000000002 ***`\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> -2 log-likelihood value at convergence = 33425.11`\
`#> AIC = 33437.11`\
`#> BIC = 33480.37`

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

The matrices of the model, with the dynr labels. The transition matrix
carries the intercept and adds the slope, and there is no process noise
(`Q` is zero).

\
`mx_matrices`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`  ``MxMat``(`\
`    values ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``0``, ``1``, ``1``)``, nrow ``=`` ``2``)``,`\
`    labels ``=`` ``"fixed"``,`\
`    name ``=`` ``"A"`\
`  ``)``,`\
`  ``MxMat``(``values ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``0``)``, labels ``=`` ``"fixed"``, name ``=`` ``"B"``)``,`\
`  ``MxMat``(`\
`    values ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``0``)``, nrow ``=`` ``1``)``,`\
`    labels ``=`` ``"fixed"``,`\
`    name ``=`` ``"C"`\
`  ``)``,`\
`  ``MxMat``(``values ``=`` ``0``, labels ``=`` ``"fixed"``, name ``=`` ``"D"``)``,`\
`  ``MxMat``(`\
`    values ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``0``, nrow ``=`` ``2``, ncol ``=`` ``2``)``,`\
`    labels ``=`` ``"fixed"``,`\
`    name ``=`` ``"Q"`\
`  ``)``,`\
`  ``MxMat``(``values ``=`` ``theta``, labels ``=`` ``"theta"``, name ``=`` ``"R"``)``,`\
`  ``MxMat``(``values ``=`` ``mu0``, labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_0"``, ``"mu0_1"``)``, name ``=`` ``"mu0"``)``,`\
`  ``MxMat``(`\
`    values ``=`` ``sigma0``,`\
`    labels ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`      data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`        ``"sigma0_00"``, ``"sigma0_01"``,`\
`        ``"sigma0_01"``, ``"sigma0_11"`\
`      ``)``,`\
`      nrow ``=`` ``2`\
`    ``)``,`\
`    name ``=`` ``"sigma0"`\
`  ``)`\
`)`

``` r
fit_mx_ssm <- FitMxSSM(
  data = data,
  observed = "y1",
  matrices = mx_matrices,
  name = "lin_growth"
)
#> 
Beginning initial fit attempt
Fit attempt 0, fit=33425.1055806257, new current best! (was 33432.008312779)
                                                                            
#> 
#>  Solution found!  Final fit=33425.106 (started at 33432.008)  (1 attempt(s): 1 valid, 0 errors)
summary(fit_mx_ssm)
#> Summary of lin_growth 
#>  
#> free parameters:
#>        name matrix row col  Estimate  Std.Error A
#> 1     theta      R   1   1 0.5936771 0.01083943  
#> 2     mu0_0    mu0   1   1 0.6201697 0.03299043  
#> 3     mu0_1    mu0   2   1 0.9860179 0.01783248  
#> 4 sigma0_00 sigma0   1   1 1.8202777 0.06914596  
#> 5 sigma0_01 sigma0   2   1 0.5573863 0.02815680  
#> 6 sigma0_11 sigma0   2   2 0.5765853 0.02014077  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:              6                   9994              33425.11
#>    Saturated:             NA                     NA                    NA
#> Independence:             NA                     NA                    NA
#> Number of observations/statistics: 10000/10000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:       13437.11               33437.11                 33437.11
#> BIC:      -58623.04               33480.37                 33461.30
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-11-02 05:51:22 
#> Wall clock time: 17.71462 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

### lavaan: Structural Equation Model

The same model in SEM form, fitted with `lavaan`, should give nearly the
same estimates. In SEM, the five occasions are five observed variables
in one row per individual (wide format). The intercept factor has
loadings of 1 on all occasions, and the slope factor has loadings 0, 1,
2, 3, and 4.

#### Wide Data Set

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) with
`long = FALSE` returns the data in wide format.

\
`data_wide`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``, long ``=`` ``FALSE``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data_wide``)`\
`#>   id        y1_0       y1_1        y1_2        y1_3       y1_4`\
`#> 1  1  1.72356430  2.1487019  3.38497364  6.22542609  5.9939649`\
`#> 2  2  1.95186051  2.6140623  5.60107919  5.59791920  7.6881040`\
`#> 3  3  0.07909854 -1.1410633 -0.02940765  0.07103948 -0.4998955`\
`#> 4  4 -0.14827866 -0.6478676 -0.33440620 -0.29683127 -0.6950588`\
`#> 5  5  1.14625267  2.3844799  0.39101270 -0.76159166  1.6750324`\
`#> 6  6  1.02604168  2.1852479  2.55168739  2.13223414  4.1499727`

#### Model

The parameter labels match the dynr names. Giving the five error
variances the same label, `theta`, constrains them to be equal, as in
the state space model.

\
`model`` ``<-`` ``"`\
`  # factor loadings`\
`  eta0 =~ 1 * y1_0 + 1 * y1_1 + 1 * y1_2 + 1 * y1_3 + 1 * y1_4`\
`  eta1 =~ 0 * y1_0 + 1 * y1_1 + 2 * y1_2 + 3 * y1_3 + 4 * y1_4`\
`  # means of latent variables`\
`  eta0 ~ mu0_0 * 1`\
`  eta1 ~ mu0_1 * 1`\
`  # variances and covariances of latent variables`\
`  eta0 ~~ sigma0_00 * eta0`\
`  eta0 ~~ sigma0_01 * eta1`\
`  eta1 ~~ sigma0_11 * eta1`\
`  # constrain error variance theta to be equal`\
`  y1_0 ~~ theta * y1_0`\
`  y1_1 ~~ theta * y1_1`\
`  y1_2 ~~ theta * y1_2`\
`  y1_3 ~~ theta * y1_3`\
`  y1_4 ~~ theta * y1_4`\
`"`

\
`fit_lavaan`` ``<-`` ``lavaan``::`[`growth`](https://rdrr.io/pkg/lavaan/man/growth.html)`(`\
`  model ``=`` ``model``,`\
`  data ``=`` ``data_wide`\
`)`\
`lavaan``::``summary``(``fit_lavaan``)`\
`#> lavaan 0.7-3 ended normally after 30 iterations`\
`#> `\
`#>   Estimator                                         ML`\
`#>   Optimization method                           NLMINB`\
`#>   Number of model parameters                        10`\
`#>   Number of equality constraints                     4`\
`#> `\
`#>   Number of observations                          2000`\
`#> `\
`#> Model Test User Model:`\
`#>                                                       `\
`#>   Test statistic                                22.264`\
`#>   Degrees of freedom                                14`\
`#>   P-value (Chi-square)                           0.073`\
`#>                                                       `\
`#>   Browne's residual (NT model-based) test             `\
`#>   Test statistic                                22.030`\
`#>   Degrees of freedom                                14`\
`#>   P-value (Chi-square)                           0.078`\
`#> `\
`#> Parameter Estimates:`\
`#> `\
`#>   Standard errors                             Standard`\
`#>   Information                                 Expected`\
`#>   Information saturated (h1) model          Structured`\
`#> `\
`#> Latent Variables:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>   eta0 =~                                             `\
`#>     y1_0              1.000                           `\
`#>     y1_1              1.000                           `\
`#>     y1_2              1.000                           `\
`#>     y1_3              1.000                           `\
`#>     y1_4              1.000                           `\
`#>   eta1 =~                                             `\
`#>     y1_0              0.000                           `\
`#>     y1_1              1.000                           `\
`#>     y1_2              2.000                           `\
`#>     y1_3              3.000                           `\
`#>     y1_4              4.000                           `\
`#> `\
`#> Covariances:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>   eta0 ~~                                             `\
`#>     eta1    (s0_0)    0.557    0.028   19.794    0.000`\
`#> `\
`#> Intercepts:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>     eta0    (m0_0)    0.620    0.033   18.800    0.000`\
`#>     eta1    (m0_1)    0.986    0.018   55.295    0.000`\
`#> `\
`#> Variances:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>     eta0    (s0_0)    1.820    0.069   26.330    0.000`\
`#>     eta1    (s0_1)    0.577    0.020   28.629    0.000`\
`#>    .y1_0    (thet)    0.594    0.011   54.772    0.000`\
`#>    .y1_1    (thet)    0.594    0.011   54.772    0.000`\
`#>    .y1_2    (thet)    0.594    0.011   54.772    0.000`\
`#>    .y1_3    (thet)    0.594    0.011   54.772    0.000`\
`#>    .y1_4    (thet)    0.594    0.011   54.772    0.000`

### OpenMx: Structural Equation Model

The same SEM in `OpenMx`, written with
[`OpenMx::mxPath()`](https://rdrr.io/pkg/OpenMx/man/mxPath.html) on the
wide data. One-headed arrows are loadings and regressions, two-headed
arrows are variances and covariances, and arrows from `"one"` are means.
A path that is not written is fixed at zero, so the zero residual
variances, covariances, and intercepts of the lavaan model need no lines
here. The labels are the dynr names, and the starting values are the
population values.

\
`waves`` ``<-`` ``0``:``(``time`` ``-`` ``1``)`\
`manifest`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y1_"``, ``waves``)`\
`mx_sem`` ``<-`` ``OpenMx``::`[`mxModel`](https://rdrr.io/pkg/OpenMx/man/mxModel.html)`(`\
`  name ``=`` ``"lin_growth_sem"``,`\
`  type ``=`` ``"RAM"``,`\
`  manifestVars ``=`` ``manifest``,`\
`  latentVars ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta0"``, ``"eta1"``)``,`\
`  ``OpenMx``::`[`mxData`](https://rdrr.io/pkg/OpenMx/man/mxData.html)`(``observed ``=`` ``data_wide``[``, ``manifest``]``, type ``=`` ``"raw"``)``,`\
`  ``# factor loadings`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``"eta0"``, to ``=`` ``manifest``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``"eta1"``, to ``=`` ``manifest``, free ``=`` ``FALSE``, values ``=`` ``waves``)``,`\
`  ``# means of latent variables`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(`\
`    from ``=`` ``"one"``,`\
`    to ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta0"``, ``"eta1"``)``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_0"``, ``"mu0_1"``)``,`\
`    values ``=`` ``mu0`\
`  ``)``,`\
`  ``# variances and covariances of latent variables`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(`\
`    from ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta0"``, ``"eta0"``, ``"eta1"``)``,`\
`    to ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta0"``, ``"eta1"``, ``"eta1"``)``,`\
`    arrows ``=`` ``2``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"sigma0_00"``, ``"sigma0_01"``, ``"sigma0_11"``)``,`\
`    values ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``sigma0``[``1``, ``1``]``, ``sigma0``[``1``, ``2``]``, ``sigma0``[``2``, ``2``]``)`\
`  ``)``,`\
`  ``# error variance theta, equal over time`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``manifest``, arrows ``=`` ``2``, labels ``=`` ``"theta"``, values ``=`` ``theta``)`\
`)`

\
`fit_mx_sem`` ``<-`` ``OpenMx``::`[`mxRun`](https://rdrr.io/pkg/OpenMx/man/mxRun.html)`(`\
`  model ``=`` ``mx_sem``,`\
`  silent ``=`` ``TRUE`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_mx_sem``)`\
`#> Summary of lin_growth_sem `\
`#>  `\
`#> free parameters:`\
`#>        name matrix  row  col  Estimate  Std.Error A`\
`#> 1     theta      S y1_0 y1_0 0.5936771 0.01083906  `\
`#> 2 sigma0_00      S eta0 eta0 1.8202777 0.06913399  `\
`#> 3 sigma0_01      S eta0 eta1 0.5573863 0.02815924  `\
`#> 4 sigma0_11      S eta1 eta1 0.5765853 0.02013960  `\
`#> 5     mu0_0      M    1 eta0 0.6201697 0.03298763  `\
`#> 6     mu0_1      M    1 eta1 0.9860179 0.01783185  `\
`#> `\
`#> Model Statistics: `\
`#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)`\
`#>        Model:              6                   9994              33425.11`\
`#>    Saturated:             20                   9980                    NA`\
`#> Independence:             10                   9990                    NA`\
`#> Number of observations/statistics: 2000/10000`\
`#> `\
`#> Information Criteria: `\
`#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted`\
`#> AIC:       13437.11               33437.11                 33437.15`\
`#> BIC:      -42538.31               33470.71                 33451.65`\
`#> To get additional fit indices, see help(mxRefModels)`\
`#> timestamp: 2026-11-02 05:52:02 `\
`#> Wall clock time: 0.01413393 secs `\
`#> optimizer:  SLSQP `\
`#> OpenMx version number: 2.22.11 `\
`#> Need help?  See help(mxSummary)`

### lme4: Multilevel Model

The linear growth curve model is also a multilevel (mixed-effects) model
of the long data: the occasions are nested in individuals, the fixed
intercept and slope are the means \\(\mu\_{\eta\_{0}}\\) and
\\(\mu\_{\eta\_{1}}\\) (`mu0_0` and `mu0_1`), the random intercept and
slope have the covariance matrix \\(\boldsymbol{\Sigma}\_{0}\\)
(`sigma0_00`, `sigma0_01`, and `sigma0_11`), and the residual variance
is \\(\theta\\).
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) fits it with a
random intercept and a random slope of `time` for each `id` (Bates et
al., 2015). `REML = FALSE` gives maximum likelihood estimates, as in the
other methods.

\
`fit_lme4`` ``<-`` ``lme4``::`[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)`(`\
`  formula ``=`` ``y1`` ``~`` ``time`` ``+`` ``(``time`` ``|`` ``id``)``,`\
`  data ``=`` ``data``,`\
`  REML ``=`` ``FALSE`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_lme4``)`\
`#> Linear mixed model fit by maximum likelihood  ['lmerMod']`\
`#> Formula: y1 ~ time + (time | id)`\
`#>    Data: data`\
`#> `\
`#>       AIC       BIC    logLik -2*log(L)  df.resid `\
`#>   33437.1   33480.4  -16712.6   33425.1      9994 `\
`#> `\
`#> Scaled residuals: `\
`#>      Min       1Q   Median       3Q      Max `\
`#> -2.93296 -0.55042  0.01077  0.54126  3.06952 `\
`#> `\
`#> Random effects:`\
`#>  Groups   Name        Variance Std.Dev. Corr `\
`#>  id       (Intercept) 1.8203   1.3492        `\
`#>           time        0.5766   0.7593   0.54 `\
`#>  Residual             0.5937   0.7705        `\
`#> Number of obs: 10000, groups:  id, 2000`\
`#> `\
`#> Fixed effects:`\
`#>             Estimate Std. Error t value`\
`#> (Intercept)  0.62017    0.03299    18.8`\
`#> time         0.98602    0.01783    55.3`\
`#> `\
`#> Correlation of Fixed Effects:`\
`#>      (Intr)`\
`#> time 0.373`

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
`population`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  mu0_0 ``=`` ``mu0``[``1``]``,`\
`  mu0_1 ``=`` ``mu0``[``2``]``,`\
`  sigma0_00 ``=`` ``sigma0``[``1``, ``1``]``,`\
`  sigma0_01 ``=`` ``sigma0``[``1``, ``2``]``,`\
`  sigma0_11 ``=`` ``sigma0``[``2``, ``2``]``,`\
`  theta ``=`` ``theta`\
`)`

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

| Parameter | Population | Estimate |    SE |   Bias |  2.5% | 97.5% | Covers |
|:----------|-----------:|---------:|------:|-------:|------:|------:|:-------|
| mu0_0     |      0.615 |    0.620 | 0.033 |  0.005 | 0.555 | 0.685 | TRUE   |
| mu0_1     |      1.006 |    0.986 | 0.018 | -0.020 | 0.951 | 1.021 | TRUE   |
| sigma0_00 |      1.932 |    1.820 | 0.069 | -0.112 | 1.685 | 1.956 | TRUE   |
| sigma0_01 |      0.618 |    0.557 | 0.028 | -0.061 | 0.502 | 0.613 | FALSE  |
| sigma0_11 |      0.587 |    0.577 | 0.020 | -0.010 | 0.537 | 0.616 | TRUE   |
| theta     |      0.600 |    0.594 | 0.011 | -0.006 | 0.572 | 0.615 | TRUE   |

Population values and dynr estimates (one simulated data set) {.table}

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

| Parameter | Population | Estimate |    SE |   Bias |  2.5% | 97.5% | Covers |
|:----------|-----------:|---------:|------:|-------:|------:|------:|:-------|
| mu0_0     |      0.615 |    0.620 | 0.033 |  0.005 | 0.556 | 0.685 | TRUE   |
| mu0_1     |      1.006 |    0.986 | 0.018 | -0.020 | 0.951 | 1.021 | TRUE   |
| sigma0_00 |      1.932 |    1.820 | 0.069 | -0.112 | 1.685 | 1.956 | TRUE   |
| sigma0_01 |      0.618 |    0.557 | 0.028 | -0.061 | 0.502 | 0.613 | FALSE  |
| sigma0_11 |      0.587 |    0.577 | 0.020 | -0.010 | 0.537 | 0.616 | TRUE   |
| theta     |      0.600 |    0.594 | 0.011 | -0.006 | 0.572 | 0.615 | TRUE   |

Population values and OpenMx state space estimates (one simulated data
set) {.table}

### lavaan: Structural Equation Model

\
`pe`` ``<-`` ``lavaan``::`[`parameterEstimates`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html)`(``fit_lavaan``)`\
`# a label shared by several paths is one parameter`\
`pe`` ``<-`` ``pe``[``pe``$``label`` ``!=`` ``""`` ``&`` ``!`[`duplicated`](https://rdrr.io/r/base/duplicated.html)`(``pe``$``label``)``, ``]`\
`lavaan_est`` ``<-`` ``stats``::`[`setNames`](https://rdrr.io/r/stats/setNames.html)`(``object ``=`` ``pe``$``est``, nm ``=`` ``pe``$``label``)`\
`lavaan_se`` ``<-`` ``stats``::`[`setNames`](https://rdrr.io/r/stats/setNames.html)`(``object ``=`` ``pe``$``se``, nm ``=`` ``pe``$``label``)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` ``EstimatesTable``(`\
`    population ``=`` ``population``,`\
`    estimate ``=`` ``lavaan_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    se ``=`` ``lavaan_se``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and lavaan estimates (one simulated data set)"`\
`)`

| Parameter | Population | Estimate |    SE |   Bias |  2.5% | 97.5% | Covers |
|:----------|-----------:|---------:|------:|-------:|------:|------:|:-------|
| mu0_0     |      0.615 |    0.620 | 0.033 |  0.005 | 0.556 | 0.685 | TRUE   |
| mu0_1     |      1.006 |    0.986 | 0.018 | -0.020 | 0.951 | 1.021 | TRUE   |
| sigma0_00 |      1.932 |    1.820 | 0.069 | -0.112 | 1.685 | 1.956 | TRUE   |
| sigma0_01 |      0.618 |    0.557 | 0.028 | -0.061 | 0.502 | 0.613 | FALSE  |
| sigma0_11 |      0.587 |    0.577 | 0.020 | -0.010 | 0.537 | 0.616 | TRUE   |
| theta     |      0.600 |    0.594 | 0.011 | -0.006 | 0.572 | 0.615 | TRUE   |

Population values and lavaan estimates (one simulated data set) {.table}

### OpenMx: Structural Equation Model

\
`mx_sem_est`` ``<-`` ``OpenMx``::`[`omxGetParameters`](https://rdrr.io/pkg/OpenMx/man/omxGetParameters.html)`(``fit_mx_sem``)`\
`mx_sem_se`` ``<-`` ``fit_mx_sem``$``output``$``standardErrors``[``, ``1``]`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` ``EstimatesTable``(`\
`    population ``=`` ``population``,`\
`    estimate ``=`` ``mx_sem_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    se ``=`` ``mx_sem_se``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and OpenMx SEM estimates (one simulated data set)"`\
`)`

| Parameter | Population | Estimate |    SE |   Bias |  2.5% | 97.5% | Covers |
|:----------|-----------:|---------:|------:|-------:|------:|------:|:-------|
| mu0_0     |      0.615 |    0.620 | 0.033 |  0.005 | 0.556 | 0.685 | TRUE   |
| mu0_1     |      1.006 |    0.986 | 0.018 | -0.020 | 0.951 | 1.021 | TRUE   |
| sigma0_00 |      1.932 |    1.820 | 0.069 | -0.112 | 1.685 | 1.956 | TRUE   |
| sigma0_01 |      0.618 |    0.557 | 0.028 | -0.061 | 0.502 | 0.613 | FALSE  |
| sigma0_11 |      0.587 |    0.577 | 0.020 | -0.010 | 0.537 | 0.616 | TRUE   |
| theta     |      0.600 |    0.594 | 0.011 | -0.006 | 0.572 | 0.615 | TRUE   |

Population values and OpenMx SEM estimates (one simulated data set)
{.table}

### lme4: Multilevel Model

lme4 gives standard errors for the fixed effects only. For the
variances, the interval is the profile likelihood interval of the
standard deviation, squared; lme4 has no interval for a covariance.

\
`vc`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``lme4``::`[`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html)`(``fit_lme4``)``)`\
`random`` ``<-`` ``vc``[``vc``$``grp`` ``==`` ``"id"``, ``]`\
`lme4_est`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  mu0_0 ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``lme4``::`[`fixef`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)`(``fit_lme4``)``[``"(Intercept)"``]``)``,`\
`  mu0_1 ``=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``lme4``::`[`fixef`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)`(``fit_lme4``)``[``"time"``]``)``,`\
`  sigma0_00 ``=`` ``random``$``vcov``[``random``$``var1`` ``==`` ``"(Intercept)"`` ``&`` `[`is.na`](https://rdrr.io/r/base/NA.html)`(``random``$``var2``)``]``,`\
`  sigma0_01 ``=`` ``random``$``vcov``[``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``random``$``var2``)``]``,`\
`  sigma0_11 ``=`` ``random``$``vcov``[``random``$``var1`` ``==`` ``"time"`` ``&`` `[`is.na`](https://rdrr.io/r/base/NA.html)`(``random``$``var2``)``]``,`\
`  theta ``=`` ``vc``$``vcov``[``vc``$``grp`` ``==`` ``"Residual"``]`\
`)`\
`lme4_se`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(`[`diag`](https://rdrr.io/r/base/diag.html)`(`[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``stats``::`[`vcov`](https://rdrr.io/r/stats/vcov.html)`(``fit_lme4``)``)``)``)``,`\
`  `[`rep`](https://rdrr.io/r/base/rep.html)`(``x ``=`` ``NA``, times ``=`` ``4``)`\
`)`\
[`names`](https://rdrr.io/r/base/names.html)`(``lme4_se``)`` ``<-`` `[`names`](https://rdrr.io/r/base/names.html)`(``lme4_est``)`\
`profile_ci`` ``<-`` ``stats``::`[`confint`](https://rdrr.io/r/stats/confint.html)`(`\
`  ``fit_lme4``,`\
`  method ``=`` ``"profile"``,`\
`  oldNames ``=`` ``FALSE`\
`)`\
`#> Warning: 'oldNames' is deprecated. Please use 'signames' instead.`\
`#> Computing profile confidence intervals ...`\
`z`` ``<-`` ``stats``::`[`qnorm`](https://rdrr.io/r/stats/Normal.html)`(``p ``=`` ``0.975``)`\
`lme4_ll`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``lme4_est``[``1``:``2``]`` ``-`` ``z`` ``*`` ``lme4_se``[``1``:``2``]``,`\
`  sigma0_00 ``=`` ``profile_ci``[``"sd_(Intercept)|id"``, ``1``]``^``2``,`\
`  sigma0_01 ``=`` ``NA``,`\
`  sigma0_11 ``=`` ``profile_ci``[``"sd_time|id"``, ``1``]``^``2``,`\
`  theta ``=`` ``profile_ci``[``"sigma"``, ``1``]``^``2`\
`)`\
`lme4_ul`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``lme4_est``[``1``:``2``]`` ``+`` ``z`` ``*`` ``lme4_se``[``1``:``2``]``,`\
`  sigma0_00 ``=`` ``profile_ci``[``"sd_(Intercept)|id"``, ``2``]``^``2``,`\
`  sigma0_01 ``=`` ``NA``,`\
`  sigma0_11 ``=`` ``profile_ci``[``"sd_time|id"``, ``2``]``^``2``,`\
`  theta ``=`` ``profile_ci``[``"sigma"``, ``2``]``^``2`\
`)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`  x ``=`` ``EstimatesTable``(`\
`    population ``=`` ``population``,`\
`    estimate ``=`` ``lme4_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    se ``=`` ``lme4_se``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    ll ``=`` ``lme4_ll``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``,`\
`    ul ``=`` ``lme4_ul``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` `[`paste`](https://rdrr.io/r/base/paste.html)`(`\
`    ``"Population values and lme4 estimates (one simulated data set);"``,`\
`    ``"intervals for the variances: profile likelihood"`\
`  ``)`\
`)`

| Parameter | Population | Estimate |    SE |   Bias |  2.5% | 97.5% | Covers |
|:----------|-----------:|---------:|------:|-------:|------:|------:|:-------|
| mu0_0     |      0.615 |    0.620 | 0.033 |  0.005 | 0.556 | 0.685 | TRUE   |
| mu0_1     |      1.006 |    0.986 | 0.018 | -0.020 | 0.951 | 1.021 | TRUE   |
| sigma0_00 |      1.932 |    1.820 |    NA | -0.112 | 1.690 | 1.961 | TRUE   |
| sigma0_01 |      0.618 |    0.557 |    NA | -0.061 |    NA |    NA | NA     |
| sigma0_11 |      0.587 |    0.577 |    NA | -0.010 | 0.539 | 0.618 | TRUE   |
| theta     |      0.600 |    0.594 |    NA | -0.006 | 0.573 | 0.615 | TRUE   |

Population values and lme4 estimates (one simulated data set); intervals
for the variances: profile likelihood {.table}

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
``     `lavaan`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``lavaan_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``)``,`\
``     `OpenMx (SEM)`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``mx_sem_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``)``,`\
``     `lme4`  ```=`` `[`unname`](https://rdrr.io/r/base/unname.html)`(``lme4_est``[`[`names`](https://rdrr.io/r/base/names.html)`(``population``)``]``)``,`\
`    check.names ``=`` ``FALSE`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and the estimates of every method"`\
`)`

| Parameter | Population |  dynr | OpenMx (SSM) | lavaan | OpenMx (SEM) |  lme4 |
|:----------|-----------:|------:|-------------:|-------:|-------------:|------:|
| mu0_0     |      0.615 | 0.620 |        0.620 |  0.620 |        0.620 | 0.620 |
| mu0_1     |      1.006 | 0.986 |        0.986 |  0.986 |        0.986 | 0.986 |
| sigma0_00 |      1.932 | 1.820 |        1.820 |  1.820 |        1.820 | 1.820 |
| sigma0_01 |      0.618 | 0.557 |        0.557 |  0.557 |        0.557 | 0.557 |
| sigma0_11 |      0.587 | 0.577 |        0.577 |  0.577 |        0.577 | 0.577 |
| theta     |      0.600 | 0.594 |        0.594 |  0.594 |        0.594 | 0.594 |

Population values and the estimates of every method {.table
style="width:99%;"}

## References

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
*67*(1). <https://doi.org/10.18637/jss.v067.i01>

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

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software*, *48*(2).
<https://doi.org/10.18637/jss.v048.i02>
