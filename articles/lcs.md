# The Bivariate Dual Change Score Model

## Overview

Latent change score models describe change itself (McArdle, 2009). Take
vocabulary and reading comprehension measured once a year in children.
From one occasion to the next, each child’s vocabulary changes. In the
dual change score model, the change has two sources: a constant amount
each child gains every year (the slope), and an amount that depends on
where the child stood the year before (the self-feedback). In the
bivariate model, the change in one variable can also depend on the
previous level of the other variable (the coupling), for example, if
children with a larger vocabulary gain more in reading comprehension.

The dual change score model is usually fitted as a structural equation
model (SEM). It can also be written as a state space model: the true
levels of the two variables are latent states, and the slopes are latent
states that never change.

This tutorial does three things:

1.  It writes the bivariate dual change score model as a state space
    model.
2.  It simulates data from the model with
    [`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md).
3.  It fits the model to the simulated data as a state space model with
    `dynr` and `OpenMx` and as a structural equation model (SEM) with
    `lavaan` and `OpenMx`, and compares every set of estimates with the
    values used to simulate the data.

No prior knowledge of state space models is needed.

## Model

The measurement model is given by \\\[\left( \begin{array}{c} Y\_{1\_{i,
t}} \\\\ Y\_{2\_{i, t}} \\\\ \end{array} \right) = \left(
\begin{array}{cccc} 1 & 0 & 0 & 0 \\\\ 0 & 1 & 0 & 0 \\\\ \end{array}
\right) \left( \begin{array}{c} \eta\_{1\_{i, t}} \\\\ \eta\_{2\_{i, t}}
\\\\ \eta\_{3\_{i, t}} \\\\ \eta\_{4\_{i, t}} \\\\ \end{array} \right) +
\left( \begin{array}{c} \varepsilon\_{1\_{i, t}} \\\\
\varepsilon\_{2\_{i, t}} \\\\ \end{array} \right), \quad \mathrm{with}
\quad \left( \begin{array}{c} \varepsilon\_{1\_{i, t}} \\\\
\varepsilon\_{2\_{i, t}} \\\\ \end{array} \right) \sim \mathcal{N}
\left( \mathbf{0}, \left( \begin{array}{cc} \theta\_{11} & 0 \\\\ 0 &
\theta\_{22} \\\\ \end{array} \right) \right)\\\] where \\(Y\_{1\_{i,
t}}\\) and \\(Y\_{2\_{i, t}}\\) are the two observed variables at time
\\(t\\) and individual \\(i\\), \\(\eta\_{1\_{i, t}}\\) and
\\(\eta\_{2\_{i, t}}\\) are their true levels, \\(\eta\_{3\_{i, t}}\\)
and \\(\eta\_{4\_{i, t}}\\) are their slopes, and \\(\theta\_{11}\\) and
\\(\theta\_{22}\\) are the measurement error variances.

In words, each observed score is the true level plus measurement error.
The slopes have loadings of 0: they affect the scores only through the
change in the true levels.

The change in the true level of the first variable from occasion \\(t -
1\\) to occasion \\(t\\) is \\\[\Delta \eta\_{1\_{i, t}} = \eta\_{3\_{i,
t - 1}} + b\_{1} \eta\_{1\_{i, t - 1}} + g\_{12} \eta\_{2\_{i, t - 1}}
,\\\] and similarly for the second variable, \\(\Delta \eta\_{2\_{i, t}}
= \eta\_{4\_{i, t - 1}} + b\_{2} \eta\_{2\_{i, t - 1}} + g\_{21}
\eta\_{1\_{i, t - 1}}\\). The true level at occasion \\(t\\) is the
previous true level plus the change, \\(\eta\_{1\_{i, t}} =
\eta\_{1\_{i, t - 1}} + \Delta \eta\_{1\_{i, t}}\\).

Writing the change into the level gives the dynamic structure \\\[\left(
\begin{array}{c} \eta\_{1\_{i, t}} \\\\ \eta\_{2\_{i, t}} \\\\
\eta\_{3\_{i, t}} \\\\ \eta\_{4\_{i, t}} \\\\ \end{array} \right) =
\left( \begin{array}{cccc} 1 + b\_{1} & g\_{12} & 1 & 0 \\\\ g\_{21} &
1 + b\_{2} & 0 & 1 \\\\ 0 & 0 & 1 & 0 \\\\ 0 & 0 & 0 & 1 \\\\
\end{array} \right) \left( \begin{array}{c} \eta\_{1\_{i, t - 1}} \\\\
\eta\_{2\_{i, t - 1}} \\\\ \eta\_{3\_{i, t - 1}} \\\\ \eta\_{4\_{i, t -
1}} \\\\ \end{array} \right) .\\\]

In words, \\(b\_{1}\\) and \\(b\_{2}\\) are the self-feedback effects: a
negative value means that individuals who stand higher gain less, so the
levels level off over time. \\(g\_{12}\\) is the coupling from the
second variable to the first: a positive value means that a higher level
of the second variable is followed by a larger gain in the first.
\\(g\_{21}\\) is the reverse. The slopes stay the same over time. As in
most applications of the model, there is no process noise: once the
initial levels and the slopes are set, the true trajectories are fixed.

The initial condition holds the means, variances, and covariances of the
initial true levels and the slopes, \\\[\boldsymbol{\eta}\_{i, 0} \sim
\mathcal{N} \left( \boldsymbol{\mu}\_{0}, \boldsymbol{\Sigma}\_{0}
\right) .\\\]

For example, a positive covariance between the initial level and the
slope of the same variable means that individuals who start higher also
tend to gain more each year.

### Diagrams

![Path diagram of the dual change score
model.](fig-vignettes-lcs-path.png)

Path diagram of the dual change score model.

The bivariate dual change score model as a path diagram, for the first
three occasions (\\(t = 0, 1, 2\\)). The true level at each occasion,
\\(\eta\_{1}\\) or \\(\eta\_{2}\\), is the previous level plus the
latent change \\(\Delta \eta\_{1}\\) or \\(\Delta \eta\_{2}\\) (the
arrows with a loading of 1). Each change is the slope (\\(\eta\_{3}\\)
or \\(\eta\_{4}\\)), plus the self-feedback effect (\\(b\_{1}\\) or
\\(b\_{2}\\)) of the previous level of the same variable, plus the
coupling effect (\\(g\_{12}\\) or \\(g\_{21}\\)) of the previous level
of the other variable. \\(\varepsilon\_{1}\\) and \\(\varepsilon\_{2}\\)
are the measurement errors, with variances \\(\theta\_{11}\\) and
\\(\theta\_{22}\\). A loop on a variable is its variance and a curved
two-headed arrow a covariance. The initial levels and the slopes have
variances \\(\sigma\_{11}\\) to \\(\sigma\_{44}\\) and all covary; to
keep the diagram readable, only the covariances \\(\sigma\_{12}\\),
\\(\sigma\_{13}\\), \\(\sigma\_{24}\\), and \\(\sigma\_{34}\\) are
drawn.

![State space form of the dual change score
model.](fig-vignettes-lcs-ssm.png)

State space form of the dual change score model.

The same model in state space form, for the first three occasions, with
every element of the state drawn. The latent changes are not states:
writing each change into the level gives the arrows \\(1 + b\_{1}\\),
\\(1 + b\_{2}\\), \\(g\_{12}\\), and \\(g\_{21}\\) between the levels
\\(\eta\_{1, t}\\) and \\(\eta\_{2, t}\\), and the arrows with a loading
of 1 from the slopes. The slopes are states with a transition of 1, so
\\(\eta\_{3, t} = \eta\_{3}\\) and \\(\eta\_{4, t} = \eta\_{4}\\) for
all \\(t\\). There is no process noise. At \\(t = 0\\), the initial
condition \\(\boldsymbol{\eta}\_{0} \sim \mathcal{N} \left(
\boldsymbol{\mu}\_{0}, \boldsymbol{\Sigma}\_{0} \right)\\) gives all the
variances and covariances of the levels and the slopes at \\(t = 0\\).
Below the diagram, all the latent variables at one occasion form the
state vector \\(\boldsymbol{\eta}\_{t}\\), and all the observed
variables form \\(\mathbf{y}\_{t}\\): \\(\boldsymbol{\beta}\\) carries
the state from one occasion to the next, \\(\boldsymbol{\Lambda}\\) maps
the state to the observed variables, and \\(\boldsymbol{\Theta}\\) is
the covariance matrix of the measurement error
\\(\boldsymbol{\varepsilon}\_{t}\\).
[`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md)
takes `beta`, `lambda`, `mu0`, and the Cholesky factors of
\\(\boldsymbol{\Theta}\\) and \\(\boldsymbol{\Sigma}\_{0}\\) (`theta_l`
and `sigma0_l`); `psi_l` is a matrix of zeros.

## Data Generation

### Notation

Let \\(t = 6\\) be the number of time points and \\(n = 2000\\) be the
number of individuals.

Let the self-feedback effects be \\(b\_{1} = -0.2\\) and \\(b\_{2} =
-0.3\\), and the coupling effects be \\(g\_{12} = 0.1\\) and \\(g\_{21}
= 0.2\\), so that

\\\[\boldsymbol{\beta} = \left( \begin{array}{cccc} 0.8 & 0.1 & 1 & 0
\\\\ 0.2 & 0.7 & 0 & 1 \\\\ 0 & 0 & 1 & 0 \\\\ 0 & 0 & 0 & 1 \\\\
\end{array} \right) .\\\]

Let the initial condition be

\\\[\boldsymbol{\mu}\_{0} = \left( \begin{array}{c} 2 \\\\ 1 \\\\ 0.5
\\\\ 0.4 \\\\ \end{array} \right) \quad \mathrm{and,} \quad
\boldsymbol{\Sigma}\_{0} = \left( \begin{array}{cccc} 1 & 0.3 & 0.15 &
0.05 \\\\ 0.3 & 1 & 0.05 & 0.15 \\\\ 0.15 & 0.05 & 0.25 & 0.05 \\\\ 0.05
& 0.15 & 0.05 & 0.25 \\\\ \end{array} \right) ,\\\]

and the measurement error variances be \\(\theta\_{11} = 0.2\\) and
\\(\theta\_{22} = 0.2\\).

### R Function Arguments

The simulation functions take covariance matrices as lower triangular
Cholesky factors (the arguments ending in `_l`). There is no process
noise, so `psi` and `psi_l` are zero matrices.

\
`n`\
`#> [1] 2000`\
`time`\
`#> [1] 6`\
`mu0`\
`#> [1] 2.0 1.0 0.5 0.4`\
`sigma0`\
`#>      [,1] [,2] [,3] [,4]`\
`#> [1,] 1.00 0.30 0.15 0.05`\
`#> [2,] 0.30 1.00 0.05 0.15`\
`#> [3,] 0.15 0.05 0.25 0.05`\
`#> [4,] 0.05 0.15 0.05 0.25`\
`sigma0_l`` ``# sigma0_l <- t(chol(sigma0))`\
`#>      [,1]        [,2]       [,3]      [,4]`\
`#> [1,] 1.00 0.000000000 0.00000000 0.0000000`\
`#> [2,] 0.30 0.953939201 0.00000000 0.0000000`\
`#> [3,] 0.15 0.005241424 0.47694080 0.0000000`\
`#> [4,] 0.05 0.141518453 0.08755435 0.4688355`\
`alpha`\
`#> [1] 0 0 0 0`\
`beta`\
`#>      [,1] [,2] [,3] [,4]`\
`#> [1,]  0.8  0.1    1    0`\
`#> [2,]  0.2  0.7    0    1`\
`#> [3,]  0.0  0.0    1    0`\
`#> [4,]  0.0  0.0    0    1`\
`psi`\
`#>      [,1] [,2] [,3] [,4]`\
`#> [1,]    0    0    0    0`\
`#> [2,]    0    0    0    0`\
`#> [3,]    0    0    0    0`\
`#> [4,]    0    0    0    0`\
`psi_l`\
`#>      [,1] [,2] [,3] [,4]`\
`#> [1,]    0    0    0    0`\
`#> [2,]    0    0    0    0`\
`#> [3,]    0    0    0    0`\
`#> [4,]    0    0    0    0`\
`nu`\
`#> [1] 0 0`\
`lambda`\
`#>      [,1] [,2] [,3] [,4]`\
`#> [1,]    1    0    0    0`\
`#> [2,]    0    1    0    0`\
`theta`\
`#>      [,1] [,2]`\
`#> [1,]  0.2  0.0`\
`#> [2,]  0.0  0.2`\
`theta_l`` ``# theta_l <- t(chol(theta))`\
`#>           [,1]      [,2]`\
`#> [1,] 0.4472136 0.0000000`\
`#> [2,] 0.0000000 0.4472136`

### Visualizing the Dynamics Without Measurement Error (n = 5 with Different Initial Condition)

Without measurement error, each individual’s true levels follow smooth
curves that bend as the self-feedback slows the gains.

![](fig-vignettes-lcs-no-error-lcs-1.png)![](fig-vignettes-lcs-no-error-lcs-2.png)

### Using the `SimSSMFixed` Function from the `simStateSpace` Package to Simulate Data

Now we add the measurement error and simulate the full data set.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the data in long format, one row per individual and time point, with the
columns `id`, `time`, `y1`, and `y2`.

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
`#>   id time       y1           y2`\
`#> 1  1    0 1.239124 -0.034741096`\
`#> 2  1    1 2.372287 -0.050275519`\
`#> 3  1    2 2.550748 -0.363167146`\
`#> 4  1    3 3.417814 -0.164406246`\
`#> 5  1    4 3.685087  0.007406512`\
`#> 6  1    5 3.473986 -0.348501993`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``data``)`\
`#>        id              time           y1               y2         `\
`#>  Min.   :   1.0   Min.   :0.0   Min.   :-4.200   Min.   :-5.4249  `\
`#>  1st Qu.: 500.8   1st Qu.:1.0   1st Qu.: 1.429   1st Qu.: 0.8064  `\
`#>  Median :1000.5   Median :2.5   Median : 2.446   Median : 1.8475  `\
`#>  Mean   :1000.5   Mean   :2.5   Mean   : 2.515   Mean   : 1.9371  `\
`#>  3rd Qu.:1500.2   3rd Qu.:4.0   3rd Qu.: 3.516   3rd Qu.: 2.9450  `\
`#>  Max.   :2000.0   Max.   :5.0   Max.   : 9.998   Max.   : 8.7711`\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-lcs-error-lcs-1.png)![](fig-vignettes-lcs-error-lcs-2.png)

## Model Fitting

The same simulated data are fitted as a state space model with `dynr`
and `OpenMx` and as a structural equation model (SEM) with `lavaan` and
`OpenMx`. All of them estimate the parameters by maximum likelihood. The
section Population Values and Estimates compares their estimates with
the population values.

### dynr

We fit the model with `dynr`, which estimates state space models by
maximum likelihood using the Kalman filter.

A `dynr` model is built from pieces that match the equations above: the
data, the initial condition, the measurement model, the dynamic model,
and the noise. Each parameter gets a name (for example, `g_1_2` for
\\(g\_{12}\\)), and the same name used in two places means one shared
parameter. The starting values are set to the population values for
convenience; with real data, they would be rough guesses.

#### Prepare Data

Tell `dynr` which columns hold the individual, the time, and the
observed variables.

\
`dynr_data`` ``<-`` ``dynr``::`[`dynr.data`](https://rdrr.io/pkg/dynr/man/dynr.data.html)`(`\
`  dataframe ``=`` ``data``,`\
`  id ``=`` ``"id"``,`\
`  time ``=`` ``"time"``,`\
`  observed ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``)`\
`)`\
`#> Registered S3 method overwritten by 'car':`\
`#>   method           from`\
`#>   na.action.merMod lme4`

#### Prepare Initial Condition

The means, variances, and covariances of the initial true levels and the
slopes are all estimated.

\
`dynr_initial`` ``<-`` ``dynr``::`[`prep.initial`](https://rdrr.io/pkg/dynr/man/prep.initial.html)`(`\
`  values.inistate ``=`` ``mu0``,`\
`  params.inistate ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_1_1"``, ``"mu0_2_1"``, ``"mu0_3_1"``, ``"mu0_4_1"``)``,`\
`  values.inicov ``=`` ``sigma0``,`\
`  params.inicov ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"sigma0_1_1"``, ``"sigma0_2_1"``, ``"sigma0_3_1"``, ``"sigma0_4_1"``,`\
`      ``"sigma0_2_1"``, ``"sigma0_2_2"``, ``"sigma0_3_2"``, ``"sigma0_4_2"``,`\
`      ``"sigma0_3_1"``, ``"sigma0_3_2"``, ``"sigma0_3_3"``, ``"sigma0_4_3"``,`\
`      ``"sigma0_4_1"``, ``"sigma0_4_2"``, ``"sigma0_4_3"``, ``"sigma0_4_4"`\
`    ``)``,`\
`    nrow ``=`` ``4`\
`  ``)`\
`)`

#### Prepare Measurement Model

The loadings are fixed.

\
`dynr_measurement`` ``<-`` ``dynr``::`[`prep.measurement`](https://rdrr.io/pkg/dynr/man/prep.measurement.html)`(`\
`  values.load ``=`` ``lambda``,`\
`  params.load ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``2``, ncol ``=`` ``4``)``,`\
`  state.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"eta_1"``, ``"eta_2"``, ``"eta_3"``, ``"eta_4"``)``,`\
`  obs.names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``)`\
`)`

#### Prepare Dynamic Process

The first two formulas read as “the new true level is the previous level
plus the change”, with the change made of the slope, the self-feedback,
and the coupling. The slopes carry over unchanged.

\
`dynr_dynamics`` ``<-`` ``dynr``::`[`prep.formulaDynamics`](https://rdrr.io/pkg/dynr/man/prep.formulaDynamics.html)`(`\
`  formula ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    ``eta_1`` ``~`` ``eta_1`` ``+`` ``eta_3`` ``+`` ``b_1`` ``*`` ``eta_1`` ``+`` ``g_1_2`` ``*`` ``eta_2``,`\
`    ``eta_2`` ``~`` ``eta_2`` ``+`` ``eta_4`` ``+`` ``b_2`` ``*`` ``eta_2`` ``+`` ``g_2_1`` ``*`` ``eta_1``,`\
`    ``eta_3`` ``~`` ``eta_3``,`\
`    ``eta_4`` ``~`` ``eta_4`\
`  ``)``,`\
`  startval ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    b_1 ``=`` ``b_1``, b_2 ``=`` ``b_2``,`\
`    g_1_2 ``=`` ``g_1_2``, g_2_1 ``=`` ``g_2_1`\
`  ``)``,`\
`  isContinuousTime ``=`` ``FALSE`\
`)`

#### Prepare Process Noise

The process noise is fixed at zero, and the two measurement error
variances are estimated.

\
`dynr_noise`` ``<-`` ``dynr``::`[`prep.noise`](https://rdrr.io/pkg/dynr/man/prep.noise.html)`(`\
`  values.latent ``=`` ``psi``,`\
`  params.latent ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``4``, ncol ``=`` ``4``)``,`\
`  values.observed ``=`` ``theta``,`\
`  params.observed ``=`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`    data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"theta_1_1"``, ``"fixed"``,`\
`      ``"fixed"``, ``"theta_2_2"`\
`    ``)``,`\
`    nrow ``=`` ``2`\
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
`  outfile ``=`` ``"lcs.c"`\
`)`

The plot below shows the model equations with the parameter names.

![](fig-vignettes-lcs-model-lcs-1.png)

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
`#> Original fitted parameters:  -0.2111052 -0.29736 0.1117813 0.1928005 -1.627613 `\
`#> -1.624232 2.0052 0.9912427 0.5111888 0.4009262 0.02999263 0.3074729 0.1540763 `\
`#> 0.06532563 -0.0973672 0.003555894 0.1489334 -1.37127 0.1642432 -1.483718 `\
`#> `\
`#> Transformed fitted parameters:  -0.2111052 -0.29736 0.1117813 0.1928005 `\
`#> 0.1963977 0.1970629 2.0052 0.9912427 0.5111888 0.4009262 1.030447 0.3168345 `\
`#> 0.1587674 0.0673146 1.004641 0.05204266 0.1558132 0.2782583 0.05253442 `\
`#> 0.2581595 `\
`#> `\
`#> Doing end processing`\
`#> Successful trial`\
`#> Total Time: 1.054222 `\
`#> Backend Time: 1.047882`

#### Summary

The summary lists each estimate with its standard error, \\(t\\) value,
confidence interval, and \\(p\\) value.

\
[`summary`](https://rdrr.io/r/base/summary.html)`(``results``)`\
`#> Coefficients:`\
`#>             Estimate Std. Error t value  ci.lower  ci.upper`\
`#> b_1        -0.211105   0.009469 -22.295 -0.229663 -0.192547`\
`#> b_2        -0.297360   0.008052 -36.931 -0.313141 -0.281579`\
`#> g_1_2       0.111781   0.008077  13.840  0.095951  0.127611`\
`#> g_2_1       0.192801   0.009444  20.416  0.174292  0.211310`\
`#> theta_1_1   0.196398   0.003113  63.082  0.190296  0.202500`\
`#> theta_2_2   0.197063   0.003137  62.823  0.190915  0.203211`\
`#> mu0_1_1     2.005200   0.024303  82.507  1.957566  2.052834`\
`#> mu0_2_1     0.991243   0.024135  41.071  0.943939  1.038546`\
`#> mu0_3_1     0.511189   0.020627  24.782  0.470760  0.551618`\
`#> mu0_4_1     0.400926   0.020259  19.790  0.361218  0.440634`\
`#> sigma0_1_1  1.030447   0.036772  28.023  0.958375  1.102519`\
`#> sigma0_2_1  0.316834   0.026613  11.905  0.264674  0.368995`\
`#> sigma0_3_1  0.158767   0.014193  11.186  0.130949  0.186586`\
`#> sigma0_4_1  0.067315   0.014578   4.617  0.038742  0.095887`\
`#> sigma0_2_2  1.004641   0.036331  27.653  0.933434  1.075847`\
`#> sigma0_3_2  0.052043   0.014124   3.685  0.024359  0.079726`\
`#> sigma0_4_2  0.155813   0.013392  11.635  0.129565  0.182061`\
`#> sigma0_3_3  0.278258   0.012877  21.610  0.253021  0.303496`\
`#> sigma0_4_3  0.052534   0.008431   6.231  0.036010  0.069059`\
`#> sigma0_4_4  0.258160   0.010700  24.128  0.237189  0.279130`\
`#>                        Pr(>|t|)    `\
`#> b_1        < 0.0000000000000002 ***`\
`#> b_2        < 0.0000000000000002 ***`\
`#> g_1_2      < 0.0000000000000002 ***`\
`#> g_2_1      < 0.0000000000000002 ***`\
`#> theta_1_1  < 0.0000000000000002 ***`\
`#> theta_2_2  < 0.0000000000000002 ***`\
`#> mu0_1_1    < 0.0000000000000002 ***`\
`#> mu0_2_1    < 0.0000000000000002 ***`\
`#> mu0_3_1    < 0.0000000000000002 ***`\
`#> mu0_4_1    < 0.0000000000000002 ***`\
`#> sigma0_1_1 < 0.0000000000000002 ***`\
`#> sigma0_2_1 < 0.0000000000000002 ***`\
`#> sigma0_3_1 < 0.0000000000000002 ***`\
`#> sigma0_4_1           0.00000196 ***`\
`#> sigma0_2_2 < 0.0000000000000002 ***`\
`#> sigma0_3_2             0.000115 ***`\
`#> sigma0_4_2 < 0.0000000000000002 ***`\
`#> sigma0_3_3 < 0.0000000000000002 ***`\
`#> sigma0_4_3 < 0.0000000000000002 ***`\
`#> sigma0_4_4 < 0.0000000000000002 ***`\
`#> ---`\
`#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`\
`#> `\
`#> -2 log-likelihood value at convergence = 51858.31`\
`#> AIC = 51898.31`\
`#> BIC = 52046.16`

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
is a fixed part (each level carries over and adds its slope) plus the
self-feedback and coupling effects, so `A` is an algebra, the sum of the
two. There is no process noise (`Q` is zero).

\
`feedback`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``0``, nrow ``=`` ``4``, ncol ``=`` ``4``)`\
`feedback``[``1``:``2``, ``1``:``2``]`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``b_1``, ``g_2_1``, ``g_1_2``, ``b_2``)`\
`labels_feedback`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(``data ``=`` ``"fixed"``, nrow ``=`` ``4``, ncol ``=`` ``4``)`\
`labels_feedback``[``1``:``2``, ``1``:``2``]`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"b_1"``, ``"g_2_1"``, ``"g_1_2"``, ``"b_2"``)`\
`labels_sigma0`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"sigma0_1_1"``, ``"sigma0_2_1"``, ``"sigma0_3_1"``, ``"sigma0_4_1"``,`\
`    ``"sigma0_2_1"``, ``"sigma0_2_2"``, ``"sigma0_3_2"``, ``"sigma0_4_2"``,`\
`    ``"sigma0_3_1"``, ``"sigma0_3_2"``, ``"sigma0_3_3"``, ``"sigma0_4_3"``,`\
`    ``"sigma0_4_1"``, ``"sigma0_4_2"``, ``"sigma0_4_3"``, ``"sigma0_4_4"`\
`  ``)``,`\
`  nrow ``=`` ``4`\
`)`\
`mx_matrices`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`  ``MxMat``(``values ``=`` ``beta`` ``-`` ``feedback``, labels ``=`` ``"fixed"``, name ``=`` ``"carry"``)``,`\
`  ``MxMat``(``values ``=`` ``feedback``, labels ``=`` ``labels_feedback``, name ``=`` ``"feedback"``)``,`\
`  ``OpenMx``::`[`mxAlgebra`](https://rdrr.io/pkg/OpenMx/man/mxAlgebra.html)`(``expression ``=`` ``carry`` ``+`` ``feedback``, name ``=`` ``"A"``)``,`\
`  ``MxMat``(``values ``=`` ``alpha``, labels ``=`` ``"fixed"``, name ``=`` ``"B"``)``,`\
`  ``MxMat``(``values ``=`` ``lambda``, labels ``=`` ``"fixed"``, name ``=`` ``"C"``)``,`\
`  ``MxMat``(``values ``=`` ``nu``, labels ``=`` ``"fixed"``, name ``=`` ``"D"``)``,`\
`  ``MxMat``(``values ``=`` ``psi``, labels ``=`` ``"fixed"``, name ``=`` ``"Q"``)``,`\
`  ``MxMat``(`\
`    values ``=`` ``theta``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"theta_1_1"``, ``"fixed"``, ``"fixed"``, ``"theta_2_2"``)``,`\
`    name ``=`` ``"R"`\
`  ``)``,`\
`  ``MxMat``(`\
`    values ``=`` ``mu0``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_1_1"``, ``"mu0_2_1"``, ``"mu0_3_1"``, ``"mu0_4_1"``)``,`\
`    name ``=`` ``"mu0"`\
`  ``)``,`\
`  ``MxMat``(``values ``=`` ``sigma0``, labels ``=`` ``labels_sigma0``, name ``=`` ``"sigma0"``)`\
`)`

``` r
fit_mx_ssm <- FitMxSSM(
  data = data,
  observed = c("y1", "y2"),
  matrices = mx_matrices,
  name = "lcs"
)
#> 
Beginning initial fit attempt
Fit attempt 0, fit=51858.3099900912, new current best! (was 51872.695914374)
                                                                            
#> 
#>  Solution found!  Final fit=51858.31 (started at 51872.696)  (1 attempt(s): 1 valid, 0 errors)
summary(fit_mx_ssm)
#> Summary of lcs 
#>  
#> free parameters:
#>          name   matrix row col    Estimate   Std.Error A
#> 1         b_1 feedback   1   1 -0.21096858 0.009462328  
#> 2       g_2_1 feedback   2   1  0.19289389 0.009440144  
#> 3       g_1_2 feedback   1   2  0.11179913 0.008071380  
#> 4         b_2 feedback   2   2 -0.29733010 0.008047631  
#> 5   theta_1_1        R   1   1  0.19640814 0.003113182  
#> 6   theta_2_2        R   2   2  0.19707365 0.003136681  
#> 7     mu0_1_1      mu0   1   1  2.00533053 0.024301722  
#> 8     mu0_2_1      mu0   2   1  0.99136005 0.024121585  
#> 9     mu0_3_1      mu0   3   1  0.51095027 0.020625726  
#> 10    mu0_4_1      mu0   4   1  0.40072153 0.020266705  
#> 11 sigma0_1_1   sigma0   1   1  1.03052285 0.036723098  
#> 12 sigma0_2_1   sigma0   2   1  0.31692473 0.026607487  
#> 13 sigma0_3_1   sigma0   3   1  0.15875743 0.014194767  
#> 14 sigma0_4_1   sigma0   4   1  0.06730386 0.014581311  
#> 15 sigma0_2_2   sigma0   2   2  1.00468893 0.036313000  
#> 16 sigma0_3_2   sigma0   3   2  0.05206261 0.014123809  
#> 17 sigma0_4_2   sigma0   4   2  0.15580819 0.013390717  
#> 18 sigma0_3_3   sigma0   3   3  0.27814602 0.012862091  
#> 19 sigma0_4_3   sigma0   4   3  0.05246000 0.008425488  
#> 20 sigma0_4_4   sigma0   4   4  0.25812078 0.010696107  
#> 
#> Model Statistics: 
#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
#>        Model:             20                  23980              51858.31
#>    Saturated:             NA                     NA                    NA
#> Independence:             NA                     NA                    NA
#> Number of observations/statistics: 12000/24000
#> 
#> Information Criteria: 
#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
#> AIC:        3898.31               51898.31                 51898.38
#> BIC:     -173377.72               52046.16                 51982.61
#> To get additional fit indices, see help(mxRefModels)
#> timestamp: 2026-11-02 05:49:49 
#> Wall clock time: 50.37064 secs 
#> optimizer:  SLSQP 
#> OpenMx version number: 2.22.11 
#> Need help?  See help(mxSummary)
```

### lavaan: Structural Equation Model

The same model in SEM form, fitted with `lavaan`, should give nearly the
same estimates. In SEM, each variable at each occasion is a separate
observed variable in one row per individual (wide format). Each observed
score has a true-level factor; each true level after the first is the
previous true level plus a change factor; and each change factor is the
slope plus the self-feedback and coupling effects.

#### Wide Data Set

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) with
`long = FALSE` returns the data in wide format.

\
`data_wide`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``, long ``=`` ``FALSE``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data_wide``)`\
`#>   id     y1_0       y2_0     y1_1        y2_1     y1_2        y2_2     y1_3`\
`#> 1  1 1.239124 -0.0347411 2.372287 -0.05027552 2.550748 -0.36316715 3.417814`\
`#> 2  2 1.414760 -1.7045943 1.880544 -0.62495355 2.246019 -0.23016644 1.492737`\
`#> 3  3 1.076995  1.6913778 2.687536  2.38251121 3.147829  3.29079426 3.582296`\
`#> 4  4 3.371977  0.2882310 2.850104  1.22590392 2.501115  1.20157000 2.255159`\
`#> 5  5 2.953999  0.3501444 1.791110 -0.73608084 1.715968  0.05072464 1.850961`\
`#> 6  6 1.076737 -1.5194131 0.735993  1.07098420 1.650132  1.57132134 2.070286`\
`#>         y2_3     y1_4         y2_4     y1_5       y2_5`\
`#> 1 -0.1644062 3.685087  0.007406512 3.473986 -0.3485020`\
`#> 2 -0.7540751 1.117000  0.347719984 1.811812  0.3486327`\
`#> 3  4.0998228 4.878402  4.993348339 4.269622  4.8999756`\
`#> 4  1.7383942 2.659687  1.991782688 2.388665  2.2080090`\
`#> 5 -0.8198842 1.300108 -1.526556346 1.172956 -0.7951067`\
`#> 6  1.7259835 1.860267  1.958970444 1.452103  2.4188041`

#### Model

The model has the same structure for every occasion, so we write it with
[`paste0()`](https://rdrr.io/r/base/paste.html). The parameter labels
match the dynr names. Giving the self-feedback and coupling effects and
the measurement error variances the same label on every occasion
constrains them to be equal over time, as in the state space model.

\
`waves`` ``<-`` ``0``:``(``time`` ``-`` ``1``)`\
`now`` ``<-`` ``waves``[``-``1``]`\
`before`` ``<-`` ``waves``[``-``time``]`\
`model`` ``<-`` `[`paste`](https://rdrr.io/r/base/paste.html)`(`\
`  `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"# true levels"``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l1_"``, ``waves``, ``" =~ 1 * y1_"``, ``waves``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l2_"``, ``waves``, ``" =~ 1 * y2_"``, ``waves``)``,`\
`    ``"# each true level is the previous level plus the change"``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l1_"``, ``now``, ``" ~ 1 * l1_"``, ``before``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l2_"``, ``now``, ``" ~ 1 * l2_"``, ``before``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l1_"``, ``now``, ``" ~~ 0 * l1_"``, ``now``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l2_"``, ``now``, ``" ~~ 0 * l2_"``, ``now``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d1_"``, ``now``, ``" =~ 1 * l1_"``, ``now``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d2_"``, ``now``, ``" =~ 1 * l2_"``, ``now``)``,`\
`    ``"# the change: self-feedback and coupling"``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d1_"``, ``now``, ``" ~ b_1 * l1_"``, ``before``, ``" + g_1_2 * l2_"``, ``before``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d2_"``, ``now``, ``" ~ b_2 * l2_"``, ``before``, ``" + g_2_1 * l1_"``, ``before``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d1_"``, ``now``, ``" ~~ 0 * d1_"``, ``now``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d2_"``, ``now``, ``" ~~ 0 * d2_"``, ``now``)``,`\
`    ``"# the change: slopes"``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"s1 =~ "``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"1 * d1_"``, ``now``, collapse ``=`` ``" + "``)``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"s2 =~ "``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"1 * d2_"``, ``now``, collapse ``=`` ``" + "``)``)``,`\
`    ``"# initial true levels and slopes: means"``,`\
`    ``"l1_0 ~ mu0_1_1 * 1"``,`\
`    ``"l2_0 ~ mu0_2_1 * 1"``,`\
`    ``"s1 ~ mu0_3_1 * 1"``,`\
`    ``"s2 ~ mu0_4_1 * 1"``,`\
`    ``"# initial true levels and slopes: variances and covariances"``,`\
`    ``"l1_0 ~~ sigma0_1_1 * l1_0"``,`\
`    ``"l2_0 ~~ sigma0_2_2 * l2_0"``,`\
`    ``"s1 ~~ sigma0_3_3 * s1"``,`\
`    ``"s2 ~~ sigma0_4_4 * s2"``,`\
`    ``"l1_0 ~~ sigma0_2_1 * l2_0"``,`\
`    ``"l1_0 ~~ sigma0_3_1 * s1"``,`\
`    ``"l1_0 ~~ sigma0_4_1 * s2"``,`\
`    ``"l2_0 ~~ sigma0_3_2 * s1"``,`\
`    ``"l2_0 ~~ sigma0_4_2 * s2"``,`\
`    ``"s1 ~~ sigma0_4_3 * s2"``,`\
`    ``"# measurement error variances, equal over time; no observed intercepts"``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y1_"``, ``waves``, ``" ~~ theta_1_1 * y1_"``, ``waves``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y2_"``, ``waves``, ``" ~~ theta_2_2 * y2_"``, ``waves``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y1_"``, ``waves``, ``" ~ 0 * 1"``)``,`\
`    `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y2_"``, ``waves``, ``" ~ 0 * 1"``)`\
`  ``)``,`\
`  collapse ``=`` ``"\n"`\
`)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``model``)`\
`#> # true levels`\
`#> l1_0 =~ 1 * y1_0`\
`#> l1_1 =~ 1 * y1_1`\
`#> l1_2 =~ 1 * y1_2`\
`#> l1_3 =~ 1 * y1_3`\
`#> l1_4 =~ 1 * y1_4`\
`#> l1_5 =~ 1 * y1_5`\
`#> l2_0 =~ 1 * y2_0`\
`#> l2_1 =~ 1 * y2_1`\
`#> l2_2 =~ 1 * y2_2`\
`#> l2_3 =~ 1 * y2_3`\
`#> l2_4 =~ 1 * y2_4`\
`#> l2_5 =~ 1 * y2_5`\
`#> # each true level is the previous level plus the change`\
`#> l1_1 ~ 1 * l1_0`\
`#> l1_2 ~ 1 * l1_1`\
`#> l1_3 ~ 1 * l1_2`\
`#> l1_4 ~ 1 * l1_3`\
`#> l1_5 ~ 1 * l1_4`\
`#> l2_1 ~ 1 * l2_0`\
`#> l2_2 ~ 1 * l2_1`\
`#> l2_3 ~ 1 * l2_2`\
`#> l2_4 ~ 1 * l2_3`\
`#> l2_5 ~ 1 * l2_4`\
`#> l1_1 ~~ 0 * l1_1`\
`#> l1_2 ~~ 0 * l1_2`\
`#> l1_3 ~~ 0 * l1_3`\
`#> l1_4 ~~ 0 * l1_4`\
`#> l1_5 ~~ 0 * l1_5`\
`#> l2_1 ~~ 0 * l2_1`\
`#> l2_2 ~~ 0 * l2_2`\
`#> l2_3 ~~ 0 * l2_3`\
`#> l2_4 ~~ 0 * l2_4`\
`#> l2_5 ~~ 0 * l2_5`\
`#> d1_1 =~ 1 * l1_1`\
`#> d1_2 =~ 1 * l1_2`\
`#> d1_3 =~ 1 * l1_3`\
`#> d1_4 =~ 1 * l1_4`\
`#> d1_5 =~ 1 * l1_5`\
`#> d2_1 =~ 1 * l2_1`\
`#> d2_2 =~ 1 * l2_2`\
`#> d2_3 =~ 1 * l2_3`\
`#> d2_4 =~ 1 * l2_4`\
`#> d2_5 =~ 1 * l2_5`\
`#> # the change: self-feedback and coupling`\
`#> d1_1 ~ b_1 * l1_0 + g_1_2 * l2_0`\
`#> d1_2 ~ b_1 * l1_1 + g_1_2 * l2_1`\
`#> d1_3 ~ b_1 * l1_2 + g_1_2 * l2_2`\
`#> d1_4 ~ b_1 * l1_3 + g_1_2 * l2_3`\
`#> d1_5 ~ b_1 * l1_4 + g_1_2 * l2_4`\
`#> d2_1 ~ b_2 * l2_0 + g_2_1 * l1_0`\
`#> d2_2 ~ b_2 * l2_1 + g_2_1 * l1_1`\
`#> d2_3 ~ b_2 * l2_2 + g_2_1 * l1_2`\
`#> d2_4 ~ b_2 * l2_3 + g_2_1 * l1_3`\
`#> d2_5 ~ b_2 * l2_4 + g_2_1 * l1_4`\
`#> d1_1 ~~ 0 * d1_1`\
`#> d1_2 ~~ 0 * d1_2`\
`#> d1_3 ~~ 0 * d1_3`\
`#> d1_4 ~~ 0 * d1_4`\
`#> d1_5 ~~ 0 * d1_5`\
`#> d2_1 ~~ 0 * d2_1`\
`#> d2_2 ~~ 0 * d2_2`\
`#> d2_3 ~~ 0 * d2_3`\
`#> d2_4 ~~ 0 * d2_4`\
`#> d2_5 ~~ 0 * d2_5`\
`#> # the change: slopes`\
`#> s1 =~ 1 * d1_1 + 1 * d1_2 + 1 * d1_3 + 1 * d1_4 + 1 * d1_5`\
`#> s2 =~ 1 * d2_1 + 1 * d2_2 + 1 * d2_3 + 1 * d2_4 + 1 * d2_5`\
`#> # initial true levels and slopes: means`\
`#> l1_0 ~ mu0_1_1 * 1`\
`#> l2_0 ~ mu0_2_1 * 1`\
`#> s1 ~ mu0_3_1 * 1`\
`#> s2 ~ mu0_4_1 * 1`\
`#> # initial true levels and slopes: variances and covariances`\
`#> l1_0 ~~ sigma0_1_1 * l1_0`\
`#> l2_0 ~~ sigma0_2_2 * l2_0`\
`#> s1 ~~ sigma0_3_3 * s1`\
`#> s2 ~~ sigma0_4_4 * s2`\
`#> l1_0 ~~ sigma0_2_1 * l2_0`\
`#> l1_0 ~~ sigma0_3_1 * s1`\
`#> l1_0 ~~ sigma0_4_1 * s2`\
`#> l2_0 ~~ sigma0_3_2 * s1`\
`#> l2_0 ~~ sigma0_4_2 * s2`\
`#> s1 ~~ sigma0_4_3 * s2`\
`#> # measurement error variances, equal over time; no observed intercepts`\
`#> y1_0 ~~ theta_1_1 * y1_0`\
`#> y1_1 ~~ theta_1_1 * y1_1`\
`#> y1_2 ~~ theta_1_1 * y1_2`\
`#> y1_3 ~~ theta_1_1 * y1_3`\
`#> y1_4 ~~ theta_1_1 * y1_4`\
`#> y1_5 ~~ theta_1_1 * y1_5`\
`#> y2_0 ~~ theta_2_2 * y2_0`\
`#> y2_1 ~~ theta_2_2 * y2_1`\
`#> y2_2 ~~ theta_2_2 * y2_2`\
`#> y2_3 ~~ theta_2_2 * y2_3`\
`#> y2_4 ~~ theta_2_2 * y2_4`\
`#> y2_5 ~~ theta_2_2 * y2_5`\
`#> y1_0 ~ 0 * 1`\
`#> y1_1 ~ 0 * 1`\
`#> y1_2 ~ 0 * 1`\
`#> y1_3 ~ 0 * 1`\
`#> y1_4 ~ 0 * 1`\
`#> y1_5 ~ 0 * 1`\
`#> y2_0 ~ 0 * 1`\
`#> y2_1 ~ 0 * 1`\
`#> y2_2 ~ 0 * 1`\
`#> y2_3 ~ 0 * 1`\
`#> y2_4 ~ 0 * 1`\
`#> y2_5 ~ 0 * 1`

[`lavaan::lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html) adds
nothing to the model by itself, so every parameter above is exactly the
one written.

\
`fit_lavaan`` ``<-`` ``lavaan``::`[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html)`(`\
`  model ``=`` ``model``,`\
`  data ``=`` ``data_wide`\
`)`\
`lavaan``::``summary``(``fit_lavaan``)`\
`#> lavaan 0.7-3 ended normally after 48 iterations`\
`#> `\
`#>   Estimator                                         ML`\
`#>   Optimization method                           NLMINB`\
`#>   Number of model parameters                        46`\
`#>   Number of equality constraints                    26`\
`#> `\
`#>   Number of observations                          2000`\
`#> `\
`#> Model Test User Model:`\
`#>                                                       `\
`#>   Test statistic                                64.866`\
`#>   Degrees of freedom                                70`\
`#>   P-value (Chi-square)                           0.651`\
`#>                                                       `\
`#>   Browne's residual (NT model-based) test             `\
`#>   Test statistic                                64.525`\
`#>   Degrees of freedom                                70`\
`#>   P-value (Chi-square)                           0.662`\
`#> `\
`#> Parameter Estimates:`\
`#> `\
`#>   Standard errors                             Standard`\
`#>   Information                                 Expected`\
`#>   Information saturated (h1) model          Structured`\
`#> `\
`#> Latent Variables:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>   l1_0 =~                                             `\
`#>     y1_0              1.000                           `\
`#>   l1_1 =~                                             `\
`#>     y1_1              1.000                           `\
`#>   l1_2 =~                                             `\
`#>     y1_2              1.000                           `\
`#>   l1_3 =~                                             `\
`#>     y1_3              1.000                           `\
`#>   l1_4 =~                                             `\
`#>     y1_4              1.000                           `\
`#>   l1_5 =~                                             `\
`#>     y1_5              1.000                           `\
`#>   l2_0 =~                                             `\
`#>     y2_0              1.000                           `\
`#>   l2_1 =~                                             `\
`#>     y2_1              1.000                           `\
`#>   l2_2 =~                                             `\
`#>     y2_2              1.000                           `\
`#>   l2_3 =~                                             `\
`#>     y2_3              1.000                           `\
`#>   l2_4 =~                                             `\
`#>     y2_4              1.000                           `\
`#>   l2_5 =~                                             `\
`#>     y2_5              1.000                           `\
`#>   d1_1 =~                                             `\
`#>     l1_1              1.000                           `\
`#>   d1_2 =~                                             `\
`#>     l1_2              1.000                           `\
`#>   d1_3 =~                                             `\
`#>     l1_3              1.000                           `\
`#>   d1_4 =~                                             `\
`#>     l1_4              1.000                           `\
`#>   d1_5 =~                                             `\
`#>     l1_5              1.000                           `\
`#>   d2_1 =~                                             `\
`#>     l2_1              1.000                           `\
`#>   d2_2 =~                                             `\
`#>     l2_2              1.000                           `\
`#>   d2_3 =~                                             `\
`#>     l2_3              1.000                           `\
`#>   d2_4 =~                                             `\
`#>     l2_4              1.000                           `\
`#>   d2_5 =~                                             `\
`#>     l2_5              1.000                           `\
`#>   s1 =~                                               `\
`#>     d1_1              1.000                           `\
`#>     d1_2              1.000                           `\
`#>     d1_3              1.000                           `\
`#>     d1_4              1.000                           `\
`#>     d1_5              1.000                           `\
`#>   s2 =~                                               `\
`#>     d2_1              1.000                           `\
`#>     d2_2              1.000                           `\
`#>     d2_3              1.000                           `\
`#>     d2_4              1.000                           `\
`#>     d2_5              1.000                           `\
`#> `\
`#> Regressions:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>   l1_1 ~                                              `\
`#>     l1_0              1.000                           `\
`#>   l1_2 ~                                              `\
`#>     l1_1              1.000                           `\
`#>   l1_3 ~                                              `\
`#>     l1_2              1.000                           `\
`#>   l1_4 ~                                              `\
`#>     l1_3              1.000                           `\
`#>   l1_5 ~                                              `\
`#>     l1_4              1.000                           `\
`#>   l2_1 ~                                              `\
`#>     l2_0              1.000                           `\
`#>   l2_2 ~                                              `\
`#>     l2_1              1.000                           `\
`#>   l2_3 ~                                              `\
`#>     l2_2              1.000                           `\
`#>   l2_4 ~                                              `\
`#>     l2_3              1.000                           `\
`#>   l2_5 ~                                              `\
`#>     l2_4              1.000                           `\
`#>   d1_1 ~                                              `\
`#>     l1_0     (b_1)   -0.211    0.009  -22.340    0.000`\
`#>     l2_0    (g_1_)    0.112    0.008   13.849    0.000`\
`#>   d1_2 ~                                              `\
`#>     l1_1     (b_1)   -0.211    0.009  -22.340    0.000`\
`#>     l2_1    (g_1_)    0.112    0.008   13.849    0.000`\
`#>   d1_3 ~                                              `\
`#>     l1_2     (b_1)   -0.211    0.009  -22.340    0.000`\
`#>     l2_2    (g_1_)    0.112    0.008   13.849    0.000`\
`#>   d1_4 ~                                              `\
`#>     l1_3     (b_1)   -0.211    0.009  -22.340    0.000`\
`#>     l2_3    (g_1_)    0.112    0.008   13.849    0.000`\
`#>   d1_5 ~                                              `\
`#>     l1_4     (b_1)   -0.211    0.009  -22.340    0.000`\
`#>     l2_4    (g_1_)    0.112    0.008   13.849    0.000`\
`#>   d2_1 ~                                              `\
`#>     l2_0     (b_2)   -0.297    0.008  -37.207    0.000`\
`#>     l1_0    (g_2_)    0.193    0.009   20.614    0.000`\
`#>   d2_2 ~                                              `\
`#>     l2_1     (b_2)   -0.297    0.008  -37.207    0.000`\
`#>     l1_1    (g_2_)    0.193    0.009   20.614    0.000`\
`#>   d2_3 ~                                              `\
`#>     l2_2     (b_2)   -0.297    0.008  -37.207    0.000`\
`#>     l1_2    (g_2_)    0.193    0.009   20.614    0.000`\
`#>   d2_4 ~                                              `\
`#>     l2_3     (b_2)   -0.297    0.008  -37.207    0.000`\
`#>     l1_3    (g_2_)    0.193    0.009   20.614    0.000`\
`#>   d2_5 ~                                              `\
`#>     l2_4     (b_2)   -0.297    0.008  -37.207    0.000`\
`#>     l1_4    (g_2_)    0.193    0.009   20.614    0.000`\
`#> `\
`#> Covariances:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>   l1_0 ~~                                             `\
`#>     l2_0    (s0_2)    0.317    0.027   11.910    0.000`\
`#>     s1    (s0_3_1)    0.159    0.014   11.190    0.000`\
`#>     s2    (s0_4_1)    0.067    0.015    4.619    0.000`\
`#>   l2_0 ~~                                             `\
`#>     s1    (s0_3_2)    0.052    0.014    3.687    0.000`\
`#>     s2    (s0_4_2)    0.156    0.013   11.635    0.000`\
`#>   s1 ~~                                               `\
`#>     s2    (s0_4_3)    0.052    0.008    6.241    0.000`\
`#> `\
`#> Intercepts:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>     l1_0    (m0_1)    2.005    0.024   82.518    0.000`\
`#>     l2_0    (m0_2)    0.991    0.024   41.092    0.000`\
`#>     s1      (m0_3)    0.511    0.021   24.810    0.000`\
`#>     s2      (m0_4)    0.401    0.020   19.789    0.000`\
`#>    .y1_0              0.000                           `\
`#>    .y1_1              0.000                           `\
`#>    .y1_2              0.000                           `\
`#>    .y1_3              0.000                           `\
`#>    .y1_4              0.000                           `\
`#>    .y1_5              0.000                           `\
`#>    .y2_0              0.000                           `\
`#>    .y2_1              0.000                           `\
`#>    .y2_2              0.000                           `\
`#>    .y2_3              0.000                           `\
`#>    .y2_4              0.000                           `\
`#>    .y2_5              0.000                           `\
`#> `\
`#> Variances:`\
`#>                    Estimate  Std.Err  z-value  P(>|z|)`\
`#>    .l1_1              0.000                           `\
`#>    .l1_2              0.000                           `\
`#>    .l1_3              0.000                           `\
`#>    .l1_4              0.000                           `\
`#>    .l1_5              0.000                           `\
`#>    .l2_1              0.000                           `\
`#>    .l2_2              0.000                           `\
`#>    .l2_3              0.000                           `\
`#>    .l2_4              0.000                           `\
`#>    .l2_5              0.000                           `\
`#>    .d1_1              0.000                           `\
`#>    .d1_2              0.000                           `\
`#>    .d1_3              0.000                           `\
`#>    .d1_4              0.000                           `\
`#>    .d1_5              0.000                           `\
`#>    .d2_1              0.000                           `\
`#>    .d2_2              0.000                           `\
`#>    .d2_3              0.000                           `\
`#>    .d2_4              0.000                           `\
`#>    .d2_5              0.000                           `\
`#>     l1_0    (s0_1)    1.031    0.037   28.069    0.000`\
`#>     l2_0    (s0_2)    1.005    0.036   27.666    0.000`\
`#>     s1      (s0_3)    0.278    0.013   21.652    0.000`\
`#>     s2      (s0_4)    0.258    0.011   24.151    0.000`\
`#>    .y1_0    (t_1_)    0.196    0.003   63.135    0.000`\
`#>    .y1_1    (t_1_)    0.196    0.003   63.135    0.000`\
`#>    .y1_2    (t_1_)    0.196    0.003   63.135    0.000`\
`#>    .y1_3    (t_1_)    0.196    0.003   63.135    0.000`\
`#>    .y1_4    (t_1_)    0.196    0.003   63.135    0.000`\
`#>    .y1_5    (t_1_)    0.196    0.003   63.135    0.000`\
`#>    .y2_0    (t_2_)    0.197    0.003   62.833    0.000`\
`#>    .y2_1    (t_2_)    0.197    0.003   62.833    0.000`\
`#>    .y2_2    (t_2_)    0.197    0.003   62.833    0.000`\
`#>    .y2_3    (t_2_)    0.197    0.003   62.833    0.000`\
`#>    .y2_4    (t_2_)    0.197    0.003   62.833    0.000`\
`#>    .y2_5    (t_2_)    0.197    0.003   62.833    0.000`

### OpenMx: Structural Equation Model

The same SEM in `OpenMx`, written with
[`OpenMx::mxPath()`](https://rdrr.io/pkg/OpenMx/man/mxPath.html) on the
wide data. One-headed arrows are loadings and regressions, two-headed
arrows are variances and covariances, and arrows from `"one"` are means.
A path that is not written is fixed at zero, so the zero residual
variances, covariances, and intercepts of the lavaan model need no lines
here. The labels are the dynr names, and the starting values are the
population values.

`connect = "unique.pairs"` adds every variance and covariance among the
initial true levels and the slopes; their labels and values are the
lower triangle of the label matrix and of
\\(\boldsymbol{\Sigma}\_{0}\\), in the same order.

\
`waves`` ``<-`` ``0``:``(``time`` ``-`` ``1``)`\
`y1`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y1_"``, ``waves``)`\
`y2`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"y2_"``, ``waves``)`\
`l1`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l1_"``, ``waves``)`\
`l2`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"l2_"``, ``waves``)`\
`d1`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d1_"``, ``waves``[``-``1``]``)`\
`d2`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"d2_"``, ``waves``[``-``1``]``)`\
`initial`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"l1_0"``, ``"l2_0"``, ``"s1"``, ``"s2"``)`\
`mx_sem`` ``<-`` ``OpenMx``::`[`mxModel`](https://rdrr.io/pkg/OpenMx/man/mxModel.html)`(`\
`  name ``=`` ``"lcs_sem"``,`\
`  type ``=`` ``"RAM"``,`\
`  manifestVars ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``y1``, ``y2``)``,`\
`  latentVars ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``l1``, ``l2``, ``d1``, ``d2``, ``"s1"``, ``"s2"``)``,`\
`  ``OpenMx``::`[`mxData`](https://rdrr.io/pkg/OpenMx/man/mxData.html)`(``observed ``=`` ``data_wide``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``y1``, ``y2``)``]``, type ``=`` ``"raw"``)``,`\
`  ``# true levels`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l1``, to ``=`` ``y1``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l2``, to ``=`` ``y2``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``# each true level is the previous level plus the change`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l1``[``-``time``]``, to ``=`` ``l1``[``-``1``]``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l2``[``-``time``]``, to ``=`` ``l2``[``-``1``]``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``d1``, to ``=`` ``l1``[``-``1``]``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``d2``, to ``=`` ``l2``[``-``1``]``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``# the change: self-feedback and coupling`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l1``[``-``time``]``, to ``=`` ``d1``, labels ``=`` ``"b_1"``, values ``=`` ``b_1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l2``[``-``time``]``, to ``=`` ``d1``, labels ``=`` ``"g_1_2"``, values ``=`` ``g_1_2``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l2``[``-``time``]``, to ``=`` ``d2``, labels ``=`` ``"b_2"``, values ``=`` ``b_2``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``l1``[``-``time``]``, to ``=`` ``d2``, labels ``=`` ``"g_2_1"``, values ``=`` ``g_2_1``)``,`\
`  ``# the change: slopes`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``"s1"``, to ``=`` ``d1``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(``from ``=`` ``"s2"``, to ``=`` ``d2``, free ``=`` ``FALSE``, values ``=`` ``1``)``,`\
`  ``# initial true levels and slopes: means`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(`\
`    from ``=`` ``"one"``,`\
`    to ``=`` ``initial``,`\
`    labels ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"mu0_1_1"``, ``"mu0_2_1"``, ``"mu0_3_1"``, ``"mu0_4_1"``)``,`\
`    values ``=`` ``mu0`\
`  ``)``,`\
`  ``# initial true levels and slopes: variances and covariances`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(`\
`    from ``=`` ``initial``,`\
`    arrows ``=`` ``2``,`\
`    connect ``=`` ``"unique.pairs"``,`\
`    labels ``=`` ``labels_sigma0``[`[`lower.tri`](https://rdrr.io/r/base/lower.tri.html)`(``labels_sigma0``, diag ``=`` ``TRUE``)``]``,`\
`    values ``=`` ``sigma0``[`[`lower.tri`](https://rdrr.io/r/base/lower.tri.html)`(``sigma0``, diag ``=`` ``TRUE``)``]`\
`  ``)``,`\
`  ``# measurement error variances, equal over time`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(`\
`    from ``=`` ``y1``, arrows ``=`` ``2``, labels ``=`` ``"theta_1_1"``, values ``=`` ``theta``[``1``, ``1``]`\
`  ``)``,`\
`  ``OpenMx``::`[`mxPath`](https://rdrr.io/pkg/OpenMx/man/mxPath.html)`(`\
`    from ``=`` ``y2``, arrows ``=`` ``2``, labels ``=`` ``"theta_2_2"``, values ``=`` ``theta``[``2``, ``2``]`\
`  ``)`\
`)`

\
`fit_mx_sem`` ``<-`` ``OpenMx``::`[`mxRun`](https://rdrr.io/pkg/OpenMx/man/mxRun.html)`(`\
`  model ``=`` ``mx_sem``,`\
`  silent ``=`` ``TRUE`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_mx_sem``)`\
`#> Summary of lcs_sem `\
`#>  `\
`#> free parameters:`\
`#>          name matrix  row  col    Estimate   Std.Error A`\
`#> 1         b_1      A d1_1 l1_0 -0.21096858 0.009462393  `\
`#> 2       g_2_1      A d2_1 l1_0  0.19289389 0.009438994  `\
`#> 3       g_1_2      A d1_1 l2_0  0.11179913 0.008071393  `\
`#> 4         b_2      A d2_1 l2_0 -0.29733010 0.008047536  `\
`#> 5   theta_1_1      S y1_0 y1_0  0.19640814 0.003113180  `\
`#> 6   theta_2_2      S y2_0 y2_0  0.19707365 0.003136680  `\
`#> 7  sigma0_1_1      S l1_0 l1_0  1.03052285 0.036715998  `\
`#> 8  sigma0_2_1      S l1_0 l2_0  0.31692473 0.026608563  `\
`#> 9  sigma0_2_2      S l2_0 l2_0  1.00468893 0.036312085  `\
`#> 10 sigma0_3_1      S l1_0   s1  0.15875743 0.014194454  `\
`#> 11 sigma0_3_2      S l2_0   s1  0.05206261 0.014123396  `\
`#> 12 sigma0_3_3      S   s1   s1  0.27814602 0.012861056  `\
`#> 13 sigma0_4_1      S l1_0   s2  0.06730386 0.014579721  `\
`#> 14 sigma0_4_2      S l2_0   s2  0.15580819 0.013391110  `\
`#> 15 sigma0_4_3      S   s1   s2  0.05246000 0.008425212  `\
`#> 16 sigma0_4_4      S   s2   s2  0.25812078 0.010696046  `\
`#> 17    mu0_1_1      M    1 l1_0  2.00533053 0.024302200  `\
`#> 18    mu0_2_1      M    1 l2_0  0.99136005 0.024125775  `\
`#> 19    mu0_3_1      M    1   s1  0.51095027 0.020625935  `\
`#> 20    mu0_4_1      M    1   s2  0.40072153 0.020261966  `\
`#> `\
`#> Model Statistics: `\
`#>                |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)`\
`#>        Model:             20                  23980              51858.31`\
`#>    Saturated:             90                  23910                    NA`\
`#> Independence:             24                  23976                    NA`\
`#> Number of observations/statistics: 2000/24000`\
`#> `\
`#> Information Criteria: `\
`#>       |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted`\
`#> AIC:        3898.31               51898.31                 51898.73`\
`#> BIC:     -130411.33               52010.33                 51946.79`\
`#> To get additional fit indices, see help(mxRefModels)`\
`#> timestamp: 2026-11-02 05:50:35 `\
`#> Wall clock time: 0.08748722 secs `\
`#> optimizer:  SLSQP `\
`#> OpenMx version number: 2.22.11 `\
`#> Need help?  See help(mxSummary)`

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
`  mu0_1_1 ``=`` ``mu0``[``1``]``,`\
`  mu0_2_1 ``=`` ``mu0``[``2``]``,`\
`  mu0_3_1 ``=`` ``mu0``[``3``]``,`\
`  mu0_4_1 ``=`` ``mu0``[``4``]``,`\
`  sigma0_1_1 ``=`` ``sigma0``[``1``, ``1``]``,`\
`  sigma0_2_1 ``=`` ``sigma0``[``2``, ``1``]``,`\
`  sigma0_3_1 ``=`` ``sigma0``[``3``, ``1``]``,`\
`  sigma0_4_1 ``=`` ``sigma0``[``4``, ``1``]``,`\
`  sigma0_2_2 ``=`` ``sigma0``[``2``, ``2``]``,`\
`  sigma0_3_2 ``=`` ``sigma0``[``3``, ``2``]``,`\
`  sigma0_4_2 ``=`` ``sigma0``[``4``, ``2``]``,`\
`  sigma0_3_3 ``=`` ``sigma0``[``3``, ``3``]``,`\
`  sigma0_4_3 ``=`` ``sigma0``[``4``, ``3``]``,`\
`  sigma0_4_4 ``=`` ``sigma0``[``4``, ``4``]``,`\
`  b_1 ``=`` ``b_1``,`\
`  b_2 ``=`` ``b_2``,`\
`  g_1_2 ``=`` ``g_1_2``,`\
`  g_2_1 ``=`` ``g_2_1``,`\
`  theta_1_1 ``=`` ``theta``[``1``, ``1``]``,`\
`  theta_2_2 ``=`` ``theta``[``2``, ``2``]`\
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

| Parameter  | Population | Estimate |    SE |   Bias |   2.5% |  97.5% | Covers |
|:-----------|-----------:|---------:|------:|-------:|-------:|-------:|:-------|
| mu0_1_1    |       2.00 |    2.005 | 0.024 |  0.005 |  1.958 |  2.053 | TRUE   |
| mu0_2_1    |       1.00 |    0.991 | 0.024 | -0.009 |  0.944 |  1.039 | TRUE   |
| mu0_3_1    |       0.50 |    0.511 | 0.021 |  0.011 |  0.471 |  0.552 | TRUE   |
| mu0_4_1    |       0.40 |    0.401 | 0.020 |  0.001 |  0.361 |  0.441 | TRUE   |
| sigma0_1_1 |       1.00 |    1.030 | 0.037 |  0.030 |  0.958 |  1.103 | TRUE   |
| sigma0_2_1 |       0.30 |    0.317 | 0.027 |  0.017 |  0.265 |  0.369 | TRUE   |
| sigma0_3_1 |       0.15 |    0.159 | 0.014 |  0.009 |  0.131 |  0.187 | TRUE   |
| sigma0_4_1 |       0.05 |    0.067 | 0.015 |  0.017 |  0.039 |  0.096 | TRUE   |
| sigma0_2_2 |       1.00 |    1.005 | 0.036 |  0.005 |  0.933 |  1.076 | TRUE   |
| sigma0_3_2 |       0.05 |    0.052 | 0.014 |  0.002 |  0.024 |  0.080 | TRUE   |
| sigma0_4_2 |       0.15 |    0.156 | 0.013 |  0.006 |  0.130 |  0.182 | TRUE   |
| sigma0_3_3 |       0.25 |    0.278 | 0.013 |  0.028 |  0.253 |  0.303 | FALSE  |
| sigma0_4_3 |       0.05 |    0.053 | 0.008 |  0.003 |  0.036 |  0.069 | TRUE   |
| sigma0_4_4 |       0.25 |    0.258 | 0.011 |  0.008 |  0.237 |  0.279 | TRUE   |
| b_1        |      -0.20 |   -0.211 | 0.009 | -0.011 | -0.230 | -0.193 | TRUE   |
| b_2        |      -0.30 |   -0.297 | 0.008 |  0.003 | -0.313 | -0.282 | TRUE   |
| g_1_2      |       0.10 |    0.112 | 0.008 |  0.012 |  0.096 |  0.128 | TRUE   |
| g_2_1      |       0.20 |    0.193 | 0.009 | -0.007 |  0.174 |  0.211 | TRUE   |
| theta_1_1  |       0.20 |    0.196 | 0.003 | -0.004 |  0.190 |  0.202 | TRUE   |
| theta_2_2  |       0.20 |    0.197 | 0.003 | -0.003 |  0.191 |  0.203 | TRUE   |

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
| mu0_1_1    |       2.00 |    2.005 | 0.024 |  0.005 |  1.958 |  2.053 | TRUE   |
| mu0_2_1    |       1.00 |    0.991 | 0.024 | -0.009 |  0.944 |  1.039 | TRUE   |
| mu0_3_1    |       0.50 |    0.511 | 0.021 |  0.011 |  0.471 |  0.551 | TRUE   |
| mu0_4_1    |       0.40 |    0.401 | 0.020 |  0.001 |  0.361 |  0.440 | TRUE   |
| sigma0_1_1 |       1.00 |    1.031 | 0.037 |  0.031 |  0.959 |  1.102 | TRUE   |
| sigma0_2_1 |       0.30 |    0.317 | 0.027 |  0.017 |  0.265 |  0.369 | TRUE   |
| sigma0_3_1 |       0.15 |    0.159 | 0.014 |  0.009 |  0.131 |  0.187 | TRUE   |
| sigma0_4_1 |       0.05 |    0.067 | 0.015 |  0.017 |  0.039 |  0.096 | TRUE   |
| sigma0_2_2 |       1.00 |    1.005 | 0.036 |  0.005 |  0.934 |  1.076 | TRUE   |
| sigma0_3_2 |       0.05 |    0.052 | 0.014 |  0.002 |  0.024 |  0.080 | TRUE   |
| sigma0_4_2 |       0.15 |    0.156 | 0.013 |  0.006 |  0.130 |  0.182 | TRUE   |
| sigma0_3_3 |       0.25 |    0.278 | 0.013 |  0.028 |  0.253 |  0.303 | FALSE  |
| sigma0_4_3 |       0.05 |    0.052 | 0.008 |  0.002 |  0.036 |  0.069 | TRUE   |
| sigma0_4_4 |       0.25 |    0.258 | 0.011 |  0.008 |  0.237 |  0.279 | TRUE   |
| b_1        |      -0.20 |   -0.211 | 0.009 | -0.011 | -0.230 | -0.192 | TRUE   |
| b_2        |      -0.30 |   -0.297 | 0.008 |  0.003 | -0.313 | -0.282 | TRUE   |
| g_1_2      |       0.10 |    0.112 | 0.008 |  0.012 |  0.096 |  0.128 | TRUE   |
| g_2_1      |       0.20 |    0.193 | 0.009 | -0.007 |  0.174 |  0.211 | TRUE   |
| theta_1_1  |       0.20 |    0.196 | 0.003 | -0.004 |  0.190 |  0.203 | TRUE   |
| theta_2_2  |       0.20 |    0.197 | 0.003 | -0.003 |  0.191 |  0.203 | TRUE   |

Population values and OpenMx state space estimates (one simulated data
set) {.table style="width:99%;"}

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

| Parameter  | Population | Estimate |    SE |   Bias |   2.5% |  97.5% | Covers |
|:-----------|-----------:|---------:|------:|-------:|-------:|-------:|:-------|
| mu0_1_1    |       2.00 |    2.005 | 0.024 |  0.005 |  1.958 |  2.053 | TRUE   |
| mu0_2_1    |       1.00 |    0.991 | 0.024 | -0.009 |  0.944 |  1.039 | TRUE   |
| mu0_3_1    |       0.50 |    0.511 | 0.021 |  0.011 |  0.471 |  0.551 | TRUE   |
| mu0_4_1    |       0.40 |    0.401 | 0.020 |  0.001 |  0.361 |  0.440 | TRUE   |
| sigma0_1_1 |       1.00 |    1.031 | 0.037 |  0.031 |  0.959 |  1.102 | TRUE   |
| sigma0_2_1 |       0.30 |    0.317 | 0.027 |  0.017 |  0.265 |  0.369 | TRUE   |
| sigma0_3_1 |       0.15 |    0.159 | 0.014 |  0.009 |  0.131 |  0.187 | TRUE   |
| sigma0_4_1 |       0.05 |    0.067 | 0.015 |  0.017 |  0.039 |  0.096 | TRUE   |
| sigma0_2_2 |       1.00 |    1.005 | 0.036 |  0.005 |  0.934 |  1.076 | TRUE   |
| sigma0_3_2 |       0.05 |    0.052 | 0.014 |  0.002 |  0.024 |  0.080 | TRUE   |
| sigma0_4_2 |       0.15 |    0.156 | 0.013 |  0.006 |  0.130 |  0.182 | TRUE   |
| sigma0_3_3 |       0.25 |    0.278 | 0.013 |  0.028 |  0.253 |  0.303 | FALSE  |
| sigma0_4_3 |       0.05 |    0.052 | 0.008 |  0.002 |  0.036 |  0.069 | TRUE   |
| sigma0_4_4 |       0.25 |    0.258 | 0.011 |  0.008 |  0.237 |  0.279 | TRUE   |
| b_1        |      -0.20 |   -0.211 | 0.009 | -0.011 | -0.229 | -0.192 | TRUE   |
| b_2        |      -0.30 |   -0.297 | 0.008 |  0.003 | -0.313 | -0.282 | TRUE   |
| g_1_2      |       0.10 |    0.112 | 0.008 |  0.012 |  0.096 |  0.128 | TRUE   |
| g_2_1      |       0.20 |    0.193 | 0.009 | -0.007 |  0.175 |  0.211 | TRUE   |
| theta_1_1  |       0.20 |    0.196 | 0.003 | -0.004 |  0.190 |  0.203 | TRUE   |
| theta_2_2  |       0.20 |    0.197 | 0.003 | -0.003 |  0.191 |  0.203 | TRUE   |

Population values and lavaan estimates (one simulated data set) {.table
style="width:99%;"}

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

| Parameter  | Population | Estimate |    SE |   Bias |   2.5% |  97.5% | Covers |
|:-----------|-----------:|---------:|------:|-------:|-------:|-------:|:-------|
| mu0_1_1    |       2.00 |    2.005 | 0.024 |  0.005 |  1.958 |  2.053 | TRUE   |
| mu0_2_1    |       1.00 |    0.991 | 0.024 | -0.009 |  0.944 |  1.039 | TRUE   |
| mu0_3_1    |       0.50 |    0.511 | 0.021 |  0.011 |  0.471 |  0.551 | TRUE   |
| mu0_4_1    |       0.40 |    0.401 | 0.020 |  0.001 |  0.361 |  0.440 | TRUE   |
| sigma0_1_1 |       1.00 |    1.031 | 0.037 |  0.031 |  0.959 |  1.102 | TRUE   |
| sigma0_2_1 |       0.30 |    0.317 | 0.027 |  0.017 |  0.265 |  0.369 | TRUE   |
| sigma0_3_1 |       0.15 |    0.159 | 0.014 |  0.009 |  0.131 |  0.187 | TRUE   |
| sigma0_4_1 |       0.05 |    0.067 | 0.015 |  0.017 |  0.039 |  0.096 | TRUE   |
| sigma0_2_2 |       1.00 |    1.005 | 0.036 |  0.005 |  0.934 |  1.076 | TRUE   |
| sigma0_3_2 |       0.05 |    0.052 | 0.014 |  0.002 |  0.024 |  0.080 | TRUE   |
| sigma0_4_2 |       0.15 |    0.156 | 0.013 |  0.006 |  0.130 |  0.182 | TRUE   |
| sigma0_3_3 |       0.25 |    0.278 | 0.013 |  0.028 |  0.253 |  0.303 | FALSE  |
| sigma0_4_3 |       0.05 |    0.052 | 0.008 |  0.002 |  0.036 |  0.069 | TRUE   |
| sigma0_4_4 |       0.25 |    0.258 | 0.011 |  0.008 |  0.237 |  0.279 | TRUE   |
| b_1        |      -0.20 |   -0.211 | 0.009 | -0.011 | -0.230 | -0.192 | TRUE   |
| b_2        |      -0.30 |   -0.297 | 0.008 |  0.003 | -0.313 | -0.282 | TRUE   |
| g_1_2      |       0.10 |    0.112 | 0.008 |  0.012 |  0.096 |  0.128 | TRUE   |
| g_2_1      |       0.20 |    0.193 | 0.009 | -0.007 |  0.174 |  0.211 | TRUE   |
| theta_1_1  |       0.20 |    0.196 | 0.003 | -0.004 |  0.190 |  0.203 | TRUE   |
| theta_2_2  |       0.20 |    0.197 | 0.003 | -0.003 |  0.191 |  0.203 | TRUE   |

Population values and OpenMx SEM estimates (one simulated data set)
{.table style="width:99%;"}

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
`    check.names ``=`` ``FALSE`\
`  ``)``,`\
`  digits ``=`` ``3``,`\
`  caption ``=`` ``"Population values and the estimates of every method"`\
`)`

| Parameter  | Population |   dynr | OpenMx (SSM) | lavaan | OpenMx (SEM) |
|:-----------|-----------:|-------:|-------------:|-------:|-------------:|
| mu0_1_1    |       2.00 |  2.005 |        2.005 |  2.005 |        2.005 |
| mu0_2_1    |       1.00 |  0.991 |        0.991 |  0.991 |        0.991 |
| mu0_3_1    |       0.50 |  0.511 |        0.511 |  0.511 |        0.511 |
| mu0_4_1    |       0.40 |  0.401 |        0.401 |  0.401 |        0.401 |
| sigma0_1_1 |       1.00 |  1.030 |        1.031 |  1.031 |        1.031 |
| sigma0_2_1 |       0.30 |  0.317 |        0.317 |  0.317 |        0.317 |
| sigma0_3_1 |       0.15 |  0.159 |        0.159 |  0.159 |        0.159 |
| sigma0_4_1 |       0.05 |  0.067 |        0.067 |  0.067 |        0.067 |
| sigma0_2_2 |       1.00 |  1.005 |        1.005 |  1.005 |        1.005 |
| sigma0_3_2 |       0.05 |  0.052 |        0.052 |  0.052 |        0.052 |
| sigma0_4_2 |       0.15 |  0.156 |        0.156 |  0.156 |        0.156 |
| sigma0_3_3 |       0.25 |  0.278 |        0.278 |  0.278 |        0.278 |
| sigma0_4_3 |       0.05 |  0.053 |        0.052 |  0.052 |        0.052 |
| sigma0_4_4 |       0.25 |  0.258 |        0.258 |  0.258 |        0.258 |
| b_1        |      -0.20 | -0.211 |       -0.211 | -0.211 |       -0.211 |
| b_2        |      -0.30 | -0.297 |       -0.297 | -0.297 |       -0.297 |
| g_1_2      |       0.10 |  0.112 |        0.112 |  0.112 |        0.112 |
| g_2_1      |       0.20 |  0.193 |        0.193 |  0.193 |        0.193 |
| theta_1_1  |       0.20 |  0.196 |        0.196 |  0.196 |        0.196 |
| theta_2_2  |       0.20 |  0.197 |        0.197 |  0.197 |        0.197 |

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

McArdle, J. J. (2009). Latent variable modeling of differences and
changes with longitudinal data. *Annual Review of Psychology*, *60*(1),
577–605. <https://doi.org/10.1146/annurev.psych.60.110707.163612>

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
