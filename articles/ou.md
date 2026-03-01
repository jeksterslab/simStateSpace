# The Ornstein-Uhlenbeck Model

## Model

The measurement model is given by
``` math
\begin{equation}
  \mathbf{y}_{i, t}
  =
  \boldsymbol{\nu}
  +
  \boldsymbol{\Lambda}
  \boldsymbol{\eta}_{i, t}
  +
  \boldsymbol{\varepsilon}_{i, t},
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
\end{equation}
```
where $`\mathbf{y}_{i, t}`$, $`\boldsymbol{\eta}_{i, t}`$, and
$`\boldsymbol{\varepsilon}_{i, t}`$ are random variables and
$`\boldsymbol{\nu}`$, $`\boldsymbol{\Lambda}`$, and
$`\boldsymbol{\Theta}`$ are model parameters. $`\mathbf{y}_{i, t}`$
represents a vector of observed random variables,
$`\boldsymbol{\eta}_{i, t}`$ a vector of latent random variables, and
$`\boldsymbol{\varepsilon}_{i, t}`$ a vector of random measurement
errors, at time $`t`$ and individual $`i`$. $`\boldsymbol{\nu}`$ denotes
a vector of intercepts, $`\boldsymbol{\Lambda}`$ a matrix of factor
loadings, and $`\boldsymbol{\Theta}`$ the covariance matrix of
$`\boldsymbol{\varepsilon}`$.

An alternative representation of the measurement error is given by
``` math
\begin{equation}
  \boldsymbol{\varepsilon}_{i, t}
  =
  \boldsymbol{\Theta}^{\frac{1}{2}}
  \mathbf{z}_{i, t},
  \quad
  \mathrm{with}
  \quad
  \mathbf{z}_{i, t}
  \sim
  \mathcal{N}
  \left(
  \mathbf{0},
  \mathbf{I}
  \right)
\end{equation}
```
where $`\mathbf{z}_{i, t}`$ is a vector of independent standard normal
random variables and
$`\left( \boldsymbol{\Theta}^{\frac{1}{2}} \right) \left( \boldsymbol{\Theta}^{\frac{1}{2}} \right)^{\prime} = \boldsymbol{\Theta}`$
.

The dynamic structure is given by
``` math
\begin{equation}
  \mathrm{d} \boldsymbol{\eta}_{i, t}
  =
  \boldsymbol{\Phi}
  \left(
  \boldsymbol{\eta}_{i, t}
  -
  \boldsymbol{\mu}
  \right)
  \mathrm{d}t
  +
  \boldsymbol{\Sigma}^{\frac{1}{2}}
  \mathrm{d}
  \mathbf{W}_{i, t}
\end{equation}
```
where $`\boldsymbol{\mu}`$ is the long-term mean or equilibrium level,
$`\boldsymbol{\Phi}`$ is the rate of mean reversion, determining how
quickly the variable returns to its mean, $`\boldsymbol{\Sigma}`$ is the
matrix of volatility or randomness in the process, and
$`\mathrm{d}\boldsymbol{W}`$ is a Wiener process or Brownian motion,
which represents random fluctuations.

## Data Generation

### Notation

Let $`t = 1000`$ be the number of time points and $`n = 100`$ be the
number of individuals.

Let the measurement model intecept vector $`\boldsymbol{\nu}`$ be given
by

``` math
\begin{equation}
\boldsymbol{\nu}
=
\left(
\begin{array}{c}
  0 \\
  0 \\
  0 \\
\end{array}
\right) .
\end{equation}
```

Let the factor loadings matrix $`\boldsymbol{\Lambda}`$ be given by

``` math
\begin{equation}
\boldsymbol{\Lambda}
=
\left(
\begin{array}{ccc}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1 \\
\end{array}
\right) .
\end{equation}
```

Let the measurement error covariance matrix $`\boldsymbol{\Theta}`$ be
given by

``` math
\begin{equation}
\boldsymbol{\Theta}
=
\left(
\begin{array}{ccc}
  0.2 & 0 & 0 \\
  0 & 0.2 & 0 \\
  0 & 0 & 0.2 \\
\end{array}
\right) .
\end{equation}
```

Let the initial condition $`\boldsymbol{\eta}_{0}`$ be given by

``` math
\begin{equation}
\boldsymbol{\eta}_{0} \sim \mathcal{N} \left( \boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}, \boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0} \right)
\end{equation}
```

``` math
\begin{equation}
\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}
=
\left(
\begin{array}{c}
  0 \\
  0 \\
  0 \\
\end{array}
\right)
\end{equation}
```

``` math
\begin{equation}
\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}
=
\left(
\begin{array}{ccc}
  0.3425148 & 0.3296023 & 0.0343817 \\
  0.3296023 & 0.5664625 & 0.2545955 \\
  0.0343817 & 0.2545955 & 0.2999909 \\
\end{array}
\right) .
\end{equation}
```

Let the long-term mean vector $`\boldsymbol{\mu}`$ be given by

``` math
\begin{equation}
\boldsymbol{\mu}
=
\left(
\begin{array}{c}
  0 \\
  0 \\
  0 \\
\end{array}
\right) .
\end{equation}
```

Let the rate of mean reversion matrix $`\boldsymbol{\Phi}`$ be given by

``` math
\begin{equation}
\boldsymbol{\Phi}
=
\left(
\begin{array}{ccc}
  -0.357 & 0 & 0 \\
  0.771 & -0.511 & 0 \\
  -0.45 & 0.729 & -0.693 \\
\end{array}
\right) .
\end{equation}
```

Let the dynamic process noise covariance matrix $`\boldsymbol{\Sigma}`$
be given by

``` math
\begin{equation}
\boldsymbol{\Sigma}
=
\left(
\begin{array}{ccc}
  0.2445556 & 0.0220159 & -0.0500476 \\
  0.0220159 & 0.070678 & 0.0153946 \\
  -0.0500476 & 0.0153946 & 0.0755306 \\
\end{array}
\right) .
\end{equation}
```

Let $`\Delta t = 0.1`$.

### R Function Arguments

``` r

n
#> [1] 100
time
#> [1] 1000
delta_t
#> [1] 0.1
mu0
#> [1] 0 0 0
sigma0
#>           [,1]      [,2]      [,3]
#> [1,] 0.3425148 0.3296023 0.0343817
#> [2,] 0.3296023 0.5664625 0.2545955
#> [3,] 0.0343817 0.2545955 0.2999909
sigma0_l # sigma0_l <- t(chol(sigma0))
#>            [,1]      [,2]      [,3]
#> [1,] 0.58524763 0.0000000 0.0000000
#> [2,] 0.56318429 0.4992855 0.0000000
#> [3,] 0.05874726 0.4436540 0.3157701
mu
#> [1] 0 0 0
phi
#>        [,1]   [,2]   [,3]
#> [1,] -0.357  0.000  0.000
#> [2,]  0.771 -0.511  0.000
#> [3,] -0.450  0.729 -0.693
sigma
#>             [,1]       [,2]        [,3]
#> [1,]  0.24455556 0.02201587 -0.05004762
#> [2,]  0.02201587 0.07067800  0.01539456
#> [3,] -0.05004762 0.01539456  0.07553061
sigma_l # sigma_l <- t(chol(sigma))
#>             [,1]      [,2]     [,3]
#> [1,]  0.49452559 0.0000000 0.000000
#> [2,]  0.04451917 0.2620993 0.000000
#> [3,] -0.10120330 0.0759256 0.243975
nu
#> [1] 0 0 0
lambda
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
theta
#>      [,1] [,2] [,3]
#> [1,]  0.2  0.0  0.0
#> [2,]  0.0  0.2  0.0
#> [3,]  0.0  0.0  0.2
theta_l # theta_l <- t(chol(theta))
#>           [,1]      [,2]      [,3]
#> [1,] 0.4472136 0.0000000 0.0000000
#> [2,] 0.0000000 0.4472136 0.0000000
#> [3,] 0.0000000 0.0000000 0.4472136
```

### Visualizing the Dynamics Without Measurement Error and Process Noise (n = 5 with Different Initial Condition)

![](fig-vignettes-ou-no-error-ou-1.png)![](fig-vignettes-ou-no-error-ou-2.png)![](fig-vignettes-ou-no-error-ou-3.png)

### Using the `SimSSMOUFixed` Function from the `simStateSpace` Package to Simulate Data

``` r

library(simStateSpace)
sim <- SimSSMOUFixed(
  n = n,
  time = time,
  delta_t = delta_t,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  mu = mu,
  phi = phi,
  sigma_l = sigma_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l,
  type = 0
)
data <- as.data.frame(sim)
head(data)
#>   id time          y1         y2         y3
#> 1  1  0.0  0.45047365 -1.2185550  0.6708638
#> 2  1  0.1 -0.84190461  0.1242594  0.1813493
#> 3  1  0.2  0.47289651 -0.2398710  0.5995107
#> 4  1  0.3  0.04089817 -0.6547109  0.4674929
#> 5  1  0.4 -1.37222946 -0.2010512 -0.3341885
#> 6  1  0.5 -0.62410106  0.5262251 -0.1892906
summary(data)
#>        id              time             y1                  y2           
#>  Min.   :  1.00   Min.   : 0.00   Min.   :-3.387237   Min.   :-3.561919  
#>  1st Qu.: 25.75   1st Qu.:24.98   1st Qu.:-0.497505   1st Qu.:-0.593119  
#>  Median : 50.50   Median :49.95   Median : 0.002308   Median : 0.002867  
#>  Mean   : 50.50   Mean   :49.95   Mean   : 0.002420   Mean   : 0.001311  
#>  3rd Qu.: 75.25   3rd Qu.:74.92   3rd Qu.: 0.507538   3rd Qu.: 0.600379  
#>  Max.   :100.00   Max.   :99.90   Max.   : 3.381345   Max.   : 3.415473  
#>        y3            
#>  Min.   :-2.9091510  
#>  1st Qu.:-0.4779403  
#>  Median : 0.0006008  
#>  Mean   : 0.0019631  
#>  3rd Qu.: 0.4815267  
#>  Max.   : 3.0074722
plot(sim)
```

![](fig-vignettes-ou-error-ou-1.png)![](fig-vignettes-ou-error-ou-2.png)![](fig-vignettes-ou-error-ou-3.png)

## Model Fitting

### Prepare Data

``` r

dynr_data <- dynr::dynr.data(
  dataframe = data,
  id = "id",
  time = "time",
  observed = c("y1", "y2", "y3")
)
```

### Prepare Initial Condition

``` r

dynr_initial <- dynr::prep.initial(
  values.inistate = mu0,
  params.inistate = c("mu0_1_1", "mu0_2_1", "mu0_3_1"),
  values.inicov = sigma0,
  params.inicov = matrix(
    data = c(
      "sigma0_1_1", "sigma0_2_1", "sigma0_3_1",
      "sigma0_2_1", "sigma0_2_2", "sigma0_3_2",
      "sigma0_3_1", "sigma0_3_2", "sigma0_3_3"
    ),
    nrow = 3
  )
)
```

### Prepare Measurement Model

``` r

dynr_measurement <- dynr::prep.measurement(
  values.load = diag(3),
  params.load = matrix(data = "fixed", nrow = 3, ncol = 3),
  state.names = c("eta_1", "eta_2", "eta_3"),
  obs.names = c("y1", "y2", "y3")
)
```

### Prepare Dynamic Process

``` r

dynr_dynamics <- dynr::prep.formulaDynamics(
  formula = list(  
    eta_1 ~ (phi_1_1 * (eta_1 - mu_1_1)) + (phi_1_2 * (eta_2 - mu_2_1)) + (phi_1_3 * (eta_3 - mu_3_1)),
    eta_2 ~ (phi_2_1 * (eta_1 - mu_1_1)) + (phi_2_2 * (eta_2 - mu_2_1)) + (phi_2_3 * (eta_3 - mu_3_1)),
    eta_3 ~ (phi_3_1 * (eta_1 - mu_1_1)) + (phi_3_2 * (eta_2 - mu_2_1)) + (phi_3_3 * (eta_3 - mu_3_1))
  ),
  startval = c(
    mu_1_1 = mu[1], mu_2_1 = mu[2], mu_3_1 = mu[3],
    phi_1_1 = phi[1, 1], phi_1_2 = phi[1, 2], phi_1_3 = phi[1, 3],
    phi_2_1 = phi[2, 1], phi_2_2 = phi[2, 2], phi_2_3 = phi[2, 3],
    phi_3_1 = phi[3, 1], phi_3_2 = phi[3, 2], phi_3_3 = phi[3, 3]
  ),
  isContinuousTime = TRUE
)
```

### Prepare Process Noise

``` r

dynr_noise <- dynr::prep.noise(
  values.latent = sigma,
  params.latent = matrix(
    data = c(
      "sigma_1_1", "sigma_2_1", "sigma_3_1",
      "sigma_2_1", "sigma_2_2", "sigma_3_2",
      "sigma_3_1", "sigma_3_2", "sigma_3_3"
    ),
    nrow = 3
  ),
  values.observed = theta,
  params.observed = matrix(
    data = c(
      "theta_1_1", "fixed", "fixed",
      "fixed", "theta_2_2", "fixed",
      "fixed", "fixed", "theta_3_3"
    ),
    nrow = 3
  )
)
```

### Prepare the Model

``` r

model <- dynr::dynr.model(
  data = dynr_data,
  initial = dynr_initial,
  measurement = dynr_measurement,
  dynamics = dynr_dynamics,
  noise = dynr_noise,
  outfile = "ou.c"
)
```

Add lower and upper bounds to aid in the optimization.

``` r

model$lb[
  c(
    "phi_1_1",
    "phi_1_2",
    "phi_1_3",
    "phi_2_1",
    "phi_2_2",
    "phi_2_3",
    "phi_3_1",
    "phi_3_2",
    "phi_3_3"
  )
] <- -1.5
model$ub[
  c(
    "phi_1_1",
    "phi_1_2",
    "phi_1_3",
    "phi_2_1",
    "phi_2_2",
    "phi_2_3",
    "phi_3_1",
    "phi_3_2",
    "phi_3_3"
  )
] <- +1.5
```

![](fig-vignettes-ou-model-ou-1.png)

### Fit the Model

``` r

results <- dynr::dynr.cook(
  model,
  debug_flag = TRUE,
  verbose = FALSE
)
#> [1] "Get ready!!!!"
#> using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0’
#> Optimization function called.
#> Starting Hessian calculation ...
#> Finished Hessian calculation.
#> Original exit flag:  3 
#> Modified exit flag:  3 
#> Optimization terminated successfully: ftol_rel or ftol_abs was reached. 
#> Original fitted parameters:  -0.001520639 -0.003533444 0.001420943 -0.3877294 
#> 0.04975107 -0.0506851 0.7574043 -0.5019982 -0.001077431 -0.448398 0.7166367 
#> -0.6816658 -1.397795 0.08919241 -0.2065966 -2.670757 0.2874013 -2.880175 
#> -1.615792 -1.611789 -1.602981 0.01208934 -0.03468598 -0.01986509 -1.022329 
#> 0.9899225 0.1954094 -1.492165 0.8768935 -2.109573 
#> 
#> Transformed fitted parameters:  -0.001520639 -0.003533444 0.001420943 
#> -0.3877294 0.04975107 -0.0506851 0.7574043 -0.5019982 -0.001077431 -0.448398 
#> 0.7166367 -0.6816658 0.2471414 0.02204313 -0.05105857 0.07116589 0.01533408 
#> 0.07238933 0.1987332 0.1995303 0.2012956 0.01208934 -0.03468598 -0.01986509 
#> 0.3597561 0.3561307 0.07029972 0.577427 0.2667917 0.3079507 
#> 
#> Doing end processing
#> Successful trial
#> Total Time: 3.002938 
#> Backend Time: 3.002721
```

## Summary

``` r

summary(results)
#> Coefficients:
#>             Estimate Std. Error t value  ci.lower  ci.upper Pr(>|t|)    
#> mu_1_1     -0.001521   0.014097  -0.108 -0.029150  0.026109   0.4570    
#> mu_2_1     -0.003533   0.022601  -0.156 -0.047830  0.040763   0.4379    
#> mu_3_1      0.001421   0.014367   0.099 -0.026737  0.029579   0.4606    
#> phi_1_1    -0.387729   0.040455  -9.584 -0.467021 -0.308438   <2e-16 ***
#> phi_1_2     0.049751   0.035494   1.402 -0.019816  0.119318   0.0805 .  
#> phi_1_3    -0.050685   0.027111  -1.870 -0.103822  0.002452   0.0308 *  
#> phi_2_1     0.757404   0.024717  30.644  0.708961  0.805848   <2e-16 ***
#> phi_2_2    -0.501998   0.021768 -23.061 -0.544663 -0.459333   <2e-16 ***
#> phi_2_3    -0.001077   0.016405  -0.066 -0.033230  0.031075   0.4738    
#> phi_3_1    -0.448398   0.025940 -17.286 -0.499239 -0.397557   <2e-16 ***
#> phi_3_2     0.716637   0.022902  31.291  0.671749  0.761524   <2e-16 ***
#> phi_3_3    -0.681666   0.017440 -39.086 -0.715848 -0.647484   <2e-16 ***
#> sigma_1_1   0.247141   0.007205  34.301  0.233020  0.261263   <2e-16 ***
#> sigma_2_1   0.022043   0.002758   7.993  0.016638  0.027448   <2e-16 ***
#> sigma_3_1  -0.051059   0.002901 -17.598 -0.056745 -0.045372   <2e-16 ***
#> sigma_2_2   0.071166   0.001955  36.404  0.067334  0.074997   <2e-16 ***
#> sigma_3_2   0.015334   0.001389  11.038  0.012611  0.018057   <2e-16 ***
#> sigma_3_3   0.072389   0.002103  34.422  0.068268  0.076511   <2e-16 ***
#> theta_1_1   0.198733   0.001181 168.257  0.196418  0.201048   <2e-16 ***
#> theta_2_2   0.199530   0.001002 199.104  0.197566  0.201494   <2e-16 ***
#> theta_3_3   0.201296   0.001017 197.979  0.199303  0.203288   <2e-16 ***
#> mu0_1_1     0.012089   0.066673   0.181 -0.118588  0.142767   0.4281    
#> mu0_2_1    -0.034686   0.082438  -0.421 -0.196262  0.126890   0.3370    
#> mu0_3_1    -0.019865   0.060762  -0.327 -0.138957  0.099227   0.3719    
#> sigma0_1_1  0.359756   0.062604   5.747  0.237054  0.482458   <2e-16 ***
#> sigma0_2_1  0.356131   0.062098   5.735  0.234421  0.477841   <2e-16 ***
#> sigma0_3_1  0.070300   0.040435   1.739 -0.008950  0.149550   0.0411 *  
#> sigma0_2_2  0.577427   0.095535   6.044  0.390182  0.764672   <2e-16 ***
#> sigma0_3_2  0.266792   0.056900   4.689  0.155269  0.378314   <2e-16 ***
#> sigma0_3_3  0.307951   0.049621   6.206  0.210696  0.405205   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> -2 log-likelihood value at convergence = 429162.30
#> AIC = 429222.30
#> BIC = 429507.68
```

### Parameter Estimates

``` r

mu_hat
#> [1] -0.001520639 -0.003533444  0.001420943
phi_hat
#>            [,1]        [,2]         [,3]
#> [1,] -0.3877294  0.04975107 -0.050685104
#> [2,]  0.7574043 -0.50199815 -0.001077431
#> [3,] -0.4483980  0.71663668 -0.681665775
sigma_hat
#>             [,1]       [,2]        [,3]
#> [1,]  0.24714136 0.02204313 -0.05105857
#> [2,]  0.02204313 0.07116589  0.01533408
#> [3,] -0.05105857 0.01533408  0.07238933
theta_hat
#>           [,1]      [,2]      [,3]
#> [1,] 0.1987332 0.0000000 0.0000000
#> [2,] 0.0000000 0.1995303 0.0000000
#> [3,] 0.0000000 0.0000000 0.2012956
mu0_hat
#> [1]  0.01208934 -0.03468598 -0.01986509
sigma0_hat
#>            [,1]      [,2]       [,3]
#> [1,] 0.35975612 0.3561307 0.07029972
#> [2,] 0.35613067 0.5774270 0.26679167
#> [3,] 0.07029972 0.2667917 0.30795070
beta_var1_hat <- expm::expm(phi_hat)
beta_var1_hat
#>            [,1]       [,2]        [,3]
#> [1,]  0.6951836 0.02135575 -0.03008751
#> [2,]  0.4900706 0.61426512 -0.01204410
#> [3,] -0.1042629 0.39412620  0.50932993
```

## References

Chow, S.-M., Ho, M. R., Hamaker, E. L., & Dolan, C. V. (2010).
Equivalence and differences between structural equation modeling and
state-space modeling techniques. *Structural Equation Modeling: A
Multidisciplinary Journal*, *17*(2), 303–332.
<https://doi.org/10.1080/10705511003661553>

Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
Continuous-time dynamic models: Connections to structural equation
models and other discrete-time models. In R. H. Hoyle (Ed.), *Handbook
of structural equation modeling* (2nd ed.). The Guilford Press.

Deboeck, P. R., & Preacher, K. J. (2015). No need to be discrete: A
method for continuous time mediation analysis. *Structural Equation
Modeling: A Multidisciplinary Journal*, *23*(1), 61–75.
<https://doi.org/10.1080/10705511.2014.973960>

Oravecz, Z., Tuerlinckx, F., & Vandekerckhove, J. (2011). A hierarchical
latent stochastic differential equation model for affective dynamics.
*Psychological Methods*, *16*(4), 468–490.
<https://doi.org/10.1037/a0024375>

Ou, L., Hunter, M. D., & Chow, S.-M. (2019). What’s for dynr: A package
for linear and nonlinear dynamic modeling in R. *The R Journal*,
*11*(1), 91. <https://doi.org/10.32614/rj-2019-012>

R Core Team. (2025). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>

Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the
brownian motion. *Physical Review*, *36*(5), 823–841.
<https://doi.org/10.1103/physrev.36.823>
