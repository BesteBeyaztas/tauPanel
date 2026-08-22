# tauPanel

## Fast Tau Estimation for Fixed- and Random-Effects Panel Data Models

`tauPanel` is an R package for robust estimation of balanced linear panel data models using fast Tau estimators. The package provides implementations for both **fixed-effects (FE)** and **random-effects (RE)** specifications, together with prediction methods and simulation tools for generating clean and contaminated panel data.

The package is designed for panel-data settings in which classical least-squares and generalized least-squares procedures may be sensitive to vertical outliers, leverage points, or combined contamination.

The main procedures are:

- `fast_tau_fe()` — fast Tau estimation for fixed-effects panel models;
- `fast_tau_re()` — fast Tau estimation for random-effects panel models;
- `fast_tau_panel()` — a unified interface for FE and RE estimation;
- `predict()` — structural-mean prediction from fitted fast Tau models;
- `generate_panel_data()` — simulation of clean and contaminated panel data;
- `fast_tau_control()` — default computational settings used by the fast Tau procedures.

---

## Methodological overview

Consider the balanced panel-data model

$$
y_{it}=x_{it}^{\top}\beta+\alpha_i+\varepsilon_{it}, \qquad i=1,\ldots,N,\quad t=1,\ldots,T,
$$

where $x_{it}$ is a vector of regressors, $\beta$ is the common regression parameter, $\alpha_i$ is an individual-specific effect, and $\varepsilon_{it}$ is the idiosyncratic error.

### Fixed effects

For the proposed FE procedure, the default transformation is first differencing,

$$
\Delta y_{it} = \Delta x_{it}^{\top}\beta + \Delta\varepsilon_{it}.
$$

This removes the time-invariant individual effect exactly. Compared with conventional within demeaning, first differencing also limits the propagation of a contaminated unit-time observation to its adjacent transformed observations.

The resulting transformed regression is estimated using a fast Tau procedure with cluster-diverse randomized starting subsamples and iterative reweighted least-squares refinement.

### Random effects

The proposed RE estimator uses a two-stage robust construction:

1. a first-difference fast Tau pilot fit;
2. robust estimation of scale-based variance components from pilot residuals;
3. residual-based observation weights;
4. robust quasi-demeaning of the panel;
5. a final fast Tau fit to the robustly transformed observations.

The RE procedure therefore robustifies both coefficient estimation and the transformation used before the final regression fit.

---

## Installation

The development version of `tauPanel` can be installed directly from GitHub.

Using `remotes`:

```r
install.packages("remotes")
remotes::install_github("BesteBeyaztas/tauPanel")
```

Alternatively, using `devtools`:

```r
install.packages("devtools")
devtools::install_github("BesteBeyaztas/tauPanel")
```

Then load the package with

```r
library(tauPanel)
```

The package requires the `robustbase` package. Required dependencies are normally installed automatically when `tauPanel` is installed from GitHub.

---

## Package manual

A complete PDF reference manual containing the documentation for all exported functions is available in this repository:

**[tauPanel Reference Manual (PDF)](tauPanel_1.0.0.pdf)**

The same documentation is also available from within R. For example,

```r
?fast_tau_panel
?fast_tau_fe
?fast_tau_re
?fast_tau_control
?generate_panel_data
```

and, for prediction methods,

```r
?predict.fast_tau_fe_fit
?predict.fast_tau_re_fit
```

---

## Quick start

### Fixed-effects model

Generate a small FE panel:

```r
set.seed(123)

dat_fe <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  effect_model = "fe",
  contam_scheme = "none"
)
```

Fit the proposed fast Tau FE estimator:

```r
fit_fe <- fast_tau_fe(
  Y = dat_fe$Y,
  X = dat_fe$X
)

fit_fe$beta
```

The default FE implementation uses first differences.

The same estimator can be fitted through the unified interface:

```r
fit_fe <- fast_tau_panel(
  Y = dat_fe$Y,
  X = dat_fe$X,
  effect = "fe"
)
```

---

## Random-effects model

Generate a random-effects panel:

```r
set.seed(123)

dat_re <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  effect_model = "re",
  contam_scheme = "none"
)
```

Fit the proposed fast Tau RE estimator:

```r
fit_re <- fast_tau_re(
  Y = dat_re$Y,
  X = dat_re$X
)

fit_re$beta
```

The fitted object also contains the robust scale-based variance-component quantities and the corresponding quasi-demeaning parameter:

```r
fit_re$sigma2_a
fit_re$sigma2_e
fit_re$theta
```

The unified interface can also be used:

```r
fit_re <- fast_tau_panel(
  Y = dat_re$Y,
  X = dat_re$X,
  effect = "re"
)
```

---

## Prediction

Both FE and RE fits support the standard R `predict()` interface.

For panel-structured predictor data,

```r
pred_fe <- predict(
  fit_fe,
  newdata = dat_fe$X
)

dim(pred_fe)
```

and similarly for the RE estimator,

```r
pred_re <- predict(
  fit_re,
  newdata = dat_re$X
)
```

Predictions are structural-mean predictions of the form

\[
\widehat{\mu}=x^{\top}\widehat{\beta}.
\]

They do not add an estimated or predicted individual-specific effect. This makes the prediction methods directly applicable to new units for which only predictor values are available.

A single predictor vector can also be supplied:

```r
x_new <- c(0.5, -0.2, 1.1, 0.7)

predict(
  fit_fe,
  newdata = x_new
)
```

For several observations, use a matrix or data frame with one row per observation:

```r
X_new <- matrix(
  rnorm(20),
  nrow = 5,
  ncol = 4
)

predict(
  fit_re,
  newdata = X_new
)
```

---

## Simulating contaminated panel data

`generate_panel_data()` can generate both clean and contaminated balanced panels.

The available contamination schemes are:

```text
"none"
"random_vertical"
"block_vertical"
"random_leverage"
"block_leverage"
"random_vertical_leverage"
```

For example, a fixed-effects panel with 10% random vertical contamination can be generated by

```r
set.seed(123)

dat_vertical <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  effect_model = "fe",
  contam_scheme = "random_vertical",
  cont = 0.10
)
```

The contaminated cells and affected units are returned explicitly:

```r
dat_vertical$outlier_cells
dat_vertical$outlier_units
```

A bad-leverage design can be generated by

```r
set.seed(123)

dat_leverage <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  effect_model = "re",
  contam_scheme = "random_leverage",
  cont = 0.10
)
```

Under leverage contamination, the regressor vector observed by the estimator is contaminated while the response is generated from the corresponding uncontaminated regressor vector. The resulting observations are therefore bad-leverage points rather than high-leverage observations that remain on the regression surface.

Joint vertical and leverage contamination can be generated using

```r
dat_joint <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  effect_model = "fe",
  contam_scheme = "random_vertical_leverage",
  cont = 0.10
)
```

---

## Error distributions

The data generator supports several idiosyncratic-error distributions:

```r
distribution = "normal"
distribution = "t"
distribution = "chi"
```

For example,

```r
dat_t <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  distribution = "t",
  effect_model = "fe",
  contam_scheme = "none"
)
```

The chi-squared option can be used for asymmetric-error experiments:

```r
dat_chi <- generate_panel_data(
  N_ = 100,
  T_ = 5,
  K = 4,
  distribution = "chi",
  chi_df = 1,
  effect_model = "fe",
  contam_scheme = "none"
)
```

---

## Computational controls

The default computational configuration can be inspected using

```r
fast_tau_control()
```

which returns

```r
$N
[1] 2000

$kk
[1] 2

$tt
[1] 20

$rr
[1] 2

$approximate
[1] 0

$cluster_start
[1] TRUE

$m_cluster
NULL

$seed
NULL

$fe_transform
[1] "first_difference"
```

The principal numerical controls can also be supplied directly to the fitting functions. For example,

```r
fit_fe <- fast_tau_fe(
  Y = dat_fe$Y,
  X = dat_fe$X,
  N = 2000,
  kk = 2,
  tt = 20,
  cluster_start = TRUE
)
```

With `cluster_start = TRUE`, randomized starting subsamples are constructed from distinct panel units, with at most one transformed observation contributed by each selected unit.

---

## Main functions

| Function | Purpose |
|---|---|
| `fast_tau_panel()` | Unified interface for FE and RE fast Tau estimation |
| `fast_tau_fe()` | Fast Tau estimator for fixed-effects panel models |
| `fast_tau_re()` | Fast Tau estimator for random-effects panel models |
| `fast_tau_control()` | Default computational settings |
| `generate_panel_data()` | Balanced panel-data generator with optional contamination |
| `predict()` | Structural-mean prediction from fitted FE and RE models |

---

## Model specification

The choice between FE and RE is a modeling decision and is **not** made automatically by `tauPanel`.

Use the FE specification when individual-specific effects may be correlated with the regressor history. The RE specification requires the corresponding mean-independence assumption between the individual effects and regressors.

For example,

```r
fit <- fast_tau_panel(
  Y = Y,
  X = X,
  effect = "fe"
)
```

or

```r
fit <- fast_tau_panel(
  Y = Y,
  X = X,
  effect = "re"
)
```

depending on the assumed panel structure.

---

## License

`tauPanel` is distributed under the **GPL-3** license.

---
