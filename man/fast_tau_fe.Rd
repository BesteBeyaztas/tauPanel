\name{fast_tau_fe}
\alias{fast_tau_fe}

\title{Fast Tau Estimation for Fixed-Effects Panel Models}

\description{
Fits the proposed robust fast Tau estimator for a balanced
fixed-effects panel regression model. By default, the estimator is
applied to first-differenced observations, although ordinary
within-mean and within-median transformations are also available.
}

\usage{
fast_tau_fe(
  Y,
  X,
  transform = c("first_difference", "within_mean", "within_median"),
  N = 2000,
  kk = 2,
  tt = 20,
  rr = 2,
  approximate = 0,
  cluster_start = TRUE,
  m_cluster = NULL,
  seed = NULL
)
}

\arguments{

  \item{Y}{
A numeric \eqn{n \times T} matrix containing the panel responses,
where \eqn{n} is the number of cross-sectional units and \eqn{T} is
the number of time periods. The panel must be balanced.
}

  \item{X}{
A numeric \eqn{n \times T \times p} array containing the regressors,
where \eqn{p} is the number of regression coefficients. The first
two dimensions must agree with the dimensions of \code{Y}.
}

  \item{transform}{
Character string specifying the transformation used to eliminate
the individual fixed effects. Available choices are
\code{"first_difference"}, \code{"within_mean"}, and
\code{"within_median"}. The default is
\code{"first_difference"}, which corresponds to the proposed
fixed-effects estimator.
}

  \item{N}{
Number of randomized starting subsamples used by the fast-Tau
algorithm. The default is \code{2000}. This quantity corresponds
to the number of random starts, denoted by \eqn{S} in the
methodological description, and should not be confused with the
number of cross-sectional units.
}

  \item{kk}{
Number of preliminary iteratively reweighted least-squares (IRLS)
updates performed for each random or diagnostic starting candidate.
The default is \code{2}, corresponding to \eqn{\kappa = 2}.
}

  \item{tt}{
Number of the most promising candidate solutions retained for
subsequent refinement. The default is \code{20}, corresponding to
the shortlist size \eqn{\ell = 20}.
}

  \item{rr}{
Number of fixed-point scale updates used when an approximate
residual-scale calculation is requested. The default is \code{2}.
This argument is relevant when \code{approximate} is nonzero.
}

  \item{approximate}{
Indicator controlling calculation of the residual M-scale during
candidate evaluation. The default \code{0} uses the full scale
calculation implemented by the fitting routine. A nonzero value
uses the shorter approximation based on \code{rr} scale updates.
}

  \item{cluster_start}{
Logical value indicating whether cluster-diverse starting
subsamples are used. The default is \code{TRUE}. When enabled,
no cross-sectional unit contributes more than one transformed
observation to a starting subsample.
}

  \item{m_cluster}{
Number of distinct cross-sectional units used in each
cluster-diverse starting subsample. If \code{NULL}, the number
of units is set equal to the number of regression coefficients,
yielding a cluster-diverse elemental start. Larger values produce
overdetermined cluster-diverse starting subsets, subject to the
number of available units.
}

  \item{seed}{
Optional integer seed for the randomized starting-subsample
generation. The default is \code{NULL}. The seed is used locally
within the fast-Tau computation.
}

}

\details{

Consider the balanced fixed-effects panel regression model
\deqn{
y_{it} = x_{it}^{\mathsf T}\beta + \alpha_i + \varepsilon_{it},
\qquad
i=1,\ldots,n,\quad t=1,\ldots,T,
}
where \eqn{\beta} is the vector of regression coefficients,
\eqn{\alpha_i} is an unobserved time-invariant individual effect,
and \eqn{\varepsilon_{it}} is the idiosyncratic error.

The proposed fixed-effects estimator uses first differencing by
default. Defining
\eqn{\Delta y_{it}=y_{it}-y_{i,t-1}} and
\eqn{\Delta x_{it}=x_{it}-x_{i,t-1}} gives
\deqn{
\Delta y_{it}
=
\Delta x_{it}^{\mathsf T}\beta
+
\Delta\varepsilon_{it},
\qquad t=2,\ldots,T,
}
so that the individual effect \eqn{\alpha_i} is removed exactly.

First differencing is particularly useful from a robustness
perspective. A contaminated observation at an interior time point
can affect at most the two adjacent first differences, whereas
ordinary within demeaning can transmit the effect of a contaminated
observation to every transformed observation belonging to the same
cross-sectional unit. Moreover, when the level idiosyncratic errors
are identically distributed and independent within a unit, their
differences are marginally symmetric about zero even when the
original error distribution is asymmetric. Adjacent differences are,
however, dependent because they share a common level error.

After transformation, the regression coefficients are estimated by
minimizing a Tau-scale of the transformed residuals. For a candidate
coefficient vector \eqn{\beta}, let
\eqn{e_j(\beta)} denote the corresponding transformed residuals.
The Tau-scale has the general form
\deqn{
\tau^2(\beta)
=
s^2(\beta)
\left\{
\frac{1}{m b_2}
\sum_{j=1}^{m}
\rho_2\left(
\frac{e_j(\beta)}{s(\beta)}
\right)
\right\},
}
where \eqn{s(\beta)} is a robust M-scale determined by a bounded
loss function \eqn{\rho_1}, and \eqn{\rho_2} is a second bounded
loss function used to control efficiency. The estimator seeks a
coefficient vector having a small Tau-scale.

The implementation uses the optimal bounded loss family and the
standard fast-Tau tuning adopted in the accompanying methodology.
The underlying calibration combines a high-breakdown M-scale with
approximately 95 percent Gaussian efficiency in the classical
independent-error regression setting. Because first differencing
induces within-unit dependence, this percentage is a tuning
calibration and should not be interpreted as an exact 95 percent
efficiency statement for the panel estimator.

The Tau objective is non-convex. To explore its relevant local
solutions efficiently, the routine uses multiple starting
subsamples, a small number of preliminary IRLS updates, candidate
screening and shortlisting, followed by more extensive refinement
of the most promising solutions. With \code{cluster_start = TRUE},
starting subsamples are cluster-diverse: distinct units are sampled
and at most one transformed observation is selected from each unit.
This respects the grouping of the transformed panel observations
and reduces the possibility that a starting candidate is dominated
by several dependent observations from one unit.

The optional transformations have the following interpretations:

\itemize{
  \item \code{"first_difference"} uses
  \eqn{y_{it}-y_{i,t-1}} and
  \eqn{x_{it}-x_{i,t-1}}. This is the default and the transformation
  used for the proposed fixed-effects estimator.

  \item \code{"within_mean"} subtracts the ordinary unit-specific
  mean from the response and each regressor.

  \item \code{"within_median"} subtracts the unit-specific median
  from the response and each regressor, providing a robust
  alternative to mean centering.
}

The last two choices are provided as alternative transformation
options. The methodology proposed for the fixed-effects estimator
uses \code{"first_difference"}.

Before applying the fast-Tau algorithm, the transformed regressor
columns are internally rescaled to improve numerical conditioning.
The estimated coefficients are subsequently transformed back to
the original regressor scale, so this rescaling does not alter the
regression model represented by the returned coefficient estimates.
}

\value{

An object of classes \code{"fast_tau_fe_fit"} and
\code{"fast_tau_panel_fit"}. It is a list containing:

\item{beta}{
Numeric vector containing the estimated fixed-effects regression
coefficients on the original regressor scale.
}

\item{fit}{
The underlying fast-Tau fit object. This contains the low-level
optimization results, including the robust coefficient and scale
quantities used by the fitting algorithm.
}

\item{effect}{
Character string equal to \code{"fixed"}, identifying the fitted
model as a fixed-effects specification.
}

\item{transform}{
Character string identifying the transformation used in fitting.
}

\item{transformed_data}{
A list containing the transformed response vector \code{y},
transformed regressor matrix \code{x}, and cross-sectional unit
labels \code{id} supplied to the fast-Tau regression engine.
}

}
\examples{

## Generate a small fixed-effects panel
dat <- generate_panel_data(
  N_ = 60,
  T_ = 4,
  K = 4,
  effect_model = "fe",
  contam_scheme = "none"
)

## Fit the proposed fixed-effects estimator.
fit <- fast_tau_fe(
  Y = dat$Y,
  X = dat$X,
  N = 100,
  tt = 5
)

## Estimated regression coefficients
fit$beta

## Transformation used
fit$transform

## Structural-mean predictions for the supplied predictors
pred <- predict(fit, newdata = dat$X)
dim(pred)
}
