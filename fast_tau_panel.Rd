\name{fast_tau_panel}
\alias{fast_tau_panel}

\title{Fast Tau Estimation for Panel Data Models}

\description{
Provides a unified interface for fitting the proposed fast Tau estimator
under either a fixed-effects or random-effects panel-data specification.
Depending on \code{effect}, the function dispatches to
\code{\link{fast_tau_fe}} or \code{\link{fast_tau_re}}.
}

\usage{
fast_tau_panel(
  Y,
  X,
  effect = c("fe", "re"),
  ...
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
two dimensions must agree with those of \code{Y}.
}

  \item{effect}{
Character string specifying the panel-data specification. Use
\code{"fe"} for the fixed-effects fast Tau estimator or
\code{"re"} for the random-effects fast Tau estimator.
The default is \code{"fe"}.
}

  \item{...}{
Additional arguments passed to the selected fitting function.
For \code{effect = "fe"}, these arguments are passed to
\code{\link{fast_tau_fe}}; for \code{effect = "re"}, they are
passed to \code{\link{fast_tau_re}}. These may include computational
controls such as the number of random starts, preliminary IRLS
iterations, shortlist size, cluster-diverse starting options, and
the random seed.
}

}

\details{

Consider the balanced linear panel-data model
\deqn{
y_{it}
=
x_{it}^{\mathsf T}\beta
+
\alpha_i
+
\varepsilon_{it},
\qquad
i=1,\ldots,n,\quad
t=1,\ldots,T,
}
where \eqn{\beta} is the common regression coefficient vector,
\eqn{\alpha_i} is an unobserved individual-specific effect, and
\eqn{\varepsilon_{it}} is the idiosyncratic error.

The appropriate treatment of \eqn{\alpha_i} depends on the assumed
relationship between the individual effect and the regressors. Under
a fixed-effects specification, \eqn{\alpha_i} may be correlated with
the regressor history. Under a random-effects specification, the
individual effect is assumed to be mean independent of the regressors,
together with the usual variance-component assumptions.

The function \code{fast_tau_panel()} does not estimate a separate
model of its own. It is a convenience dispatcher:

\itemize{

  \item If \code{effect = "fe"}, the call is passed to
  \code{\link{fast_tau_fe}}. The default proposed fixed-effects
  estimator uses first differencing,
  \deqn{
  \Delta y_{it}
  =
  \Delta x_{it}^{\mathsf T}\beta
  +
  \Delta\varepsilon_{it},
  }
  which eliminates the time-invariant individual effect exactly.
  First differencing also limits propagation of a contaminated
  unit-time observation to at most two adjacent transformed
  observations.

  \item If \code{effect = "re"}, the call is passed to
  \code{\link{fast_tau_re}}. The proposed random-effects procedure
  uses a two-stage construction. A difference-based fast-Tau pilot
  is first obtained and used to construct robust variance-component
  estimates and residual-based weights. The original panel is then
  robustly quasi-demeaned and a second fast-Tau fit produces the
  final regression coefficient estimate.
}

Both procedures use the same underlying regression Tau-estimation
principle. For transformed residuals \eqn{e_j(\beta)}, the estimator
minimizes a Tau-scale of the general form
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
where \eqn{s(\beta)} is a robust M-scale based on a bounded loss
function \eqn{\rho_1}, while a second bounded loss function
\eqn{\rho_2} controls the efficiency component of the criterion.

The Tau objective is non-convex, so both the fixed- and random-effects
procedures employ a fast multi-start strategy. Cluster-diverse
starting subsamples, preliminary iteratively reweighted least-squares
(IRLS) updates, candidate screening, shortlisting, and final
refinement are used to obtain a high-quality numerical solution
without exhaustively searching all local minima.

The loss-function calibration used by the package combines a
high-breakdown M-scale with approximately 95 percent Gaussian
efficiency in the classical independent-error regression setting.
This percentage is a tuning calibration and should not be interpreted
as an exact efficiency statement for transformed panel data, where
within-unit dependence may be present.
}

\value{

If \code{effect = "fe"}, an object of classes
\code{"fast_tau_fe_fit"} and \code{"fast_tau_panel_fit"} is returned.
The object is the result of \code{\link{fast_tau_fe}}.

If \code{effect = "re"}, an object of classes
\code{"fast_tau_re_fit"} and \code{"fast_tau_panel_fit"} is returned.
The object is the result of \code{\link{fast_tau_re}}.

In both cases, the returned object contains at least the estimated
regression coefficient vector \code{beta}, together with fitting
information specific to the selected panel specification.
}
\examples{

## Generate a small fixed-effects panel
dat_fe <- generate_panel_data(
  N_ = 60,
  T_ = 4,
  K = 4,
  effect_model = "fe",
  contam_scheme = "none"
)

## Fit through the unified interface.
fit_fe <- fast_tau_panel(
  Y = dat_fe$Y,
  X = dat_fe$X,
  effect = "fe",
  N = 100,
  tt = 5
)

fit_fe$beta


## Generate a small random-effects panel
dat_re <- generate_panel_data(
  N_ = 60,
  T_ = 4,
  K = 4,
  effect_model = "re",
  contam_scheme = "none"
)

fit_re <- fast_tau_panel(
  Y = dat_re$Y,
  X = dat_re$X,
  effect = "re",
  N = 100,
  tt = 5
)

fit_re$beta


## Structural-mean predictions
pred_fe <- predict(fit_fe, newdata = dat_fe$X)
pred_re <- predict(fit_re, newdata = dat_re$X)
}
