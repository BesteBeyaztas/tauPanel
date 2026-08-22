\name{fast_tau_re}
\alias{fast_tau_re}

\title{Fast Tau Estimation for Random-Effects Panel Models}

\description{
Fits the proposed robust fast Tau estimator for a balanced
random-effects panel regression model. The procedure combines a
first-difference fast-Tau pilot fit, robust scale-based variance-component
estimation, residual-based observation weights, robust quasi-demeaning,
and a final fast-Tau regression.
}

\usage{
fast_tau_re(
  Y,
  X,
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
where \eqn{n} is the number of cross-sectional units and \eqn{T}
is the number of time periods. The panel must be balanced.
}

  \item{X}{
A numeric \eqn{n \times T \times p} array containing the regressors,
where \eqn{p} is the number of regression coefficients. The first
two dimensions must agree with those of \code{Y}.
}

  \item{N}{
Number of randomized starting subsamples used by each application
of the fast-Tau algorithm. The default is \code{2000}. This quantity
corresponds to the number of random starts, denoted by \eqn{S} in
the methodological description, and should not be confused with the
number of cross-sectional units.
}

  \item{kk}{
Number of preliminary iteratively reweighted least-squares (IRLS)
updates applied to each random or diagnostic starting candidate.
The default is \code{2}, corresponding to \eqn{\kappa = 2}.
}

  \item{tt}{
Number of the most promising candidate solutions retained for
subsequent refinement. The default is \code{20}, corresponding to
the shortlist size \eqn{\ell = 20}.
}

  \item{rr}{
Number of fixed-point scale updates used when an approximate
residual M-scale calculation is requested. The default is \code{2}.
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
at most one transformed observation from any selected
cross-sectional unit enters a starting subsample.
}

  \item{m_cluster}{
Number of distinct cross-sectional units used to form each
cluster-diverse starting subsample. If \code{NULL}, the number
of selected units is set equal to the number of regression
coefficients, yielding a cluster-diverse elemental start.
Larger values give overdetermined cluster-diverse starts,
subject to the number of available units.
}

  \item{seed}{
Optional integer seed used for randomized starting-subsample
generation. The default is \code{NULL}. The seed is handled
locally within the fast-Tau computation.
}

}

\details{

Consider the balanced random-effects panel regression model
\deqn{
y_{it}
=
x_{it}^{\mathsf T}\beta
+
\alpha_i
+
\varepsilon_{it},
\qquad
i=1,\ldots,n,\quad t=1,\ldots,T,
}
where \eqn{\beta} is the common regression coefficient vector,
\eqn{\alpha_i} is an individual-specific random effect, and
\eqn{\varepsilon_{it}} is the idiosyncratic error.

The random-effects specification requires the individual effects to
satisfy the corresponding mean-independence condition with respect
to the regressor history. Under the classical variance-components
model,
\deqn{
\operatorname{Var}(\alpha_i)=\sigma_\alpha^2,
\qquad
\operatorname{Var}(\varepsilon_{it})=\sigma_\varepsilon^2,
}
the conventional quasi-demeaning parameter is
\deqn{
\theta
=
1-
\left\{
\frac{\sigma_\varepsilon^2}
{\sigma_\varepsilon^2+T\sigma_\alpha^2}
\right\}^{1/2}.
}
The proposed estimator replaces the classical preliminary
variance-component and unit-centering steps by robust counterparts.

The procedure consists of two fast-Tau regression stages connected
by robust variance-component estimation and robust quasi-demeaning.

\itemize{

  \item \strong{Pilot fast-Tau fit.}
  The individual effects are first eliminated by first differencing,
  \deqn{
  \Delta y_{it}
  =
  \Delta x_{it}^{\mathsf T}\beta
  +
  \Delta\varepsilon_{it}.
  }
  The fast-Tau estimator is applied to this first-differenced sample,
  producing a robust pilot coefficient vector
  \eqn{\widehat\beta_\tau^{(0)}}.

  \item \strong{Robust scale-based variance components.}
  Pilot residuals are reconstructed on the original level scale,
  \deqn{
  r_{it}^{(0)}
  =
  y_{it}-x_{it}^{\mathsf T}\widehat\beta_\tau^{(0)}.
  }
  Pairwise residual differences and sums are formed over distinct
  time points. Let \eqn{q_D} and \eqn{q_S} denote normal-calibrated
  \eqn{Q_n} robust scales of the pooled residual differences and
  residual sums, respectively. The robust component estimates are
  \deqn{
  \widehat\sigma_\varepsilon^2
  =
  \max\left\{
  \frac{q_D^2}{2},\delta
  \right\},
  }
  and
  \deqn{
  \widehat\sigma_\alpha^2
  =
  \max\left\{
  \frac{q_S^2-q_D^2}{4},\delta
  \right\},
  }
  where a small positive lower bound is used for numerical
  stability.

  \item \strong{Robust quasi-demeaning.}
  The estimated quasi-demeaning parameter is
  \deqn{
  \widehat\theta
  =
  1-
  \left\{
  \frac{\widehat\sigma_\varepsilon^2}
  {\widehat\sigma_\varepsilon^2+
  T\widehat\sigma_\alpha^2}
  \right\}^{1/2},
  }
  restricted numerically to the interval \eqn{[0,1]}.

  To reduce the sensitivity of unit-specific centering to isolated
  contamination, the pilot residuals are centered within each unit
  and converted to bounded, redescending Tukey-bisquare weights.
  Using weights \eqn{\varpi_{it}}, weighted unit means are formed as
  \deqn{
  \bar y_{i,\varpi}
  =
  \frac{\sum_t \varpi_{it}y_{it}}
       {\sum_t \varpi_{it}},
  \qquad
  \bar x_{i,\varpi}
  =
  \frac{\sum_t \varpi_{it}x_{it}}
       {\sum_t \varpi_{it}}.
  }
  The robustly quasi-demeaned observations are then
  \deqn{
  y_{it}^{\dagger}
  =
  y_{it}
  -
  \widehat\theta\bar y_{i,\varpi},
  \qquad
  x_{it}^{\dagger}
  =
  x_{it}
  -
  \widehat\theta\bar x_{i,\varpi}.
  }

  \item \strong{Final fast-Tau fit.}
  The fast-Tau estimator is applied a second time to
  \eqn{(y_{it}^{\dagger},x_{it}^{\dagger})}. The resulting
  coefficient vector is returned as the final random-effects
  fast-Tau estimate.
}

The same multi-start fast-Tau settings are used in both the pilot
and final regression stages. The Tau estimator minimizes a robust
Tau-scale constructed from an M-scale and a second bounded loss
function. This permits the robustness and efficiency components of
the regression criterion to be calibrated separately.

Because the Tau objective is non-convex, the numerical algorithm
uses multiple cluster-diverse starting values, a small number of
preliminary IRLS refinements, screening and shortlisting, followed
by more extensive refinement of promising candidates. The default
configuration uses 2000 random starts, two preliminary IRLS updates,
and a shortlist of 20 candidates.

The regressor matrix supplied to each fast-Tau stage is internally
rescaled columnwise to improve numerical conditioning. Coefficients
are transformed back to the original regressor scale before being
returned.

Under the Gaussian random-effects model, the normal calibration of
the \eqn{Q_n} scales makes the quantities
\eqn{\widehat\sigma_\varepsilon^2} and
\eqn{\widehat\sigma_\alpha^2} target the corresponding classical
variance components. Under non-Gaussian distributions, squared
robust scales need not equal second moments. In that setting, these
quantities should be interpreted as robust, normal-calibrated
scale-based analogues of the variance components rather than
necessarily as consistent estimators of the population variances.
}

\value{

An object of classes \code{"fast_tau_re_fit"} and
\code{"fast_tau_panel_fit"}. It is a list containing:

\item{beta}{
Numeric vector containing the final fast-Tau estimates of the
regression coefficients on the original regressor scale.
}

\item{sigma2_a}{
Robust scale-based estimate of the individual-effect component
\eqn{\sigma_\alpha^2}. Under the Gaussian RE model this targets
the classical individual-effect variance component.
}

\item{sigma2_e}{
Robust scale-based estimate of the idiosyncratic component
\eqn{\sigma_\varepsilon^2}. Under the Gaussian RE model this
targets the classical idiosyncratic-error variance component.
}

\item{theta}{
Estimated robust quasi-demeaning parameter,
\deqn{
\widehat\theta
=
1-
\left\{
\frac{\widehat\sigma_\varepsilon^2}
{\widehat\sigma_\varepsilon^2+
T\widehat\sigma_\alpha^2}
\right\}^{1/2}.
}
}

\item{fit}{
A list containing additional fitting information:

\itemize{
  \item \code{final}: the final fast-Tau fit obtained from the
  robustly quasi-demeaned panel;

  \item \code{preliminary}: the first-difference fast-Tau pilot fit;

  \item \code{beta_initial}: the coefficient vector obtained from
  the pilot stage;

  \item \code{transformed_data}: the robustly quasi-demeaned data
  supplied to the final fast-Tau fit. This object includes the
  transformed response and regressors, unit identifiers,
  observation-specific weights, weighted unit means, the robust
  residual scale, and the quasi-demeaning parameter.
}
}

\item{effect}{
Character string equal to \code{"random"}, identifying the fitted
model as a random-effects specification.
}

}
\examples{

## Generate a small random-effects panel
dat <- generate_panel_data(
  N_ = 60,
  T_ = 4,
  K = 4,
  effect_model = "re",
  contam_scheme = "none"
)

## Fit the proposed random-effects estimator.
fit <- fast_tau_re(
  Y = dat$Y,
  X = dat$X,
  N = 100,
  tt = 5
)

## Final regression coefficients
fit$beta

## Robust scale-based variance components
fit$sigma2_a
fit$sigma2_e

## Robust quasi-demeaning parameter
fit$theta

## Structural-mean predictions
pred <- predict(
  fit,
  newdata = dat$X
)

dim(pred)
}
