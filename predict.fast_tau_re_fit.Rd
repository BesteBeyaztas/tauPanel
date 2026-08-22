\name{predict.fast_tau_re_fit}
\alias{predict.fast_tau_re_fit}

\title{Predict from a Fast Tau Random-Effects Model}

\description{
Computes marginal structural-mean predictions from a fitted fast Tau
random-effects panel model. Predictions are based on the estimated common
regression component, \eqn{x^{\mathsf T}\widehat\beta_{\tau,\mathrm{RE}}},
and do not include a predicted individual-specific random effect.
}

\usage{
\method{predict}{fast_tau_re_fit}(object, newdata, ...)
}

\arguments{

  \item{object}{
A fitted object of class \code{"fast_tau_re_fit"}, typically returned
by \code{\link{fast_tau_re}} or by \code{\link{fast_tau_panel}} with
\code{effect = "re"}.
}

  \item{newdata}{
New predictor values on the original regressor scale.

The argument may be supplied in one of three forms:

\itemize{
  \item a numeric vector of length \eqn{p}, representing one
  predictor vector;

  \item a numeric matrix or data frame with \eqn{p} columns,
  representing multiple predictor vectors;

  \item a numeric \eqn{N \times T \times p} array, representing
  panel-structured predictor values for \eqn{N} units and
  \eqn{T} time periods.
}

In every case, \eqn{p} must equal the number of regression
coefficients in the fitted model.
}

  \item{...}{
Additional arguments included for compatibility with the generic
\code{\link[stats]{predict}} interface. They are currently ignored.
}

}

\details{

Consider the random-effects panel regression model
\deqn{
y_{it}
=
x_{it}^{\mathsf T}\beta
+
\alpha_i
+
\varepsilon_{it},
}
where \eqn{\alpha_i} is the individual-specific random effect and
\eqn{\varepsilon_{it}} is the idiosyncratic error.

The proposed random-effects estimator yields an estimate
\eqn{\widehat\beta_{\tau,\mathrm{RE}}} of the common regression
coefficient vector. For a supplied predictor vector \eqn{x}, this
method returns
\deqn{
\widehat\mu(x)
=
x^{\mathsf T}\widehat\beta_{\tau,\mathrm{RE}}.
}

Thus, the method predicts the structural or marginal regression
component rather than an individual-specific conditional response.

The random-effects specification assumes that the individual effect
has zero mean conditional on the relevant regressor information.
For a previously unobserved unit with no response history, the
realized random effect cannot be estimated from the new regressors
alone. Its corresponding prediction is therefore based on
\deqn{
E(y_{\mathrm{new}}\mid x_{\mathrm{new}})
=
x_{\mathrm{new}}^{\mathsf T}\beta,
}
and the fitted counterpart is
\deqn{
\widehat y_{\mathrm{new}}
=
x_{\mathrm{new}}^{\mathsf T}
\widehat\beta_{\tau,\mathrm{RE}}.
}

For a unit already observed in the fitting sample, an
individual-specific conditional prediction can additionally use a
prediction of its random effect. In the accompanying methodology,
the random effect is predicted using the plug-in shrinkage form
\deqn{
\widehat\alpha_i
=
\widehat\lambda\,\bar r_i,
}
where
\deqn{
\bar r_i
=
\frac{1}{T}
\sum_{t=1}^{T}
\left(
y_{it}
-
x_{it}^{\mathsf T}\widehat\beta_{\tau,\mathrm{RE}}
\right)
}
and
\deqn{
\widehat\lambda
=
\frac{\widehat\sigma_\alpha^2}
{\widehat\sigma_\alpha^2+
 \widehat\sigma_\varepsilon^2/T}.
}
Equivalently,
\deqn{
\widehat\lambda
=
1-(1-\widehat\theta)^2,
}
where \eqn{\widehat\theta} is the robust quasi-demeaning parameter
returned by \code{\link{fast_tau_re}}.

A corresponding conditional fitted value for an observed unit is
therefore
\deqn{
\widehat y_{it}
=
x_{it}^{\mathsf T}
\widehat\beta_{\tau,\mathrm{RE}}
+
\widehat\alpha_i.
}

The present \code{predict()} method does not perform this
individual-effect prediction, because it receives only
\code{newdata} and the fitted model object and does not use the
response history required to construct \eqn{\bar r_i}. It therefore
returns the structural component only.

The predictor values supplied through \code{newdata} must be on the
original regressor scale. Users should not apply the robust
quasi-demeaning transformation before prediction. That transformation
is an estimation device used internally to obtain
\eqn{\widehat\beta_{\tau,\mathrm{RE}}}; predictions from new
regressors are evaluated on the original scale.

The shape of the returned result follows the shape of
\code{newdata}.
}

\value{

Predicted structural means based on
\eqn{x^{\mathsf T}\widehat\beta_{\tau,\mathrm{RE}}}. The form of
the returned object depends on \code{newdata}:

\itemize{

  \item If \code{newdata} is a numeric vector of length \eqn{p},
  a single numeric prediction is returned.

  \item If \code{newdata} is a numeric matrix or data frame with
  \eqn{m} rows and \eqn{p} columns, a numeric vector of length
  \eqn{m} is returned.

  \item If \code{newdata} is an
  \eqn{N \times T \times p} numeric array, an
  \eqn{N \times T} numeric matrix of predictions is returned.
}

For array input, available dimension names for the first two
dimensions are preserved in the returned prediction matrix.
}
\examples{

## Generate a random-effects panel

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

## Panel-structured structural predictions
pred_panel <- predict(
  fit,
  newdata = dat$X
)

dim(pred_panel)


## Prediction for one new predictor vector
x_new <- c(0.5, -0.2, 1.1, 0.7)

predict(
  fit,
  newdata = x_new
)


## Predictions for several new predictor vectors
X_new <- matrix(
  rnorm(20),
  nrow = 5,
  ncol = 4
)

predict(
  fit,
  newdata = X_new
)
}
