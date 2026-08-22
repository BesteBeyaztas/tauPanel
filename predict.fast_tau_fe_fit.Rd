\name{predict.fast_tau_fe_fit}
\alias{predict.fast_tau_fe_fit}

\title{Predict from a Fast Tau Fixed-Effects Model}

\description{
Computes structural-mean predictions from a fitted fast Tau
fixed-effects panel model. Predictions are obtained from the estimated
common regression component, \eqn{x^{\mathsf T}\widehat\beta_\tau},
and do not include an estimated individual fixed effect.
}

\usage{
\method{predict}{fast_tau_fe_fit}(object, newdata, ...)
}

\arguments{

  \item{object}{
A fitted object of class \code{"fast_tau_fe_fit"}, typically returned
by \code{\link{fast_tau_fe}} or by \code{\link{fast_tau_panel}} with
\code{effect = "fe"}.
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

In all cases, \eqn{p} must equal the number of coefficients in the
fitted model.
}

  \item{...}{
Additional arguments for compatibility with the generic
\code{\link[stats]{predict}} interface. They are currently ignored.
}

}

\details{

For the fixed-effects panel model
\deqn{
y_{it}
=
x_{it}^{\mathsf T}\beta
+
\alpha_i
+
\varepsilon_{it},
}
the fitted fast Tau estimator provides an estimate
\eqn{\widehat\beta_\tau} of the common regression coefficient vector.
For a supplied predictor vector \eqn{x}, this method returns the
structural component
\deqn{
\widehat\mu(x)
=
x^{\mathsf T}\widehat\beta_\tau.
}

The prediction is evaluated using the regressors on their
\emph{original scale}. Users should therefore supply the original
predictor values in \code{newdata}; first differences, within
transformations, or other fixed-effects transformations should not
be applied before calling \code{predict()}.

This distinction is important because the transformation used during
estimation serves to eliminate the individual-specific effect and
estimate the common slope vector. Once
\eqn{\widehat\beta_\tau} has been estimated, structural prediction
for new predictor values is naturally made on the original scale.

For an individual already observed in the fitting sample, a
conditional fitted value can in principle be written as
\deqn{
\widehat y_{it}
=
x_{it}^{\mathsf T}\widehat\beta_\tau
+
\widehat\alpha_i,
}
where, for the proposed fixed-effects procedure, the accompanying
methodology estimates the individual effect using the unitwise
residual median,
\deqn{
\widehat\alpha_i
=
\operatorname{median}_{t}
\left\{
y_{it}
-
x_{it}^{\mathsf T}\widehat\beta_\tau
\right\}.
}
Such a quantity requires observed response history for that
individual.

The present method is intentionally different: it uses only
\code{newdata} and the estimated coefficient vector stored in
\code{object}. Consequently, it returns
\eqn{x^{\mathsf T}\widehat\beta_\tau} only and does not estimate or
add an individual-specific fixed effect. This makes the method
appropriate for structural prediction and for predictor sets
belonging to previously unobserved units.

The shape of the returned object follows the shape of
\code{newdata}. No additional panel transformation is performed
during prediction.
}

\value{

Predicted structural means based on
\eqn{x^{\mathsf T}\widehat\beta_\tau}. The returned object depends
on the form of \code{newdata}:

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
dimensions are preserved in the prediction matrix.
}
\examples{

## Generate a fixed-effects panel

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

## Panel-structured predictions
pred_panel <- predict(
  fit,
  newdata = dat$X
)

dim(pred_panel)


## Prediction for one predictor vector
x_new <- c(0.5, -0.2, 1.1, 0.7)

predict(
  fit,
  newdata = x_new
)


## Predictions for several predictor vectors
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
