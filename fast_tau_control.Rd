\name{fast_tau_control}
\alias{fast_tau_control}

\title{Default Control Parameters for Fast Tau Panel Estimation}

\description{
Returns the default computational and transformation settings used by the
fast Tau estimators for fixed- and random-effects panel data models.
The defaults correspond to the numerical configuration adopted for the
proposed estimators.
}

\usage{
fast_tau_control()
}

\details{
The proposed panel estimator is based on regression Tau-estimation, in which
the coefficient vector is obtained by minimizing a Tau-scale of the transformed
regression residuals. The Tau-scale combines a robust M-scale with a second
bounded loss function, thereby allowing robustness and efficiency to be
controlled separately.

Because the Tau objective is non-convex and may contain multiple local minima,
the estimator is computed using a multi-start fast-Tau strategy. Random
starting subsets are partially refined by iteratively reweighted least squares
(IRLS), promising candidates are retained in a shortlist, and only this
shortlist is subjected to more extensive refinement.

For panel data, the default implementation uses cluster-diverse starting
subsamples. Distinct cross-sectional units are selected and at most one
transformed observation is drawn from each selected unit. This construction
reduces the possibility that a starting subsample is dominated by several
dependent observations from the same individual. When \code{m_cluster = NULL},
the number of selected units is set equal to the number of regression
coefficients, producing a cluster-diverse elemental start.

The default computational settings correspond to 2000 random starts, two
preliminary IRLS updates per candidate, and a shortlist containing the 20 most
promising candidates. In the notation of the fast-Tau algorithm, these
quantities correspond respectively to \eqn{S = 2000}, \eqn{\kappa = 2}, and
\eqn{\ell = 20}. Thus, the component \code{N} returned by this function denotes
the number of random starts and should not be confused with the conventional
use of \eqn{N} for the number of cross-sectional units in a panel.

For fixed-effects estimation, the default transformation is first
differencing. If
\deqn{y_{it} = x_{it}^{\mathsf T}\beta + \alpha_i + \varepsilon_{it},}
then first differencing gives
\deqn{
\Delta y_{it}
=
\Delta x_{it}^{\mathsf T}\beta
+
\Delta\varepsilon_{it},
}
and therefore eliminates the time-invariant individual effect
\eqn{\alpha_i}. In addition to this standard identification property, first
differencing has a robustness advantage in the proposed procedure: contamination
of a single unit-time observation affects at most two adjacent transformed
observations, whereas ordinary within demeaning can propagate its effect
through the entire transformed block of that unit.

The argument \code{approximate} controls how the residual M-scale is evaluated
during candidate comparison. With the default \code{approximate = 0}, the
algorithm uses the full M-scale calculation implemented by the fitting routine.
When a nonzero value is supplied, a shorter scale approximation is used, with
\code{rr} determining the number of scale-update iterations.
}

\value{
A named list containing the following components:

\item{N}{
Number of randomized starting subsamples used by the fast-Tau algorithm.
The default is \code{2000}. This corresponds to \eqn{S} in the algorithmic
notation of the accompanying methodology, and not to the number of
cross-sectional units.
}

\item{kk}{
Number of preliminary IRLS updates applied to each random or diagnostic
starting candidate before candidate screening and shortlisting.
The default is \code{2}, corresponding to \eqn{\kappa = 2}.
}

\item{tt}{
Number of best candidate solutions retained for subsequent refinement.
The default is \code{20}, corresponding to a shortlist size
\eqn{\ell = 20}.
}

\item{rr}{
Number of fixed-point scale-update iterations used when an approximate
M-scale calculation is requested. The default is \code{2}. This component
is relevant when \code{approximate} is nonzero.
}

\item{approximate}{
Indicator controlling the calculation of the residual M-scale during the
candidate-evaluation stage. The default \code{0} uses the full scale
calculation implemented by the fitting routine. A nonzero value requests
the shorter approximation based on \code{rr} scale updates.
}

\item{cluster_start}{
Logical value indicating whether cluster-diverse starting subsamples are
used. The default is \code{TRUE}. Under cluster-diverse sampling, each
starting subsample contains at most one transformed observation from any
selected cross-sectional unit.
}

\item{m_cluster}{
Number of distinct cross-sectional units used to construct each
cluster-diverse starting subsample. The default is \code{NULL}, in which
case the fitting routine uses the number of regression coefficients,
producing an elemental cluster-diverse starting subset. Values larger than
the number of coefficients produce overdetermined cluster-diverse starts,
subject to the number of available units.
}

\item{seed}{
Optional random-number seed used for the randomized starting subsamples.
The default is \code{NULL}, so no fixed local seed is imposed.
}

\item{fe_transform}{
Character string specifying the default transformation for the
fixed-effects estimator. The default is \code{"first_difference"}.
}
}
\examples{
## Obtain the default control settings
control <- fast_tau_control()
control

## Inspect individual settings
control$N
control$kk
control$tt
control$fe_transform

## Control values are returned as an ordinary list and may be modified
## before being supplied to the corresponding fitting function.
control$N <- 1000
control$tt <- 10
}
