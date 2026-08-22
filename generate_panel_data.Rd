\name{generate_panel_data}
\alias{generate_panel_data}

\title{Generate Balanced Panel Data with Optional Contamination}

\description{
Generates balanced linear panel data under either a random-effects or
correlated fixed-effects specification. The function supports Gaussian,
Student-\eqn{t}, and chi-squared idiosyncratic errors together with several
cellwise and blockwise vertical and leverage contamination mechanisms.
}

\usage{
generate_panel_data(
  N_ = 100,
  T_ = 3,
  K = 4,
  distribution = "normal",
  beta = c(-3.4, 1.5, 2.8, -1.2),
  sigma2_a = 1,
  sigma2_e = 1,
  effect_model = c("re", "fe"),
  alpha_corr = 0.6,
  gamma_alpha = NULL,
  contam_scheme = c(
    "none",
    "random_vertical",
    "block_vertical",
    "random_leverage",
    "block_leverage",
    "random_vertical_leverage"
  ),
  eps_mean = 10,
  eps_std = 1,
  cont = 0.05,
  mu_leverage = 5,
  std_leverage = 1,
  block_fraction = 0.5,
  chi_df = NULL
)
}

\arguments{

  \item{N_}{
Number of cross-sectional units. The default is \code{100}.
}

  \item{T_}{
Number of time periods observed for each cross-sectional unit.
The default is \code{3}. The generated panel is balanced.
}

  \item{K}{
Number of regressors. The default is \code{4}. The length of
\code{beta} must equal \code{K}.
}

  \item{distribution}{
Character string specifying the idiosyncratic-error distribution.
Available choices are \code{"normal"}, \code{"t"}, and \code{"chi"}.

For \code{"normal"}, errors are generated from a centered Gaussian
distribution with variance \code{sigma2_e}. For \code{"t"}, errors
are generated from a Student-\eqn{t} distribution with 5 degrees of
freedom. For \code{"chi"}, errors are generated from an uncentered
chi-squared distribution with \code{chi_df} degrees of freedom.
}

  \item{beta}{
Numeric vector of regression coefficients. Its length must equal
\code{K}. The default is
\code{c(-3.4, 1.5, 2.8, -1.2)}.
}

  \item{sigma2_a}{
Variance parameter controlling the magnitude of the individual-specific
effects. The default is \code{1}.
}

  \item{sigma2_e}{
Idiosyncratic-error variance used when
\code{distribution = "normal"}. The default is \code{1}.
}

  \item{effect_model}{
Character string specifying the individual-effect design.
Use \code{"re"} for random effects or \code{"fe"} for correlated
fixed effects. The default is \code{"re"}.
}

  \item{alpha_corr}{
Target correlation parameter used in the fixed-effects
Mundlak-style construction. The default is \code{0.6}.
Internally, the value is restricted to the interval
\eqn{[-0.99,0.99]}. This argument is ignored when
\code{effect_model = "re"}.
}

  \item{gamma_alpha}{
Optional numeric vector of length \code{K} defining the direction
of the Mundlak index used to generate correlated fixed effects.
The vector is normalized internally to have Euclidean norm one.

If \code{NULL} and \code{K = 4}, the default direction is
\code{c(1.0, -0.7, 0.5, 0.3)}. For other values of \code{K},
a decreasing sequence from \code{1} to \code{0.4} is used.
This argument is relevant only when \code{effect_model = "fe"}.
}

  \item{contam_scheme}{
Character string specifying the contamination mechanism. Available
choices are:

\code{"none"},
\code{"random_vertical"},
\code{"block_vertical"},
\code{"random_leverage"},
\code{"block_leverage"}, and
\code{"random_vertical_leverage"}.

See \code{Details} for the precise construction.
}

  \item{eps_mean}{
Mean of the Gaussian distribution used to replace idiosyncratic
errors under vertical contamination. The default is \code{10}.
}

  \item{eps_std}{
Standard deviation of the Gaussian distribution used to replace
idiosyncratic errors under vertical contamination. The default is
\code{1}.
}

  \item{cont}{
Target fraction of the \code{N_ T_} unit-time observations to
contaminate. The default is \code{0.05}. For random cellwise
contamination, approximately
\code{round(cont * N_ * T_)} observations are selected without
replacement. Under block contamination, the achieved fraction may
differ slightly because complete time blocks are contaminated.
}

  \item{mu_leverage}{
Mean of each regressor component used to generate leverage
contamination. The default is \code{5}.
}

  \item{std_leverage}{
Standard deviation of each regressor component used to generate
leverage contamination. The default is \code{1}.
}

  \item{block_fraction}{
Fraction of the time periods included in the consecutive
contamination block for each selected unit under
\code{"block_vertical"} and \code{"block_leverage"}.
The block length is
\eqn{\max\{1,\lceil f T\rceil\}} where \eqn{f} is specified by \code{block_fraction}.
The default is \code{0.5}.
}

  \item{chi_df}{
Degrees of freedom of the chi-squared idiosyncratic-error
distribution when \code{distribution = "chi"}.

If \code{NULL}, the function uses \code{2} for the fixed-effects
design and \code{1} for the random-effects design. A user-supplied
value must be a positive finite scalar.
}

}

\details{

The generated data follow the balanced static panel model
\deqn{
y_{it}
=
x_{it}^{\mathsf T}\beta
+
\alpha_i
+
\varepsilon_{it},
\qquad
i=1,\ldots,N,\quad
t=1,\ldots,T,
}
where \eqn{x_{it}} is a \eqn{K}-dimensional regressor vector,
\eqn{\beta} is the common slope vector, \eqn{\alpha_i} is the
individual-specific effect, and \eqn{\varepsilon_{it}} is the
idiosyncratic error.

\strong{Regressor generation}

The regressor array has dimensions
\eqn{N \times T \times K}. Regressors are initially generated
independently from standard normal distributions, except for the
first regressor, which is replaced by
\deqn{
x_{it1}=w_{it}-2,
\qquad
w_{it}\sim\chi^2_2.
}
Thus, the first regressor is right-skewed but centered at zero,
whereas the remaining regressors are standard Gaussian.

No intercept column is included in the generated design.

\strong{Random-effects design}

When \code{effect_model = "re"}, the individual effects are generated
independently of the regressors as
\deqn{
\alpha_i\sim N(0,\sigma_\alpha^2),
}
where \code{sigma2_a} specifies \eqn{\sigma_\alpha^2}. This corresponds
to the standard random-effects requirement that the individual effect
be independent of, and hence mean independent of, the regressor
history in the data-generating mechanism.

\strong{Fixed-effects design}

When \code{effect_model = "fe"}, the individual effects are explicitly
correlated with the regressor history through a Mundlak-style
construction.

Let
\deqn{
\bar{x}_i
=
\frac{1}{T}\sum_{t=1}^{T} x_{it}
}
denote the vector of unit-specific regressor means. After normalizing
\code{gamma_alpha} to unit Euclidean norm, define the standardized
Mundlak index
\deqn{
z_i
=
\operatorname{standardize}
\left(
\gamma_\alpha^{\mathsf T}\bar x_i
\right).
}

An independent Gaussian vector is generated, orthogonalized with
respect to \eqn{z}, and standardized to produce \eqn{\eta_i}.
The fixed effect is then constructed as
\deqn{
\alpha_i
=
\sigma_\alpha
\left\{
\rho z_i
+
\sqrt{1-\rho^2}\eta_i
\right\},
}
where
\eqn{\rho} is the value of \code{alpha_corr} restricted to
\eqn{[-0.99,0.99]}.

The orthogonalization makes the realized dependence between the
individual effects and the Mundlak index closely controlled by
\code{alpha_corr}, while \code{sigma2_a} determines the scale of
the individual effects. This construction is used to generate an
FE setting in which the individual effects are intentionally
correlated with the regressors.

\strong{Idiosyncratic errors}

Before contamination, the idiosyncratic errors are generated
according to \code{distribution}:

\itemize{
  \item \code{"normal"}:
  \deqn{
\varepsilon_{it}\sim N(0,\sigma_\varepsilon^2),
}
where \eqn{\sigma_\varepsilon^2} is specified by
\code{sigma2_e}.

  \item \code{"t"}:
  \deqn{
  \varepsilon_{it}\sim t_5.
  }
  The current implementation does not rescale the Student-\eqn{t}
  draws using \code{sigma2_e}.

  \item \code{"chi"}:
\deqn{
\varepsilon_{it}\sim\chi^2_{\nu}.
}
The degrees of freedom \eqn{\nu} are specified by
\code{chi_df}. The chi-squared errors are deliberately
uncentered. The current implementation does not rescale these
draws using \code{sigma2_e}.
}

\strong{Contamination schemes}

Contamination is applied after the clean regressors and individual
effects have been generated.

\itemize{

  \item \code{"none"}:
  no observations are contaminated.

  \item \code{"random_vertical"}:
  approximately a fraction \code{cont} of unit-time cells is sampled
  uniformly without replacement. At each selected cell, the original
  idiosyncratic error is replaced by an independent draw from a
  Gaussian distribution with mean \code{eps_mean} and standard deviation
  \code{eps_std}.

  \item \code{"block_vertical"}:
  selected cross-sectional units receive vertical contamination over
  a randomly positioned consecutive block of time periods. The block
  length is determined by \code{block_fraction}. Errors in those cells
  are replaced by draws from a Gaussian distribution with mean \code{eps_mean}
  and standard deviation \code{eps_std}.

  \item \code{"random_leverage"}:
  approximately a fraction \code{cont} of unit-time cells is sampled
  uniformly without replacement. At every selected cell, all
  \eqn{K} observed regressor components are replaced independently by
  draws from independent Gaussian draws with mean \code{mu_leverage}
  and standard deviation \code{std_leverage}.

  \item \code{"block_leverage"}:
  leverage contamination is introduced over randomly positioned
  consecutive blocks within selected cross-sectional units. All
  regressor components in a contaminated cell are replaced by draws
  from independent Gaussian draws with mean \code{mu_leverage}
  and standard deviation \code{std_leverage}.

  \item \code{"random_vertical_leverage"}:
  the same randomly selected unit-time cells receive both forms of
  contamination. Their idiosyncratic errors are replaced by draws from
  \eqn{N(\mathrm{eps\_mean}, \mathrm{eps\_std}^2)}, and their observed
  regressor vectors are replaced by Gaussian leverage draws.
}

For all schemes involving leverage contamination, the response is
generated from the \emph{uncontaminated} regressor vector rather than
from the contaminated value subsequently observed by the estimator.
Consequently, the leverage and combined contamination schemes create
bad leverage observations: the contaminated carrier does not remain
on the regression surface.

The structural component returned as \code{MU} is therefore
\deqn{
\mu_{it}
=
x_{it,\mathrm{clean}}^{\mathsf T}\beta
}
for leverage-contaminated designs. Under pure vertical contamination
or no contamination, the observed and response-generating regressors
coincide.

Finally, responses are generated as
\deqn{
y_{it}
=
\mu_{it}
+
\alpha_i
+
\varepsilon_{it}.
}
}

\value{

A list containing:

\item{X}{
Numeric array of dimensions
\eqn{N \times T \times K} containing the regressors observed by
the estimator. Under leverage contamination, selected entries differ
from the clean regressors used to generate the corresponding responses.
}

\item{Y}{
Numeric \eqn{N \times T} matrix containing the generated panel
responses.
}

\item{MU}{
Numeric \eqn{N \times T} matrix containing the structural regression
component \eqn{x_{it}^{\mathsf T}\beta} used to generate the response.
For leverage-contaminated cells, this is calculated from the
uncontaminated regressor vector.
}

\item{alpha}{
Numeric \eqn{N \times 1} matrix containing the generated
individual-specific effects.
}

\item{epsilon}{
Numeric \eqn{N \times T} matrix containing the final idiosyncratic
errors after any vertical contamination has been applied.
}

\item{outlier_cells}{
Integer vector identifying contaminated unit-time cells. Indices are
stored in unit-major order, so that cell \eqn{(i,t)} corresponds to
\code{(i - 1) * T_ + t}. An empty integer vector is returned when no
contamination is applied.
}

\item{outlier_units}{
Sorted integer vector identifying cross-sectional units containing
at least one contaminated cell. An empty integer vector is returned
when no contamination is applied.
}

\item{gamma_alpha}{
Normalized coefficient vector defining the Mundlak index in the
fixed-effects design. For the random-effects design this is typically
\code{NULL}.
}

\item{effect_model}{
Character string identifying the generated individual-effect
specification: \code{"fe"} or \code{"re"}.
}

\item{chi_df}{
Degrees of freedom used for chi-squared idiosyncratic errors.
The value is returned even when a different error distribution is used.
}

}
\examples{

## Clean random-effects panel
dat_re <- generate_panel_data(
  N_ = 50,
  T_ = 4,
  K = 4,
  effect_model = "re",
  contam_scheme = "none"
)

dim(dat_re$X)
dim(dat_re$Y)
dat_re$effect_model


## Correlated fixed-effects panel
dat_fe <- generate_panel_data(
  N_ = 50,
  T_ = 4,
  K = 4,
  effect_model = "fe",
  alpha_corr = 0.6,
  contam_scheme = "none"
)

dat_fe$gamma_alpha


## Random vertical contamination
dat_vertical <- generate_panel_data(
  N_ = 50,
  T_ = 5,
  K = 4,
  effect_model = "fe",
  contam_scheme = "random_vertical",
  cont = 0.10
)

length(dat_vertical$outlier_cells)
dat_vertical$outlier_units


## Bad-leverage contamination
dat_leverage <- generate_panel_data(
  N_ = 50,
  T_ = 5,
  K = 4,
  effect_model = "re",
  contam_scheme = "random_leverage",
  cont = 0.10
)

length(dat_leverage$outlier_cells)
}
