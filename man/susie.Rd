% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susie_ctwas_all.R
\name{susie}
\alias{susie}
\alias{susie_suff_stat}
\title{Sum of Single Effects (SuSiE) Regression}
\usage{
susie(
  X,
  Y,
  L = min(10, ncol(X)),
  scaled_prior_variance = 0.2,
  residual_variance = NULL,
  prior_weights = NULL,
  null_weight = NULL,
  standardize = TRUE,
  intercept = TRUE,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "simple"),
  check_null_threshold = 0,
  prior_tol = 1e-09,
  residual_variance_upperbound = Inf,
  s_init = NULL,
  coverage = 0.95,
  min_abs_corr = 0.5,
  compute_univariate_zscore = FALSE,
  na.rm = FALSE,
  max_iter = 100,
  tol = 0.001,
  verbose = FALSE,
  track_fit = FALSE,
  residual_variance_lowerbound = var(drop(Y))/10000,
  refine = FALSE
)

susie_suff_stat(
  bhat,
  shat,
  R,
  n,
  var_y,
  XtX,
  Xty,
  yty,
  maf = NULL,
  maf_thresh = 0,
  L = 10,
  scaled_prior_variance = 0.2,
  residual_variance = NULL,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "simple"),
  check_null_threshold = 0,
  prior_tol = 1e-09,
  r_tol = 1e-08,
  prior_weights = NULL,
  null_weight = NULL,
  standardize = TRUE,
  max_iter = 100,
  s_init = NULL,
  intercept_value = 0,
  coverage = 0.95,
  min_abs_corr = 0.5,
  tol = 0.001,
  verbose = FALSE,
  track_fit = FALSE,
  check_input = FALSE,
  refine = FALSE
)
}
\arguments{
\item{X}{An n by p matrix of covariates.}

\item{Y}{The observed responses, a vector of length n.}

\item{L}{Number of components (nonzero coefficients) in the susie
regression model. If L is larger than the number of covariates, p,
L is set to p.}

\item{scaled_prior_variance}{The scaled prior variance. This is
either a scalar or a vector of length \code{L}. The prior variance
of each non-zero element of b is set to \code{var(Y) *
  scaled_prior_variance}. If \code{estimate_prior_variance = TRUE},
this provides initial estimates of the prior variances.}

\item{residual_variance}{Variance of the residual. If
\code{estimate_residual_variance = TRUE}, this value provides the
initial estimate of the residual variance. By default, it is
\code{var(Y)}.}

\item{prior_weights}{A vector of length p, in which each entry
gives the prior probability that corresponding column of X has a
nonzero effect on the outcome, Y.}

\item{null_weight}{Prior probability of no effect (a number between
0 and 1, and cannot be exactly 1).}

\item{standardize}{If \code{standardize = TRUE}, standardize the
columns of X (or XtX and Xty) to unit variance prior to
fitting. Note that \code{scaled_prior_variance} specifies the prior
on the coefficients of X \emph{after} standardization (if it is
performed). If you do not standardize, you may need to think more
carefully about specifying \code{scaled_prior_variance}. Whatever
your choice, the coefficients returned by \code{coef} are given for
\code{X} on the original input scale. Any column of \code{X} that
has zero variance is not standardized.}

\item{intercept}{If \code{intercept = TRUE}, the intercept is
fitted; it \code{intercept = FALSE}, the intercept is set to
zero. Setting \code{intercept = FALSE} is generally not
recommended.}

\item{estimate_residual_variance}{If
\code{estimate_residual_variance = TRUE}, the residual variance is
estimated, using \code{residual_variance} as an initial value. If
\code{estimate_residual_variance = FALSE}, the residual variance is
fixed to the value supplied by \code{residual_variance}.}

\item{estimate_prior_variance}{If \code{estimate_prior_variance =
  TRUE}, the prior variance is estimated (this is a separate
parameter for each of the L effects). If provided,
\code{scaled_prior_variance} is then used as an initial value for
the optimization. When \code{estimate_prior_variance = FALSE}, the
prior variance for each of the L effects is determined by the
value supplied to \code{scaled_prior_variance}.}

\item{estimate_prior_method}{The method used for estimating prior
variance. When \code{estimate_prior_method = "simple"} is used, the
likelihood at the specified prior variance is compared to the
likelihood at a variance of zero, and the setting with the larger
likelihood is retained.}

\item{check_null_threshold}{When the prior variance is estimated,
compare the estimate with the null, and set the prior variance to
zero unless the log-likelihood using the estimate is larger by this
threshold amount. For example, if you set
\code{check_null_threshold = 0.1}, this will "nudge" the estimate
towards zero when the difference in log-likelihoods is small. A
note of caution that setting this to a value greater than zero may
lead the IBSS fitting procedure to occasionally decrease the ELBO.}

\item{prior_tol}{When the prior variance is estimated, compare the
estimated value to \code{prior_tol} at the end of the computation,
and exclude a single effect from PIP computation if the estimated
prior variance is smaller than this tolerance value.}

\item{residual_variance_upperbound}{Upper limit on the estimated
residual variance. It is only relevant when
\code{estimate_residual_variance = TRUE}.}

\item{s_init}{A previous susie fit with which to initialize.}

\item{coverage}{A number between 0 and 1 specifying the
\dQuote{coverage} of the estimated confidence sets.}

\item{min_abs_corr}{Minimum absolute correlation allowed in a
credible set. The default, 0.5, corresponds to a squared
correlation of 0.25, which is a commonly used threshold for
genotype data in genetic studies.}

\item{compute_univariate_zscore}{If \code{compute_univariate_zscore
  = TRUE}, the univariate regression z-scores are outputted for each
variable.}

\item{na.rm}{Drop any missing values in Y from both X and Y.}

\item{max_iter}{Maximum number of IBSS iterations to perform.}

\item{tol}{A small, non-negative number specifying the convergence
tolerance for the IBSS fitting procedure. The fitting procedure
will halt when the difference in the variational lower bound, or
\dQuote{ELBO} (the objective function to be maximized), is
less than \code{tol}.}

\item{verbose}{If \code{verbose = TRUE}, the algorithm's progress,
and a summary of the optimization settings, are printed to the
console.}

\item{track_fit}{If \code{track_fit = TRUE}, \code{trace}
is also returned containing detailed information about the
estimates at each iteration of the IBSS fitting procedure.}

\item{residual_variance_lowerbound}{Lower limit on the estimated
residual variance. It is only relevant when
\code{estimate_residual_variance = TRUE}.}

\item{refine}{If \code{refine = TRUE}, we use a procedure to help
SuSiE get out of local optimum.}

\item{bhat}{A p-vector of estimated effects.}

\item{shat}{A p-vector of standard errors.}

\item{R}{A p by p symmetric, positive semidefinite matrix. It
can be \eqn{X'X}, the covariance matrix \eqn{X'X/(n-1)}, or a
correlation matrix. It should be estimated from the same samples
used to compute \code{bhat} and \code{shat}. Using an out-of-sample
matrix may produce unreliable results.}

\item{n}{The sample size.}

\item{var_y}{The sample variance of y, defined as \eqn{y'y/(n-1)}.
When the sample variance cannot be provided, the coefficients
(returned from \code{coef}) are computed on the "standardized" X, y
scale.}

\item{XtX}{A p by p matrix \eqn{X'X} in which the columns of X
are centered to have mean zero.}

\item{Xty}{A p-vector \eqn{X'y} in which y and the columns of X are
centered to have mean zero.}

\item{yty}{A scalar \eqn{y'y} in which y is centered to have mean
zero.}

\item{maf}{Minor allele frequency; to be used along with
\code{maf_thresh} to filter input summary statistics.}

\item{maf_thresh}{Variants having a minor allele frequency smaller
than this threshold are not used.}

\item{r_tol}{Tolerance level for eigenvalue check of positive
semidefinite matrix of R.}

\item{intercept_value}{The intercept. (The intercept cannot be
estimated from centered summary data.) This setting will be used by
\code{coef} to assign an intercept value, mainly for consistency
with \code{susie}. Set to \code{NULL} if you want \code{coef} not
to include an intercept term (and so only a p-vector is returned).}

\item{check_input}{If \code{check_input = TRUE},
\code{susie_suff_stat} performs additional checks on \code{XtX} and
\code{Xty}. The checks are: (1) check that \code{XtX} is positive
semidefinite; (2) check that \code{Xty} is in the space spanned by
the non-zero eigenvectors of \code{XtX}.

@keywords internal}
}
\value{
A \code{"susie"} object with some or all of the following
elements:

\item{alpha}{An L by p matrix of posterior inclusion probabilites.}

\item{mu}{An L by p matrix of posterior means, conditional on
inclusion.}

\item{mu2}{An L by p matrix of posterior second moments,
conditional on inclusion.}

\item{Xr}{A vector of length n, equal to \code{X \%*\% colSums(alpha
  * mu)}.}

\item{lbf}{log-Bayes Factor for each single effect.}

\item{lbf_variable}{log-Bayes Factor for each variable and single effect.}

\item{intercept}{Intercept (fixed or estimated).}

\item{sigma2}{Residual variance (fixed or estimated).}

\item{V}{Prior variance of the non-zero elements of b, equal to
\code{scaled_prior_variance * var(Y)}.}

\item{elbo}{The value of the variational lower bound, or
\dQuote{ELBO} (objective function to be maximized), achieved at
each iteration of the IBSS fitting procedure.}

\item{fitted}{Vector of length n containing the fitted values of
the outcome.}

\item{sets}{Credible sets estimated from model fit; see
\code{\link{susie_get_cs}} for details.}

\item{pip}{A vector of length p giving the (marginal) posterior
inclusion probabilities for all p covariates.}

\item{z}{A vector of univariate z-scores.}

\item{niter}{Number of IBSS iterations that were performed.}

\item{converged}{\code{TRUE} or \code{FALSE} indicating whether
the IBSS converged to a solution within the chosen tolerance
level.}

\code{susie_suff_stat} returns also outputs:

\item{XtXr}{A p-vector of \code{t(X)} times the fitted values,
\code{X \%*\% colSums(alpha*mu)}.}
}
\description{
Performs Bayesian multiple linear regression of Y on
X; that is, this function fits the regression model \eqn{Y = \sum_l
  X b_{l=1}^L + e}, where elements of e are \emph{i.i.d.} normal with
zero mean and variance \code{residual_variance}, and
\eqn{\sum_{l=1}^L b_l} is a vector of length p representing the
effects to be estimated. The \dQuote{susie assumption} is that each
\eqn{b_l} has exactly one non-zero element. The prior on the
non-zero element is normal with zero mean and variance \code{var(Y)
  * scaled_prior_variance}. The model is fitted using the
\dQuote{Iterative Bayesian Stepwise Selection} (IBSS) algorithm.
See also \code{\link{susie_trendfilter}} for applying susie to
non-parametric regression, particularly changepoint problems.
}
\details{
\code{susie_suff_stat} performs sum of single-effect
linear regression with summary statistics. The required summary
data are either: \code{bhat}, \code{shat}, the p by p symmetric,
positive semidefinite correlation (or covariance) matrix \code{R},
the sample size \code{n}, and the variance of y; or the p by p
matrix \eqn{X'X}, the p-vector \eqn{X'y}, the sum of squares
\eqn{y'y}, and the sample size \code{n}. The summary statistics
should come from the same individuals. Both the columns of X and
the vector y should be centered to have mean zero before computing
these summary statistics; you may also want to scale each column of
X and y to have variance 1 (see examples).
}
\references{
G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
new approach to variable selection in regression, with application
to genetic fine-mapping. \emph{Journal of the Royal Statistical
Society, Series B} \url{https://doi.org/10.1101/501114}.
}
\seealso{
\code{\link{susie_rss}}
}
\keyword{internal}
