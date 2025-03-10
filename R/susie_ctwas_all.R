#' @title Compute sufficient statistics for susie.
#' 
#' @param X An n by p matrix of covariates.
#' 
#' @param y An n vector.
#' 
#' @param standardize Logical flag indicating whether to standardize
#'   columns of X to unit variance prior to fitting.
#' 
#' @return A list of sufficient statistics.
#' 
#' @importFrom methods as
#'
#' @examples
#' data(N2finemapping)
#' ss = compute_ss(N2finemapping$X, N2finemapping$Y[,1])
#' 
compute_ss = function(X, y, standardize = TRUE) {
  y = y - mean(y)
  is.sparse = !is.matrix(X)
  X = set_X_attributes(as.matrix(X),center=TRUE,scale = standardize)
  X = t((t(X) - attr(X,"scaled:center"))/attr(X,"scaled:scale"))
  XtX = crossprod(X)
  if(is.sparse)
    XtX = as(XtX,"dgCMatrix")
  Xty = c(y %*% X)
  n = length(y)
  yty = sum(y^2)
  return(list(XtX = XtX,Xty = Xty,yty = yty,n = n))
}
# @title Get objective function from data and susie fit object.
# @param data A flash data object.
# @param f A flash fit object.
# @keywords internal
get_objective = function (X, Y, s) {
  return(Eloglik(X,Y,s) - sum(s$KL))
}

# Expected loglikelihood for a susie fit.
Eloglik = function (X, Y, s) {
  n = nrow(X)
  return(-(n/2) * log(2*pi*s$sigma2) - (1/(2*s$sigma2)) * get_ER2(X,Y,s))
}

# Expected squared residuals.
# s$Xr is column sum of Xr_L
get_ER2 = function (X, Y, s) {
  Xr_L = compute_MXt(s$alpha * s$mu,X) # L by N matrix
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  return(sum((Y - s$Xr)^2) - sum(Xr_L^2) + sum(attr(X,"d") * t(postb2)))
}

# @title posterior expected loglikelihood for a single effect regression
# @param X an n by p matrix of covariates
# @param Y an n vector of regression outcome
# @param s2 the residual variance
# @param Eb the posterior mean of b (p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik = function (X, Y, s2, Eb, Eb2) {
  n = nrow(X)
  return(-0.5*n*log(2*pi*s2) - 0.5/s2*(sum(Y*Y)
                                       - 2*sum(Y*compute_Xb(X,Eb))
                                       + sum(attr(X,"d") * Eb2)))
}
# @title Get objective function from data and susie fit object.
# @param R p by p corelation matrix
# @param z length p vector
# @param s a susie fit object
get_objective_rss = function (R, z, s) 
  Eloglik_rss(s$sigma2, R,z,s) - sum(s$KL)

# @title expected loglikelihood for a susie fit
Eloglik_rss = function (sigma2, R, z, s) {
  d = sigma2 * attr(R,"eigen")$values + attr(R,"lambda")
  if(attr(R,"lambda") == 0)
    result = -(sum(d != 0)/2) * log(2*pi*sigma2) - 0.5*get_ER2_rss(sigma2, R,z,s)
  else
    result = -(length(z)/2)*log(2*pi) - 0.5*sum(log(d)) - 0.5*get_ER2_rss(sigma2, R,z,s)
  return(result)
}

# @title expected squared residuals
# @importFrom Matrix diag
get_ER2_rss = function (sigma2, R, z, s) {
  d = sigma2 * attr(R,"eigen")$values + attr(R,"lambda")
  Dinv = 1/d
  Dinv[is.infinite(Dinv)] = 0
  SinvR = attr(R,"eigen")$vectors %*%
          ((Dinv*attr(R,"eigen")$values) * t(attr(R,"eigen")$vectors))
  Utz = crossprod(attr(R,"eigen")$vectors,z)
  zSinvz = sum(Utz * (Dinv * Utz))

  Z = s$alpha * s$mu
  if(attr(R,"lambda") == 0)
    RSinvR = R/sigma2
  else
    RSinvR = R %*% SinvR
  RZ2 = sum((Z%*%RSinvR) * Z)

  zbar = colSums(Z)
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  return(zSinvz - 2*sum((SinvR %*% z) * zbar)
         + sum(zbar * (RSinvR %*% zbar))
         - RZ2 + sum(diag(RSinvR) * t(postb2)))
}

# @title posterior expected loglikelihood for a single effect regression
# @param R a p by p symmetric and positive semidefinite correlation matrix.
# @param Sigma residual_var * R + lambda I
# @param r residuals
# @param Eb the posterior mean of b (p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik_rss = function (R, Sigma, r, Ez, Ez2) {
  eigenS = attr(Sigma,'eigenS')
  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  rR = R %*% r
  SinvEz = eigenS$vectors %*% (Dinv * crossprod(eigenS$vectors, Ez))
  return(-0.5*(-2*sum(rR*SinvEz) +
               sum(attr(Sigma,"RjSinvRj") * as.vector(Ez2))))
}
# @title Get objective function from data and susie fit object.
# @param XtX a p by p matrix, X'X
# @param Xty a p vector, X'y,
# @param s a susie fit object
# @param yty a scaler, y'y, where y is centered to have mean 0
# @param n sample size
get_objective_ss = function (XtX, Xty, s, yty, n)
  Eloglik_ss(XtX,Xty,s,yty,n) - sum(s$KL)

# Expected loglikelihood for a susie fit.
Eloglik_ss = function (XtX, Xty, s, yty, n)
  -n/2*log(2*pi*s$sigma2) - 1/(2*s$sigma2) * get_ER2_ss(XtX,Xty,s,yty)

# Expected squared residuals.
get_ER2_ss = function (XtX, Xty, s, yty) {
  B = s$alpha * s$mu
  XB2 = sum((B %*% XtX) * B)
  betabar = colSums(B)
  d = attr(XtX,"d")
  postb2 = s$alpha * s$mu2 # Posterior second moment.
  return(yty - 2*sum(betabar * Xty) + sum(betabar * (XtX %*% betabar)) -
         XB2 + sum(d * t(postb2)))
}

# @title posterior expected loglikelihood for a single effect regression
# @param dXtX a p vector of diagonal elements of XtX
# @param Xty a p vector
# @param s2 the residual variance
# @param Eb the posterior mean of b (p vector) (alpha * mu)
# @param Eb2 the posterior second moment of b (p vector) (alpha * mu2)
SER_posterior_e_loglik_ss = function (dXtX, Xty, s2, Eb, Eb2)
  -0.5/s2 * (-2*sum(Eb*Xty) + sum(dXtX * as.vector(Eb2)))
# @title Estimate residual variance
# @param X an n by p matrix of covariantes
# @param Y an n vector of data
# @param s a susie fit
estimate_residual_variance = function (X, Y, s) {
  n = nrow(X)
  return((1/n)*get_ER2(X,Y,s))
}

# @title Estimate residual variance for summary statistics
# @param XtX a p by p matrix
# @param Xty a p vector
# @param s a susie fit
# @param yty a scaler, Y'Y, where Y is centered to have mean 0
# @param n sample size
estimate_residual_variance_ss = function (XtX, Xty, s, yty, n)
  (1/n)*get_ER2_ss(XtX,Xty,s,yty)
#' @title Initialize a susie object using regression coefficients
#'
#' @param coef_index An L-vector containing the the indices of the
#'   nonzero coefficients.
#'
#' @param coef_value An L-vector containing initial coefficient
#' estimates.
#'
#' @param p A scalar giving the number of variables.
#'
#' @return A list with elements \code{alpha}, \code{mu} and \code{mu2}
#'   to be used by \code{susie}.
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[sample(1:1000,4)] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#'
#' # Initialize susie to ground-truth coefficients.
#' s = susie_init_coef(which(beta != 0),beta[beta != 0],length(beta))
#' res = susie(X,y,L = 10,s_init=s)
#'
susie_init_coef = function (coef_index, coef_value, p) {
  L = length(coef_index)
  if (L <= 0)
    stop("Need at least one non-zero effect")
  if (!all(coef_value != 0))
    stop("Input coef_value must be non-zero for all its elements")
  if (L != length(coef_value))
    stop("Inputs coef_index and coef_value must of the same length")
  if (max(coef_index) > p)
    stop("Input coef_index exceeds the boundary of p")
  alpha = matrix(0,nrow = L,ncol = p)
  mu = matrix(0,nrow = L,ncol = p)
  for(i in 1:L){
    alpha[i,coef_index[i]] = 1
    mu[i,coef_index[i]] = coef_value[i]
  }
  out = list(alpha = alpha,mu = mu,mu2 = mu*mu)
  class(out) <- c("susie","list")
  return(out)
}

# Set default susie initialization.
init_setup = function (n, p, L, scaled_prior_variance, residual_variance,
                       prior_weights, null_weight, varY, standardize) {
  if (!is.numeric(scaled_prior_variance) || scaled_prior_variance < 0)
    stop("Scaled prior variance should be positive number")
  if (scaled_prior_variance > 1 && standardize)
    stop("Scaled prior variance should be no greater than 1 when ",
         "standardize = TRUE")
  if(is.null(residual_variance))
    residual_variance = varY
  if(is.null(prior_weights))
    prior_weights = rep(1/p,p)
  else
    prior_weights = prior_weights / sum(prior_weights)
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p
  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Xr     = rep(0,n),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           V      = scaled_prior_variance*varY,
           pi     = prior_weights)
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

# Update a susie fit object in order to initialize susie model.
init_finalize = function (s, X = NULL, Xr = NULL) {
  # different form susieR
  # if(length(s$V) == 1)
  #   s$V = rep(s$V,nrow(s$alpha))

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("Residual variance sigma2 must be positive (is your var(Y) zero?)")

  # check prior variance
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")

  # different form susieR
  # if (nrow(s$alpha) != length(s$V))
  #   stop("Input prior variance V must have length of nrow of alpha in ",
  #        "input object")

  # Update Xr.
  if (!missing(Xr))
    s$Xr = Xr
  if (!missing(X))
    s$Xr = compute_Xb(X,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}
# Set default susie initialization.
init_setup_rss = function(p, L, prior_variance, residual_variance,
                          prior_weights, null_weight) {
  if (!is.numeric(prior_variance) || prior_variance < 0)
    stop("Prior variance should be positive number.")
  if(!is.null(residual_variance) &&
     (residual_variance > 1 | residual_variance < 0))
    stop("Residual variance should be a scaler between 0 and 1")
  if (is.null(residual_variance))
    residual_variance = 1
  if (is.null(prior_weights))
    prior_weights = rep(1/p,p)
  if(length(prior_weights) != p)
    stop("Prior weights must have length p")
  if (p < L)
    L = p
  s = list(alpha  = matrix(1/p,nrow = L,ncol = p),
           mu     = matrix(0,nrow = L,ncol = p),
           mu2    = matrix(0,nrow = L,ncol = p),
           Rz     = rep(0,p),
           KL     = rep(as.numeric(NA),L),
           lbf    = rep(as.numeric(NA),L),
           lbf_variable = matrix(as.numeric(NA),L,p),
           sigma2 = residual_variance,
           V      = prior_variance,
           pi     = prior_weights)
  if (is.null(null_weight))
    s$null_index = 0
  else
    s$null_index = p
  class(s) = "susie"
  return(s)
}

# Update a susie fit object in order to initialize susie model.
init_finalize_rss = function (s, R = NULL, Rz = NULL) {
  # different from susieR
  # if(length(s$V) == 1)
  #   s$V = rep(s$V,nrow(s$alpha))

  # Check sigma2.
  if (!is.numeric(s$sigma2))
    stop("Input residual variance sigma2 must be numeric")

  # Avoid problems with dimension if input is a 1 x 1 matrix.
  s$sigma2 = as.numeric(s$sigma2)
  if (length(s$sigma2) != 1)
    stop("Input residual variance sigma2 must be a scalar")
  if (s$sigma2 <= 0)
    stop("residual variance sigma2 must be positive (is your var(Y) zero?)")

  # Check prior variance.
  if (!is.numeric(s$V))
    stop("Input prior variance must be numeric")
  if (!all(s$V >= 0))
    stop("prior variance must be non-negative")
  if (!all(dim(s$mu) == dim(s$mu2)))
    stop("dimension of mu and mu2 in input object do not match")
  if (!all(dim(s$mu) == dim(s$alpha)))
    stop("dimension of mu and alpha in input object do not match")
  # different from susieR
  # if (nrow(s$alpha) != length(s$V))
  #   stop("Input prior variance V must have length of nrow of alpha in ",
  #        "input object")

  # Update Rz.
  if (!missing(Rz))
    s$Rz = Rz
  if (!missing(R))
    s$Rz = compute_Xb(R,colSums(s$mu * s$alpha))

  # Reset KL and lbf.
  s$KL = rep(as.numeric(NA),nrow(s$alpha))
  s$lbf = rep(as.numeric(NA),nrow(s$alpha))
  class(s) = "susie"
  return(s)
}
#' @title Extract regression coefficients from susie fit
#' 
#' @param object A susie fit.
#'
#' @param \dots Additional arguments passed to the generic \code{coef}
#'   method.
#' 
#' @return A p+1 vector, the first element being an intercept, and the
#'   remaining p elements being estimated regression coefficients.
#'
#' @importFrom stats coef
#'
#' @method coef susie
#' 
coef.susie = function (object, ...) {
  s = object
  return(c(s$intercept,colSums(s$alpha*s$mu)/s$X_column_scale_factors))
}

#' @title Predict outcomes or extract coefficients from susie fit.
#' 
#' @param object A susie fit.
#' 
#' @param newx A new value for X at which to do predictions.
#' 
#' @param type The type of output. For \code{type = "response"},
#'   predicted or fitted outcomes are returned; for \code{type =
#'   "coefficients"}, the estimated coefficients are returned.
#'
#' @param \dots Other arguments used by generic predict function. These
#'   extra arguments are not used here.
#' 
#' @return For \code{type = "response"}, predicted or fitted outcomes
#'   are returned; for \code{type = "coefficients"}, the estimated
#'   coefficients are returned.
#' 
#' @importFrom stats coef
#'
#' @method predict susie
#' 
predict.susie = function (object, newx = NULL,
                          type = c("response","coefficients"), ...) {
  s = object
  type = match.arg(type)
  if (type == "coefficients") {
    if (!missing(newx))
      stop("Do not supply newx when predicting coefficients")
    return(coef(s))
  }
  if (missing(newx))
    return(s$fitted)
  return(drop(s$intercept + newx %*% coef(s)[-1]))
}
# @title sets the attributes for the R matrix
# @param R a p by p LD matrix
# @param r_tol tolerance level for eigen value check of positive
#   semidefinite matrix of R.
# @return R with attribute e.g., attr(R, 'eigenR') is the eigen
#   decomposition of R.
set_R_attributes = function (R, r_tol) {
  if (is.null(attr(R,"eigen")))
    eigenR = eigen(R,symmetric = TRUE)
  else
    eigenR = attr(R,"eigen")

  # Drop small eigenvalues.
  eigenR$values[abs(eigenR$values) < r_tol] = 0
  if(any(eigenR$values < 0)) {
    min_lambda = min(eigenR$values)
    eigenR$values[eigenR$values < 0] = 0
    warning(paste0("The input correlation matrix has negative eigenvalues ",
                   "(smallest one is ", min_lambda, "). The correlation ",
                   "matrix is adjusted such that these negative eigenvalues ",
                   "are now zeros. You can ignore this message, only if you ",
                   "believe the negative eigenvalue is result of numerical ",
                   "rounding errors."))
  }
  res = eigenR$vectors %*% (t(eigenR$vectors) * eigenR$values)

  attr(res,"eigen") = eigenR
  attr(res,"d") = diag(res)
  attr(res,"scaled:scale") = rep(1,length = nrow(R))
  return(res)
}
# @title sets three attributes for matrix X
# @param X an n by p data matrix that can be either a trend filtering
#   matrix or a regular dense/sparse matrix
# @param center boolean indicating centered by column means or not
# @param scale boolean indicating scaled by column standard deviations or not
# @return X with three attributes e.g. `attr(X, 'scaled:center') is a
# p vector of column means of X if center=TRUE, a p vector of zeros
# otherwise. 'attr(X, 'scaled:scale') is a p vector of columan standard
# deviations of X if scale=TRUE, a p vector of 1s otherwise. 'attr(X,
# 'd') is a p vector of column sums of X.standardized^2,' where
# X.standardized is the matrix X centered by attr(X, 'scaled:center')
# and scaled by attr(X, 'scaled:scale').
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
set_X_attributes = function(X, center = TRUE, scale = TRUE) {
    
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
      
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))
    X.std = (t(X) - cm)/csd
    
    # Set three attributes for X.
    attr(X,"d") = rowSums(X.std * X.std)
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

# @title computes column standard deviations for any type of matrix
# @details This should give the same result as matrixStats::colSds(X),
#   but allows for sparse matrices as well as dense ones.
# @param X an n by p matrix of any type, e.g. sparse, dense.
# @return a p vector of column standard deviations.
#
#' @importFrom Matrix colSums
compute_colSds = function(X) {
  n = nrow(X)
  return(sqrt((colSums(X^2)/n - (colSums(X)/n)^2)*(n/(n-1))))
}
#' @rdname single_effect_regression
#'
#' @title Bayesian single-effect linear regression
#' 
#' @description These methods fit the regression model \eqn{y = Xb +
#'   e}, where elements of e are \emph{i.i.d.}  \eqn{N(0,s^2)}, and b is
#'   a p-vector of effects to be estimated. The assumption is that b has
#'   exactly one non-zero element, with all elements equally likely to
#'   be non-zero. The prior on the coefficient of the non-zero element
#'   is \eqn{N(0,V)}.
#'
#' @details \code{single_effect_regression_ss} performs single-effect
#' linear regression with summary data, in which only the statistcs
#' \eqn{X^Ty} and diagonal elements of \eqn{X^TX} are provided to the
#' method.
#' 
#' \code{single_effect_regression_rss} performs single-effect linear
#' regression with z scores. That is, this function fits the
#' regression model \eqn{z = R*b + e}, where e is \eqn{N(0,Sigma)},
#' \eqn{Sigma = residual_var*R + lambda*I}, and the b is a p-vector of
#' effects to be estimated. The assumption is that b has exactly one
#' non-zero element, with all elements equally likely to be non-zero.
#' The prior on the non-zero element is \eqn{N(0,V)}. The required
#' summary data are the p-vector \code{z} and the p by p matrix
#' \code{Sigma}. The summary statistics should come from the same
#' individuals.
#' 
#' @param Y An n-vector.
#' 
#' @param X An n by p matrix of covariates.
#' 
#' @param V A scalar giving the (initial) prior variance
#' 
#' @param residual_variance The residual variance.
#' 
#' @param prior_weights A p-vector of prior weights.
#' 
#' @param optimize_V The optimization method to use for fitting the
#'   prior variance.
#' 
#' @param check_null_threshold Scalar specifying threshold on the
#'   log-scale to compare likelihood between current estimate and zero
#'   the null.
#' 
#' @return A list with the following elements:
#' 
#' \item{alpha}{Vector of posterior inclusion probabilities;
#'   \code{alpha[i]} is posterior probability that the ith coefficient
#'   is non-zero.}
#' 
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#' 
#' \item{mu2}{Vector of posterior second moments (conditional on
#'   inclusion).}
#' 
#' \item{lbf}{Vector of log-Bayes factors for each variable.}
#' 
#' \item{lbf_model}{Log-Bayes factor for the single effect regression.}
#'
#' \code{single_effect_regression} and \code{single_effect_regression_ss}
#' additionally output:
#' 
#' \item{V}{Prior variance (after optimization if \code{optimize_V !=
#'   "none"}).}
#' 
#' \item{loglik}{The log-likelihood, \eqn{\log p(y | X, V)}.}
#'
#' @importFrom stats uniroot
#' @importFrom stats optim
#' @importFrom Matrix colSums
#'
#' @keywords internal
#' 
single_effect_regression =
  function (Y, X, V, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  optimize_V = match.arg(optimize_V)
  Xty = compute_Xty(X,Y)
  betahat = (1/attr(X,"d")) * Xty
  shat2 = residual_variance/attr(X,"d")
  if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X),ncol(X))
  if (optimize_V != "EM" && optimize_V != "none")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
        alpha = NULL,post_mean2 = NULL,V_init = V,
        check_null_threshold = check_null_threshold)

  # log(bf) for each SNP
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 
  maxlbf = max(lbf)

  # w is proportional to BF, but subtract max for numerical stability.
  w = exp(lbf - maxlbf)
  
  # Posterior prob for each SNP.
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  post_var = (1/V + attr(X,"d")/residual_variance)^(-1) # Posterior variance.
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # Second moment.
  
  # BF for single effect model.
  lbf_model = maxlbf + log(weighted_sum_w)
  loglik = lbf_model + sum(dnorm(Y,0,sqrt(residual_variance),log = TRUE))

  if(optimize_V == "EM")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                alpha,post_mean2,
                                check_null_threshold = check_null_threshold)

  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              lbf_model = lbf_model,V = V,loglik = loglik))
}

# Estimate prior variance.
est_V_uniroot = function (betahat, shat2, prior_weights) {
  V.u = uniroot(negloglik.grad.logscale,c(-10,10),extendInt = "upX",
                betahat = betahat,shat2 = shat2,prior_weights = prior_weights)
  return(exp(V.u$root))
}

optimize_prior_variance = function (optimize_V, betahat, shat2, prior_weights,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL, check_null_threshold = 0) {
  V = V_init
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
      lV = optim(par = log(max(c(betahat^2-shat2,1),na.rm = TRUE)),
          fn = neg.loglik.logscale,betahat = betahat,shat2 = shat2,
          prior_weights = prior_weights,method = "Brent",lower = -30,
          upper = 15)$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik.logscale(lV, betahat = betahat,shat2 = shat2,prior_weights = prior_weights) > 
         neg.loglik.logscale(log(V), betahat = betahat,
                             shat2 = shat2,prior_weights = prior_weights)){
        lV = log(V)
      }
      V = exp(lV)
    } else if (optimize_V == "uniroot")
      V = est_V_uniroot(betahat,shat2,prior_weights)
    else if (optimize_V == "EM")
      V = sum(alpha * post_mean2)
    else
     stop("Invalid option for optimize_V method")
  }
  
  # Set V exactly 0 if that beats the numerical value by
  # check_null_threshold in loglik. check_null_threshold = 0.1 is
  # exp(0.1) = 1.1 on likelihood scale; it means that for parsimony
  # reasons we set estiate of V to zero, if its numerical estimate is
  # only "negligibly" different from zero. We use a likelihood ratio
  # of exp(check_null_threshold) to define "negligible" in this
  # context. This is fairly modest condition compared to, say, a
  # formal LRT with p-value 0.05. But the idea is to be lenient to
  # non-zeros estimates unless they are indeed small enough to be
  # neglible. See more intuition at
  # https://stephens999.github.io/fiveMinuteStats/LR_and_BF.html
  if (loglik(0,betahat,shat2,prior_weights) +
      check_null_threshold >= loglik(V,betahat,shat2,prior_weights))
    V = 0
  return(V)
}

# In these functions, s2 represents residual_variance, and shat2 is an
# estimate of it.

# The log likelihood function for SER model (based on summary data
# betahat, shat2) as a function of prior variance V.
# 
#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik = function (V, betahat, shat2, prior_weights) {

  #log(bf) for each SNP
  lbf = dnorm(betahat,0,sqrt(V+shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 

  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) # w = BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w) + maxlbf)
}

neg.loglik.logscale = function(lV,betahat,shat2,prior_weights)
  -loglik(exp(lV),betahat,shat2,prior_weights)

#' @importFrom Matrix colSums
#' @importFrom stats dnorm
loglik.grad = function(V, betahat, shat2, prior_weights) {

  # log(bf) for each SNP.
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 

  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) # w = BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  return(sum(alpha * lbf.grad(V,shat2,betahat^2/shat2)))
}

# Define loglikelihood and gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale = function (lV, betahat, shat2, prior_weights)
  -exp(lV) * loglik.grad(exp(lV),betahat,shat2,prior_weights)

# Vector of gradients of logBF_j for each j, with respect to prior
# variance V.
lbf.grad = function (V, shat2, T2) {
  l = 0.5*(1/(V + shat2)) * ((shat2/(V + shat2))*T2 - 1)
  l[is.nan(l)] = 0
  return(l)
}

lbf = function (V, shat2, T2) {
  l = 0.5*log(shat2/(V + shat2)) + 0.5*T2*(V/(V + shat2))
  l[is.nan(l)] = 0
  return(l)
}
#' @rdname single_effect_regression
#' 
#' @param z A p-vector of z scores.
#' 
#' @param Sigma \code{residual_var*R + lambda*I}
#' 
#' @keywords internal
#' 
single_effect_regression_rss =
  function (z, Sigma, V = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  p = length(z)
  shat2 = 1/attr(Sigma,"RjSinvRj")
  if (is.null(prior_weights))
    prior_weights = rep(1/p,p)

  if (optimize_V != "EM" && optimize_V != "none") 
    V = optimize_prior_variance_rss(optimize_V,z,Sigma,prior_weights,
                                    alpha = NULL,post_mean2 = NULL,V_init = V,
                                    check_null_threshold=check_null_threshold)

  lbf = sapply(1:p, function(j)
    -0.5 * log(1 + (V/shat2[j])) +
     0.5 * (V/(1 + (V/shat2[j]))) * sum(attr(Sigma,"SinvRj")[,j] * z)^2
  )

  # Deal with special case of infinite shat2 (e.g., happens if X does not
  # vary).
  lbf[is.infinite(shat2)] = 0 

  # w is proportional to BF, but subtract max for numerical stability.
  maxlbf = max(lbf)
  w = exp(lbf-maxlbf)
  
  # posterior prob on each SNP
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w

  post_var = (attr(Sigma,"RjSinvRj") + 1/V)^(-1) # Posterior variance.
  post_mean = sapply(1:p,function(j) (post_var[j]) *
              sum(attr(Sigma,"SinvRj")[,j] * z))
  post_mean2 = post_var + post_mean^2 # Second moment.
  lbf_model = maxlbf + log(weighted_sum_w) # Analogue of loglik in the
                                           # non-summary case.

  if (optimize_V=="EM") 
    V = optimize_prior_variance_rss(optimize_V,z,Sigma,prior_weights,
        alpha,post_mean2,check_null_threshold = check_null_threshold)
  
  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              V = V,lbf_model = lbf_model))
}

loglik_rss = function (V, z, Sigma, prior_weights) {
  p = length(z)
  shat2 = 1/attr(Sigma,"RjSinvRj")

  # log(bf) for each SNP.
  lbf = sapply(1:p,function (j)
    -0.5 * log(1 + (V/shat2[j])) +
     0.5 * (V/(1 + (V/shat2[j]))) * sum(attr(Sigma,"SinvRj")[,j] * z)^2)
  
  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 

  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w = BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w) + maxlbf)
}

neg.loglik_z.logscale_rss = function (lV, z, Sigma, prior_weights)
  -loglik_rss(exp(lV),z,Sigma,prior_weights)

optimize_prior_variance_rss = function (optimize_V, z, Sigma, prior_weights,
                                        alpha = NULL, post_mean2 = NULL,
                                        V_init = NULL,
                                        check_null_threshold = 0) {
  V = V_init
  if (optimize_V != "simple") {
    if(optimize_V == "optim") {
      lV = optim(par = log(max(c((colSums(attr(Sigma,"SinvRj") * z)^2) -
                       (1/attr(Sigma,"RjSinvRj")),1e-6),na.rm = TRUE)),
                 fn = neg.loglik_z.logscale_rss,z = z,Sigma = Sigma,
                 prior_weights = prior_weights,method = "Brent",
                 lower = -30,upper = 15)$par
      ## if the estimated one is worse than current one, don't change it.
      if(neg.loglik_z.logscale_rss(lV, z = z,Sigma = Sigma,prior_weights = prior_weights) > 
         neg.loglik_z.logscale_rss(log(V), z = z,Sigma = Sigma,prior_weights = prior_weights)){
        lV = log(V)
      }
      V = exp(lV)
    } else if (optimize_V == "EM")
      V = sum(alpha * post_mean2)
    else
      stop("Invalid option for optimize_V")
  }
  
  # Set V exactly 0 if that beats the numerical value. By
  # check_null_threshold in loglik. check_null_threshold = 0.1 is
  # exp(0.1) = 1.1 on likelihood scale; it means that for parsimony
  # reasons we set estiate of V to zero, if its numerical estimate is
  # only "negligibly" different from zero. We use a likelihood ratio
  # of exp(check_null_threshold) to define "negligible" in this
  # context. This is fairly modest condition compared to, say, a
  # formal LRT with p-value 0.05. But the idea is to be lenient to
  # non-zeros estimates unless they are indeed small enough to be
  # neglible. See more intuition at
  # https://stephens999.github.io/fiveMinuteStats/LR_and_BF.html
  if (loglik_rss(0,z,Sigma,prior_weights) + check_null_threshold >=
      loglik_rss(V,z,Sigma,prior_weights))
    V = 0
  return(V)
}
#' @rdname single_effect_regression
#' 
#' @param Xty A p-vector.
#' 
#' @param dXtX A p-vector containing the diagonal elements of
#'   \code{crossprod(X)}.
#' 
#' @importFrom stats uniroot
#' @importFrom stats optim
#'
#' @keywords internal
#' 
single_effect_regression_ss =
  function (Xty, dXtX, V = 1, residual_variance = 1, prior_weights = NULL,
            optimize_V = c("none", "optim", "uniroot", "EM", "simple"),
            check_null_threshold = 0) {
  optimize_V = match.arg(optimize_V)
  betahat = (1/dXtX) * Xty
  shat2 = residual_variance/dXtX
  if (is.null(prior_weights))
    prior_weights = rep(1/length(dXtX),length(dXtX))

  if (optimize_V != "EM" && optimize_V != "none")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,
                                alpha = NULL,post_mean2 = NULL,V_init = V,
                                check_null_threshold = check_null_threshold)

  # log(bf) for each SNP.
  lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
        dnorm(betahat,0,sqrt(shat2),log = TRUE)

  # Deal with special case of infinite shat2 (e.g., happens if X does
  # not vary).
  lbf[is.infinite(shat2)] = 0 

  # w is proportional to BF, but subtract max for numerical stability
  # posterior prob on each SNP.
  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) 
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w

  post_var = (1/V + dXtX/residual_variance)^(-1) # Posterior variance.
  post_mean = (1/residual_variance) * post_var * Xty
  post_mean2 = post_var + post_mean^2 # Second moment.
  lbf_model = maxlbf + log(weighted_sum_w) # Analogue of loglik in the
                                           # non-summary case.

  if (optimize_V == "EM")
    V = optimize_prior_variance(optimize_V,betahat,shat2,prior_weights,alpha,
          post_mean2,check_null_threshold = check_null_threshold)
  return(list(alpha = alpha,mu = post_mean,mu2 = post_mean2,lbf = lbf,
              V = V,lbf_model = lbf_model))
}


# @title Computes standardized.X %*% b using sparse multiplication trick
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @param b a p vector
# @return an n vector
# 
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_Xb = function (X, b) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  
  # Scale Xb.
  if (!is.null(attr(X,"matrix.type")))

    # When X is a trend filtering matrix.
    scaled.Xb = compute_tf_Xb(attr(X,"order"),b/csd)
  else
      
    # When X is an ordinary sparse/dense matrix.
    scaled.Xb = tcrossprod(X,t(b/csd))
  
  # Center Xb.
  Xb = scaled.Xb - sum(cm*b/csd)
  return(as.numeric(Xb))
}

# @title Computes t(standardized.X) %*% y using sparse multiplication trick
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @param y an n vector
# @return a p vector
# 
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty = function (X, y) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  ytX = crossprod(y,X)
  
  # Scale Xty.
  if (!is.null(attr(X,"matrix.type")))

    # When X is a trend filtering matrix.
    scaled.Xty = compute_tf_Xty(attr(X,"order"),y)/csd
  else

    # When X is an ordinary sparse/dense matrix.
    scaled.Xty = t(ytX/csd)
  
  # Center Xty.
  centered.scaled.Xty = scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

# @title Computes M %* %t(standardized.X) using sparse multiplication trick
# @param M a L by p matrix
# @param X an n by p unstandardized matrix with three attributes:
# attr(X,"scaled:center"), attr(X,"scaled:scale") and attr(X,"d")
# @return a L by n matrix
# 
#' @importFrom Matrix t
compute_MXt = function (M, X) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  
  if (!is.null(attr(X,"matrix.type")))

    # When X is a trend filtering matrix.
    return(as.matrix(t(apply(M,1,function(b) compute_Xb(X,b)))))
  else

    # When X is an ordinary sparse/dense matrix.
    return(as.matrix(tcrossprod(M,sweep(X,2,csd,"/")) -
                     drop(tcrossprod(M,t(cm/csd)))))

  # This should be the same as
  #
  #   t(apply(M, 1, function(b) compute_Xb(X, b))))
  #
  # as well as
  #
  #   M %*% (t(X)/csd) - drop(tcrossprod(M,t(cm/csd)))
  #
  # but should be more memory-efficient.
}
#' @rdname susie
#'
#' @title Sum of Single Effects (SuSiE) Regression
#'
#' @description Performs Bayesian multiple linear regression of Y on
#'   X; that is, this function fits the regression model \eqn{Y = \sum_l
#'   X b_{l=1}^L + e}, where elements of e are \emph{i.i.d.} normal with
#'   zero mean and variance \code{residual_variance}, and
#'   \eqn{\sum_{l=1}^L b_l} is a vector of length p representing the
#'   effects to be estimated. The \dQuote{susie assumption} is that each
#'   \eqn{b_l} has exactly one non-zero element. The prior on the
#'   non-zero element is normal with zero mean and variance \code{var(Y)
#'   * scaled_prior_variance}. The model is fitted using the
#'   \dQuote{Iterative Bayesian Stepwise Selection} (IBSS) algorithm.
#'   See also \code{\link{susie_trendfilter}} for applying susie to
#'   non-parametric regression, particularly changepoint problems.
#'
#' @details \code{susie_suff_stat} performs sum of single-effect
#' linear regression with summary statistics. The required summary
#' data are either: \code{bhat}, \code{shat}, the p by p symmetric,
#' positive semidefinite correlation (or covariance) matrix \code{R},
#' the sample size \code{n}, and the variance of y; or the p by p
#' matrix \eqn{X'X}, the p-vector \eqn{X'y}, the sum of squares
#' \eqn{y'y}, and the sample size \code{n}. The summary statistics
#' should come from the same individuals. Both the columns of X and
#' the vector y should be centered to have mean zero before computing
#' these summary statistics; you may also want to scale each column of
#' X and y to have variance 1 (see examples).
#'
#' @param X An n by p matrix of covariates.
#'
#' @param Y The observed responses, a vector of length n.
#'
#' @param L Number of components (nonzero coefficients) in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param scaled_prior_variance The scaled prior variance. This is
#'   either a scalar or a vector of length \code{L}. The prior variance
#'   of each non-zero element of b is set to \code{var(Y) *
#'   scaled_prior_variance}. If \code{estimate_prior_variance = TRUE},
#'   this provides initial estimates of the prior variances.
#'
#' @param residual_variance Variance of the residual. If
#'   \code{estimate_residual_variance = TRUE}, this value provides the
#'   initial estimate of the residual variance. By default, it is
#'   \code{var(Y)}.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that corresponding column of X has a
#'   nonzero effect on the outcome, Y.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X (or XtX and Xty) to unit variance prior to
#'   fitting. Note that \code{scaled_prior_variance} specifies the prior
#'   on the coefficients of X \emph{after} standardization (if it is
#'   performed). If you do not standardize, you may need to think more
#'   carefully about specifying \code{scaled_prior_variance}. Whatever
#'   your choice, the coefficients returned by \code{coef} are given for
#'   \code{X} on the original input scale. Any column of \code{X} that
#'   has zero variance is not standardized.
#'
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; it \code{intercept = FALSE}, the intercept is set to
#'   zero. Setting \code{intercept = FALSE} is generally not
#'   recommended.
#'
#' @param estimate_residual_variance If
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated, using \code{residual_variance} as an initial value. If
#'   \code{estimate_residual_variance = FALSE}, the residual variance is
#'   fixed to the value supplied by \code{residual_variance}.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L effects). If provided,
#'   \code{scaled_prior_variance} is then used as an initial value for
#'   the optimization. When \code{estimate_prior_variance = FALSE}, the
#'   prior variance for each of the L effects is determined by the
#'   value supplied to \code{scaled_prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance. When \code{estimate_prior_method = "simple"} is used, the
#'   likelihood at the specified prior variance is compared to the
#'   likelihood at a variance of zero, and the setting with the larger
#'   likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated
#'   prior variance is smaller than this tolerance value.
#'
#' @param residual_variance_upperbound Upper limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param s_init A previous susie fit with which to initialize.
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param compute_univariate_zscore If \code{compute_univariate_zscore
#'   = TRUE}, the univariate regression z-scores are outputted for each
#'   variable.
#'
#' @param na.rm Drop any missing values in Y from both X and Y.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param track_fit If \code{track_fit = TRUE}, \code{trace}
#'   is also returned containing detailed information about the
#'   estimates at each iteration of the IBSS fitting procedure.
#'
#' @param residual_variance_lowerbound Lower limit on the estimated
#'   residual variance. It is only relevant when
#'   \code{estimate_residual_variance = TRUE}.
#'
#' @param refine If \code{refine = TRUE}, we use a procedure to help
#'   SuSiE get out of local optimum.
#'
#' @return A \code{"susie"} object with some or all of the following
#'   elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An L by p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{An L by p matrix of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{Xr}{A vector of length n, equal to \code{X \%*\% colSums(alpha
#'   * mu)}.}
#'
#' \item{lbf}{log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{log-Bayes Factor for each variable and single effect.}
#'
#' \item{intercept}{Intercept (fixed or estimated).}
#'
#' \item{sigma2}{Residual variance (fixed or estimated).}
#'
#' \item{V}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(Y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{fitted}{Vector of length n containing the fitted values of
#'   the outcome.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{susie_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{z}{A vector of univariate z-scores.}
#'
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' \code{susie_suff_stat} returns also outputs:
#'
#' \item{XtXr}{A p-vector of \code{t(X)} times the fitted values,
#'   \code{X \%*\% colSums(alpha*mu)}.}
#'
#' @references
#'
#' G. Wang, A. Sarkar, P. Carbonetto and M. Stephens (2020). A simple
#'   new approach to variable selection in regression, with application
#'   to genetic fine-mapping. \emph{Journal of the Royal Statistical
#'   Society, Series B} \url{https://doi.org/10.1101/501114}.
#'
#' @seealso \code{\link{susie_rss}}
#'
#' @examples
#' # susie example.
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' res1 = susie(X,y,L = 10)
#' plot(beta,coef(res1)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#' plot(y,predict(res1))
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' # susie_suff_stat example.
#' input_ss = compute_ss(X,y,standardize = TRUE)
#' res2 = with(input_ss,
#'             susie_suff_stat(XtX = XtX,Xty = Xty,yty = yty,n = n,L = 10))
#' plot(coef(res1)[-1],coef(res2)[-1])
#' abline(a = 0,b = 1,col = "skyblue",lty = "dashed")
#'
#' @importFrom stats var
#' @importFrom utils modifyList
#'
susie <- function (X,Y,L = min(10,ncol(X)),
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
                   prior_tol = 1e-9,
                   residual_variance_upperbound = Inf,
                   s_init = NULL,
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   compute_univariate_zscore = FALSE,
                   na.rm = FALSE,
                   max_iter = 100,
                   tol = 1e-3,
                   verbose = FALSE,
                   track_fit = FALSE,
                   residual_variance_lowerbound = var(drop(Y))/1e4,
                   refine = FALSE) {

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)

  # Check input X.
  if (!(is.double(X) & is.matrix(X)) &
      !inherits(X,"CsparseMatrix") &
      is.null(attr(X,"matrix.type")))
    stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
         "a trend filtering matrix")
  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight) && is.null(attr(X,"matrix.type"))) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(X) * (1 - null_weight),ncol(X)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    X = cbind(X,0)
  }
  if (any(is.na(X)))
    stop("Input X must not contain missing values")
  if (any(is.na(Y))) {
    if (na.rm) {
      samples_kept = which(!is.na(Y))
      Y = Y[samples_kept]
      X = X[samples_kept,]
    } else
      stop("Input Y must not contain missing values")
  }

  # Check input Y.
  p = ncol(X)
  n = nrow(X)
  mean_y = mean(Y)

  # Center and scale input.
  if (intercept)
    Y = Y - mean_y
  X = set_X_attributes(X,center = intercept,scale = standardize)

  # Initialize susie fit.
  s = init_setup(n,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,as.numeric(var(Y)),standardize)
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init, verbose=FALSE)
    # Then prune or expand
    s_init = susie_prune_single_effects(s_init, L, s$V, verbose)
    s = modifyList(s,s_init)
    s = init_finalize(s,X = X)
  } else {
    s = init_finalize(s)
  }
  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect(X,Y,s,estimate_prior_variance,estimate_prior_method,
                           check_null_threshold)
    if (verbose)
      print(paste0("objective:",get_objective(X,Y,s)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective(X,Y,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      s$sigma2 = pmax(residual_variance_lowerbound,
                      estimate_residual_variance(X,Y,s))
      if (s$sigma2 > residual_variance_upperbound)
        s$sigma2 = residual_variance_upperbound
      if (verbose)
        print(paste0("objective:",get_objective(X,Y,s)))
    }
  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)]
  s$elbo = elbo
  s$niter = i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  if (intercept) {

    # Estimate unshrunk intercept.
    s$intercept = mean_y - sum(attr(X,"scaled:center") *
      (colSums(s$alpha * s$mu)/attr(X,"scaled:scale")))
    s$fitted = s$Xr + mean_y
  } else {
    s$intercept = 0
    s$fitted = s$Xr
  }
  s$fitted = drop(s$fitted)
  names(s$fitted) = `if`(is.null(names(Y)), rownames(X), names(Y))

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s,coverage = coverage,X = X,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }
  # different from susieR, with null weight, this seems doesn't work
  # if (!is.null(colnames(X))) {
  #   variable_names = colnames(X)
  #   if (!is.null(null_weight))
  #     variable_names = c("null", variable_names)
  #   colnames(s$alpha) = variable_names
  #   colnames(s$mu) = variable_names
  #   colnames(s$mu2) = variable_names
  #   colnames(s$lbf_variable) = variable_names
  #   names(s$pip) = variable_names
  # }
  # report z-scores from univariate regression.
  if (compute_univariate_zscore) {
    if (!is.null(null_weight) && null_weight != 0)
      X = X[,1:(ncol(X) - 1)]
    s$z = calc_z(X,Y,center = intercept,scale = standardize)
  }

  # For prediction.
  s$X_column_scale_factors = attr(X,"scaled:scale")

  if(refine){
    if(!is.null(null_weight) && null_weight!=0 && !compute_univariate_zscore){
      ## if null_weight is specified, and the extra 0 column is not removed from compute_univariate_zscore,
      ## we remove it here
      X = X[,1:(ncol(X) - 1)]
    }
    conti = TRUE
    while(conti){
      m = list()
      for(cs in 1:length(s$sets$cs)){
        if(!missing(s_init) && !is.null(s_init)){
          warning('The given s_init is not used in refinement.')
        }
        pw = rep(1, ncol(X))
        pw[s$sets$cs[[cs]]] = 0
        s2 = susie(X,Y,L = L,
                   scaled_prior_variance = scaled_prior_variance,residual_variance = residual_variance,
                   prior_weights = pw,s_init = NULL,
                   null_weight = null_weight,standardize = standardize,intercept = intercept,
                   estimate_residual_variance = estimate_residual_variance,
                   estimate_prior_variance = estimate_prior_variance,
                   estimate_prior_method = estimate_prior_method,
                   check_null_threshold = check_null_threshold,
                   prior_tol = prior_tol,coverage = coverage,
                   residual_variance_upperbound = residual_variance_upperbound,
                   min_abs_corr = min_abs_corr,compute_univariate_zscore = FALSE,
                   na.rm = na.rm,max_iter = max_iter,tol = tol,
                   verbose = FALSE,track_fit = FALSE,residual_variance_lowerbound = var(drop(Y))/1e4,
                   refine = FALSE)
        sinit2 = s2[c('alpha', 'mu', 'mu2')]
        class(sinit2) = 'susie'
        s3 = susie(X,Y,L = L,
                   scaled_prior_variance = scaled_prior_variance,residual_variance = residual_variance,
                   prior_weights = NULL, s_init = sinit2,
                   null_weight = null_weight,standardize = standardize,intercept = intercept,
                   estimate_residual_variance = estimate_residual_variance,
                   estimate_prior_variance = estimate_prior_variance,
                   estimate_prior_method = estimate_prior_method,
                   check_null_threshold = check_null_threshold,
                   prior_tol = prior_tol,coverage = coverage,
                   residual_variance_upperbound = residual_variance_upperbound,
                   min_abs_corr = min_abs_corr,compute_univariate_zscore = FALSE,
                   na.rm = na.rm,max_iter = max_iter,tol = tol,
                   verbose = FALSE,track_fit = FALSE,residual_variance_lowerbound = var(drop(Y))/1e4,
                   refine = FALSE)
        m = c(m, list(s3))
      }
      elbo = sapply(m, function(x) susie_get_objective(x))
      if((max(elbo) - susie_get_objective(s)) <= 0){
        conti=FALSE
      }else{
        s = m[[which.max(elbo)]]
      }
    }
  }
  return(s)
}
#' @title Sum of Single Effects (SuSiE) Regression using summary statistics
#'
#' @description \code{susie_rss} performs sum of single-effect linear
#'   regression with z scores; all posterior calculations are for
#'   z-scores. This function fits the regression model \eqn{z = \sum_l
#'   R*b_l + e}, where e is \eqn{N(0,R)} and \eqn{\sum_l b_l} is a
#'   p-vector of effects to be estimated. The required summary data are
#'   the p by p correlation matrix, \code{R}, and the p-vector \code{z}.
#'
#' @param z A p-vector of z scores.
#'
#' @param R A p by p symmetric, positive semidefinite correlation
#' matrix.
#'
#' @param maf Minor allele frequency; to be used along with
#'   \code{maf_thresh} to filter input summary statistics.
#'
#' @param maf_thresh Variants having a minor allele frequency smaller
#'   than this threshold are not used.
#'
#' @param z_ld_weight This feature is not recommended. The weights
#'   assigned to the z scores in the LD matrix. If \code{z_ld_weight >
#'   0}, the LD matrix used in the model is \code{cov2cor((1-w)*R +
#'   w*tcrossprod(z))}, where \code{w = z_ld_weight}.
#'
#' @param L Number of components (nonzero coefficients) in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param prior_variance The prior variance. It is either a scalar or
#'   a vector of length L.
#'
#' @param residual_variance Variance of the residual.
#'   If it is not specified, we set it to 1.
#'
#' @param prior_weights A vector of length p, in which each entry
#'   gives the prior probability that SNP j has non-zero effect.
#'
#' @param null_weight Prior probability of no effect (a number between
#'   0 and 1, and cannot be exactly 1).
#'
#' @param estimate_residual_variance The residual variance is
#'   fixed to the value supplied by \code{residual_variance}. We don't
#'   estimate residual variance in susie_rss.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L effects). If provided,
#'   \code{prior_variance} is then used as an initial value for
#'   the optimization. When \code{estimate_prior_variance = FALSE}, the
#'   prior variance for each of the L effects is determined by the
#'   value supplied to \code{prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior
#'   variance. When \code{estimate_prior_method = "simple"} is used, the
#'   likelihood at the specified prior variance is compared to the
#'   likelihood at a variance of zero, and the setting with the larger
#'   likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimate with the null, and set the prior variance to
#'   zero unless the log-likelihood using the estimate is larger by this
#'   threshold amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated
#'   prior variance is smaller than this tolerance value.
#'
#' @param max_iter Maximum number of IBSS iterations to perform.
#'
#' @param s_init A previous susie fit with which to initialize.
#'
#' @param intercept_value The intercept. (The intercept cannot be
#'   estimated from centered summary data.) This setting will be used by
#'   \code{coef} to assign an intercept value, mainly for consistency
#'   with \code{susie}. Set to \code{NULL} if you want \code{coef} not
#'   to include an intercept term (and so only a p-vector is returned).
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#'   and a summary of the optimization settings, are printed to the
#'   console.
#'
#' @param track_fit If \code{track_fit = TRUE}, \code{trace}
#'   is also returned containing detailed information about the
#'   estimates at each iteration of the IBSS fitting procedure.
#'
#' @param check_R If \code{check_R = TRUE}, check that \code{R} is
#'   positive semidefinite.
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param refine If \code{refine = TRUE}, we use a procedure to help
#'   SuSiE get out of local optimum.
#'
#' @return A \code{"susie"} object with some or all of the following
#'   elements:
#'
#' \item{alpha}{An L by p matrix of posterior inclusion probabilites.}
#'
#' \item{mu}{An L by p matrix of posterior means, conditional on
#'   inclusion.}
#'
#' \item{mu2}{An L by p matrix of posterior second moments,
#'   conditional on inclusion.}
#'
#' \item{lbf}{log-Bayes Factor for each single effect.}
#'
#' \item{lbf_variable}{log-Bayes Factor for each variable and single effect.}
#'
#' \item{intercept}{Fixed Intercept.}
#'
#' \item{sigma2}{Fixed Residual variance.}
#'
#' \item{V}{Prior variance of the non-zero elements of b, equal to
#'   \code{scaled_prior_variance * var(Y)}.}
#'
#' \item{elbo}{The value of the variational lower bound, or
#'   \dQuote{ELBO} (objective function to be maximized), achieved at
#'   each iteration of the IBSS fitting procedure.}
#'
#' \item{fitted}{Vector of length n containing the fitted values of
#'   the outcome.}
#'
#' \item{sets}{Credible sets estimated from model fit; see
#'   \code{\link{susie_get_cs}} for details.}
#'
#' \item{pip}{A vector of length p giving the (marginal) posterior
#'   inclusion probabilities for all p covariates.}
#'
#' \item{niter}{Number of IBSS iterations that were performed.}
#'
#' \item{converged}{\code{TRUE} or \code{FALSE} indicating whether
#'   the IBSS converged to a solution within the chosen tolerance
#'   level.}
#'
#' \item{Rr}{An p-vector of \code{t(X)} times fitted values, \code{X
#'   \%*\% colSums(alpha*mu)}.}
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#'
#' input_ss <- compute_ss(X,y,standardize = TRUE)
#' ss   <- univariate_regression(X,y)
#' R    <- with(input_ss,cov2cor(XtX))
#' zhat <- with(ss,betahat/sebetahat)
#' res  <- susie_rss(zhat,R,L = 10)
#'
susie_rss = function (z, R, maf = NULL, maf_thresh = 0, z_ld_weight = 0,
                      L = 10, prior_variance = 50, residual_variance = NULL,
                      prior_weights = NULL, null_weight = NULL,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = TRUE,
                      estimate_prior_method = c("optim", "EM", "simple"),
                      check_null_threshold = 0, prior_tol = 1e-9,
                      max_iter = 100, s_init = NULL, intercept_value = 0,
                      coverage = 0.95, min_abs_corr = 0.5,
                      tol = 1e-03, verbose = FALSE, track_fit = FALSE,
                      check_R = FALSE, r_tol = 1e-08, refine = FALSE) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (", nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))
  if (!is_symmetric_matrix(R))
    stop("R is not a symmetric matrix")
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix, or a sparse matrix")

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(z))
      stop(paste0("The length of maf does not agree with expected ",length(z)))
    id = which(maf > maf_thresh)
    R = R[id,id]
    z = z[id]
  }

  if (any(is.infinite(z)))
    stop("z contains infinite value")

  # Check for NAs in R.
  if (any(is.na(R)))
    stop("R matrix contains missing values")

  # Replace NAs in z with zeros.
  if (any(is.na(z))) {
    warning("NA values in z-scores are replaced with 0")
    z[is.na(z)] = 0
  }

  if (check_R && any(attr(R,"eigen")$values < -r_tol)){
    semi_pd = check_semi_pd(R,r_tol)
    if (!semi_pd$status){
      stop(paste0("The correlation matrix (",nrow(R)," by ",ncol(R),
                  ") is not a positive semidefinite matrix. The smallest ",
                  "eigenvalue is ",min(semi_pd$eigenvalues),"."))
    }
  }

  # Modify R as needed.
  if (z_ld_weight > 0) {
    warning('From version 0.11.0, the non-zero z_ld_weight is no longer recommended.')
    R = muffled_cov2cor((1-z_ld_weight)*R + z_ld_weight * tcrossprod(z))
    R = (R + t(R))/2
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R) * (1 - null_weight),ncol(R)),null_weight)
    else
      prior_weights = c(prior_weights * (1 - null_weight),null_weight)
    R = cbind(rbind(R,0),0)
    z = c(z,0)
    # different from susieR
    prior_variance = cbind(prior_variance, 0)
  }

  if(estimate_residual_variance){
    warning("SuSiE-RSS no longer estimates residual variance, since we found it didn't help.")
    estimate_residual_variance = FALSE
  }

  if (!is.null(residual_variance) &&
      (residual_variance > 1 | residual_variance < 0))
    stop("Residual variance should be a scalar between 0 and 1")
  if (is.null(residual_variance))
    residual_variance = 1

  s = susie_suff_stat(XtX = R, Xty = z, n = length(z), yty = length(z)-1,
                      L = L, scaled_prior_variance = prior_variance,
                      residual_variance = residual_variance,
                      estimate_residual_variance = FALSE,
                      estimate_prior_variance = estimate_prior_variance,
                      estimate_prior_method = estimate_prior_method,
                      check_null_threshold = check_null_threshold, prior_tol = prior_tol,
                      r_tol = r_tol, prior_weights = prior_weights,
                      null_weight = NULL, standardize = FALSE,
                      max_iter = max_iter, s_init = s_init, intercept_value = intercept_value,
                      coverage = coverage, min_abs_corr = min_abs_corr,
                      tol = tol, verbose = verbose, track_fit = track_fit, check_input = FALSE,
                      refine = refine)
  s$fitted = s$Xtfitted
  s$Rr = s$XtXr
  s$Xtfitted = s$XtXr = NULL

  # different from susieR, this is a bug of susie_rss, make it consistent to susie
  if (!is.null(null_weight)) {
    s$null_index <- length(z)
    s$pip <- s$pip[1:(length(z)-1)]
  }
  return(s)
}

# Performs sum of single-effect (SuSiE) linear regression with z
# scores (with lambda). The summary data required are the p by p
# correlation matrix R, the p vector z. The summary stats should come
# from the same individuals. This function fits the regression model z
# = sum_l Rb_l + e, where e is N(0,residual_variance * R + lambda I)
# and the sum_l b_l is a p vector of effects to be estimated. The
# assumption is that each b_l has exactly one non-zero element, with
# all elements equally likely to be non-zero. The prior on the
# non-zero element is N(0,var = prior_variance).
#
#' @importFrom stats optimize
susie_rss_lambda = function(z, R, maf = NULL, maf_thresh = 0,
                            L = 10, lambda = 0,
                            prior_variance = 50, residual_variance = NULL,
                            r_tol = 1e-08, prior_weights = NULL,
                            null_weight = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            estimate_prior_method = c("optim", "EM", "simple"),
                            check_null_threshold = 0, prior_tol = 1e-9,
                            max_iter = 100, s_init = NULL, intercept_value = 0,
                            coverage = 0.95, min_abs_corr = 0.5,
                            tol = 1e-3, verbose = FALSE, track_fit = FALSE,
                            check_R = TRUE, check_z = FALSE) {

  # Check input R.
  if (nrow(R) != length(z))
    stop(paste0("The dimension of correlation matrix (",nrow(R)," by ",
                ncol(R),") does not agree with expected (",length(z)," by ",
                length(z),")"))
  if (!is_symmetric_matrix(R))
    stop("R is not a symmetric matrix")
  if (!(is.double(R) & is.matrix(R)) & !inherits(R,"CsparseMatrix"))
    stop("Input R must be a double-precision matrix or a sparse matrix")

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(z))
      stop(paste0("The length of maf does not agree with expected ",length(z)))
    id = which(maf > maf_thresh)
    R = R[id,id]
    z = z[id]
  }

  if (any(is.infinite(z)))
    stop("z contains infinite values")

  # Check for NAs in R.
  if (any(is.na(R)))
    stop("R matrix contains missing values")

  # Replace NAs in z with zero.
  if (any(is.na(z))) {
    warning("NA values in z-scores are replaced with 0")
    z[is.na(z)] = 0
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (missing(prior_weights))
      prior_weights = c(rep(1/ncol(R)*(1-null_weight),ncol(R)),null_weight)
    else
      prior_weights = c(prior_weights * (1-null_weight),null_weight)
    R = cbind(rbind(R,0),0)
    z = c(z,0)
  }

  # Eigen decomposition for R, fileter on eigenvalues.
  p = ncol(R)
  attr(R,"eigen") = eigen(R,symmetric = TRUE)
  if (check_R && any(attr(R,"eigen")$values < -r_tol))
    stop(paste0("The correlation matrix (",nrow(R)," by ",ncol(R),
                "is not a positive semidefinite matrix. ",
                "The smallest eigenvalue is ",min(attr(R,"eigen")$values),
                ". You can bypass this by \"check_R = FALSE\" which instead ",
                "sets negative eigenvalues to 0 to allow for continued ",
                "computations."))

  # Check whether z in space spanned by the non-zero eigenvectors of R.
  if (check_z) {
    proj = check_projection(R,z)
    if (!proj$status)
      warning("Input z does not lie in the space of non-zero eigenvectors ",
              "of R.")
    else
      message("Input z is in space spanned by the non-zero eigenvectors of ",
              "R.")
  }
  R = set_R_attributes(R,r_tol)

  if (lambda == 'estimate'){
    colspace = which(attr(R,"eigen")$values > 0)
    if(length(colspace) == length(z)){
      lambda = 0
    }else{
      znull = crossprod(attr(R,"eigen")$vectors[,-colspace], z) # U2^T z
      lambda = sum(znull^2)/length(znull)
    }
  }

  # Initialize susie fit.
  s = init_setup_rss(p,L,prior_variance,residual_variance,prior_weights,
                     null_weight)
  if (!missing(s_init) && !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init, verbose=FALSE)
    # Then prune or expand
    s_init = susie_prune_single_effects(s_init, L, s$V, verbose)
    s = modifyList(s,s_init)
    s = init_finalize_rss(s,R = R)
  } else
    s = init_finalize_rss(s)

  s$sigma2 = s$sigma2 - lambda
  estimate_prior_method = match.arg(estimate_prior_method)

  # Intialize elbo to NA.
  elbo = rep(NA,max_iter+1)
  elbo[1] = -Inf;
  tracking = list()

  attr(R,"lambda") = lambda
  Sigma = update_Sigma(R,s$sigma2,z)  # sigma2*R + lambda*I

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect_rss(R,z,s,Sigma,estimate_prior_variance,
                               estimate_prior_method,check_null_threshold)
    if (verbose)
      print(paste0("before estimate sigma2 objective:",
                   get_objective_rss(R,z,s)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance, before the update.
    elbo[i+1] = get_objective_rss(R,z,s)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      if (lambda == 0) {
        est_sigma2 = (1/sum(attr(R,"eigen")$values != 0))*get_ER2_rss(1,R,z,s)
        if (est_sigma2 < 0)
          stop("Estimating residual variance failed: the estimated value ",
               "is negative")
        if (est_sigma2 > 1)
          est_sigma2 = 1
      } else {
        est_sigma2 = optimize(Eloglik_rss, interval = c(1e-4, 1-lambda),
                              R = R, z = z, s = s, maximum = TRUE)$maximum
        if(Eloglik_rss(est_sigma2, R, z, s) < Eloglik_rss(1-lambda, R, z, s)){
          est_sigma2 = 1-lambda
        }
      }
      s$sigma2 = est_sigma2

      if (verbose)
        print(paste0("after estimate sigma2 objective:",
                     get_objective_rss(R,z,s)))
      Sigma = update_Sigma(R,s$sigma2,z)
    }
  }

  # Remove first (infinite) entry, and trailing NAs.
  elbo = elbo[2:(i+1)]
  s$elbo = elbo
  s$niter = i
  s$lambda = lambda

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$fitted = s$Rz

  s$X_column_scale_factors = attr(R,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(R)
    s$sets = susie_get_cs(s,coverage = coverage,Xcorr = R,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }
  if (!is.null(names(z))) {
    variable_names = names(z)
    if (!is.null(null_weight))
      variable_names = c("null",variable_names)
    colnames(s$alpha) = variable_names
    colnames(s$mu) = variable_names
    colnames(s$mu2) = variable_names
    colnames(s$lbf_variable) = variable_names
    names(s$pip) = variable_names
  }

  return(s)
}

update_Sigma = function (R, sigma2, z) {
  Sigma = sigma2*R + attr(R,"lambda") * diag(length(z))
  eigenS = attr(R,"eigen")
  eigenS$values = sigma2 * eigenS$values + attr(R,"lambda")

  Dinv = 1/(eigenS$values)
  Dinv[is.infinite(Dinv)] = 0
  attr(Sigma,"eigenS") = eigenS

  # Sigma^(-1) R_j = U (sigma2 D + lambda)^(-1) D U^T e_j
  attr(Sigma,"SinvRj") = eigenS$vectors %*% (Dinv * attr(R,"eigen")$values *
                           t(eigenS$vectors))

  if (attr(R,"lambda") == 0)
    attr(Sigma,"RjSinvRj") = attr(R,"d")/sigma2
  else {
    tmp = t(eigenS$vectors)
    attr(Sigma,"RjSinvRj") =
      colSums(tmp * (Dinv*(attr(R,"eigen")$values^2) * tmp))
  }

  return(Sigma)
}

#' @rdname susie
#'
#' @param bhat A p-vector of estimated effects.
#'
#' @param shat A p-vector of standard errors.
#'
#' @param R A p by p symmetric, positive semidefinite matrix. It
#'   can be \eqn{X'X}, the covariance matrix \eqn{X'X/(n-1)}, or a
#'   correlation matrix. It should be estimated from the same samples
#'   used to compute \code{bhat} and \code{shat}. Using an out-of-sample
#'   matrix may produce unreliable results.
#'
#' @param n The sample size.
#'
#' @param var_y The sample variance of y, defined as \eqn{y'y/(n-1)}.
#'   When the sample variance cannot be provided, the coefficients
#'   (returned from \code{coef}) are computed on the "standardized" X, y
#'   scale.
#'
#' @param XtX A p by p matrix \eqn{X'X} in which the columns of X
#'   are centered to have mean zero.
#'
#' @param Xty A p-vector \eqn{X'y} in which y and the columns of X are
#'   centered to have mean zero.
#'
#' @param yty A scalar \eqn{y'y} in which y is centered to have mean
#'   zero.
#'
#' @param maf Minor allele frequency; to be used along with
#'   \code{maf_thresh} to filter input summary statistics.
#'
#' @param maf_thresh Variants having a minor allele frequency smaller
#'   than this threshold are not used.
#'
#' @param r_tol Tolerance level for eigenvalue check of positive
#'   semidefinite matrix of R.
#'
#' @param intercept_value The intercept. (The intercept cannot be
#'   estimated from centered summary data.) This setting will be used by
#'   \code{coef} to assign an intercept value, mainly for consistency
#'   with \code{susie}. Set to \code{NULL} if you want \code{coef} not
#'   to include an intercept term (and so only a p-vector is returned).
#'
#' @param check_input If \code{check_input = TRUE},
#'   \code{susie_suff_stat} performs additional checks on \code{XtX} and
#'   \code{Xty}. The checks are: (1) check that \code{XtX} is positive
#'   semidefinite; (2) check that \code{Xty} is in the space spanned by
#'   the non-zero eigenvectors of \code{XtX}.
#'
susie_suff_stat = function (bhat, shat, R, n, var_y, XtX, Xty, yty,
                            maf = NULL, maf_thresh = 0, L = 10,
                            scaled_prior_variance = 0.2,
                            residual_variance = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE,
                            estimate_prior_method = c("optim","EM","simple"),
                            check_null_threshold = 0, prior_tol = 1e-9,
                            r_tol = 1e-08, prior_weights = NULL,
                            null_weight = NULL, standardize = TRUE,
                            max_iter = 100, s_init = NULL,
                            intercept_value = 0, coverage = 0.95,
                            min_abs_corr = 0.5, tol = 1e-3, verbose = FALSE,
                            track_fit = FALSE, check_input = FALSE, refine = FALSE) {

  # Process input estimate_prior_method.
  estimate_prior_method = match.arg(estimate_prior_method)

  if (missing(n))
    stop("n must be provided")

  # Check sufficient statistics.
  missing_bhat = c(missing(bhat), missing(shat), missing(R))
  missing_XtX = c(missing(XtX), missing(Xty), missing(yty))

  if (all(missing_bhat) & all(missing_XtX))
    stop("Please provide either all of bhat, shat, R, n, var_y or all of ",
         "XtX, Xty, yty, n")
  if (any(missing_bhat) & any(missing_XtX))
    stop("Please provide either all of bhat, shat, R, n, var_y or all of ",
         "XtX, Xty, yty, n")
  if (all(missing_bhat) & any(missing_XtX))
    stop("Please provide all of XtX, Xty, yty, n")
  if (all(missing_XtX) & any(missing_bhat))
    stop("Please provide all of bhat, shat, R, n, var_y")
  if ((!any(missing_XtX)) & (!all(missing_bhat)))
    warning("Only using information from XtX, Xty, yty, n")
  if (!any(missing_bhat)) {
    if (!all(missing_XtX))
      warning("Only using information from bhat, shat, R, n, var_y")

    # Compute XtX, Xty, yty from bhat, shat, R, n, var_y.
    if (length(shat) == 1)
      shat = rep(shat,length(bhat))
    if (length(bhat) != length(shat))
      stop("The length of bhat does not agree with length of shat")
    if (anyNA(bhat) || anyNA(shat))
      stop("The input summary statistics have missing values")
    if (any(shat == 0))
      stop("shat contains one or more zeros")

    that = bhat/shat
    that[is.na(that)] = 0
    R2 = that^2/(that^2 + n-2)
    sigma2 = (n-1)*(1-R2)/(n-2)

    # Convert any input R to correlation matrix.
    # If R has 0 colums and rows, cov2cor produces NaN and warning.
    X0 = diag(R) == 0
    R = muffled_cov2cor(R)

    # Change the columns and rows with NaN to 0.
    if (sum(X0) > 0)
      R[X0,] = R[,X0] = 0

    if (missing(var_y)) {
      XtX = (n-1)*R
      Xty = sqrt(sigma2) * sqrt(n-1) * that
      var_y = 1
    } else {
      XtXdiag = var_y * sigma2/(shat^2)
      Xty = that * var_y * sigma2/shat
      XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
    }
    yty = var_y * (n-1)
  }

  # Check input XtX.
  if (ncol(XtX) != length(Xty))
    stop(paste0("The dimension of XtX (",nrow(XtX)," by ",ncol(XtX),
                ") does not agree with expected (",length(Xty)," by ",
                length(Xty),")"))
  if (!is_symmetric_matrix(XtX))
    stop("XtX is not a symmetric matrix")

  # MAF filter.
  if (!is.null(maf)) {
    if (length(maf) != length(Xty))
      stop(paste("The length of maf does not agree with expected",length(Xty)))
    id = which(maf > maf_thresh)
    XtX = XtX[id,id]
    Xty = Xty[id]
  }

  if (any(is.infinite(Xty)))
    stop("Xty contains infinite values")
  if (!(is.double(XtX) & is.matrix(XtX)) & !inherits(XtX,"CsparseMatrix"))
    stop("Input X must be a double-precision matrix, or a sparse matrix")
  if (any(is.na(XtX)))
    stop("XtX matrix contains NAs")

  if (check_input) {

    # Check whether XtX is positive semidefinite.
    semi_pd = check_semi_pd(XtX,r_tol)
    if (!semi_pd$status)
      stop("XtX is not a positive semidefinite matrix")

    # Check whether Xty in space spanned by the non-zero eigenvectors of XtX
    proj = check_projection(semi_pd$matrix,Xty)
    if (!proj$status)
      warning("Xty does not lie in the space of the non-zero eigenvectors ",
              "of XtX")
  }

  if (is.numeric(null_weight) && null_weight == 0)
    null_weight = NULL
  if (!is.null(null_weight)) {
    if (!is.numeric(null_weight))
      stop("Null weight must be numeric")
    if (null_weight < 0 || null_weight >= 1)
      stop("Null weight must be between 0 and 1")
    if (is.null(prior_weights))
      prior_weights = c(rep(1/ncol(XtX)*(1-null_weight),ncol(XtX)),null_weight)
    else
      prior_weights = c(prior_weights*(1 - null_weight),null_weight)
    XtX = cbind(rbind(XtX,0),0)
    Xty = c(Xty,0)
  }

  p = ncol(XtX)

  if (standardize) {
    dXtX = diag(XtX)
    csd = sqrt(dXtX/(n-1))
    csd[csd == 0] = 1
    XtX = (1/csd) * t((1/csd) * XtX)
    Xty = (1/csd) * Xty
  } else
    csd = rep(1, length = p)
  attr(XtX,"d") = diag(XtX)
  attr(XtX,"scaled:scale") = csd

  # Initialize susie fit.
  s = init_setup(0,p,L,scaled_prior_variance,residual_variance,prior_weights,
                 null_weight,yty/(n-1),standardize)
  s$Xr = NULL
  s$XtXr = rep(0,p)

  if (!missing(s_init)&& !is.null(s_init)) {
    if (!inherits(s_init,"susie"))
      stop("s_init should be a susie object")
    if (max(s_init$alpha) > 1 || min(s_init$alpha) < 0)
      stop("s_init$alpha has invalid values outside range [0,1]; please ",
           "check your input")
    # First, remove effects with s_init$V = 0
    s_init = susie_prune_single_effects(s_init, verbose=FALSE)
    # Then prune or expand
    s_init = susie_prune_single_effects(s_init, L, s$V, verbose)
    s = modifyList(s,s_init)
    s = init_finalize(s,X = XtX)
    s$XtXr = s$Xr
    s$Xr = NULL
  } else
    s = init_finalize(s)

  # Initialize elbo to NA.
  elbo = rep(as.numeric(NA),max_iter + 1)
  elbo[1] = -Inf;
  tracking = list()

  for (i in 1:max_iter) {
    if (track_fit)
      tracking[[i]] = susie_slim(s)
    s = update_each_effect_ss(XtX,Xty,s,estimate_prior_variance,
                              estimate_prior_method,check_null_threshold)

    if (verbose)
      print(paste0("objective:",get_objective_ss(XtX,Xty,s,yty,n)))

    # Compute objective before updating residual variance because part
    # of the objective s$kl has already been computed under the
    # residual variance before the update.
    elbo[i+1] = get_objective_ss(XtX,Xty,s,yty,n)
    if ((elbo[i+1] - elbo[i]) < tol) {
      s$converged = TRUE
      break
    }
    if (estimate_residual_variance) {
      est_sigma2 = estimate_residual_variance_ss(XtX,Xty,s,yty,n)
      if (est_sigma2 < 0)
        stop("Estimating residual variance failed: the estimated value ",
             "is negative")
      s$sigma2 = est_sigma2
      if (verbose)
        print(paste0("objective:",get_objective_ss(XtX,Xty,s,yty,n)))
    }
  }
  elbo = elbo[2:(i+1)] # Remove first (infinite) entry, and trailing NAs.
  s$elbo = elbo
  s$niter = i

  if (is.null(s$converged)) {
    warning(paste("IBSS algorithm did not converge in",max_iter,"iterations!"))
    s$converged = FALSE
  }

  s$intercept = intercept_value
  s$Xtfitted = s$XtXr

  s$X_column_scale_factors = attr(XtX,"scaled:scale")

  if (track_fit)
    s$trace = tracking

  # SuSiE CS and PIP.
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    R = muffled_cov2cor(XtX)
    s$sets = susie_get_cs(s,coverage = coverage,Xcorr = R,
                          min_abs_corr = min_abs_corr)
    s$pip = susie_get_pip(s,prune_by_cs = FALSE,prior_tol = prior_tol)
  }

  if(refine){
    if(!is.null(null_weight) && null_weight!=0){
      ## if null_weight is specified
      ## we remove the extra 0 column
      XtX = XtX[1:(ncol(XtX)-1), 1:(ncol(XtX)-1)]
      Xty = Xty[1:(ncol(XtX)-1)]
    }
    conti = TRUE
    while(conti){
      m = list()
      for(cs in 1:length(s$sets$cs)){
        if(!missing(s_init) && !is.null(s_init)){
          warning('The given s_init is not used in refinement.')
        }
        pw = rep(1, ncol(XtX))
        pw[s$sets$cs[[cs]]] = 0
        s2 = susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                             prior_weights = pw, s_init = NULL,
                             scaled_prior_variance = scaled_prior_variance,
                             residual_variance = residual_variance,
                             estimate_residual_variance = estimate_residual_variance,
                             estimate_prior_variance = estimate_prior_variance,
                             estimate_prior_method = estimate_prior_method,
                             check_null_threshold = check_null_threshold, prior_tol = prior_tol,
                             r_tol = r_tol, max_iter = max_iter,
                             null_weight = NULL, standardize = standardize,
                             intercept_value = intercept_value, coverage = coverage,
                             min_abs_corr = min_abs_corr, tol = tol, verbose = FALSE,
                             track_fit = FALSE, check_input = FALSE, refine = FALSE)
        sinit2 = s2[c('alpha', 'mu', 'mu2')]
        class(sinit2) = 'susie'
        s3 = susie_suff_stat(XtX = XtX, Xty = Xty, yty = yty, n = n, L = L,
                             prior_weights = NULL, s_init = sinit2,
                             scaled_prior_variance = scaled_prior_variance,
                             residual_variance = residual_variance,
                             estimate_residual_variance = estimate_residual_variance,
                             estimate_prior_variance = estimate_prior_variance,
                             estimate_prior_method = estimate_prior_method,
                             check_null_threshold = check_null_threshold, prior_tol = prior_tol,
                             r_tol = r_tol, max_iter = max_iter,
                             null_weight = NULL, standardize = standardize,
                             intercept_value = intercept_value, coverage = coverage,
                             min_abs_corr = min_abs_corr, tol = tol, verbose = FALSE,
                             track_fit = FALSE, check_input = FALSE, refine = FALSE)
        m = c(m, list(s3))
      }
      elbo = sapply(m, function(x) susie_get_objective(x))
      if((max(elbo) - susie_get_objective(s)) <= 0){
        conti=FALSE
      }else{
        s = m[[which.max(elbo)]]
      }
    }
  }

  return(s)
}

# @title Check whether A is positive semidefinite
# @param A a symmetric matrix
# @return a list of result:
# \item{matrix}{The matrix with eigen decomposition}
# \item{status}{whether A is positive semidefinite}
# \item{eigenvalues}{eigenvalues of A truncated by r_tol}
check_semi_pd = function (A, tol) {
  attr(A,"eigen") = eigen(A,symmetric = TRUE)
  eigenvalues = attr(A,"eigen")$values
  eigenvalues[abs(eigenvalues) < tol] = 0
  return(list(matrix = A,
              status = !any(eigenvalues < 0),
              eigenvalues = eigenvalues))
}

# @title Check whether b is in space spanned by the non-zero eigenvectors
#   of A
# @param A a p by p matrix
# @param b a length p vector
# @return a list of result:
# \item{status}{whether b in space spanned by the non-zero
#  eigenvectors of A}
# \item{msg}{msg gives the difference between the projected b and b if
#   status is FALSE}
check_projection = function (A, b) {
  if (is.null(attr(A,"eigen")))
    attr(A,"eigen") = eigen(A,symmetric = TRUE)
  B = attr(A,"eigen")$vectors[,attr(A,"eigen")$values > .Machine$double.eps]
  msg = all.equal(as.vector(B %*% crossprod(B,b)),as.vector(b),check.names = FALSE)
  if (!is.character(msg))
    return(list(status = TRUE,msg = NA))
  else
    return(list(status = FALSE,msg = msg))
}
#' @rdname susie_get_methods
#'
#' @title Inferences From Fitted SuSiE Model
#'
#' @description These functions access basic properties or draw
#'   inferences from a fitted susie model.
#'
#' @param res A susie fit, typically an output from
#'   \code{\link{susie}} or one of its variants. For
#'   \code{susie_get_pip} and \code{susie_get_cs}, this may instead be
#'   the posterior inclusion probability matrix, \code{alpha}.
#'
#' @param last_only If \code{last_only = FALSE}, return the ELBO from
#'   all iterations; otherwise return the ELBO from the last iteration
#'   only.
#'
#' @param warning_tol Warn if ELBO is decreasing by this
#'   tolerance level.
#'
#' @return \code{susie_get_objective} returns the evidence lower bound
#' (ELBO) achieved by the fitted susie model and, optionally, at each
#' iteration of the IBSS fitting procedure.
#'
#' \code{susie_get_residual_variance} returns the (estimated or
#' fixed) residual variance parameter.
#'
#' \code{susie_get_prior_variance} returns the (estimated or fixed)
#' prior variance parameters.
#'
#' \code{susie_get_posterior_mean} returns the posterior mean for the
#' regression coefficients of the fitted susie model.
#'
#' \code{susie_get_posterior_sd} returns the posterior standard
#' deviation for coefficients of the fitted susie model.
#'
#' \code{susie_get_niter} returns the number of model fitting
#' iterations performed.
#'
#' \code{susie_get_pip} returns a vector containing the posterior
#' inclusion probabilities (PIPs) for all variables.
#'
#' \code{susie_get_lfsr} returns a vector containing the average lfsr
#' across variables for each single-effect, weighted by the posterior
#' inclusion probability (alpha).
#'
#' \code{susie_get_posterior_samples} returns a list containing the
#' effect sizes samples and causal status with two components: \code{b},
#' an \code{num_variables} x \code{num_samples} matrix of effect
#' sizes; \code{gamma}, an \code{num_variables} x \code{num_samples}
#' matrix of causal status random draws.
#'
#' \code{susie_get_cs} returns credible sets (CSs) from a susie fit,
#' as well as summaries of correlation among the variables included in
#' each CS. If desired, one can filter out CSs that do not meet a
#' specified \dQuote{purity} threshold; to do this, either \code{X} or
#' \code{Xcorr} must be supplied. It returns a list with the following
#' elements:
#'
#' \item{cs}{A list in which each list element is a vector containing
#'   the indices of the variables in the CS.}
#'
#' \item{coverage}{The nominal coverage specified for each CS.}
#'
#' \item{purity}{If \code{X} or \code{Xcorr} iis provided), the
#'   purity of each CS.}
#'
#' \item{cs_index}{If \code{X} or \code{Xcorr} is provided) the index
#'   (number between 1 and L) of each reported CS in the supplied susie
#'   fit.}
#'
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' s = susie(X,y,L = 10)
#' susie_get_objective(s)
#' susie_get_objective(s, last_only=FALSE)
#' susie_get_residual_variance(s)
#' susie_get_prior_variance(s)
#' susie_get_posterior_mean(s)
#' susie_get_posterior_sd(s)
#' susie_get_niter(s)
#' susie_get_pip(s)
#' susie_get_lfsr(s)
#'
susie_get_objective = function (res, last_only = TRUE, warning_tol = 1e-6) {
  if (!all(diff(res$elbo) >= (-1*warning_tol)))
    warning("Objective is decreasing")
  if (last_only)
    return(res$elbo[length(res$elbo)])
  else
    return(res$elbo)
}

#' @rdname susie_get_methods
#'
susie_get_posterior_mean = function (res, prior_tol = 1e-9) {

  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else
    include_idx = 1:nrow(res$alpha)

  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0)
    return(colSums((res$alpha*res$mu)[include_idx,,drop=FALSE])/
           res$X_column_scale_factors)
  else
    return(numeric(ncol(res$mu)))
}

#' @rdname susie_get_methods
#'
susie_get_posterior_sd = function (res, prior_tol = 1e-9) {

  # Drop the single-effects with estimated prior of zero.
  if (is.numeric(res$V))
    include_idx = which(res$V > prior_tol)
  else
    include_idx = 1:nrow(res$alpha)

  # Now extract relevant rows from alpha matrix.
  if (length(include_idx) > 0)
    return(sqrt(colSums((res$alpha * res$mu2 -
                         (res$alpha*res$mu)^2)[include_idx,,drop=FALSE]))/
           (res$X_column_scale_factors))
  else
    return(numeric(ncol(res$mu)))
}

#' @rdname susie_get_methods
#'
susie_get_niter = function (res)
  res$niter

#' @rdname susie_get_methods
#'
susie_get_prior_variance = function (res)
  res$V

#' @rdname susie_get_methods
#'
susie_get_residual_variance = function (res)
  res$sigma2

#' @rdname susie_get_methods
#'
#' @importFrom stats pnorm
#'
susie_get_lfsr = function (res) {
  pos_prob = pnorm(0,mean = t(res$mu),sd = sqrt(res$mu2 - res$mu^2))
  neg_prob = 1 - pos_prob
  return(1 - rowSums(res$alpha * t(pmax(pos_prob,neg_prob))))
}

#' @rdname susie_get_methods
#'
#' @param susie_fit A susie fit, an output from \code{\link{susie}}.
#'
#' @param num_samples The number of draws from the posterior
#'   distribution.
#'
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#'
susie_get_posterior_samples <- function (susie_fit, num_samples) {

  # Remove effects having estimated prior variance equals zero.
  if (is.numeric(susie_fit$V))
    include_idx = which(susie_fit$V > 1e-9)
  else
    include_idx = 1:nrow(susie_fit$alpha)

  posterior_mean = sweep(susie_fit$mu,2,susie_fit$X_column_scale_factors,"/")
  posterior_sd = sweep(sqrt(susie_fit$mu2 - (susie_fit$mu)^2),2,
                       susie_fit$X_column_scale_factors,"/")

  pip = susie_fit$alpha
  L = nrow(pip)
  num_snps = ncol(pip)
  b_samples = matrix(as.numeric(NA),num_snps,num_samples)
  gamma_samples = matrix(as.numeric(NA),num_snps,num_samples)
  for (sample_i in 1:num_samples) {
    b = 0
    if (length(include_idx) > 0) {
      for (l in include_idx) {
        gamma_l = rmultinom(1,1,pip[l,])
        effect_size = rnorm(1,mean = posterior_mean[l,which(gamma_l != 0)],
                            sd = posterior_sd[l,which(gamma_l != 0)])
        b_l = gamma_l * effect_size
        b = b + b_l
      }
    }
    b_samples[, sample_i] = b
    gamma_samples[, sample_i] = as.numeric(b != 0)
  }
  return(list(b = b_samples,gamma = gamma_samples))
}

#' @rdname susie_get_methods
#'
#' @param X n by p matrix of values of the p variables (covariates) in
#'   n samples. When provided, correlation between variables will be
#'   computed and used to remove CSs whose minimum correlation among
#'   variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables
#'   (covariates). When provided, it will be used to remove CSs whose
#'   minimum correlation among variables is smaller than
#'   \code{min_abs_corr}.
#'
#' @param coverage A number between 0 and 1 specifying desired
#'   coverage of each CS.
#'
#' @param min_abs_corr A "purity" threshold for the CS. Any CS that
#'   contains a pair of variables with correlation less than this
#'   threshold will be filtered out and not reported.
#'
#' @param dedup If \code{dedup = TRUE}, remove duplicate CSs.
#'
#' @param squared If \code{squared = TRUE}, report min, mean and
#' median of squared correlation instead of the absolute correlation.
#'
susie_get_cs = function (res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE) {
  if (!is.null(X) && !is.null(Xcorr))
    stop("Only one of X or Xcorr should be specified")
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr))
    stop("Xcorr matrix must be symmetric")
  null_index = 0
  include_idx = rep(TRUE,nrow(res$alpha))
  if (!is.null(res$null_index)) null_index = res$null_index
  # different from susieR
  # if (is.numeric(res$V)) include_idx = res$V > 1e-9

  # L x P binary matrix.
  status = in_CS(res$alpha,coverage)

  # L-list of CS positions.
  cs = lapply(1:nrow(status),function(i) which(status[i,]!=0))
  claimed_coverage = sapply(1:length(cs),
                            function (i) sum(res$alpha[i,][cs[[i]]]))
  include_idx = include_idx * (lapply(cs,length) > 0)

  # FIXME: see issue 21
  # https://github.com/stephenslab/susieR/issues/21
  if (dedup)
    include_idx = include_idx * (!duplicated(cs))
  include_idx = as.logical(include_idx)
  if (sum(include_idx) == 0)
    return(list(cs = NULL,
                coverage = NULL,
                requested_coverage = coverage))
  cs = cs[include_idx]
  claimed_coverage = claimed_coverage[include_idx]

  # Compute and filter by "purity".
  if (is.null(Xcorr) && is.null(X)) {
    names(cs) = paste0("L",which(include_idx))
    return(list(cs = cs,
                coverage = claimed_coverage,
                requested_coverage = coverage))
  } else {
    purity = data.frame(do.call(rbind,lapply(1:length(cs),function (i) {
              if (null_index > 0 && null_index %in% cs[[i]])
                c(-9,-9,-9)
              else
                get_purity(cs[[i]],X,Xcorr,squared)
             })))
    if (squared)
      colnames(purity) = c("min.sq.corr","mean.sq.corr","median.sq.corr")
    else
      colnames(purity) = c("min.abs.corr","mean.abs.corr","median.abs.corr")
    threshold = ifelse(squared,min_abs_corr^2,min_abs_corr)
    is_pure = which(purity[,1] >= threshold)
    if (length(is_pure) > 0) {
      cs = cs[is_pure]
      purity = purity[is_pure,]
      row_names = paste0("L",which(include_idx)[is_pure])
      names(cs) = row_names
      rownames(purity) = row_names

      # Re-order CS list and purity rows based on purity.
      ordering = order(purity[,1],decreasing = TRUE)
      return(list(cs       = cs[ordering],
                  purity   = purity[ordering,],
                  cs_index = which(include_idx)[is_pure[ordering]],
                  coverage = claimed_coverage[ordering],
                  requested_coverage=coverage))
    } else
      return(list(cs = NULL,coverage = NULL, requested_coverage = coverage))
  }
}

#' @title Get correlations between CS, using variable with maximum PIP from each CS
#'
#' @param X n by p matrix of values of the p variables (covariates) in
#'   n samples. When provided, correlation between variables will be
#'   computed and used to remove CSs whose minimum correlation among
#'   variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables
#'   (covariates). When provided, it will be used to remove CSs whose
#'   minimum correlation among variables is smaller than
#'   \code{min_abs_corr}.
#'
#' @keywords internal
#'
get_cs_correlation = function (res, X = NULL, Xcorr = NULL, max = FALSE) {
  if (is.null(res$sets$cs) || length(res$sets$cs) == 1) return(NA)
  if (!is.null(X) && !is.null(Xcorr))
    stop("Only one of X or Xcorr should be specified")
  if (is.null(Xcorr) && is.null(X))
    stop("One of X or Xcorr must be specified")
  if (!is.null(Xcorr) && !is_symmetric_matrix(Xcorr))
    stop("Xcorr matrix must be symmetric")
  # Get index for the best PIP per CS
  max_pip_idx = sapply(res$sets$cs, function(cs) cs[which.max(res$pip[cs])])
  if (is.null(Xcorr)) {
    X_sub = X[,max_pip_idx]
    cs_corr = muffled_corr(as.matrix(X_sub))
  } else {
    cs_corr = Xcorr[max_pip_idx, max_pip_idx]
  }
  if (max) {
    cs_corr = max(abs(cs_corr[upper.tri(cs_corr)]))
  }
  return(cs_corr)
}

#' @rdname susie_get_methods
#'
#' @param prune_by_cs Whether or not to ignore single effects not in
#'   a reported CS when calculating PIP.
#'
#' @param prior_tol Filter out effects having estimated prior variance
#'   smaller than this threshold.
#'
susie_get_pip = function (res, prune_by_cs = FALSE, prior_tol = 1e-9) {

  if (inherits(res,"susie")) {

    # Drop null weight columns.
    if (!is.null(res$null_index) && res$null_index > 0)
      res$alpha = res$alpha[,-res$null_index,drop=FALSE]

    # Drop the single-effects with estimated prior of zero.
    # different from susieR, will affect value of small PIPs
    # if (is.numeric(res$V))
    #   include_idx = which(res$V > prior_tol)
    # else
      include_idx = 1:nrow(res$alpha)

    # Only consider variables in reported CS.
    # This is not what we do in the SuSiE paper.
    # So by default prune_by_cs = FALSE means we do not run the
    # following code.
    if (!is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = intersect(include_idx,res$sets$cs_index)
    if (is.null(res$sets$cs_index) && prune_by_cs)
      include_idx = numeric(0)

    # now extract relevant rows from alpha matrix
    if (length(include_idx) > 0)
      res = res$alpha[include_idx,,drop = FALSE]
    else
      res = matrix(0,1,ncol(res$alpha))
  }

  return(as.vector(1 - apply(1 - res,2,prod)))
}

# Find how many variables in the CS.
# x is a probability vector.
#' @keywords internal
n_in_CS_x = function (x, coverage = 0.9)
  sum(cumsum(sort(x,decreasing = TRUE)) < coverage) + 1

# Return binary vector indicating if each point is in CS.
# x is a probability vector.
#' @keywords internal
in_CS_x = function (x, coverage = 0.9) {
  n = n_in_CS_x(x,coverage)
  o = order(x,decreasing = TRUE)
  result = rep(0,length(x))
  result[o[1:n]] = 1
  return(result)
}

# Returns an l-by-p binary matrix indicating which variables are in
# susie credible sets.
#' @keywords internal
in_CS = function (res, coverage = 0.9) {
  if (inherits(res,"susie"))
    res = res$alpha
  return(t(apply(res,1,function(x) in_CS_x(x,coverage))))
}

#' @keywords internal
n_in_CS = function(res, coverage = 0.9) {
  if (inherits(res,"susie"))
    res = res$alpha
  return(apply(res,1,function(x) n_in_CS_x(x,coverage)))
}

# Subsample and compute min, mean, median and max abs corr.
#
#' @importFrom stats median
get_purity = function(pos, X, Xcorr, squared = FALSE, n = 100) {
  if (length(pos) == 1)
    c(1,1,1)
  else {
    if (length(pos) > n)
      pos = sample(pos, n)
    if (is.null(Xcorr)) {
      X_sub = X[,pos]
      if (length(pos) > n) {

        # Remove identical columns.
        pos_rm = sapply(1:ncol(X_sub),
                       function(i) all(abs(X_sub[,i] - mean(X_sub[,i])) <
                                       .Machine$double.eps^0.5))
        if (length(pos_rm))
          X_sub = X_sub[,-pos_rm]
      }
      value = abs(muffled_corr(as.matrix(X_sub)))
    } else
      value = abs(Xcorr[pos,pos])
    if (squared)
      value = value^2
    return(c(min(value,na.rm = TRUE),
             mean(value,na.rm = TRUE),
             median(value,na.rm = TRUE)))
  }
}

# Correlation function with specified warning muffled.
#
#' @importFrom stats cor
muffled_corr = function (x)
  withCallingHandlers(cor(x),
                      warning = function(w) {
                        if (grepl("the standard deviation is zero",w$message))
                          invokeRestart("muffleWarning")
                      })

# cov2cor function with specified warning muffled.
#
#' @importFrom stats cov2cor
#' @keywords internal
muffled_cov2cor = function (x)
  withCallingHandlers(cov2cor(x),
    warning = function(w) {
      if (grepl("had 0 or NA entries; non-finite result is doubtful",
                w$message))
          invokeRestart("muffleWarning")
      })

# Check for symmetric matrix.
#' @keywords internal
is_symmetric_matrix = function (x) {
  res = isSymmetric(x)
  if (!res)
    res = isSymmetric(unname(x))
  return(res)
}

# Compute standard error for regression coef.
# S = (X'X)^-1 \Sigma
#' @keywords internal
calc_stderr = function (X, residuals)
  sqrt(diag(sum(residuals^2)/(nrow(X) - 2) * chol2inv(chol(t(X) %*% X))))

# Return residuals of Y after removing the linear effects of the susie
# model.
#
#' @importFrom stats coef
#' @keywords internal
get_R = function (X, Y, s)
  Y - X %*% coef(s)

# Slim the result of fitted susie model.
#' @keywords internal
susie_slim = function (res)
  list(alpha = res$alpha,niter = res$niter,V = res$V,sigma2 = res$sigma2)

# Prune single effects to given number L in susie model object.
#' @keywords internal
susie_prune_single_effects = function (s,L = 0,V = NULL,verbose = FALSE) {
  num_effects = nrow(s$alpha)
  if (L == 0) {

    # Filtering will be based on non-zero elements in s$V.
    if (!is.null(s$V))
      L = length(which(s$V > 0))
    else
      L = num_effects
  }
  if (L == num_effects) {
    s$sets = NULL
    return(s)
  }
  if (!is.null(s$sets$cs_index))
    effects_rank = c(s$sets$cs_index,setdiff(1:num_effects,s$sets$cs_index))
  else
    effects_rank = 1:num_effects
  if (verbose)
    warning(paste("Specified number of effects L =",L,
                  "does not match the number of effects",num_effects,
                  "in input SuSiE model. It will be",
                  ifelse(L < num_effects,"pruned","expanded"),"to have",L,
                  "effects."))
  if (L < num_effects) {
    for (n in c("alpha","mu","mu2","lbf_variable"))
      if (!is.null(s[[n]]))
        s[[n]] = s[[n]][effects_rank,][1:L,]
    for (n in c("KL","lbf","V"))
      if (!is.null(s[[n]]))
        s[[n]] = s[[n]][effects_rank][1:L]
  } else {
    s$alpha = rbind(s$alpha[effects_rank,],
                    matrix(1/ncol(s$alpha),L - num_effects,ncol(s$alpha)))
    for (n in c("mu","mu2","lbf_variable"))
      if (!is.null(s[[n]]))
        s[[n]] = rbind(s[[n]][effects_rank,],
                       matrix(0,L - num_effects,ncol(s[[n]])))
    for (n in c("KL", "lbf"))
      if (!is.null(s[[n]]))
        s[[n]] = c(s[[n]][effects_rank],rep(NA, L-num_effects))
    if (!is.null(V)) {
      if (length(V) > 1)
        V[1:num_effects] = s$V[effects_rank]
      else V = rep(V,L)
    }
    s$V = V
  }
  s$sets = NULL
  return(s)
}
#' @title Perform Univariate Linear Regression Separately for Columns of X
#' 
#' @description This function performs the univariate linear
#' regression y ~ x separately for each column x of X. Each regression
#' is implemented using \code{.lm.fit()}. The estimated effect size
#' and stardard error for each variable are outputted.
#' 
#' @param X n by p matrix of regressors.
#' 
#' @param y n-vector of response variables.
#' 
#' @param Z Optional n by k matrix of covariates to be included in all
#'   regresions. If Z is not \code{NULL}, the linear effects of
#'   covariates are removed from y first, and the resulting residuals
#'   are used in place of y.
#' 
#' @param center If \code{center = TRUE}, center X, y and Z.
#' 
#' @param scale If \code{scale = TRUE}, scale X, y and Z.
#' 
#' @param return_residuals Whether or not to output the residuals if Z
#' is not \code{NULL}.
#'
#' @return A list with two vectors containing the least-squares
#'   estimates of the coefficients (\code{betahat}) and their standard
#'   errors (\code{sebetahat}). Optionally, and only when a matrix of
#'   covariates \code{Z} is provided, a third vector \code{residuals}
#'   containing the residuals is returned.
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' X = scale(X,center = TRUE,scale = TRUE)
#' y = drop(X %*% beta + rnorm(n))
#' res = univariate_regression(X,y)
#' plot(res$betahat/res$sebetahat)
#' 
#' @importFrom stats lm
#' @importFrom stats .lm.fit
#' @importFrom stats coef
#' @importFrom stats summary.lm
#'
univariate_regression = function (X, y, Z = NULL, center = TRUE,
                                  scale = FALSE, return_residuals = FALSE) {
  y_na = which(is.na(y))
  if (length(y_na)) {
    X = X[-y_na,]
    y = y[-y_na]
  }
  if (center) {
    y = y - mean(y)
    X = scale(X,center = TRUE,scale = scale)
  } else 
    X = scale(X,center = FALSE,scale = scale)
  X[is.nan(X)] = 0
  if (!is.null(Z)) {
    if (center)
      Z = scale(Z,center = TRUE,scale = scale)
    y = .lm.fit(Z,y)$residuals
  }
  output = try(do.call(rbind,
                       lapply(1:ncol(X), function (i) {
                         g = .lm.fit(cbind(1,X[,i]),y)
                         return(c(coef(g)[2],calc_stderr(cbind(1,X[,i]),
                                                         g$residuals)[2]))
                       })),
               silent = TRUE)
  
  # Exception occurs, fall back to a safer but slower calculation.
  if (inherits(output,"try-error")) {
    output = matrix(0,ncol(X),2)
    for (i in 1:ncol(X)) {
      fit = summary(lm(y ~ X[,i]))$coef
      if (nrow(fit) == 2)
        output[i,] = as.vector(summary(lm(y ~ X[,i]))$coef[2,1:2])
      else
        output[i,] = c(0,0)
    }
  }
  if (return_residuals && !is.null(Z)) 
    return(list(betahat = output[,1],sebetahat = output[,2],residuals = y))
  else
    return(list(betahat = output[,1],sebetahat = output[,2]))
}

# Computes the z-scores (t-statistics) for association between Y and
# each column of X.
calc_z = function (X, Y, center = FALSE, scale = FALSE) {
  univariate_z = function(X,Y,center,scale) {
    out = univariate_regression(X,Y,center = center,scale = scale)
    return(out$betahat/out$sebetahat)
  }
  if (is.null(dim(Y)))
    return(univariate_z(X,Y,center,scale))
  else
    return(do.call(cbind,lapply(1:ncol(Y),
                                function(i) univariate_z(X,Y[,i],
                                                         center = center,
                                                         scale = scale))))
}
# @title update each effect once
# @param X an n by p matrix of regressor variables
# @param Y an n vector of response variable
# @param s a SuSiE fit
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
update_each_effect = function (X, Y, s, estimate_prior_variance = FALSE,
                               estimate_prior_method = "optim",
                               check_null_threshold) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  L = nrow(s$alpha)
  if (L > 0)
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$Xr = s$Xr - compute_Xb(X,s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      R = Y - s$Xr
      #different from susieR s$V[l]
      res = single_effect_regression(R,X,s$V[l,],s$sigma2,s$pi,
                                     estimate_prior_method,
                                     check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l,]      = res$V  #different from susieR
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$loglik +
        SER_posterior_e_loglik(X,R,s$sigma2,res$alpha * res$mu,
                               res$alpha * res$mu2)

      s$Xr = s$Xr + compute_Xb(X,s$alpha[l,] * s$mu[l,])
    }
  return(s)
}


# @title Update each effect once.
# @param R a p by p symmetric and positive semidefinite correlation matrix.
# @param z a p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param Sigma sigma2*R + lambda I
# @param estimate_prior_variance boolean indicating whether to estimate
#   prior variance
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
#
#' @importFrom Matrix diag
update_each_effect_rss = function (R, z, s_init, Sigma,
                                   estimate_prior_variance = FALSE,
                                   estimate_prior_method = "optim",
                                   check_null_threshold = 0) {

  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  if(L > 0)
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$Rz = s$Rz - R %*% (s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      r = z - s$Rz
      #different from susieR s$V[l]
      res = single_effect_regression_rss(as.vector(r),Sigma,s$V[l,],s$pi,
              estimate_prior_method,check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l,]      = res$V  #different from susieR
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_rss(R,Sigma,r,res$alpha * res$mu,
                                   res$alpha * res$mu2)
      s$Rz = s$Rz + R %*% (s$alpha[l,] * s$mu[l,])
    }
  s$Rz = unname(as.matrix(s$Rz))
  return(s)
}
# @title update each effect once
# @param XtX a p by p matrix, X'X
# @param Xty a p vector
# @param s_init a list with elements sigma2, V, alpha, mu, Xr
# @param estimate_prior_variance boolean indicating whether to
#   estimate prior variance
# @param estimate_prior_method The method used for estimating prior
#   variance, 'optim' or 'EM'.
# @param check_null_threshold float a threshold on the log scale to
#   compare likelihood between current estimate and zero the null
#
#' @importFrom Matrix diag
update_each_effect_ss = function (XtX, Xty, s_init,
                                  estimate_prior_variance = FALSE,
                                  estimate_prior_method = "optim",
                                  check_null_threshold = 0) {
  if (!estimate_prior_variance)
    estimate_prior_method = "none"

  # Repeat for each effect to update.
  s = s_init
  L = nrow(s$alpha)
  if (L > 0) {
    for (l in 1:L) {

      # Remove lth effect from fitted values.
      s$XtXr = s$XtXr - XtX %*% (s$alpha[l,] * s$mu[l,])

      # Compute residuals.
      XtR = Xty - s$XtXr
      #different from susieR s$V[l]
      res = single_effect_regression_ss(as.matrix(XtR),attr(XtX,"d"),s$V[l,],
              s$sigma2,s$pi,estimate_prior_method,check_null_threshold)

      # Update the variational estimate of the posterior mean.
      s$mu[l,]    = res$mu
      s$alpha[l,] = res$alpha
      s$mu2[l,]   = res$mu2
      s$V[l,]      = res$V  #different from susieR
      s$lbf[l]    = res$lbf_model
      s$lbf_variable[l,] = res$lbf
      s$KL[l]     = -res$lbf_model +
        SER_posterior_e_loglik_ss(attr(XtX,"d"),XtR,s$sigma2,
                                  res$alpha * res$mu,res$alpha * res$mu2)

      s$XtXr = s$XtXr + XtX %*% (s$alpha[l,] * s$mu[l,])
    }
  }
  s$XtXr = unname(as.matrix(s$XtXr))
  return(s)
}
