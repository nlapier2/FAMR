# a simplified re-implementation of cTWAS which does not require TWAS-style
#   prediction models as input and places fewer restrictions on input data.
# designed to be used for mendelian randomization.


# given susie results, update PIP (alpha) and
#   posterior mean (post_mu) and variance (tau2) of each variable
update_alpha_tau2 = function(res, n_vars, use_vb=F) {
  alpha = list()
  tau2 = list()
  post_mu = list()
  pips = as.numeric(res$pip)
  if(use_vb) {
    mu2 = as.numeric(colSums(t(res$s)))
  } else {
    mu2 = as.numeric(colSums(res$mu2))
  }
  post_means = as.numeric(colSums(res$alpha * res$mu))
  pips = sapply(pips, function(x) max(x, 1e-10))  # prevent underflow issues
  mu2 = sapply(mu2, function(x) max(x, 1e-10))
  post_means = sapply(post_means, function(x) max(x, 1e-10))
  prev = 1
  for(k in names(n_vars)) {
    end = prev + n_vars[[k]] - 1
    alpha[[k]] = pips[prev:end]
    tau2[[k]] = mu2[prev:end]
    post_mu[[k]] = post_means[prev:end]
    prev = end + 1
  }
  return(list('alpha' = alpha, 'tau2' = tau2, 'post_mu' = post_mu))
}


# given updated alpha and tau2, update priors for nonzero effect and variance
update_pi_sigma2 = function(alpha, tau2) {
  new_pi = list()
  new_sigma2 = list()
  for(k in names(alpha)) {
    ak = alpha[[k]]  # pips in class k
    sum_ak = sum(alpha[[k]])
    new_pi[[k]] = sum_ak / length(ak)
    new_sigma2[[k]] = sum(ak * tau2[[k]]) / sum_ak
  }
  return(list('pi' = new_pi, 'sigma2' = new_sigma2))
}


# run susie given prior parameters and return results
run_susie_with_priors = function(pred_vars=c(), y=c(), zscores=c(), R=c(),
                                 prior_pi=c(), prior_sigma2=c(), L=10, use_vb=F) {
  prior_pi_plain = as.numeric(unlist(prior_pi))
  prior_sigma2_plain = as.numeric(unlist(prior_sigma2))
  # prevent underflow issues
  prior_pi_plain = sapply(prior_pi_plain, function(x) max(x, 1e-10))
  prior_sigma2_plain = sapply(prior_sigma2_plain, function(x) max(x, 1e-10))
  prior_sigma2_plain = matrix(rep(prior_sigma2_plain, each = L), nrow=L)
  nw = max(0, 1 - sum(prior_pi_plain))  # null weight
  prior_pi_plain = prior_pi_plain / (1 - nw)

  # run susie
  z_plain = as.numeric(unlist(zscores))
  res = susie_rss(z_plain, R, L=L, # n=n,
                        prior_weights = prior_pi_plain,
                        prior_variance = prior_sigma2_plain,
                        null_weight = nw,
                        estimate_prior_variance=F)
  return(res)
}


# run E-M algorithm to learn per-variable and group priors
estimate_priors_EM = function(pred_vars=c(), y=c(), zscores=c(), R=c(),
                              L=1, n_iter=30, use_vb=F) {
  # initialize priors
  prior_pi = list()
  prior_sigma2 = list()
  n_vars = list()
  if(length(zscores) > 0) {  # sumstats mode
    M = dim(R)[1]  # total number of vars
    for(k in names(zscores)) {
      m_k = length(zscores[[k]])
      prior_pi[[k]] = rep(1 / M, m_k)
      prior_sigma2[[k]] = rep(50, m_k)
      n_vars[[k]] = m_k
    }
  } 

  for(iter in 1:n_iter) {
    # iteratively run susie with current priors, update per-variable params
    #   (alpha and tau2), and update group priors (pi and sigma2)
    message('Running prior estimation iteration ', iter)
    res = run_susie_with_priors(pred_vars=pred_vars, y=y, zscores=zscores, R=R,
                                     prior_pi=prior_pi, prior_sigma2=prior_sigma2,
                                     L=L, use_vb=use_vb)
    new_params = update_alpha_tau2(res, n_vars, use_vb=use_vb)
    new_group_priors = update_pi_sigma2(new_params$alpha, new_params$tau2)
    for(k in names(new_group_priors$pi)) {
      prior_pi[[k]] = rep(new_group_priors$pi[[k]], n_vars[[k]])
    }
    for(k in names(new_group_priors$pi)) {
      prior_sigma2[[k]] = rep(new_group_priors$sigma2[[k]], n_vars[[k]])
    }
    message('Current pi:')
    message(new_group_priors$pi)
    message('Current sigma2:')
    message(new_group_priors$sigma2)
  }
  
  # allow equal priors between all variables if n_iter=0
  if(n_iter == 0) {
    for(k in names(prior_pi)) {
      prior_pi[[k]] = prior_pi[[k]] * 0 + 1
      prior_sigma2[[k]] = prior_sigma2[[k]] * 0 + 50
    }
  }
  return(list('pi' = prior_pi, 'sigma2' = prior_sigma2, 'n_vars' = n_vars))
}


# run final susie(-rss) regression to determine PIPs for each variable
get_results = function(pred_vars=c(), y=c(), zscores=c(), R=c(), priors=c(), L=10, use_vb=F) {
  susieres = run_susie_with_priors(pred_vars=pred_vars, y=y, zscores=zscores, R=R,
                                   prior_pi=priors$pi, prior_sigma2=priors$sigma2,
                                   L=L, use_vb=use_vb)
  res_df = list('pips' = list(), 'posterior_var' = list(), 'posterior_mean' = list())
  res_df$susieres = susieres
  res_df$zscores = zscores
  res_df$R = R
  res_df$x_pred = if(length(R) == 0) pred_vars$exposures else c()
  res_df$priors = priors
  res_df$L = L
  final_params  = update_alpha_tau2(susieres, priors$n_vars, use_vb=use_vb)
  for(k in names(final_params$alpha)) {
    res_df$pips[[k]] = final_params$alpha[[k]]
    res_df$posterior_var[[k]] = final_params$tau2[[k]]
    res_df$posterior_mean[[k]] = final_params$post_mu[[k]]
    # retain variable names
    varnames = if(length(R) > 0) names(zscores[[k]]) else colnames(pred_vars[[k]])
    names(res_df$pips[[k]]) = varnames
    names(res_df$posterior_var[[k]]) = varnames
    names(res_df$posterior_mean[[k]]) = varnames
  }
  res_df$sorted_pips = sort(unlist(res_df$pips), decreasing=T)
  return(res_df)
}
