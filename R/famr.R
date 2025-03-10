# Main user-facing FAMR script

# main function for users to run.
# a wrapper function that runs all parts of the FAMR-Susie method.
run_famr_susie = function(in_dir, ld_dir, y_gwas_file='NONE', 
                          fa_method='gfa', num_samples=10000, n_iter=30, 
                          susieL=30, prune_thresh=0.5, fa_prune_thresh=0.1, 
                          final_prune_thresh=0.1, annihilate=FALSE,
                          idcol=1, betacol=2, secol=3, header=TRUE) {
  
  # read outcome gwas, if appropriate
  y_gwas = read_y_gwas(y_gwas_file, idcol, betacol, secol, header)
  
  # read data from input files, add sumstats for outcome
  sumstats = gather_sumstats(in_dir, ld_dir, prune_thresh=prune_thresh,
                             fa_prune_thresh=fa_prune_thresh, y_gwas=y_gwas)
  
  # generate factors
  factor_ss = generate_factors(fa_method, sumstats$ss_for_fa, num_samples,
                               full_ss=sumstats)
  
  # project out factors from exposures and outcome if requested by user
  if(annihilate) {
    sumstats = annihilate_factors(sumstats, factor_ss)
  }
  
  # precompute Z-scores and LD between phenotypes required for susie_rss
  dat = learn_wts_precomp_merge(sumstats, factor_ss, N=num_samples, 
                                prune_thresh=final_prune_thresh)
  
  # run ctwas-style framework and return
  n_expo = ncol(dat$sumstats$betas) - factor_ss$n_factors
  famr_res = run_modified_ctwas(dat, susieL, n_iter, n_expo, 
                                num_samples, annih=annihilate)
  return(famr_res)
}


# read outcome gwas sumstats if provided
read_y_gwas = function(y_gwas_file, idcol, betacol, secol, header=T) {
  message('Reading outcome GWAS...')
  if(y_gwas_file == 'NONE')  return(c())
  if(endsWith(y_gwas_file, '.vcf.gz') || endsWith(y_gwas_file, '.vcf')) {
    gwas = read_vcf(y_gwas_file)
    gwas = gwas %>% dplyr::select(ID, BETA, SE)  # remove unnecessary fields
  } else {
    gwas = fread(y_gwas_file, sep='\t', header=header,
                 select=c(idcol, betacol, secol))
  }
  colnames(gwas)  = c('names_y', 'betas_y', 'stderrs_y')
  gwas$zscores_y = gwas$betas_y / gwas$stderrs_y
  return(gwas)
}


# project out factors from exposures and outcome if requested by user
annihilate_factors = function(sumstats, factor_ss) {
  if(factor_ss$n_factors == 0)  return(sumstats)
  sumstats$orig_betas_y = sumstats$betas_y
  sumstats$orig_zscores_y = sumstats$zscores_y
  xpredfss = as.matrix(lm(sumstats$betas ~ factor_ss$betas)$fitted.values)
  ypredfss = as.numeric(lm(sumstats$betas_y ~ factor_ss$betas)$fitted.values)
  sumstats$betas = sumstats$betas - xpredfss
  sumstats$betas_y = sumstats$betas_y - ypredfss
  sumstats$Z = sumstats$betas / sumstats$stderrs
  sumstats$zscores_y = sumstats$betas_y / sumstats$stderrs_y
  sumstats$annih_y = ypredfss
  return(sumstats)
}


# learn regularization parameters for susie-rss weight learning
learn_lambda = function(sumstats, ld, oracle_mode=F) {
  lambda = rep(0, ncol(as.matrix(sumstats$betas)))
  if(oracle_mode)  return(lambda)
  z_gx = as.matrix(sumstats$betas / sumstats$stderrs)
  if(ncol(z_gx) == 1)  z_gx = t(z_gx)
  # downsample large LD matrices to save time on eigendecomposition
  if(nrow(ld) > 2000) {
    samp_idx = sort(sample(1:nrow(ld), 2000, replace=F))
    ld = ld[samp_idx, samp_idx]
    z_gx = z_gx[samp_idx, ]
  }
  # adapted from: https://github.com/stephenslab/susieR/blob/master/R/susie_rss_lambda.R
  R_adj = susieR:::set_R_attributes(as.matrix(ld), 1e-8)
  colspace = which(attr(R_adj,"eigen")$values > 0)
  if (length(colspace) == nrow(z_gx)) {
    lambda = rep(0, ncol(z_gx)) 
  } else {
    znull = as.matrix(apply(z_gx, 2, function(x) crossprod(attr(R_adj,"eigen")$vectors[,-colspace], x)))
    if(ncol(znull) == 1)  znull = t(znull)
    lambda = apply(znull, 2, function(x) sum(x^2) / length(x))
    lambda = sapply(lambda, function(x) max(0, min(1, x)))
  }
  return(lambda)
}


# learn weights of SNP effects on exposures
learn_weights = function(sumstats, ld, lambda, nome=F, N=10000, L=10, 
                         prune_min=0.1, prune_max=0.5) {
  z_gx = as.matrix(sumstats$betas / sumstats$stderrs)
  if(ncol(z_gx) == 1)  z_gx = t(z_gx)
  ld = as.matrix(ld)
  weights = t(as.matrix(sumstats$betas))
  if(nome)  return(weights)

  # impute lambda for factors as weighted linear combination of loaded exposures
  if(length(lambda) < ncol(z_gx)) {
    if(nrow(z_gx) == 1) {
      newvals = rep(mean(lambda), ncol(z_gx) - length(lambda))
    } else {
      l_lam = length(lambda)
      tmp = t(abs(cor(z_gx)[1:l_lam,(l_lam+1):ncol(z_gx)]))
      newvals = as.numeric(tmp %*% lambda) / rowSums(tmp)
    }
    lambda = c(lambda, newvals)
  }
  
  # adaptively set number of possible causal SNPs based on lambda
  Lvec = ceiling(10 * (1 - lambda) + 1e-10)
  Lvec = sapply(Lvec, function(x) max(3, min(x, L)))
  
  # LD prune all SNPs based on prune_max (the most liberal prune threshold)
  if(length(sumstats) > 0) {
    idx = ld_prune_famr(sumstats, ld, prune_max)
    sumstats = subset_sumstats(sumstats, idx)
    ld = ld[idx, idx, drop=F]
    cnames = colnames(ld)
    colnames(ld) = NULL
    R_adj = susieR:::set_R_attributes(ld, 1e-8)
    z_gx = as.matrix(sumstats$betas / sumstats$stderrs)
    if(ncol(z_gx) == 1)  z_gx = t(z_gx)
  }

  # perform per-exposure "pseudo-prune" adaptively based on lambda
  # e.g. if lambda=1 the LD is diagonal, so strict prune
  prunes = 1 - lambda
  prunes = sapply(prunes, function(x) max(prune_min, min(prune_max, x)))
  for(i in 1:length(prunes)) {
    if(prunes[i] < prune_max) {
      idx = ld_prune_famr(sumstats, ld, prunes[i])
      z_gx[-idx, i] = 0  # "pseudo-prune" by setting Z=0 for those SNPs for this trait
    }
  }
  
  # learn weights with susie for each exposure
  susie_rss_res = lapply(1:ncol(z_gx), function(x) susieR::susie_rss(z_gx[,x], 
                         (1 - lambda[x]) * ld + lambda[x] * diag(nrow(ld)), N, L=Lvec[x]))
  weights = t(sapply(susie_rss_res, function(x) colSums(x$alpha * x$mu)))
  # set weights to zero for SNPs not in credible sets
  all_cs_snps = sapply(susie_rss_res, function(x) sort(as.numeric(unlist((x)$sets$cs))))
  for(i in 1:nrow(weights)) {
    if(length(all_cs_snps[[i]]) == 0) {
      weights[i, ] = weights[i, ] * 0
    } else {
      weights[i, -all_cs_snps[[i]]] <- 0
    }
  }
  if(ncol(t(weights)) == 1) weights = t(weights)
  return(t(weights))
}


# precompute, for each independent locus, the numerator and denominator of the
#   formula for imputing the exposure-outcome z-score (zxy)
precompute_zxy = function(weights, z_gy, ld, orig_z_gy=c(), n_factors=0) {
  z_gy = as.matrix(z_gy)
  if(n_factors == 0 || length(orig_z_gy) == 0) {
    numerator = t(t(weights) %*% z_gy)
  } else {  
    # factor z_scores computed using annihilated (factor-associated) part of y
    orig_z_gy = as.matrix(orig_z_gy)
    wts_expo = weights[, 1:(ncol(weights) - n_factors)]
    wts_fac = weights[, (ncol(weights) - n_factors + 1):ncol(weights)]
    numerator = c(t(t(wts_expo) %*% z_gy), t(t(wts_fac) %*% orig_z_gy))
  }
  Rgx = ld %*% weights
  wRw = t(weights) %*% ld %*% weights
  return(list('numerator' = numerator, 'wRw' = wRw, 'Rgx' = Rgx))
}


# gather summary statistics across all loci.
# can merge in gwas for the outcome (y) if provided separately (in y_gwas)
gather_sumstats = function(in_dir, ld_dir, y_gwas=c(), oracle_mode=F, 
                           prune_thresh=1.0, fa_prune_thresh=1.0) {
  message('Gathering summary statistics data...')
  fa_prune_thresh = min(fa_prune_thresh, prune_thresh)
  dat = list('lambda' = c(), 'sumstats' = list(), 'ss_for_fa' = list(), 'ld' = list())
  ss_fnames = list.files(path = in_dir, pattern = '*.rds', full.names = T)
  ld_fnames = list.files(path = ld_dir, pattern = '*.rds', full.names = T)

  # loop through loci, LD pruning SNPs, possibly intersecting with y_gwas, and
  #   appending the SNPs that pass to the overall summary statistics
  for(f_idx in 1:length(ss_fnames)) {
    message('Gathering from locus ', f_idx)
    sumstats = readRDS(ss_fnames[f_idx])
    if(!is.null(nrow(sumstats$pos)))  sumstats$pos = sumstats$pos$ID
    ld = Matrix(as.matrix(readRDS(ld_fnames[f_idx])))

    # remove SNPs not in outcome gwas, if appropriate
    if(length(y_gwas) > 0) {
      with_y = merge_outcome(sumstats, ld, y_gwas)
      if(length(with_y$ld) == 0)  next
      sumstats = with_y$sumstats
      ld = with_y$ld
    }
    
    if(length(ld) == 0)  next  # if no SNPs remaining, iterate loop

    # learn "lambda" regularization parameter
    sumstats$lambda = learn_lambda(sumstats, ld, oracle_mode = oracle_mode)
    sumstats$weights = sumstats$betas  # weights will be learned later
    
    # perform LD pruning
    dat = prune_and_merge(dat, sumstats, ld, prune_thresh, fa_prune_thresh, 
                          f_idx, oracle_mode)
  }
  dat$betas = as.matrix(dat$betas)
  dat$stderrs = as.matrix(dat$stderrs)
  dat$weights = as.matrix(dat$weights)
  dat$ss_for_fa$betas = as.matrix(dat$ss_for_fa$betas)
  dat$ss_for_fa$stderrs = as.matrix(dat$ss_for_fa$stderrs)
  dat$ss_for_fa$weights = as.matrix(dat$ss_for_fa$weights)
  return(dat)
}


# equivalent of gather_sumstats, but takes data directly as input rather than
#   reading it in from files
gather_sumstats_from_dat = function(all_ss, all_ld, y_gwas=c(), oracle_mode=F, 
                                    prune_thresh=1.0, fa_prune_thresh=1.0) {
  message('Gathering summary statistics data...')
  fa_prune_thresh = min(fa_prune_thresh, prune_thresh)
  dat = list('lambda' = c(), 'sumstats' = list(), 'ss_for_fa' = list(), 'ld' = list())
  
  # loop through loci, LD pruning SNPs, possibly intersecting with y_gwas, and
  #   appending the SNPs that pass to the overall summary statistics
  for(f_idx in 1:max(all_ss$locus_idx)) {
    message('Gathering from locus ', f_idx)
    idx = all_ss$locus_idx == f_idx
    sumstats = subset_sumstats(all_ss, idx)
    ld = all_ld[[f_idx]]
    
    # remove SNPs not in outcome gwas, if appropriate
    if(length(y_gwas) > 0) {
      with_y = merge_outcome(sumstats, ld, y_gwas)
      if(length(with_y$ld) == 0)  next
      sumstats = with_y$sumstats
      ld = with_y$ld
    }
    
    if(length(ld) == 0)  next  # if no SNPs remaining, iterate loop
    
    # learn "lambda" regularization parameter
    sumstats$lambda = learn_lambda(sumstats, ld, oracle_mode = oracle_mode)
    sumstats$weights = sumstats$betas  # weights will be learned later
    
    # perform LD pruning
    dat = prune_and_merge(dat, sumstats, ld, prune_thresh, fa_prune_thresh, 
                          f_idx, oracle_mode)
  }
  dat$betas = as.matrix(dat$betas)
  dat$stderrs = as.matrix(dat$stderrs)
  dat$weights = as.matrix(dat$weights)
  dat$ss_for_fa$betas = as.matrix(dat$ss_for_fa$betas)
  dat$ss_for_fa$stderrs = as.matrix(dat$ss_for_fa$stderrs)
  dat$ss_for_fa$weights = as.matrix(dat$ss_for_fa$weights)
  return(dat)
}



# generate factors for FAMR if requested by user
generate_factors = function(fa_method, sumstats, N=10000, given_factors='NONE',
                            full_ss=c()) {
  factor_sumstats = list('n_factors' = 0)
  if(is.na(fa_method) || tolower(fa_method) == 'none')  return(factor_sumstats)
  message('Generating factors with: ', fa_method)
  # generate appropriate factors
  if(given_factors != 'NONE' && file.exists(given_factors)) {  # use pre-specified factors
    # for given factor SNPs set to true value, for others set to 0
    factors = readRDS(given_factors)
    causal_idx = match(factors$pos, sumstats$pos)
    locus_wts = factors$betas[!is.na(causal_idx),]  # only causal SNPs in this locus
    causal_idx = causal_idx[!is.na(causal_idx)]
    factor_wts = matrix(0, nrow(sumstats$betas), ncol(factors$betas))
    factor_wts[causal_idx,] = locus_wts
  } else if(fa_method == 'gfa') {
    gfares = run_gfa_full(sumstats$betas, sumstats$stderrs, N)
    if(length(gfares$L_hat) > 0 && length(full_ss) > 0
       && nrow(full_ss$betas) > nrow(sumstats$betas)) {
      # impute GFA L_hat for full set of SNPs, if applicable
      if(is.null(dim(gfares$F_hat)))  gfares$F_hat = as.matrix(gfares$F_hat)
      factor_wts = as.matrix(GFA:::loadings_gls(full_ss$betas, full_ss$stderrs,
                                      cor(full_ss$betas), gfares$F_hat)$L)
    } else {
      factor_wts = gfares$L_hat
    }
    factor_sumstats$gfares = gfares
  } else {
    return(factor_sumstats)
  }
  if(length(factor_wts) == 0 || max(abs(factor_wts)) == 0) {
    # return no factors if method returns empty or all-zero factors
    message('Number of factors detected: 0')
    return(list('n_factors' = 0))
  } else if(is.null(ncol(factor_wts))) {  # convert vector to matrix
    factor_wts = as.matrix(factor_wts)
  }
  colnames(factor_wts) = paste0('factor', 1:ncol(factor_wts))
  factor_sumstats$betas = factor_wts
  factor_sumstats$weights = factor_wts

  # also provide stderrs, zscores, snp id, etc
  factor_sumstats$stderrs = factor_sumstats$betas * 0 + (1 / sqrt(N))
  factor_sumstats$Z = factor_sumstats$betas / factor_sumstats$stderrs
  if(fa_method == 'gfa' && length(full_ss) > 0 && nrow(full_ss$betas) > nrow(sumstats$betas)) {
    factor_sumstats$pos = full_ss$pos
    factor_sumstats$corrs = cor(factor_sumstats$betas, full_ss$betas)
  } else {
    factor_sumstats$pos = sumstats$pos
    factor_sumstats$corrs = cor(factor_sumstats$betas, sumstats$betas)
  }
  factor_sumstats$n_factors = max(0, ncol(factor_sumstats$betas))
  message('Number of factors detected: ', factor_sumstats$n_factors)
  return(factor_sumstats)
}


# impute Z and R values needed for FAMR, including for factors
learn_wts_precomp_merge = function(all_sumstats, factor_ss, N=10000, 
                                   nome=F, prune_thresh=1.0) {
  message('Computing Z-scores and LD for exposures and factors...')
  res = list('sumstats' = list(), 'ld' = c(),
             'numerator' = 0, 'wRw' = 0, 'Rgx' = c(), 'factor_corrs' = c())
  lam_idx = 0
  for(f_idx in 1:max(all_sumstats$locus_idx)) {
    idx = all_sumstats$locus_idx == f_idx
    sumstats = subset_sumstats(all_sumstats, idx)
    if(nrow(sumstats$betas) == 0)  next
    ld = all_sumstats$ld[[f_idx]]
    lam_idx = lam_idx + 1
    lambda = all_sumstats$lambda[lam_idx, ]
    if(factor_ss$n_factors > 0) {  # merge in factor summary statistics if given
      factor_ss_sub = subset_sumstats(factor_ss, idx)
      sumstats = merge_sumstats(sumstats, factor_ss_sub, by_row=F)
    }
    
    # learn SNP weights on exposures and factors, precompute components of Z and R for this locus 
    sumstats$weights = learn_weights(sumstats, ld, lambda, N=N, nome=nome,
                                     prune_max=1.0, prune_min=0.1)
    precomp = precompute_zxy(sumstats$weights, sumstats$zscores_y, ld, 
                             orig_z_gy=sumstats$orig_zscores_y, 
                             n_factors=factor_ss$n_factors)
    res$numerator = res$numerator + precomp$numerator
    res$wRw = res$wRw + precomp$wRw
    
    # filter SNPs for final regression / variable selection procedure, 
    #   then merge into overall sumstats
    p_idx = ld_prune_famr(sumstats, ld, prune_thresh)
    message(length(p_idx), ' SNPs passed filters in this locus.')
    if(length(p_idx != 0)) {
      sumstats = subset_sumstats(sumstats, p_idx)
      res$sumstats = merge_sumstats(res$sumstats, sumstats)
      ld = ld[p_idx, p_idx, drop=F]
      res$ld = Matrix(merge_ld(res$ld, ld))
      res$Rgx = rbind(res$Rgx, precomp$Rgx[p_idx, ])
    }
  }
  res$wRw[res$wRw == 0] = 1e-10  # prevent errors from dividing by 0 later
  
  if(factor_ss$n_factors > 0)  res$factor_corrs = factor_ss$corrs
  message(nrow(res$sumstats$betas), ' total variants passed FAMR filters.')
  return(res)
}


# run zscore and ld imputation, EM estimation, and regression on all data
run_modified_ctwas = function(dat, L, n_iter, n_expo, n_samp, annih) {
  message('Running FAMR...')
  # impute exposure z-scores on outcome and ld with SNPs using precomputed data
  R = impute_exposure_corrs(dat$ld, t(as.matrix(dat$sumstats$weights)),
                                  wRw = dat$wRw, Rgx = t(dat$Rgx))
  
  # run a separate regression with factors only, if factors were inferred
  factor_susie_res = list()
  z_factors = c()
  if(length(dat$factor_corrs) > 0) {
    KplusJ = length(as.numeric(dat$numerator))
    f_idx = (n_expo+1):KplusJ  # indices of factors
    R_factors = R[f_idx, f_idx, drop=F]
    z_factors = dat$numerator[f_idx] / sqrt(diag(dat$wRw))[f_idx]
    factor_susie_res = susieR::susie_rss(z_factors, as.matrix(R_factors), n_samp, L=L)
    if(annih) {  # annihilator mode: separate exposures from factors
      non_f_idx = c(1:n_expo, (KplusJ+1):nrow(R))  # non-factor indices
      R = R[non_f_idx, non_f_idx, drop=F]
      dat$numerator = dat$numerator[1:n_expo]
      dat$wRw = dat$wRw[1:n_expo, 1:n_expo]
    }
  }
  
  R = Matrix(R, sparse=TRUE)
  expo_effects = dat$numerator / sqrt(diag(dat$wRw))
  effects = list('exposures' = expo_effects, 'snps' = dat$sumstats$zscores_y)
  names(effects$exposures) = names(dat$sumstats$Z)[1:length(effects$exposures)]
  names(effects$snps) = dat$sumstats$pos
  # calculate priors w/ EM algorithm
  priors = estimate_priors_EM(zscores=effects, R=R, L=L, n_iter=n_iter)
  M_expo = length(expo_effects)
  if(M_expo > n_expo && n_iter > 0)  priors$pi$exposures[(n_expo+1):M_expo] = priors$pi$exposures[1] * 1000

  # run final variable selection given priors
  famr_res = get_results(zscores=effects, R=R, priors=priors, L=L)
  
  # if appropriate, set variables for results of separate factor variable selection
  if(length(dat$factor_corrs) > 0) {
    famr_res$factor_susie_res = factor_susie_res
    famr_res$zscores$factors = z_factors
    famr_res$pips$factors = factor_susie_res$pip
    famr_res$pip$factors = factor_susie_res$pip
    names(famr_res$zscores$factors) = rownames(dat$factor_corrs)
    names(famr_res$pips$factors) = rownames(dat$factor_corrs)
    names(famr_res$pip$factors) = rownames(dat$factor_corrs)
    famr_res$factor_corrs = dat$factor_corrs
  }
  return(famr_res)
}
