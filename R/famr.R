# Main user-facing FAMR script

### COMMAND LINE OPTION PARSING

if (sys.nframe() == 0) {
  option_list = list(
    make_option(c("-i", "--in_dir"), type="character", default='./',
                help="Directory of summary statistics files to use", metavar="character"),
    make_option(c("--ld_dir"), type="character", default='./',
                help="Directory of LD matrices for loci", metavar="character"),
    make_option(c("--precomp_dat"), type="character", default='NONE', metavar="character",
                help="Pass in a precomputed dataframe with sumstats, ld, factors"),
    make_option(c("--y_gwas_file"), type="character", default='NONE',
                help="File with GWAS results for outcome phenotype.", metavar="character"),
    make_option(c("--fa_method"), type="character", default='gfa',
                help="Factor analysis method to use", metavar="character"),
    make_option(c("-N", "--num_samples"), type="numeric", default=10000,
                help="Sample size of exposure data.", metavar="numeric"),
    make_option(c("--susieL"), type="numeric", default=30,
                help="SuSiE 'L' setting (max number of effect variables)", metavar="numeric"),
    make_option(c("--n_iter"), type="numeric", default=30,
                help="Number of iterations to run FAMR EM loop", metavar="numeric"),
    make_option(c("--idcol"), type="numeric", default=1,
                help="Column number with SNP IDs in y_gwas_file.", metavar="numeric"),
    make_option(c("--betacol"), type="numeric", default=2,
                help="Column number with SNP betas in y_gwas_file.", metavar="numeric"),
    make_option(c("--secol"), type="numeric", default=3,
                help="Column number with SNP std errors in y_gwas_file.", metavar="numeric"),
    make_option(c("--no_header"), type="logical", action="store_true", default=FALSE,
                help="Use if sumstat file has no header.", metavar="logical"),
    make_option(c("--out"), type="character", default='res_',
                help="Prefix of output and temporary files", metavar="character"),
    make_option(c("--prune_thresh"), type="numeric", default=0.5, metavar="numeric",
                help="Threshold for LD pruning (default = 0.5)"),
    make_option(c("--fa_prune_thresh"), type="numeric", default=0.1, metavar="numeric",
                help="Threshold for LD pruning for factor analysis (default = 0.1)"),
    make_option(c("--final_prune_thresh"), type="numeric", default=0.1, metavar="numeric",
                help="Threshold for LD pruning for final regressions (default = 0.1)"),
    make_option(c("--nome"), type="logical", action="store_true", default=FALSE,
                help="Assume no measurement error in SNP effect estimates"),
    make_option(c("--annihilate_factors"), type="logical", action="store_true", default=FALSE,
                help="Project factors out of X & Y instead of including as covariates"),
    make_option(c("--given_factors"), type="character", default='NONE',
                help="Optional file with prespecified factors", metavar="character")
  );

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  out = opt$out
}


# read outcome gwas sumstats if provided
read_y_gwas = function(y_gwas_file, idcol, betacol, secol, header=T) {
  print('Reading outcome GWAS...')
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


# merge outcome (y) summary statistics into sumstats object,
#   keeping only the snps available for both exposures and outcome
merge_outcome = function(sumstats, ld, y_gwas) {
  if(is.null(nrow(sumstats$pos))) {  # if pos is only a vector of IDs
    y_gwas = y_gwas[names_y %in% sumstats$pos, ]
    idx = which(sumstats$pos %in% y_gwas$names_y)
  } else {  # if pos is a list/df with an ID row
    y_gwas = y_gwas[names_y %in% sumstats$pos$ID, ]
    idx = which(sumstats$pos$ID %in% y_gwas$names_y)
  }
  if(length(idx) == 0)  return(list('sumstats' = c(), 'ld' = c()))
  sumstats = subset_sumstats(sumstats, idx)
  sumstats$zscores_y = y_gwas$zscores_y
  sumstats$betas_y = y_gwas$betas_y
  sumstats$stderrs_y = y_gwas$stderrs_y
  ld = ld[idx, idx, drop=F]
  return(list('sumstats' = sumstats, 'ld' = ld))
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


# simple wrapper function to write results (for use with simulation script only)
write_famr_res = function(out, famr_res, K, fa_method, for_sim=F) {
  # construct method_name which identifies the FA method used
  method_name = ifelse(toupper(fa_method) != 'NONE' && !is.na(fa_method),
                       paste0('famr_', fa_method), 'famr')

  # save full results in RDS file
  saveRDS(famr_res, paste0(out, method_name, '.rds'))

  # write simpler results text file for evaluating simulations, if requested
  if(for_sim) {
    expo_pips = t(famr_res$pips$exposures[1:K])  # exclude factor pips
    expo_post_means = t(famr_res$posterior_mean$exposures[1:K])
    expo_post_stderr = t(sqrt(famr_res$posterior_var$exposures[1:K]))
    res_vars = c(expo_post_means, expo_post_stderr, expo_pips)
    write.table(t(res_vars), paste0(out, method_name, '.txt'),
                row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}


# helper function for LD pruning SNPs
ld_prune_famr = function(sumstats, ld, prune_thresh) {
  if(length(ld) == 1)  return(c(1))  # if only one SNP, no pruning
  # Sort betas in descending order, reshuffle LD matrix
  max_abs_betas <- apply(abs(sumstats$betas), 1, max)
  sorted_ss <- data.frame(max_abs_betas = max_abs_betas, ord = 1:nrow(sumstats$betas))
  sorted_ss <- sorted_ss[order(-sorted_ss$max_abs_betas), ]
  ld <- abs(ld[sorted_ss$ord, sorted_ss$ord])

  # Prune SNPs based on LD threshold in the upper triangle of sorted LD matrix
  to_prune <- rep(FALSE, ncol(ld))
  for (idx in 2:ncol(ld)) {
    if (to_prune[idx]) next  # Skip already pruned SNPs
    vals <- ld[1:(idx-1), idx]  # upper triangle
    if (any(vals[!to_prune[1:(idx-1)]] > prune_thresh)) {
      to_prune[idx] <- TRUE
    }
  }
  keep_idx <- sorted_ss$ord[!to_prune]  # original order
  return(sort(keep_idx))
}


# wrapper function to LD prune and merge summary statistics
prune_and_merge = function(dat, sumstats, ld, prune_thresh, fa_prune_thresh, 
                           f_idx, oracle_mode) {
  if(!oracle_mode) {  # don't prune in oracle mode since all SNPs are truly causal
    idx = ld_prune_famr(sumstats, ld, prune_thresh)
    sumstats = subset_sumstats(sumstats, idx)
    ld = ld[idx, idx, drop=F]
    if(fa_prune_thresh < prune_thresh) {
      # if applicable, keep a more restricted set of SNPs for FA methods to use
      fa_idx = ld_prune_famr(sumstats, ld, fa_prune_thresh)
      ss_for_fa = subset_sumstats(sumstats, fa_idx)
    } else {
      ss_for_fa = sumstats
    }
  }
  # keep these indices and merge locus sumstats into overall sumstats
  sumstats$locus_idx = rep(f_idx, nrow(sumstats$betas))
  ss_for_fa$locus_idx = rep(f_idx, nrow(ss_for_fa$betas))
  dat = merge_sumstats(dat, sumstats)
  dat$ss_for_fa = merge_sumstats(dat$ss_for_fa, ss_for_fa)
  dat$ld[[f_idx]] = ld
  return(dat)
}


# gather summary statistics across all loci.
# can merge in gwas for the outcome (y) if provided separately (in y_gwas)
gather_sumstats = function(in_dir, ld_dir, y_gwas=c(), oracle_mode=F, 
                           prune_thresh=1.0, fa_prune_thresh=1.0) {
  print('Gathering summary statistics data...')
  fa_prune_thresh = min(fa_prune_thresh, prune_thresh)
  dat = list('lambda' = c(), 'sumstats' = list(), 'ss_for_fa' = list(), 'ld' = list())
  ss_fnames = list.files(path = in_dir, pattern = '*.rds', full.names = T)
  ld_fnames = list.files(path = ld_dir, pattern = '*.rds', full.names = T)

  # loop through loci, LD pruning SNPs, possibly intersecting with y_gwas, and
  #   appending the SNPs that pass to the overall summary statistics
  for(f_idx in 1:length(ss_fnames)) {
    print(paste0('Gathering from locus ', f_idx))
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
  print('Gathering summary statistics data...')
  fa_prune_thresh = min(fa_prune_thresh, prune_thresh)
  dat = list('lambda' = c(), 'sumstats' = list(), 'ss_for_fa' = list(), 'ld' = list())
  
  # loop through loci, LD pruning SNPs, possibly intersecting with y_gwas, and
  #   appending the SNPs that pass to the overall summary statistics
  for(f_idx in 1:max(all_ss$locus_idx)) {
    print(paste0('Gathering from locus ', f_idx))
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
  print(paste0('Generating factors with: ', fa_method))
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
    print('Number of factors detected: 0')
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
  print(paste0('Number of factors detected: ', factor_sumstats$n_factors))
  return(factor_sumstats)
}


# impute Z and R values needed for FAMR, including for factors
learn_wts_precomp_merge = function(all_sumstats, factor_ss, N=10000, 
                                   nome=F, prune_thresh=1.0) {
  print('Computing Z-scores and LD for exposures and factors...')
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
    print(paste0(length(p_idx), ' SNPs passed filters in this locus.'))
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
  print(paste0(nrow(res$sumstats$betas), ' total variants passed FAMR filters.'))
  return(res)
}


# run zscore and ld imputation, EM estimation, and regression on all data
run_famr_rss_all = function(dat, L, n_iter, n_expo, n_samp, annih) {
  print('Running FAMR...')
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


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # read outcome gwas, if appropriate
  y_gwas = read_y_gwas(opt$y_gwas_file, opt$idcol, opt$betacol, opt$secol, !opt$no_header)

  # read data from input files, add sumstats for outcome
  sumstats = gather_sumstats(opt$in_dir, opt$ld_dir, prune_thresh=opt$prune_thresh,
                             fa_prune_thresh=opt$fa_prune_thresh, y_gwas=y_gwas)

  # generate factors
  factor_ss = generate_factors(opt$fa_method, sumstats$ss_for_fa, opt$num_samples,
                               opt$given_factors, full_ss=sumstats)

  # project out factors from exposures and outcome if requested by user
  if(opt$annihilate_factors) {
    sumstats = annihilate_factors(sumstats, factor_ss)
  }

  # precompute Z-scores and LD between phenotypes required for susie_rss
  dat = learn_wts_precomp_merge(sumstats, factor_ss, N=opt$num_samples, 
                                nome=opt$nome, prune_thresh=opt$final_prune_thresh)

  # run famr and save results
  n_expo = ncol(dat$sumstats$betas) - factor_ss$n_factors
  famr_res = run_famr_rss_all(dat, opt$susieL, opt$n_iter, n_expo, 
                              opt$num_samples, annih=opt$annihilate_factors)
  saveRDS(famr_res, opt$out)
}
