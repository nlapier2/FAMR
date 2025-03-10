# a simplified re-implementation of cTWAS which does not require TWAS-style
#   prediction models as input and places fewer restrictions on input data.
# designed to be used for mendelian randomization.

if (sys.nframe() == 0) {
  # only runs if script is run directly, as opposed to being sourced
  option_list = list(
    make_option(c("--pgen_dir"), type="character", metavar="character",
                help="File with plink pgen files for exposures"),
    make_option(c("--outcome"), type="character", metavar="character",
                help="File with outcome trait values, one per line"),
    make_option(c("--sumstat_dir"), type="character", metavar="character",
                help="File with summary statistics for exposures"),
    make_option(c("--ld_dir"), type="character",  metavar="character",
                help="File with LD matrix for SNPs"),
    make_option(c("--gwas_fname"), type="character",  metavar="character",
                help="File with summary statistics for outcome"),
    make_option(c("--out"), type="character", default='ctwas_results.rds',
                help="Output RDS file name containing results", metavar="character"),
    make_option(c("--idcol"), type="numeric", default=3,
                help="Column with SNP IDs", metavar="numeric"),
    make_option(c("--a1col"), type="numeric", default=4,
                help="Column with SNP minor allele", metavar="numeric"),
    make_option(c("--a2col"), type="numeric", default=5,
                help="Column with SNP major allele", metavar="numeric"),
    make_option(c("--zcol"), type="numeric", default=11,
                help="Column with Z scores", metavar="numeric"),
    make_option(c("--L"), type="numeric", default=5,
                help="Susie L setting", metavar="numeric"),
    make_option(c("--threshold"), type="numeric", default=5.2,
                help="Filter Z-scores below this threshold", metavar="numeric"),
    make_option(c("--thin"), type="numeric", default=0.1, metavar="numeric",
                help="Fraction of SNPs to keep during thinning process"),
    make_option(c("--n_iter"), type="numeric", default=30, metavar="numeric",
                help="Num of iterations to run EM algorithm to fit params"),
    make_option(c("--no_header"), type="logical", action="store_true", default=FALSE,
                help="No header in sumstat/gwas files", metavar="logical"),
    make_option(c("--use_varbvs"), type="logical", action="store_true", default=FALSE,
                help="Use varbvs instead of susie (only individual-level)", metavar="logical")
  );

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}


# filter and thin SNPs for a locus
thresh_and_thin = function(locus_sumstats, thresh, thin) {
  best = apply(abs(locus_sumstats$Z), 1, max)
  pass_thresh = which(best > thresh)
  n_pass = length(pass_thresh)
  pass_thin = sort(sample(1:n_pass, ceiling(n_pass*thin), replace=F))
  indices = pass_thresh[pass_thin]
  return(indices)
}


# merge sumstats for a locus onto overall sumstats
simple_merge_sumstats = function(all_sumstats, locus_sumstats) {
  if(length(all_sumstats) == 0) {
    return(locus_sumstats)
  }
  all_sumstats$Z = rbind(all_sumstats$Z, locus_sumstats$Z)
  all_sumstats$betas = rbind(all_sumstats$betas, locus_sumstats$betas)
  all_sumstats$stderrs = rbind(all_sumstats$stderrs, locus_sumstats$stderrs)
  all_sumstats$pos = rbind(all_sumstats$pos, locus_sumstats$pos)
  return(all_sumstats)
}


# merge LD for a locus onto overall LD matrix
merge_ld = function(all_ld, locus_ld) {
  if(length(all_ld) == 0) {
    return(locus_ld)
  }
  if(length(all_ld) == 1) {
    alen = 1
  } else {
    alen = dim(all_ld)[1]
  }
  if(length(locus_ld) == 1) {
    llen = 1
  } else {
    llen = dim(locus_ld)[1]
  }
  newdim = alen + llen
  all_ld = rbind(all_ld, matrix(0, llen, alen))
  all_ld = cbind(all_ld, matrix(0, newdim, llen))
  start = newdim - llen + 1
  all_ld[start:newdim, start:newdim] = locus_ld
  return(all_ld)
}


# read genotypes in as well as the SNP names
read_genos = function(pgen_fname) {
  pvar_fname = paste0(strsplit(pgen_fname, 'pgen')[[1]], 'pvar')
  pvar = NewPvar(pvar_fname)
  pgen = NewPgen(pgen_fname, pvar=pvar)
  count = GetVariantCt(pvar)
  genos = ReadList(pgen, 1:count)
  var_ids = sapply(1:count, function(x) GetVariantId(pvar, x))
  if(is.null(nrow(genos))) {  # if one-dimensional
    genos = as.matrix(genos)  # change vector to matrix
  }
  colnames(genos) = var_ids
  return(genos)
}


# read in genotypes and SNP weights and filter
read_and_filter_genos = function(pgen_dir, sumstat_dir, thresh, thin) {
  all_genos = c()
  all_sumstats = list()
  pgen_fnames = list.files(path=pgen_dir, pattern='*.pgen', full.names=T)
  sumstat_fnames = list.files(path=sumstat_dir, pattern='*.rds', full.names=T)

  for(i in 1:length(pgen_fnames)) {
    genos = read_genos(pgen_fnames[i])
    sumstats = readRDS(sumstat_fnames[i])

    # apply threshold and thinning procedures, filter SNPs accordingly
    indices = thresh_and_thin(sumstats, thresh, thin)
    print(paste0(length(indices), ' SNPs passed filters in ', sumstat_fnames[i]))
    if(length(indices) == 0) {
      next
    }
    sumstats$betas = sumstats$betas[indices, ]
    sumstats$pos = sumstats$pos[indices, ]
    genos = genos[, indices]

    # now merge the SNPs retained from this locus onto overall genos & sumstats
    all_genos = cbind(all_genos, genos)
    all_sumstats = simple_merge_sumstats(all_sumstats, sumstats)
  }
  return(list('genos' = as.matrix(all_genos), 'sumstats' = all_sumstats))
}


# read in summary statistics and LD matrices and filter
read_sumstats_ld = function(sumstat_dir, ld_dir, thresh, thin) {
  all_sumstats = list()
  all_ld = list()
  sumstat_fnames = list.files(path=sumstat_dir, pattern='*.rds', full.names=T)
  ld_fnames = list.files(path=ld_dir, pattern='*.ld', full.names=T)

  for(i in 1:length(sumstat_fnames)) {
    sumstats = readRDS(sumstat_fnames[i])
    ld = as.matrix(fread(ld_fnames[i]))

    # apply threshold and thinning procedures, filter SNPs accordingly
    indices = thresh_and_thin(sumstats, thresh, thin)
    print(paste0(length(indices), ' SNPs passed filters in ', sumstat_fnames[i]))
    if(length(indices) == 0) {
      next
    }
    sumstats$betas = sumstats$betas[indices, ]
    sumstats$pos = sumstats$pos[indices, ]
    ld = ld[indices, indices]

    # now merge the SNPs retained from this locus onto overall sumstats & LD
    all_sumstats = simple_merge_sumstats(all_sumstats, sumstats)
    all_ld = merge_ld(all_ld, ld)
  }
  return(list('sumstats' = all_sumstats, 'ld' = all_ld))
}


# read in z-scores for SNPs on the outcome trait
read_z_snp = function(opt, sumstats) {
  z_snp = fread(opt$gwas_fname, sep='\t', header=!opt$no_header,
                select=c(opt$idcol, opt$a1col, opt$a2col, opt$zcol))
  names(z_snp) = c('id', 'A1', 'A2', 'z')
  z_snp = z_snp[z_snp$id %in% sumstats$pos$ID, ]
  ids = z_snp$id
  z_snp = as.numeric(z_snp$z)
  names(z_snp) = ids
  return(z_snp)
}


# impute exposure z scores from z_snp
impute_exposure_zscores = function(weights, z_snp, ld) {
  expo_weighted = t(weights %*% z_snp)
  variance = weights %*% ld %*% t(weights)
  z_expo = expo_weighted / sqrt(diag(variance))
  names(z_expo) = colnames(z_expo)
  return(z_expo)
}


# impute exposure betas from b_snp
impute_exposure_betas = function(weights, b_snp, ld) {
  expo_weighted = t(weights %*% b_snp)
  variance = weights %*% ld %*% t(weights)
  b_expo = expo_weighted / diag(variance)
  names(b_expo) = colnames(b_expo)
  return(b_expo)
}


# impute correlations of exposures with SNPs and with each other
impute_exposure_corrs = function(ld, weights, wRw=c(), Rgx=c()) {
  # initialize bottom corner of R as SNP ld matrix
  weights = as.matrix(weights)
  M = sum(dim(weights))  # total number of variables
  M_expo = M - dim(ld)[1]  # number of exposures
  if(typeof(ld) == "S4") {  # sparse matrix
    R = Matrix(0, M, M)
  } else {
    R = matrix(0, M, M)
  }
  R[(M_expo+1):M, (M_expo+1):M] = ld

  # compute exposure-exposure correlations according to cTWAS formula
  if(length(wRw) == 0)  wRw = weights %*% ld %*% t(weights)
  invwts = diag(wRw)
  Rij = wRw / sqrt(invwts %*% t(invwts))
  R[1:M_expo, 1:M_expo] = Rij

  # compute exposure-SNP correlations according to cTWAS formula
  if(length(Rgx) == 0)  Rgx = t(ld %*% t(weights))
  Rgx = apply(Rgx, 2, function(x) x / sqrt(invwts))
  R[1:M_expo, (M_expo+1):M] = Rgx

  # make symmetric and unit diagonal
  diag(R) = 1
  R[lower.tri(R)] = t(R)[lower.tri(R)]
  return(R)
}


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
  if(length(zscores) > 0) {  # sumstats mode
    z_plain = as.numeric(unlist(zscores))
    res = susie_rss(z_plain, R, L=L, # n=n,
                          prior_weights = prior_pi_plain,
                          prior_variance = prior_sigma2_plain,
                          null_weight = nw,
                          estimate_prior_variance=F)
  } else {  # individual level mode
    X = do.call(cbind, pred_vars)  # cbind all matrices in pred_vars
    # scale the variance unless this is the first iteration
    if(max(prior_sigma2_plain) != 0.2 || min(prior_sigma2_plain) != 0.2) {
      prior_sigma2_plain = prior_sigma2_plain / var(y)
    }
    if(use_vb) {
      logodds = log10(prior_pi_plain / (1 - prior_pi_plain))
      logodds_l = matrix(logodds, length(logodds), L)
      sa_l = rep(0.01, L)
      res = varbvs(X, y, Z=NULL, family='gaussian', logodds=logodds_l, sa=sa_l)
    } else {
      res = susie(X, y, L=L, # n=n,
                        prior_weights = prior_pi_plain,
                        scaled_prior_variance = prior_sigma2_plain,
                        null_weight = nw,
                        standardize = T,
                        estimate_prior_variance=F)
    }
  }
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
  } else {  # individual-level mode
    M = dim(pred_vars$exposures)[2] + dim(pred_vars$snps)[2]  # total number of vars
    if('factors' %in% names(pred_vars))  M = M + dim(pred_vars$factors)[2]
    for(k in names(pred_vars)) {
      m_k = ncol(pred_vars[[k]])
      prior_pi[[k]] = rep(1 / M, m_k)
      prior_sigma2[[k]] = rep(0.2, m_k)
      n_vars[[k]] = m_k
    }
  }

  for(iter in 1:n_iter) {
    # iteratively run susie with current priors, update per-variable params
    #   (alpha and tau2), and update group priors (pi and sigma2)
    print(paste0('Running prior estimation iteration ', iter))
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
    print('Current pi:')
    print(new_group_priors$pi)
    print('Current sigma2:')
    print(new_group_priors$sigma2)
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


# wrapper function to run summary statistics mode
run_sumstats_mode = function(opt) {
  # read zscores and ld
  dat = read_sumstats_ld(opt$sumstat_dir, opt$ld_dir, opt$threshold, opt$thin)
  sumstats = dat$sumstats
  ld = dat$ld

  # read/generate z_scores for SNPs and exposures on outcome from GWAS data
  # also harmoize to make sure same SNPs in z_snp, sumstats, & R, and same order
  z_snp = read_z_snp(opt, sumstats)
  idx = which(sumstats$pos$ID %in% names(z_snp))
  sumstats$pos = sumstats$pos[idx, ]
  sumstats$betas = sumstats$betas[idx, ]
  ld = ld[idx, idx]
  print(paste0(length(idx), ' total variants passed filters and matched GWAS file.'))
  weights = t(as.matrix(sumstats$betas))
  z_expo = impute_exposure_zscores(weights, z_snp, ld)
  z_all = list('exposures' = z_expo, 'snps' = z_snp)  # currently two classes

  # compute SNP-exposure correlations (currently only works for two classes)
  R = impute_exposure_corrs(ld, weights)

  # calculate priors w/ EM algorithm, then run cTWAS-rss
  priors = estimate_priors_EM(zscores=z_all, R=R, L=1, n_iter=opt$n_iter)
  print('Finished fitting prior. Computing final results...')
  ctwas_res = get_results(zscores=z_all, R=R, priors=priors, L=opt$L)
  print(ctwas_res$sorted_pips[1:20])
  saveRDS(ctwas_res, opt$out)
}


# wrapper function to run individual-level mode
run_individual_mode = function(opt) {
  # read and filter genotypes and weights; read outcome
  dat = read_and_filter_genos(opt$pgen_dir, opt$sumstat_dir, opt$threshold, opt$thin)
  print(paste0(ncol(dat$genos), ' total variants passed filters.'))
  y = read.table(opt$outcome, header=F)$V1

  # currently impute X with simple OLS weights
  weights = as.matrix(dat$sumstats$betas)
  x_pred = dat$genos %*% weights
  colnames(x_pred) = names(dat$sumstats$betas)
  pred_vars = list('exposures' = x_pred, 'snps' = dat$genos)

  # calculate priors w/ EM algorithm, then run cTWAS
  priors = estimate_priors_EM(pred_vars=pred_vars, y=y, L=1, n_iter=opt$n_iter,
                              use_vb=opt$use_varbvs)
  print('Finished fitting prior. Computing final results...')
  ctwas_res = get_results(pred_vars=pred_vars, y=y, priors=priors, L=opt$L,
                          use_vb=opt$use_varbvs)
  print(ctwas_res$sorted_pips[1:20])
  saveRDS(ctwas_res, opt$out)
}


if (sys.nframe() == 0) {  # if running from shell or Rscript (not sourcing)
  # check whether running individual level mode or summary mode,
  #    then run the appropriate mode
  individual = FALSE
  summary = FALSE
  if(is.null(opt$sumstat_dir)) {
    stop('--sumstat_dir is required.')
  }
  if(!is.null(opt$pgen_dir) && !is.null(opt$outcome)) {
    individual = TRUE
  }
  if(!is.null(opt$ld_dir) && !is.null(opt$gwas_fname)) {
    summary = TRUE
  }
  if(!individual && !summary) {
    stop(paste0('Must provide pgen_dir and outcome for individual level mode,',
                ' or ld_dir and gwas_fname for summary statistics mode.'))
  } else if(individual && summary) {
    stop('Ambiguous input: both individual level and summary statics input provided.')
  } else if(individual) {
    run_individual_mode(opt)
  } else {
    run_sumstats_mode(opt)
  }
}
