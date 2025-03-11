#' Factor-Augmented Mendelian Randomization with SuSiE (FAMR-Susie)
#' 
#' @description
#' This is the main user-facing function that runs all parts of the 
#' FAMR-Susie method, described in the manuscript
#' "Factor-based methods for robust high-dimensional Mendelian randomization".
#' FAMR-Susie is intended for use in "high-dimensional" MR, where there are
#' many genetically-correlated exposures, and simultaneously learns which of 
#' these exposures have a causal effect on the "outcome" trait.
#' FAMR-Susie uses factor analysis to infer and correct for unobserved
#' confounders, and uses sparse polygenic modeling and Bayesian regression to
#' alleviate other challenges such as other violations of MR assumptions and
#' sparsity. For more details on the method, please see the manuscript.
#' 
#' @details
#' The famr_susie method expects summary statistics and 
#' linkage disequilibrium (LD) matrices as input. The user is expected to have
#' split their data into independent loci before running famr_susie.
#' Thus famr_susie expects a directory (in_dir) of summary statistics and 
#' another (ld_dir) for LD matrices, with one file each per locus.
#' The data for the outcome phenotype may be included in the summary statistics
#' files or supplied separately via the y_gwas_file argument.
#' For the formats of these files, see the documentation for gather_sumstats
#' and read_y_gwas, respectively.
#' The remaining arguments are explained below.
#' 
#' 
#' @param in_dir A directory with input rds files containing summary statistics
#' for the exposures, with one file per locus. Each file should be a list object
#' with three matrices with SNPs as rows and exposure traits as columns. These
#' are "betas" (effect size estimates of SNPs on exposures), "stderrs" 
#' (standard errors of the betas), and "Z" (the z-scores, i.e. betas/stderrs).
#' Optionally can also include corresponding vectors betas_y, stderrs_y, and
#' zscores_y for the outcome trait ("y"); alternatively, those can be read in
#' separately via the read_y_gwas function (see below).
#' Finally, it has a vector 'pos' with the SNP names.
#' See gather_sumstats for more information.
#' 
#' @param ld_dir A directory with input rds files containing LD matrices
#' corresponding to the SNPs in the summary statistics files, with one file
#' per locus. These should simply be numerical nSNP-by-nSNP matrices.
#' 
#' @param y_gwas_file An optional file path to the file with the outcome
#' summary statistics. Can either be a TSV file with columns for the SNP name,
#' effect size estimates, and standard errors (see idcol, betacol, secol),
#' or a VCF file in the MRC IEU GWAS format, described here:
#' https://github.com/MRCIEU/gwas-vcf-specification
#' See also read_vcf in utils.R.
#' 
#' @param fa_method The factor analysis (FA) method to be used. This defaults to
#' 'gfa' and generally should not be changed, as FAMR-Susie currently only
#' recognizes 'gfa'. If the user wishes to modify the code to allow different FA
#' methods, they can specify those different methods here. Currently, if the
#' user puts in something other than 'gfa', no factors are inferred.
#' This can potentially be useful if the user wants to assess FAMR-Susie both
#' with and without factors.
#' 
#' @param num_samples The number of samples that were used to generate the
#' exposure summary statistics. If this varies between the exposures, it can be
#' something like the average -- it does not have to be exact.
#' 
#' @param n_iter The number of Expectation Maximization iterations to fit the
#' cTWAS-style algorithm to learn prior effect probabilities for SNPs and
#' exposures/factors. Default: 30.
#' 
#' @param susieL The maximum number of effects (credible sets) allowed in the
#' susie_rss runs. Default: 30.
#' 
#' @param prune_thresh The absolute value correlation threshold used in the LD
#' clumping step performed by FAMR-Susie before any further analysis. This is 
#' done to improve computational speed and stabilize results. We generally
#' recommend a value between 0.3 and 0.9; the default is 0.5.
#' 
#' @param fa_prune_thresh A stricter clumping threshold to generate a set of 
#' SNPs for input to factor analysis with GFA. Default: 0.1. It is not generally
#' recommended to change this value much.
#' 
#' @param final_prune_thresh A stricted clumping threshold to generate a set of
#' SNPs to model in the cTWAS-style framework / susie_rss runs. By default it is
#' 0.1, the same as fa_prune_thresh.
#' 
#' @param annihilate Boolean; whether to regress factor summary statistics out
#' of exposure and outcome summary statistics after the FA stage and before
#' polygenic modeling and Bayesian regression. Setting this to TRUE will likely
#' yield a conservative result. Default: FALSE.
#' 
#' @param idcol Column of y_gwas_file with the SNP IDs (names). Default is 1,
#' i.e. the first column. Not used if a VCF is provided.
#' 
#' @param betacol Column of y_gwas_file with the SNP effect sizes (betas). 
#' Default is 2, i.e. the second column. Not used if a VCF is provided.
#' 
#' @param secol Column of y_gwas_file with the SNP effect standard errors. 
#' Default is 3, i.e. the third column. Not used if a VCF is provided.
#' 
#' @param header Boolean indicating whether the y_gwas_file has a header line,
#' e.g. with column names. Default=TRUE. Not used if a VCF is provided.
#' 
#' 
#' 
#' @returns A list with some or all of the following elements:
#'
#' \describe{
#' \item{pips}{A list of the learned PIPs for each exposure, factor, and SNP.}
#' 
#' \item{zscores}{A list of the marginal association z-scores learned for each
#' exposure, factor, and SNP with the outcome.}
#' 
#' \item{R}{A matrix of the estimated correlations between each pair of 
#' variables (exposures, factors, and SNPs).}
#' 
#' \item{susieres}{The full susie_rss results object. See the susieR 
#' documentation for more details.}
#' 
#' \item{priors}{A list containing the learned prior effect probabilities and
#' prior effect variances learned for each class of variables (i.e. SNPs and
#' exposures/factors).}
#' 
#' \item{posterior_mean}{A list containing the posterior means learned for 
#' each variable by susie_rss.}
#' 
#' \item{posterior_var}{A list containing the posterior variances learned for 
#' each variable by susie_rss.}
#' 
#' \item{factor_corrs}{A matrix of correlations between the learned factor
#' summary statistics and those of the exposures. Deprecated in favor of the
#' R matrix.}
#' }
#' 
#' @export
famr_susie = function(in_dir, ld_dir, y_gwas_file='NONE', fa_method='gfa', 
                      num_samples=10000, n_iter=30, susieL=30, 
                      prune_thresh=0.5, fa_prune_thresh=0.1, 
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


#' Read in outcome trait summary statistics
#' 
#' @description
#' Optional helper function to read in summary statistics for the outcome trait,
#' used if they are not already provided in the summary statistics files read
#' by in_dir (see famr_susie function documentation).
#' 
#' @param y_gwas_file File path to the file with the outcome
#' summary statistics. Can either be a TSV file with columns for the SNP name,
#' effect size estimates, and standard errors (see idcol, betacol, secol),
#' or a VCF file in the MRC IEU GWAS format, described here:
#' https://github.com/MRCIEU/gwas-vcf-specification
#' See also read_vcf in utils.R.
#' 
#' @param idcol Column of y_gwas_file with the SNP IDs (names). Not used if a 
#' VCF is provided.
#' 
#' @param betacol Column of y_gwas_file with the SNP effect sizes (betas). 
#' Not used if a VCF is provided.
#' 
#' @param secol Column of y_gwas_file with the SNP effect standard errors. 
#' Not used if a VCF is provided.
#' 
#' @param header Boolean indicating whether the y_gwas_file has a header line,
#' e.g. with column names. Default=TRUE. Not used if a VCF is provided.
#' 
#' @returns Returns a list with four items: names_y, betas_y, stderrs_y, and
#' zscores_y, giving the IDs, estimated effect sizes, standard errors, and 
#' z-scores of the SNPs on the outcome trait, respectively.
#' 
#' @importFrom data.table fread
#' @importFrom dplyr %>% select
#' @importFrom rlang .data
#' 
#' @export
read_y_gwas = function(y_gwas_file, idcol, betacol, secol, header=T) {
  message('Reading outcome GWAS...')
  if(y_gwas_file == 'NONE')  return(c())
  if(endsWith(y_gwas_file, '.vcf.gz') || endsWith(y_gwas_file, '.vcf')) {
    gwas = read_vcf(y_gwas_file)
    gwas = gwas %>% dplyr::select(.data$ID, .data$BETA, .data$SE)
  } else {
    gwas = fread(y_gwas_file, sep='\t', header=header,
                 select=c(idcol, betacol, secol))
  }
  colnames(gwas)  = c('names_y', 'betas_y', 'stderrs_y')
  gwas$zscores_y = gwas$betas_y / gwas$stderrs_y
  return(gwas)
}


#' Regress factors out of exposures and outcome
#' 
#' @description
#' This is an optional method that regresses factor summary statistics out
#' of exposure and outcome summary statistics after the factor analysis stage 
#' and before polygenic modeling and Bayesian regression. 
#' This will likely yield a conservative result. By default, famr_susie does
#' not run this. See gather_sumstats for more details.
#' 
#' @param sumstats A list object, usually read by gather_sumstats. The list has
#' three matrices with SNPs as rows and exposure traits as columns. These
#' are "betas" (effect size estimates of SNPs on exposures), "stderrs" 
#' (standard errors of the betas), and "Z" (the z-scores, i.e. betas/stderrs).
#' It also has three corresponding vectors betas_y, stderrs_y, and
#' zscores_y for the outcome trait ("y"). It also has a vector 'pos' with
#' the SNP names.
#' 
#' @param factor_ss Like sumstats, this is a list object with betas, stderrs, 
#' and Z matrices, except for the factors instead of the exposures. Usually
#' generated by generate_factors. It must also have an integer n_factors
#' giving the number of factors. It must have the same set of SNPs as sumstats.
#' 
#' @returns A modified version of the sumstats object. The summary statistics for
#' the exposures and outcome will have had the factors projected (regressed) out
#' of them. It adds the fields orig_betas_y and orig_zscores_y, which are the
#' values those fields had before annihilation, and annih_y, which is the
#' portion of the outcome that was projected out.
#' Thus, orig_betas_y = betas_y + annih_y in the returned object.
#' 
#' @export
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


#' Learn regularization parameters for Susie-RSS weight learning
#' 
#' @description
#' This function learns a "lambda" parameter to regularize the LD matrix of SNPs
#' in case of a mismatch in populations between an LD reference panel and the
#' population the summary statistics were gathered in. This improves the 
#' stability of weights learned by Susie-RSS in learn_weights. The lambdas
#' will be between 0 and 1 and there is one learned for each trait in each 
#' independent locus. lambda=0 is no regularization, lambda=1 sets R to the 
#' identity matrix. 
#' 
#' @param sumstats A list object, usually read by gather_sumstats. The list has
#' three matrices with SNPs as rows and exposure traits as columns. These
#' are "betas" (effect size estimates of SNPs on exposures), "stderrs" 
#' (standard errors of the betas), and "Z" (the z-scores, i.e. betas/stderrs).
#' See gather_sumstats for more details.
#' 
#' @param ld numeric LD matrix corresponding to the SNPs in sumstats
#' 
#' @param oracle_mode DEPRECATED. used only for internal testing. in this case,
#' it automatically returns 0 for the lambda parameters, i.e. no regularization.
#' 
#' @returns a numerical vector of learned lambda parameters for each trait for 
#' this locus. each is between 0 and 1.
#' 
#' @export
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


#' Learn weights of SNP effects on exposures and factors
#' 
#' @description
#' Learns the "weights" (joint effect sizes) of SNPs on each exposure and
#' factor, using Susie-RSS. Uses the lambda parameters learned by learn_lambda
#' for regularization. Adjusts other Susie parameters adaptively according to
#' lambda, i.e. the "L" (number of effects) parameter and clumping threshold.
#' 
#' @param sumstats A list object, usually read by gather_sumstats. The list has
#' three matrices with SNPs as rows and exposure traits as columns. These
#' are "betas" (effect size estimates of SNPs on exposures), "stderrs" 
#' (standard errors of the betas), and "Z" (the z-scores, i.e. betas/stderrs).
#' See gather_sumstats for more details.
#' 
#' @param ld numeric LD matrix corresponding to the SNPs in sumstats
#' 
#' @param lambda numeric vector of lambda parameters, one per trait, learned by
#' learne_lambda.
#' 
#' @param nome DEPRECATED. used only for internal testing. in this case,
#' it automatically returns the "betas" as the weights, as it assumes a no
#' measurement error (nome) model.
#' 
#' @param N number of samples used to generate the exposure data, passed as an
#' argument to susie-rss.
#' 
#' @param L maximum number of effects (per-trait) susie-rss can model
#' 
#' @param prune_min this sets the minimum LD clumping threshold allowed to be
#' applied for any trait, after adjusting for the lambda parameter.
#' 
#' @param prune_max this sets the maximum LD clumping threshold allowed to be
#' applied for any trait, after adjusting for the lambda parameter.
#' 
#' @returns a numeric matrix of weights learned for each SNP on each exposure
#' and factor, with one row per SNP and one column per factor.
#' 
#' @export
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


#' Precompute values for Z-scores and R for this locus
#' 
#' @description
#' FAMR-Susie builds a Z-score for each exposure and factor with the outcome, as
#' well as a matrix of their correlations with one another (R). In general, 
#' these would be constructed using all SNPs across the genome with an effect on
#' any exposure or factor. However, that can yield a very large R matrix. 
#' Instead, we can "precompute" some values needed to compute Z and R separately
#' for each independent locus, which is what this method does, and then build
#' Z and R later using the precomputed values from each locus.
#' 
#' @param weights a numeric matrix of "weights" of each SNP on each exposure and
#' factor in this locus, usually learned by learn_weights, with one row per SNP
#' and one column per exposure/factor.
#' 
#' @param z_gy a vector of SNP Z-scores on the outcome trait
#' 
#' @param ld the LD matrix corresponding to the SNPs in weights
#' 
#' @param orig_z_gy SNP z-scores on the outcome trait prior to annihilation, if
#' applicable. see annihilate_factors.
#' 
#' @param n_factors number of inferred factors
#' 
#' @returns a list with three values, "numerator", "wRw", and "Rgx". numerator
#' is the numerator of the z-score for this locus. wRw is a normalizing matrix
#' of associations between traits, used for computing both Z and R. "Rgx" is
#' a matrix of SNP-trait correlations, with one row per SNP and one column per
#' exposure/factor.
#' 
#' @export
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


#' Gather summary statistics for exposures
#' 
#' @description
#' Given directories containing summary statistics and LD matrices for one or
#' more independent loci, reads in this data and merges it with data for the
#' outcome trait, if applicable. Also performs LD clumping to reduce 
#' computational cost.
#' 
#' @param in_dir A directory with input rds files containing summary statistics
#' for the exposures, with one file per locus. Each file should be a list object
#' with three matrices with SNPs as rows and exposure traits as columns. These
#' are "betas" (effect size estimates of SNPs on exposures), "stderrs" 
#' (standard errors of the betas), and "Z" (the z-scores, i.e. betas/stderrs).
#' Optionally can also include corresponding vectors betas_y, stderrs_y, and
#' zscores_y for the outcome trait ("y"); alternatively, those can be read in
#' separately via the read_y_gwas function.
#' Finally, it has a vector 'pos' with the SNP names.
#' 
#' @param ld_dir A directory with input rds files containing LD matrices
#' corresponding to the SNPs in the summary statistics files, with one file
#' per locus. These should simply be numerical nSNP-by-nSNP matrices.
#' 
#' @param y_gwas A list with four items: names_y, betas_y, stderrs_y, and
#' zscores_y, giving the IDs, estimated effect sizes, standard errors, and 
#' z-scores of the SNPs on the outcome trait, respectively. Usually produced
#' by the read_y_gwas function. If data for the outcome is included in the files
#' in in_dir, leave this blank (default).
#' 
#' @param oracle_mode DEPRECATED. Used for internal testing purposes.
#' 
#' @param prune_thresh LD clumping threshold; a numeric value between 0 and 1.
#' Default is 0.5.
#' 
#' @param fa_prune_thresh Stricter LD clumping threshold used to generate a set
#' of SNPs to be input to factor analysis. Default: 0.1.
#' 
#' @returns A list object with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{weights}{Currently this is the same as betas; later used to store
#' learned weights of SNPs on exposures and factors.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{betas_y}{A vector of estimated effect sizes of SNPs on the outcome.}
#' 
#' \item{stderrs_y}{A vector of standard errors for betas_y.}
#' 
#' \item{zscores_y}{A vector of z-scores of SNPs on the outcome, equal to
#' betas_y divided by stderrs_y.}
#' 
#' \item{pos}{A vector of SNP names.}
#' 
#' \item{locus_idx}{A vector matching SNPs to which locus they come from, 
#' indexed in the order the loci are read from in_dir.}
#' 
#' \item{lambda}{A matrix of lambda values learned by learn_lambda for each 
#' exposure in each locus, where each row of the matrix corresponds to a locus
#' and the columns are the exposures.}
#' 
#' \item{ss_for_fa}{A list object with the "betas", "stderrs", "Z", "weights", 
#' and "pos" fields, but only for the SNPs to be input to factor analysis.}
#' }
#' 
#' @importFrom Matrix Matrix
#' 
#' @export
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


#' Gather summary statistics for exposures from R objects
#' 
#' @description
#' Similar to gather_sumstats, but takes R objects as input instead of 
#' directories containing files. In other words, it assumes the data has
#' already been read in, so we just need to merge in the outcome data, 
#' perform LD clumping, and learn lambda regularization parameters.
#' 
#' @param all_ss A list object with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{weights}{Currently this is the same as betas; later used to store
#' learned weights of SNPs on exposures and factors.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{pos}{A vector of SNP names.}
#' 
#' \item{locus_idx}{A vector matching SNPs to which locus they come from, 
#' indexed in the order the loci are read from in_dir.}
#' }
#' 
#' @param all_ld A list mapping locus indices (see above) to the corresponding
#' LD matrices for that locus. These should be numerical nSNP-by-nSNP matrices.
#' 
#' @param y_gwas A list with four items: names_y, betas_y, stderrs_y, and
#' zscores_y, giving the IDs, estimated effect sizes, standard errors, and 
#' z-scores of the SNPs on the outcome trait, respectively. Usually produced
#' by the read_y_gwas function. If data for the outcome is included in the files
#' in in_dir, leave this blank (default).
#' 
#' @param oracle_mode DEPRECATED. Used for internal testing purposes.
#' 
#' @param prune_thresh LD clumping threshold; a numeric value between 0 and 1.
#' Default is 0.5.
#' 
#' @param fa_prune_thresh Stricter LD clumping threshold used to generate a set
#' of SNPs to be input to factor analysis. Default: 0.1.
#' 
#' @returns A list object with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{betas_y}{A vector of estimated effect sizes of SNPs on the outcome.}
#' 
#' \item{stderrs_y}{A vector of standard errors for betas_y.}
#' 
#' \item{zscores_y}{A vector of z-scores of SNPs on the outcome, equal to
#' betas_y divided by stderrs_y.}
#' 
#' \item{pos}{A vector of SNP names.}
#' 
#' \item{locus_idx}{A vector matching SNPs to which locus they come from, 
#' indexed in the order the loci are read from in_dir.}
#' 
#' \item{lambda}{A matrix of lambda values learned by learn_lambda for each 
#' exposure in each locus, where each row of the matrix corresponds to a locus
#' and the columns are the exposures.}
#' 
#' \item{ss_for_fa}{A list object with the "betas", "stderrs", "Z", "weights", 
#' and "pos" fields, but only for the SNPs to be input to factor analysis.}
#' }
#' 
#' @importFrom Matrix Matrix
#' 
#' @export
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



#' Generate factors for FAMR-Susie
#' 
#' @description
#' Generate factors given the summary statistics read in by gather_sumstats.
#' Since the factors are usually learned with a more strictly LD-clumped set
#' of summary statistics, this function also imputes the SNP-factor effect size
#' estimates for the rest of the SNPs.
#' 
#' @param fa_method A string providing the factor analysis method to use.
#' Currently this only recognizes 'gfa', but this could be extended to recognize
#' other options in the future. If anything other than 'gfa' is passed in, the
#' method assumes no factors should be generated.
#' 
#' @param sumstats Input summary statistics to use for factor analysis, usually
#' from gather_sumstats. A list object with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{pos}{A vector of SNP names.}
#' }
#' 
#' @param N Sample size of the studies used to generate GWAS summary statistics.
#' 
#' @param given_factors DEPRECATED. A file path to factor summary statistics 
#' that have been pre-generated by the user. Used for internal testing purposes.
#' 
#' @param full_ss Similar to sumstats, but contains the full set of SNPs to be
#' used in further analysis, i.e. if prune_thresh is greater than 
#' fa_prune_thresh in gather_sumstats.
#' 
#' @returns A list object containing inferred summary statistics for factors,
#' with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{pos}{A vector of SNP names.}
#' 
#' \item{n_factors}{Number of inferred factors.}
#' }
#' 
#' @export
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


#' Learn weights for SNPs, compute values for Z and R, merge data across loci
#' 
#' @description
#' Essentially this is a large wrapper function that takes the summary 
#' statistics for exposures and factors and prepares the data necessary to
#' impute Z-scores and the R matrix used in the cTWAS/Susie regression step.
#' Merges sumstats for exposures and factors, learns weights on both exposures 
#' and factors, precomputes values needed to impute entries of the Z-scores 
#' vector and R matrix corresponding to exposures and factors, and LD clumps 
#' and merges all of this data into a single sumstats object and LD object.
#' 
#' @param all_sumstats A list object, usually generated by gather_sumstats, 
#' with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{betas_y}{A vector of estimated effect sizes of SNPs on the outcome.}
#' 
#' \item{stderrs_y}{A vector of standard errors for betas_y.}
#' 
#' \item{zscores_y}{A vector of z-scores of SNPs on the outcome, equal to
#' betas_y divided by stderrs_y.}
#' 
#' \item{pos}{A vector of SNP names.}
#' 
#' \item{locus_idx}{A vector matching SNPs to which locus they come from, 
#' indexed in the order the loci are read from in_dir.}
#' 
#' \item{lambda}{A matrix of lambda values learned by learn_lambda for each 
#' exposure in each locus, where each row of the matrix corresponds to a locus
#' and the columns are the exposures.}
#' }
#' 
#' @param factor_ss A list object containing inferred summary statistics for 
#' factors, usually from generate_factors, with the following elements:
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' 
#' \item{pos}{A vector of SNP names.}
#' 
#' \item{n_factors}{Number of inferred factors.}
#' }
#' 
#' @param N Number of samples that summary statistics were generated using.
#' 
#' @param nome DEPRECATED. Boolean value indicating whether to assume a "no
#' measurement error" model. Default: FALSE. Used for internal testing purposes.
#' 
#' @param prune_thresh Numeric value between 0 and 1 indicating a final LD
#' clumping threshold to be applied to the data before merging the loci.
#' SNPs removed in this step will still be used to learn weights and compute
#' Z and R, but will not be included in the final cTWAS/Susie regression.
#' 
#' @returns A list object with the following elements:
#' 
#' \describe{
#' \item{sumstats}{A list object containing summary statistics for all exposures
#' and factors; see the all_sumstats parameter for format.}
#' 
#' \item{ld}{A sparse Matrix object with the LD for the SNPs in the sumstats 
#' object.}
#' 
#' \item{numerator}{A numeric vector representing the numerator of the Z-scores
#' for each exposure and factor.}
#' 
#' \item{wRw}{A numeric matrix of associations between traits (exposures and 
#' factors), used for computing both Z and R}
#' 
#' \item{Rgx}{A numeric matrix of inferred correlations between SNPs and 
#' traits (exposures and outcomes), with one row per SNP and one column per 
#' exposure/factor.}
#' }
#' 
#' @importFrom Matrix Matrix
#' 
#' @export
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


#' Run cTWAS-style regression framework
#' 
#' @description
#' Given the output from learn_wts_precomp_merge, computes the final Z-score
#' vector and R matrix, uses a cTWAS-style framework to infer priors for each
#' "class" of variables (SNPs and exposures/factors), and then generates the
#' final PIPs for each variable using Susie-RSS. Also runs a separate
#' factor-only regression to isolate factor effects from exposure effects.
#' 
#' @param dat A list object, usually from learn_wts_precomp_merge, with the
#' following elements:
#' 
#' \describe{
#' \item{sumstats}{A list object containing summary statistics for all exposures
#' and factors; see the all_sumstats parameter of learn_wts_precomp_merge 
#' for the format.}
#' 
#' \item{ld}{A sparse Matrix object with the LD for the SNPs in the sumstats 
#' object.}
#' 
#' \item{numerator}{A numeric vector representing the numerator of the Z-scores
#' for each exposure and factor.}
#' 
#' \item{wRw}{A numeric matrix of associations between traits (exposures and 
#' factors), used for computing both Z and R}
#' 
#' \item{Rgx}{A numeric matrix of inferred correlations between SNPs and 
#' traits (exposures and outcomes), with one row per SNP and one column per 
#' exposure/factor.}
#' }
#' 
#' @param L a numeric value representing the maximum number of effects (more
#' precisely, credible sets) Susie-RSS can model on the outcome
#' 
#' @param n_iter number of cTWAS expectation-maximization (EM) iterations to 
#' run for fitting prior parameters; usually 10-30 is good.
#' 
#' @param n_expo number of exposures being analyzed
#' 
#' @param n_samp number of samples used to generate exposure summary statistics
#' 
#' @param annih boolean indicating whether the annihilate flag was used
#' (whether factors were regressed out of exposures and outcome). See
#' annihilate_factors for more details.
#' 
#' @returns A list with some or all of the following elements:
#'
#' \describe{
#' \item{pips}{A list of the learned PIPs for each exposure, factor, and SNP.}
#' 
#' \item{zscores}{A list of the marginal association z-scores learned for each
#' exposure, factor, and SNP with the outcome.}
#' 
#' \item{R}{A matrix of the estimated correlations between each pair of 
#' variables (exposures, factors, and SNPs).}
#' 
#' \item{susieres}{The full susie_rss results object. See the susieR 
#' documentation for more details.}
#' 
#' \item{priors}{A list containing the learned prior effect probabilities and
#' prior effect variances learned for each class of variables (i.e. SNPs and
#' exposures/factors).}
#' 
#' \item{posterior_mean}{A list containing the posterior means learned for 
#' each variable by susie_rss.}
#' 
#' \item{posterior_var}{A list containing the posterior variances learned for 
#' each variable by susie_rss.}
#' 
#' \item{factor_corrs}{A matrix of correlations between the learned factor
#' summary statistics and those of the exposures. Deprecated in favor of the
#' R matrix.}
#' }
#' 
#' @importFrom Matrix Matrix
#' 
#' @export
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
