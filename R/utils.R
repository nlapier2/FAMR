# Provides various utility functions used by multiple scripts, 
#   and/or intended to be called internally only.

#' Subset summary statistics objects
#' 
#' @description
#' Subset a summary statistics object ("ss") according to the indices provided
#' by idx, i.e. retain only the SNPs at the positions given by idx.
#' 
#' @param ss A summary statistics list object, containing some or all of the
#' following fields (only betas, stderrs, weights, Z required):
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
#' \item{annih_y}{A vector containing the portion of betas_y projected out
#' by factors; see annihilate_factors in famr.R for more details}
#' 
#' \item{orig_betas_y}{A vector containing the original betas_y if the latter
#' was modified by annihilating factors; see annihilate_factors in famr.R 
#' for more details}
#' 
#' \item{orig_zscores_y}{Similar to orig_betas_y but for zscores_y}
#' }
#' 
#' @param idx A vector of SNP indices to retain from ss
#' 
#' @returns A modified version of ss with only the SNPs specified by idx kept.
#' 
subset_sumstats = function(ss, idx) {
  ss$betas = ss$betas[idx, , drop=F]
  ss$stderrs = ss$stderrs[idx, , drop=F]
  ss$weights = ss$weights[idx, , drop=F]
  ss$Z = ss$Z[idx, , drop=F]
  if(is.null(nrow(ss$pos))) {
    ss$pos = ss$pos[idx]
  } else {
    ss$pos = ss$pos[idx, , drop=F]
  }
  if('zscores_y' %in% names(ss)) {
    ss$zscores_y = ss$zscores_y[idx]
    ss$betas_y = ss$betas_y[idx]
    ss$stderrs_y = ss$stderrs_y[idx]
  }
  if('locus_idx' %in% names(ss)) {
    ss$locus_idx = ss$locus_idx[idx]
  }
  if('annih_y' %in% names(ss)) {
    ss$annih_y = ss$annih_y[idx]
    ss$orig_betas_y = ss$orig_betas_y[idx]
    ss$orig_zscores_y = ss$orig_zscores_y[idx]
  }
  return(ss)
}


#' Merge sumstats for a locus into overall sumstats
#' 
#' @description
#' This function is meant to be called from several functions in famr.R that
#' loop through loci. Given an object containing summary statistics merged from
#' all loci so far (all_sumstats) and summary statistics from a new locus,
#' merges these two and returns a new summary statistics object.
#' Also used to merge exposure and factor summary statistics.
#' 
#' @param all_sumstats A summary statistics list object, containing some or all 
#' of the following fields (only betas, stderrs, weights, Z required):
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
#' \item{lambda}{A matrix of lambda values learned by learn_lambda in famr.R for 
#' each exposure in each locus, where each row of the matrix corresponds to a 
#' locus and the columns are the exposures.}
#' }
#' 
#' @param locus_sumstats Like all_sumstats, but for the current locus only.
#' 
#' @param by_row Boolean indicating whether the sumstats should be merged by row
#' or by column. The latter is used when merging exposure and factor summary
#' statistics rather than merging different loci.
#' 
#' @returns A list object in the same format as all_sumstats, with the merged
#' summary statistics from all_sumstats and locus_sumstats.
#' 
merge_sumstats = function(all_sumstats, locus_sumstats, by_row=T) {
  if(length(all_sumstats) == 0 || !('betas' %in% names(all_sumstats))) {
    return(locus_sumstats)
  }
  if(by_row) {
    all_sumstats$Z = rbind(all_sumstats$Z, locus_sumstats$Z)
    all_sumstats$betas = rbind(all_sumstats$betas, locus_sumstats$betas)
    all_sumstats$stderrs = rbind(all_sumstats$stderrs, locus_sumstats$stderrs)
    all_sumstats$weights = rbind(all_sumstats$weights, locus_sumstats$weights)
  } else {
    all_sumstats$Z = cbind(all_sumstats$Z, locus_sumstats$Z)
    all_sumstats$betas = cbind(all_sumstats$betas, locus_sumstats$betas)
    all_sumstats$stderrs = cbind(all_sumstats$stderrs, locus_sumstats$stderrs)
    all_sumstats$weights = cbind(all_sumstats$weights, locus_sumstats$weights)
  }
  if(is.null(nrow(all_sumstats$pos)) && by_row) {
    all_sumstats$pos = c(all_sumstats$pos, locus_sumstats$pos)
  } else if(by_row) {
    all_sumstats$pos = rbind(all_sumstats$pos, locus_sumstats$pos)
  }
  if('zscores_y' %in% names(all_sumstats)) {
    all_sumstats$zscores_y = c(all_sumstats$zscores_y, locus_sumstats$zscores_y)
    all_sumstats$betas_y = c(all_sumstats$betas_y, locus_sumstats$betas_y)
    all_sumstats$stderrs_y = c(all_sumstats$stderrs_y, locus_sumstats$stderrs_y)
  }
  if('locus_idx' %in% names(all_sumstats)) {
    all_sumstats$locus_idx = c(all_sumstats$locus_idx, locus_sumstats$locus_idx)
  }
  if('lambda' %in% names(all_sumstats) && by_row) {
    all_sumstats$lambda = rbind(all_sumstats$lambda, locus_sumstats$lambda)
  }
  return(all_sumstats)
}


#' Merge LD for a locus into overall LD matrix
#' 
#' @description
#' This function is meant to be called from several functions in famr.R that
#' loop through loci. Given an LD matrix (all_ld) representing all merged loci
#' so far and an LD matrix for the current locus, merges the two together
#' into a new LD matrix and returns that. LD between different loci is assumed
#' to be zero.
#' 
#' @param all_ld A numeric matrix or sparse Matrix object containing the linkage
#' disequilibrium (LD) values between each pair of SNPs from all previous loci.
#' 
#' @param locus_ld A numeric matrix or sparse Matrix object containing the 
#' LD values between each pair of SNPs in the current locus.
#' 
#' @returns A numeric matrix or sparse Matrix object containing the SNPs of both
#' all_ld and locus_ld.
#' 
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


#' Read summary statistics from a VCF file
#' 
#' @description
#' Given a VCF file containing summary statistics for a trait, read in and
#' preprocess the summary statistics, i.e. performing NA, MAF, QC filtering
#' and more. The VCF file must be in the MRC IEU GWAS format, described here:
#' https://github.com/MRCIEU/gwas-vcf-specification
#' 
#' @param fname VCF file name
#' 
#' @returns A data.frame object with one row per SNP and the following column
#' fields (omitting those not used by other methods):
#' 
#' \describe{
#' \item{CHR}{Chromosome the SNP is in.}
#' 
#' \item{ID}{ID (name) of the SNP.}
#' 
#' \item{BETA}{Estimated marginal effect size of the SNP on the trait.}
#' 
#' \item{SE}{Standard errors of the BETAs.}
#' 
#' \item{Z}{Z-scores for the SNPs on the trait.}
#' 
#' \item{logp}{Log p-values for the association between the SNP and trait.}
#' }
#' 
#' @importFrom data.table fread setnames
#' @importFrom dplyr %>% filter
#' @importFrom stats na.omit
#' @importFrom rlang .data
#' 
read_vcf = function(fname) {
  # read gwas file, filter for NA, MAF, multi-character allele, duplicates, etc
  gwas = fread(paste0("zgrep -v '^##' ", fname))
  setnames(gwas, '#CHROM', 'CHR')

  # Filter NAs, non-passing QC SNPs, reverse-complements, MHC region, and multi-allelic SNPs
  gwas = gwas %>% na.omit() %>%
    filter(.data$FILTER == 'PASS',
           nchar(.data$REF) == 1, nchar(.data$ALT) == 1,
           !(.data$INFO == "ReverseComplementedAlleles"),
           !(.data$CHR == 6 & .data$POS>=25000000 & .data$POS<=36000000))


  # extract allele frequency (AF) from INFO field, filter at MAF < 0.01
  if(grepl('AF', gwas[['INFO']][1])) {
    AF = as.numeric(sapply(gwas$INFO, function(x) unlist(strsplit(x, '='))[2]))
    gwas = gwas[AF > 0.01 & AF < 0.99, ]  # MAF filter at 0.01
  }

  # gather the betas and stderrs from the last column using the FORMAT column
  #   to tell us which columns tell us those
  format_split = unlist(strsplit(gwas[1,]$FORMAT, ':'))
  beta_col = which(format_split == 'ES')
  se_col = which(format_split == 'SE')
  logp_col = which(format_split == 'LP')
  ss_col = gwas[[names(gwas)[length(gwas)]]]
  ss_split = strsplit(ss_col, ':')
  gwas$BETA = as.numeric(lapply(ss_split, function(x) x[beta_col]))
  gwas$SE = as.numeric(lapply(ss_split, function(x) x[se_col]))
  gwas$logp = as.numeric(lapply(ss_split, function(x) x[logp_col]))
  gwas$Z = gwas$BETA / gwas$SE
  gwas = gwas %>% na.omit()
  return(gwas)
}


#' Run GFA for factor analysis
#' 
#' @description
#' Wrapper function to run GFA to infer factors for FAMR, while catching errors
#' that happen when no factors are inferred and removing factors strongly 
#' associated with only one factor.
#' 
#' @param x_betas Numeric matrix of marginal effect size estimates ("betas") of
#' SNPs on each exposure, with rows corresponding to SNPs and columns 
#' corresponding to exposures.
#' 
#' @param x_stderrs Numeric matrix of standard errors of the x_betas.
#' 
#' @param N DEPRECATED. Sample size used to generate the exposures.
#' 
#' @returns A GFA object; see GFA for more details. For FAMR, the most relevant
#' fields are L_hat and F_hat, which contain the loadings and factors inferred
#' by GFA. L_hat will be of dimension nSNPs-by-nFactors, and t(F_hat) will be of
#' dimension nFactors-by-nExposures.
#' 
run_gfa_full = function(x_betas, x_stderrs, N) {
  gfa_factors = tryCatch({
    gfares = GFA::gfa_fit(B_hat = x_betas, S = x_stderrs)
    return(gfa_factor_prune_full(x_betas, gfares))
  }, error = function(e) {
    message('Error in GFA (may simply be no factors identified): ', e)
    return(c())
  })
  return(gfa_factors)
}


#' Prune GFA factors associated with only one exposure
#' 
#' @description
#' GFA occasionally infers factors that are very strongly associated with only
#' one exposure and not with others, so that they are essentially a proxy for
#' that exposure. This helper function removes those factors heuristically, by
#' removing factors that do not have an R^2 above some value ("thresh") with at
#' least "min_count" exposures.
#' 
#' @param x_betas Numeric matrix of marginal effect size estimates ("betas") of
#' SNPs on each exposure, with rows corresponding to SNPs and columns 
#' corresponding to exposures.
#' 
#' @param gfares GFA results object, usually generated by run_gfa_full.
#' 
#' @param thresh The R^2 threshold with exposures that the factors in gfares
#' must meet to avoid being removed. Default: 0.1.
#' 
#' @param min_count The number of exposures a factor must be sufficiently
#' correlated with to avoid being removed. Default: 2.
#' 
#' @returns A GFA object (see GFA for more details); a modified version of 
#' gfares that contains only the factors that passed the filter.
#' 
gfa_factor_prune_full = function(x_betas, gfares, thresh = 0.1, min_count = 2) {
  r2_matrix = cor(x_betas, gfares$L_hat)^2
  keep_indices = which(apply(r2_matrix, 2, function(x) sum(x > thresh)) >= min_count)
  gfares$L_hat = gfares$L_hat[, keep_indices]
  gfares$F_hat = gfares$F_hat[, keep_indices]
  gfares$gfa_pve$pve = gfares$gfa_pve$pve[, keep_indices]
  return(gfares)
}


#' Form correlation matrix of all variables
#' 
#' @description
#' Generates a correlation matrix between all pairs of variables, including
#' SNPs, exposures, and factors, given the weights learned by learn_weights in f
#' amr.R and the LD between the SNPs. Can optionally take in wRw and Rgx
#' (see below) if those have already been computed.
#' 
#' @param ld A numeric matrix or sparse Matrix object with the linkage
#' disequilibrium (LD) values between each pair of SNPs, for the SNPs in the
#' weights parameter/
#' 
#' @param weights A numeric matrix of SNP weights on the exposures and factors,
#' with rows corresponding to SNPs and columns to expsosures and factors.
#' 
#' @param wRw A numeric matrix of associations between traits (exposures and 
#' factors), used for computing both Z and R. See precompute_zxy in famr.R
#' for more details.
#' 
#' @param Rgx A numeric matrix of inferred correlations between SNPs and 
#' traits (exposures and outcomes), with one row per SNP and one column per 
#' exposure/factor.See precompute_zxy in famr.R for more details.
#' 
#' @returns A numeric matrix or sparse Matrix object containing the inferred 
#' correlations between all SNPs, exposures, and factors.
#' 
#' @importFrom Matrix Matrix
#' 
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


#' Write FAMR-Susie results to a file
#' 
#' @description
#' Writes FAMR-Susie results to a file, and optionally a simple results summary
#' to an additional text file. Mainly intended to be used internally for
#' simulations.
#' 
#' @param out Base name for output files
#' 
#' @param famr_res An FAMR-Susie results object; see the famr_susie method in
#' famr.R for more details.
#' 
#' @param K The number of exposures in the analysis.
#' 
#' @param fa_method Which factor analysis method was used to infer factors
#' for FAMR-Susie.
#' 
#' @param for_sim A boolen indicating whether to generate an additional text
#' file output for use with internal simulations.
#' 
#' @returns nothing
#' 
#' @importFrom utils write.table
#' 
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


#' LD clump SNPs in a locus
#' 
#' @description
#' A simple helper function to LD clump SNPs in a locus. Given a clumping
#' threshold and a summary statistics object and LD matrix, returns a list of
#' SNP indices to keep after clumping. SNPs are sorted in descending order by
#' highest Z-score with any exposure, and SNPs are marked for removal when they
#' have high enough LD with a SNP with a stronger peak Z-score.
#' The kept indices are usually passed to subset_sumstats to actually 
#' remove the SNPs.
#' 
#' @param sumstats A summary statistics list object, containing the
#' following fields (and possibly others):
#' 
#' \describe{
#' \item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
#' rows corresponding to SNPs and columns corresponding to exposures.}
#' 
#' \item{stderrs}{A matrix of standard errors for the betas.}
#' 
#' \item{Z}{A matrix of z-scores for the exposures, equal to betas divided by 
#' stderrs.}
#' }
#' 
#' @param ld A numeric matrix of LD values between SNPs in sumstats
#' 
#' @param prune_thresh A numeric value between 0 and 1 indicating the LD
#' clumping threshold.
#' 
#' @returns A vector of SNP indices marked to be retained after LD clumping.
#' 
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


#' Wrapper function to clump and merge SNPs
#' 
#' @description
#' A simple wrapper function that performs LD clumping for SNPs in a locus,
#' then subsets those SNPs, then merges the SNPs for this locus into an object
#' collecting summary statistics over all loci. Used by gather_sumstats to read
#' individual locus files into a single object in memory.
#' 
#' @param dat A list object formed from merging all previous loci, 
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
#' 
#' \item{ss_for_fa}{A list object with the "betas", "stderrs", "Z", "weights", 
#' and "pos" fields, but only for the SNPs to be input to factor analysis.}
#' }
#' 
#' @param sumstats A list object in the same format as dat, but corresponding
#' to the current locus
#' 
#' @param ld A numertic LD matrix between the SNPs in sumstats (current locus)
#' 
#' @param prune_thresh A numeric value between 0 and 1 indicating the LD
#' clumping threshold to be used by ld_prune_famr
#' 
#' @param fa_prune_thresh A numeric value between 0 and 1 indicating the LD
#' clumping threshold to be used to generate a set of SNPs for input to factor
#' analysis (FA)
#' 
#' @param f_idx The current file index, used to indicate which locus the SNPs in
#' sumstats come from
#' 
#' @param oracle_mode DEPRECATED. Used for internal testing only.
#' 
#' @returns A list object in the same format as dat, containing the SNPs in both
#' dat and those from sumstats that were not removed during LD clumping.
#' 
prune_and_merge = function(dat, sumstats, ld, prune_thresh, fa_prune_thresh, 
                           f_idx, oracle_mode=F) {
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


#' Merge outcome and exposure summary statistics
#' 
#' @description
#' Given exposure summary statistics and outcome summary statistics,
#' intersects the SNPs, then merges them into a single object. Used if
#' read_y_gwas is run, i.e. if outcome summary statistics are not already
#' included in the exposure summary statistics files.
#' 
#' @param sumstats A list object with the following elements, i.e. summary
#' statistics for exposures:
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
#' }
#' 
#' @param ld A numeric LD matrix between the SNPs in sumstats.
#' 
#' @param y_gwas A list object with the following elements, i.e. summary
#' statistics for the outcome:
#' 
#' \describe{
#' \item{names_y}{A vector of SNP IDs (names)}
#' 
#' \item{betas_y}{A vector of estimated effect sizes of SNPs on the outcome.}
#' 
#' \item{stderrs_y}{A vector of standard errors for betas_y.}
#' 
#' \item{zscores_y}{A vector of z-scores of SNPs on the outcome, equal to
#' betas_y divided by stderrs_y.}
#' }
#' 
#' @returns A list object with the following elements:
#' 
#' \describe{
#' \item{sumstats}{A list with the fields of both sumstats and y_gwas (except
#' names_y, which is merged into pos), containing the merged summary statistics}
#' 
#' \item{ld}{A numeric LD matrix for the SNPs in the new sumstats object; 
#' the same as the input LD matrix except without the SNPs that are unmeasured
#' in y_gwas.}
#' }
#' 
merge_outcome = function(sumstats, ld, y_gwas) {
  if(is.null(nrow(sumstats$pos))) {  # if pos is only a vector of IDs
    y_gwas = y_gwas[y_gwas$names_y %in% sumstats$pos, ]
    idx = which(sumstats$pos %in% y_gwas$names_y)
  } else {  # if pos is a list/df with an ID row
    y_gwas = y_gwas[y_gwas$names_y %in% sumstats$pos$ID, ]
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
