# Provides various utility functions used by multiple scripts, or intended to
#    be called internally only.

# helper function to subset summary statistics based on provided indices
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


# merge sumstats for a locus onto overall sumstats
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


# read in VCF file, clean data, and preprocess for use in real data analysis
read_vcf = function(fname) {
  # read gwas file, filter for NA, MAF, multi-character allele, duplicates, etc
  gwas = fread(paste0("zgrep -v '^##' ", fname))
  setnames(gwas, '#CHROM', 'CHR')

  # Filter NAs, non-passing QC SNPs, reverse-complements, MHC region, and multi-allelic SNPs
  gwas = gwas %>% na.omit() %>%
    filter(FILTER == 'PASS',
           nchar(REF) == 1, nchar(ALT) == 1,
           !(INFO == "ReverseComplementedAlleles"),
           !(CHR == 6 & POS>=25000000 & POS<=36000000))


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


# helper function to heuristically prune excess factors for GFA: prune factors
#   that don't explain at least "thresh" variance of at least "min_count" exposures
gfa_factor_prune_full = function(x_betas, gfares, thresh = 0.1, min_count = 2) {
  r2_matrix = cor(x_betas, gfares$L_hat)^2
  keep_indices = which(apply(r2_matrix, 2, function(x) sum(x > thresh)) >= min_count)
  gfares$L_hat = gfares$L_hat[, keep_indices]
  gfares$F_hat = gfares$F_hat[, keep_indices]
  gfares$gfa_pve$pve = gfares$gfa_pve$pve[, keep_indices]
  return(gfares)
}


# run gfa to estimate confounders, and return the full GFA object
run_gfa_full = function(x_betas, x_stderrs, N) {
  gfa_factors = tryCatch({
    gfares = gfa_fit(B_hat = x_betas, S = x_stderrs)
    return(gfa_factor_prune_full(x_betas, gfares))
  }, error = function(e) {
    message('Error in GFA (may simply be no factors identified): ', e)
    return(c())
  })
  return(gfa_factors)
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

