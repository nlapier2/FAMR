# Provides various utility functions used by multiple scripts

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
    print(paste0('Error in GFA (may simply be no factors identified): ', e))
    return(c())
  })
  return(gfa_factors)
}
