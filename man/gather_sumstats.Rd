% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/famr.R
\name{gather_sumstats}
\alias{gather_sumstats}
\title{Gather summary statistics for exposures}
\usage{
gather_sumstats(
  in_dir,
  ld_dir,
  y_gwas = c(),
  oracle_mode = F,
  prune_thresh = 1,
  fa_prune_thresh = 1
)
}
\arguments{
\item{in_dir}{A directory with input rds files containing summary statistics
for the exposures, with one file per locus. Each file should be a list object
with three matrices with SNPs as rows and exposure traits as columns. These
are "betas" (effect size estimates of SNPs on exposures), "stderrs"
(standard errors of the betas), and "Z" (the z-scores, i.e. betas/stderrs).
Optionally can also include corresponding vectors betas_y, stderrs_y, and
zscores_y for the outcome trait ("y"); alternatively, those can be read in
separately via the read_y_gwas function.
Finally, it has a vector 'pos' with the SNP names.}

\item{ld_dir}{A directory with input rds files containing LD matrices
corresponding to the SNPs in the summary statistics files, with one file
per locus. These should simply be numerical nSNP-by-nSNP matrices.}

\item{y_gwas}{A list with four items: names_y, betas_y, stderrs_y, and
zscores_y, giving the IDs, estimated effect sizes, standard errors, and
z-scores of the SNPs on the outcome trait, respectively. Usually produced
by the read_y_gwas function. If data for the outcome is included in the files
in in_dir, leave this blank (default).}

\item{oracle_mode}{DEPRECATED. Used for internal testing purposes.}

\item{prune_thresh}{LD clumping threshold; a numeric value between 0 and 1.
Default is 0.5.}

\item{fa_prune_thresh}{Stricter LD clumping threshold used to generate a set
of SNPs to be input to factor analysis. Default: 0.1.}
}
\value{
A list object with the following elements:

\describe{
\item{betas}{A matrix of estimated effect sizes of SNPs on exposures, with
rows corresponding to SNPs and columns corresponding to exposures.}

\item{weights}{Currently this is the same as betas; later used to store
learned weights of SNPs on exposures and factors.}

\item{stderrs}{A matrix of standard errors for the betas.}

\item{Z}{A matrix of z-scores for the exposures, equal to betas divided by
stderrs.}

\item{betas_y}{A vector of estimated effect sizes of SNPs on the outcome.}

\item{stderrs_y}{A vector of standard errors for betas_y.}

\item{zscores_y}{A vector of z-scores of SNPs on the outcome, equal to
betas_y divided by stderrs_y.}

\item{pos}{A vector of SNP names.}

\item{locus_idx}{A vector matching SNPs to which locus they come from,
indexed in the order the loci are read from in_dir.}

\item{lambda}{A matrix of lambda values learned by learn_lambda for each
exposure in each locus, where each row of the matrix corresponds to a locus
and the columns are the exposures.}

\item{ss_for_fa}{A list object with the "betas", "stderrs", "Z", "weights",
and "pos" fields, but only for the SNPs to be input to factor analysis.}
}
}
\description{
Given directories containing summary statistics and LD matrices for one or
more independent loci, reads in this data and merges it with data for the
outcome trait, if applicable. Also performs LD clumping to reduce
computational cost.
}
