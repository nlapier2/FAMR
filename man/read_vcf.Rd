% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{read_vcf}
\alias{read_vcf}
\title{Read summary statistics from a VCF file}
\usage{
read_vcf(fname)
}
\arguments{
\item{fname}{VCF file name}
}
\value{
A data.frame object with one row per SNP and the following column
fields (omitting those not used by other methods):

\describe{
\item{CHR}{Chromosome the SNP is in.}

\item{ID}{ID (name) of the SNP.}

\item{BETA}{Estimated marginal effect size of the SNP on the trait.}

\item{SE}{Standard errors of the BETAs.}

\item{Z}{Z-scores for the SNPs on the trait.}

\item{logp}{Log p-values for the association between the SNP and trait.}
}
}
\description{
Given a VCF file containing summary statistics for a trait, read in and
preprocess the summary statistics, i.e. performing NA, MAF, QC filtering
and more. The VCF file must be in the MRC IEU GWAS format, described here:
https://github.com/MRCIEU/gwas-vcf-specification
}
