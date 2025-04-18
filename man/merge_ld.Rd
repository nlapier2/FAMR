% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{merge_ld}
\alias{merge_ld}
\title{Merge LD for a locus into overall LD matrix}
\usage{
merge_ld(all_ld, locus_ld)
}
\arguments{
\item{all_ld}{A numeric matrix or sparse Matrix object containing the linkage
disequilibrium (LD) values between each pair of SNPs from all previous loci.}

\item{locus_ld}{A numeric matrix or sparse Matrix object containing the
LD values between each pair of SNPs in the current locus.}
}
\value{
A numeric matrix or sparse Matrix object containing the SNPs of both
all_ld and locus_ld.
}
\description{
This function is meant to be called from several functions in famr.R that
loop through loci. Given an LD matrix (all_ld) representing all merged loci
so far and an LD matrix for the current locus, merges the two together
into a new LD matrix and returns that. LD between different loci is assumed
to be zero.
}
\keyword{internal}
