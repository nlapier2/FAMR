% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/famr.R
\name{precompute_zxy}
\alias{precompute_zxy}
\title{Precompute values for Z-scores and R for this locus}
\usage{
precompute_zxy(weights, z_gy, ld, orig_z_gy = c(), n_factors = 0)
}
\arguments{
\item{weights}{a numeric matrix of "weights" of each SNP on each exposure and
factor in this locus, usually learned by learn_weights, with one row per SNP
and one column per exposure/factor.}

\item{z_gy}{a vector of SNP Z-scores on the outcome trait}

\item{ld}{the LD matrix corresponding to the SNPs in weights}

\item{orig_z_gy}{SNP z-scores on the outcome trait prior to annihilation, if
applicable. see annihilate_factors.}

\item{n_factors}{number of inferred factors}
}
\value{
a list with three values, "numerator", "wRw", and "Rgx". numerator
is the numerator of the z-score for this locus. wRw is a normalizing matrix
of associations between traits, used for computing both Z and R. "Rgx" is
a matrix of SNP-trait correlations, with one row per SNP and one column per
exposure/factor.
}
\description{
FAMR-Susie builds a Z-score for each exposure and factor with the outcome, as
well as a matrix of their correlations with one another (R). In general,
these would be constructed using all SNPs across the genome with an effect on
any exposure or factor. However, that can yield a very large R matrix.
Instead, we can "precompute" some values needed to compute Z and R separately
for each independent locus, which is what this method does, and then build
Z and R later using the precomputed values from each locus.
}
