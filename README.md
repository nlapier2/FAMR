# Factor-Augmented Mendelian Randomization with Susie (FAMR-Susie)

The `FAMR` package implements the FAMR-Susie method [^1] 
described in the manuscript
"Factor-based methods for robust high-dimensional Mendelian randomization" [^2].
FAMR-Susie is intended for use in "high-dimensional" MR, where there are
many genetically-correlated exposures, and simultaneously learns which of 
these exposures have a causal effect on the "outcome" trait.
For example, one application would be to simultaneously infer which of many
metabolites measured in a biobank dataset have a causal effect on heart disease.
FAMR-Susie uses factor analysis to infer and correct for unobserved
confounders, builds sparse polygenic models of exposures and inferred factors,
and then uses sparse Bayesian regression to infer which of the exposures and 
factors are causal.
The input is GWAS summary statistics for the exposures and outcome trait of 
interest, and the primary output is posterior inclusion probabilities (PIPs) 
for the exposures and inferred factors, indicating the estimated probability 
that each of those traits has a causal effect on the outcome.

We recommend interested users read our manuscript [^2] for more details on
our method.
Additionally, the sparse polygenic modeling and Bayesian regression stages
are based around modified versions of the
[cTWAS](https://xinhe-lab.github.io/ctwas/)
and [susieR](https://stephenslab.github.io/susieR/index.html)
packages, so users may find it helpful to look at those for more context.

Please see below for installation instructions and a link to a simple vignette.
We recommend running that vignette to get a feel for using the package and
the input and output format.


### Installation

The simplest way to install this package is to open an R shell and run

```
library(devtools)
devtools::install_github('nlapier2/FAMR')
```

Then you can run

```
library(FAMR)
```

to load the package.


### Quick Start

We recommend the following vignette as a quick introduction to the package
and its input and output formats:

(FORTHCOMING)


### Citation

If you use the package, please cite our manuscript:

(FORTHCOMING)


### Footnotes

[^1]: Factor-Augmented Mendelian Randomization (FAMR) is a general framework
described in the manuscript, and FAMR-Susie is a specific method, but the 
package is simply called "FAMR" for simplicity.

[^2]: Citation (FORTHCOMING)