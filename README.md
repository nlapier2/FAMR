# Factor-Augmented Mendelian Randomization with Susie (FAMR-Susie)

This is the repository for the FAMR package, which implements the 
FAMR-Susie method described in the manuscript
"Factor-based methods for robust high-dimensional Mendelian randomization".
Factor-Augmented Mendelian Randomization (FAMR) is a general framework
described in the manuscript, and FAMR-Susie is a specific method, but the 
package is simply called "FAMR" for simplicity.
FAMR-Susie is intended for use in "high-dimensional" MR, where there are
many genetically-correlated exposures, and simultaneously learns which of 
these exposures have a causal effect on the "outcome" trait.
FAMR-Susie uses factor analysis to infer and correct for unobserved
confounders, and uses sparse polygenic modeling and Bayesian regression to
alleviate other challenges such as other violations of MR assumptions and 
sparsity.


This repository contains code modified from the susieR and cTWAS packages,
with permission from the authors.

susieR can be found at https://github.com/stephenslab/susieR

cTWAS can be found at https://github.com/xinhe-lab/ctwas/