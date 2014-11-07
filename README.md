cFDR-common-controls
====================

This repository contains R functions to compute the 'Conditional False Discovery Rate', a method to analyse genome-wide association studies (GWAS) for related diseases originally proposed by Andreasson et al (2013), in the case where the two GWAS have shared controls. These functions widen the scope of the existing technique and enable improved power by allowing larger overall control groups.

Functions are also included to compute an upper bound on the overall false discovery rate amongst SNPs for which cFDR < x; which, counterintuitively, is not x. This provides a robust value at which the false discovery rate can be controlled using this technique.

Each file contains R code for a specific purpose:

- cfdr Compute cFDR values from a set of pairs of p values
- exp_quantile Evaluate the conditional expected quantile of a p value in the shared control design.
- exp_quantile_point Evaluate the conditional expected quantile of a p value in the shared control design conditional on a p value bound for a second phenotype
- fdr_bound Computation of overall false discovery rate given a cFDR cutoff
- fit.em  fit a mixture of two normals with means 0 and variances 1 and 1+sigma^2 to a set of Z scores
- simulation.r Simulation to demonstrate adjustment of p values
