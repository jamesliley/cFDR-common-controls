cFDR-common-controls
====================

This repository contains R functions to compute the 'Conditional False Discovery Rate', a method to analyse genome-wide association studies (GWAS) for related diseases originally proposed by Andreasson et al (2013), in the case where the two GWAS have shared controls. These functions widen the scope of the existing technique and enable improved power by allowing larger overall control groups.

Functions are also included to compute an upper bound on the overall false discovery rate amongst SNPs for which cFDR < x, which is almost always higher than x. This provides a robust value at which the false discovery rate can be controlled using this technique.
