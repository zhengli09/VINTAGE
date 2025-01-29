# VINTAGE: Variant-set test INtegrative TWAS for GEne-based analysis

## An alternative framework for transcriptome-wide association studies to detect and decipher gene-trait associations

![scheme](https://github.com/zhengli09/VINTAGE/blob/master/docs/schematic.png)

<p align="justify">

Integrative analysis of genome-wide association studies (GWAS) with gene
expression mapping studies holds the promise to unravel the molecular
mechanisms underlying disease etiology. Here, we present VINTAGE, an
alternative statistical framework for such integrative analysis to identify
and decipher gene-trait associations. VINTAGE explicitly quantifies and
tests the proportion of genetic effects on a trait potentially mediated
through gene expression using a local genetic correlation test, and
further leverages such information to guide the integration of gene
expression mapping study towards gene association mapping in GWAS
through a genetic variance test. The explicit quantification of local
genetic correlation in VINTAGE allows its gene association test to unify
two seemingly unrelated methods, SKAT and TWAS, into the same analytic
framework and include both as special cases, thus achieving robust
performance across a range of scenarios. The framework of VINTAGE offers
a clearer understanding of TWAS false signals when gene expression is not
relevant to SNP-trait associations. VINTAGE also stands out as the only
effective method in assessing potential gene mediation effects, revealing
that most genes lack detectable mediation effects, which explains the power
advantage of VINTAGE and SKAT over TWAS.

</p>

## Installation

VINTAGE is implemented as an R package with underlying efficient C++
code interfaced through Rcpp and RcppArmadillo. VINTAGE depends on a few
other R packages. Please refer to the package
[DESCRIPTION](https://github.com/zhengli09/VINTAGE/blob/master/DESCRIPTION)
file for details. Dependent packages are supposed to be automatically
installed while installing VINTAGE. The installation of VINTAGE has been
tested on Ubuntu 20.04.6.

Install the VINTAGE R package maintained in github through the
`devtools` package. It takes a few minutes to install the package
depending on the number of dependencies that also need to be installed
in the process.

``` r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/VINTAGE")
library(VINTAGE)
?VINTAGE
?run_vintage
```

## How to use `VINTAGE` (v0.1.008)

Check our [vignettes](https://zhengli09.github.io/VINTAGE-analysis/).

## Citing the work

If you find the `VINTAGE` package or any of the source code in this
repository useful for your work, please cite:

> Li, Z., Gao, B, Zhou, X. VINTAGE: An alternative framework for
> transcriptome-wide association studies to detect and decipher
> gene-trait associations.

Visit our [group website](https://xiangzhou.github.io/) for more
statistical tools on analyzing genetics and genomics data.
