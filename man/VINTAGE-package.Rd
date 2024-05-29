\name{VINTAGE-package}
\alias{VINTAGE-package}
\alias{VINTAGE}
\docType{package}
\title{
  VINTAGE: Variant-set test INtegrative TWAS for GEne-based analysis
}
\description{
  This package implements the VINTAGE method for an integrative analysis of
  gene expression mapping studies and genome-wide association studies to 
  identify and decipher gene-trait associations. VINTAGE consists of two tests:
  a gene-wise genetic variance test that examines the joint effects of all
  genetic variants in a gene region on a trait of interest; and a gene-wise
  genetic correlation test that examines the potential mediation effects of 
  gene expression on the trait.
}
\details{
  Refer to the core function "run_vintage" for performing the integrative
  analysis of a gene expression study and GWAS. The function streamlines the 
  analysis for one chromosome at a time and takes as input the eQTL summary
  statistics, GWAS summary statistics, and an LD reference panel. The main
  outputs are p values from the gene-wise genetic variance test and gene-wise
  genetic correlation test.
}
\section{Package Content}{\packageIndices{VINTAGE}}
\author{
  Zheng Li, Boran Gao, Xiang Zhou
}
\section{Maintainer}{
  Zheng Li <zlisph@umich.edu>
}
\references{
  VINTAGE: A unified framework integrating gene expression mapping studies with
  genome-wide association studies for detecting and deciphering gene-trait 
  associations
}
\keyword{package}
