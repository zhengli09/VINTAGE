% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VINTAGE.R
\name{run_vintage}
\alias{run_vintage}
\title{Run VINTAGE}
\usage{
run_vintage(
  eqtl,
  gwas,
  refpath,
  gene_info,
  outfile,
  ref_type = c("plink"),
  ambiguity = TRUE,
  maf = 0.05,
  window = 1e+05,
  B = 1e+06,
  dofilter = TRUE,
  pg_cutoff = 0.05/20000,
  ncores = 1,
  maxIterEM = 1e+05,
  tolEM = 0.001,
  save_profile = 1,
  seed = 0,
  debug = FALSE
)
}
\arguments{
\item{eqtl}{data.frame. eQTL summary statistics that minimally should contain
the gene id, variant id, variant chromosome, variant genomic position, 
effect allele, other allele, z scores, and sample size.}

\item{gwas}{data.frame. GWAS summary statistics that minimally should contain
the variant id, variant chromosome, variant genomic position, effect allele, 
other allele, z scores, and sample size.}

\item{refpath}{character. Path to an LD reference panel for calculating the 
SNP-SNP correlations. If the reference panel is in plink format, please 
specify the path to the ".bed" file.}

\item{gene_info}{data.frame. Gene annotations that minimally should contain
the id, chromosome, transcription start site (TSS), and transcription end 
site (TES) of each gene.}

\item{outfile}{character. Path to an output file.}

\item{ref_type}{character. Format of the LD reference panel file. As of 
May 28, 2024, VINTAGE only supports plink file format.}

\item{ambiguity}{logical. Whether to filter out strand-ambiguous variants.}

\item{maf}{numeric. Filter out variants with minor allele frequencies below 
this threshold in the LD reference panel.}

\item{window}{integer. Buffer region size flanking the gene TSS and TES.}

\item{B}{integer. Number of simulations/permutations used to evaluate p values
in the gene-wise genetic variance test.}

\item{dofilter}{logical. Whether to filter out LD mis-matched variants
between the LD reference panel and eQTL/GWAS summary statistics}

\item{pg_cutoff}{numeric. Perform gene-wise genetic correlation test if the 
p value of the gene-wise genetic variance test passed this threshold.}

\item{ncores}{integer. Number of CPU cores to use.}

\item{maxIterEM}{integer. Maximum number of iterations for the EM algorithm.}

\item{tolEM}{numeric. Tolerance for likelihood difference to claim convergence.}

\item{save_profile}{integer. Evaluate likelihood every "save_profile" iterations.}

\item{seed}{integer. Random seed.}

\item{debug}{logical. Whether to show intermediate quantities for debugging.}
}
\value{
The function does not have a return value. Instead, it directly 
  writes the results to the specified output file. The results consist of
  the following information for each gene:
\item{gene}{Gene id.}
\item{start}{Genomic position of gene TSS.}
\item{end}{Genomic position of gene TES.}
\item{window}{Buffer region size flanking the gene TSS and TES.}
\item{p}{Number of SNPs included in the analysis.}
\item{k}{Number of eigenvectors/eigenvalues used in the analysis.}
\item{s1}{Estimated shrinkage parameter between the LD reference panel and 
  eQTL summary statistics.}
\item{s2}{Estimated shrinkage parameter between the LD reference panel and
  GWAS summary statistics.}
\item{h1_init}{Initial estimate of gene expression cis-heritability.}
\item{h2_init}{Initial estimate of trait cis-heritability.}
\item{h1}{Final estimate of gene expression cis-heritability.}
\item{h2}{Final estimate of trait cis-heritability.}
\item{r}{Local genetic correlation estimate.}
\item{sigma1_sq}{Residual variance estimate from gene expression.}
\item{sigma2_sq}{Residual variance estimate from trait.}
\item{vin_var_pval00-10}{P values from each individual genetic variance test 
  constructed under a particular weight. The weight "00" corresponds to SKAT 
  whereas the weight "10" corresponds to TWAS. The other weights are sequential
  and are in-between these two extreme scenarios.}
\item{vin_var_test}{P value from VINTAGE's gene-wise genetic variance test
  which aggregates p values from all individual tests through the Cauchy
  combination approach.}
\item{vin_corr_test}{P value from VINTAGE's gene-wise genetic correlation
  test. An NA would be reported if the p value from the genetic variance
  test did pass the specified threshold, pg_cutoff.}
\item{skat}{P value from SKAT under an unweighted and linear kernel.}
\item{time}{Elapsed time (secs) for analyzing each gene.}
}
\description{
Perform gene-wise genetic variance tests and gene-wise genetic correlation 
tests for genes on one chromosome.
}
