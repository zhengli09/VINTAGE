# Author: Zheng Li
# Date: 2023-05-16

#' Run VINTAGE
#'
#' Perform gene-wise genetic variance tests and gene-wise genetic correlation 
#' tests for genes on one chromosome.
#' 
#' @param eqtl data.frame. eQTL summary statistics that minimally should contain
#'   the gene id, variant id, variant chromosome, variant genomic position, 
#'   effect allele, other allele, z scores, and sample size.
#' @param gwas data.frame. GWAS summary statistics that minimally should contain
#'   the variant id, variant chromosome, variant genomic position, effect allele, 
#'   other allele, z scores, and sample size.
#' @param refpath character. Path to an LD reference panel for calculating the 
#'   SNP-SNP correlations. If the reference panel is in plink format, please 
#'   specify the path to the ".bed" file.
#' @param gene_info data.frame. Gene annotations that minimally should contain
#'   the id, chromosome, transcription start site (TSS), and transcription end 
#'   site (TES) of each gene.
#' @param outfile character. Path to an output file.
#' @param ref_type character. Format of the LD reference panel file. As of 
#'   May 28, 2024, VINTAGE only supports plink file format.
#' @param ambiguity logical. Whether to filter out strand-ambiguous variants.
#' @param maf numeric. Filter out variants with minor allele frequencies below 
#'   this threshold in the LD reference panel.
#' @param window integer. Buffer region size flanking the gene TSS and TES.
#' @param B integer. Number of simulations/permutations used to evaluate p values
#'   in the gene-wise genetic variance test.
#' @param dofilter logical. Whether to filter out LD mis-matched variants
#'   between the LD reference panel and eQTL/GWAS summary statistics
#' @param pg_cutoff numeric. Perform gene-wise genetic correlation test if the 
#'   p value of the gene-wise genetic variance test passed this threshold.
#' @param ncores integer. Number of CPU cores to use.
#' @param maxIterEM integer. Maximum number of iterations for the EM algorithm.
#' @param tolEM numeric. Tolerance for likelihood difference to claim convergence.
#' @param save_profile integer. Evaluate likelihood every "save_profile" iterations.
#' @param seed integer. Random seed.
#' @param debug logical. Whether to show intermediate quantities for debugging.
#' 
#' @return The function does not have a return value. Instead, it directly 
#'   writes the results to the specified output file. The results consist of
#'   the following information for each gene:
#' \item{gene}{Gene id.}
#' \item{start}{Genomic position of gene TSS.}
#' \item{end}{Genomic position of gene TES.}
#' \item{window}{Buffer region size flanking the gene TSS and TES.}
#' \item{p}{Number of SNPs included in the analysis.}
#' \item{k}{Number of eigenvectors/eigenvalues used in the analysis.}
#' \item{s1}{Estimated shrinkage parameter between the LD reference panel and 
#'   eQTL summary statistics.}
#' \item{s2}{Estimated shrinkage parameter between the LD reference panel and
#'   GWAS summary statistics.}
#' \item{h1_init}{Initial estimate of gene expression cis-heritability.}
#' \item{h2_init}{Initial estimate of trait cis-heritability.}
#' \item{h1}{Final estimate of gene expression cis-heritability.}
#' \item{h2}{Final estimate of trait cis-heritability.}
#' \item{r}{Local genetic correlation estimate.}
#' \item{sigma1_sq}{Residual variance estimate from gene expression.}
#' \item{sigma2_sq}{Residual variance estimate from trait.}
#' \item{vin_var_pval00-10}{P values from each individual genetic variance test 
#'   constructed under a particular weight. The weight "00" corresponds to SKAT 
#'   whereas the weight "10" corresponds to TWAS. The other weights are sequential
#'   and are in-between these two extreme scenarios.}
#' \item{vin_var_test}{P value from VINTAGE's gene-wise genetic variance test
#'   which aggregates p values from all individual tests through the Cauchy
#'   combination approach.}
#' \item{vin_corr_test}{P value from VINTAGE's gene-wise genetic correlation
#'   test. An NA would be reported if the p value from the genetic variance
#'   test did pass the specified threshold, pg_cutoff.}
#' \item{skat}{P value from SKAT under an unweighted and linear kernel.}
#' \item{time}{Elapsed time (secs) for analyzing each gene.}
#' @export
#' 
run_vintage <- function(eqtl, gwas, refpath, gene_info, outfile, 
  ref_type = c("plink"), ambiguity = TRUE, maf = 0.05, window = 1e5, B = 1e6, 
  dofilter = TRUE, pg_cutoff = 0.05 / 20000, ncores = 1, maxIterEM = 1e5, 
  tolEM = 1e-3, save_profile = 1, seed = 0, debug = FALSE)
{
  ref_type <- match.arg(ref_type)
  
  # 1.load LD reference panel
  cat("***** Load LD reference data in", ref_type, "format *****\n")
  if(ref_type == "plink"){
    if(!file.exists(bigsnpr::sub_bed(refpath, ".bk"))){
      rds <- bigsnpr::snp_readBed(refpath)
    } else{
      rds <- bigsnpr::sub_bed(refpath, ".rds")
    }
    ldref <- bigsnpr::snp_attach(rds)
    ldref$map$marker.ID <- with(ldref$map, paste(chromosome, physical.pos,
      allele2, allele1, sep = ":"))
    # filter by minor allele frequency
    mafs <- bigsnpr::snp_MAF(ldref$genotypes)
    idx_keep <- mafs >= maf
    # message(sprintf(paste0("%.0f/%.0f(%.0f%%) variants excluded in LD ",
    #   "reference by MAF cutoff"), sum(!idx_keep), nrow(ldref$map), 
    #   sum(!idx_keep) / nrow(ldref$map) * 100))
    rdssub <- bigsnpr::snp_subset(ldref, ind.col = which(idx_keep))
    ldref <- bigsnpr::snp_attach(rdssub)
  }

  # 2.handle allele ambiguity
  if(ambiguity){
    cat("***** Exclude strand ambiguous variants (eQTL, GWAS, LD) *****\n")
    cat("(eQTL) ")
    amb_idx <- .find_ambiguity(eqtl$effect_allele, eqtl$other_allele)
    eqtl <- eqtl[!amb_idx, ]
    cat("(GWAS) ")
    amb_idx <- .find_ambiguity(gwas$effect_allele, gwas$other_allele)
    gwas <- gwas[!amb_idx, ]
    cat("(LD) ")
    amb_idx <- .find_ambiguity(ldref$map$allele1, ldref$map$allele2)
    if(ref_type == "plink"){
      rdssubsub <- bigsnpr::snp_subset(ldref, ind.col = which(!amb_idx))
      ldref <- bigsnpr::snp_attach(rdssubsub)
    }
  }
  
  # 3.analyze one gene at a time
  # 3.0.prepare output file
  header <- c("gene", "start", "end", "window", "p", "k", "s1", "s2", 
    "h1_init", "h2_init", "h1", "h2", "r", "sigma1_sq", "sigma2_sq",
    paste0("vin_var_pval", c("00", "01", "02", "03", "04", "05", "06", "07", 
      "08", "09", "10")), "vin_var_test", "vin_corr_test", "skat", "time")
  write.table(t(header), file = outfile, col.names = F, row.names = F,
    quote = F)
  
  verbose <- parallel::mclapply(1:nrow(gene_info), function(i){
    start_time <- Sys.time()
    set.seed(seed + i)
    gene <- gene_info[i, "gene_id"]
    start <- gene_info[i, "start"]
    end <- gene_info[i, "end"]
    cat("(", i, "/", nrow(gene_info), ") Handling gene: ", gene, " (", start, 
      "-", end, ")\n", sep = "")
    
    # 3.1.subset eQTL and GWAS data
    eqtl_gene <- eqtl[eqtl$gene == gene, , drop = F]
    eqtl_gene <- subset(eqtl_gene, pos >= (start - window) & 
        pos <= (end + window))
    gwas_gene <- subset(gwas, pos >= (start - window) & 
        pos <= (end + window))
    if(nrow(eqtl_gene) == 0 | nrow(gwas_gene) == 0){
      warning("No cis-SNPs remianed after merging datasets")
      return(NULL)
    }
    
    # 3.2.merge three datasets (use alleles in refld as reference)
    cat("  - Merge three datasets (LD, eQTL, and GWAS)\n")
    if(ref_type == "plink"){
      snpref <- subset(ldref$map, physical.pos >= (start - window) &
          physical.pos <= (end + window))[, "marker.ID", drop = F]
      colnames(snpref) <- "variant"
    }
    eqtl_matched <- .match_snp(snpref, eqtl_gene)
    if(nrow(eqtl_matched) == 0){
      warning("No cis-SNPs remianed after merging datasets")
      return(NULL)
    }
    snpref <- subset(snpref, variant %in% eqtl_matched$variant)
    gwas_matched <- .match_snp(snpref, gwas_gene)
    if(nrow(gwas_matched) == 0){
      warning("No cis-SNPs remianed after merging datasets")
      return(NULL)
    }
    snpref <- subset(snpref, variant %in% gwas_matched$variant)
    snpref <- snpref$variant
    eqtl_matched <- eqtl_matched[match(snpref, eqtl_matched$variant), , drop = F]
    gwas_matched <- gwas_matched[match(snpref, gwas_matched$variant), , drop = F]
    n1 <- median(eqtl_matched$N)
    n2 <- median(gwas_matched$N)
    n <- c(n1, n2)
    
    # 3.3.prepare model input
    eqtl_matched$zscore_adj <- with(eqtl_matched, zscore * sqrt((N - 1) /
        (N - 2 + zscore^2)))
    gwas_matched$zscore_adj <- with(gwas_matched, zscore * sqrt((N - 1) /
        (N - 2 + zscore^2)))
    if(ref_type == "plink"){
      R <- cor(ldref$genotypes[, match(snpref, ldref$map$marker.ID), drop = F])
    }
    if(dofilter){
      cat("  - Filter out variants that have LD mismatch\n")
      filter1 <- susieR::kriging_rss(eqtl_matched$zscore, R, n1)$conditional_dist
      filter2 <- susieR::kriging_rss(gwas_matched$zscore, R, n2)$conditional_dist
      filter_idx <- filter1$logLR > 2 | filter2$logLR > 2 |
        abs(filter1$z_std_diff) > qnorm(0.995) | 
        abs(filter2$z_std_diff) > qnorm(0.995)
      # message(sprintf("    %.0f%% variants filtered due to LD mismatch", 
      #   signif(mean(filter_idx), digits = 4) * 100))
      eqtl_matched <- eqtl_matched[!filter_idx, , drop = F]
      gwas_matched <- gwas_matched[!filter_idx, , drop = F]
      R <- R[!filter_idx, !filter_idx, drop = F]
    }
    colnames(R) <- rownames(R) <- eqtl_matched$variant
    p <- nrow(eqtl_matched)
    if(p == 0){
      warning("No cis-SNPs remianed after LD mismatch filtering")
      return(NULL)
    }
    # get around LAPACK error code 1 with eigen
    svdR <- tryCatch({
      svd(R)
    }, error = function(e){
      eigR <- eigen(R)
      list(u = eigR$vectors, d = eigR$values, v = eigR$vectors)
    })
    k <- which.min(abs(cumsum(svdR$d) / sum(svdR$d) - 0.99))
    s1 <- susieR::estimate_s_rss(eqtl_matched$zscore, R, n1)
    s2 <- susieR::estimate_s_rss(gwas_matched$zscore, R, n2)
    lambdas1 <- (1 - s1) * svdR$d + s1
    lambdas2 <- (1 - s2) * svdR$d + s2
    
    # 3.4.initial values
    h1_init <- .init_h2(t(svdR$u) %*% eqtl_matched$zscore_adj, lambdas1, n1, p, k)
    h2_init <- .init_h2(t(svdR$u) %*% gwas_matched$zscore_adj, lambdas2, n2, p, k)
    init_D <- matrix(c(h1_init / p, 0, 0, h2_init / p), 2, 2)
    
    # 3.5.run VINTAGE
    # under the null (sigma_beta2_sq = rho = 0)
    null_g <- vintage(z1 = eqtl_matched$zscore_adj, z2 = gwas_matched$zscore_adj, 
      Q = svdR$u, lambdas1 = lambdas1, lambdas2 = lambdas2, n = n, p = p, k = k, 
      init_D = init_D, maxIterEM = maxIterEM, tolEM = tolEM, scenario = "null", 
      B = B, save_profile = save_profile, debug = debug)
    pw <- null_g$var_test$pvalues[, 1]
    pg <- ACAT::ACAT(pw)
    pg <- ifelse(pg == 0, 1 / B, pg)

    # under the alternative
    alt <- vintage(z1 = eqtl_matched$zscore_adj, z2 = gwas_matched$zscore_adj, 
      Q = svdR$u, lambdas1 = lambdas1, lambdas2 = lambdas2, n = n, p = p, k = k, 
      init_D = init_D, maxIterEM = maxIterEM, tolEM = tolEM, scenario = "alt", 
      B = B, save_profile = save_profile, debug = debug)

    # hypothesis testing for H0:r=0
    if(pg < pg_cutoff){
      null_r <- vintage(z1 = eqtl_matched$zscore_adj, z2 = gwas_matched$zscore_adj, 
        Q = svdR$u, lambdas1 = lambdas1, lambdas2 = lambdas2, n = n, p = p, k = k, 
        init_D = init_D, maxIterEM = maxIterEM, tolEM = tolEM, scenario = "r0", 
        B = B, save_profile = save_profile, debug = debug)
      pr <- pchisq(null_r$corr_test$Tr, df = 1, lower.tail = F)
    } else{
      pr <- NA
    }

    # # 3.6.run SKAT
    Q <- sum(gwas_matched$zscore_adj^2) / 2
    W <- svdR$u %*% diag(lambdas2) %*% t(svdR$v)
    skat <- SKAT::Get_Davies_PVal(Q, W)$p.value

    # 3.7.collect output
    end_time <- Sys.time()
    elapsed_time <- as.double(end_time-start_time, units = "secs")
    out <- c(gene, start, end, window, p, k, s1, s2, h1_init, h2_init, 
      alt$h1_sq, alt$h2_sq, alt$r, alt$sigma1_sq, alt$sigma2_sq, pw, pg, 
      pr, skat, elapsed_time)
    if(!null_g$fail & !alt$fail){
      write.table(t(out), file = outfile, col.names = F, row.names = F, 
        quote = F, append = T) 
    }
  }, mc.cores = ncores)
  
  # 4.remove intermediate files
  if(ref_type == "plink"){
    file.remove(rdssub)
    file.remove(gsub(".rds", ".bk", rdssub))
    if(ambiguity){
      file.remove(rdssubsub)
      file.remove(gsub(".rds", ".bk", rdssubsub))
    }
  }
  
  # "done"
  return(invisible(NULL))
}

.find_ambiguity <- function(allele1, allele2)
{
  amb_idx <- paste0(allele1, allele2) %in% c("AT", "TA", "CG", "GC")
  cat(sprintf("%.0f/%.0f (%.0f%%) ambiguous variants identified",
    sum(amb_idx), length(amb_idx), sum(amb_idx) / length(amb_idx) * 100),
    "\n", sep = "")
  amb_idx
}

.flip_strand <- function(allele)
{
  dplyr::case_when(
    allele == "A" ~ "T",
    allele == "C" ~ "G",
    allele == "T" ~ "A",
    allele == "G" ~ "C",
    TRUE ~ NA_character_
  )
}

.match_snp <- function(snpref, sumstat)
{
  sumstat$reversed <- FALSE
  sumstat$flipped <- FALSE
  
  # reversed allele order
  snprev <- sumstat
  snprev$variant <- with(sumstat, paste(chr, pos, 
    effect_allele, other_allele, sep = ":"))
  snprev$effect_allele <- sumstat$other_allele
  snprev$other_allele <- sumstat$effect_allele
  snprev$zscore <- -snprev$zscore
  snprev$reversed <- TRUE
  sumstat <- rbind(sumstat, snprev)
  
  # flipped alleles
  snpflip <- sumstat
  snpflip$effect_allele <- .flip_strand(snpflip$effect_allele)
  snpflip$other_allele <- .flip_strand(snpflip$other_allele)
  snpflip$variant <- with(snpflip, paste(chr, pos, 
    other_allele, effect_allele, sep = ":"))
  snpflip$flipped <- TRUE
  sumstat <- rbind(sumstat, snpflip)
  
  matched <- merge(sumstat, snpref, by = "variant")
  # message(sprintf("%.0f variants excluded, %.0f reversed, and %.0f flipped",
  #   nrow(sumstat) / 4 - nrow(matched), sum(matched$reversed), 
  #   sum(matched$flipped)))
  matched <- subset(matched, select = -c(reversed, flipped))
}

#' Find an initial value for the heritability
#' @noRd
.init_h2 <- function(Qtz, lambdas, n, p, k)
{
  Qtz <- Qtz[1:k]
  lambdas <- lambdas[1:k]
  neglogllk <- function(sigma_beta2_sq, Qtz, lambdas, n, k)
  {
    logllk <- -k * log(sigma_beta2_sq) - 
      sum(log(n * lambdas + 1.0 / sigma_beta2_sq)) +
      n * sum(Qtz * Qtz / (n * lambdas + 1.0 / sigma_beta2_sq))
    -logllk
  }
  opt_res <- optim(par = 0, fn = neglogllk, method = "Brent",
    Qtz = Qtz, lambdas = lambdas, n = n, k = k, lower = 0, upper = 1)
  opt_res$par * p
}

