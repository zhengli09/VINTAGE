# Author: Zheng Li
# Date: 2023-05-16

#' Run VINTAGE
#'
#' VINTAGE analyzes genes on one chromosome at a time
#' 
#' @param eqtl eQTL summary statistics
#' @param gwas GWAS summary statistics
#' @param refpath Reference panel for calculating SNP-SNP correlation
#' @param gene_info Gene annotations
#' @param outfile Ouput file name
#' @param ref_type Format of reference panel
#' @param ambiguity Whether to filter out strand ambiguous variants
#' @param maf Filter out variants with minor allele frequency below this value
#' @param window Integer. Buffer region size flanking the gene TSS and TES
#' @param B Number of simulations to evaluate the simulation-based p-value
#' @param dofilter Whether to filter out LD mismatched variants
#' @param ncores Number of cores to use
#' @param maxIterEM Maximum number of iterations for the EM algorithm
#' @param tolEM Tolerance for likelihood difference to claim convergence
#' @param save_profile Evaluate likelihood 
#' @param seed Random seed
run_vintage <- function(eqtl, gwas, refpath, gene_info, outfile, 
  ref_type = c("plink", "vcf"), ambiguity = TRUE, maf = 0.05, 
  window = 50000, B = 1e7, dofilter = TRUE, ncores = 1,
  maxIterEM = 1e5, tolEM = 1e-5, save_profile = 1, seed = 0)
{
  set.seed(seed)
  ref_type <- match.arg(ref_type)
  
  # 1.load LD reference panel
  cat("***** Load LD reference data in", ref_type, "format *****\n")
  if(ref_type == "plink"){
    if(!file.exists(paste0(bigsnpr::sub_bed(refpath), ".bk"))){
      rds <- bigsnpr::snp_readBed(refpath)
    } else{
      rds <- paste0(bigsnpr::sub_bed(refpath), ".rds")
    }
    ldref <- bigsnpr::snp_attach(rds)
    ldref$map$marker.ID <- with(ldref$map, paste(chromosome, physical.pos,
      allele2, allele1, sep = ":"))
    # filter by minor allele frequency
    mafs <- bigsnpr::snp_MAF(ldref$genotypes)
    idx_keep <- mafs >= maf
    message(sprintf(paste0("%.0f/%.0f(%.0f%%) variants excluded in LD ",
      "reference by MAF cutoff"), sum(!idx_keep), nrow(ldref$map), 
      sum(!idx_keep) / nrow(ldref$map) * 100))
    rdssub <- bigsnpr::snp_subset(ldref, ind.col = which(idx_keep))
    ldref <- bigsnpr::snp_attach(rdssub)
  }
  
  # 2.handle allele ambiguity
  if(ambiguity){
    cat("***** Exclude strand ambiguous variants (eQTL, GWAS, LD) *****\n")
    amb_idx <- find_ambiguity(eqtl$effect_allele, eqtl$other_allele)
    eqtl <- eqtl[!amb_idx, ]
    amb_idx <- find_ambiguity(gwas$effect_allele, gwas$other_allele)
    gwas <- gwas[!amb_idx, ]
    amb_idx <- find_ambiguity(ldref$map$allele1, ldref$map$allele2)
    if(ref_type == "plink"){
      rdssub <- bigsnpr::snp_subset(ldref, ind.col = which(!amb_idx))
      ldref <- bigsnpr::snp_attach(rdssub)
    }
  }
  
  # 3.analyze one gene at a time
  # 3.0.prepare output file
  header <- c("gene", "start", "end", "window", "p", "s1", "s2", "h1_init",
    "h2_init", "h1", "h2", "r", "sigma1_sq", "sigma2_sq", "null_converged",
    "null_fix_D", "alt_converged", "alt_fix_D", "pval10", "pval91", "pval82",
    "pval73", "pval64", "pval55", "pval46", "pval37", "pval28", "pval19", 
    "pval01", "VINTAGE", "skat", "time")
  write.table(t(header), file = outfile, col.names = F, row.names = F, 
    quote = F)
  
  verbose <- pbmcapply::pbmclapply(1:nrow(gene_info), function(i){
    gene <- gene_info[i, "gene_id"]
    start <- gene_info[i, "start"]
    end <- gene_info[i, "end"]
    cat("- Handling gene", i, ": ", gene, "(", start, "-", end, ")\n", sep = "")
    
    # 3.1.subset eQTL and GWAS data
    eqtl_gene <- eqtl[eqtl$gene == gene, ]
    eqtl_gene <- subset(eqtl_gene, pos >= (start - window) & 
        pos <= (end + window))
    gwas_gene <- subset(gwas, pos >= (start - window) & 
        pos <= (end + window))
    if(nrow(eqtl_gene) == 0 | nrow(gwas_gene) == 0){
      warning("No cis-SNPs remianed after merging datasets")
      return(NULL)
    }
    
    # 3.2.merge three datasets (use alleles in refld as reference)
    if(ref_type == "plink"){
      snpref <- ldref$map[, "marker.ID", drop = F]
      colnames(snpref) <- "variant"
    }
    eqtl_matched <- match_snp(snpref, eqtl_gene)
    if(nrow(eqtl_matched) == 0){
      warning("No cis-SNPs remianed after merging datasets")
      return(NULL)
    }
    snpref <- subset(snpref, variant %in% eqtl_matched$variant)
    gwas_matched <- match_snp(snpref, gwas_gene)
    if(nrow(gwas_matched) == 0){
      warning("No cis-SNPs remianed after merging datasets")
      return(NULL)
    }
    snpref <- subset(snpref, variant %in% gwas_matched$variant)
    snpref <- snpref$variant
    eqtl_matched <- eqtl_matched[match(snpref, eqtl_matched$variant), ]
    gwas_matched <- gwas_matched[match(snpref, gwas_matched$variant), ]
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
      cat("-- Filter out variants that have LD mismatch\n")
      filter1 <- susieR::kriging_rss(eqtl_matched$zscore, R, n1)
      filter2 <- susieR::kriging_rss(gwas_matched$zscore, R, n2)
      filter_idx <- filter1$conditional_dist$logLR > 2 | 
        filter2$conditional_dist$logLR > 2
      message(sprintf("%.0f variants filtered due to LD mismatch", 
        sum(filter_idx)))
      eqtl_matched <- eqtl_matched[!filter_idx, ]
      gwas_matched <- gwas_matched[!filter_idx, ]
      R <- R[!filter_idx, !filter_idx]
    }
    p <- nrow(eqtl_matched)
    if(p == 0){
      warning("No cis-SNPs remianed after LD mismatch filtering")
      return(NULL)
    }
    svdR <- svd(R)
    s1 <- susieR::estimate_s_rss(eqtl_matched$zscore, R, n1, 
      method = "null-pseudomle")
    s2 <- susieR::estimate_s_rss(gwas_matched$zscore, R, n2,
      method = "null-pseudomle")
    lambdas1 <- (1 - s1) * svdR$d + s1
    lambdas2 <- (1 - s2) * svdR$d + s2
    
    # 3.4.initial values
    l2_h1 <- apply((svdR$u %*% diag(lambdas1) %*% t(svdR$v))^2, 1, sum)
    l2_h2 <- apply((svdR$u %*% diag(lambdas2) %*% t(svdR$v))^2, 1, sum)
    h1_init <- tryCatch({bigsnpr::snp_ldsc(ld_score = l2_h1, ld_size = p, 
      chi2 = eqtl_matched$zscore_adj^2, sample_size = n1, blocks = NULL)["h2"]
    }, error = function(e) 0.01)
    h2_init <- tryCatch({bigsnpr::snp_ldsc(ld_score = l2_h2, ld_size = p, 
      chi2 = gwas_matched$zscore_adj^2, sample_size = n2, blocks = NULL)["h2"]
    }, error = function(e) 0.0001)
    if(h1_init < 0 | h1_init > 1) h1_init <- 0.01
    if(h2_init < 0 | h2_init > 1) h2_init <- 0.0001
    init_D <- matrix(c(h1_init, 0, 0, h2_init), 2, 2)
    
    # 3.5.run VINTAGE
    # under the null (sigma_beta2_sq = rho = 0)
    null <- vintage(eqtl_matched$zscore_adj, gwas_matched$zscore_adj, svdR$u, 
      lambdas1, lambdas2, n, init_D, maxIterEM = maxIterEM, tolEM = tolEM, 
      save_profile = save_profile, scenario = "null", B = B)
    wts_pvalues <- null$test_h2$pvalues[, 1]
    pvalue <- ACAT::ACAT(wts_pvalues)
    # under the alternative
    alt <- vintage(eqtl_matched$zscore_adj, gwas_matched$zscore_adj, svdR$u, 
      lambdas1, lambdas2, n, init_D, maxIterEM = maxIterEM, tolEM = tolEM, 
      save_profile = save_profile, scenario = "alt", B = B)

    # 3.6.run SKAT
    skat <- CompQuadForm::davies(sum(gwas_matched$zscore_adj^2), 
      lambdas2, acc = 1e-10)$Qq
    skat <- ifelse(skat > 1, 1, skat)
    
    # 3.7.collect output
    out <- c(gene, start, end, window, p, s1, s2, h1_init, h2_init, alt$h1_sq,
      alt$h2_sq, alt$r, alt$sigma1_sq, alt$sigma2_sq, null$is_converged,
      null$is_D_fixed, alt$is_converged, alt$is_D_fixed, wts_pvalues, pvalue,
      skat, null$elapsed_time)
    write.table(t(out), file = outfile, col.names = F, row.names = F, 
      quote = F, append = T)
  }, mc.cores = ncores)
  
  # 5.remove intermediate files
  files_remove <- list.files(dirname(refpath), pattern = paste0(
    bigsnpr::sub_bed(basename(refpath)), ".+rds|", 
    bigsnpr::sub_bed(basename(refpath)), ".+bk"), full.names = T)
  file.remove(files_remove)
  "done"
}

find_ambiguity <- function(allele1, allele2)
{
  amb_idx <- paste0(allele1, allele2) %in% c("AT", "TA", "CG", "GC")
  message(sprintf("%.0f/%.0f (%.0f%%) ambiguous variants identified",
    sum(amb_idx), length(amb_idx), sum(amb_idx) / length(amb_idx) * 100))
  amb_idx
}

flip_strand <- function(allele)
{
  dplyr::case_when(
    allele == "A" ~ "T",
    allele == "C" ~ "G",
    allele == "T" ~ "A",
    allele == "G" ~ "C",
    TRUE ~ NA_character_
  )
}

match_snp <- function(snpref, sumstat)
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
  snpflip$effect_allele <- flip_strand(snpflip$effect_allele)
  snpflip$other_allele <- flip_strand(snpflip$other_allele)
  snpflip$variant <- with(snpflip, paste(chr, pos, 
    other_allele, effect_allele, sep = ":"))
  snpflip$flipped <- TRUE
  sumstat <- rbind(sumstat, snpflip)
  
  matched <- merge(sumstat, snpref, by = "variant")
  message(sprintf("%.0f variants excluded, %.0f reversed, and %.0f flipped",
    nrow(sumstat) / 4 - nrow(matched), sum(matched$reversed), 
    sum(matched$flipped)))
  matched <- subset(matched, select = -c(reversed, flipped))
}
