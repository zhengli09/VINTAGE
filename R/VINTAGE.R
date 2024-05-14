# Author: Zheng Li
# Date: 2023-05-16
# VINTAGE: A unified framework integrating gene expression mapping studies 
# with genome-wide association studies for detecting and deciphering gene-trait 
# associations

#' Run VINTAGE
#'
#' VINTAGE analyzes genes on one chromosome at a time
#' 
#' @param eqtl eQTL summary statistics
#' @param gwas GWAS summary statistics
#' @param refpath Reference panel for calculating SNP-SNP correlation
#' @param gene_info Gene annotations
#' @param outfile Output file name
#' @param ref_type Format of reference panel
#' @param ambiguity Whether to filter out strand ambiguous variants
#' @param maf Filter out variants with minor allele frequency below this value
#' @param window Integer. Buffer region size flanking the gene TSS and TES
#' @param B Number of simulations to evaluate the simulation-based p-value
#' @param dofilter Whether to filter out LD mismatched variants
#' @param pg_cutoff Perform genetic correlation test if the p-value of the
#'   genetic variance test passed pg_cutoff
#' @param ncores Number of cores to use
#' @param maxIterEM Maximum number of iterations for the EM algorithm
#' @param tolEM Tolerance for likelihood difference to claim convergence
#' @param save_profile Evaluate likelihood 
#' @param seed Random seed
run_vintage <- function(eqtl, gwas, refpath, gene_info, outfile, 
  ref_type = c("plink"), ambiguity = TRUE, maf = 0.05, window = 1e5, B = 1e6, 
  dofilter = TRUE, pg_cutoff = 0.05 / 20000, ncores = 1, maxIterEM = 1e5, 
  tolEM = 1e-3, save_profile = 1, seed = 0)
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
    message(sprintf(paste0("%.0f/%.0f(%.0f%%) variants excluded in LD ",
      "reference by MAF cutoff"), sum(!idx_keep), nrow(ldref$map), 
      sum(!idx_keep) / nrow(ldref$map) * 100))
    rdssub <- bigsnpr::snp_subset(ldref, ind.col = which(idx_keep))
    ldref <- bigsnpr::snp_attach(rdssub)
  }

  # 2.handle allele ambiguity
  if(ambiguity){
    cat("***** Exclude strand ambiguous variants (eQTL, GWAS, LD) *****\n")
    amb_idx <- .find_ambiguity(eqtl$effect_allele, eqtl$other_allele)
    eqtl <- eqtl[!amb_idx, ]
    amb_idx <- .find_ambiguity(gwas$effect_allele, gwas$other_allele)
    gwas <- gwas[!amb_idx, ]
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
  
  verbose <- pbmcapply::pbmclapply(1:nrow(gene_info), function(i){
    set.seed(seed + i)
    gene <- gene_info[i, "gene_id"]
    start <- gene_info[i, "start"]
    end <- gene_info[i, "end"]
    cat("- Handling gene", i, ": ", gene, "(", start, "-", end, ")\n", sep = "")
    
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
      cat("-- Filter out variants that have LD mismatch\n")
      filter1 <- susieR::kriging_rss(eqtl_matched$zscore, R, n1)$conditional_dist
      filter2 <- susieR::kriging_rss(gwas_matched$zscore, R, n2)$conditional_dist
      filter_idx <- filter1$logLR > 2 | filter2$logLR > 2 |
        abs(filter1$z_std_diff) > qnorm(0.995) | 
        abs(filter2$z_std_diff) > qnorm(0.995)
      message(sprintf("%.0f variants filtered due to LD mismatch", 
        sum(filter_idx)))
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
    null_g <- vintage(eqtl_matched$zscore_adj, gwas_matched$zscore_adj, svdR$u, 
      lambdas1, lambdas2, n, p, k, init_D, maxIterEM = maxIterEM, tolEM = tolEM, 
      save_profile = save_profile, scenario = "null", B = B)
    pw <- null_g$var_test$pvalues[, 1]
    pg <- ACAT::ACAT(pw)
    pg <- ifelse(pg == 0, 1 / B, pg)

    # under the alternative
    alt <- vintage(eqtl_matched$zscore_adj, gwas_matched$zscore_adj, svdR$u, 
      lambdas1, lambdas2, n, p, k, init_D, maxIterEM = maxIterEM, tolEM = tolEM, 
      save_profile = save_profile, scenario = "alt", B = B)

    # hypothesis testing for H0:r=0
    if(pg < pg_cutoff){
      null_r <- vintage(eqtl_matched$zscore_adj, gwas_matched$zscore_adj, 
        svdR$u, lambdas1, lambdas2, n, p, k, init_D, maxIterEM = maxIterEM, 
        tolEM = tolEM, save_profile = save_profile, scenario = "r0", B = B)
      pr <- pchisq(null_r$corr_test$Tr, df = 1, lower.tail = F)
    } else{
      pr <- NA
    }

    # # 3.6.run SKAT
    Q <- sum(gwas_matched$zscore_adj^2) / 2
    W <- svdR$u %*% diag(lambdas2) %*% t(svdR$v)
    skat <- SKAT::Get_Davies_PVal(Q, W)$p.value

    # 3.7.collect output
    out <- c(gene, start, end, window, p, k, s1, s2, h1_init, h2_init, 
      alt$h1_sq, alt$h2_sq, alt$r, alt$sigma1_sq, alt$sigma2_sq, pw, pg, 
      pr, skat, null_g$elapsed_time)
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
  "done"
}

.find_ambiguity <- function(allele1, allele2)
{
  amb_idx <- paste0(allele1, allele2) %in% c("AT", "TA", "CG", "GC")
  message(sprintf("%.0f/%.0f (%.0f%%) ambiguous variants identified",
    sum(amb_idx), length(amb_idx), sum(amb_idx) / length(amb_idx) * 100))
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
  message(sprintf("%.0f variants excluded, %.0f reversed, and %.0f flipped",
    nrow(sumstat) / 4 - nrow(matched), sum(matched$reversed), 
    sum(matched$flipped)))
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

#' Run MESuSiE to estimate the number of shared SNP signals
#' and the number of unique SNP signals
#' @noRd
.mesusie <- function(ss_eqtl, ss_gwas, R1, R2, L)
{
  ss_eqtl$Beta <- ss_eqtl$Z / sqrt(ss_eqtl$N)
  ss_gwas$Beta <- ss_gwas$Z / sqrt(ss_gwas$N)
  ss_eqtl$Se <- 1 / sqrt(ss_eqtl$N)
  ss_gwas$Se <- 1 / sqrt(ss_gwas$N)
  ss_list <- list(eqtl = ss_eqtl, gwas = ss_gwas)
  R_list <- list(eqtl = R1, gwas = R2)
  mesusie <- MESuSiE::meSuSie_core(R_list, ss_list, L = L)
  n_causal_eqtl <- sum(mesusie$cs$cs_category == "eqtl")
  n_causal_gwas <- sum(mesusie$cs$cs_category == "gwas")
  n_causal_both <- sum(mesusie$cs$cs_category == "eqtl_gwas")
  out <- c(n_causal_eqtl, n_causal_gwas, n_causal_both)
}

#' Run TWMR and revTWMR to check for reverse mediation effects
#' Implementation follows:
#' https://github.com/eleporcu/TWMR/blob/master/MR.R
#' https://github.com/eleporcu/revTWMR/blob/main/revTWMR.R
#' @noRd
.twmr <- function(ldref, ss1, ss2, do_clump = TRUE, constrain_size = TRUE)
{
  stopifnot(all(ss1$variant == ss2$variant))
  snps <- ss1$variant
  ss1$beta <- ss1$zscore / sqrt(ss1$N)
  ss2$beta <- ss2$zscore / sqrt(ss2$N)
  prop_valid <- NA
  
  if(constrain_size){
    # only use instruments with larger effects on the mediator than on outcome
    snps_constrain <- which(abs(ss1$beta) > abs(ss2$beta))
    prop_valid <- mean(abs(ss1$beta) > abs(ss2$beta))
    if(length(snps_constrain) == 0) return(c(NA, prop_valid))
    snps <- snps[snps_constrain]
    ss1 <- ss1[snps_constrain, ]
    ss2 <- ss2[snps_constrain, ]
  }
  if(do_clump){
    ind.col <- match(snps, ldref$map$marker.ID)
    ref_gene_path <- bigsnpr::snp_subset(ldref, ind.col = ind.col)
    ref_gene <- bigsnpr::snp_attach(ref_gene_path)
    infos.chr <- rep(1, ncol(ref_gene$genotypes))
    snps_clumped <- bigsnpr::snp_clumping(ref_gene$genotypes, infos.chr, 
      S = abs(ss1$zscore), thr.r2 = 0.01)
    file.remove(ref_gene_path)
    file.remove(gsub(".rds", ".bk", ref_gene_path))
    if(length(snps_clumped) == 0) return(c(NA, prop_valid))
    snps <- snps[snps_clumped]
    ss1 <- ss1[snps_clumped, ]
    ss2 <- ss2[snps_clumped, ]
  }
  alpha <- sum(ss1$beta * ss2$beta) / sum(ss1$beta * ss1$beta)
  c(alpha, prop_valid)
}

