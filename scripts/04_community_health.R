#!/usr/bin/env Rscript
# =============================================================================
# 04_community_health.R
# Community health analysis for VEG-05 microbiome (OSD-766)
#
# Three pillars of community health:
#   1. Diversity: alpha (Observed, Shannon, Simpson) + beta (Bray-Curtis, Jaccard)
#      with PERMANOVA for flight × light effects
#   2. Dysbiosis index: normalized distance to ground-control centroid
#   3. Differential abundance: ALDEx2 (replacing ANCOMBC) for flight and light effects
#
# Preprocessing:
#   - 16S: Remove chloroplast + mitochondria ASVs (host contamination)
#   - ITS: Use as-is (clean fungal community)
#   - Keep all samples (per user decision on low-depth 16S leaf)
#
# Input:  results/microbiome/{16S,ITS}/phyloseq_{amplicon}.rds
#         data/metadata/sample_metadata_microbiome.csv
#
# Output: results/microbiome/community_health/
#           alpha_diversity_{amplicon}.tsv
#           beta_diversity_permanova_{amplicon}.tsv
#           dysbiosis_index_{amplicon}.tsv
#           aldex2_results_{amplicon}_{contrast}.tsv
#           phyloseq_16S_filtered.rds  (chloroplast/mitochondria removed)
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ALDEx2)
  library(ggplot2)
  library(dplyr)
})

# --- Configuration ----------------------------------------------------------
REPO_ROOT  <- "/mnt/shared-workspace/veg05-integrated-omics"
OUT_DIR    <- file.path(REPO_ROOT, "results/microbiome/community_health")
PDF_DIR    <- "/workspace/community_pdfs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

# Compartments for community health analysis
COMMUNITY_COMPS <- c("leaf", "AdvRoot", "root", "wick", "soil")
PERMANOVA_PERMS <- 999
ALDEX_MC <- 128  # Monte Carlo instances for ALDEx2

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# --- Helper: load and preprocess phyloseq objects ---------------------------
load_phyloseq <- function(amplicon) {
  ps <- readRDS(file.path(REPO_ROOT, "results/microbiome", amplicon,
                           paste0("phyloseq_", amplicon, ".rds")))
  
  if (amplicon == "16S") {
    # Remove chloroplast and mitochondria (host contamination)
    ps <- subset_taxa(ps, !grepl("Chloroplast", Order) & !grepl("Mitochondria", Family))
    log_msg("16S: After removing chloroplast/mitochondria: ", ntaxa(ps), " ASVs")
    
    # Save filtered phyloseq for downstream use
    saveRDS(ps, file.path(REPO_ROOT, "results/microbiome/16S/phyloseq_16S_filtered.rds"))
  }
  
  # Remove samples with 0 reads
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  log_msg(amplicon, ": ", nsamples(ps), " samples, ", ntaxa(ps), " ASVs, ",
          sum(otu_table(ps)), " total reads")
  return(ps)
}

# --- 1. Alpha diversity -----------------------------------------------------
compute_alpha <- function(ps, amplicon) {
  log_msg("Computing alpha diversity for ", amplicon, "...")
  
  # Observed ASV richness
  # Note: estimate_richness converts sample names via make.names (hyphens -> dots)
  # Use sample_names() to get original names
  rich <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
  rich$sample_id <- sample_names(ps)
  
  # Add metadata
  sd <- data.frame(sample_data(ps))
  sd$sample_id <- rownames(sd)
  alpha <- merge(rich, sd, by = "sample_id")
  
  # Add read depth
  alpha$read_depth <- sample_sums(ps)
  
  # Save
  write.table(alpha, file.path(OUT_DIR, paste0("alpha_diversity_", amplicon, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Summary by compartment × flight × light (for integration compartments)
  log_msg("  Alpha diversity summary (integration compartments):")
  for (comp in c("leaf", "AdvRoot")) {
    sub <- alpha[alpha$compartment == comp, ]
    if (nrow(sub) > 0) {
      log_msg("    ", comp, " (n=", nrow(sub), "):")
      for (metric in c("Observed", "Shannon")) {
        vals <- tapply(sub[[metric]], paste(sub$flight, sub$light, sep = "_"), 
                       function(x) sprintf("%.1f ± %.1f", mean(x), sd(x)))
        log_msg("      ", metric, ": ", paste(names(vals), vals, sep = "=", collapse = ", "))
      }
    }
  }
  
  return(alpha)
}

# --- 2. Beta diversity + PERMANOVA ------------------------------------------
compute_beta <- function(ps, amplicon) {
  log_msg("Computing beta diversity for ", amplicon, "...")
  
  # Bray-Curtis dissimilarity
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  
  # Relative abundance for Bray-Curtis
  otu_rel <- sweep(otu, 1, rowSums(otu), "/")
  otu_rel[is.na(otu_rel)] <- 0
  
  bray_dist <- vegdist(otu_rel, method = "bray")
  
  sd <- data.frame(sample_data(ps))
  
  # PERMANOVA: flight × light, stratified by compartment
  permanova_results <- list()
  
  for (comp in COMMUNITY_COMPS) {
    comp_samples <- rownames(sd)[sd$compartment == comp]
    if (length(comp_samples) < 4) next
    
    # Subset distance matrix
    comp_dist <- as.dist(as.matrix(bray_dist)[comp_samples, comp_samples])
    comp_sd <- sd[comp_samples, ]
    
    # Check if we have both flight levels
    if (length(unique(comp_sd$flight)) < 2) next
    
    # PERMANOVA with flight and light (if available)
    if (length(unique(comp_sd$light)) > 1 && !all(is.na(comp_sd$light))) {
      # Replace NA light with "NA"
      comp_sd$light[is.na(comp_sd$light)] <- "NA"
      set.seed(42)
      perm <- adonis2(comp_dist ~ flight * light, data = comp_sd,
                      permutations = PERMANOVA_PERMS)
    } else {
      set.seed(42)
      perm <- adonis2(comp_dist ~ flight, data = comp_sd,
                      permutations = PERMANOVA_PERMS)
    }
    
    # Extract results
    for (term in rownames(perm)) {
      permanova_results[[length(permanova_results) + 1]] <- data.frame(
        amplicon = amplicon,
        compartment = comp,
        term = term,
        df = perm[term, "Df"],
        sum_sq = perm[term, "SumOfSqs"],
        r2 = perm[term, "R2"],
        f_stat = perm[term, "F"],
        p_value = perm[term, "Pr(>F)"],
        n_samples = length(comp_samples)
      )
    }
  }
  
  permanova_df <- do.call(rbind, permanova_results)
  write.table(permanova_df, file.path(OUT_DIR, paste0("beta_diversity_permanova_", amplicon, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  log_msg("  PERMANOVA results saved (", nrow(permanova_df), " terms)")
  
  # Print key results for integration compartments
  for (comp in c("leaf", "AdvRoot")) {
    comp_res <- permanova_df[permanova_df$compartment == comp, ]
    if (nrow(comp_res) > 0) {
      log_msg("    ", comp, " PERMANOVA:")
      for (i in 1:nrow(comp_res)) {
        r <- comp_res[i, ]
        log_msg("      ", r$term, ": R²=", round(r$r2, 3), 
                ", F=", round(r$f_stat, 2), ", p=", signif(r$p_value, 3))
      }
    }
  }
  
  return(permanova_df)
}

# --- 3. Dysbiosis index -----------------------------------------------------
compute_dysbiosis <- function(ps, amplicon) {
  log_msg("Computing dysbiosis index for ", amplicon, "...")
  
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  
  # Relative abundance
  otu_rel <- sweep(otu, 1, rowSums(otu), "/")
  otu_rel[is.na(otu_rel)] <- 0
  
  bray_dist <- vegdist(otu_rel, method = "bray")
  dist_mat <- as.matrix(bray_dist)
  
  sd <- data.frame(sample_data(ps))
  
  dysbiosis_results <- list()
  
  for (comp in COMMUNITY_COMPS) {
    comp_samples <- rownames(sd)[sd$compartment == comp]
    if (length(comp_samples) < 4) next
    
    comp_sd <- sd[comp_samples, ]
    
    # Find ground control samples
    ground_samples <- comp_samples[comp_sd$flight == "Ground"]
    flight_samples <- comp_samples[comp_sd$flight == "Flight"]
    
    if (length(ground_samples) < 2 || length(flight_samples) < 1) next
    
    # Ground centroid: mean pairwise distance among ground samples
    ground_dist <- dist_mat[ground_samples, ground_samples]
    ground_within <- mean(ground_dist[upper.tri(ground_dist)])
    
    if (ground_within == 0) ground_within <- NA  # Avoid division by zero
    
    # For each flight sample: distance to ground centroid
    # (mean distance to all ground samples)
    for (fs in flight_samples) {
      dist_to_ground <- mean(dist_mat[fs, ground_samples])
      dysbiosis_score <- dist_to_ground / ground_within
      
      dysbiosis_results[[length(dysbiosis_results) + 1]] <- data.frame(
        amplicon = amplicon,
        sample_id = fs,
        compartment = comp,
        flight = "Flight",
        dist_to_ground_centroid = dist_to_ground,
        ground_within_group_dist = ground_within,
        dysbiosis_index = dysbiosis_score,
        light = comp_sd[fs, "light"]
      )
    }
    
    # Also compute for ground samples (distance to ground centroid, should be ~1)
    for (gs in ground_samples) {
      other_ground <- setdiff(ground_samples, gs)
      if (length(other_ground) < 1) next
      dist_to_centroid <- mean(dist_mat[gs, other_ground])
      dysbiosis_score <- dist_to_centroid / ground_within
      
      dysbiosis_results[[length(dysbiosis_results) + 1]] <- data.frame(
        amplicon = amplicon,
        sample_id = gs,
        compartment = comp,
        flight = "Ground",
        dist_to_ground_centroid = dist_to_centroid,
        ground_within_group_dist = ground_within,
        dysbiosis_index = dysbiosis_score,
        light = comp_sd[gs, "light"]
      )
    }
  }
  
  dysbiosis_df <- do.call(rbind, dysbiosis_results)
  write.table(dysbiosis_df, file.path(OUT_DIR, paste0("dysbiosis_index_", amplicon, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  log_msg("  Dysbiosis index saved (", nrow(dysbiosis_df), " samples)")
  
  # Print summary for integration compartments
  for (comp in c("leaf", "AdvRoot")) {
    sub <- dysbiosis_df[dysbiosis_df$compartment == comp, ]
    if (nrow(sub) > 0) {
      flight_scores <- sub$dysbiosis_index[sub$flight == "Flight"]
      ground_scores <- sub$dysbiosis_index[sub$flight == "Ground"]
      log_msg("    ", comp, ": Flight dysbiosis=", round(mean(flight_scores, na.rm=TRUE), 2),
              " ± ", round(sd(flight_scores, na.rm=TRUE), 2),
              " | Ground=", round(mean(ground_scores, na.rm=TRUE), 2),
              " ± ", round(sd(ground_scores, na.rm=TRUE), 2))
    }
  }
  
  return(dysbiosis_df)
}

# --- 4. ALDEx2 differential abundance ---------------------------------------
run_aldex2 <- function(ps, amplicon, compartment, contrast_var, group_levels) {
  log_msg("  ALDEx2: ", amplicon, " / ", compartment, " / ", contrast_var)
  
  # Subset to compartment
  sd <- data.frame(sample_data(ps))
  comp_samples <- rownames(sd)[sd$compartment == compartment]
  if (length(comp_samples) < 4) return(NULL)
  
  ps_comp <- prune_samples(comp_samples, ps)
  sd_comp <- data.frame(sample_data(ps_comp))
  
  # Subset to the two groups being compared
  sd_comp[[contrast_var]] <- as.character(sd_comp[[contrast_var]])
  keep_samples <- rownames(sd_comp)[sd_comp[[contrast_var]] %in% group_levels]
  if (length(keep_samples) < 4) return(NULL)
  
  ps_comp <- prune_samples(keep_samples, ps_comp)
  sd_comp <- data.frame(sample_data(ps_comp))
  
  # Get count matrix (taxa as rows) - ensure numeric
  otu <- as(otu_table(ps_comp), "matrix")
  if (!taxa_are_rows(ps_comp)) otu <- t(otu)
  storage.mode(otu) <- "integer"  # Ensure integer counts
  
  # Filter low-prevalence taxa (present in >= 10% of samples)
  prev_cutoff <- ceiling(ncol(otu) * 0.10)
  keep_taxa <- rowSums(otu > 0) >= prev_cutoff
  otu_filt <- otu[keep_taxa, ]
  log_msg("    Taxa after prevalence filter: ", nrow(otu_filt), " / ", nrow(otu))
  
  if (nrow(otu_filt) < 2) return(NULL)
  
  # Group vector (ALDEx2 requires character vector, not factor)
  groups <- as.character(sd_comp[[contrast_var]])
  
  # Run ALDEx2
  set.seed(42)
  clr <- aldex.clr(otu_filt, groups, mc.samples = ALDEX_MC)
  
  # Welch's t-test + Wilcoxon rank-sum (for 2-group comparison)
  welch <- aldex.ttest(clr, paired.test = FALSE)
  
  # Effect size
  effect <- aldex.effect(clr)
  
  # Combine results
  results <- data.frame(
    asv_id = rownames(otu_filt),
    compartment = compartment,
    contrast = paste0(contrast_var, "_", group_levels[1], "_vs_", group_levels[2]),
    welch_p = welch$we.eBH,    # Welch's t BH-adjusted p-value
    wilcox_p = welch$wi.eBH,   # Wilcoxon BH-adjusted p-value
    effect_size = effect$effect,
    diff_btw = effect$diff.btw,
    diff_win = effect$diff.win,
    overlap = effect$overlap
  )
  
  # Significance: BH-adjusted p < 0.05 and |effect| > 1
  results$significant <- results$welch_p < 0.05 & abs(results$effect_size) > 1
  
  # Add taxonomy
  tax <- data.frame(tax_table(ps_comp))
  tax$asv_id <- rownames(tax)
  results <- merge(results, tax, by = "asv_id", all.x = TRUE)
  
  return(results)
}

# =============================================================================
# RUN ANALYSIS
# =============================================================================
log_msg("Starting community health analysis for VEG-05 microbiome")
log_msg("R version: ", R.version.string)

all_alpha <- list()
all_beta <- list()
all_dysbiosis <- list()
all_aldex <- list()

for (amplicon in c("16S", "ITS")) {
  log_msg("\n========== ", amplicon, " ==========")
  
  ps <- load_phyloseq(amplicon)
  
  # 1. Alpha diversity
  all_alpha[[amplicon]] <- compute_alpha(ps, amplicon)
  
  # 2. Beta diversity + PERMANOVA
  all_beta[[amplicon]] <- compute_beta(ps, amplicon)
  
  # 3. Dysbiosis index
  all_dysbiosis[[amplicon]] <- compute_dysbiosis(ps, amplicon)
  
  # 4. ALDEx2 differential abundance (for integration compartments)
  for (comp in c("leaf", "AdvRoot")) {
    # Flight effect
    aldex_flight <- run_aldex2(ps, amplicon, comp, "flight", c("Flight", "Ground"))
    if (!is.null(aldex_flight)) {
      outfile <- file.path(OUT_DIR, paste0("aldex2_results_", amplicon, "_", comp, "_flight.tsv"))
      write.table(aldex_flight, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
      n_sig <- sum(aldex_flight$significant)
      log_msg("    Significant taxa: ", n_sig, " / ", nrow(aldex_flight))
      all_aldex[[paste0(amplicon, "_", comp, "_flight")]] <- aldex_flight
    }
    
    # Light effect (only for samples with light info)
    aldex_light <- run_aldex2(ps, amplicon, comp, "light", c("Red", "Blue"))
    if (!is.null(aldex_light)) {
      outfile <- file.path(OUT_DIR, paste0("aldex2_results_", amplicon, "_", comp, "_light.tsv"))
      write.table(aldex_light, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
      n_sig <- sum(aldex_light$significant)
      log_msg("    Significant taxa: ", n_sig, " / ", nrow(aldex_light))
      all_aldex[[paste0(amplicon, "_", comp, "_light")]] <- aldex_light
    }
  }
}

# --- Save combined alpha diversity table ---
alpha_all <- do.call(rbind, all_alpha)
write.table(alpha_all, file.path(OUT_DIR, "alpha_diversity_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Copy PDFs
pdfs <- list.files(PDF_DIR, full.names = TRUE)
for (p in pdfs) file.copy(p, file.path(OUT_DIR, basename(p)), overwrite = TRUE)

log_msg("\n========================================")
log_msg("COMMUNITY HEALTH ANALYSIS COMPLETE")
log_msg("Outputs saved to: ", OUT_DIR)
log_msg("========================================")
