#!/usr/bin/env Rscript
# =============================================================================
# 06_wgcna.R
# Weighted Gene Co-expression Network Analysis for VEG-05 RNA-seq
#
# For each tissue (Leaf, AdvRoot):
#   1. Select top N most variable genes from VST matrix
#   2. Choose soft threshold (power) for scale-free topology
#   3. Build co-expression network, detect modules
#   4. Correlate modules with traits (flight, light, dysbiosis index)
#   5. Identify hub genes per module
#   6. Export module assignments, ME-trait correlations, hub genes
#
# Input:  results/rnaseq/vst_{tissue}.tsv
#         data/metadata/sample_metadata_rnaseq.csv
#         results/microbiome/community_health/dysbiosis_index_{amplicon}.tsv
#
# Output: results/rnaseq/wgcna_{tissue}/
#           module_assignments.tsv
#           module_trait_correlations.tsv
#           hub_genes.tsv
#           soft_threshold_{tissue}.pdf
#           module_colors_{tissue}.pdf
#           wgcna_{tissue}.rds
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
})

# Allow multi-threading
allowWGCNAThreads()

# --- Configuration ----------------------------------------------------------
REPO_ROOT  <- "/mnt/shared-workspace/veg05-integrated-omics"
VST_DIR    <- file.path(REPO_ROOT, "results/rnaseq")
META_FILE  <- file.path(REPO_ROOT, "data/metadata/sample_metadata_rnaseq.csv")
DYS_DIR    <- file.path(REPO_ROOT, "results/microbiome/community_health")
OUT_BASE   <- file.path(REPO_ROOT, "results/rnaseq")
PDF_DIR    <- "/workspace/wgcna_pdfs"
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters
TOP_N_GENES     <- 5000
MIN_MODULE_SIZE <- 30
MERGE_CUT_H     <- 0.25
SOFT_POWER_R2   <- 0.8
TISSUES         <- c("Leaf", "Adv-Root")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# --- Load dysbiosis index for module-trait correlation ----------------------
load_dysbiosis <- function() {
  dys_list <- list()
  for (amp in c("16S", "ITS")) {
    f <- file.path(DYS_DIR, paste0("dysbiosis_index_", amp, ".tsv"))
    if (file.exists(f)) {
      dys <- read.csv(f, sep = "\t", stringsAsFactors = FALSE)
      dys_list[[amp]] <- dys
    }
  }
  return(dys_list)
}

# --- Match dysbiosis to RNA-seq samples -------------------------------------
match_dysbiosis_to_rnaseq <- function(rnaseq_samples, dys_data) {
  # RNA-seq sample IDs look like: VEG-05-Flt-SN05-Leaf-Red_L008
  # Microbiome sample IDs look like: VEG-05F-SN05-leaf-tom-red-rich_S107_L001
  # Need to match by plant + tissue + flight + light
  
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
  
  dys_matched <- data.frame(
    sample_id = rnaseq_samples,
    dysbiosis_16S = NA_real_,
    dysbiosis_ITS = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(rnaseq_samples)) {
    sid <- rnaseq_samples[i]
    row <- meta[meta$sample_id == sid, ]
    if (nrow(row) == 0) next
    
    flight <- row$flight  # "Flight" or "Ground"
    plant <- row$plant    # "SN05"
    light <- row$light    # "Red" or "Blue"
    tissue <- row$tissue  # "Leaf" or "Adv-Root"
    
    # Map tissue names
    comp <- ifelse(tissue == "Leaf", "leaf", "AdvRoot")
    
    for (amp in c("16S", "ITS")) {
      if (is.null(dys_data[[amp]])) next
      dys <- dys_data[[amp]]
      # Match: same compartment, flight, and plant
      # Normalize plant ID: RNA-seq uses "SN01", microbiome uses "SN-01" for Ground
      plant_norm <- gsub("-", "", plant)  # "SN01" from either "SN01" or "SN-01"
      matches <- dys[dys$compartment == comp & 
                      dys$flight == flight &
                      grepl(plant_norm, gsub("-", "", dys$sample_id)), ]
      if (nrow(matches) > 0) {
        dys_matched$dysbiosis_16S[i] <- if (amp == "16S") matches$dysbiosis_index[1] else dys_matched$dysbiosis_16S[i]
        dys_matched$dysbiosis_ITS[i] <- if (amp == "ITS") matches$dysbiosis_index[1] else dys_matched$dysbiosis_ITS[i]
      }
    }
  }
  
  return(dys_matched)
}

# --- WGCNA per tissue -------------------------------------------------------
run_wgcna <- function(tissue) {
  log_msg("========== WGCNA: ", tissue, " ==========")
  
  out_dir <- file.path(OUT_BASE, paste0("wgcna_", tissue))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load VST matrix
  vst_file <- file.path(VST_DIR, paste0("vst_", tissue, ".tsv"))
  vst <- read.csv(vst_file, sep = "\t", row.names = 1, check.names = FALSE)
  vst <- as.matrix(vst)
  log_msg("VST matrix: ", nrow(vst), " genes x ", ncol(vst), " samples")
  
  # Load metadata
  meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_id
  meta <- meta[colnames(vst), ]
  
  # Select top N most variable genes (by MAD)
  gene_mad <- apply(vst, 1, mad)
  top_genes <- names(sort(gene_mad, decreasing = TRUE))[1:min(TOP_N_GENES, nrow(vst))]
  datExpr <- t(vst[top_genes, ])
  log_msg("Selected top ", length(top_genes), " variable genes")
  
  # Check for genes/samples with too many missing values
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      datExpr <- datExpr[, gsg$goodGenes]
      log_msg("Removed ", sum(!gsg$goodGenes), " bad genes")
    }
    if (sum(!gsg$goodSamples) > 0) {
      datExpr <- datExpr[gsg$goodSamples, ]
      log_msg("Removed ", sum(!gsg$goodSamples), " bad samples")
    }
  }
  
  # --- Soft threshold selection ---
  log_msg("Selecting soft threshold...")
  powers <- c(1:10, seq(12, 20, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                           networkType = "signed")
  
  # Choose power: first power where R² >= threshold
  power_table <- sft$fitIndices
  power_table$R2 <- -sign(power_table$slope) * power_table$SFT.R.sq
  chosen_power <- power_table$Power[power_table$R2 >= SOFT_POWER_R2][1]
  if (is.na(chosen_power)) {
    chosen_power <- power_table$Power[which.max(power_table$R2)]
  }
  log_msg("Chosen soft power: ", chosen_power, " (R²=", 
          round(power_table$R2[power_table$Power == chosen_power], 3), ")")
  
  # Plot soft threshold selection
  pdf(file.path(PDF_DIR, paste0("soft_threshold_", tissue, ".pdf")), width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(power_table$Power, power_table$R2,
       xlab = "Soft Power", ylab = "Scale-free R²",
       main = paste(tissue, "- Soft Threshold"))
  abline(h = SOFT_POWER_R2, col = "red")
  abline(v = chosen_power, col = "blue", lty = 2)
  plot(power_table$Power, power_table$mean.k.,
       xlab = "Soft Power", ylab = "Mean Connectivity",
       main = paste(tissue, "- Mean Connectivity"))
  dev.off()
  
  # --- Build network ---
  log_msg("Building co-expression network (power=", chosen_power, ")...")
  net <- blockwiseModules(
    datExpr,
    power = chosen_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = MIN_MODULE_SIZE,
    mergeCutHeight = MERGE_CUT_H,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    verbose = 3
  )
  
  moduleColors <- labels2colors(net$colors)
  n_modules <- length(unique(net$colors))
  n_grey <- sum(net$colors == 0)
  log_msg("Modules detected: ", n_modules, " (grey/unassigned: ", n_grey, " genes)")
  
  # Plot module dendrogram
  pdf(file.path(PDF_DIR, paste0("module_colors_", tissue, ".pdf")), width = 12, height = 6)
  plotDendroAndColors(net$dendrograms[[1]],
                      moduleColors[net$blockGenes[[1]]],
                      "Module Colors",
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05,
                      main = paste(tissue, "- Gene Dendrogram and Module Colors"))
  dev.off()
  
  # --- Module eigengenes ---
  MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs)
  
  # --- Module-trait correlations ---
  # Traits: flight (binary), light (binary), dysbiosis indices
  traits <- data.frame(
    flight = as.numeric(factor(meta$flight, levels = c("Ground", "Flight"))) - 1,
    light = as.numeric(factor(meta$light, levels = c("Red", "Blue"))) - 1,
    row.names = colnames(vst)
  )
  
  # Add dysbiosis indices
  dys_data <- load_dysbiosis()
  dys_matched <- match_dysbiosis_to_rnaseq(colnames(vst), dys_data)
  traits$dysbiosis_16S <- dys_matched$dysbiosis_16S
  traits$dysbiosis_ITS <- dys_matched$dysbiosis_ITS
  
  # Match MEs to traits (same sample order)
  MEs <- MEs[rownames(traits), ]
  
  # Correlate MEs with traits
  moduleTraitCor <- cor(MEs, traits, use = "pairwise.complete.obs")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  # Format results
  mt_results <- data.frame(
    tissue = tissue,
    module = colnames(MEs),
    module_color = gsub("ME", "", colnames(MEs)),
    n_genes = as.numeric(table(moduleColors)[gsub("ME", "", colnames(MEs))]),
    stringsAsFactors = FALSE
  )
  
  for (trait in colnames(traits)) {
    mt_results[[paste0("cor_", trait)]] <- moduleTraitCor[, trait]
    mt_results[[paste0("p_", trait)]] <- moduleTraitPvalue[, trait]
    mt_results[[paste0("padj_", trait)]] <- p.adjust(moduleTraitPvalue[, trait], method = "BH")
  }
  
  write.table(mt_results, file.path(out_dir, "module_trait_correlations.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  log_msg("Module-trait correlations saved")
  
  # Print significant correlations
  for (trait in c("flight", "light", "dysbiosis_16S", "dysbiosis_ITS")) {
    padj_col <- paste0("padj_", trait)
    cor_col <- paste0("cor_", trait)
    sig <- mt_results[!is.na(mt_results[[padj_col]]) & mt_results[[padj_col]] < 0.05, ]
    if (nrow(sig) > 0) {
      log_msg("  Significant module-", trait, " correlations (padj<0.05):")
      for (i in 1:nrow(sig)) {
        log_msg("    ", sig$module_color[i], " (n=", sig$n_genes[i], "): r=",
                round(sig[[cor_col]][i], 2), ", padj=", signif(sig[[padj_col]][i], 3))
      }
    }
  }
  
  # --- Module assignments ---
  module_assign <- data.frame(
    gene_id = colnames(datExpr),
    module_color = moduleColors,
    module_label = net$colors,
    tissue = tissue,
    stringsAsFactors = FALSE
  )
  # Add gene significance (correlation with each trait)
  for (trait in colnames(traits)) {
    GS <- cor(datExpr, traits[[trait]], use = "pairwise.complete.obs")
    module_assign[[paste0("GS_", trait)]] <- GS[, 1]
  }
  
  write.table(module_assign, file.path(out_dir, "module_assignments.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # --- Hub genes (top 10 per module by module membership) ---
  log_msg("Identifying hub genes...")
  hub_list <- list()
  for (color in unique(moduleColors)) {
    if (color == "grey") next
    module_genes <- colnames(datExpr)[moduleColors == color]
    if (length(module_genes) < 5) next
    
    # Module membership = correlation with module eigengene
    ME <- MEs[[paste0("ME", color)]]
    MM <- cor(datExpr[, module_genes, drop = FALSE], ME, use = "pairwise.complete.obs")
    
    # Hub genes: top 10 by |MM|
    top_idx <- order(abs(MM[, 1]), decreasing = TRUE)[1:min(10, length(module_genes))]
    
    for (idx in top_idx) {
      hub_list[[length(hub_list) + 1]] <- data.frame(
        tissue = tissue,
        module_color = color,
        gene_id = module_genes[idx],
        module_membership = MM[idx, 1],
        n_genes_in_module = length(module_genes)
      )
    }
  }
  
  hub_df <- do.call(rbind, hub_list)
  write.table(hub_df, file.path(out_dir, "hub_genes.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("Hub genes saved: ", nrow(hub_df), " genes across ", length(unique(hub_df$module_color)), " modules")
  
  # --- Save WGCNA object ---
  wgcna_obj <- list(
    datExpr = datExpr,
    moduleColors = moduleColors,
    MEs = MEs,
    traits = traits,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    softPower = chosen_power,
    net = net
  )
  saveRDS(wgcna_obj, file.path(out_dir, paste0("wgcna_", tissue, ".rds")))
  
  # Copy PDFs
  for (p in list.files(PDF_DIR, pattern = tissue, full.names = TRUE)) {
    file.copy(p, file.path(out_dir, basename(p)), overwrite = TRUE)
  }
  
  return(mt_results)
}

# =============================================================================
# RUN WGCNA
# =============================================================================
log_msg("Starting WGCNA for VEG-05 RNA-seq")
log_msg("R version: ", R.version.string)
log_msg("WGCNA version: ", as.character(packageVersion("WGCNA")))

all_mt <- do.call(rbind, lapply(TISSUES, run_wgcna))

log_msg("\n========================================")
log_msg("WGCNA COMPLETE")
log_msg("========================================")
print(all_mt[, c("tissue", "module_color", "n_genes", "cor_flight", "padj_flight", "cor_dysbiosis_16S", "padj_dysbiosis_16S")])
