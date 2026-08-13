#!/usr/bin/env Rscript
# =============================================================================
# 03_rnaseq_deseq2.R
# Differential expression analysis for VEG-05 RNA-seq (OSD-767)
#
# Uses RSEM unnormalized counts with DESeq2:
#   1. Load counts + metadata, filter low-count genes
#   2. Create DESeq2 dataset with design ~ flight * light
#   3. Run DESeq2 pipeline per tissue (Leaf, AdvRoot)
#   4. Extract 5 contrasts per tissue with LFC shrinkage (apeglm)
#   5. VST transformation for downstream visualization/integration
#   6. Export results tables, VST matrix, DESeq2 objects
#
# Input:  data/rnaseq/GLDS-709_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv
#         data/metadata/sample_metadata_rnaseq.csv
#
# Output: results/rnaseq/
#           dds_{tissue}.rds              — DESeq2 object
#           vst_{tissue}.tsv              — VST-transformed expression matrix
#           degs_{tissue}_{contrast}.tsv  — DE results per contrast
#           degs_summary.tsv              — Summary of DEG counts
#           qc_sample_distances_{tissue}.pdf — Sample distance heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

# --- Configuration ----------------------------------------------------------
REPO_ROOT  <- "/mnt/shared-workspace/veg05-integrated-omics"
COUNTS_FILE <- file.path(REPO_ROOT, "data/rnaseq/GLDS-709_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv")
META_FILE  <- file.path(REPO_ROOT, "data/metadata/sample_metadata_rnaseq.csv")
OUT_DIR    <- file.path(REPO_ROOT, "results/rnaseq")
PDF_DIR    <- "/workspace/rnaseq_pdfs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters from config
MIN_COUNTS   <- 10    # Minimum count in a sample
MIN_SAMPLES  <- 4     # Minimum samples passing count threshold
PADJ_THRESH  <- 0.05
LFC_THRESH   <- 1.0   # |log2FC| threshold
TISSUES      <- c("Leaf", "Adv-Root")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# --- Load data --------------------------------------------------------------
log_msg("Loading RSEM counts...")
counts <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
# Clean gene IDs: remove "gene-" prefix
rownames(counts) <- gsub("^gene-", "", rownames(counts))
# Round to integers (RSEM outputs fractional counts)
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"
log_msg("Counts matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta$flight <- factor(meta$flight, levels = c("Ground", "Flight"))
meta$light <- factor(meta$light, levels = c("Red", "Blue"))
meta$tissue <- factor(meta$tissue)
log_msg("Metadata: ", nrow(meta), " samples")

# Verify sample matching
stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), ]

# --- DESeq2 analysis per tissue ---------------------------------------------
run_deseq2 <- function(tissue_name) {
  log_msg("========== DESeq2: ", tissue_name, " ==========")
  
  # Subset to tissue
  samples <- rownames(meta)[meta$tissue == tissue_name]
  counts_sub <- counts[, samples, drop = FALSE]
  meta_sub <- meta[samples, , drop = FALSE]
  log_msg("Samples: ", length(samples))
  
  # Print group sizes
  log_msg("Group sizes:")
  print(table(meta_sub$flight, meta_sub$light))
  
  # Filter low-count genes: keep genes with >= MIN_COUNTS in >= MIN_SAMPLES
  keep <- rowSums(counts_sub >= MIN_COUNTS) >= MIN_SAMPLES
  counts_filt <- counts_sub[keep, , drop = FALSE]
  log_msg("Genes after filtering: ", nrow(counts_filt), " (removed ", nrow(counts_sub) - nrow(counts_filt), ")")
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filt,
    colData = meta_sub,
    design = ~ flight * light
  )
  
  # Run DESeq2
  log_msg("Running DESeq2...")
  dds <- DESeq(dds, parallel = FALSE)
  log_msg("DESeq2 complete")
  
  # --- Extract results for each contrast ---
  contrasts <- list(
    flight_vs_ground = c("flight", "Flight", "Ground"),
    red_vs_blue = c("light", "Red", "Blue"),
    flt_red_vs_gnd_red = c("flight", "Flight", "Ground"),  # subset to Red
    flt_blue_vs_gnd_blue = c("flight", "Flight", "Ground"),  # subset to Blue
    interaction = NULL  # interaction term
  )
  
  degs_all <- list()
  
  for (contrast_name in names(contrasts)) {
    log_msg("  Contrast: ", contrast_name)
    
    if (contrast_name == "interaction") {
      # Interaction: flightFlight:lightBlue
      res <- results(dds, name = "flightFlight.lightBlue", alpha = PADJ_THRESH)
    } else if (contrast_name %in% c("flt_red_vs_gnd_red", "flt_blue_vs_gnd_blue")) {
      # Subset to specific light condition
      light_val <- ifelse(contrast_name == "flt_red_vs_gnd_red", "Red", "Blue")
      dds_sub <- dds[, dds$light == light_val]
      dds_sub <- DESeqDataSetFromMatrix(
        countData = counts(dds_sub),
        colData = colData(dds_sub),
        design = ~ flight
      )
      dds_sub <- DESeq(dds_sub, parallel = FALSE)
      res <- results(dds_sub, contrast = c("flight", "Flight", "Ground"), alpha = PADJ_THRESH)
    } else {
      coef <- contrasts[[contrast_name]]
      res <- results(dds, contrast = coef, alpha = PADJ_THRESH)
    }
    
    # LFC shrinkage: use apeglm with coef when possible, ashr for contrast-based
    if (contrast_name == "interaction") {
      res_shr <- lfcShrink(dds, coef = "flightFlight.lightBlue", type = "apeglm")
    } else if (contrast_name == "flight_vs_ground") {
      res_shr <- lfcShrink(dds, coef = "flight_Flight_vs_Ground", type = "apeglm")
    } else if (contrast_name == "red_vs_blue") {
      # Red vs Blue = negative of Blue_vs_Red coef
      res_shr <- lfcShrink(dds, coef = "light_Blue_vs_Red", type = "apeglm")
      res_shr$log2FoldChange <- -res_shr$log2FoldChange
    } else if (contrast_name %in% c("flt_red_vs_gnd_red", "flt_blue_vs_gnd_blue")) {
      # For subset analyses, use ashr (apeglm needs the full model)
      light_val <- ifelse(contrast_name == "flt_red_vs_gnd_red", "Red", "Blue")
      dds_sub <- dds[, dds$light == light_val]
      dds_sub <- DESeqDataSetFromMatrix(
        countData = counts(dds_sub),
        colData = colData(dds_sub),
        design = ~ flight
      )
      dds_sub <- DESeq(dds_sub, parallel = FALSE)
      res_shr <- lfcShrink(dds_sub, contrast = c("flight", "Flight", "Ground"), type = "ashr")
    } else {
      res_shr <- lfcShrink(dds, contrast = contrasts[[contrast_name]], type = "ashr")
    }
    
    # Create results dataframe
    degs <- as.data.frame(res_shr)
    degs$gene_id <- rownames(degs)
    degs$contrast <- contrast_name
    degs$tissue <- tissue_name
    
    # Classify significance
    degs$significant <- degs$padj < PADJ_THRESH & abs(degs$log2FoldChange) >= LFC_THRESH
    degs$significant[is.na(degs$significant)] <- FALSE
    
    # Reorder columns
    # Note: shrunk results may not have 'stat' column
    available_cols <- intersect(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), colnames(degs))
    degs <- degs[, c("gene_id", "tissue", "contrast", available_cols, "significant")]
    
    # Sort by padj
    degs <- degs[order(degs$padj), ]
    
    # Save
    outfile <- file.path(OUT_DIR, paste0("degs_", tissue_name, "_", contrast_name, ".tsv"))
    write.table(degs, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
    
    n_sig <- sum(degs$significant)
    n_up <- sum(degs$significant & degs$log2FoldChange > 0)
    n_down <- sum(degs$significant & degs$log2FoldChange < 0)
    log_msg("    DEGs: ", n_sig, " (up=", n_up, ", down=", n_down, ")")
    
    degs_all[[contrast_name]] <- data.frame(
      tissue = tissue_name,
      contrast = contrast_name,
      n_total = nrow(degs),
      n_sig = n_sig,
      n_up = n_up,
      n_down = n_down
    )
  }
  
  # --- VST transformation ---
  log_msg("Computing VST transformation...")
  vsd <- vst(dds, blind = TRUE)
  vst_mat <- assay(vsd)
  colnames(vst_mat) <- colnames(dds)
  write.table(as.data.frame(vst_mat), file.path(OUT_DIR, paste0("vst_", tissue_name, ".tsv")),
              sep = "\t", quote = FALSE)
  log_msg("VST matrix saved: ", nrow(vst_mat), " genes x ", ncol(vst_mat), " samples")
  
  # --- Sample distance QC plot ---
  sample_dist <- dist(t(vst_mat))
  sample_dist_mat <- as.matrix(sample_dist)
  rownames(sample_dist_mat) <- paste(meta_sub$flight, meta_sub$light, sep = "_")
  colnames(sample_dist_mat) <- rownames(sample_dist_mat)
  
  pdf(file.path(PDF_DIR, paste0("qc_sample_distances_", tissue_name, ".pdf")), width = 10, height = 9)
  col_fun <- colorRamp2(c(0, max(sample_dist_mat)/2, max(sample_dist_mat)),
                         c("#4ECDC4", "#FAF9F3", "#FF6B6B"))
  ha <- HeatmapAnnotation(
    Flight = meta_sub$flight,
    Light = meta_sub$light,
    col = list(Flight = c(Flight = "#FF6B6B", Ground = "#4ECDC4"),
               Light = c(Red = "#E9ED4C", Blue = "#0279EE"))
  )
  draw(Heatmap(sample_dist_mat, name = "Distance",
               top_annotation = ha,
               col = col_fun,
               show_row_names = TRUE, show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 7),
               column_names_gp = gpar(fontsize = 7)))
  dev.off()
  
  # Save DESeq2 object
  saveRDS(dds, file.path(OUT_DIR, paste0("dds_", tissue_name, ".rds")))
  saveRDS(vsd, file.path(OUT_DIR, paste0("vsd_", tissue_name, ".rds")))
  
  return(do.call(rbind, degs_all))
}

# --- Run for both tissues ---------------------------------------------------
log_msg("Starting DESeq2 analysis for VEG-05 RNA-seq")
log_msg("R version: ", R.version.string)
log_msg("DESeq2 version: ", as.character(packageVersion("DESeq2")))

all_summary <- do.call(rbind, lapply(TISSUES, run_deseq2))

# --- Save summary ---
write.table(all_summary, file.path(OUT_DIR, "degs_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Copy PDFs to results
pdfs <- list.files(PDF_DIR, full.names = TRUE)
for (p in pdfs) file.copy(p, file.path(OUT_DIR, basename(p)), overwrite = TRUE)

log_msg("========================================")
log_msg("DESeq2 PIPELINE COMPLETE")
log_msg("Summary:")
print(all_summary)
log_msg("Outputs saved to: ", OUT_DIR)
log_msg("========================================")
