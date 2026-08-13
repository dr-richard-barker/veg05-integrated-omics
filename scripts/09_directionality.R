#!/usr/bin/env Rscript
# ==============================================================================
# Step 9: Directionality Inference
#
# 1. GO enrichment of WGCNA modules (Biological Process)
#    Uses custom tomato GO annotations from Ensembl Plants BioMart
#    (org.Sl.eg.db not available for Bioconductor 3.20)
# 2. Directionality classification: host-driven vs microbe-driven vs co-regulated
#    - Host-driven:   module correlates with flight/light but NOT dysbiosis
#    - Microbe-driven: module correlates with dysbiosis but NOT flight/light
#    - Co-regulated:  module correlates with both environmental and microbial signals
# 3. Module-trait correlation heatmaps with directionality labels
#
# Inputs:
#   results/rnaseq/wgcna_Leaf/module_assignments.tsv
#   results/rnaseq/wgcna_Leaf/module_trait_correlations.tsv
#   results/rnaseq/wgcna_Adv-Root/ (same files)
#   results/integration/networks/module_taxon_correlations.tsv
#   data/rnaseq/tomato_go_mapping.tsv
#
# Outputs:
#   results/rnaseq/go_enrichment/go_enrichment_Leaf.tsv
#   results/rnaseq/go_enrichment/go_enrichment_AdvRoot.tsv
#   results/integration/networks/module_directionality.tsv
#   results/integration/networks/module_directionality_Leaf.pdf
#   results/integration/networks/module_directionality_AdvRoot.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
})

# --- Paths ---
REPO_ROOT <- "/mnt/shared-workspace/veg05-integrated-omics"
OUT_DIR   <- file.path(REPO_ROOT, "results/integration/networks")
GO_DIR    <- file.path(REPO_ROOT, "results/rnaseq/go_enrichment")
PDF_DIR   <- "/workspace/directionality_pdfs"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GO_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

# ==============================================================================
# Part 1: GO Enrichment of WGCNA Modules
# ==============================================================================
log_msg("=== Part 1: GO Enrichment of WGCNA Modules ===")

# Load GO mapping (from Ensembl Plants BioMart)
go_mapping <- read.csv(file.path(REPO_ROOT, "data/rnaseq/tomato_go_mapping.tsv"),
                       sep = "\t", stringsAsFactors = FALSE)
term2gene <- unique(go_mapping[, c("go_id", "gene_base")])
colnames(term2gene) <- c("term", "gene")
term2name <- unique(go_mapping[, c("go_id", "go_name")])
colnames(term2name) <- c("term", "name")

run_go_enrichment <- function(tissue, wgcna_dir) {
  log_msg("  GO enrichment for ", tissue)

  module_assignments <- read.csv(file.path(wgcna_dir, "module_assignments.tsv"),
                                  sep = "\t", stringsAsFactors = FALSE)

  # Strip version suffix from gene IDs (e.g. Solyc02g090330.3 -> Solyc02g090330)
  if ("gene_id" %in% colnames(module_assignments)) {
    module_assignments$gene_base <- sapply(module_assignments$gene_id,
      function(x) sub("\\.[0-9]+$", "", as.character(x)))
  } else {
    module_assignments$gene_base <- sapply(module_assignments[, 1],
      function(x) sub("\\.[0-9]+$", "", as.character(x)))
  }

  # Resolve module column name (may be "module" or "module_color")
  if (!"module" %in% colnames(module_assignments)) {
    if ("module_color" %in% colnames(module_assignments)) {
      module_assignments$module <- module_assignments$module_color
    } else {
      mod_col <- grep("module|color|label", colnames(module_assignments),
                      ignore.case = TRUE, value = TRUE)
      if (length(mod_col) > 0) module_assignments$module <- module_assignments[[mod_col[1]]]
    }
  }

  universe <- unique(module_assignments$gene_base[module_assignments$module != "grey"])
  log_msg("    Universe size: ", length(universe))

  all_enrich <- list()
  modules <- unique(module_assignments$module[
    module_assignments$module != "grey" & !is.na(module_assignments$module)])

  for (mod in modules) {
    gene_set <- module_assignments$gene_base[
      module_assignments$module == mod & !is.na(module_assignments$module)]
    if (length(gene_set) < 10) next

    tryCatch({
      ego <- enricher(gene = gene_set,
                      universe = universe,
                      TERM2GENE = term2gene,
                      TERM2NAME = term2name,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 5,
                      maxGSSize = 500)

      if (nrow(as.data.frame(ego)) > 0) {
        res_df <- as.data.frame(ego)
        res_df$tissue <- tissue
        res_df$module <- mod
        res_df$n_genes_in_module <- length(gene_set)
        all_enrich[[length(all_enrich) + 1]] <- res_df
        n_sig <- sum(res_df$p.adjust < 0.05)
        log_msg("    Module ", mod, " (n=", length(gene_set), "): ",
                n_sig, " significant GO terms")
      }
    }, error = function(e) {
      log_msg("    Module ", mod, ": error - ", conditionMessage(e))
    })
  }

  if (length(all_enrich) > 0) {
    return(do.call(rbind, all_enrich))
  } else {
    return(data.frame())
  }
}

# Run for both tissues
for (tissue in c("Leaf", "Adv-Root")) {
  wgcna_dir <- file.path(REPO_ROOT, "results/rnaseq", paste0("wgcna_", tissue))
  enrich_res <- run_go_enrichment(tissue, wgcna_dir)

  if (nrow(enrich_res) > 0) {
    tissue_clean <- gsub("-", "", tissue)
    write.table(enrich_res, file.path(GO_DIR, paste0("go_enrichment_", tissue_clean, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    log_msg("  Saved GO enrichment for ", tissue, ": ", nrow(enrich_res), " terms")

    sig_go <- enrich_res[enrich_res$p.adjust < 0.05, ]
    if (nrow(sig_go) > 0) {
      log_msg("  Top GO terms (", tissue, "):")
      sig_go <- sig_go[order(sig_go$p.adjust), ]
      for (i in 1:min(15, nrow(sig_go))) {
        log_msg("    ", sig_go$module[i], " -> ", sig_go$Description[i],
                " (p.adj=", signif(sig_go$p.adjust[i], 3),
                ", ", sig_go$Count[i], " genes)")
      }
    }
  } else {
    log_msg("  No GO enrichment results for ", tissue)
  }
}

# ==============================================================================
# Part 2: Directionality Inference
# ==============================================================================
log_msg("\n=== Part 2: Directionality Inference ===")

infer_directionality <- function(tissue) {
  log_msg("  Inferring directionality for ", tissue)

  mtc <- read.csv(file.path(REPO_ROOT, "results/rnaseq",
                            paste0("wgcna_", tissue),
                            "module_trait_correlations.tsv"),
                  sep = "\t", stringsAsFactors = FALSE)

  net_file <- file.path(OUT_DIR, "module_taxon_correlations.tsv")
  if (!file.exists(net_file)) return(data.frame())
  net <- read.csv(net_file, sep = "\t", stringsAsFactors = FALSE)
  net_tissue <- net[net$tissue == tissue, ]

  results <- list()
  for (i in 1:nrow(mtc)) {
    mod <- mtc$module_color[i]

    env_sig   <- (mtc$padj_flight[i] < 0.05) || (mtc$padj_light[i] < 0.05)
    env_rho   <- max(abs(mtc$cor_flight[i]), abs(mtc$cor_light[i]))
    micro_sig <- (mtc$padj_dysbiosis_16S[i] < 0.05) ||
                 (mtc$padj_dysbiosis_ITS[i] < 0.05)
    micro_rho <- max(abs(mtc$cor_dysbiosis_16S[i]),
                     abs(mtc$cor_dysbiosis_ITS[i]), na.rm = TRUE)

    mod_net <- net_tissue[net_tissue$module == mod & !is.na(net_tissue$padj), ]
    n_sig_taxa   <- sum(mod_net$padj < 0.05)
    max_taxon_rho <- if (nrow(mod_net) > 0) max(abs(mod_net$rho)) else 0

    if (mod == "grey") {
      direction <- "unassigned"
    } else if (env_sig && !micro_sig) {
      direction <- "host_driven"
    } else if (!env_sig && micro_sig) {
      direction <- "microbe_driven"
    } else if (env_sig && micro_sig) {
      direction <- "co_regulated"
    } else {
      direction <- "unresponsive"
    }

    results[[length(results) + 1]] <- data.frame(
      tissue = tissue,
      module = mod,
      n_genes = mtc$n_genes[i],
      cor_flight = mtc$cor_flight[i],
      padj_flight = mtc$padj_flight[i],
      cor_light = mtc$cor_light[i],
      padj_light = mtc$padj_light[i],
      cor_dysbiosis_16S = mtc$cor_dysbiosis_16S[i],
      padj_dysbiosis_16S = mtc$padj_dysbiosis_16S[i],
      cor_dysbiosis_ITS = mtc$cor_dysbiosis_ITS[i],
      padj_dysbiosis_ITS = mtc$padj_dysbiosis_ITS[i],
      n_sig_taxa_correlations = n_sig_taxa,
      max_taxon_rho = max_taxon_rho,
      directionality = direction
    )
  }

  return(do.call(rbind, results))
}

all_directions <- list()
for (tissue in c("Leaf", "Adv-Root")) {
  dir_res <- infer_directionality(tissue)
  all_directions[[tissue]] <- dir_res

  log_msg("\n  Directionality summary for ", tissue, ":")
  dir_table <- table(dir_res$directionality)
  for (d in names(dir_table)) {
    log_msg("    ", d, ": ", dir_table[d])
  }

  interesting <- dir_res[dir_res$directionality %in% c("microbe_driven", "co_regulated"), ]
  if (nrow(interesting) > 0) {
    log_msg("  Key modules:")
    for (i in 1:nrow(interesting)) {
      log_msg("    ", interesting$module[i], " (", interesting$directionality[i],
              ", n=", interesting$n_genes[i], ")")
    }
  }
}

direction_df <- do.call(rbind, all_directions)
write.table(direction_df, file.path(OUT_DIR, "module_directionality.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# Part 3: Module-Directionality Heatmaps
# ==============================================================================
log_msg("\n=== Part 3: Module-Directionality Heatmaps ===")

for (tissue in c("Leaf", "Adv-Root")) {
  tissue_clean <- gsub("-", "", tissue)
  mtc <- read.csv(file.path(REPO_ROOT, "results/rnaseq",
                            paste0("wgcna_", tissue),
                            "module_trait_correlations.tsv"),
                  sep = "\t", stringsAsFactors = FALSE)
  dir_res <- direction_df[direction_df$tissue == tissue, ]

  mtc$directionality <- dir_res$directionality[match(mtc$module_color, dir_res$module)]

  traits  <- c("flight", "light", "dysbiosis_16S", "dysbiosis_ITS")
  cor_cols  <- paste0("cor_", traits)
  padj_cols <- paste0("padj_", traits)

  cor_mat <- as.matrix(mtc[, cor_cols])
  rownames(cor_mat) <- mtc$module_color
  colnames(cor_mat) <- c("Flight", "Light", "16S Dysbiosis", "ITS Dysbiosis")

  padj_mat <- as.matrix(mtc[, padj_cols])
  star_mat <- matrix("", nrow(padj_mat), ncol(padj_mat))
  star_mat[padj_mat < 0.05]  <- "*"
  star_mat[padj_mat < 0.01]  <- "**"
  star_mat[padj_mat < 0.001] <- "***"

  pdf(file.path(PDF_DIR, paste0("module_directionality_", tissue_clean, ".pdf")),
      width = 10, height = 8)
  breaks <- seq(-1, 1, length.out = 101)
  col_palette <- colorRampPalette(c("blue", "white", "red"))(100)

  old_par <- par(mar = c(8, 6, 4, 10))
  image(1:ncol(cor_mat), 1:nrow(cor_mat), t(cor_mat),
        col = col_palette, breaks = breaks,
        xlab = "", ylab = "", axes = FALSE,
        main = paste0(tissue, " Module-Trait Correlations + Directionality"))
  axis(1, at = 1:ncol(cor_mat), labels = colnames(cor_mat), las = 2, cex.axis = 0.8)
  axis(2, at = 1:nrow(cor_mat), labels = rownames(cor_mat), las = 2, cex.axis = 0.7)
  for (i in 1:nrow(cor_mat)) {
    for (j in 1:ncol(cor_mat)) {
      text(j, i, star_mat[i, j], cex = 0.8)
    }
  }
  mtext(mtc$directionality, side = 4, at = 1:nrow(cor_mat), las = 2, cex = 0.6)
  par(old_par)
  dev.off()

  log_msg("  Saved directionality heatmap for ", tissue)
}

# Copy PDFs to results
for (p in list.files(PDF_DIR, full.names = TRUE)) {
  system(paste0("cp '", p, "' '", file.path(OUT_DIR, basename(p)), "'"))
}

# ==============================================================================
# Summary
# ==============================================================================
log_msg("\n========================================")
log_msg("DIRECTIONALITY INFERENCE COMPLETE")
log_msg("Outputs:")
log_msg("  GO enrichment: ", GO_DIR)
log_msg("  Module directionality: ", file.path(OUT_DIR, "module_directionality.tsv"))
log_msg("  Directionality heatmaps: ", OUT_DIR)
log_msg("========================================")
