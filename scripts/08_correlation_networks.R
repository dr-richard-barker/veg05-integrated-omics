#!/usr/bin/env Rscript
# ==============================================================================
# Step 8: Bipartite Correlation Networks
#
# Build bipartite module-taxon correlation networks by correlating WGCNA
# module eigengenes with genus-level taxa abundances (relative abundance)
# for matched RNA-seq + microbiome samples.
#
# Uses Spearman correlation with BH-FDR correction.
# Prevalence filter: taxa present in >20% of matched samples.
#
# Inputs:
#   results/rnaseq/wgcna_Leaf/wgcna_Leaf.rds
#   results/rnaseq/wgcna_Adv-Root/wgcna_Adv-Root.rds
#   results/microbiome/16S/phyloseq_16S_filtered.rds
#   results/microbiome/ITS/phyloseq_ITS.rds
#   data/metadata/sample_metadata_rnaseq.csv
#   data/metadata/sample_metadata_microbiome.csv
#   results/microbiome/16S/taxonomy_16S.tsv
#   results/microbiome/ITS/taxonomy_ITS.tsv
#
# Outputs:
#   results/integration/networks/module_taxon_correlations.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(WGCNA)
})

# Allow multi-threading for WGCNA
enableWGCNAThreads()

# --- Paths ---
REPO_ROOT <- "/mnt/shared-workspace/veg05-integrated-omics"
OUT_DIR   <- file.path(REPO_ROOT, "results/integration/networks")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

# ==============================================================================
# Load shared data
# ==============================================================================
log_msg("Loading metadata and taxonomy...")

meta_rna <- read.csv(file.path(REPO_ROOT, "data/metadata/sample_metadata_rnaseq.csv"),
                     stringsAsFactors = FALSE)
meta_mb  <- read.csv(file.path(REPO_ROOT, "data/metadata/sample_metadata_microbiome.csv"),
                     stringsAsFactors = FALSE)

normalize_key <- function(x) gsub("-", "", x)

tax_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/16S/taxonomy_16S.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)
tax_its <- read.csv(file.path(REPO_ROOT, "results/microbiome/ITS/taxonomy_ITS.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)

# Load phyloseq objects (filtered — chloroplast/mito removed for 16S)
ps_16s <- readRDS(file.path(REPO_ROOT, "results/microbiome/16S/phyloseq_16S_filtered.rds"))
ps_its <- readRDS(file.path(REPO_ROOT, "results/microbiome/ITS/phyloseq_ITS.rds"))

# ==============================================================================
# Build bipartite network for one tissue
# ==============================================================================
build_bipartite_network <- function(tissue, wgcna_dir, ps_16s, ps_its, meta_rna, meta_mb) {
  log_msg("  Building network for ", tissue)

  comp_name <- ifelse(tissue == "Leaf", "leaf", "AdvRoot")

  # Load WGCNA results
  wgcna <- readRDS(file.path(wgcna_dir, paste0("wgcna_", tissue, ".rds")))
  module_assignments <- read.csv(file.path(wgcna_dir, "module_assignments.tsv"),
                                  sep = "\t", stringsAsFactors = FALSE)

  MEs <- wgcna$MEs
  me_cols <- colnames(MEs)

  # Get RNA-seq samples for this tissue
  tissue_samples <- meta_rna[meta_rna$tissue == tissue, ]
  if (!is.null(rownames(MEs))) {
    me_sample_names <- rownames(MEs)
  } else {
    me_sample_names <- tissue_samples$sample_id[1:nrow(MEs)]
  }

  # Build match keys: Flight_SN05_Red (tissue-independent)
  me_match_keys <- sapply(me_sample_names, function(s) {
    parts <- strsplit(s, "_")[[1]]
    flight <- ifelse(grepl("Flt", s), "Flight", "Ground")
    plant <- gsub(".*SN", "SN", parts[1])
    plant <- gsub("-.*", "", plant)
    light <- ifelse(grepl("Red", s), "Red", "Blue")
    paste(flight, plant, light, sep = "_")
  })

  # --- 16S correlations ---
  sd_16s <- sample_data(ps_16s)
  keep_16s <- which(sd_16s$compartment == comp_name)
  if (length(keep_16s) == 0) {
    keep_16s <- which(grepl(comp_name, sd_16s$compartment, ignore.case = TRUE))
  }
  ps_16s_tissue <- prune_samples(sample_names(ps_16s)[keep_16s], ps_16s)

  ps_16s_ra <- transform_sample_counts(ps_16s_tissue, function(x) x / sum(x))
  otu_16s <- as(otu_table(ps_16s_ra), "matrix")
  if (taxa_are_rows(ps_16s_ra)) otu_16s <- t(otu_16s)

  mb_samples_16s <- sample_names(ps_16s_ra)
  mb_meta_16s <- as(sample_data(ps_16s_ra), "data.frame")

  mb_match_keys_16s <- sapply(mb_samples_16s, function(s) {
    flight <- mb_meta_16s[s, "flight"]
    plant <- normalize_key(mb_meta_16s[s, "plant"])
    light <- mb_meta_16s[s, "light"]
    if (is.na(light) || light == "") light <- "NA"
    paste(flight, plant, light, sep = "_")
  })

  common_16s <- intersect(me_match_keys, mb_match_keys_16s)
  log_msg("    16S common samples: ", length(common_16s))

  if (length(common_16s) < 5) {
    log_msg("    WARNING: Too few common 16S samples, skipping 16S correlations")
    cor_16s_df <- data.frame()
  } else {
    me_ordered <- MEs[match(common_16s, me_match_keys), , drop = FALSE]
    otu_16s_ordered <- otu_16s[match(common_16s, mb_match_keys_16s), , drop = FALSE]

    # Prevalence filter: present in >20% of matched samples
    prev_16s <- colMeans(otu_16s_ordered > 0)
    otu_16s_filt <- otu_16s_ordered[, prev_16s >= 0.2, drop = FALSE]
    log_msg("    16S ASVs after prevalence filter (>20%): ", ncol(otu_16s_filt))

    cor_16s_list <- list()
    for (me in colnames(me_ordered)) {
      me_color <- gsub("^ME", "", me)
      for (asv in colnames(otu_16s_filt)) {
        ct <- cor.test(me_ordered[, me], otu_16s_filt[, asv], method = "spearman")
        cor_16s_list[[length(cor_16s_list) + 1]] <- data.frame(
          tissue = tissue,
          amplicon = "16S",
          module = me_color,
          asv_id = asv,
          rho = ct$estimate,
          p_value = ct$p.value,
          n = sum(!is.na(me_ordered[, me]) & !is.na(otu_16s_filt[, asv]))
        )
      }
    }
    cor_16s_df <- do.call(rbind, cor_16s_list)
    cor_16s_df$padj <- p.adjust(cor_16s_df$p_value, method = "BH")
  }

  # --- ITS correlations ---
  sd_its <- sample_data(ps_its)
  keep_its <- which(sd_its$compartment == comp_name)
  if (length(keep_its) == 0) {
    keep_its <- which(grepl(comp_name, sd_its$compartment, ignore.case = TRUE))
  }
  ps_its_tissue <- prune_samples(sample_names(ps_its)[keep_its], ps_its)

  ps_its_ra <- transform_sample_counts(ps_its_tissue, function(x) x / sum(x))
  otu_its <- as(otu_table(ps_its_ra), "matrix")
  if (taxa_are_rows(ps_its_ra)) otu_its <- t(otu_its)

  mb_samples_its <- sample_names(ps_its_ra)
  mb_meta_its <- as(sample_data(ps_its_ra), "data.frame")

  mb_match_keys_its <- sapply(mb_samples_its, function(s) {
    flight <- mb_meta_its[s, "flight"]
    plant <- normalize_key(mb_meta_its[s, "plant"])
    light <- mb_meta_its[s, "light"]
    if (is.na(light) || light == "") light <- "NA"
    paste(flight, plant, light, sep = "_")
  })

  common_its <- intersect(me_match_keys, mb_match_keys_its)
  log_msg("    ITS common samples: ", length(common_its))

  if (length(common_its) < 5) {
    log_msg("    WARNING: Too few common ITS samples, skipping ITS correlations")
    cor_its_df <- data.frame()
  } else {
    me_ordered <- MEs[match(common_its, me_match_keys), , drop = FALSE]
    otu_its_ordered <- otu_its[match(common_its, mb_match_keys_its), , drop = FALSE]

    prev_its <- colMeans(otu_its_ordered > 0)
    otu_its_filt <- otu_its_ordered[, prev_its >= 0.2, drop = FALSE]
    log_msg("    ITS ASVs after prevalence filter (>20%): ", ncol(otu_its_filt))

    cor_its_list <- list()
    for (me in colnames(me_ordered)) {
      me_color <- gsub("^ME", "", me)
      for (asv in colnames(otu_its_filt)) {
        ct <- cor.test(me_ordered[, me], otu_its_filt[, asv], method = "spearman")
        cor_its_list[[length(cor_its_list) + 1]] <- data.frame(
          tissue = tissue,
          amplicon = "ITS",
          module = me_color,
          asv_id = asv,
          rho = ct$estimate,
          p_value = ct$p.value,
          n = sum(!is.na(me_ordered[, me]) & !is.na(otu_its_filt[, asv]))
        )
      }
    }
    cor_its_df <- do.call(rbind, cor_its_list)
    cor_its_df$padj <- p.adjust(cor_its_df$p_value, method = "BH")
  }

  # Combine and add genus annotation
  all_cor <- rbind(cor_16s_df, cor_its_df)

  if (nrow(cor_16s_df) > 0) {
    tax_lookup_16s <- tax_16s
    rownames(tax_lookup_16s) <- tax_lookup_16s$ASV_ID
    all_cor$Genus <- sapply(all_cor$asv_id, function(a) {
      if (a %in% rownames(tax_lookup_16s)) {
        g <- tax_lookup_16s[a, "Genus"]
        if (is.na(g)) return(tax_lookup_16s[a, "Family"])
        return(g)
      }
      return(NA)
    })
  }

  return(all_cor)
}

# ==============================================================================
# Build networks for both tissues
# ==============================================================================
log_msg("\n=== Building Bipartite Module-Taxon Networks ===")

all_networks <- list()
for (tissue in c("Leaf", "Adv-Root")) {
  wgcna_dir <- file.path(REPO_ROOT, "results/rnaseq", paste0("wgcna_", tissue))
  net <- build_bipartite_network(tissue, wgcna_dir, ps_16s, ps_its, meta_rna, meta_mb)
  all_networks[[tissue]] <- net
}

network_df <- do.call(rbind, all_networks)
write.table(network_df, file.path(OUT_DIR, "module_taxon_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# Summary
# ==============================================================================
sig_net <- network_df[!is.na(network_df$padj) & network_df$padj < 0.05, ]
log_msg("\n  Significant module-taxon correlations (padj<0.05): ", nrow(sig_net))
if (nrow(sig_net) > 0) {
  log_msg("  By tissue:")
  for (t in unique(sig_net$tissue)) {
    log_msg("    ", t, ": ", sum(sig_net$tissue == t))
  }
  log_msg("  By amplicon:")
  for (a in unique(sig_net$amplicon)) {
    log_msg("    ", a, ": ", sum(sig_net$amplicon == a))
  }
  sig_net <- sig_net[order(abs(sig_net$rho), decreasing = TRUE), ]
  log_msg("\n  Top 20 module-taxon correlations:")
  for (i in 1:min(20, nrow(sig_net))) {
    log_msg("    ", sig_net$tissue[i], " ", sig_net$amplicon[i],
            " ", sig_net$module[i], " ~ ", sig_net$asv_id[i],
            " (", sig_net$Genus[i], "): rho=", round(sig_net$rho[i], 2),
            ", padj=", signif(sig_net$padj[i], 3))
  }
}

if (nrow(sig_net) < 10) {
  log_msg("\n  Top unadjusted correlations (p<0.01):")
  top_unadj <- network_df[!is.na(network_df$p_value) & network_df$p_value < 0.01, ]
  top_unadj <- top_unadj[order(abs(top_unadj$rho), decreasing = TRUE), ]
  for (i in 1:min(20, nrow(top_unadj))) {
    log_msg("    ", top_unadj$tissue[i], " ", top_unadj$amplicon[i],
            " ", top_unadj$module[i], " ~ ", top_unadj$asv_id[i],
            " (", top_unadj$Genus[i], "): rho=", round(top_unadj$rho[i], 2),
            ", p=", signif(top_unadj$p_value[i], 3))
  }
}

log_msg("\n========================================")
log_msg("BIPARTITE CORRELATION NETWORKS COMPLETE")
log_msg("Output: ", file.path(OUT_DIR, "module_taxon_correlations.tsv"))
log_msg("========================================")
