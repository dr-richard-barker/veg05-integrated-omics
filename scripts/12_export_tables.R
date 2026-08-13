#!/usr/bin/env Rscript
# ==============================================================================
# Step 12: Export Supplementary Tables
#
#   Table S1:  Sample metadata (RNA-seq)
#   Table S2:  DEG summary (all contrasts, both tissues)
#   Table S3:  Alpha diversity (all samples, all metrics)
#   Table S4:  PERMANOVA results (16S + ITS)
#   Table S5:  Module directionality classification
#   Table S6:  FAPROTAX functional report (16S)
#   Table S7:  FUNGuild guild assignments (ITS)
#   Table S8:  MOFA+ factor-trait correlations
#   Table S9:  Module-taxon correlations (significant only)
#   Table S10: GO enrichment summary (significant terms per module)
#
# Inputs:
#   data/metadata/sample_metadata_rnaseq.csv
#   results/rnaseq/degs_summary.tsv
#   results/microbiome/community_health/alpha_diversity_all.tsv
#   results/microbiome/community_health/beta_diversity_permanova_16S.tsv
#   results/microbiome/community_health/beta_diversity_permanova_ITS.tsv
#   results/integration/networks/module_directionality.tsv
#   results/microbiome/functional/faprotax_report_16S.tsv
#   results/microbiome/functional/funguild_guilds_ITS_manual.tsv
#   results/integration/factor_trait_correlations.tsv
#   results/integration/networks/module_taxon_correlations.tsv
#   results/rnaseq/go_enrichment/go_enrichment_Leaf.tsv
#   results/rnaseq/go_enrichment/go_enrichment_AdvRoot.tsv
#
# Outputs:
#   results/supplementary_tables/Table_S1_sample_metadata_rnaseq.csv
#   results/supplementary_tables/Table_S2_deg_summary.csv
#   results/supplementary_tables/Table_S3_alpha_diversity.csv
#   results/supplementary_tables/Table_S4_permanova.csv
#   results/supplementary_tables/Table_S5_module_directionality.csv
#   results/supplementary_tables/Table_S6_faprotax.csv
#   results/supplementary_tables/Table_S7_funguild.csv
#   results/supplementary_tables/Table_S8_mofa_correlations.csv
#   results/supplementary_tables/Table_S9_module_taxon_correlations.csv
#   results/supplementary_tables/Table_S10_go_enrichment.csv
# ==============================================================================

# --- Paths ---
REPO_ROOT <- "/mnt/shared-workspace/veg05-integrated-omics"
TBL_DIR   <- file.path(REPO_ROOT, "results/supplementary_tables")

dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""), "\n", sep = "")
}

log_msg("Exporting supplementary tables...")

# --- Table S1: Sample metadata (RNA-seq) ---
log_msg("  Table S1: Sample metadata")
meta_rna <- read.csv(file.path(REPO_ROOT, "data/metadata/sample_metadata_rnaseq.csv"),
                     stringsAsFactors = FALSE)
write.csv(meta_rna, file.path(TBL_DIR, "Table_S1_sample_metadata_rnaseq.csv"), row.names = FALSE)

# --- Table S2: DEG summary ---
log_msg("  Table S2: DEG summary")
deg_summary <- read.csv(file.path(REPO_ROOT, "results/rnaseq/degs_summary.tsv"),
                        sep = "\t", stringsAsFactors = FALSE)
write.csv(deg_summary, file.path(TBL_DIR, "Table_S2_deg_summary.csv"), row.names = FALSE)

# --- Table S3: Alpha diversity ---
log_msg("  Table S3: Alpha diversity")
alpha_all <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/alpha_diversity_all.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
write.csv(alpha_all, file.path(TBL_DIR, "Table_S3_alpha_diversity.csv"), row.names = FALSE)

# --- Table S4: PERMANOVA ---
log_msg("  Table S4: PERMANOVA")
permanova_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/beta_diversity_permanova_16S.tsv"),
                           sep = "\t", stringsAsFactors = FALSE)
permanova_its <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/beta_diversity_permanova_ITS.tsv"),
                           sep = "\t", stringsAsFactors = FALSE)
permanova_all <- rbind(
  cbind(amplicon = "16S", permanova_16s),
  cbind(amplicon = "ITS", permanova_its)
)
write.csv(permanova_all, file.path(TBL_DIR, "Table_S4_permanova.csv"), row.names = FALSE)

# --- Table S5: Module directionality ---
log_msg("  Table S5: Module directionality")
directionality <- read.csv(file.path(REPO_ROOT, "results/integration/networks/module_directionality.tsv"),
                           sep = "\t", stringsAsFactors = FALSE)
write.csv(directionality, file.path(TBL_DIR, "Table_S5_module_directionality.csv"), row.names = FALSE)

# --- Table S6: FAPROTAX report ---
log_msg("  Table S6: FAPROTAX")
faprotax_report <- read.csv(file.path(REPO_ROOT, "results/microbiome/functional/faprotax_report_16S.tsv"),
                            sep = "\t", stringsAsFactors = FALSE)
# Filter out parsing artifact (taxon misclassified as function)
faprotax_report <- faprotax_report[!grepl("\\*", faprotax_report[, 1]), ]
write.csv(faprotax_report, file.path(TBL_DIR, "Table_S6_faprotax.csv"), row.names = FALSE)

# --- Table S7: FUNGuild ---
log_msg("  Table S7: FUNGuild")
funguild <- read.csv(file.path(REPO_ROOT, "results/microbiome/functional/funguild_guilds_ITS_manual.tsv"),
                     sep = "\t", stringsAsFactors = FALSE)
write.csv(funguild, file.path(TBL_DIR, "Table_S7_funguild.csv"), row.names = FALSE)

# --- Table S8: MOFA factor-trait correlations ---
log_msg("  Table S8: MOFA correlations")
mofa_cor <- read.csv(file.path(REPO_ROOT, "results/integration/factor_trait_correlations.tsv"),
                     sep = "\t", stringsAsFactors = FALSE)
write.csv(mofa_cor, file.path(TBL_DIR, "Table_S8_mofa_correlations.csv"), row.names = FALSE)

# --- Table S9: Module-taxon correlations (significant only) ---
log_msg("  Table S9: Module-taxon correlations")
net_df <- read.csv(file.path(REPO_ROOT, "results/integration/networks/module_taxon_correlations.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
sig_net_out <- net_df[!is.na(net_df$padj) & net_df$padj < 0.05, ]
write.csv(sig_net_out, file.path(TBL_DIR, "Table_S9_module_taxon_correlations.csv"), row.names = FALSE)

# --- Table S10: GO enrichment summary ---
log_msg("  Table S10: GO enrichment")
go_leaf <- read.csv(file.path(REPO_ROOT, "results/rnaseq/go_enrichment/go_enrichment_Leaf.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)
go_root <- read.csv(file.path(REPO_ROOT, "results/rnaseq/go_enrichment/go_enrichment_AdvRoot.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)
go_all <- rbind(go_leaf, go_root)
go_sig <- go_all[go_all$p.adjust < 0.05, ]
write.csv(go_sig, file.path(TBL_DIR, "Table_S10_go_enrichment.csv"), row.names = FALSE)

log_msg("\n========================================")
log_msg("SUPPLEMENTARY TABLES COMPLETE")
log_msg("Saved 10 tables to: ", TBL_DIR)
log_msg("========================================")
