#!/usr/bin/env Rscript
# ==============================================================================
# Step 10: Main Manuscript Figures (npj Microgravity)
#
#   Fig 2: Microbiome community composition (alpha diversity + dysbiosis)
#   Fig 3: Differential expression summary (DEG bar plot + volcano plots)
#   Fig 4: WGCNA module-trait correlations + directionality
#   Fig 5: MOFA+ integration (variance explained + factor-trait correlations)
#   Fig 6: Module-taxon bipartite network (significant correlations)
#   Fig 7: FAPROTAX functional summary (16S)
#
# All figures saved as SVG + PNG + PDF.
#
# Inputs:
#   results/microbiome/community_health/alpha_diversity_16S.tsv
#   results/microbiome/community_health/dysbiosis_index_16S.tsv
#   results/microbiome/community_health/dysbiosis_index_ITS.tsv
#   results/rnaseq/degs_summary.tsv
#   results/rnaseq/degs_Leaf_flight_vs_ground.tsv
#   results/rnaseq/degs_Adv-Root_flight_vs_ground.tsv
#   results/integration/networks/module_directionality.tsv
#   results/integration/factor_variance_explained.tsv
#   results/integration/factor_trait_correlations.tsv
#   results/integration/networks/module_taxon_correlations.tsv
#   results/microbiome/functional/faprotax_report_16S.tsv
#
# Outputs:
#   results/figures/fig2a_alpha_diversity_16S.{svg,png}
#   results/figures/fig2b_dysbiosis_16S.{svg,png}
#   results/figures/fig2b_dysbiosis_ITS.{svg,png}
#   results/figures/fig3_deg_summary.{svg,png}
#   results/figures/fig3b_volcano_leaf_flight.{svg,png}
#   results/figures/fig3c_volcano_advroot_flight.{svg,png}
#   results/figures/fig4_module_traits_Leaf.{svg,png}
#   results/figures/fig4_module_traits_AdvRoot.{svg,png}
#   results/figures/fig5a_mofa_variance.{svg,png}
#   results/figures/fig5b_mofa_correlations.{svg,png}
#   results/figures/fig6_module_taxon_network.{svg,png}
#   results/figures/fig7_faprotax.{svg,png}
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- Paths ---
REPO_ROOT <- "/mnt/shared-workspace/veg05-integrated-omics"
FIG_DIR   <- file.path(REPO_ROOT, "results/figures")
PDF_DIR   <- "/workspace/figure_pdfs"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Settings ---
theme_set(theme_bw(base_size = 11))

# Color palette (Phylo)
FLIGHT_COLORS <- c("Flight" = "#E9ED4C", "Ground" = "#75A025")
LIGHT_COLORS  <- c("Red" = "#FF9400", "Blue" = "#0279EE")

log_msg <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""), "\n", sep = "")
}

# Helper: save figure in all three formats
save_fig <- function(plot, name, w, h) {
  ggsave(file.path(PDF_DIR, paste0(name, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".svg")),  plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")),  plot, width = w, height = h, dpi = 300)
}

# ==============================================================================
# Fig 2: Microbiome Community Composition
# ==============================================================================
log_msg("Generating Fig 2: Microbiome community composition")

alpha_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/alpha_diversity_16S.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
dys_16s   <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/dysbiosis_index_16S.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
dys_its   <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/dysbiosis_index_ITS.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)

# Fig 2A: Alpha diversity (Observed) for 16S leaf and AdvRoot
alpha_16s$group <- paste(alpha_16s$flight, alpha_16s$light, sep = "_")
alpha_16s_main  <- alpha_16s[alpha_16s$compartment %in% c("leaf", "AdvRoot"), ]

fig2a <- ggplot(alpha_16s_main, aes(x = group, y = Observed, fill = flight)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_wrap(~compartment, scales = "free_y") +
  scale_fill_manual(values = FLIGHT_COLORS) +
  labs(x = "", y = "Observed ASVs", title = "16S Alpha Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")
save_fig(fig2a, "fig2a_alpha_diversity_16S", 8, 5)

# Fig 2B: Dysbiosis index (16S)
dys_16s_main <- dys_16s[dys_16s$compartment %in% c("leaf", "AdvRoot"), ]
fig2b <- ggplot(dys_16s_main, aes(x = flight, y = dysbiosis_index, fill = flight)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_wrap(~compartment, scales = "free_y") +
  scale_fill_manual(values = FLIGHT_COLORS) +
  labs(x = "", y = "Dysbiosis Index", title = "16S Dysbiosis Index") +
  theme(legend.position = "bottom")
save_fig(fig2b, "fig2b_dysbiosis_16S", 6, 5)

# Fig 2B (ITS): Dysbiosis index (ITS)
dys_its_main <- dys_its[dys_its$compartment %in% c("leaf", "AdvRoot"), ]
fig2b_its <- ggplot(dys_its_main, aes(x = flight, y = dysbiosis_index, fill = flight)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_wrap(~compartment, scales = "free_y") +
  scale_fill_manual(values = FLIGHT_COLORS) +
  labs(x = "", y = "Dysbiosis Index", title = "ITS Dysbiosis Index") +
  theme(legend.position = "bottom")
save_fig(fig2b_its, "fig2b_dysbiosis_ITS", 6, 5)

log_msg("  Saved Fig 2A-B")

# ==============================================================================
# Fig 3: Differential Expression Summary
# ==============================================================================
log_msg("Generating Fig 3: DEG summary")

deg_summary <- read.csv(file.path(REPO_ROOT, "results/rnaseq/degs_summary.tsv"),
                        sep = "\t", stringsAsFactors = FALSE)

deg_long <- deg_summary %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "n_degs") %>%
  mutate(direction = gsub("n_", "", direction))

fig3 <- ggplot(deg_long, aes(x = contrast, y = n_degs, fill = direction)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~tissue, scales = "free_x") +
  scale_fill_manual(values = c("up" = "#E9ED4C", "down" = "#0279EE")) +
  labs(x = "", y = "Number of DEGs (padj<0.05, |lfc|>=1)",
       title = "Differentially Expressed Genes", fill = "Direction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")
save_fig(fig3, "fig3_deg_summary", 10, 6)

# Volcano: Leaf flight_vs_ground
degs_leaf <- read.csv(file.path(REPO_ROOT, "results/rnaseq/degs_Leaf_flight_vs_ground.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
degs_leaf$significance <- ifelse(degs_leaf$padj < 0.05 & abs(degs_leaf$log2FoldChange) >= 1,
                                  "Significant", "NS")
fig3b <- ggplot(degs_leaf, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Significant" = "#FF9400", "NS" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (Flight vs Ground)", y = "-log10(padj)",
       title = "Leaf: Flight vs Ground") +
  theme(legend.position = "bottom")
save_fig(fig3b, "fig3b_volcano_leaf_flight", 7, 6)

# Volcano: AdvRoot flight_vs_ground
degs_root <- read.csv(file.path(REPO_ROOT, "results/rnaseq/degs_Adv-Root_flight_vs_ground.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
degs_root$significance <- ifelse(degs_root$padj < 0.05 & abs(degs_root$log2FoldChange) >= 1,
                                  "Significant", "NS")
fig3c <- ggplot(degs_root, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Significant" = "#FF9400", "NS" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  labs(x = "log2 Fold Change (Flight vs Ground)", y = "-log10(padj)",
       title = "Adv-Root: Flight vs Ground") +
  theme(legend.position = "bottom")
save_fig(fig3c, "fig3c_volcano_advroot_flight", 7, 6)

log_msg("  Saved Fig 3 (DEG summary + volcano plots)")

# ==============================================================================
# Fig 4: WGCNA Module-Trait Correlations + Directionality
# ==============================================================================
log_msg("Generating Fig 4: WGCNA module-trait correlations")

directionality <- read.csv(file.path(REPO_ROOT, "results/integration/networks/module_directionality.tsv"),
                           sep = "\t", stringsAsFactors = FALSE)

for (tissue in c("Leaf", "Adv-Root")) {
  dir_sub <- directionality[directionality$tissue == tissue, ]
  dir_sub <- dir_sub[dir_sub$module != "grey", ]

  traits  <- c("flight", "light", "dysbiosis_16S", "dysbiosis_ITS")
  cor_cols  <- paste0("cor_", traits)
  padj_cols <- paste0("padj_", traits)

  cor_mat <- as.matrix(dir_sub[, cor_cols])
  rownames(cor_mat) <- dir_sub$module
  colnames(cor_mat) <- c("Flight", "Light", "16S\nDysbiosis", "ITS\nDysbiosis")

  padj_mat <- as.matrix(dir_sub[, padj_cols])
  sig_labels <- matrix("", nrow(padj_mat), ncol(padj_mat))
  sig_labels[padj_mat < 0.05]  <- "*"
  sig_labels[padj_mat < 0.01]  <- "**"
  sig_labels[padj_mat < 0.001] <- "***"

  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Module", "Trait", "Correlation")
  sig_df <- as.data.frame(as.table(sig_labels))
  colnames(sig_df) <- c("Module", "Trait", "Sig")
  cor_df <- merge(cor_df, sig_df, by = c("Module", "Trait"))

  dir_lookup <- setNames(dir_sub$directionality, dir_sub$module)
  cor_df$Directionality <- dir_lookup[cor_df$Module]

  tissue_clean <- gsub("-", "", tissue)

  fig4 <- ggplot(cor_df, aes(x = Trait, y = Module, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Sig), size = 3, vjust = 0.8) +
    scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(x = "", y = "WGCNA Module",
         title = paste0(tissue, ": Module-Trait Correlations"),
         fill = "Spearman\nrho") +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          panel.grid = element_blank())

  h <- max(5, nrow(cor_mat) * 0.4)
  save_fig(fig4, paste0("fig4_module_traits_", tissue_clean), 6, h)
}

log_msg("  Saved Fig 4 (module-trait heatmaps)")

# ==============================================================================
# Fig 5: MOFA+ Integration
# ==============================================================================
log_msg("Generating Fig 5: MOFA+ integration")

ve_df  <- read.csv(file.path(REPO_ROOT, "results/integration/factor_variance_explained.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
cor_df <- read.csv(file.path(REPO_ROOT, "results/integration/factor_trait_correlations.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)

# Fig 5A: Variance explained
ve_long <- ve_df %>%
  pivot_longer(cols = starts_with("r2_"), names_to = "View", values_to = "R2") %>%
  mutate(View = gsub("r2_", "", View),
         View = recode(View, "transcriptome" = "Transcriptome",
                       "microbiome_16S" = "16S", "microbiome_ITS" = "ITS"),
         View = factor(View, levels = c("Transcriptome", "16S", "ITS")))

fig5a <- ggplot(ve_long, aes(x = factor, y = R2, fill = View)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Transcriptome" = "#0279EE", "16S" = "#75A025", "ITS" = "#FD9BED")) +
  labs(x = "", y = "Variance Explained (%)", title = "MOFA+ Variance Explained") +
  theme(legend.position = "bottom")
save_fig(fig5a, "fig5a_mofa_variance", 8, 5)

# Fig 5B: Factor-trait correlations
cor_df$significant <- cor_df$padj < 0.05
fig5b <- ggplot(cor_df, aes(x = factor, y = trait, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(significant, "*", "")), size = 5) +
  scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = "", y = "", title = "MOFA+ Factor-Trait Correlations",
       fill = "Spearman\nrho") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank())
save_fig(fig5b, "fig5b_mofa_correlations", 7, 5)

log_msg("  Saved Fig 5 (MOFA+ integration)")

# ==============================================================================
# Fig 6: Module-Taxon Network (key correlations)
# ==============================================================================
log_msg("Generating Fig 6: Module-taxon network")

net_df <- read.csv(file.path(REPO_ROOT, "results/integration/networks/module_taxon_correlations.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
sig_net <- net_df[!is.na(net_df$padj) & net_df$padj < 0.05, ]

if (nrow(sig_net) > 0) {
  fig6 <- ggplot(sig_net, aes(x = module, y = Genus)) +
    geom_point(aes(size = abs(rho), color = rho), alpha = 0.7) +
    facet_grid(amplicon ~ tissue, scales = "free", space = "free") +
    scale_color_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                          midpoint = 0, limits = c(-1, 1)) +
    scale_size_continuous(range = c(2, 8), name = "|rho|") +
    labs(x = "WGCNA Module", y = "Taxon (Genus)",
         title = "Significant Module-Taxon Correlations (padj<0.05)",
         color = "Spearman\nrho") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          panel.grid.minor = element_blank())
  save_fig(fig6, "fig6_module_taxon_network", 10, 8)
  log_msg("  Saved Fig 6 (module-taxon network, ", nrow(sig_net), " edges)")
} else {
  log_msg("  No significant correlations for Fig 6")
}

# ==============================================================================
# Fig 7: FAPROTAX Functional Summary
# ==============================================================================
log_msg("Generating Fig 7: FAPROTAX functional summary")

faprotax_report <- read.csv(file.path(REPO_ROOT, "results/microbiome/functional/faprotax_report_16S.tsv"),
                            sep = "\t", stringsAsFactors = FALSE)

# Filter out the misclassified taxon-as-function (parsing artifact)
faprotax_report <- faprotax_report[!grepl("\\*", faprotax_report[, 1]), ]
faprotax_top <- faprotax_report[order(faprotax_report$total_reads, decreasing = TRUE), ][1:15, ]
faprotax_top[, 1] <- factor(faprotax_top[, 1], levels = rev(faprotax_top[, 1]))

fig7 <- ggplot(faprotax_top, aes(x = total_reads, y = .data[[colnames(faprotax_top)[1]]], fill = total_reads)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#75A025", high = "#E9ED4C") +
  labs(x = "Total Reads", y = "Function",
       title = "FAPROTAX Functional Prediction (16S)") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))
save_fig(fig7, "fig7_faprotax", 8, 6)

log_msg("  Saved Fig 7 (FAPROTAX)")

# ==============================================================================
# Copy PDFs to results
# ==============================================================================
for (p in list.files(PDF_DIR, full.names = TRUE)) {
  system(paste0("cp '", p, "' '", file.path(FIG_DIR, basename(p)), "'"))
}

log_msg("\n========================================")
log_msg("MAIN FIGURES COMPLETE")
log_msg("Figures saved to: ", FIG_DIR)
log_msg("========================================")
