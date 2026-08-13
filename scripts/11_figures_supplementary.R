#!/usr/bin/env Rscript
# ==============================================================================
# Step 11: Supplementary Figures
#
#   Fig S1: DADA2 read retention through pipeline (16S + ITS)
#   Fig S2: Alpha diversity all metrics x all compartments
#   Fig S3: Bray-Curtis PCoA for all compartments (16S + ITS)
#   Fig S4: Stacked bar charts of top 15 genera per compartment (16S)
#   Fig S5: Full module-taxon correlation heatmap (all significant)
#   Fig S6: MOFA+ factor weight dot plots (top features per factor per view)
#   Fig S7: WGCNA soft power scale-free topology fit curves
#
# Uses only base R + ggplot2/dplyr/tidyr (no phyloseq/vegan/WGCNA required).
# Bray-Curtis dissimilarity and PCoA computed manually.
# Power curve recomputed from datExpr stored in WGCNA RDS.
#
# Inputs:
#   results/microbiome/16S/track_stats_16S.tsv
#   results/microbiome/ITS/track_stats_ITS.tsv
#   results/microbiome/community_health/alpha_diversity_all.tsv
#   results/microbiome/16S/asv_table_16S.tsv
#   results/microbiome/ITS/asv_table_ITS.tsv
#   results/microbiome/16S/taxonomy_16S.tsv
#   results/microbiome/ITS/taxonomy_ITS.tsv
#   data/metadata/sample_metadata_microbiome.csv
#   results/integration/networks/module_taxon_correlations.tsv
#   results/integration/factor_weights_top.tsv
#   results/rnaseq/wgcna_Leaf/wgcna_Leaf.rds
#   results/rnaseq/wgcna_Adv-Root/wgcna_Adv-Root.rds
#
# Outputs:
#   results/figures/figS1_dada2_read_retention.{svg,png}
#   results/figures/figS2_alpha_diversity_all.{svg,png}
#   results/figures/figS3_pcoa_bray_16S.{svg,png}
#   results/figures/figS3_pcoa_bray_ITS.{svg,png}
#   results/figures/figS4_stacked_bar_16S.{svg,png}
#   results/figures/figS5_module_taxon_heatmap.{svg,png}
#   results/figures/figS6_mofa_factor_weights.{svg,png}
#   results/figures/figS7_wgcna_power_curve.{svg,png}
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# --- Paths ---
REPO_ROOT <- "/mnt/shared-workspace/veg05-integrated-omics"
FIG_DIR   <- file.path(REPO_ROOT, "results/figures")
PDF_DIR   <- "/workspace/supp_figure_pdfs"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Settings ---
theme_set(theme_bw(base_size = 11))

FLIGHT_COLORS <- c("Flight" = "#E9ED4C", "Ground" = "#75A025")
LIGHT_COLORS  <- c("Red" = "#FF9400", "Blue" = "#0279EE")

log_msg <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""), "\n", sep = "")
}

save_fig <- function(plot, name, w, h) {
  ggsave(file.path(PDF_DIR, paste0(name, ".pdf")), plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".svg")),  plot, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")),  plot, width = w, height = h, dpi = 300)
}

# ==============================================================================
# Fig S1: DADA2 Read Retention
# ==============================================================================
log_msg("Generating Fig S1: DADA2 read retention")

track_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/16S/track_stats_16S.tsv"),
                      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
track_its <- read.csv(file.path(REPO_ROOT, "results/microbiome/ITS/track_stats_ITS.tsv"),
                      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Add amplicon label and combine
track_16s$amplicon <- "16S"
track_its$amplicon <- "ITS"
track_all <- rbind(track_16s, track_its)

# Reshape to long format
steps <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track_long <- track_all %>%
  pivot_longer(cols = all_of(steps), names_to = "step", values_to = "reads") %>%
  mutate(step = factor(step, levels = steps),
         amplicon = factor(amplicon, levels = c("16S", "ITS")))

# Summarize: median reads at each step
track_summary <- track_long %>%
  group_by(amplicon, step) %>%
  summarise(median_reads = median(reads, na.rm = TRUE),
            q25 = quantile(reads, 0.25, na.rm = TRUE),
            q75 = quantile(reads, 0.75, na.rm = TRUE),
            .groups = "drop")

figS1 <- ggplot(track_long, aes(x = step, y = reads, fill = amplicon)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c("16S" = "#75A025", "ITS" = "#FD9BED")) +
  scale_y_log10() +
  labs(x = "DADA2 Pipeline Step", y = "Reads (log scale)",
       title = "DADA2 Read Retention Through Pipeline",
       fill = "Amplicon") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")
save_fig(figS1, "figS1_dada2_read_retention", 8, 5)

log_msg("  Saved Fig S1")

# ==============================================================================
# Fig S2: Alpha Diversity All Metrics x All Compartments
# ==============================================================================
log_msg("Generating Fig S2: Alpha diversity all metrics")

alpha_all <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/alpha_diversity_all.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)

# Filter to plant-relevant compartments
plant_comps <- c("leaf", "AdvRoot", "root", "wick", "soil")
alpha_all <- alpha_all[alpha_all$compartment %in% plant_comps, ]
alpha_all$compartment <- factor(alpha_all$compartment, levels = plant_comps)
alpha_all$flight <- factor(alpha_all$flight, levels = c("Flight", "Ground"))

# Reshape to long for faceting
alpha_long <- alpha_all %>%
  select(sample_id, Observed, Shannon, Simpson, amplicon, compartment, flight, light) %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson")))

figS2 <- ggplot(alpha_long, aes(x = compartment, y = value, fill = flight)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_grid(metric ~ amplicon, scales = "free_y") +
  scale_fill_manual(values = FLIGHT_COLORS) +
  labs(x = "Compartment", y = "Diversity Value",
       title = "Alpha Diversity Across Compartments",
       fill = "Flight") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        legend.position = "bottom")
save_fig(figS2, "figS2_alpha_diversity_all", 10, 8)

log_msg("  Saved Fig S2")

# ==============================================================================
# Fig S3: Bray-Curtis PCoA for All Compartments
# ==============================================================================
log_msg("Generating Fig S3: Bray-Curtis PCoA")

# Helper: compute Bray-Curtis dissimilarity matrix (samples in rows, taxa in cols)
bray_curtis <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      num <- sum(abs(mat[i, ] - mat[j, ]))
      den <- sum(mat[i, ] + mat[j, ])
      if (den > 0) d[i, j] <- num / den else d[i, j] <- 1
      d[j, i] <- d[i, j]
    }
  }
  return(as.dist(d))
}

# Load microbiome metadata
meta_mb <- read.csv(file.path(REPO_ROOT, "data/metadata/sample_metadata_microbiome.csv"),
                    stringsAsFactors = FALSE)

run_pcoa <- function(amplicon_label, asv_file, out_name) {
  asv <- read.csv(asv_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  asv_ids <- asv[, 1]
  count_data <- t(asv[, -1])
  colnames(count_data) <- asv_ids
  # count_data: rows = samples, cols = ASVs

  # Get metadata for these samples
  sample_names <- rownames(count_data)
  meta_sub <- meta_mb[match(sample_names, meta_mb$sample_id), ]

  # Filter to plant-relevant compartments with enough samples
  plant_comps <- c("leaf", "AdvRoot", "root", "wick", "soil")
  keep <- meta_sub$compartment %in% plant_comps
  count_data <- count_data[keep, ]
  meta_sub <- meta_sub[keep, ]

  # Remove samples with zero total reads
  row_sums <- rowSums(count_data)
  count_data <- count_data[row_sums > 0, ]
  meta_sub <- meta_sub[row_sums > 0, ]

  if (nrow(count_data) < 5) return(NULL)

  # Relative abundance transform
  ra <- sweep(count_data, 1, rowSums(count_data), "/")

  # Compute Bray-Curtis and PCoA
  bc_dist <- bray_curtis(ra)
  pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
  var_explained <- pcoa_res$eig[1:2] / sum(pcoa_res$eig) * 100

  pcoa_df <- data.frame(
    PC1 = pcoa_res$points[, 1],
    PC2 = pcoa_res$points[, 2],
    compartment = factor(meta_sub$compartment, levels = plant_comps),
    flight = factor(meta_sub$flight, levels = c("Flight", "Ground"))
  )

  COMP_COLORS <- c("leaf" = "#75A025", "AdvRoot" = "#FD9BED", "root" = "#FF9400",
                   "wick" = "#0279EE", "soil" = "#E9ED4C")

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = compartment, shape = flight)) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_ellipse(aes(group = compartment), level = 0.95, linewidth = 0.5, show.legend = FALSE) +
    scale_color_manual(values = COMP_COLORS) +
    labs(x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
         title = paste0("Bray-Curtis PCoA (", amplicon_label, ")"),
         color = "Compartment", shape = "Flight") +
    theme(legend.position = "right")

  save_fig(p, out_name, 9, 6)
  log_msg("  Saved ", out_name, " (", nrow(pcoa_df), " samples)")
}

run_pcoa("16S", file.path(REPO_ROOT, "results/microbiome/16S/asv_table_16S.tsv"),
         "figS3_pcoa_bray_16S")
run_pcoa("ITS", file.path(REPO_ROOT, "results/microbiome/ITS/asv_table_ITS.tsv"),
         "figS3_pcoa_bray_ITS")

# ==============================================================================
# Fig S4: Stacked Bar Charts of Top 15 Genera per Compartment (16S)
# ==============================================================================
log_msg("Generating Fig S4: Stacked bar charts (16S genera)")

asv_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/16S/asv_table_16S.tsv"),
                    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
tax_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/16S/taxonomy_16S.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)

asv_ids <- asv_16s[, 1]
count_data <- asv_16s[, -1]
rownames(count_data) <- asv_ids

# Map ASV to genus
tax_lookup <- tax_16s
rownames(tax_lookup) <- tax_lookup$ASV_ID
genera <- sapply(asv_ids, function(a) {
  g <- tax_lookup[a, "Genus"]
  if (is.na(g)) {
    f <- tax_lookup[a, "Family"]
    if (is.na(f)) return("Unclassified")
    return(paste0("Unclassified ", f))
  }
  return(g)
})

# Get sample metadata
sample_names <- colnames(count_data)
meta_sub <- meta_mb[match(sample_names, meta_mb$sample_id), ]

# Filter to plant-relevant compartments
plant_comps <- c("leaf", "AdvRoot", "root", "wick", "soil")
keep <- meta_sub$compartment %in% plant_comps
count_data <- count_data[, keep]
meta_sub <- meta_sub[keep, ]

# Convert to matrix and compute relative abundance per sample
count_mat <- as.matrix(count_data)
# Remove samples with zero total reads
col_sums <- colSums(count_mat)
count_mat <- count_mat[, col_sums > 0, drop = FALSE]
meta_sub <- meta_sub[col_sums > 0, , drop = FALSE]
ra <- sweep(count_mat, 2, colSums(count_mat), "/")

# Aggregate by genus using rowsum (robust, no length mismatch issues)
ra_genus <- base::rowsum(ra, group = genera)
# ra_genus: rows = genera, cols = samples

# For each compartment, compute mean relative abundance
genus_names <- rownames(ra_genus)
comp_cols <- list()
for (comp in plant_comps) {
  comp_samples <- meta_sub$sample_id[meta_sub$compartment == comp]
  comp_samples <- intersect(comp_samples, colnames(ra_genus))
  if (length(comp_samples) > 0) {
    comp_cols[[comp]] <- as.numeric(rowMeans(ra_genus[, comp_samples, drop = FALSE]))
  } else {
    comp_cols[[comp]] <- rep(0, length(genus_names))
  }
}
# Build data.frame with explicit numeric columns
comp_df <- data.frame(
  genus = genus_names,
  leaf    = as.numeric(comp_cols[["leaf"]]),
  AdvRoot = as.numeric(comp_cols[["AdvRoot"]]),
  root    = as.numeric(comp_cols[["root"]]),
  wick    = as.numeric(comp_cols[["wick"]]),
  soil    = as.numeric(comp_cols[["soil"]]),
  stringsAsFactors = FALSE
)

# Get top 15 genera overall (by mean across all compartments)
comp_df$mean_all <- rowMeans(comp_df[, plant_comps, drop = FALSE])
top_genera <- comp_df[order(comp_df$mean_all, decreasing = TRUE), ][1:15, ]
other_abun <- comp_df[!(comp_df$genus %in% top_genera$genus), ]
other_row <- data.frame(
  genus = "Other",
  leaf    = sum(other_abun$leaf),
  AdvRoot = sum(other_abun$AdvRoot),
  root    = sum(other_abun$root),
  wick    = sum(other_abun$wick),
  soil    = sum(other_abun$soil),
  mean_all = sum(other_abun$mean_all),
  stringsAsFactors = FALSE
)
bar_data <- rbind(top_genera, other_row)

# Reshape for plotting
bar_long <- bar_data %>%
  pivot_longer(cols = all_of(plant_comps), names_to = "compartment", values_to = "abundance") %>%
  mutate(compartment = factor(compartment, levels = plant_comps),
         genus = factor(genus, levels = rev(bar_data$genus)))

figS4 <- ggplot(bar_long, aes(x = compartment, y = abundance, fill = genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_viridis_d() +
  labs(x = "Compartment", y = "Mean Relative Abundance",
       title = "Top 15 Genera by Compartment (16S)",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 7))
save_fig(figS4, "figS4_stacked_bar_16S", 10, 6)

log_msg("  Saved Fig S4")

# ==============================================================================
# Fig S5: Full Module-Taxon Correlation Heatmap
# ==============================================================================
log_msg("Generating Fig S5: Full module-taxon heatmap")

net_df <- read.csv(file.path(REPO_ROOT, "results/integration/networks/module_taxon_correlations.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)

# Show all correlations with p < 0.05 (not just padj < 0.05)
sig_net <- net_df[!is.na(net_df$p_value) & net_df$p_value < 0.05, ]

if (nrow(sig_net) > 0) {
  # Create label combining module and genus
  sig_net$module_genus <- paste0(sig_net$module, " | ", sig_net$Genus)

  figS5 <- ggplot(sig_net, aes(x = module, y = Genus, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_grid(amplicon ~ tissue, scales = "free", space = "free") +
    scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#E9ED4C",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(x = "WGCNA Module", y = "Taxon (Genus)",
         title = "Module-Taxon Correlations (p<0.05)",
         fill = "Spearman\nrho") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank())
  save_fig(figS5, "figS5_module_taxon_heatmap", 12, 10)
  log_msg("  Saved Fig S5 (", nrow(sig_net), " correlations)")
} else {
  log_msg("  No correlations for Fig S5")
}

# ==============================================================================
# Fig S6: MOFA+ Factor Weight Dot Plots
# ==============================================================================
log_msg("Generating Fig S6: MOFA+ factor weights")

weights_df <- read.csv(file.path(REPO_ROOT, "results/integration/factor_weights_top.tsv"),
                       sep = "\t", stringsAsFactors = FALSE)

# Top 10 features per factor per view
figS6 <- ggplot(weights_df, aes(x = reorder(feature, abs(weight)), y = weight, fill = view)) +
  geom_col(alpha = 0.8) +
  coord_flip() +
  facet_grid(view ~ factor, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("transcriptome" = "#0279EE",
                               "microbiome_16S" = "#75A025",
                               "microbiome_ITS" = "#FD9BED")) +
  labs(x = "Feature", y = "Weight",
       title = "MOFA+ Top Feature Weights per Factor",
       fill = "View") +
  theme(axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8),
        legend.position = "bottom")
save_fig(figS6, "figS6_mofa_factor_weights", 14, 8)

log_msg("  Saved Fig S6")

# ==============================================================================
# Fig S7: WGCNA Soft Power Scale-Free Topology Fit Curves
# ==============================================================================
log_msg("Generating Fig S7: WGCNA power curves")

compute_power_curve <- function(datExpr, powers = 1:20) {
  # Compute correlation matrix
  cor_mat <- cor(datExpr, use = "pairwise.complete.obs")

  sft <- data.frame(power = powers, sft_r2 = NA, mean_k = NA)
  for (p in powers) {
    # Adjacency
    adj <- abs(cor_mat)^p
    # Connectivity (excluding self)
    k <- colSums(adj) - 1
    # Scale-free topology fit
    hist_info <- hist(k, breaks = 50, plot = FALSE)
    k_bins <- hist_info$mids
    p_k <- hist_info$counts
    # Filter non-zero bins
    valid <- p_k > 0 & k_bins > 0
    if (sum(valid) >= 3) {
      fit <- lm(log10(p_k[valid]) ~ log10(k_bins[valid]))
      sft$sft_r2[sft$power == p] <- summary(fit)$r.squared
    }
    sft$mean_k[sft$power == p] <- mean(k)
  }
  return(sft)
}

power_curves <- list()
for (tissue in c("Leaf", "Adv-Root")) {
  wgcna_file <- file.path(REPO_ROOT, "results/rnaseq",
                          paste0("wgcna_", tissue),
                          paste0("wgcna_", tissue, ".rds"))
  if (file.exists(wgcna_file)) {
    w <- readRDS(wgcna_file)
    sft <- compute_power_curve(w$datExpr)
    sft$tissue <- tissue
    sft$chosen_power <- w$softPower
    power_curves[[tissue]] <- sft
    log_msg("  ", tissue, ": chosen power = ", w$softPower)
  }
}

if (length(power_curves) > 0) {
  power_df <- do.call(rbind, power_curves)
  power_df$tissue <- factor(power_df$tissue, levels = c("Leaf", "Adv-Root"))

  figS7 <- ggplot(power_df, aes(x = power, y = sft_r2, color = tissue)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(data = power_df[!duplicated(power_df$tissue), ],
               aes(xintercept = chosen_power, color = tissue),
               linetype = "dotted", linewidth = 0.5) +
    scale_color_manual(values = c("Leaf" = "#75A025", "Adv-Root" = "#FD9BED")) +
    labs(x = "Soft Power", y = "Scale-Free Topology Fit (R\u00b2)",
         title = "WGCNA Soft Power Selection",
         color = "Tissue") +
    theme(legend.position = "bottom")
  save_fig(figS7, "figS7_wgcna_power_curve", 7, 5)
  log_msg("  Saved Fig S7")
} else {
  log_msg("  No WGCNA data for Fig S7")
}

# ==============================================================================
# Copy PDFs to results
# ==============================================================================
for (p in list.files(PDF_DIR, full.names = TRUE)) {
  system(paste0("cp '", p, "' '", file.path(FIG_DIR, basename(p)), "'"))
}

log_msg("\n========================================")
log_msg("SUPPLEMENTARY FIGURES COMPLETE")
log_msg("Figures saved to: ", FIG_DIR)
log_msg("========================================")
