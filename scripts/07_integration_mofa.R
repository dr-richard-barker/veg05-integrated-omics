#!/usr/bin/env Rscript
# =============================================================================
# 07_integration_mofa.R
# MOFA+ unsupervised multi-omics integration for VEG-05
#
# Integrates three views for Leaf tissue (best sample overlap):
#   1. Transcriptome: top 2000 variable genes (VST)
#   2. Microbiome 16S: ASV relative abundances (bacterial, chloro/mito removed)
#   3. Microbiome ITS: ASV relative abundances (fungal)
#
# MOFA+ learns latent factors that capture coordinated variation across views,
# then factors are correlated with traits (flight, light, dysbiosis).
#
# Input:  results/rnaseq/vst_Leaf.tsv
#         results/microbiome/16S/phyloseq_16S_filtered.rds
#         results/microbiome/ITS/phyloseq_ITS.rds
#         data/metadata/sample_metadata_rnaseq.csv
#         data/metadata/sample_crosswalk.csv
#
# Output: results/integration/
#           mofa_model.rds
#           factor_variance_explained.tsv
#           factor_trait_correlations.tsv
#           factor_weights_top.tsv
#           mofa_factor_values.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(phyloseq)
  library(ggplot2)
})
# Point reticulate to the Python with mofapy2 installed
reticulate::use_python("/workspace/.venv/bin/python")

# --- Configuration ----------------------------------------------------------
REPO_ROOT  <- "/mnt/shared-workspace/veg05-integrated-omics"
OUT_DIR    <- file.path(REPO_ROOT, "results/integration")
PDF_DIR    <- "/workspace/mofa_pdfs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N_GENES <- 2000
N_FACTORS   <- 5  # Reduced from 15 due to small sample size (n=21)
TISSUE      <- "Leaf"  # Focus on leaf (best sample overlap)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# --- Load data --------------------------------------------------------------
log_msg("Loading data for MOFA+ integration (", TISSUE, ")")

# 1. Transcriptome: VST matrix
vst <- read.csv(file.path(REPO_ROOT, "results/rnaseq", paste0("vst_", TISSUE, ".tsv")),
                sep = "\t", row.names = 1, check.names = FALSE)
vst <- as.matrix(vst)
log_msg("Transcriptome: ", nrow(vst), " genes x ", ncol(vst), " samples")

# 2. 16S microbiome: filtered phyloseq (chloro/mito removed)
ps16 <- readRDS(file.path(REPO_ROOT, "results/microbiome/16S/phyloseq_16S_filtered.rds"))

# 3. ITS microbiome
psITS <- readRDS(file.path(REPO_ROOT, "results/microbiome/ITS/phyloseq_ITS.rds"))
psITS <- prune_samples(sample_sums(psITS) > 0, psITS)

# --- Build sample crosswalk -------------------------------------------------
# RNA-seq sample IDs: VEG-05-Flt-SN05-Leaf-Red_L008
# Microbiome sample IDs: VEG-05F-SN05-leaf-tom-red-rich_S107_L001
# Need to match by plant + flight + light + tissue

meta_rna <- read.csv(file.path(REPO_ROOT, "data/metadata/sample_metadata_rnaseq.csv"),
                      stringsAsFactors = FALSE)
meta_rna <- meta_rna[meta_rna$tissue == TISSUE, ]

# Build matching key
meta_rna$match_key <- paste(meta_rna$flight, meta_rna$plant, meta_rna$light, sep = "_")

# For microbiome, build matching key from sample metadata
sd16 <- data.frame(sample_data(ps16))
sd16$match_key <- paste(sd16$flight, sd16$plant, sd16$light, sep = "_")
sd16$match_key <- gsub("-", "", sd16$match_key)  # Normalize SN-01 -> SN01

# Same for ITS
sdITS <- data.frame(sample_data(psITS))
sdITS$match_key <- paste(sdITS$flight, sdITS$plant, sdITS$light, sep = "_")
sdITS$match_key <- gsub("-", "", sdITS$match_key)

# Normalize RNA-seq match key too
meta_rna$match_key <- gsub("-", "", meta_rna$match_key)

# Find samples present in all three views
common_keys <- intersect(
  meta_rna$match_key,
  intersect(sd16$match_key, sdITS$match_key)
)
log_msg("Common samples (all 3 views): ", length(common_keys))

# Build matched sample sets
rna_samples <- sapply(common_keys, function(k) meta_rna$sample_id[meta_rna$match_key == k][1])
names(rna_samples) <- common_keys

# For microbiome, pick one sample per key (first match)
get_microbiome_sample <- function(sd, ps, key) {
  matches <- rownames(sd)[sd$match_key == key]
  if (length(matches) == 0) return(NULL)
  # Pick the one with most reads
  reads <- sample_sums(ps)[matches]
  names(which.max(reads))
}

s16_samples <- sapply(common_keys, function(k) get_microbiome_sample(sd16, ps16, k))
sITS_samples <- sapply(common_keys, function(k) get_microbiome_sample(sdITS, psITS, k))

log_msg("Matched: ", length(rna_samples), " RNA-seq, ", 
        sum(!sapply(s16_samples, is.null)), " 16S, ",
        sum(!sapply(sITS_samples, is.null)), " ITS")

# --- Prepare view matrices --------------------------------------------------
# Use common sample names (the match_key) for MOFA+
sample_names <- common_keys

# View 1: Transcriptome (top variable genes)
gene_var <- apply(vst, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(TOP_N_GENES, nrow(vst))]
rna_mat <- vst[top_genes, rna_samples, drop = FALSE]
colnames(rna_mat) <- sample_names
# Center and scale
rna_mat <- t(scale(t(rna_mat)))
log_msg("Transcriptome view: ", nrow(rna_mat), " genes x ", ncol(rna_mat), " samples")

# View 2: 16S microbiome (relative abundance, CLR-transformed)
ps16_sub <- prune_samples(na.omit(s16_samples), ps16)
otu16 <- as(otu_table(ps16_sub), "matrix")
if (!taxa_are_rows(ps16_sub)) otu16 <- t(otu16)
# Use the matched sample names
colnames(otu16) <- names(s16_samples)[!is.na(s16_samples)]
# Relative abundance
otu16_rel <- sweep(otu16, 2, colSums(otu16), "/")
otu16_rel[is.na(otu16_rel)] <- 0
# CLR transform with pseudocount
otu16_clr <- log2(otu16_rel + 1e-6) - rowMeans(log2(otu16_rel + 1e-6))
# Scale
otu16_clr <- t(scale(t(otu16_clr)))
log_msg("16S view: ", nrow(otu16_clr), " ASVs x ", ncol(otu16_clr), " samples")

# View 3: ITS microbiome
psITS_sub <- prune_samples(na.omit(sITS_samples), psITS)
otuITS <- as(otu_table(psITS_sub), "matrix")
if (!taxa_are_rows(psITS_sub)) otuITS <- t(otuITS)
colnames(otuITS) <- names(sITS_samples)[!is.na(sITS_samples)]
otuITS_rel <- sweep(otuITS, 2, colSums(otuITS), "/")
otuITS_rel[is.na(otuITS_rel)] <- 0
otuITS_clr <- log2(otuITS_rel + 1e-6) - rowMeans(log2(otuITS_rel + 1e-6))
otuITS_clr <- t(scale(t(otuITS_clr)))
log_msg("ITS view: ", nrow(otuITS_clr), " ASVs x ", ncol(otuITS_clr), " samples")

# --- Create MOFA object -----------------------------------------------------
log_msg("Creating MOFA+ object...")

# Ensure all views have the same samples
common_samples <- intersect(colnames(rna_mat), 
                            intersect(colnames(otu16_clr), colnames(otuITS_clr)))
log_msg("Final common samples: ", length(common_samples))

rna_mat <- rna_mat[, common_samples, drop = FALSE]
otu16_clr <- otu16_clr[, common_samples, drop = FALSE]
otuITS_clr <- otuITS_clr[, common_samples, drop = FALSE]

# Create data list
data_list <- list(
  transcriptome = rna_mat,
  microbiome_16S = otu16_clr,
  microbiome_ITS = otuITS_clr
)

# Create MOFA object
mofa <- create_mofa(data_list)
# Get defaults and modify
train_opts <- get_default_training_options(mofa)
train_opts$maxiter <- 1000
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
train_opts$verbose <- FALSE

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- N_FACTORS

mofa <- prepare_mofa(
  mofa,
  model_options = model_opts,
  training_options = train_opts
)

# --- Train MOFA ---
log_msg("Training MOFA+ model...")
# HDF5 requires random-access writes — save to /workspace first, then copy to results
mofa_hdf5_local <- "/workspace/mofa_model.hdf5"
if (file.exists(mofa_hdf5_local)) file.remove(mofa_hdf5_local)
mofa <- run_mofa(mofa, outfile = mofa_hdf5_local)
system(paste0("cp ", mofa_hdf5_local, " ", file.path(OUT_DIR, "mofa_model.hdf5")))
log_msg("Training complete")

# --- Extract results --------------------------------------------------------
# Variance explained per factor per view
ve <- calculate_variance_explained(mofa)
# r2_per_factor can be a list (per view) or nested differently; handle robustly
r2pf <- ve$r2_per_factor
cat("Variance explained structure:\n")
str(r2pf, max.level = 2)

# Build data frame robustly
# r2_per_factor is a relistable list; actual data in $group1 as matrix (factors x views)
r2mat <- r2pf$group1
if (is.null(r2mat)) r2mat <- r2pf[[1]]
colnames(r2mat) <- views_names(mofa)
rownames(r2mat) <- paste0("Factor", 1:nrow(r2mat))

ve_df <- data.frame(
  factor = rownames(r2mat),
  r2_transcriptome = r2mat[, "transcriptome"],
  r2_microbiome_16S = r2mat[, "microbiome_16S"],
  r2_microbiome_ITS = r2mat[, "microbiome_ITS"]
)
ve_df$r2_total <- rowMeans(ve_df[, -1], na.rm = TRUE)
write.table(ve_df, file.path(OUT_DIR, "factor_variance_explained.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

log_msg("Variance explained per factor:")
for (i in 1:min(N_FACTORS, nrow(ve_df))) {
  log_msg("  Factor", i, ": RNA=", round(ve_df$r2_transcriptome[i]*100, 1), "%, ",
          "16S=", round(ve_df$r2_microbiome_16S[i]*100, 1), "%, ",
          "ITS=", round(ve_df$r2_microbiome_ITS[i]*100, 1), "%")
}

# Factor values (samples x factors)
factors <- get_factors(mofa)
factor_df <- as.data.frame(factors[[1]])
factor_df$sample_id <- rownames(factor_df)

# Add traits
meta_match <- meta_rna[match(common_samples, meta_rna$match_key), ]
factor_df$flight <- meta_match$flight
factor_df$light <- meta_match$light

# Add dysbiosis indices
dys_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/dysbiosis_index_16S.tsv"), sep="\t")
dys_its <- read.csv(file.path(REPO_ROOT, "results/microbiome/community_health/dysbiosis_index_ITS.tsv"), sep="\t")

for (i in 1:nrow(factor_df)) {
  key <- common_samples[i]
  # Find dysbiosis by matching
  plant <- meta_match$plant[i]
  flight <- meta_match$flight[i]
  plant_norm <- gsub("-", "", plant)
  
  d16 <- dys_16s[dys_16s$compartment == "leaf" & dys_16s$flight == flight &
                  grepl(plant_norm, gsub("-", "", dys_16s$sample_id)), ]
  if (nrow(d16) > 0) factor_df$dysbiosis_16S[i] <- d16$dysbiosis_index[1] else factor_df$dysbiosis_16S[i] <- NA
  
  dITS <- dys_its[dys_its$compartment == "leaf" & dys_its$flight == flight &
                   grepl(plant_norm, gsub("-", "", dys_its$sample_id)), ]
  if (nrow(dITS) > 0) factor_df$dysbiosis_ITS[i] <- dITS$dysbiosis_index[1] else factor_df$dysbiosis_ITS[i] <- NA
}

write.table(factor_df, file.path(OUT_DIR, "mofa_factor_values.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Factor-trait correlations
trait_cols <- c("flight", "light", "dysbiosis_16S", "dysbiosis_ITS")
factor_cols <- paste0("Factor", 1:N_FACTORS)

cor_results <- list()
for (fc in factor_cols) {
  for (tc in trait_cols) {
    if (tc %in% c("flight", "light")) {
      # Binary trait
      trait_vec <- as.numeric(factor(factor_df[[tc]], levels = c("Ground", "Flight"))) - 1
      if (tc == "light") trait_vec <- as.numeric(factor(factor_df[[tc]], levels = c("Red", "Blue"))) - 1
    } else {
      trait_vec <- factor_df[[tc]]
    }
    
    factor_vec <- factor_df[[fc]]
    valid <- !is.na(trait_vec) & !is.na(factor_vec)
    if (sum(valid) < 3) next
    
    r <- cor(factor_vec[valid], trait_vec[valid], method = "spearman")
    p <- cor.test(factor_vec[valid], trait_vec[valid], method = "spearman")$p.value
    
    cor_results[[length(cor_results) + 1]] <- data.frame(
      factor = fc,
      trait = tc,
      rho = r,
      p_value = p,
      n = sum(valid)
    )
  }
}

cor_df <- do.call(rbind, cor_results)
cor_df$padj <- p.adjust(cor_df$p_value, method = "BH")
write.table(cor_df, file.path(OUT_DIR, "factor_trait_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

log_msg("\nSignificant factor-trait correlations (padj<0.05):")
sig <- cor_df[cor_df$padj < 0.05, ]
if (nrow(sig) > 0) {
  for (i in 1:nrow(sig)) {
    log_msg("  ", sig$factor[i], " ~ ", sig$trait[i], ": rho=", round(sig$rho[i], 2),
            ", padj=", signif(sig$padj[i], 3))
  }
} else {
  log_msg("  None at padj<0.05")
  log_msg("  Top unadjusted:")
  cor_df <- cor_df[order(cor_df$p_value), ]
  for (i in 1:min(10, nrow(cor_df))) {
    log_msg("  ", cor_df$factor[i], " ~ ", cor_df$trait[i], ": rho=", round(cor_df$rho[i], 2),
            ", p=", signif(cor_df$p_value[i], 3), ", padj=", signif(cor_df$padj[i], 3))
  }
}

# --- Top feature weights per factor ---
log_msg("Extracting top feature weights...")
weights_list <- list()
for (view in names(data_list)) {
  w <- get_weights(mofa, view = view)
  # w is a list: $group1 = matrix (features x factors)
  wmat <- w$group1
  if (is.null(wmat)) wmat <- w[[1]]
  for (f in 1:N_FACTORS) {
    if (f > ncol(wmat)) next
    wf <- wmat[, f]
    wf <- wf[!is.na(wf)]
    if (length(wf) == 0) next
    top_pos <- names(sort(wf, decreasing = TRUE))[1:min(10, length(wf))]
    top_neg <- names(sort(wf, decreasing = FALSE))[1:min(10, length(wf))]

    for (feat in c(top_pos, top_neg)) {
      weights_list[[length(weights_list) + 1]] <- data.frame(
        view = view,
        factor = paste0("Factor", f),
        feature = feat,
        weight = wf[feat]
      )
    }
  }
}

weights_df <- do.call(rbind, weights_list)
weights_df <- weights_df[order(abs(weights_df$weight), decreasing = TRUE), ]
write.table(weights_df, file.path(OUT_DIR, "factor_weights_top.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Save model ---
saveRDS(mofa, file.path(OUT_DIR, "mofa_model.rds"))

# --- Plots ---
pdf(file.path(PDF_DIR, "mofa_variance_explained.pdf"), width = 10, height = 6)
plot_variance_explained(mofa, plot_total = TRUE)
dev.off()

pdf(file.path(PDF_DIR, "mofa_factor_correlation.pdf"), width = 8, height = 7)
plot_factor_cor(mofa)
dev.off()

# Copy PDFs (use system cp — R file.copy fails on S3-backed filesystem)
for (p in list.files(PDF_DIR, full.names = TRUE)) {
  system(paste0("cp '", p, "' '", file.path(OUT_DIR, basename(p)), "'"))
}

log_msg("\n========================================")
log_msg("MOFA+ INTEGRATION COMPLETE")
log_msg("Samples: ", length(common_samples))
log_msg("Views: transcriptome (", nrow(rna_mat), " genes), 16S (", nrow(otu16_clr), 
        " ASVs), ITS (", nrow(otuITS_clr), " ASVs)")
log_msg("Factors: ", N_FACTORS)
log_msg("Outputs saved to: ", OUT_DIR)
log_msg("========================================")
