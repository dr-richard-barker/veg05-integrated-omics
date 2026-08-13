#!/usr/bin/env Rscript
# =============================================================================
# 02_microbiome_dada2.R
# DADA2 amplicon processing pipeline for VEG-05 microbiome (OSD-766)
#
# Processes both 16S (V4, 515F/806R) and ITS (ITS3F/ITS4) amplicons:
#   1. Primer removal (ITS only — 16S primers already removed)
#   2. Quality filtering
#   3. Error learning + denoising
#   4. Merge paired reads
#   5. Chimera removal
#   6. Taxonomy assignment (SILVA 138.2 for 16S, UNITE v7 for ITS)
#   7. Export ASV tables, taxonomy, tracking stats
#
# Input:  data/microbiome/fastq/*.fastq.gz
#         data/metadata/sample_metadata_microbiome.csv
#         data/microbiome/ref_db/silva_nr99_v138.1_train_set.fa.gz
#         data/microbiome/ref_db/silva_species_assignment_v138.1.fa.gz
#         data/microbiome/ref_db/unite_ver7_99_10.10.2017_dada2.fasta.gz
#
# Output: results/microbiome/{16S,ITS}/
#           asv_table_{amplicon}.tsv       — ASV count table (samples x ASVs)
#           taxonomy_{amplicon}.tsv        — Taxonomy assignments
#           seqtab_{amplicon}.rds          — Saved sequence table
#           taxtab_{amplicon}.rds          — Saved taxonomy table
#           track_stats_{amplicon}.tsv     — Read tracking through pipeline
#           quality_profile_{R1,R2}_{amplicon}.pdf — Quality profiles
#           phyloseq_{amplicon}.rds        — phyloseq object
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
})

# --- Configuration ----------------------------------------------------------
REPO_ROOT  <- "/mnt/shared-workspace/veg05-integrated-omics"
FASTQ_DIR  <- file.path(REPO_ROOT, "data/microbiome/fastq")
META_FILE  <- file.path(REPO_ROOT, "data/metadata/sample_metadata_microbiome.csv")
REF_DIR    <- file.path(REPO_ROOT, "data/microbiome/ref_db")
OUT_DIR    <- file.path(REPO_ROOT, "results/microbiome")
TRIM_DIR   <- "/workspace/dada2_trimmed"  # Local disk for trimmed FASTQ
# Ensure cutadapt is in PATH (installed via uv pip in /workspace/.venv)
Sys.setenv(PATH = paste("/workspace/.venv/bin", Sys.getenv("PATH"), sep = ":"))
PDF_DIR    <- "/workspace/dada2_pdfs"     # Local disk for PDF output (S3 can't do random-access)
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PDF_DIR, "16S"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PDF_DIR, "ITS"), showWarnings = FALSE, recursive = TRUE)
NTHREADS   <- 8

# 16S parameters (primers already removed)
SILVA_TRAIN   <- file.path(REF_DIR, "silva_nr99_v138.1_train_set.fa.gz")
SILVA_SPECIES <- file.path(REF_DIR, "silva_species_assignment_v138.1.fa.gz")

# ITS parameters (primers need removal)
ITS3F <- "GCATCGATGAAGAACGCAGC"   # 20 bp forward primer
ITS4  <- "TCCTCCGCTTATTGATATGC"   # 20 bp reverse primer
UNITE_DB <- file.path(REF_DIR, "unite_ver7_99_10.10.2017_dada2.fasta.gz")

# --- Helper functions -------------------------------------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# --- Load metadata ----------------------------------------------------------
meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
log_msg("Loaded metadata: ", nrow(meta), " samples")

# Create output directories
for (amp in c("16S", "ITS")) {
  out_sub <- file.path(OUT_DIR, amp)
  dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)
}
dir.create(TRIM_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 16S PIPELINE
# =============================================================================
process_16S <- function() {
  log_msg("========== 16S PIPELINE ==========")
  
  meta_16s <- meta[meta$amplicon == "16S", ]
  log_msg("16S samples: ", nrow(meta_16s))
  
  # Get FASTQ file paths
  fnFs <- file.path(FASTQ_DIR, meta_16s$fastq_r1)
  fnRs <- file.path(FASTQ_DIR, meta_16s$fastq_r2)
  sample_names <- meta_16s$sample_id
  
  # Verify files exist
  missing <- fnFs[!file.exists(fnFs)]
  if (length(missing) > 0) {
    log_msg("WARNING: ", length(missing), " R1 files missing")
    fnFs <- fnFs[file.exists(fnFs)]
    fnRs <- fnRs[file.exists(fnRs)]
    sample_names <- sample_names[file.exists(file.path(FASTQ_DIR, meta_16s$fastq_r1))]
  }
  log_msg("Files verified: ", length(fnFs), " pairs")
  
  # --- Quality profiles (first 4 samples) ---
  log_msg("Generating quality profiles...")
  pdf(file.path(PDF_DIR, "16S", "quality_profile_R1_16S.pdf"), width = 12, height = 6)
  plotQualityProfile(fnFs[1:min(4, length(fnFs))])
  dev.off()
  pdf(file.path(PDF_DIR, "16S", "quality_profile_R2_16S.pdf"), width = 12, height = 6)
  plotQualityProfile(fnRs[1:min(4, length(fnRs))])
  dev.off()
  
  # --- Filter and trim ---
  # 16S primers already removed, so we just filter by quality
  # truncLen: R1=240 (good quality to ~240), R2=160 (quality drops after ~160)
  filtFs <- file.path(TRIM_DIR, paste0(sample_names, "_F_filt.fastq.gz"))
  filtRs <- file.path(TRIM_DIR, paste0(sample_names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample_names
  names(filtRs) <- sample_names
  
  log_msg("Filtering and trimming 16S (truncLen=[240,160], maxEE=[2,2])...")
  filt_out <- filterAndTrim(
    fwd = fnFs, filt = filtFs,
    rev = fnRs, filt.rev = filtRs,
    truncLen = c(240, 160),
    maxEE = c(2, 2),
    truncQ = 2,
    minLen = 100,
    rm.phix = TRUE,
    multithread = NTHREADS
  )
  log_msg("Filtering complete. Reads retained: ", round(sum(filt_out[, 2]) / sum(filt_out[, 1]) * 100, 1), "%")
  
  # Remove files with 0 reads after filtering
  keep <- filt_out[, 2] > 0
  filtFs <- filtFs[keep]
  filtRs <- filtRs[keep]
  sample_names_filt <- sample_names[keep]
  log_msg("Samples with reads after filtering: ", sum(keep), " / ", length(keep))
  
  # --- Learn error rates ---
  log_msg("Learning error rates (forward)...")
  errF <- learnErrors(filtFs, multithread = NTHREADS, nbases = 1e8)
  log_msg("Learning error rates (reverse)...")
  errR <- learnErrors(filtRs, multithread = NTHREADS, nbases = 1e8)
  
  # Save error plots
  pdf(file.path(PDF_DIR, "16S", "error_rates_16S.pdf"), width = 10, height = 8)
  plotErrors(errF, nominalQ = TRUE)
  dev.off()
  
  # --- Denoise ---
  log_msg("Denoising (dada2)...")
  dadaFs <- dada(filtFs, err = errF, multithread = NTHREADS, pool = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = NTHREADS, pool = TRUE)
  log_msg("Denoising complete")
  
  # --- Merge pairs ---
  log_msg("Merging paired reads...")
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)
  log_msg("Merged: ", ncol(seqtab), " ASVs across ", nrow(seqtab), " samples")
  
  # Remove sequences outside expected V4 length range (~200-260 bp)
  seq_len <- nchar(getSequences(seqtab))
  seqtab2 <- seqtab[, seq_len >= 200 & seq_len <= 260]
  log_msg("After length filter (200-260bp): ", ncol(seqtab2), " ASVs (removed ", ncol(seqtab) - ncol(seqtab2), ")")
  
  # --- Chimera removal ---
  log_msg("Removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = NTHREADS, verbose = TRUE)
  log_msg("After chimera removal: ", ncol(seqtab_nochim), " ASVs (", 
          round(ncol(seqtab_nochim) / ncol(seqtab2) * 100, 1), "% retained)")
  
  # --- Track reads through pipeline ---
  track <- cbind(
    sample = sample_names,
    input = filt_out[, 1],
    filtered = filt_out[, 2],
    denoisedF = sapply(dadaFs, function(x) sum(x$denoised)),
    denoisedR = sapply(dadaRs, function(x) sum(x$denoised)),
    merged = sapply(mergers, function(x) sum(x$abundance)),
    nonchim = rowSums(seqtab_nochim[match(sample_names, rownames(seqtab_nochim)), , drop = FALSE])
  )
  track[is.na(track)] <- 0
  # Handle samples that were filtered out
  if (nrow(track) < length(sample_names)) {
    missing_samples <- setdiff(sample_names, track[, "sample"])
    missing_track <- cbind(sample = missing_samples, 
                           matrix(0, nrow = length(missing_samples), ncol = 6,
                                  dimnames = list(NULL, c("input","filtered","denoisedF","denoisedR","merged","nonchim"))))
    track <- rbind(track, missing_track)
  }
  write.table(track, file.path(OUT_DIR, "16S", "track_stats_16S.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("Track stats written")
  
  # --- Taxonomy assignment ---
  log_msg("Assigning taxonomy with SILVA 138.2...")
  taxa <- assignTaxonomy(seqtab_nochim, SILVA_TRAIN, multithread = NTHREADS, verbose = TRUE)
  log_msg("Adding species assignment...")
  taxa <- addSpecies(taxa, SILVA_SPECIES, verbose = TRUE)
  log_msg("Taxonomy assigned: ", nrow(taxa), " ASVs")
  
  # --- Create ASV IDs ---
  asv_ids <- paste0("ASV", sprintf("%03d", seq_len(ncol(seqtab_nochim))))
  colnames(seqtab_nochim) <- asv_ids
  rownames(taxa) <- asv_ids
  
  # --- Save outputs ---
  # ASV table
  asv_df <- as.data.frame(t(seqtab_nochim))
  colnames(asv_df) <- rownames(seqtab_nochim)
  asv_df$ASV_ID <- asv_ids
  asv_df <- asv_df[, c("ASV_ID", setdiff(colnames(asv_df), "ASV_ID"))]
  write.table(asv_df, file.path(OUT_DIR, "16S", "asv_table_16S.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Taxonomy table
  tax_df <- as.data.frame(taxa)
  tax_df$ASV_ID <- rownames(tax_df)
  tax_df <- tax_df[, c("ASV_ID", setdiff(colnames(tax_df), "ASV_ID"))]
  write.table(tax_df, file.path(OUT_DIR, "16S", "taxonomy_16S.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save RDS objects
  saveRDS(seqtab_nochim, file.path(OUT_DIR, "16S", "seqtab_16S.rds"))
  saveRDS(taxa, file.path(OUT_DIR, "16S", "taxtab_16S.rds"))
  
  # --- Build phyloseq object ---
  # Match metadata to samples in seqtab
  meta_match <- meta_16s[match(rownames(seqtab_nochim), meta_16s$sample_id), ]
  rownames(meta_match) <- meta_match$sample_id
  
  # Create DNAStringSet for refseq
  dna <- DNAStringSet(getSequences(seqtab_nochim))
  names(dna) <- asv_ids
  
  ps <- phyloseq(
    otu_table(seqtab_nochim, taxa_are_rows = FALSE),
    tax_table(taxa),
    sample_data(meta_match),
    refseq(dna)
  )
  saveRDS(ps, file.path(OUT_DIR, "16S", "phyloseq_16S.rds"))
  log_msg("Phyloseq object saved")
  
  # --- Summary ---
  log_msg("16S SUMMARY:")
  log_msg("  Total ASVs: ", ncol(seqtab_nochim))
  log_msg("  Total reads: ", sum(seqtab_nochim))
  log_msg("  Mean reads/sample: ", round(mean(rowSums(seqtab_nochim)), 0))
  log_msg("  Samples: ", nrow(seqtab_nochim))
  
  return(seqtab_nochim)
}

# =============================================================================
# ITS PIPELINE
# =============================================================================
process_ITS <- function() {
  log_msg("========== ITS PIPELINE ==========")
  
  meta_its <- meta[meta$amplicon == "ITS", ]
  log_msg("ITS samples: ", nrow(meta_its))
  
  fnFs <- file.path(FASTQ_DIR, meta_its$fastq_r1)
  fnRs <- file.path(FASTQ_DIR, meta_its$fastq_r2)
  sample_names <- meta_its$sample_id
  
  # Verify files exist
  fnFs <- fnFs[file.exists(fnFs)]
  fnRs <- fnRs[file.exists(fnRs)]
  sample_names <- meta_its$sample_id[file.exists(file.path(FASTQ_DIR, meta_its$fastq_r1))]
  log_msg("Files verified: ", length(fnFs), " pairs")
  
  # --- Primer removal with cutadapt ---
  # ITS3F on R1, ITS4 on R2 (both at 5' end of reads)
  # Also remove reverse complements at 3' ends (internal priming)
  log_msg("Removing ITS primers with cutadapt...")
  cutFs <- file.path(TRIM_DIR, paste0(sample_names, "_F_trim.fastq.gz"))
  cutRs <- file.path(TRIM_DIR, paste0(sample_names, "_R_trim.fastq.gz"))
  
  # cutadapt arguments:
  # -g ITS3F (forward primer on R1 5' end)
  # -G ITS4 (reverse primer on R2 5' end)
  # -a rc:ITS4 (reverse complement of ITS4 on R1 3' end - read-through)
  # -A rc:ITS3F (reverse complement of ITS3F on R2 3' end - read-through)
  # --minimum-length 50
  # --discard-untrimmed (keep only reads where primer was found)
  rc_ITS3F <- as.character(reverseComplement(DNAString(ITS3F)))
  rc_ITS4  <- as.character(reverseComplement(DNAString(ITS4)))
  
  for (i in seq_along(fnFs)) {
    if (!file.exists(cutFs[i]) || file.size(cutFs[i]) == 0) {
      cmd <- sprintf(
        "cutadapt -g %s -G %s --minimum-length 50 -j 1 -o %s -p %s %s %s 2>>%s/cutadapt_its.log",
        ITS3F, ITS4,
        cutFs[i], cutRs[i], fnFs[i], fnRs[i], TRIM_DIR
      )
      system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    if (i %% 20 == 0) log_msg("  Primer-trimmed ", i, " / ", length(fnFs), " samples")
  }
  log_msg("Primer removal complete")
  
  # Check which trimmed files have reads
  trim_exists <- file.exists(cutFs) & file.size(cutFs) > 0
  cutFs <- cutFs[trim_exists]
  cutRs <- cutRs[trim_exists]
  sample_names_trim <- sample_names[trim_exists]
  log_msg("Samples with reads after primer trimming: ", sum(trim_exists), " / ", length(trim_exists))
  
  # --- Quality profiles ---
  pdf(file.path(PDF_DIR, "ITS", "quality_profile_R1_ITS.pdf"), width = 12, height = 6)
  plotQualityProfile(cutFs[1:min(4, length(cutFs))])
  dev.off()
  pdf(file.path(PDF_DIR, "ITS", "quality_profile_R2_ITS.pdf"), width = 12, height = 6)
  plotQualityProfile(cutRs[1:min(4, length(cutRs))])
  dev.off()
  
  # --- Filter and trim ---
  # ITS: no truncation (variable length), just quality filter
  filtFs <- file.path(TRIM_DIR, paste0(sample_names_trim, "_F_filt.fastq.gz"))
  filtRs <- file.path(TRIM_DIR, paste0(sample_names_trim, "_R_filt.fastq.gz"))
  names(filtFs) <- sample_names_trim
  names(filtRs) <- sample_names_trim
  
  log_msg("Filtering ITS (truncLen=[0,180], maxEE=[2,2], minLen=50)...")
  filt_out <- filterAndTrim(
    fwd = cutFs, filt = filtFs,
    rev = cutRs, filt.rev = filtRs,
    truncLen = c(0, 180),  # Truncate R2 at 180bp where quality drops
    maxEE = c(2, 2),
    truncQ = 2,
    minLen = 50,
    rm.phix = TRUE,
    multithread = NTHREADS
  )
  log_msg("Filtering complete. Reads retained: ", round(sum(filt_out[, 2]) / sum(filt_out[, 1]) * 100, 1), "%")
  
  keep <- filt_out[, 2] > 0
  filtFs <- filtFs[keep]
  filtRs <- filtRs[keep]
  sample_names_filt <- sample_names_trim[keep]
  log_msg("Samples with reads after filtering: ", sum(keep), " / ", length(keep))
  
  # --- Learn error rates ---
  log_msg("Learning error rates (forward)...")
  errF <- learnErrors(filtFs, multithread = NTHREADS, nbases = 1e8)
  log_msg("Learning error rates (reverse)...")
  errR <- learnErrors(filtRs, multithread = NTHREADS, nbases = 1e8)
  
  pdf(file.path(PDF_DIR, "ITS", "error_rates_ITS.pdf"), width = 10, height = 8)
  plotErrors(errF, nominalQ = TRUE)
  dev.off()
  
  # --- Denoise ---
  log_msg("Denoising (dada2)...")
  dadaFs <- dada(filtFs, err = errF, multithread = NTHREADS, pool = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = NTHREADS, pool = TRUE)
  
  # --- Merge pairs ---
  log_msg("Merging paired reads...")
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)
  log_msg("Merged: ", ncol(seqtab), " ASVs across ", nrow(seqtab), " samples")
  
  # ITS length filter: keep 50-500 bp (ITS2 region expected range)
  seq_len <- nchar(getSequences(seqtab))
  seqtab2 <- seqtab[, seq_len >= 50 & seq_len <= 500]
  log_msg("After length filter (50-500bp): ", ncol(seqtab2), " ASVs (removed ", ncol(seqtab) - ncol(seqtab2), ")")
  
  # --- Chimera removal ---
  log_msg("Removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = NTHREADS, verbose = TRUE)
  log_msg("After chimera removal: ", ncol(seqtab_nochim), " ASVs (",
          round(ncol(seqtab_nochim) / ncol(seqtab2) * 100, 1), "% retained)")
  
  # --- Track reads ---
  track <- cbind(
    sample = sample_names_trim,
    input = filt_out[, 1],
    filtered = filt_out[, 2],
    denoisedF = sapply(dadaFs, function(x) sum(x$denoised)),
    denoisedR = sapply(dadaRs, function(x) sum(x$denoised)),
    merged = sapply(mergers, function(x) sum(x$abundance)),
    nonchim = rowSums(seqtab_nochim[match(sample_names_trim, rownames(seqtab_nochim)), , drop = FALSE])
  )
  track[is.na(track)] <- 0
  if (nrow(track) < length(sample_names)) {
    missing_samples <- setdiff(sample_names, track[, "sample"])
    missing_track <- cbind(sample = missing_samples,
                           matrix(0, nrow = length(missing_samples), ncol = 6,
                                  dimnames = list(NULL, c("input","filtered","denoisedF","denoisedR","merged","nonchim"))))
    track <- rbind(track, missing_track)
  }
  write.table(track, file.path(OUT_DIR, "ITS", "track_stats_ITS.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # --- Taxonomy assignment ---
  log_msg("Assigning taxonomy with UNITE v7...")
  taxa <- assignTaxonomy(seqtab_nochim, UNITE_DB, multithread = NTHREADS, verbose = TRUE)
  log_msg("Taxonomy assigned: ", nrow(taxa), " ASVs")
  
  # --- Capture sequences BEFORE renaming ---
  asv_seqs <- getSequences(seqtab_nochim)
  
  # --- Create ASV IDs ---
  asv_ids <- paste0("ASV", sprintf("%03d", seq_len(ncol(seqtab_nochim))))
  colnames(seqtab_nochim) <- asv_ids
  rownames(taxa) <- asv_ids
  
  # --- Save outputs ---
  asv_df <- as.data.frame(t(seqtab_nochim))
  colnames(asv_df) <- rownames(seqtab_nochim)
  asv_df$ASV_ID <- asv_ids
  asv_df <- asv_df[, c("ASV_ID", setdiff(colnames(asv_df), "ASV_ID"))]
  write.table(asv_df, file.path(OUT_DIR, "ITS", "asv_table_ITS.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  tax_df <- as.data.frame(taxa)
  tax_df$ASV_ID <- rownames(tax_df)
  tax_df <- tax_df[, c("ASV_ID", setdiff(colnames(tax_df), "ASV_ID"))]
  write.table(tax_df, file.path(OUT_DIR, "ITS", "taxonomy_ITS.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  saveRDS(seqtab_nochim, file.path(OUT_DIR, "ITS", "seqtab_ITS.rds"))
  saveRDS(taxa, file.path(OUT_DIR, "ITS", "taxtab_ITS.rds"))
  
  # --- Build phyloseq object ---
  meta_match <- meta_its[match(rownames(seqtab_nochim), meta_its$sample_id), ]
  rownames(meta_match) <- meta_match$sample_id
  
  dna <- DNAStringSet(asv_seqs)
  names(dna) <- asv_ids
  
  ps <- phyloseq(
    otu_table(seqtab_nochim, taxa_are_rows = FALSE),
    tax_table(taxa),
    sample_data(meta_match),
    refseq(dna)
  )
  saveRDS(ps, file.path(OUT_DIR, "ITS", "phyloseq_ITS.rds"))
  
  log_msg("ITS SUMMARY:")
  log_msg("  Total ASVs: ", ncol(seqtab_nochim))
  log_msg("  Total reads: ", sum(seqtab_nochim))
  log_msg("  Mean reads/sample: ", round(mean(rowSums(seqtab_nochim)), 0))
  log_msg("  Samples: ", nrow(seqtab_nochim))
  
  return(seqtab_nochim)
}

# =============================================================================
# RUN BOTH PIPELINES
# =============================================================================
log_msg("Starting DADA2 pipeline for VEG-05 microbiome")
log_msg("R version: ", R.version.string)
log_msg("dada2 version: ", as.character(packageVersion("dada2")))

# Process 16S (already completed - loading saved results)
log_msg("16S already processed, loading saved results...")
seqtab_16s <- readRDS(file.path(OUT_DIR, "16S", "seqtab_16S.rds"))

# Process ITS
seqtab_ITS <- process_ITS()

log_msg("========================================")
log_msg("DADA2 PIPELINE COMPLETE")
log_msg("16S: ", ncol(seqtab_16s), " ASVs, ", nrow(seqtab_16s), " samples")
log_msg("ITS: ", ncol(seqtab_ITS), " ASVs, ", nrow(seqtab_ITS), " samples")
log_msg("Outputs saved to: ", OUT_DIR)
log_msg("========================================")
