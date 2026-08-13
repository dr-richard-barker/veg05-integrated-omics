#!/usr/bin/env Rscript
# ==============================================================================
# Step 5: Functional Prediction
#
# 1. FAPROTAX — functional categories for 16S bacterial communities
#    (calls scripts/faprotax_collapse.py, a custom Python implementation
#     because the original collapse_table.py was unavailable for download)
# 2. FUNGuild — ecological guild assignment for ITS fungal communities
#    (attempts FUNGuild.py; falls back to manual genus-level assignment
#     if the script or database is unavailable)
#
# Inputs:
#   results/microbiome/16S/asv_table_16S.tsv
#   results/microbiome/16S/taxonomy_16S.tsv
#   results/microbiome/ITS/asv_table_ITS.tsv
#   results/microbiome/ITS/taxonomy_ITS.tsv
#   data/microbiome/ref_db/FAPROTAX.txt
#
# Outputs:
#   results/microbiome/functional/faprotax_functions_16S.tsv
#   results/microbiome/functional/faprotax_report_16S.tsv
#   results/microbiome/functional/funguild_guilds_ITS_manual.tsv
#   results/microbiome/functional/funguild_summary_ITS.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
})

# --- Paths ---
REPO_ROOT  <- "/mnt/shared-workspace/veg05-integrated-omics"
FUNC_DIR   <- file.path(REPO_ROOT, "results/microbiome/functional")
SCRIPT_DIR <- file.path(REPO_ROOT, "scripts")

dir.create(FUNC_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

# ==============================================================================
# Part 1: FAPROTAX Functional Prediction (16S)
# ==============================================================================
log_msg("=== Part 1: FAPROTAX Functional Prediction (16S) ===")

faprotax_db   <- file.path(REPO_ROOT, "data/microbiome/ref_db/FAPROTAX.txt")
collapse_py   <- file.path(SCRIPT_DIR, "faprotax_collapse.py")

# Load taxonomy
tax_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/16S/taxonomy_16S.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)

if (file.exists(faprotax_db) && file.exists(collapse_py)) {
  log_msg("  FAPROTAX database: ", faprotax_db)
  log_msg("  collapse script:   ", collapse_py)

  # Prepare 16S ASV table for FAPROTAX
  asv_16s <- read.csv(file.path(REPO_ROOT, "results/microbiome/16S/asv_table_16S.tsv"),
                      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  asv_ids <- asv_16s[, 1]
  count_data <- asv_16s[, -1]
  rownames(count_data) <- asv_ids

  # Build taxonomy string: Kingdom;Phylum;Class;Order;Family;Genus;Species
  tax_lookup <- tax_16s
  rownames(tax_lookup) <- tax_lookup$ASV_ID

  tax_strings <- sapply(asv_ids, function(a) {
    t <- tax_lookup[a, ]
    parts <- c(t$Kingdom, t$Phylum, t$Class, t$Order, t$Family, t$Genus, t$Species)
    parts[is.na(parts)] <- ""
    paste(parts, collapse = ";")
  })

  # FAPROTAX input: rows = taxa, cols = samples, first column = taxonomy
  faprotax_input <- data.frame(
    taxonomy = tax_strings,
    t(count_data),
    check.names = FALSE
  )
  rownames(faprotax_input) <- asv_ids

  faprotax_infile <- "/workspace/faprotax_input_16S.tsv"
  write.table(faprotax_input, faprotax_infile, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("  FAPROTAX input: ", nrow(faprotax_input), " ASVs, ", ncol(faprotax_input) - 1, " samples")

  # Run custom collapse script
  faprotax_outdir <- "/workspace/faprotax_output"
  dir.create(faprotax_outdir, showWarnings = FALSE)

  cmd <- paste0(
    "python3 '", collapse_py, "' ",
    "-i '", faprotax_infile, "' ",
    "-o '", file.path(faprotax_outdir, "faprotax_report.tsv"), "' ",
    "-g '", faprotax_db, "' ",
    " 2>&1"
  )
  log_msg("  Running FAPROTAX collapse...")
  faprotax_result <- system(cmd, intern = TRUE)
  cat(paste(faprotax_result, collapse = "\n"), "\n")

  # Read and save output
  faprotax_report_file <- file.path(faprotax_outdir, "faprotax_report.tsv")
  if (file.exists(faprotax_report_file)) {
    fap_res <- read.csv(faprotax_report_file, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE)
    write.table(fap_res, file.path(FUNC_DIR, "faprotax_functions_16S.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    log_msg("  FAPROTAX: ", nrow(fap_res), " functional groups predicted")

    # Build summary (total reads per function)
    if (ncol(fap_res) > 2) {
      func_col <- colnames(fap_res)[1]
      sums <- rowSums(fap_res[, 3:ncol(fap_res)], na.rm = TRUE)
      func_summary <- data.frame(
        `function` = fap_res[, func_col],
        total_reads = sums,
        check.names = FALSE
      )
      func_summary <- func_summary[order(func_summary$total_reads, decreasing = TRUE), ]
      write.table(func_summary, file.path(FUNC_DIR, "faprotax_report_16S.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      log_msg("  Top functions:")
      for (i in 1:min(15, nrow(func_summary))) {
        if (func_summary$total_reads[i] > 0) {
          log_msg("    ", func_summary[, 1][i], ": ", func_summary$total_reads[i])
        }
      }
    }
  } else {
    log_msg("  WARNING: FAPROTAX output not found")
  }
} else {
  log_msg("  WARNING: FAPROTAX files not found, skipping")
  log_msg("  FAPROTAX.txt exists: ", file.exists(faprotax_db))
  log_msg("  faprotax_collapse.py exists: ", file.exists(collapse_py))
}

# ==============================================================================
# Part 2: FUNGuild Ecological Guild Prediction (ITS)
# ==============================================================================
log_msg("\n=== Part 2: FUNGuild Ecological Guild Prediction (ITS) ===")

funguild_script <- file.path(REPO_ROOT, "data/microbiome/ref_db/FUNGuild.py")
tax_its <- read.csv(file.path(REPO_ROOT, "results/microbiome/ITS/taxonomy_ITS.tsv"),
                    sep = "\t", stringsAsFactors = FALSE)

# Load ITS ASV table
asv_its <- read.csv(file.path(REPO_ROOT, "results/microbiome/ITS/asv_table_ITS.tsv"),
                    sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
asv_ids_its <- asv_its[, 1]
count_data_its <- asv_its[, -1]
rownames(count_data_its) <- asv_ids_its
total_counts <- rowSums(count_data_its)

# Build taxonomy strings in FUNGuild format
tax_lookup_its <- tax_its
rownames(tax_lookup_its) <- tax_lookup_its$ASV_ID

tax_strings_its <- sapply(asv_ids_its, function(a) {
  t <- tax_lookup_its[a, ]
  parts <- c(
    paste0("k__", ifelse(is.na(t$Kingdom), "", t$Kingdom)),
    paste0("p__", ifelse(is.na(t$Phylum), "", t$Phylum)),
    paste0("c__", ifelse(is.na(t$Class), "", t$Class)),
    paste0("o__", ifelse(is.na(t$Order), "", t$Order)),
    paste0("f__", ifelse(is.na(t$Family), "", t$Family)),
    paste0("g__", ifelse(is.na(t$Genus), "", t$Genus)),
    paste0("s__", ifelse(is.na(t$Species), "", t$Species))
  )
  paste(parts, collapse = ";")
})

# Manual guild assignment based on genus (used as fallback or primary method)
genus_guilds <- c(
  "Fusarium" = "Plant Pathogen",
  "Penicillium" = "Saprotroph",
  "Aspergillus" = "Saprotroph",
  "Trichoderma" = "Saprotroph",
  "Cladosporium" = "Leaf Saprotroph",
  "Alternaria" = "Plant Pathogen",
  "Botrytis" = "Plant Pathogen",
  "Mucor" = "Saprotroph",
  "Rhizopus" = "Saprotroph",
  "Chaetomium" = "Saprotroph",
  "Acremonium" = "Saprotroph",
  "Verticillium" = "Plant Pathogen",
  "Colletotrichum" = "Plant Pathogen"
)

assign_manual_guilds <- function() {
  fg_manual <- data.frame(
    ASV_ID = asv_ids_its,
    taxonomy = tax_strings_its,
    abundance = total_counts,
    genus = sapply(asv_ids_its, function(a) {
      t <- tax_lookup_its[a, "Genus"]
      if (is.na(t)) return(NA)
      return(t)
    }),
    stringsAsFactors = FALSE
  )
  fg_manual$guild <- sapply(fg_manual$genus, function(g) {
    if (is.na(g)) return("Unassigned")
    if (g %in% names(genus_guilds)) return(genus_guilds[g])
    return("Unassigned")
  })
  fg_manual$confidence <- ifelse(fg_manual$guild != "Unassigned", "Probable", "NA")

  write.table(fg_manual, file.path(FUNC_DIR, "funguild_guilds_ITS_manual.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("  Manual guild assignment: ", nrow(fg_manual), " ASVs")

  guild_table <- table(fg_manual$guild)
  guild_df <- data.frame(guild = names(guild_table), count = as.numeric(guild_table))
  guild_df <- guild_df[order(guild_df$count, decreasing = TRUE), ]
  log_msg("  Guild summary:")
  for (i in 1:min(10, nrow(guild_df))) {
    log_msg("    ", guild_df$guild[i], ": ", guild_df$count[i])
  }
  write.table(guild_df, file.path(FUNC_DIR, "funguild_summary_ITS.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Try FUNGuild script first; fall back to manual
if (file.exists(funguild_script)) {
  log_msg("  FUNGuild script found: ", funguild_script)

  funguild_input <- data.frame(
    OTU = asv_ids_its,
    taxonomy = tax_strings_its,
    abundance = total_counts,
    stringsAsFactors = FALSE
  )
  funguild_infile <- "/workspace/funguild_input_ITS.tsv"
  write.table(funguild_input, funguild_infile, sep = "\t",
              row.names = FALSE, quote = FALSE)
  log_msg("  FUNGuild input: ", nrow(funguild_input), " ASVs")

  funguild_outfile <- "/workspace/funguild_output_ITS.tsv"
  # Note: -db UNITE is invalid in FUNGuild v1.1 (choices: fungi, nematode)
  # Use -db fungi
  cmd <- paste0("python3 '", funguild_script, "' -otu '", funguild_infile,
                "' -db fungi -output '", funguild_outfile, "' 2>&1")
  log_msg("  Running FUNGuild...")
  funguild_result <- system(cmd, intern = TRUE)
  cat(paste(funguild_result, collapse = "\n"), "\n")

  if (file.exists(funguild_outfile)) {
    fg_res <- read.csv(funguild_outfile, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
    write.table(fg_res, file.path(FUNC_DIR, "funguild_guilds_ITS.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    log_msg("  FUNGuild: ", nrow(fg_res), " ASVs with guild assignments")

    if ("guild" %in% colnames(fg_res)) {
      guild_table <- table(fg_res$guild)
      guild_df <- data.frame(guild = names(guild_table), count = as.numeric(guild_table))
      guild_df <- guild_df[order(guild_df$count, decreasing = TRUE), ]
      write.table(guild_df, file.path(FUNC_DIR, "funguild_summary_ITS.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  } else {
    log_msg("  WARNING: FUNGuild output not found, using manual guild assignment")
    assign_manual_guilds()
  }
} else {
  log_msg("  FUNGuild script not found, using manual genus-level assignment")
  assign_manual_guilds()
}

# ==============================================================================
# Summary
# ==============================================================================
log_msg("\n========================================")
log_msg("FUNCTIONAL PREDICTION COMPLETE")
log_msg("Outputs:")
log_msg("  FAPROTAX: ", file.path(FUNC_DIR, "faprotax_functions_16S.tsv"))
log_msg("  FAPROTAX report: ", file.path(FUNC_DIR, "faprotax_report_16S.tsv"))
log_msg("  FUNGuild: ", file.path(FUNC_DIR, "funguild_guilds_ITS_manual.tsv"))
log_msg("  FUNGuild summary: ", file.path(FUNC_DIR, "funguild_summary_ITS.tsv"))
log_msg("========================================")
