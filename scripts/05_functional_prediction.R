#!/usr/bin/env Rscript
# ==============================================================================
# Step 5: Functional Prediction
#
# 1. FAPROTAX — functional categories for 16S bacterial communities, using
#    FAPROTAX 1.2.12 and its own collapse_table.py (Louca et al. 2016;
#    https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/)
# 2. FUNGuild — ecological guild assignment for ITS fungal communities
#    (attempts FUNGuild.py; falls back to manual genus-level assignment
#     if the script or database is unavailable)
#
# Inputs:
#   results/microbiome/16S/phyloseq_16S_filtered.rds (chloroplast/mitochondria
#     removed, written by 04_community_health.R)
#   results/microbiome/ITS/asv_table_ITS.tsv
#   results/microbiome/ITS/taxonomy_ITS.tsv
#   data/microbiome/ref_db/FAPROTAX_1.2.12/{FAPROTAX.txt,collapse_table.py}
#     (unmodified FAPROTAX 1.2.12 release)
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

# FAPROTAX's own collapse_table.py matches each group's wildcard patterns
# against the whole taxonomy string and applies the database's group set
# operations. (scripts/faprotax_collapse.py, used before, matched substrings,
# skipped the set operations and read three database lines as group names.)
faprotax_dir <- file.path(REPO_ROOT, "data/microbiome/ref_db/FAPROTAX_1.2.12")
faprotax_db  <- file.path(faprotax_dir, "FAPROTAX.txt")
collapse_py  <- file.path(faprotax_dir, "collapse_table.py")
python_bin   <- Sys.getenv("PYTHON", "python3")
stopifnot(file.exists(faprotax_db), file.exists(collapse_py))

# Bacterial ASVs only: chloroplast and mitochondrial ASVs are removed, as for
# the community analyses (FAPROTAX files SILVA's mitochondria, which sit in
# Rickettsiales, under intracellular parasites)
ps16 <- readRDS(file.path(REPO_ROOT, "results/microbiome/16S/phyloseq_16S_filtered.rds"))
otu16 <- as(otu_table(ps16), "matrix")
if (!taxa_are_rows(ps16)) otu16 <- t(otu16)
tax16 <- as.data.frame(as(tax_table(ps16), "matrix"), stringsAsFactors = FALSE)
tax16 <- tax16[rownames(otu16), ]
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_strings <- apply(tax16[, ranks], 1, function(p) {
  p[is.na(p)] <- ""
  paste(p, collapse = ";")
})

work_dir <- tempfile("faprotax_")
dir.create(work_dir)
faprotax_infile <- file.path(work_dir, "faprotax_input_16S.tsv")
write.table(data.frame(taxonomy = tax_strings, otu16, check.names = FALSE),
            faprotax_infile, sep = "\t", row.names = FALSE, quote = FALSE)
log_msg("  FAPROTAX input: ", nrow(otu16), " ASVs, ", ncol(otu16), " samples")

collapsed_file <- file.path(work_dir, "functions.tsv")
g2r_file <- file.path(work_dir, "groups2records.tsv")
cmd <- sprintf(paste("'%s' '%s' -i '%s' -g '%s' -d taxonomy",
                     "--column_names_are_in first_data_line -n none",
                     "-o '%s' --out_groups2records_table '%s' -v 2>&1"),
               python_bin, collapse_py, faprotax_infile, faprotax_db, collapsed_file, g2r_file)
log_msg("  Running FAPROTAX collapse_table.py...")
faprotax_log <- system(cmd, intern = TRUE)
cat(grep("Assigned|leftovers", faprotax_log, value = TRUE), sep = "\n")
stopifnot(file.exists(collapsed_file), file.exists(g2r_file))

# Function x sample read counts
func <- read.delim(collapsed_file, check.names = FALSE, stringsAsFactors = FALSE)
colnames(func)[1] <- "function"
# Group x ASV membership (rows are records, i.e. ASV taxonomy strings, in
# input order)
g2r <- read.delim(g2r_file, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
member <- as.matrix(g2r[, -1]) != 0
stopifnot(identical(g2r[[1]], unname(tax_strings)))
n_taxa <- colSums(member)[func[["function"]]]
asv_reads <- rowSums(otu16)
representative <- sapply(func[["function"]], function(f) {
  idx <- which(member[, f])
  idx <- idx[order(asv_reads[idx], decreasing = TRUE)][seq_len(min(3, length(idx)))]
  paste(unique(tax_strings[idx]), collapse = "; ")
})
report <- data.frame(
  `function` = func[["function"]],
  n_taxa = as.integer(n_taxa),
  total_reads = rowSums(func[, -1, drop = FALSE]),
  representative_taxa = representative,
  check.names = FALSE
)
keep <- report$n_taxa > 0
report <- report[keep, ]
report <- report[order(report$total_reads, decreasing = TRUE), ]
write.table(func[keep, ], file.path(FUNC_DIR, "faprotax_functions_16S.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(report, file.path(FUNC_DIR, "faprotax_report_16S.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
log_msg("  FAPROTAX: ", nrow(report), " functions with at least one ASV; ",
        sum(rowSums(member) > 0), " of ", nrow(member), " ASVs assigned")
for (i in seq_len(min(10, nrow(report)))) {
  log_msg("    ", report[["function"]][i], ": ", report$total_reads[i], " reads, ",
          report$n_taxa[i], " ASVs")
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
