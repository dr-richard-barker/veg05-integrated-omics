# Integrated Transcriptome-Microbiome Analysis of ISS VEG-05 Tomatoes

## Plan for GitHub/Zenodo-Ready Repository + Full Analysis Run
**Target journal:** npj Microgravity
**Data sources:** OSD-767 (transcriptome, DOI: 10.26030/rrtk-h481) + OSD-766 (microbiome, DOI: 10.26030/ywjn-5e23)

---

## 1. Project Summary

The VEG-05 experiment grew dwarf tomato (*Solanum lycopersicum* cv. Red Robin) aboard the ISS (Dec 2022–Mar 2023, 101 days) under two lighting regimes: red-rich (90%R:10%B) and blue-rich (50%R:50%B), with matched ground controls at Kennedy Space Center. Two companion datasets capture the plant transcriptome (OSD-767: 36 RNA-seq samples, leaf + adventitious root) and the associated microbial communities (OSD-766: 283 amplicon samples, 16S + ITS across leaf, root, AdvRoot, wick, soil, water, swabs, fruit).

This project integrates both datasets to:
1. Characterize **community health** of the plant-associated microbiome under each treatment (Flight×Red, Flight×Blue, Ground×Red, Ground×Blue) using three complementary perspectives: diversity/composition, dysbiosis relative to ground baseline, and predicted functional capacity.
2. Identify **coordinated changes** between plant gene expression and microbial abundance via unsupervised factor analysis (MOFA+) and WGCNA module-taxon correlations.
3. Infer **directionality** of host-microbe interactions — distinguishing microbes shifting in response to plant gene expression changes (host-driven) from those potentially driving plant transcriptional responses (microbe-driven) — using gene-function annotation of correlated modules.

---

## 2. Data Inventory (Phase 1 findings)

### OSD-767 (Transcriptome)
- **36 samples**: 18 Flight + 18 Ground; 21 Leaf + 15 AdvRoot; 19 Red + 17 Blue
- **Plants**: SN01–SN12 (14 plants, unbalanced across treatments)
- **Processed data available**: RSEM unnormalized gene counts (GLDS-709_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv, 4.47 MB), ISA metadata, runsheet, QC metrics
- **Reference genome**: Ensembl Plants release 63, *S. lycopersicum* SL4.0 (GCA_000188115v5)
- **Processing pipeline**: GeneLab GLbulkRNAseq v2 (STAR 2.7.11b + RSEM 1.3.3, Nextflow 25.10.3)
- **Sample size per group (Leaf)**: Flt-Red=5, Flt-Blue=4, Gnd-Red=6, Gnd-Blue=6 (n=21, adequate)
- **Sample size per group (AdvRoot)**: Flt-Red=5, Flt-Blue=4, Gnd-Red=3, Gnd-Blue=2 (n=14, limited — Ground-Blue has only 2 replicates)

### OSD-766 (Microbiome)
- **283 samples** (142 16S + 141 ITS) across compartments:
  - **Integration-relevant**: leaf (42), AdvRoot (30) — match RNA-seq tissues
  - **Rhizosphere context**: root (20), wick (16), soil (47, 3 replicates per plant)
  - **Additional**: water (8), surface swabs (24), fruit (36)
- **Only raw FASTQ available** (1.32 GB total: 1.0 GB 16S + 0.35 GB ITS) — no ASV/BIOM/taxonomy tables. Must process from scratch with DADA2.
- **Plants**: SN01–SN12 + water controls SN003/SN006
- **Treatment structure**: Flight vs Ground × Red-rich vs Blue-rich (same as OSD-767)

### Sample matching for integration
RNA-seq and microbiome samples share the same plants (SN01–SN12), treatments (Flight/Ground × Red/Blue), and tissues (Leaf, AdvRoot). Matching key: **Flight/Ground × Plant × Tissue × Light**. This enables per-sample integration, not just group-level comparison.

---

## 3. Repository Structure (GitHub + Zenodo)

```
veg05-integrated-omics/
├── README.md                          # Project goals, data sources, pipeline overview, reproduction instructions
├── LICENSE                            # MIT (code), CC-BY-4.0 (data/results)
├── CITATION.cff                       # Citation metadata for Zenodo
├── .gitignore                         # Exclude raw data, large intermediates
├── requirements.txt                   # Python dependencies
├── environment.yml                    # Conda environment (R + Python)
├── config/
│   └── config.yaml                    # Parameters, thresholds, file paths
├── data/
│   ├── README.md                      # Data sources, DOIs, download instructions
│   ├── metadata/
│   │   ├── sample_metadata_rnaseq.csv     # Parsed from OSD-767 ISA
│   │   ├── sample_metadata_microbiome.csv # Parsed from OSD-766 ISA + FASTQ names
│   │   └── sample_crosswalk.csv           # RNA-seq ↔ microbiome sample matching
│   ├── rnaseq/
│   │   └── RSEM_Unnormalized_Counts.csv   # Downloaded from OSDR (4.47 MB)
│   └── microbiome/
│       └── README.md                  # FASTQ not in repo — download from OSDR
├── scripts/
│   ├── 00_download_data.py            # Download processed RNA-seq + microbiome FASTQ from OSDR API
│   ├── 01_parse_metadata.py           # Parse ISA metadata, build sample crosswalk
│   ├── 02_microbiome_dada2.R          # DADA2: filter, learn errors, denoise, merge, chimera removal, taxonomy (SILVA 138.2 + UNITE 10.0)
│   ├── 03_rnaseq_deseq2.R             # DESeq2: DE analysis per tissue, contrasts, LFC shrinkage, VST
│   ├── 04_community_health.R          # Alpha diversity, beta diversity, dispersion, dysbiosis index, differential abundance (ANCOM-BC2)
│   ├── 05_functional_prediction.R     # FAPROTAX (16S functional categories) + FUNGuild (ITS fungal guilds)
│   ├── 06_wgcna.R                     # WGCNA gene modules, module-trait correlations, module eigengenes
│   ├── 07_integration_mofa.R          # MOFA+ unsupervised integration (transcriptome + microbiome views)
│   ├── 08_correlation_networks.R      # Bipartite correlation networks: WGCNA modules ↔ taxa genera
│   ├── 09_directionality.R            # GO enrichment of correlated modules → host-driven vs microbe-driven classification
│   ├── 10_figures_main.R              # Main manuscript figures (npj microgravity)
│   ├── 11_figures_supplementary.R     # Supplementary figures
│   └── 12_export_tables.R             # Supplementary tables
├── results/
│   ├── microbiome/
│   │   ├── asv_table_16S.rds          # ASV count table (16S)
│   │   ├── asv_table_ITS.rds          # ASV count table (ITS)
│   │   ├── taxonomy_16S.csv           # SILVA taxonomy
│   │   ├── taxonomy_ITS.csv           # UNITE taxonomy
│   │   ├── phyloseq_16S.rds           # phyloseq object
│   │   ├── phyloseq_ITS.rds           # phyloseq object
│   │   ├── alpha_diversity.csv        # Per-sample diversity metrics
│   │   ├── beta_diversity.rds         # Distance matrices + PERMANOVA results
│   │   ├── dysbiosis_index.csv        # Per-sample dysbiosis scores
│   │   ├── differential_abundance.csv # ANCOM-BC2 results
│   │   ├── faprotax_functional.csv    # 16S predicted functions
│   │   └── funguild_guilds.csv        # ITS fungal guild assignments
│   ├── rnaseq/
│   │   ├── deseq2_dds.rds             # DESeqDataSet object
│   │   ├── deseq2_results_leaf.csv    # DE results (leaf tissue)
│   │   ├── deseq2_results_advroot.csv # DE results (advroot tissue)
│   │   ├── vst_transformed.csv        # VST normalized expression
│   │   └── qc_metrics.csv             # QC summary
│   ├── integration/
│   │   ├── mofa_model.rds             # MOFA+ model
│   │   ├── mofa_factor_values.csv     # Sample factor scores
│   │   ├── mofa_weights.csv           # Feature weights per factor
│   │   ├── wgcna_modules.csv          # Gene-module assignments
│   │   ├── wgcna_eigengenes.csv       # Module eigengenes per sample
│   │   ├── module_taxa_correlations.csv  # Module-taxon correlation matrix
│   │   ├── directionality_classification.csv  # Host-driven vs microbe-driven
│   │   └── bipartite_network.graphml  # Network for Cytoscape
│   ├── figures/
│   │   ├── main/                      # Main text figures (SVG + PNG)
│   │   └── supplementary/             # Supplementary figures
│   └── tables/
│       └── supplementary/             # Supplementary tables (CSV)
├── manuscript/
│   ├── npj_microgravity_draft.md      # Manuscript draft
│   ├── figure_captions.md
│   └── supplementary_notes.md
└── docs/
    └── analysis_plan.md               # This document (archived)
```

**Zenodo deposit contents**: All code, metadata, processed results (ASV tables, DE results, integration outputs), figures, tables. Raw FASTQ excluded (available via NASA OSDR with DOIs). README links to OSDR for raw data.

---

## 4. Analysis Pipeline

### Step 0: Data acquisition & metadata parsing
- Download RSEM counts, ISA metadata, runsheet, QC metrics from OSD-767 (~5 MB)
- Download all microbiome FASTQ from OSD-766 (~1.32 GB) via OSDR API
- Parse ISA metadata for both studies into standardized CSVs
- Build sample crosswalk: match RNA-seq ↔ microbiome by Flight/Ground × Plant × Tissue × Light
- **Output**: `sample_metadata_rnaseq.csv`, `sample_metadata_microbiome.csv`, `sample_crosswalk.csv`

### Step 1: Microbiome processing (DADA2)
- **Primer trimming**: cutadapt to remove 16S (515F/806R) and ITS (ITS1f/ITS2) primers
- **Quality filtering**: Inspect read quality profiles from a subset of samples; set truncation parameters based on quality plots (typical: truncLen for 16S ~240/160, ITS variable)
- **DADA2 pipeline** (separate for 16S and ITS):
  - Filter and trim (`filterAndTrim`)
  - Learn error rates (`learnErrors`, with multithreading)
  - Denoise (`dada`)
  - Merge paired reads (`mergePairs`) — ITS may need higher mismatch tolerance
  - Make sequence table (`makeSequenceTable`)
  - Remove chimeras (`removeBimeraDenovo`)
  - Assign taxonomy: SILVA 138.2 for 16S, UNITE 10.0 for ITS
- **ASV table construction**: Build count matrices, remove chloroplast/mitochondria contaminants (16S), assign ASV IDs
- **Phyloseq objects**: Combine ASV tables + taxonomy + sample metadata
- **Output**: ASV tables, taxonomy CSVs, phyloseq RDS objects

### Step 2: RNA-seq differential expression (DESeq2)
- Load RSEM unnormalized counts; filter low-count genes (≥10 counts in ≥4 samples)
- **Design**: Analyze each tissue separately (Leaf: n=21, AdvRoot: n=14)
  - Leaf: `~ flight * light` (4 groups: Flt-Red, Flt-Blue, Gnd-Red, Gnd-Blue, n=4–6 per group)
  - AdvRoot: `~ flight * light` (note: Gnd-Blue has only n=2 — flag as exploratory)
- **Contrasts** (per tissue):
  1. Flight vs Ground (main spaceflight effect, pooled across light)
  2. Red vs Blue (main light effect, pooled across flight)
  3. Flight:Red vs Ground:Red (spaceflight effect under red light)
  4. Flight:Blue vs Ground:Blue (spaceflight effect under blue light)
  5. Flight×Light interaction
- **LFC shrinkage**: apeglm for visualization/ranking
- **Transformation**: VST for downstream WGCNA and MOFA+
- **QC**: PCA, dispersion plots, sample distance heatmap
- **Output**: DE result tables per tissue, VST matrix, DESeqDataSet RDS

### Step 3: Community health metrics
**3a. Diversity + composition:**
- **Alpha diversity**: Observed ASVs, Shannon, Simpson, Faith's phylogenetic diversity (16S only — requires phylogenetic tree from DADA2)
- **Beta diversity**: Bray-Curtis, Jaccard, weighted/unweighted UniFrac (16S); Bray-Curtis, Jaccard (ITS)
- **PERMANOVA**: Test flight, light, compartment, and interactions (adonis2, 999 permutations)
- **Dispersion**: betadisper + permutest (homogeneity of multivariate dispersions)
- **Composition**: Relative abundance stacked bar plots at phylum/genus level; differential abundance via ANCOM-BC2 (flight vs ground, within each compartment × light)

**3b. Dysbiosis index:**
- For each flight sample: calculate Bray-Curtis distance to the centroid of its matched ground-control group (same compartment × light)
- Normalize by the average within-group ground distance
- Dysbiosis index = (distance to ground centroid) / (mean ground within-group distance)
- Values >1 = greater deviation than expected within-group variation
- Compute per compartment (leaf, AdvRoot, root, wick, soil) and per treatment

**3c. Functional prediction:**
- **16S → FAPROTAX**: Predict functional categories (aerobic chemoheterotrophy, fermentation, nitrogen fixation, nitrate reduction, plant pathogenesis, etc.). FAPROTAX is fast, well-established for environmental/plant microbiomes, and avoids the heavy PICRUSt2 installation/runtime overhead.
- **ITS → FUNGuild**: Assign fungal functional guilds (plant pathogens, endophytes, arbuscular mycorrhizal, saprotrophs, wood saprotrophs, etc.) with confidence scores
- Compare functional profiles across treatments (Flight vs Ground, Red vs Blue)
- **Output**: Functional profile tables, guild assignments, comparison statistics

### Step 4: WGCNA (gene co-expression modules)
- Input: VST-transformed expression matrix (top 5000 most variable genes)
- **Separate by tissue**: Run WGCNA on Leaf (n=21) and AdvRoot (n=14) separately
  - Leaf (n=21): adequate for WGCNA (≥15 minimum, 20+ recommended)
  - AdvRoot (n=14): below recommended minimum — run but flag results as exploratory
- Soft power selection (scale-free topology fit R² ≥ 0.8)
- Module detection (minModuleSize=30, mergeCutHeight=0.25)
- Module eigengene calculation
- Module-trait correlations: correlate eigengenes with flight, light, and interaction
- **Output**: Module assignments, eigengenes, module-trait correlations, hub genes

### Step 5: MOFA+ integration (unsupervised)
- **Views**: 
  - View 1: Transcriptome (VST expression, top 2000 most variable genes, Leaf samples)
  - View 2: Microbiome (genus-level CLR-transformed abundances, 16S, Leaf samples)
  - View 3: Microbiome (genus-level CLR-transformed abundances, ITS, Leaf samples)
- **Samples**: Matched RNA-seq + microbiome leaf samples (via crosswalk)
- **Preprocessing**: Filter low-variance features, CLR transform microbiome, center/scale
- **MOFA+ training**: 10–15 factors, Gaussian likelihood for all views
- **Factor interpretation**: 
  - Correlate factor scores with treatment variables (flight, light)
  - Identify top-weighted features per factor per view
  - Variance decomposition (R² per factor per view)
- **Output**: MOFA model, factor scores, weights, variance explained

### Step 6: Correlation networks + directionality inference
**6a. Bipartite correlation networks:**
- Correlate WGCNA module eigengenes with genus-level taxa abundances (CLR-transformed) for matched samples
- Use Spearman correlation with BH-FDR correction (padj < 0.05, |rho| > 0.5)
- Build bipartite network: modules ↔ taxa genera
- Export as GraphML for Cytoscape visualization

**6b. Directionality inference (gene-function based):**
- For each module significantly correlated with a taxon:
  - Perform GO enrichment (Biological Process) on module genes using *S. lycopersicum* GO annotations (org.Slycersicum.db or InterPro2GO mapping)
  - Classify module function:
    - **Host-driven** (plant gene changes → microbial shift): modules enriched for defense response (GO:0006952), immune system process (GO:0002376), secondary metabolism (GO:0019748), hormone signaling (JA/SA/ethylene), cell wall modification
    - **Microbe-driven** (microbial shift → plant gene response): modules enriched for nutrient transport (GO:0006810), ion homeostasis (GO:0050801), oxidative stress response (GO:0006979), water deprivation response (GO:0009414)
    - **Ambiguous**: modules with mixed or unclear enrichment
  - Assign directionality label to each significant module-taxon edge
- **Output**: Directionality classification table, annotated bipartite network

### Step 7: Figures and tables
**Main figures (npj microgravity format):**
1. **Study design overview**: Schematic of VEG-05 experiment, sample types, treatments, integration approach
2. **Transcriptomic responses**: PCA (flight/light/tissue), volcano plots (Flight vs Ground per tissue), heatmap of top DE genes
3. **Microbiome community health**: Alpha diversity boxplots (4 treatments × compartments), beta diversity PCoA (Bray-Curtis), dysbiosis index heatmap
4. **Microbiome composition**: Stacked bar plots (top genera per treatment × compartment), differential abundance dot plot
5. **Functional profiles**: FAPROTAX functional categories (16S), FUNGuild guilds (ITS), comparison across treatments
6. **Integrated analysis**: MOFA+ variance decomposition, factor scores colored by treatment, top feature weights
7. **Host-microbe correlation network**: Bipartite network (modules ↔ taxa), colored by directionality (host-driven vs microbe-driven)
8. **Directionality summary**: Sankey or alluvial diagram showing flow from gene modules → functions → taxa → direction

**Supplementary:**
- S1: Full sample metadata table
- S2: DADA2 quality profiles and filtering statistics
- S3: Full DESeq2 results (all contrasts)
- S4: Alpha diversity all metrics
- S5: PERMANOVA full results
- S6: ANCOM-BC2 differential abundance full results
- S7: WGCNA module assignments and hub genes
- S8: MOFA+ full factor weights
- S9: Full module-taxon correlation matrix
- S10: GO enrichment results per module

---

## 5. Key Methodological Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| RNA-seq starting point | GeneLab RSEM counts | Already processed with validated GeneLab pipeline (STAR+RSEM, SL4.0 genome) |
| RNA-seq DE design | `~ flight * light` per tissue | Captures main effects + interaction; tissue analyzed separately due to distinct biology |
| AdvRoot Ground-Blue | Flag as exploratory (n=2) | Too few replicates for robust inference; report but don't over-interpret |
| Microbiome processing | DADA2 from raw FASTQ | No processed ASV tables available; DADA2 is the standard for amplicon ASV inference |
| Taxonomy databases | SILVA 138.2 (16S), UNITE 10.0 (ITS) | Current standard reference databases |
| Functional prediction | FAPROTAX (16S) + FUNGuild (ITS) | Fast, reliable, well-suited for plant microbiome; PICRUSt2 noted as optional extension |
| Integration approach | MOFA+ (unsupervised) + WGCNA correlations (supervised) | Complementary: MOFA finds shared latent factors; WGCNA gives interpretable gene modules |
| Directionality | Gene-function inference via GO enrichment of correlated modules | Cross-sectional data cannot prove causation; GO-based inference is the most defensible approach |
| Compartments for integration | Leaf + AdvRoot (1:1 matched) | Only tissues with both RNA-seq and microbiome data |
| Compartments for community health | Leaf + AdvRoot + root + wick + soil | Full plant-microbiome system context |
| Multiple testing | BH-FDR (padj < 0.05) for all tests | Standard for genomics; balances discovery and rigor |
| Correlation threshold | |rho| > 0.5, padj < 0.05 | Conservative for small sample sizes |
| WGCNA gene filter | Top 5000 most variable | Standard for reducing noise while retaining signal |
| MOFA+ gene filter | Top 2000 most variable | Reduces dimensionality for factor analysis |
| Figure format | SVG (primary) + PNG (300 DPI) | Per user preference; SVG for editable manuscript figures |

---

## 6. Compute & Resource Plan

### Machine provisioning
- **DADA2 step**: Create 8-CPU / 32 GB machine via ManageMachine (DADA2 is CPU-bound for error learning and denoising)
- **All other steps**: Continue on the 8-CPU machine (DESeq2, WGCNA, MOFA+, figures all benefit from more cores)

### Time estimates (8-CPU machine)
| Step | Estimated time | Basis |
|------|---------------|-------|
| Data download (1.32 GB FASTQ + 5 MB RNA-seq) | 20–30 min | OSDR API download speed |
| Metadata parsing + crosswalk | 5 min | Script execution |
| DADA2 16S (142 samples, 8 threads) | 60–90 min | DADA2 benchmarks for ~10K reads/sample amplicon data |
| DADA2 ITS (141 samples, 8 threads) | 45–75 min | ITS reads typically shorter/fewer |
| DESeq2 (36 samples, 2 tissues) | 10 min | Standard DESeq2 runtime |
| Community health (diversity + dysbiosis + ANCOM-BC2) | 30 min | vegan + ANCOM-BC2 on 283 samples |
| FAPROTAX + FUNGuild | 15 min | Script-based, fast |
| WGCNA (2 tissues) | 20 min | WGCNA on 5K genes, 21/14 samples |
| MOFA+ (3 views, ~21 samples) | 10 min | MOFA+ typical runtime |
| Correlation networks + directionality | 20 min | Spearman + GO enrichment |
| Figures + tables | 60 min | ggplot2/ComplexHeatmap rendering |
| **Total** | **~4–5 hours** | Within 24h sandbox limit |

### Disk space
- Raw FASTQ: 1.32 GB
- DADA2 intermediates (filtered FASTQ): ~2 GB
- ASV tables + taxonomy: ~50 MB
- RNA-seq counts + DE results: ~50 MB
- Results/figures: ~200 MB
- **Total**: ~4 GB (well within disk limits)

### Memory
- DADA2: ~4–8 GB (depends on ASV count, manageable with 32 GB)
- DESeq2: ~2 GB
- WGCNA: ~4 GB
- MOFA+: ~2 GB
- **Peak**: ~8 GB (32 GB machine is sufficient)

### Package installations needed
- **R (via BiocManager)**: dada2, phyloseq, vegan, WGCNA, MOFA2, ANCOMBC, apeglm, ggprism, ggrepel, svglite
- **Python (via uv pip)**: cutadapt (primer trimming)
- **Standalone**: FAPROTAX (download script + database), FUNGuild (R/Python script)
- **Already installed**: DESeq2, clusterProfiler, ggplot2, ComplexHeatmap, biom, skbio

---

## 7. Assumptions & Limitations

1. **AdvRoot Ground-Blue (n=2)**: DESeq2 results for this group are exploratory. We will report but not draw strong conclusions.
2. **WGCNA on AdvRoot (n=14)**: Below the recommended 15-sample minimum. Results will be reported with caution.
3. **Cross-sectional data**: Directionality inference is hypothesis-generating, not proof of causation. We will explicitly state this in the manuscript.
4. **FAPROTAX vs PICRUSt2**: FAPROTAX provides functional categories (not pathway-level prediction). If reviewers request PICRUSt2, it can be run as a revision step.
5. **Sample matching**: RNA-seq and microbiome samples are from the same plants but may not be from the exact same tissue aliquots. We assume they represent the same biological condition.
6. **Contamination**: Surface swabs and water samples serve as environmental controls. We will check for cross-contamination but do not plan formal decontam analysis unless contamination is evident.
7. **Reference genome**: Using the same Ensembl Plants release 63 / SL4.0 as GeneLab for consistency.

---

## 8. Execution Order

1. Create 8-CPU machine via ManageMachine
2. Build repo scaffold (all directories, README, config, .gitignore, CITATION.cff)
3. Download data (RNA-seq processed + microbiome FASTQ)
4. Parse metadata + build crosswalk
5. Run DADA2 (16S + ITS) → ASV tables, taxonomy, phyloseq objects
6. Run DESeq2 (Leaf + AdvRoot) → DE results, VST matrix
7. Run community health (diversity, dysbiosis, ANCOM-BC2)
8. Run functional prediction (FAPROTAX + FUNGuild)
9. Run WGCNA (Leaf + AdvRoot)
10. Run MOFA+ integration
11. Run correlation networks + directionality inference
12. Generate all figures and tables
13. Write manuscript draft
14. Save everything to /mnt/results/ and organize repo
