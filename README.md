# Integrated Transcriptome-Microbiome Analysis of ISS VEG-05 Tomatoes

![Space_Tomato_Multi-Omics_Research](Space_Tomato_Multi-Omics_Research.png)

🌐 **Live site:** https://dr-richard-barker.github.io/veg05-integrated-omics/

## Overview

Complete analysis pipeline for an integrated multi-omics study of dwarf tomato (*Solanum lycopersicum* cv. Red Robin) grown aboard the International Space Station (ISS) during the VEG-05 experiment. The experiment compared two lighting regimes (red-rich and blue-rich) with matched ground controls at Kennedy Space Center.

## Data Sources

| Dataset | OSDR Accession | Description |
|---------|---------------|-------------|
| RNA-seq | OSD-767 | Host transcriptome (Leaf, Adv-Root) |
| Microbiome | OSD-766 | 16S rRNA + ITS amplicon sequencing |

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 0 | `scripts/00_download_data.py` | Download processed RNA-seq data and microbiome FASTQ from the OSDR API |
| 1 | `scripts/01_parse_metadata.py` | Parse ISA-Tab metadata; build the RNA-seq ↔ microbiome sample crosswalk |
| 2 | `scripts/02_microbiome_dada2.R` | DADA2 amplicon processing (16S + ITS) |
| 3 | `scripts/03_rnaseq_deseq2.R` | DESeq2 differential expression |
| 4 | `scripts/04_community_health.R` | Alpha/beta diversity, PERMANOVA, dysbiosis index, ALDEx2 |
| 5 | `scripts/05_functional_prediction.R` | FAPROTAX 1.2.12 (its own `collapse_table.py`) and FUNGuild |
| 6 | `scripts/06_wgcna.R` | WGCNA co-expression networks |
| 7 | `scripts/07_integration_mofa.R` | MOFA+ multi-omics integration (leaf) |
| 8 | `scripts/08_correlation_networks.R` | Module–taxon correlation networks |
| 9 | `scripts/09_directionality.R` | GO enrichment of modules; directionality inference |
| 10 | `scripts/10_figures_main.R` | Main figures 2–7 |
| 11 | `scripts/11_figures_supplementary.R` | Supplementary figures S1–S7 |
| 12 | `scripts/12_export_tables.R` | Supplementary tables S1–S10 |
| — | `build_fig1_schematic.py` | Figure 1 (study-design schematic) |
| — | `build_module_annotations_and_fig4.py`, `fix_fig4_wgcna.py`, `fix_fig6_module_taxon.py` | Module GO labels (`results/module_annotations.tsv`); redrawn Figures 4 and 6 |
| — | `build_manuscript_figures.py` | Builds `manuscript_with_figures.md` and `.docx` from `manuscript_draft.md` |
| — | `build_slides.py` | Builds the slide deck `ISS_Tomato_Multi-Omics.pptx` from the figures |
| — | `build_infographic.py` | Builds the README banner `Space_Tomato_Multi-Omics_Research.png` |

The R and Python pipeline scripts were run in an analysis sandbox and hard-code
`REPO_ROOT <- "/mnt/shared-workspace/veg05-integrated-omics"`. They read inputs from
`data/` and write to `results/` under that root; neither directory is tracked here.
The outputs used in the manuscript were copied into the top-level folders listed
below, so the scripts are not runnable from this repository without editing those paths.

## Key Results

### Microbiome
- **16S**: 348 bacterial ASVs (after chloroplast/mitochondria removal), 142 samples, 3,558,814 reads
- **ITS**: 78 fungal ASVs, 141 samples (121 with reads after processing), 156,102 reads
- **PERMANOVA (16S leaf)**: flight R²=0.28, F=8.24, p=0.001; light n.s.; flight × light R²=0.11, p=0.035 (terms added sequentially)
- **Dysbiosis (16S leaf)**: Flight 2.06±0.99 vs Ground 1.00±0.47
- **Dysbiosis (ITS leaf)**: Flight 3.69±2.78 vs Ground 1.00±0.22

### Transcriptome
- **Leaf**: 19,844 genes, 21 samples. Flight×Light interaction: 3,189 DEGs
- **Adv-Root**: 20,875 genes, 15 samples. Flight main effect (averaged over light): 1,501 DEGs (83% upregulated)
- **Key finding**: Spaceflight response under blue light (4,716 DEGs) >> red light (523 DEGs) in leaf

### Integration
- **WGCNA**: Leaf 12 modules, Adv-Root 16 modules (plus unassigned genes)
- **Key module**: Adv-Root black (169 genes) — microbe-driven (r=−0.85, padj=0.009 with ITS dysbiosis, 12 plants; r=−0.38, padj=0.32 with flight)
- **MOFA+**: 5 factors, 18 leaf samples with all three views from the same leaf. Factor 1 explains 44% of transcriptome variance and correlates with flight (ρ=−0.68, padj=0.037) and light (ρ=+0.63, padj=0.053), which are partly confounded in these samples
- **Module-taxon correlations**: 31 significant (BH padj<0.05; 28 leaf 16S, 3 root ITS), linking flight modules to Methylobacterium, Burkholderia, Azospirillum
- **GO enrichment**: 80 module–term associations (73 unique terms, Leaf), 96 (86 unique, Adv-Root). Photosynthesis, oxidative stress, ribosome

### Functional Prediction
- **FAPROTAX** (1.2.12, hyphens as word separators): 205 of 348 bacterial ASVs in 44 groups; largest by reads are chemoheterotrophy, aerobic chemoheterotrophy, ureolysis, nitrogen fixation, methylotrophy and methanol oxidation (the last two almost entirely *Methylobacterium-Methylorubrum*)
- **FUNGuild**: Saprotrophs (11 ASVs), Plant Pathogens (2 ASVs: Fusarium spp.)

## Repository Structure

```
veg05-integrated-omics/
├── scripts/                     # Analysis pipeline, steps 00–12 (see Pipeline table)
├── microbiome/                  # DADA2 taxonomy and read tracking (16S, ITS); RERUN_2026-10-15.md
│   ├── community_health/        # Alpha/beta diversity, PERMANOVA, dysbiosis, ALDEx2
│   └── functional/              # FAPROTAX and FUNGuild outputs
├── rnaseq/                      # DESeq2 DEG tables and summary, sample-distance QC; RERUN_2026-10-14.md
│   ├── go_enrichment/           # GO enrichment of WGCNA modules
│   ├── wgcna_Leaf/              # Module–trait correlations, hub genes
│   └── wgcna_Adv-Root/
├── integration/                 # MOFA+ model, factor values, variance, weights
│   └── networks/                # Module–taxon correlations, directionality
├── results/                     # module_annotations.tsv (module GO labels)
├── figures/                     # Figures 1–7 and S1–S7 (PDF, PNG, SVG)
├── supplementary_tables/        # Tables S1–S10 (CSV)
├── manuscript/latex/            # LaTeX manuscript (main.tex, references.bib) + figure copies
├── manuscript_draft.md          # Markdown manuscript (source for the two files below)
├── manuscript_with_figures.md   # Draft with figures embedded
├── manuscript_with_figures.docx
├── build_*.py, fix_*.py         # Figure and manuscript builders (see Pipeline table)
├── docs/                        # GitHub Pages site and compiled manuscript PDF
├── execution_trace/             # Original analysis plan and agent notebook
├── COMPARISON_with_SA_integration.md
├── ISS_Tomato_Multi-Omics.pptx   # summary slides, built by build_slides.py
└── environment.yml, requirements.txt, renv.lock, sessionInfo.txt
```

## Software Versions

- Original run: R 4.4.2, Bioconductor 3.20; DADA2 1.34.0, phyloseq 1.50.0, DESeq2 1.46.0, WGCNA 1.74, MOFA2 1.16.0, ALDEx2 1.38.0, clusterProfiler 4.14.0, vegan 2.7-5; Python 3.12, mofapy2 0.7.4, cutadapt
- Re-runs of October 2026 (`rnaseq/RERUN_2026-10-14.md`, `microbiome/RERUN_2026-10-15.md`): R 4.6.1; DESeq2 1.52.0, WGCNA 1.74, DADA2 1.40.0, phyloseq 1.56.0, vegan 2.7-6, ALDEx2 1.44.0, MOFA2 1.22.1; Python mofapy2 0.7.5, cutadapt 5.2
- SILVA 138.1, UNITE v7.2 QIIME release (1 Dec 2017, 99%), FAPROTAX 1.2.12
- `renv.lock` pins the R packages of the re-runs (243 packages with their dependencies; restore with `renv::restore()`); `requirements.txt` pins the Python packages

## Limitations

1. Small sample sizes (n=21 leaf, n=15 root for RNA-seq)
2. 16S leaf samples have very low bacterial reads after host DNA removal (median 89)
3. ITS dataset has uneven group representation (some flight groups n=1-2)
4. UNITE v7.2 (2017) used instead of a current release; the original run's 10 Oct 2017 build is no longer distributed, so the re-run uses the 1 Dec 2017 build
5. Two flight adventitious-root samples retained a single ITS read each
6. FUNGuild used manual genus-level assignment (script database incompatible)
7. Study is correlational; causal relationships require experimental validation

## License

MIT License. See LICENSE file for details.

## Citation

If using this code or data, please cite:
- NASA OSDR: OSD-766, OSD-767
- This repository: DOI via Zenodo — badge pre-wired to repo id 1287757020; activates on Zenodo link + release:
  [![DOI](https://zenodo.org/badge/1287757020.svg)](https://zenodo.org/badge/latestdoi/1287757020)
