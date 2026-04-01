# Supplementary Materials for ARES-SENTINEL-GBM

## Overview

This document describes the supplementary materials accompanying the manuscript:

**"ARES-SENTINEL-GBM: An Open-Source Pipeline for Drug Target Discovery in Glioblastoma with Comprehensive Multi-Omics Validation"**

---

## Supplementary Files

### Supplementary Table S1: Complete List of Differential Metabolic Reactions

**File:** `Supplementary_Table_S1.xlsx`

**Description:** Complete list of 847 metabolic reactions with differential flux between GBM and normal astrocytes.

**Columns:**
- `reaction_id`: Reaction identifier (Human-GEM nomenclature)
- `reaction_name`: Human-readable reaction name
- `flux_gbm`: Mean flux in GBM condition
- `flux_astro`: Mean flux in astrocyte condition
- `delta`: Flux difference (GBM - Astro)
- `abs_delta`: Absolute flux difference
- `p_value`: Statistical significance (t-test)
- `fdr`: False discovery rate (Benjamini-Hochberg)
- `pathway`: Associated metabolic pathway
- `gene_association`: Associated genes

---

### Supplementary Table S2: Gene Knockout Results

**File:** `Supplementary_Table_S2.xlsx`

**Description:** Complete results of in silico gene knockout simulations for 200 metabolic genes.

**Columns:**
- `gene_symbol`: HGNC gene symbol
- `gene_ensg`: Ensembl gene ID
- `knockout_type`: Single or combinatorial
- `growth_ratio`: Growth rate (KO / WT)
- `effect_classification`: LETHAL / STRONG / WEAK
- `essentiality_score`: Continuous essentiality score (0-1)
- `literature_support`: Number of PubMed citations
- `druggability_score`: Druggability assessment (0-1)

---

### Supplementary Table S3: Survival Analysis Results

**File:** `Supplementary_Table_S3.xlsx`

**Description:** Complete survival analysis for all metabolic genes in TCGA-GBM cohort.

**Columns:**
- `gene_symbol`: HGNC gene symbol
- `hazard_ratio`: Cox PH hazard ratio
- `hr_ci_lower`: 95% CI lower bound
- `hr_ci_upper`: 95% CI upper bound
- `cox_p_value`: Cox model p-value
- `km_p_value`: Log-rank test p-value
- `median_survival_high`: Median survival (high expression group), months
- `median_survival_low`: Median survival (low expression group), months
- `n_high`: Number of patients in high expression group
- `n_low`: Number of patients in low expression group

---

### Supplementary Table S4: Literature Curation of Top Targets

**File:** `Supplementary_Table_S4.xlsx`

**Description:** Manual curation of literature evidence for top 50 metabolic targets.

**Columns:**
- `gene_symbol`: HGNC gene symbol
- `pubmed_count`: Number of relevant PubMed citations
- `tcga_evidence`: Evidence from TCGA analysis (YES/NO)
- `drug_target`: Known drug target (YES/NO)
- `essentiality_data`: CRISPR/RNAi essentiality data (YES/NO)
- `gbm_specific`: GBM-specific evidence (YES/NO)
- `key_references`: Key reference PMIDs
- `clinical_trials`: Associated clinical trials (if any)
- `notes`: Additional comments

---

### Supplementary Figure S1: Pipeline Workflow Diagram

**File:** `Supplementary_Figure_S1.pdf`

**Description:** Complete workflow diagram showing all steps of the ARES-GBM pipeline from input data to final results.

**Panel A:** Data preprocessing and quality control
**Panel B:** Context-specific model reconstruction using iMAT
**Panel C:** Flux balance analysis and differential analysis
**Panel D:** Downstream analyses (knockout, survival, validation)

---

### Supplementary Figure S2: Additional Benchmark Results

**File:** `Supplementary_Figure_S2.pdf`

**Description:** Extended benchmark comparison with additional metrics.

**Panel A:** Precision-recall curves for all methods
**Panel B:** Runtime comparison across different model sizes
**Panel C:** Memory usage comparison
**Panel D:** Scalability analysis (genes vs. runtime)

---

### Supplementary Figure S3: Pathway Enrichment - Full Results

**File:** `Supplementary_Figure_S3.pdf`

**Description:** Complete pathway enrichment analysis for all Reactome pathways.

**Panel A:** Volcano plot of pathway enrichment
**Panel B:** Network visualization of enriched pathways
**Panel C:** Heatmap of pathway activity across samples

---

### Supplementary Figure S4: Additional Survival Curves

**File:** `Supplementary_Figure_S4.pdf`

**Description:** Kaplan-Meier survival curves for all top 20 metabolic targets.

Each panel shows:
- Survival curves for high vs. low expression
- Number at risk table
- Hazard ratio with 95% CI
- Log-rank p-value

---

### Supplementary Figure S5: Network Analysis - Additional Views

**File:** `Supplementary_Figure_S5.pdf`

**Description:** Additional network topology analyses.

**Panel A:** Degree distribution with power-law fit
**Panel B:** Betweenness centrality distribution
**Panel C:** Clustering coefficient distribution
**Panel D:** Network motif analysis

---

### Supplementary Figure S6: Drug-Target Interaction Network

**File:** `Supplementary_Figure_S6.pdf`

**Description:** Network visualization of identified targets with known/potential drugs.

**Nodes:** Metabolic targets (circles) and drugs (squares)
**Edges:** Known interactions (solid) and predicted interactions (dashed)
**Colors:** Drug approval status (FDA-approved, clinical trial, preclinical)

---

### Supplementary Methods

**File:** `Supplementary_Methods.pdf`

**Description:** Detailed methodological information.

**Sections:**
1. Data Preprocessing
   - Expression data normalization
   - Quality control procedures
   - Batch effect correction

2. Genome-Scale Metabolic Model
   - Human-GEM version and modifications
   - Reaction-gene mapping
   - Constraint definition

3. iMAT Algorithm Implementation
   - Mathematical formulation
   - Solver settings
   - Convergence criteria

4. Flux Balance Analysis
   - Objective function definition
   - Constraint handling
   - Solution validation

5. Statistical Analysis
   - Differential expression testing
   - Multiple testing correction
   - Cross-validation procedure
   - Bootstrap methodology

6. Literature Curation Protocol
   - Search strategy
   - Inclusion/exclusion criteria
   - Evidence scoring system

7. Software Implementation
   - Programming languages and versions
   - Dependencies
   - Performance optimization

---

### Supplementary Data S1: Example Input Files

**File:** `Supplementary_Data_S1.zip`

**Description:** Example input files for running ARES-GBM.

**Contents:**
- `expression_matrix.csv`: Example gene expression matrix
- `survival_data.csv`: Example clinical survival data
- `config.yaml`: Example configuration file
- `run_analysis.sh`: Example batch script

---

### Supplementary Data S2: Complete Results Archive

**File:** `Supplementary_Data_S2.zip`

**Description:** Complete results from all analyses.

**Contents:**
- All flux analysis results (CSV)
- All knockout simulation results (CSV)
- All survival analysis results (CSV)
- All figure source data (CSV)
- Model files (SBML)

---

## File Formats

| Format | Description |
|--------|-------------|
| `.xlsx` | Microsoft Excel (Office Open XML) |
| `.pdf` | Portable Document Format (high resolution) |
| `.csv` | Comma-separated values (UTF-8) |
| `.zip` | Compressed archive |
| `.sbml` | Systems Biology Markup Language |

---

## Data Availability

All supplementary files are available at:
- **GitHub:** https://github.com/ARES-Sentinel-GBM/-gbm-research
- **Zenodo:** [DOI will be assigned upon publication]

---

## Contact

For questions about supplementary materials:
- **Email:** admin@ares-bio.com
- **GitHub Issues:** https://github.com/ARES-Sentinel-GBM/-gbm-research/issues
