# Data Availability Statement

## Manuscript Title
ARES-SENTINEL-GBM: An Open-Source Pipeline for Drug Target Discovery in Glioblastoma with Comprehensive Multi-Omics Validation

---

## Data Sources

### 1. Gene Expression Data

**Source:** The Cancer Genome Atlas (TCGA) - Glioblastoma Multiforme (GBM)

**Accession:** TCGA-GBM

**URL:** https://portal.gdc.cancer.gov/

**Data Type:** RNA-seq gene expression (HTSeq - TPM normalized)

**Samples:** 
- Tumor samples: n = 156
- Normal brain tissue: n = 148 (from GTEx)

**Access:** Open access (no controlled access required)

**Citation:**
> The Cancer Genome Atlas Research Network. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. Nature. 2008;455(7216):1061-1068.

---

### 2. Clinical/Survival Data

**Source:** The Cancer Genome Atlas (TCGA) - Clinical Data

**Accession:** TCGA-GBM Clinical

**URL:** https://portal.gdc.cancer.gov/

**Data Type:** Clinical annotations including:
- Overall survival time (days)
- Vital status (alive/deceased)
- Age at diagnosis
- Gender
- WHO tumor grade
- Treatment information

**Access:** Open access

---

### 3. Genome-Scale Metabolic Model

**Source:** Human-GEM

**Version:** 1.18.0

**URL:** https://github.com/SysBioChalmers/human-gem

**Repository:** SysBioChalmers/human-gem

**License:** CC-BY 4.0

**Citation:**
> Wang H, Robinson JL, Kocabas P, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. Genome Biology. 2023;24(1):1-20.

---

### 4. Pathway Annotations

**Source:** Reactome Pathway Database

**URL:** https://reactome.org/

**Version:** v84 (2024)

**License:** CC-BY 4.0

**Citation:**
> Gillespie M, Jassal B, Stephan R, et al. The reactome pathway knowledgebase. Nucleic Acids Research. 2022;50(D1):D687-D692.

---

### 5. Drug-Target Information

**Source:** DrugBank

**URL:** https://www.drugbank.ca/

**Version:** 5.1.11 (2024)

**License:** Academic use permitted with attribution

**Citation:**
> Wishart DS, Feunang YD, Guo AC, et al. DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Research. 2018;46(D1):D1074-D1082.

---

## Generated Data

### A. Processed Expression Matrix

**Description:** Normalized and batch-corrected gene expression matrix

**Format:** CSV

**Location:** GitHub repository (`data/` directory)

**DOI:** [Will be assigned via Zenodo upon publication]

---

### B. Differential Flux Results

**Description:** Complete list of differential metabolic reactions between GBM and normal conditions

**Format:** CSV / Excel

**Location:** Supplementary Table S1

**DOI:** [Will be assigned via Zenodo upon publication]

---

### C. Gene Knockout Results

**Description:** In silico knockout simulation results for 200+ metabolic genes

**Format:** CSV / Excel

**Location:** Supplementary Table S2

**DOI:** [Will be assigned via Zenodo upon publication]

---

### D. Survival Analysis Results

**Description:** Cox proportional hazards and Kaplan-Meier results for all metabolic genes

**Format:** CSV / Excel

**Location:** Supplementary Table S3

**DOI:** [Will be assigned via Zenodo upon publication]

---

### E. Literature Curation Database

**Description:** Manually curated evidence for top metabolic targets

**Format:** Excel

**Location:** Supplementary Table S4

**DOI:** [Will be assigned via Zenodo upon publication]

---

## Code Availability

### ARES-GBM Pipeline

**Repository:** GitHub

**URL:** https://github.com/ARES-Sentinel-GBM/-gbm-research

**Version:** v1.0.0

**License:** MIT License (Open Source)

**DOI:** [Will be assigned via Zenodo upon publication]

**Archive:** All code will be archived on Zenodo upon manuscript acceptance

**Docker Image:**
- **Docker Hub:** [Will be uploaded upon acceptance]
- **Command:** `docker pull aresbio/ares-gbm:latest`

---

## Data Access Timeline

| Data Type | Availability |
|-----------|--------------|
| Source code | ✅ Available now (GitHub) |
| Example data | ✅ Available now (GitHub) |
| Processed data | 🟡 Upon publication (Zenodo) |
| Complete results | 🟡 Upon publication (Supplementary) |
| Docker image | 🟡 Upon publication (Docker Hub) |

---

## Data Citation

For citing the data and code associated with this manuscript:

**Software:**
> Smith J, Rossi M, Bianchi M, Verdi G. ARES-SENTINEL-GBM v1.0.0. GitHub. 2026. https://github.com/ARES-Sentinel-GBM/-gbm-research

**Data (Zenodo):**
> Smith J, Rossi M, Bianchi M, Verdi G. Data for: ARES-SENTINEL-GBM. Zenodo. 2026. [DOI will be assigned]

---

## Contact

**Data Requests:** admin@ares-bio.com

**Technical Support:** https://github.com/ARES-Sentinel-GBM/-gbm-research/issues

---

## Compliance

This manuscript complies with:
- ✅ **FAIR Data Principles** (Findable, Accessible, Interoperable, Reusable)
- ✅ **TCGA Data Use Policies** (Open access data only)
- ✅ **Bioinformatics Data Availability Policy**
- ✅ **FORCE11 Software Citation Principles**

---

## Ethical Considerations

- No primary data collection was performed
- All data are from publicly available sources (TCGA, GTEx)
- No patient identifiers are included
- IRB approval was not required for this study
