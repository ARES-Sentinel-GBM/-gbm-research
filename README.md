# ARES-Sentinel-GBM Research Tools 🧬

**Free Python scripts for glioblastoma metabolic analysis.**

A collection of Python scripts for genome-scale metabolic modeling (GMM) of glioblastoma (GBM). Part of the ARES-Sentinel-GBM project.

> ⚠️ **Note**: This is the **research edition** (Python scripts only).  
> For the full web-based platform with authentication, results management, and cloud execution, visit: **[ares-bio.com](https://ares-bio.com)**

---

## 💻 Free vs SaaS Versions

| Feature | **This Version (Free)** | **SaaS Platform** |
|---------|------------------------|-------------------|
| **Scripts Python** | ✅ Included | ✅ Included |
| **Web Interface** | ❌ Command-line only | ✅ Streamlit app |
| **Authentication** | ❌ Local execution | ✅ User accounts |
| **Results Storage** | ❌ Local files | ✅ Cloud storage |
| **Job Queue** | ❌ Manual execution | ✅ Background jobs |
| **Data Upload** | ❌ Manual file handling | ✅ Encrypted upload |
| **Support** | ❌ Community (GitHub) | ✅ Email support |
| **Cost** | **Free** (MIT License) | From €49/month |

**For students and academic researchers**: This free version is fully functional for thesis and research projects.

**For labs and companies**: The SaaS platform provides team management, data privacy, and priority support.

👉 **Try the SaaS platform**: [https://ares-bio.com](https://ares-bio.com)

---

## 📦 What's Included

- **`flux_analysis.py`** — Differential flux analysis between GBM and astrocytes
- **`gene_ko.py`** — Simulated gene knockout effects on cell growth
- **`survival_analysis.py`** — Survival analysis based on gene expression
- **`utils.py`** — Helper functions for data loading and visualization

---

## 🚀 Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Run an Analysis

```bash
# Flux analysis
python scripts/flux_analysis.py --high 75 --low 25 --delta 0.01

# Gene knockout
python scripts/gene_ko.py --genes RRM1 RRM2 TYMS GLS

# Survival analysis
python scripts/survival_analysis.py --gene GLS --percentile 50
```

### 3. Check Results

Results are saved to `results/` directory.

---

## 📋 Requirements

- Python 3.8+
- See `requirements.txt` for full list

**Core dependencies:**
- `cobrapy` — Metabolic modeling
- `pandas`, `numpy` — Data analysis
- `lifelines` — Survival analysis
- `plotly` — Interactive visualizations

---

## 📊 Example Usage

### Flux Analysis

Compare metabolic fluxes between GBM and normal astrocytes:

```bash
python scripts/flux_analysis.py \
    --gbm-expr data/gse4412_gbm_ensg.csv \
    --astro-expr data/gse4412_astro_ensg.csv \
    --high 75 \
    --low 25 \
    --delta 0.01 \
    --output results/flux_results.csv
```

### Gene Knockout

Simulate single and combinatorial gene knockouts:

```bash
python scripts/gene_ko.py \
    --genes RRM1 RRM2 TYMS GLS LDHA HK2 \
    --combos "LDHA+RRM2" "GLS+RRM2" \
    --output results/ko_results.csv
```

### Survival Analysis

Kaplan-Meier survival curves based on gene expression:

```bash
python scripts/survival_analysis.py \
    --gene GLS \
    --expr data/gse4412_gbm_ensg.csv \
    --survival data/survival_cox_data.csv \
    --percentile 50 \
    --output results/survival_results.csv
```

---

## 📁 Directory Structure

```
ares-gbm-research/
├── README.md
├── requirements.txt
├── scripts/
│   ├── flux_analysis.py
│   ├── gene_ko.py
│   ├── survival_analysis.py
│   └── utils.py
├── notebooks/
│   └── demo_analysis.ipynb
├── data/
│   ├── gse4412_gbm_ensg.csv       # GBM expression (example)
│   ├── gse4412_astro_ensg.csv     # Astrocyte expression (example)
│   └── survival_cox_data.csv      # Clinical survival data (example)
└── results/                       # Output directory
```

---

## 🔬 Methods

These tools implement the iMAT algorithm for context-specific metabolic modeling:

1. **Gene expression integration** — Transcriptomics data mapped to metabolic reactions
2. **Flux Balance Analysis** — Constraint-based optimization of metabolic networks
3. **Differential analysis** — Identify reactions with significantly different fluxes
4. **Knockout simulation** — Predict effects of gene deletions on growth
5. **Survival analysis** — Cox proportional hazards and Kaplan-Meier curves

**Genome-scale model:** Human-GEM (available at [github.com/SysBioChalmers/human-gem](https://github.com/SysBioChalmers/human-gem))

---

## 📚 Citation

If you use these tools in your research, please cite:

> Giovanni S. **ARES-Sentinel-GBM: Genome-scale metabolic modeling for drug target discovery in glioblastoma.** *bioRxiv*, 2026.  
> GitHub: [github.com/ARES-Sentinel-GBM/ARES-Sentinel-GBM](https://github.com/ARES-Sentinel-GBM/ARES-Sentinel-GBM)

---

## 📄 License

**MIT License** — Free for academic and research use.

See [LICENSE](LICENSE) for details.

---

## 🤝 Support

- **Issues:** Open an issue on GitHub
- **Email:** admin@ares-bio.com
- **Full SaaS Platform:** [ares-bio.com](https://ares-bio.com) (advanced features, cloud execution)

---

**Built for the glioblastoma research community** 🎗️
