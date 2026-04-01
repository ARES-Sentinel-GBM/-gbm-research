# ARES-Sentinel-GBM Research Tools 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)
[![GitHub issues](https://img.shields.io/github/issues/ARES-Sentinel-GBM/-gbm-research.svg)](https://github.com/ARES-Sentinel-GBM/-gbm-research/issues)

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

### Option 1: Local Installation

#### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

#### 2. Run an Analysis

```bash
# Flux analysis
python scripts/flux_analysis.py --high 75 --low 25 --delta 0.01

# Gene knockout
python scripts/gene_ko.py --genes RRM1 RRM2 TYMS GLS

# Survival analysis
python scripts/survival_analysis.py --gene GLS --percentile 50
```

#### 3. Check Results

Results are saved to `results/` directory.

---

### Option 2: Docker (Recommended)

#### Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop) installed

#### Quick Start

```bash
# Build the image
docker build -t ares-gbm .

# Run flux analysis
docker run --rm -v $(pwd)/results:/app/results ares-gbm \
    python scripts/flux_analysis.py --high 75 --low 25 --delta 0.01

# Run gene knockout
docker run --rm -v $(pwd)/results:/app/results ares-gbm \
    python scripts/gene_ko.py --genes RRM2 GLS LDHA

# Run survival analysis
docker run --rm -v $(pwd)/results:/app/results ares-gbm \
    python scripts/survival_analysis.py --gene GLS --expr data/expression_example.csv --survival data/survival_example.csv
```

#### Using Docker Compose

```bash
# Run analysis
docker compose run ares-gbm python scripts/flux_analysis.py

# Start Jupyter notebook for interactive analysis
docker compose --profile jupyter up jupyter
```

Then open `http://localhost:8888` in your browser.

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
├── README.md                 # This file
├── LICENSE                   # MIT License
├── CITATION.md               # How to cite this work
├── requirements.txt          # Python dependencies
├── setup.py                  # Package setup
├── Dockerfile                # Docker image definition
├── docker-compose.yml        # Docker Compose configuration
├── .dockerignore             # Docker ignore rules
│
├── scripts/
│   ├── flux_analysis.py      # Differential flux analysis
│   ├── gene_ko.py            # Gene knockout simulation
│   ├── survival_analysis.py  # Kaplan-Meier survival analysis
│   ├── utils.py              # Helper functions
│   └── manuscript_figures.py # Generate manuscript figures
│
├── manuscript/
│   ├── main.tex              # LaTeX manuscript for Biovirt
│   ├── references.bib        # Bibliography database
│   ├── build.bat             # Windows build script
│   ├── build.sh              # Linux/Mac build script
│   └── README.md             # Manuscript instructions
│
├── notebooks/
│   └── demo_analysis.ipynb   # Interactive demo notebook
│
├── data/
│   ├── expression_example.csv    # Example gene expression
│   └── survival_example.csv      # Example survival data
│
└── results/                      # Output directory (generated)
    ├── figure1_benchmark.pdf
    ├── figure2_validation.pdf
    ├── figure3_pathways.pdf
    ├── figure4_network.pdf
    └── figure5_literature.pdf
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

```bibtex
@article{smith2026ares,
  title={ARES-Sentinel-GBM: Genome-scale metabolic modeling for drug target discovery in glioblastoma},
  author={Smith, John and Rossi, Maria and Bianchi, Marco and Verdi, Giuseppe},
  journal={bioRxiv},
  year={2026},
  url={https://github.com/ARES-Sentinel-GBM/-gbm-research}
}
```

**Software citation:**
> Giovanni S. **ARES-Sentinel-GBM: Genome-scale metabolic modeling for drug target discovery in glioblastoma.** *GitHub*, 2026. https://github.com/ARES-Sentinel-GBM/-gbm-research

---

## 📄 License

**MIT License** — Free for academic and research use.

See [LICENSE](LICENSE) for details.

---

## 🔧 Troubleshooting

### Common Issues

#### `cobrapy` installation fails

**Problem:** Error installing cobrapy on Windows

**Solution:**
```bash
# Install Microsoft C++ Build Tools first
# Then install with conda (recommended)
conda install -c conda-forge cobrapy

# Or use pre-built wheels
pip install cobrapy --only-binary :all:
```

#### GLPK solver not found

**Problem:** `RuntimeError: GLPK solver not available`

**Solution:**
```bash
# Windows: Install via conda
conda install -c conda-forge glpk

# Or use alternative solver
pip install scipy
# ARES-GBM will automatically use scipy's LP solver
```

#### Docker permission denied

**Problem:** `permission denied while trying to connect to Docker daemon socket`

**Solution:**
```bash
# Linux: Add user to docker group
sudo usermod -aG docker $USER
# Then log out and log back in

# Or run with sudo (not recommended for production)
sudo docker run ...
```

#### Memory error during flux analysis

**Problem:** `MemoryError` when running on large models

**Solution:**
```bash
# Use Docker with memory limit
docker run --memory="4g" --rm -v $(pwd)/results:/app/results ares-gbm ...

# Or reduce model size by filtering low-expression genes
python scripts/flux_analysis.py --min-expression 1.0 ...
```

#### LaTeX compilation fails

**Problem:** Manuscript doesn't compile with pdflatex

**Solution:**
1. Install full LaTeX distribution (TeX Live or MiKTeX)
2. Or use [Overleaf](https://overleaf.com) - upload `manuscript/` files
3. Missing packages will be auto-installed by MiKTeX

---

## 🤝 Contributing

Contributions are welcome! This is an open-source project for the GBM research community.

### How to Contribute

1. **Fork** the repository
2. **Create a branch** for your feature (`git checkout -b feature/amazing-feature`)
3. **Make your changes** and test thoroughly
4. **Commit** your changes (`git commit -m 'Add amazing feature'`)
5. **Push** to your branch (`git push origin feature/amazing-feature`)
6. **Open a Pull Request**

### Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR-USERNAME/-gbm-research.git
cd -gbm-research

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -r requirements.txt
pip install -e .

# Install development dependencies
pip install pytest black flake8 mypy
```

### Code Style

- Follow [PEP 8](https://pep8.org/) style guidelines
- Use type hints for function signatures
- Write docstrings for all public functions
- Add tests for new features

### Testing

```bash
# Run tests
pytest tests/

# Run linting
black --check scripts/
flake8 scripts/
```

### Reporting Issues

- Use [GitHub Issues](https://github.com/ARES-Sentinel-GBM/-gbm-research/issues)
- Include Python version, OS, and error traceback
- For bugs, provide minimal reproducible example

### Pull Request Guidelines

- **Features:** Add tests and documentation
- **Bug fixes:** Add test case demonstrating the fix
- **Documentation:** Update README and docstrings
- **Breaking changes:** Document migration path

---

## 🤝 Support

- **📚 Documentation:** Browse the [manuscript/](manuscript/) folder for LaTeX template
- **🐛 Issues:** Open an issue on [GitHub Issues](https://github.com/ARES-Sentinel-GBM/-gbm-research/issues)
- **📧 Email:** admin@ares-bio.com
- **🌐 Full SaaS Platform:** [ares-bio.com](https://ares-bio.com) (advanced features, cloud execution, priority support)

---

**Built for the glioblastoma research community** 🎗️

*For academic and non-commercial research use. For commercial applications, please contact admin@ares-bio.com for licensing options.*
