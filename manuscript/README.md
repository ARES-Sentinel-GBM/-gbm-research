# ARES-GBM Manuscript

LaTeX manuscript template for submission to **Bioinformatics** (Biovirt).

## Files

| File | Description |
|------|-------------|
| `main.tex` | Main manuscript document |
| `references.bib` | Bibliography database |
| `build.bat` | Windows build script |
| `build.sh` | Linux/Mac build script (to be created) |

## Requirements

To compile the manuscript, you need a LaTeX distribution:

- **Windows**: [MiKTeX](https://miktex.org/) or [TeX Live](https://tug.org/texlive/)
- **Mac**: [MacTeX](https://tug.org/mactex/)
- **Linux**: [TeX Live](https://tug.org/texlive/)

## Quick Start

### Windows

1. Install MiKTeX or TeX Live
2. Double-click `build.bat`
3. The PDF will open automatically

### Command Line (all platforms)

```bash
cd manuscript

# Compile
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex

# Open PDF
open main.pdf        # Mac
start main.pdf       # Windows
xdg-open main.pdf    # Linux
```

## Manuscript Structure

```
main.tex
├── Title and Authors
├── Abstract
├── Introduction
├── Materials and Methods
│   ├── Pipeline Overview
│   ├── Metabolic Model
│   ├── Context-Specific Reconstruction
│   ├── Flux Balance Analysis
│   ├── Gene Knockout Analysis
│   ├── Survival Analysis
│   └── Benchmark Analysis
├── Results
│   ├── Pipeline Overview
│   ├── Benchmark Performance
│   ├── Validation Results
│   ├── Pathway Enrichment
│   ├── Network Topology
│   └── Clinical Validation
├── Discussion
│   ├── Key Findings
│   ├── Comparison with Existing Tools
│   ├── Limitations
│   ├── Future Directions
│   └── Conclusions
├── Acknowledgements
├── Funding
├── Conflict of Interest
└── References
```

## Figures

The manuscript references figures from the `results/` directory:

- `figure1_benchmark.pdf` - Benchmark comparison
- `figure2_validation.pdf` - Validation results
- `figure3_pathways.pdf` - Pathway enrichment
- `figure4_network.pdf` - Network topology
- `figure5_literature.pdf` - Clinical validation

Generate figures by running:
```bash
python scripts/manuscript_figures.py
```

## Customization

### Update Author Information

Edit lines 47-52 in `main.tex`:

```latex
\author[1,2]{John Smith\thanks{Corresponding author: john.smith@email.com}}
\author[1]{Maria Rossi}
\author[2]{Marco Bianchi}
\author[1,*]{Giuseppe Verdi}
```

### Update Abstract

Edit the abstract environment (lines 60-80).

### Add/Modify Sections

Add new sections using standard LaTeX commands:

```latex
\section{New Section}
\subsection{New Subsection}
```

### Add Citations

1. Add entry to `references.bib`
2. Cite in text using `\citep{key}` or `\citet{key}`

## Troubleshooting

### Missing Packages

If you get "Package not found" errors, install missing packages:

**MiKTeX**: Use MiKTeX Console to install packages automatically
**TeX Live**: Use `tlmgr install <package-name>`

### Figure Not Found

Ensure figures are in the correct location:
```
ares-gbm-research/
├── manuscript/
│   └── main.tex
└── results/
    └── figure*.pdf
```

The paths in `main.tex` are relative: `../results/figure1_benchmark.pdf`

### Bibliography Not Showing

Run the compilation sequence again:
1. `pdflatex main.tex`
2. `bibtex main`
3. `pdflatex main.tex`
4. `pdflatex main.tex`

## Submission Checklist

- [ ] All figures generated and referenced
- [ ] All citations in references.bib
- [ ] Author affiliations correct
- [ ] Funding information updated
- [ ] Word count within limits (check journal requirements)
- [ ] Supplementary materials prepared (if needed)
- [ ] Cover letter prepared

## License

This manuscript template is provided under the MIT License.
