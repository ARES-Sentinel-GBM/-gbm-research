# ARES-GBM Research Tools

Python scripts for genome-scale metabolic modeling of glioblastoma.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run flux analysis
python scripts/flux_analysis.py --gbm-expr data/gbm.csv --astro-expr data/astro.csv

# Run gene knockout
python scripts/gene_ko.py --genes RRM1 RRM2 GLS --combos "LDHA+RRM2"

# Run survival analysis
python scripts/survival_analysis.py --gene GLS --expr data/gbm.csv --survival data/clinical.csv
```

## License

MIT License - Free for academic and research use.

## Full Documentation

See [README.md](README.md) for complete documentation.
