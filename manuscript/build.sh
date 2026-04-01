#!/bin/bash
# Build script for ARES-GBM manuscript
# Requires LaTeX distribution (TeX Live, MacTeX, etc.)

echo "============================================"
echo "ARES-GBM Manuscript Build Script"
echo "============================================"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Check if pdflatex is available
if ! command -v pdflatex &> /dev/null; then
    echo "ERROR: pdflatex not found!"
    echo "Please install TeX Live (https://tug.org/texlive/) or MacTeX (https://tug.org/mactex/)"
    exit 1
fi

echo "Compiling manuscript..."
echo ""

# First compilation (to generate aux files)
echo "[1/4] First compilation..."
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1

# Bibliography
echo "[2/4] Generating bibliography..."
bibtex main > /dev/null 2>&1

# Second compilation (to resolve references)
echo "[3/4] Second compilation..."
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1

# Final compilation (for correct page numbers)
echo "[4/4] Final compilation..."
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1

echo ""
echo "============================================"
if [ -f "main.pdf" ]; then
    echo "SUCCESS: Manuscript compiled to main.pdf"
    echo ""
    echo "Output files:"
    echo "  - main.pdf (final manuscript)"
    echo "  - main.log (compilation log)"
    echo "  - main.aux (auxiliary file)"
    echo "  - main.bbl (bibliography)"
    echo "  - main.blg (bibliography log)"
    
    # Open PDF based on OS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        open main.pdf
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        xdg-open main.pdf 2>/dev/null || echo "Open main.pdf manually"
    else
        echo "Open main.pdf manually"
    fi
else
    echo "ERROR: Compilation failed!"
    echo "Check main.log for details."
    exit 1
fi
echo "============================================"
echo ""
