@echo off
REM Build script for ARES-GBM manuscript
REM Requires LaTeX distribution (TeX Live, MiKTeX, etc.)

echo ============================================
echo ARES-GBM Manuscript Build Script
echo ============================================
echo.

cd /d "%~dp0"

REM Check if pdflatex is available
where pdflatex >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: pdflatex not found!
    echo Please install TeX Live (https://tug.org/texlive/) or MiKTeX (https://miktex.org/)
    pause
    exit /b 1
)

echo Compiling manuscript...
echo.

REM First compilation (to generate aux files)
echo [1/3] First compilation...
pdflatex -interaction=nonstopmode main.tex >nul 2>nul

REM Bibliography
echo [2/3] Generating bibliography...
bibtex main >nul 2>nul

REM Second compilation (to resolve references)
echo [3/3] Final compilation...
pdflatex -interaction=nonstopmode main.tex >nul 2>nul
pdflatex -interaction=nonstopmode main.tex >nul 2>nul

echo.
echo ============================================
if exist main.pdf (
    echo SUCCESS: Manuscript compiled to main.pdf
    echo.
    echo Output files:
    echo   - main.pdf (final manuscript)
    echo   - main.log (compilation log)
    echo   - main.aux (auxiliary file)
    echo   - main.bbl (bibliography)
    echo   - main.blg (bibliography log)
    start main.pdf
) else (
    echo ERROR: Compilation failed!
    echo Check main.log for details.
)
echo ============================================
echo.

pause
