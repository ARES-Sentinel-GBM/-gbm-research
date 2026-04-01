# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Flux Analysis

Differential flux analysis between GBM and normal astrocytes using iMAT algorithm.

Compares metabolic fluxes under two conditions to identify differential reactions
that could serve as potential drug targets.
"""
import os
import sys
import argparse
import json
from datetime import datetime
from typing import Optional

import pandas as pd
import numpy as np

try:
    import cobra
    from cobra import Model
    from cobra.io import read_sbml_model
except ImportError:
    print("Error: cobrapy is required. Install with: pip install cobrapy")
    sys.exit(1)

from utils import (
    load_expression_matrix,
    prepare_model,
    apply_imat,
    get_gene_symbol_map,
    optimize_model,
    plot_flux_comparison,
    save_plot,
    ensure_dir
)


def run_flux_analysis(
    gbm_expr_path: str,
    astro_expr_path: str,
    model_path: Optional[str] = None,
    high_percentile: float = 75,
    low_percentile: float = 25,
    delta_threshold: float = 0.01,
    solver: str = "glpk",
    output_dir: str = "results",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform differential flux analysis between GBM and astrocyte conditions.
    
    Args:
        gbm_expr_path: Path to GBM expression CSV
        astro_expr_path: Path to astrocyte expression CSV
        model_path: Path to SBML model (optional, uses demo model if not provided)
        high_percentile: Percentile for high expression threshold
        low_percentile: Percentile for low expression threshold
        delta_threshold: Minimum |delta| to include in results
        solver: COBRA solver to use
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        DataFrame with differential flux results
    """
    print("=" * 60)
    print("ARES-GBM Flux Analysis")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load expression data
    if verbose:
        print("Loading expression data...")
    
    gbm_expr = load_expression_matrix(gbm_expr_path)
    astro_expr = load_expression_matrix(astro_expr_path)
    
    # Use mean across samples if multiple columns
    if gbm_expr.shape[1] > 1:
        gbm_expr = gbm_expr.mean(axis=1)
    else:
        gbm_expr = gbm_expr.iloc[:, 0]
    
    if astro_expr.shape[1] > 1:
        astro_expr = astro_expr.mean(axis=1)
    else:
        astro_expr = astro_expr.iloc[:, 0]
    
    if verbose:
        print(f"  GBM samples: {gbm_expr.shape[0]} genes")
        print(f"  Astrocyte samples: {astro_expr.shape[0]} genes")
        print()
    
    # Load or create metabolic model
    if model_path and os.path.exists(model_path):
        if verbose:
            print(f"Loading model from: {model_path}")
        model = read_sbml_model(model_path)
    else:
        if verbose:
            print("No model provided. Using demo mode with mock data.")
            print("For real analysis, provide a Human-GEM SBML model.")
        # Create mock results for demonstration
        return _create_mock_flux_results(delta_threshold, output_dir)
    
    # Prepare model
    model = prepare_model(model, solver=solver)
    
    # Get wild-type solution
    if verbose:
        print("Computing wild-type solution...")
    wt_obj, wt_fluxes = optimize_model(model)
    if verbose:
        print(f"  WT objective: {wt_obj:.4f}")
        print()
    
    # Apply iMAT for GBM
    if verbose:
        print("Applying iMAT for GBM condition...")
    model_gbm = apply_imat(
        model,
        gbm_expr,
        high_percentile=high_percentile,
        low_percentile=low_percentile,
        verbose=verbose
    )
    obj_gbm, fluxes_gbm = optimize_model(model_gbm)
    if verbose:
        print(f"  GBM objective: {obj_gbm:.4f}")
    
    # Apply iMAT for Astrocyte
    if verbose:
        print("Applying iMAT for Astrocyte condition...")
    model_astro = apply_imat(
        model,
        astro_expr,
        high_percentile=high_percentile,
        low_percentile=low_percentile,
        verbose=verbose
    )
    obj_astro, fluxes_astro = optimize_model(model_astro)
    if verbose:
        print(f"  Astrocyte objective: {obj_astro:.4f}")
        print()
    
    # Compute differential fluxes
    if verbose:
        print("Computing differential fluxes...")
    
    df = pd.DataFrame({
        'flux_gbm': fluxes_gbm,
        'flux_astro': fluxes_astro
    }).fillna(0)
    
    df['delta'] = df['flux_gbm'] - df['flux_astro']
    df['abs_delta'] = df['delta'].abs()
    
    # Filter by threshold
    df = df[df['abs_delta'] > delta_threshold]
    df = df.sort_values('abs_delta', ascending=False)
    
    if verbose:
        print(f"  Total reactions: {len(df)}")
        print(f"  Reactions with |Δ| > {delta_threshold}: {len(df)}")
        print()
    
    # Save results
    ensure_dir(os.path.join(output_dir, 'flux_analysis.csv'))
    output_path = os.path.join(output_dir, 'flux_analysis.csv')
    df.to_csv(output_path)
    if verbose:
        print(f"Results saved to: {output_path}")
    
    # Generate plot
    if len(df) > 0:
        fig = plot_flux_comparison(
            df['flux_gbm'],
            df['flux_astro'],
            top_n=20,
            title="Top Differential Reactions (GBM vs Astrocyte)"
        )
        plot_path = os.path.join(output_dir, 'flux_analysis_top20.png')
        save_plot(fig, plot_path)
    
    # Show top results
    print()
    print("Top 10 Differential Reactions:")
    print("-" * 60)
    print(df.head(10).to_string())
    print()
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return df


def _create_mock_flux_results(
    delta_threshold: float,
    output_dir: str
) -> pd.DataFrame:
    """
    Create mock flux analysis results for demonstration.
    
    Args:
        delta_threshold: Minimum |delta| to include
        output_dir: Directory for output files
    
    Returns:
        Mock results DataFrame
    """
    np.random.seed(42)
    
    # Create mock reaction data
    reactions = [
        'R_GK', 'R_PYK', 'R_LDH', 'R_G6PDH', 'R_GLS', 'R_GLS2',
        'R_IDH', 'R_MDH', 'R_SDH', 'R_FUM', 'R_ATPS', 'R_NDH',
        'R_PFK', 'R_ALDO', 'R_TPI', 'R_GAPD', 'R_PGK', 'R_PGM',
        'R_ENO', 'R_PPCK', 'R_PC', 'R_PEPCK', 'R_FBP', 'R_TKT',
        'R_TALA', 'R_RPE', 'R_RPI', 'R_GND', 'R_PGLS', 'R_SPS'
    ]
    
    data = []
    for rxn in reactions:
        flux_gbm = np.random.uniform(-10, 10)
        flux_astro = np.random.uniform(-10, 10)
        delta = flux_gbm - flux_astro
        if abs(delta) > delta_threshold:
            data.append({
                'reaction': rxn,
                'flux_gbm': flux_gbm,
                'flux_astro': flux_astro,
                'delta': delta,
                'abs_delta': abs(delta)
            })
    
    df = pd.DataFrame(data)
    df = df.sort_values('abs_delta', ascending=False)
    df = df.set_index('reaction')
    
    # Save mock results
    ensure_dir(os.path.join(output_dir, 'flux_analysis_demo.csv'))
    output_path = os.path.join(output_dir, 'flux_analysis_demo.csv')
    df.to_csv(output_path)
    print(f"Demo results saved to: {output_path}")
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Flux Analysis - Differential flux analysis for GBM',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python flux_analysis.py --gbm-expr data/gbm.csv --astro-expr data/astro.csv
  python flux_analysis.py --gbm-expr data/gbm.csv --astro-expr data/astro.csv --high 80 --low 20
  python flux_analysis.py --gbm-expr data/gbm.csv --astro-expr data/astro.csv --model Human-GEM.xml --output my_results
        """
    )
    
    parser.add_argument('--gbm-expr', required=True, help='Path to GBM expression CSV')
    parser.add_argument('--astro-expr', required=True, help='Path to astrocyte expression CSV')
    parser.add_argument('--model', help='Path to SBML model (optional)')
    parser.add_argument('--high', type=float, default=75, help='High expression percentile (default: 75)')
    parser.add_argument('--low', type=float, default=25, help='Low expression percentile (default: 25)')
    parser.add_argument('--delta', type=float, default=0.01, help='Delta threshold (default: 0.01)')
    parser.add_argument('--solver', default='glpk', help='COBRA solver (default: glpk)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_flux_analysis(
        gbm_expr_path=args.gbm_expr,
        astro_expr_path=args.astro_expr,
        model_path=args.model,
        high_percentile=args.high,
        low_percentile=args.low,
        delta_threshold=args.delta,
        solver=args.solver,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
