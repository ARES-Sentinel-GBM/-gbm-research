# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Gene Knockout Analysis

Simulated gene knockout effects on cancer cell growth.

Supports single gene knockouts and combinatorial knockouts to identify
synthetic lethal interactions for drug target discovery.
"""
import os
import sys
import argparse
import json
from datetime import datetime
from typing import List, Optional, Dict

import pandas as pd
import numpy as np

try:
    import cobra
    from cobra import manipulation
    from cobra.io import read_sbml_model
except ImportError:
    print("Error: cobrapy is required. Install with: pip install cobrapy")
    sys.exit(1)

from utils import (
    prepare_model,
    get_gene_symbol_map,
    optimize_model,
    plot_knockout_results,
    save_plot,
    ensure_dir
)


def run_gene_knockout(
    genes: List[str],
    combos: Optional[List[str]] = None,
    model_path: Optional[str] = None,
    solver: str = "glpk",
    output_dir: str = "results",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Simulate gene knockout effects on cell growth.
    
    Args:
        genes: List of gene symbols to knock out (single KOs)
        combos: List of combinations like "GENE1+GENE2" for combinatorial KOs
        model_path: Path to SBML model (optional, uses demo model if not provided)
        solver: COBRA solver to use
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        DataFrame with knockout results
    """
    print("=" * 60)
    print("ARES-GBM Gene Knockout Analysis")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    if combos is None:
        combos = []
    
    # Load metabolic model
    if model_path and os.path.exists(model_path):
        if verbose:
            print(f"Loading model from: {model_path}")
        model = read_sbml_model(model_path)
    else:
        if verbose:
            print("No model provided. Using demo mode with mock data.")
            print("For real analysis, provide a Human-GEM SBML model.")
        return _create_mock_ko_results(genes, combos, output_dir)
    
    # Prepare model
    model = prepare_model(model, solver=solver)
    
    # Build gene symbol map
    symbol_map = get_gene_symbol_map(model)
    
    # Compute wild-type solution
    if verbose:
        print("Computing wild-type solution...")
    wt_obj, _ = optimize_model(model)
    if verbose:
        print(f"  WT objective: {wt_obj:.4f}")
        print()
    
    results = []
    
    # Single gene knockouts
    if genes:
        if verbose:
            print(f"Running {len(genes)} single gene knockouts...")
        
        for symbol in genes:
            ensg = symbol_map.get(symbol)
            if not ensg:
                if verbose:
                    print(f"  ⚠️  {symbol}: not found in model")
                continue
            
            with model:
                try:
                    manipulation.knock_out_model_genes(model, [ensg])
                    sol_obj, _ = optimize_model(model)
                    
                    if wt_obj > 0 and sol_obj is not None:
                        ratio = sol_obj / wt_obj
                    else:
                        ratio = 0
                    
                    effect = classify_effect(ratio)
                    
                    results.append({
                        'type': 'single',
                        'genes': symbol,
                        'ensg': ensg,
                        'ratio': round(ratio, 4),
                        'effect': effect
                    })
                    
                    if verbose:
                        print(f"  {symbol}: ratio={ratio:.4f} [{effect}]")
                        
                except Exception as e:
                    if verbose:
                        print(f"  {symbol}: error - {e}")
        
        if verbose:
            print()
    
    # Combinatorial knockouts
    if combos:
        if verbose:
            print(f"Running {len(combos)} combinatorial knockouts...")
        
        for combo in combos:
            symbols = [s.strip() for s in combo.split('+')]
            ensgs = [symbol_map[s] for s in symbols if symbol_map.get(s)]
            
            if not ensgs:
                if verbose:
                    print(f"  ⚠️  {combo}: no valid genes found")
                continue
            
            with model:
                try:
                    manipulation.knock_out_model_genes(model, ensgs)
                    sol_obj, _ = optimize_model(model)
                    
                    if wt_obj > 0 and sol_obj is not None:
                        ratio = sol_obj / wt_obj
                    else:
                        ratio = 0
                    
                    effect = classify_effect(ratio)
                    
                    results.append({
                        'type': 'combo',
                        'genes': combo,
                        'ensg': '+'.join(ensgs),
                        'ratio': round(ratio, 4),
                        'effect': effect
                    })
                    
                    if verbose:
                        print(f"  {combo}: ratio={ratio:.4f} [{effect}]")
                        
                except Exception as e:
                    if verbose:
                        print(f"  {combo}: error - {e}")
        
        if verbose:
            print()
    
    # Create results DataFrame
    df = pd.DataFrame(results)
    
    if len(df) > 0:
        # Sort by ratio (most lethal first)
        df = df.sort_values('ratio')
        
        # Save results
        ensure_dir(os.path.join(output_dir, 'gene_ko.csv'))
        output_path = os.path.join(output_dir, 'gene_ko.csv')
        df.to_csv(output_path, index=False)
        if verbose:
            print(f"Results saved to: {output_path}")
        
        # Generate plot
        fig = plot_knockout_results(
            results,
            title="Gene Knockout — Growth Ratio"
        )
        plot_path = os.path.join(output_dir, 'gene_ko_results.png')
        save_plot(fig, plot_path)
        
        # Summary statistics
        if verbose:
            print()
            print("Summary:")
            print(f"  Total KOs tested: {len(df)}")
            print(f"  Lethal (ratio < 0.1): {len(df[df['effect'] == 'LETHAL'])}")
            print(f"  Strong (ratio < 0.5): {len(df[df['effect'] == 'STRONG'])}")
            print(f"  Weak (ratio >= 0.5): {len(df[df['effect'] == 'WEAK'])}")
            print()
    else:
        if verbose:
            print("No results generated. Check gene symbols.")
    
    # Show top results
    if len(df) > 0:
        print("Top 10 Most Lethal Knockouts:")
        print("-" * 60)
        print(df.head(10).to_string())
        print()
    
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return df


def classify_effect(ratio: float) -> str:
    """
    Classify knockout effect based on growth ratio.
    
    Args:
        ratio: Growth ratio (KO/WT)
    
    Returns:
        Effect classification string
    """
    if ratio < 0.1:
        return "LETHAL"
    elif ratio < 0.5:
        return "STRONG"
    else:
        return "WEAK"


def _create_mock_ko_results(
    genes: List[str],
    combos: List[str],
    output_dir: str
) -> pd.DataFrame:
    """
    Create mock knockout results for demonstration.
    
    Args:
        genes: List of single genes
        combos: List of combinations
        output_dir: Directory for output files
    
    Returns:
        Mock results DataFrame
    """
    np.random.seed(42)
    
    results = []
    
    # Single KOs with varying effects
    for symbol in genes:
        ratio = np.random.uniform(0, 1)
        results.append({
            'type': 'single',
            'genes': symbol,
            'ensg': f'ENSG{np.random.randint(10000, 99999)}',
            'ratio': round(ratio, 4),
            'effect': classify_effect(ratio)
        })
    
    # Combo KOs (tend to be more lethal)
    for combo in combos:
        ratio = np.random.uniform(0, 0.7)  # More likely to be lethal
        results.append({
            'type': 'combo',
            'genes': combo,
            'ensg': '+'.join([f'ENSG{np.random.randint(10000, 99999)}' for _ in combo.split('+')]),
            'ratio': round(ratio, 4),
            'effect': classify_effect(ratio)
        })
    
    df = pd.DataFrame(results)
    df = df.sort_values('ratio')
    
    # Save mock results
    ensure_dir(os.path.join(output_dir, 'gene_ko_demo.csv'))
    output_path = os.path.join(output_dir, 'gene_ko_demo.csv')
    df.to_csv(output_path, index=False)
    print(f"Demo results saved to: {output_path}")
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Gene Knockout - Simulate gene knockout effects',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gene_ko.py --genes RRM1 RRM2 TYMS GLS LDHA
  python gene_ko.py --genes RRM1 GLS --combos "LDHA+RRM2" "GLS+RRM2"
  python gene_ko.py --genes RRM1 RRM2 --model Human-GEM.xml --output my_results
        """
    )
    
    parser.add_argument('--genes', nargs='+', required=True,
                        help='Genes to knock out (single KOs)')
    parser.add_argument('--combos', nargs='+', default=[],
                        help='Gene combinations (e.g., "GENE1+GENE2")')
    parser.add_argument('--model', help='Path to SBML model (optional)')
    parser.add_argument('--solver', default='glpk', help='COBRA solver (default: glpk)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_gene_knockout(
        genes=args.genes,
        combos=args.combos,
        model_path=args.model,
        solver=args.solver,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
