# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Survival Analysis

Kaplan-Meier survival analysis based on gene expression.

Analyzes the association between gene expression levels and patient survival
using Cox proportional hazards model and Kaplan-Meier curves.
"""
import os
import sys
import argparse
import re
from datetime import datetime
from typing import Optional, Tuple

import pandas as pd
import numpy as np

try:
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test
    from lifelines import KaplanMeierFitter
except ImportError:
    print("Error: lifelines is required. Install with: pip install lifelines")
    sys.exit(1)

from utils import (
    load_expression_matrix,
    load_survival_data,
    plot_kaplan_meier,
    save_plot,
    ensure_dir
)


def run_survival_analysis(
    gene: str,
    expr_path: str,
    survival_path: str,
    percentile: int = 50,
    output_dir: str = "results",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform survival analysis for a single gene.
    
    Args:
        gene: Gene symbol to analyze
        expr_path: Path to gene expression CSV
        survival_path: Path to clinical survival data CSV
        percentile: Percentile to split high/low expression groups
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        DataFrame with survival analysis results
    """
    print("=" * 60)
    print("ARES-GBM Survival Analysis")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load survival data
    if verbose:
        print("Loading survival data...")
    surv_df = load_survival_data(survival_path)
    if verbose:
        print(f"  Samples: {len(surv_df)}")
        print(f"  Events: {surv_df['event'].sum()}")
        print()
    
    # Load expression data
    if verbose:
        print("Loading expression data...")
    expr_df = load_expression_matrix(expr_path)
    
    # Use mean if multiple samples
    if expr_df.shape[1] > 1:
        expr_series = expr_df.mean(axis=1)
    else:
        expr_series = expr_df.iloc[:, 0]
    
    if verbose:
        print(f"  Genes: {len(expr_series)}")
        print()
    
    # Find gene in expression data
    # Try exact match first, then case-insensitive
    gene_found = None
    
    if gene in expr_series.index:
        gene_found = gene
    else:
        # Try case-insensitive
        for idx in expr_series.index:
            if idx.upper() == gene.upper():
                gene_found = idx
                break
    
    if not gene_found:
        print(f"Error: Gene '{gene}' not found in expression data.")
        print(f"Available genes: {list(expr_series.index[:10])}...")
        return pd.DataFrame()
    
    if verbose:
        print(f"Gene found: {gene_found}")
    
    # Merge expression with survival data
    common_samples = surv_df.index.intersection(expr_df.columns)
    if len(common_samples) == 0:
        print("Error: No common samples between expression and survival data.")
        return pd.DataFrame()
    
    expr_values = expr_series.loc[gene_found].reindex(common_samples).values.astype(float)
    
    # Create analysis DataFrame
    analysis_df = surv_df.loc[common_samples, ['survival_days', 'event']].copy()
    analysis_df[gene] = expr_values
    analysis_df = analysis_df.dropna()
    
    if verbose:
        print(f"Samples for analysis: {len(analysis_df)}")
        print()
    
    # Cox Proportional Hazards Model
    if verbose:
        print("Fitting Cox Proportional Hazards model...")
    
    cox_df = analysis_df[['survival_days', 'event', gene]].copy()
    cox_df[gene] = (cox_df[gene] - cox_df[gene].mean()) / cox_df[gene].std()
    
    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col='survival_days', event_col='event')
    
    hr = float(cph.hazard_ratios_.iloc[0])
    p_value = float(cph.summary['p'].iloc[0])
    ci_lower = float(cph.confidence_intervals_.loc[gene, 'lower 95%'])
    ci_upper = float(cph.confidence_intervals_.loc[gene, 'upper 95%'])
    
    if verbose:
        print(f"  Hazard Ratio: {hr:.3f}")
        print(f"  95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
        print(f"  p-value: {p_value:.4f}")
        print()
    
    # Kaplan-Meier Analysis
    if verbose:
        print(f"Computing Kaplan-Meier curves (split at {percentile}th percentile)...")
    
    cutoff = np.percentile(analysis_df[gene].values, percentile)
    
    high_expr = analysis_df[analysis_df[gene] >= cutoff].copy()
    low_expr = analysis_df[analysis_df[gene] < cutoff].copy()
    
    if verbose:
        print(f"  High expression group: n={len(high_expr)}")
        print(f"  Low expression group: n={len(low_expr)}")
    
    # Log-rank test
    lr_result = logrank_test(
        high_expr['survival_days'],
        low_expr['survival_days'],
        event_observed_A=high_expr['event'],
        event_observed_B=low_expr['event']
    )
    km_p_value = lr_result.p_value
    
    if verbose:
        print(f"  Log-rank p-value: {km_p_value:.4f}")
        print()
    
    # Generate Kaplan-Meier plot
    fig = plot_kaplan_meier(
        high_expr,
        low_expr,
        gene_name=gene,
        title=f"Kaplan-Meier — {gene} (HR={hr:.2f}, p={km_p_value:.4f})"
    )
    plot_path = os.path.join(output_dir, f'survival_{gene}.png')
    save_plot(fig, plot_path)
    
    # Prepare results DataFrame
    results = pd.DataFrame({
        'sample': list(high_expr.index) + list(low_expr.index),
        'group': ['HIGH'] * len(high_expr) + ['LOW'] * len(low_expr),
        'survival': list(high_expr['survival_days']) + list(low_expr['survival_days']),
        'event': list(high_expr['event']) + list(low_expr['event']),
        'gene': gene,
        'expression': list(high_expr[gene]) + list(low_expr[gene])
    })
    
    # Save results
    ensure_dir(os.path.join(output_dir, f'survival_{gene}.csv'))
    output_path = os.path.join(output_dir, f'survival_{gene}.csv')
    results.to_csv(output_path, index=False)
    if verbose:
        print(f"Results saved to: {output_path}")
    
    # Save summary
    summary = {
        'gene': gene,
        'hazard_ratio': hr,
        'hr_ci_lower': ci_lower,
        'hr_ci_upper': ci_upper,
        'cox_p_value': p_value,
        'km_p_value': km_p_value,
        'cutoff': cutoff,
        'n_high': len(high_expr),
        'n_low': len(low_expr),
        'n_total': len(analysis_df)
    }
    
    summary_path = os.path.join(output_dir, f'survival_{gene}_summary.json')
    import json
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    if verbose:
        print(f"Summary saved to: {summary_path}")
    
    # Print conclusion
    print()
    print("Conclusion:")
    print("-" * 60)
    if km_p_value < 0.05:
        sig_level = "***" if km_p_value < 0.001 else "**" if km_p_value < 0.01 else "*"
        print(f"  {gene} expression is SIGNIFICANTLY associated with survival (p={km_p_value:.4f}) {sig_level}")
        if hr > 1:
            print(f"  High {gene} expression → Worse prognosis (HR={hr:.2f})")
        else:
            print(f"  High {gene} expression → Better prognosis (HR={hr:.2f})")
    else:
        print(f"  {gene} expression is NOT significantly associated with survival (p={km_p_value:.4f})")
    
    print()
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Survival Analysis - Gene expression vs patient survival',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python survival_analysis.py --gene GLS --expr data/gbm.csv --survival data/survival.csv
  python survival_analysis.py --gene LDHA --expr data/gbm.csv --survival data/survival.csv --percentile 50
  python survival_analysis.py --gene RRM2 --expr data/gbm.csv --survival data/survival.csv --output my_results
        """
    )
    
    parser.add_argument('--gene', required=True, help='Gene symbol to analyze')
    parser.add_argument('--expr', required=True, help='Path to expression data CSV')
    parser.add_argument('--survival', required=True, help='Path to survival data CSV')
    parser.add_argument('--percentile', type=int, default=50, help='Percentile split (default: 50)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_survival_analysis(
        gene=args.gene,
        expr_path=args.expr,
        survival_path=args.survival,
        percentile=args.percentile,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
