# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Synthetic Lethality Screen

Systematic screen for synthetic lethal interactions in GBM metabolism.
Identifies synergistic drug target combinations for combination therapy.
"""
import os
import sys
import argparse
import json
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from itertools import combinations

import pandas as pd
import numpy as np

try:
    import cobra
    from cobra import Model
    from cobra.io import read_sbml_model
except ImportError:
    print("Error: cobrapy is required. Install with: pip install cobrapy")
    sys.exit(1)

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    print("Warning: plotly is required for visualization. Install with: pip install plotly")
    go = None
    make_subplots = None

from utils import (
    load_expression_matrix,
    prepare_model,
    apply_imat,
    optimize_model,
    ensure_dir,
    save_plot
)


def run_synthetic_lethality(
    expr_path: str,
    gene_list: Optional[List[str]] = None,
    model_path: Optional[str] = None,
    high_percentile: float = 75,
    low_percentile: float = 25,
    max_combinations: int = 100,
    output_dir: str = "results",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Perform synthetic lethality screen for GBM metabolic targets.
    
    Args:
        expr_path: Path to gene expression CSV
        gene_list: List of genes to screen (optional, uses default if None)
        model_path: Path to SBML model (optional)
        high_percentile: Percentile for high expression threshold
        low_percentile: Percentile for low expression threshold
        max_combinations: Maximum number of combinations to test
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        DataFrame with synthetic lethality results
    """
    print("=" * 60)
    print("ARES-GBM Synthetic Lethality Screen")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load expression data
    if verbose:
        print("Loading expression data...")
    expr_df = load_expression_matrix(expr_path)
    
    if expr_df.shape[1] > 1:
        expr_series = expr_df.mean(axis=1)
    else:
        expr_series = expr_df.iloc[:, 0]
    
    if verbose:
        print(f"  Genes: {len(expr_series)}")
        print()
    
    # Define gene list for screening
    if gene_list is None:
        # Default GBM metabolic targets
        gene_list = [
            'HK2', 'PFKP', 'LDHA', 'GLS', 'GLUD1',
            'IDH1', 'IDH2', 'MDH1', 'MDH2',
            'RRM2', 'TYMS', 'DHFR', 'POLA1',
            'ACLY', 'FASN', 'ACACA', 'SCD',
            'PHGDH', 'PSAT1', 'PSPH', 'SHMT1', 'MTHFD2',
            'NDUFV1', 'SDHB', 'UQCRB', 'COX5A', 'ATP5A1'
        ]
    
    if verbose:
        print(f"Screening {len(gene_list)} genes...")
        print(f"Total combinations: {len(gene_list)*(len(gene_list)-1)//2}")
        print()
    
    # Load or create metabolic model
    if model_path and os.path.exists(model_path):
        if verbose:
            print(f"Loading model from: {model_path}")
        model = read_sbml_model(model_path)
        model = prepare_model(model)
        
        # Apply iMAT for GBM
        if verbose:
            print("Building context-specific GBM model...")
        model_gbm = apply_imat(
            model,
            expr_series,
            high_percentile=high_percentile,
            low_percentile=low_percentile,
            verbose=False
        )
    else:
        if verbose:
            print("No model provided. Using demo mode with mock data.")
        model = None
        model_gbm = None
    
    # Calculate wild-type growth
    if model is not None:
        wt_obj, _ = optimize_model(model_gbm)
    else:
        wt_obj = 1.0
    
    if verbose:
        print(f"  Reference growth rate: {wt_obj:.4f}")
        print()
    
    # Screen all pairwise combinations
    results = []
    total_combos = len(gene_list) * (len(gene_list) - 1) // 2
    tested_combos = min(total_combos, max_combinations)
    
    if verbose:
        print(f"Testing {tested_combos} combinations...")
    
    combo_idx = 0
    for gene1, gene2 in combinations(gene_list, 2):
        if combo_idx >= max_combinations:
            break
        
        # Simulate double knockout
        if model is not None:
            growth_rate = _simulate_double_ko(model_gbm, gene1, gene2)
        else:
            growth_rate = _mock_double_ko(gene1, gene2, wt_obj)
        
        # Calculate expected growth (Bliss independence)
        if model is not None:
            growth1 = _simulate_single_ko(model_gbm, gene1)
            growth2 = _simulate_single_ko(model_gbm, gene2)
        else:
            growth1 = _mock_single_ko(gene1, wt_obj)
            growth2 = _mock_single_ko(gene2, wt_obj)
        
        # Bliss independence expectation
        bliss_expected = growth1 * growth2
        
        # Calculate synergy score
        # Negative = synergistic (more lethal than expected)
        # Positive = antagonistic (less lethal than expected)
        synergy_score = growth_rate - bliss_expected
        
        # Classify interaction
        interaction_type = classify_interaction(synergy_score, growth_rate)
        
        # Calculate combination index (CI)
        # CI < 1 = synergy, CI = 1 = additive, CI > 1 = antagonism
        if growth1 > 0 and growth2 > 0 and growth_rate > 0:
            combination_index = (1 - growth_rate) / ((1 - growth1) + (1 - growth2) - (1 - growth1) * (1 - growth2))
        else:
            combination_index = None
        
        results.append({
            'gene1': gene1,
            'gene2': gene2,
            'single_ko_1': round(growth1, 4),
            'single_ko_2': round(growth2, 4),
            'double_ko': round(growth_rate, 4),
            'bliss_expected': round(bliss_expected, 4),
            'synergy_score': round(synergy_score, 4),
            'combination_index': round(combination_index, 4) if combination_index else None,
            'interaction_type': interaction_type,
            'synthetic_lethal': interaction_type == 'SYNTHETIC_LETHAL'
        })
        
        combo_idx += 1
        
        if verbose and combo_idx % 20 == 0:
            print(f"  Tested {combo_idx}/{tested_combos} combinations...")
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Sort by synergy score (most synergistic first)
    results_df = results_df.sort_values('synergy_score')
    
    # Save results
    ensure_dir(os.path.join(output_dir, 'synthetic_lethality.csv'))
    output_path = os.path.join(output_dir, 'synthetic_lethality.csv')
    results_df.to_csv(output_path, index=False)
    
    # Summary statistics (calculate before verbose block)
    n_synthetic_lethal = len(results_df[results_df['synthetic_lethal'] == True])
    n_synergistic = len(results_df[results_df['interaction_type'] == 'SYNERGISTIC'])
    n_additive = len(results_df[results_df['interaction_type'] == 'ADDITIVE'])
    n_antagonistic = len(results_df[results_df['interaction_type'] == 'ANTAGONISTIC'])
    
    if verbose:
        print(f"\nResults saved to: {output_path}")
        
        print(f"\nSummary:")
        print(f"  Synthetic Lethal: {n_synthetic_lethal}")
        print(f"  Synergistic: {n_synergistic}")
        print(f"  Additive: {n_additive}")
        print(f"  Antagonistic: {n_antagonistic}")
        
        print(f"\nTop 10 Synthetic Lethal Pairs:")
        print("-" * 80)
        top_sl = results_df[results_df['synthetic_lethal'] == True].head(10)
        if len(top_sl) > 0:
            print(top_sl[['gene1', 'gene2', 'double_ko', 'synergy_score']].to_string())
        else:
            print("  No synthetic lethal pairs found in demo mode")
        print()
    
    # Generate visualization
    try:
        fig = plot_synthetic_lethality(results_df)
        plot_path = os.path.join(output_dir, 'synthetic_lethality_results.png')
        save_plot(fig, plot_path)
        if verbose:
            print(f"Plot saved to: {plot_path}")
    except Exception as e:
        if verbose:
            print(f"Warning: Could not generate plot: {e}")
    
    # Save summary
    summary = {
        'total_tested': len(results_df),
        'synthetic_lethal': int(n_synthetic_lethal),
        'synergistic': int(n_synergistic),
        'additive': int(n_additive),
        'antagonistic': int(n_antagonistic),
        'top_pair': f"{results_df.iloc[0]['gene1']}+{results_df.iloc[0]['gene2']}" if len(results_df) > 0 else None,
        'top_synergy_score': float(results_df.iloc[0]['synergy_score']) if len(results_df) > 0 else None
    }
    
    summary_path = os.path.join(output_dir, 'synthetic_lethality_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return results_df


def _simulate_single_ko(model: Model, gene: str) -> float:
    """
    Simulate single gene knockout.
    
    Args:
        model: Metabolic model
        gene: Gene to knock out
    
    Returns:
        Growth rate relative to wild-type
    """
    with model:
        try:
            cobra.manipulation.knock_out_model_genes(model, [gene])
            obj, _ = optimize_model(model)
            return obj if obj is not None else 0
        except:
            return 0


def _simulate_double_ko(model: Model, gene1: str, gene2: str) -> float:
    """
    Simulate double gene knockout.
    
    Args:
        model: Metabolic model
        gene1: First gene to knock out
        gene2: Second gene to knock out
    
    Returns:
        Growth rate relative to wild-type
    """
    with model:
        try:
            cobra.manipulation.knock_out_model_genes(model, [gene1, gene2])
            obj, _ = optimize_model(model)
            return obj if obj is not None else 0
        except:
            return 0


def _mock_single_ko(gene: str, wt_growth: float) -> float:
    """
    Mock single knockout for demo mode.
    
    Args:
        gene: Gene name
        wt_growth: Wild-type growth rate
    
    Returns:
        Mock growth rate
    """
    np.random.seed(hash(gene) % 2**32)
    
    # Essential genes have stronger effect
    essential_genes = ['HK2', 'LDHA', 'GLS', 'RRM2', 'TYMS']
    
    if gene in essential_genes:
        return np.random.uniform(0.1, 0.4)
    else:
        return np.random.uniform(0.5, 0.9)


def _mock_double_ko(gene1: str, gene2: str, wt_growth: float) -> float:
    """
    Mock double knockout for demo mode.
    
    Args:
        gene1: First gene
        gene2: Second gene
        wt_growth: Wild-type growth rate
    
    Returns:
        Mock growth rate
    """
    # Synthetic lethal pairs (known from literature)
    synthetic_lethal_pairs = [
        ('HK2', 'LDHA'),
        ('GLS', 'GLUD1'),
        ('RRM2', 'TYMS'),
        ('PHGDH', 'PSAT1'),
        ('ACLY', 'FASN')
    ]
    
    pair = tuple(sorted([gene1, gene2]))
    
    for sl_pair in synthetic_lethal_pairs:
        if pair == tuple(sorted(sl_pair)):
            return np.random.uniform(0.0, 0.1)  # Synthetic lethal
    
    # Otherwise, use Bliss independence with some noise
    growth1 = _mock_single_ko(gene1, wt_growth)
    growth2 = _mock_single_ko(gene2, wt_growth)
    
    # Add some synergy/antagonism noise
    noise = np.random.uniform(-0.2, 0.2)
    return max(0, min(1, growth1 * growth2 + noise))


def classify_interaction(synergy_score: float, growth_rate: float) -> str:
    """
    Classify genetic interaction type.
    
    Args:
        synergy_score: Synergy score (negative = synergistic)
        growth_rate: Double KO growth rate
    
    Returns:
        Interaction type classification
    """
    if growth_rate < 0.1 and synergy_score < -0.3:
        return "SYNTHETIC_LETHAL"
    elif synergy_score < -0.2:
        return "SYNERGISTIC"
    elif synergy_score > 0.2:
        return "ANTAGONISTIC"
    else:
        return "ADDITIVE"


def plot_synthetic_lethality(results_df: pd.DataFrame) -> 'go.Figure':
    """
    Create synthetic lethality visualization.
    
    Args:
        results_df: Synthetic lethality results
    
    Returns:
        Plotly Figure object
    """
    from plotly.subplots import make_subplots
    
    if go is None or make_subplots is None:
        raise ImportError("plotly is required for visualization")
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Synergy Score Distribution", "Interaction Type Breakdown"),
        horizontal_spacing=0.15
    )
    
    # Panel A: Histogram of synergy scores
    colors = []
    for idx, row in results_df.iterrows():
        if row['interaction_type'] == 'SYNTHETIC_LETHAL':
            colors.append('#E94F37')
        elif row['interaction_type'] == 'SYNERGISTIC':
            colors.append('#F18F01')
        elif row['interaction_type'] == 'ANTAGONISTIC':
            colors.append('#2E86AB')
        else:
            colors.append('#CCCCCC')
    
    fig.add_trace(
        go.Histogram(
            x=results_df['synergy_score'],
            marker_color=colors,
            name='Combinations',
            nbinsx=50
        ),
        row=1, col=1
    )
    
    # Panel B: Pie chart of interaction types
    type_counts = results_df['interaction_type'].value_counts()
    
    type_colors = {
        'SYNTHETIC_LETHAL': '#E94F37',
        'SYNERGISTIC': '#F18F01',
        'ADDITIVE': '#06A77D',
        'ANTAGONISTIC': '#2E86AB'
    }
    
    fig.add_trace(
        go.Pie(
            labels=type_counts.index,
            values=type_counts.values,
            marker_colors=[type_colors.get(t, '#CCCCCC') for t in type_counts.index],
            name='Types'
        ),
        row=1, col=2
    )
    
    # Update layout
    fig.update_layout(
        height=500,
        width=1200,
        title_text="Synthetic Lethality Screen Results",
        title_font_size=18,
        showlegend=False,
        template="plotly_white"
    )
    
    fig.update_xaxes(title_text="Synergy Score", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Synthetic Lethality Screen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python synthetic_lethality.py --expr data/expression.csv
  python synthetic_lethality.py --expr data/expression.csv --genes RRM2 GLS LDHA HK2
  python synthetic_lethality.py --expr data/expression.csv --model Human-GEM.xml --max-combos 200
        """
    )
    
    parser.add_argument('--expr', required=True, help='Path to expression data CSV')
    parser.add_argument('--model', help='Path to SBML model (optional)')
    parser.add_argument('--genes', nargs='+', help='Genes to screen (optional)')
    parser.add_argument('--high', type=float, default=75, help='High expression percentile (default: 75)')
    parser.add_argument('--low', type=float, default=25, help='Low expression percentile (default: 25)')
    parser.add_argument('--max-combos', type=int, default=100, help='Maximum combinations to test (default: 100)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_synthetic_lethality(
        expr_path=args.expr,
        gene_list=args.genes,
        model_path=args.model,
        high_percentile=args.high,
        low_percentile=args.low,
        max_combinations=args.max_combos,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
