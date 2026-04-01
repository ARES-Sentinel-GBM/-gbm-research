# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — GBM Subtype Analysis

Analyzes metabolic differences across GBM molecular subtypes:
- Proneural
- Neural
- Classical
- Mesenchymal

Based on TCGA gene expression signatures.
"""
import os
import sys
import argparse
import json
from datetime import datetime
from typing import Dict, List, Tuple, Optional

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


# TCGA-GBM subtype classifier genes (Verhaak et al., 2010)
SUBTYPE_SIGNATURES = {
    'Proneural': [
        'DLL3', 'ASCL1', 'OLIG2', 'SOX2', 'NKX2-2', 'TCF3', 'IDH1', 'PDGFRA'
    ],
    'Neural': [
        'NEFL', 'GABRA1', 'SLC12A5', 'GABRG2', 'SYT1', 'HOMER1', 'PCP4'
    ],
    'Classical': [
        'EGFR', 'ERBB2', 'STAT3', 'NKX2-1', 'SOX10', 'GPM6A', 'OLIG1'
    ],
    'Mesenchymal': [
        'NF1', 'RELB', 'CD44', 'MET', 'CLIC3', 'FAP', 'TNC', 'VEGFA'
    ]
}


def run_subtype_analysis(
    expr_path: str,
    model_path: Optional[str] = None,
    high_percentile: float = 75,
    low_percentile: float = 25,
    output_dir: str = "results",
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze metabolic differences across GBM molecular subtypes.
    
    Args:
        expr_path: Path to gene expression CSV (samples in columns)
        model_path: Path to SBML model (optional)
        high_percentile: Percentile for high expression threshold
        low_percentile: Percentile for low expression threshold
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        Tuple of (subtype_assignments, differential_results)
    """
    print("=" * 60)
    print("ARES-GBM Subtype Analysis")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load expression data
    if verbose:
        print("Loading expression data...")
    expr_df = load_expression_matrix(expr_path)
    
    if verbose:
        print(f"  Samples: {expr_df.shape[1]}")
        print(f"  Genes: {expr_df.shape[0]}")
        print()
    
    # Classify samples into subtypes
    if verbose:
        print("Classifying samples into molecular subtypes...")
    subtype_assignments = classify_subtypes(expr_df)
    
    if verbose:
        print(f"\nSubtype distribution:")
        print(subtype_assignments['subtype'].value_counts().to_string())
        print()
    
    # Build subtype-specific metabolic models
    if verbose:
        print("Building subtype-specific metabolic models...")
    
    # Load or create metabolic model
    if model_path and os.path.exists(model_path):
        model = read_sbml_model(model_path)
        model = prepare_model(model)
    else:
        if verbose:
            print("No model provided. Using demo mode with mock data.")
        model = None
    
    # Analyze each subtype
    subtype_results = {}
    all_differential = []
    
    for subtype in ['Proneural', 'Neural', 'Classical', 'Mesenchymal']:
        subtype_samples = subtype_assignments[subtype_assignments['subtype'] == subtype].index.tolist()
        
        if len(subtype_samples) == 0:
            if verbose:
                print(f"  {subtype}: No samples assigned (skipping)")
            continue
        
        if verbose:
            print(f"\n  {subtype}: {len(subtype_samples)} samples")
        
        # Get mean expression for subtype
        subtype_expr = expr_df[subtype_samples].mean(axis=1)
        
        # Build context-specific model if model available
        if model is not None:
            subtype_model = apply_imat(
                model,
                subtype_expr,
                high_percentile=high_percentile,
                low_percentile=low_percentile,
                verbose=False
            )
            obj, fluxes = optimize_model(subtype_model)
        else:
            fluxes = _create_mock_fluxes(subtype_expr)
            obj = 1.0
        
        subtype_results[subtype] = {
            'samples': subtype_samples,
            'n_samples': len(subtype_samples),
            'fluxes': fluxes,
            'objective': obj,
            'mean_expression': subtype_expr
        }
        
        if verbose:
            print(f"    Model objective: {obj:.4f}")
    
    # Differential analysis between subtypes
    if verbose:
        print("\nPerforming differential analysis...")
    
    if model is not None:
        all_differential = compare_subtypes_metabolic(subtype_results, model.reactions)
    else:
        all_differential = _create_mock_differential_results()
    
    # Convert to DataFrame
    diff_df = pd.DataFrame(all_differential)
    diff_df = diff_df.sort_values('max_fold_change', ascending=False)
    
    # Save results
    ensure_dir(os.path.join(output_dir, 'subtype_assignments.csv'))
    subtype_assignments.to_csv(os.path.join(output_dir, 'subtype_assignments.csv'))
    
    ensure_dir(os.path.join(output_dir, 'subtype_differential_analysis.csv'))
    diff_df.to_csv(os.path.join(output_dir, 'subtype_differential_analysis.csv'), index=False)
    
    if verbose:
        print(f"\nResults saved to:")
        print(f"  - subtype_assignments.csv")
        print(f"  - subtype_differential_analysis.csv")
    
    # Generate visualizations
    try:
        fig = plot_subtype_analysis(subtype_assignments, diff_df)
        plot_path = os.path.join(output_dir, 'subtype_analysis_results.png')
        save_plot(fig, plot_path)
        if verbose:
            print(f"Plot saved to: {plot_path}")
    except Exception as e:
        if verbose:
            print(f"Warning: Could not generate plot: {e}")
    
    # Save summary
    summary = {
        'total_samples': len(subtype_assignments),
        'subtype_counts': subtype_assignments['subtype'].value_counts().to_dict(),
        'differential_reactions': len(diff_df),
        'top_differential': diff_df.iloc[0]['reaction'] if len(diff_df) > 0 else None
    }
    
    summary_path = os.path.join(output_dir, 'subtype_analysis_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return subtype_assignments, diff_df


def classify_subtypes(expr_df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify samples into GBM molecular subtypes based on signature genes.
    
    Uses correlation-based classification with TCGA subtype signatures.
    
    Args:
        expr_df: Expression DataFrame (genes x samples)
    
    Returns:
        DataFrame with sample assignments
    """
    assignments = []
    
    # Normalize expression (z-score)
    means = expr_df.mean(axis=1)
    stds = expr_df.std(axis=1)
    expr_normalized = (expr_df.sub(means, axis=0)) / (stds + 1e-6)
    
    for sample in expr_df.columns:
        sample_expr = expr_normalized[sample]
        
        # Calculate correlation with each subtype signature
        correlations = {}
        for subtype, signature_genes in SUBTYPE_SIGNATURES.items():
            # Find overlapping genes
            overlapping = [g for g in signature_genes if g in sample_expr.index or g.upper() in sample_expr.index]
            
            if len(overlapping) == 0:
                correlations[subtype] = 0
                continue
            
            # Get expression of signature genes
            subtype_expr = []
            for gene in overlapping:
                if gene in sample_expr.index:
                    subtype_expr.append(sample_expr[gene])
                elif gene.upper() in sample_expr.index:
                    subtype_expr.append(sample_expr[gene.upper()])
            
            # Mean expression of signature (simple scoring)
            correlations[subtype] = np.mean(subtype_expr)
        
        # Assign to subtype with highest correlation
        best_subtype = max(correlations, key=correlations.get)
        best_score = correlations[best_subtype]
        
        assignments.append({
            'sample': sample,
            'subtype': best_subtype,
            'score': round(best_score, 4),
            'proneural_score': round(correlations.get('Proneural', 0), 4),
            'neural_score': round(correlations.get('Neural', 0), 4),
            'classical_score': round(correlations.get('Classical', 0), 4),
            'mesenchymal_score': round(correlations.get('Mesenchymal', 0), 4)
        })
    
    assignments_df = pd.DataFrame(assignments)
    assignments_df = assignments_df.set_index('sample')
    
    return assignments_df


def compare_subtypes_metabolic(
    subtype_results: Dict,
    reactions: List
) -> List[Dict]:
    """
    Compare metabolic fluxes across subtypes.
    
    Args:
        subtype_results: Dictionary of subtype-specific results
        reactions: List of reactions to compare
    
    Returns:
        List of differential reaction results
    """
    results = []
    subtypes = list(subtype_results.keys())
    
    for rxn in reactions:
        rxn_id = rxn.id
        
        # Get flux for each subtype
        fluxes = {}
        for subtype in subtypes:
            fluxes[subtype] = subtype_results[subtype]['fluxes'].get(rxn_id, 0)
        
        # Calculate statistics
        flux_values = list(fluxes.values())
        mean_flux = np.mean(flux_values)
        std_flux = np.std(flux_values)
        max_flux = max(flux_values)
        min_flux = min(flux_values)
        
        # Fold change (max / min, with pseudocount)
        if abs(min_flux) > 1e-6:
            fold_change = max_flux / min_flux if min_flux > 0 else abs(min_flux / max_flux) if max_flux > 0 else 1
        else:
            fold_change = abs(max_flux) + 1
        
        # Find which subtype has highest flux
        highest_subtype = max(fluxes, key=lambda k: abs(fluxes[k]))
        
        results.append({
            'reaction': rxn_id,
            'mean_flux': round(mean_flux, 4),
            'std_flux': round(std_flux, 4),
            'max_fold_change': round(fold_change, 3),
            'highest_subtype': highest_subtype,
            **{f'{s}_flux': round(fluxes[s], 4) for s in subtypes}
        })
    
    return results


def _create_mock_fluxes(expr_series: pd.Series) -> Dict:
    """
    Create mock flux values for demonstration.
    
    Args:
        expr_series: Expression series
    
    Returns:
        Dictionary of mock fluxes
    """
    np.random.seed(hash(str(expr_series.name if hasattr(expr_series, 'name') else 'mock')) % 2**32)
    
    # Create fluxes based on expression
    base_fluxes = {}
    for i in range(100):  # Mock 100 reactions
        rxn_id = f'R_RXN{i:04d}'
        base_fluxes[rxn_id] = np.random.uniform(-10, 10) * (expr_series.mean() / 10)
    
    return base_fluxes


def _create_mock_differential_results() -> List[Dict]:
    """
    Create mock differential analysis results.
    
    Returns:
        List of differential results
    """
    np.random.seed(42)
    
    reactions = [
        'R_GK', 'R_PYK', 'R_LDH', 'R_G6PDH', 'R_GLS', 'R_GLS2',
        'R_IDH', 'R_MDH', 'R_SDH', 'R_FUM', 'R_ATPS', 'R_NDH',
        'R_PFK', 'R_ALDO', 'R_TPI', 'R_GAPD', 'R_PGK', 'R_PGM',
        'R_ENO', 'R_PPCK', 'R_PC', 'R_PEPCK', 'R_FBP', 'R_TKT'
    ]
    
    subtypes = ['Proneural', 'Neural', 'Classical', 'Mesenchymal']
    
    results = []
    for rxn in reactions:
        fluxes = {s: np.random.uniform(-10, 10) for s in subtypes}
        flux_values = list(fluxes.values())
        
        max_flux = max(flux_values)
        min_flux = min([f for f in flux_values if abs(f) > 0.1] or [1])
        fold_change = abs(max_flux / min_flux) if abs(min_flux) > 0.1 else 1
        
        highest_subtype = max(fluxes, key=lambda k: abs(fluxes[k]))
        
        results.append({
            'reaction': rxn,
            'mean_flux': round(np.mean(flux_values), 4),
            'std_flux': round(np.std(flux_values), 4),
            'max_fold_change': round(fold_change, 3),
            'highest_subtype': highest_subtype,
            **{f'{s}_flux': round(f, 4) for s, f in fluxes.items()}
        })
    
    return results


def plot_subtype_analysis(
    assignments_df: pd.DataFrame,
    diff_df: pd.DataFrame
) -> 'go.Figure':
    """
    Create subtype analysis visualization.
    
    Args:
        assignments_df: Subtype assignments
        diff_df: Differential analysis results
    
    Returns:
        Plotly Figure object
    """
    if go is None or make_subplots is None:
        raise ImportError("plotly is required for visualization")
    
    from plotly.subplots import make_subplots
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Subtype Distribution",
            "Subtype Signature Scores",
            "Top Differential Reactions",
            "Flux Distribution by Subtype"
        ),
        vertical_spacing=0.15,
        horizontal_spacing=0.12
    )
    
    # Panel A: Subtype distribution (pie chart)
    subtype_counts = assignments_df['subtype'].value_counts()
    colors = {
        'Proneural': '#2E86AB',
        'Neural': '#06A77D',
        'Classical': '#F18F01',
        'Mesenchymal': '#A23B72'
    }
    
    fig.add_trace(
        go.Pie(
            labels=subtype_counts.index,
            values=subtype_counts.values,
            marker_colors=[colors.get(s, '#CCCCCC') for s in subtype_counts.index],
            name='Distribution'
        ),
        row=1, col=1
    )
    
    # Panel B: Signature scores (violin plot)
    for subtype in ['Proneural', 'Neural', 'Classical', 'Mesenchymal']:
        score_col = f'{subtype.lower()}_score'
        if score_col in assignments_df.columns:
            fig.add_trace(
                go.Violin(
                    y=assignments_df[score_col],
                    name=subtype,
                    marker_color=colors.get(subtype, '#CCCCCC'),
                    box_visible=True,
                    meanline_visible=True,
                    showlegend=False
                ),
                row=1, col=2
            )
    
    # Panel C: Top differential reactions (bar chart)
    top_diff = diff_df.nlargest(15, 'max_fold_change')
    
    fig.add_trace(
        go.Bar(
            x=top_diff['max_fold_change'],
            y=top_diff['reaction'],
            orientation='h',
            marker_color=top_diff['highest_subtype'].map(colors),
            name='Reactions',
            hovertemplate='<b>%{y}</b><br>FC: %{x:.2f}<br>Highest: %{customdata}<extra></extra>',
            customdata=top_diff['highest_subtype']
        ),
        row=2, col=1
    )
    
    # Panel D: Flux distribution by subtype (box plot)
    for subtype in ['Proneural', 'Neural', 'Classical', 'Mesenchymal']:
        flux_col = f'{subtype}_flux'
        if flux_col in diff_df.columns:
            fig.add_trace(
                go.Box(
                    y=diff_df[flux_col],
                    name=subtype,
                    marker_color=colors.get(subtype, '#CCCCCC'),
                    showlegend=False
                ),
                row=2, col=2
            )
    
    # Update layout
    fig.update_layout(
        height=800,
        width=1200,
        title_text="GBM Subtype Analysis Results",
        title_font_size=18,
        showlegend=True,
        template="plotly_white"
    )
    
    fig.update_xaxes(title_text="Count", row=1, col=1)
    fig.update_xaxes(title_text="Signature Score", row=1, col=2)
    fig.update_xaxes(title_text="Fold Change", row=2, col=1)
    fig.update_xaxes(title_text="Flux", row=2, col=2)
    
    fig.update_yaxes(title_text="Subtype", row=1, col=1)
    fig.update_yaxes(title_text="Reaction", row=2, col=1)
    
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Subtype Analysis - Molecular subtype classification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python subtype_analysis.py --expr data/gbm_expression.csv
  python subtype_analysis.py --expr data/gbm.csv --model Human-GEM.xml
  python subtype_analysis.py --expr data/gbm.csv --output results/subtype
        """
    )
    
    parser.add_argument('--expr', required=True, help='Path to expression data CSV')
    parser.add_argument('--model', help='Path to SBML model (optional)')
    parser.add_argument('--high', type=float, default=75, help='High expression percentile (default: 75)')
    parser.add_argument('--low', type=float, default=25, help='Low expression percentile (default: 25)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_subtype_analysis(
        expr_path=args.expr,
        model_path=args.model,
        high_percentile=args.high,
        low_percentile=args.low,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
