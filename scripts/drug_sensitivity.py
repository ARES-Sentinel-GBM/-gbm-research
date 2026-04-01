# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Drug Sensitivity Prediction

Predicts drug sensitivity based on metabolic profile and drug-target interactions.
Integrates GDSC/CTRP drug response data with metabolic model predictions.
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


def run_drug_sensitivity(
    expr_path: str,
    model_path: Optional[str] = None,
    drug_targets_path: Optional[str] = None,
    ic50_path: Optional[str] = None,
    high_percentile: float = 75,
    low_percentile: float = 25,
    output_dir: str = "results",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Predict drug sensitivity based on metabolic profile.
    
    Args:
        expr_path: Path to gene expression CSV
        model_path: Path to SBML model (optional)
        drug_targets_path: Path to drug-target database (optional)
        ic50_path: Path to IC50 data (optional)
        high_percentile: Percentile for high expression threshold
        low_percentile: Percentile for low expression threshold
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        DataFrame with drug sensitivity predictions
    """
    print("=" * 60)
    print("ARES-GBM Drug Sensitivity Prediction")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Load expression data
    if verbose:
        print("Loading expression data...")
    expr_df = load_expression_matrix(expr_path)
    
    # Use mean across samples
    if expr_df.shape[1] > 1:
        expr_series = expr_df.mean(axis=1)
    else:
        expr_series = expr_df.iloc[:, 0]
    
    if verbose:
        print(f"  Genes: {len(expr_series)}")
        print()
    
    # Load or create metabolic model
    if model_path and os.path.exists(model_path):
        if verbose:
            print(f"Loading model from: {model_path}")
        model = read_sbml_model(model_path)
    else:
        if verbose:
            print("No model provided. Using demo mode with mock data.")
        return _create_mock_drug_results(output_dir)
    
    # Prepare model
    model = prepare_model(model)
    
    # Apply iMAT for GBM
    if verbose:
        print("Applying iMAT for GBM condition...")
    model_gbm = apply_imat(
        model,
        expr_series,
        high_percentile=high_percentile,
        low_percentile=low_percentile,
        verbose=False
    )
    
    # Get wild-type solution
    wt_obj, wt_fluxes = optimize_model(model)
    gbm_obj, gbm_fluxes = optimize_model(model_gbm)
    
    if verbose:
        print(f"  WT objective: {wt_obj:.4f}")
        print(f"  GBM objective: {gbm_obj:.4f}")
        print()
    
    # Load drug-target database
    if drug_targets_path and os.path.exists(drug_targets_path):
        drug_targets = pd.read_csv(drug_targets_path)
    else:
        # Use built-in drug-target database
        drug_targets = _get_builtin_drug_targets()
        if verbose:
            print(f"Using built-in drug-target database ({len(drug_targets)} drugs)")
    
    # Load IC50 data (optional)
    if ic50_path and os.path.exists(ic50_path):
        ic50_data = pd.read_csv(ic50_path)
        has_ic50 = True
    else:
        ic50_data = None
        has_ic50 = False
    
    # Predict drug sensitivity for each drug
    results = []
    
    if verbose:
        print(f"\nPredicting sensitivity for {len(drug_targets)} drugs...")
    
    for idx, drug_row in drug_targets.iterrows():
        drug_name = drug_row['drug_name']
        targets = drug_row['targets'].split(';') if isinstance(drug_row['targets'], str) else [drug_row['targets']]
        mechanism = drug_row.get('mechanism', 'Unknown')
        clinical_status = drug_row.get('clinical_status', 'Preclinical')
        
        # Calculate target expression score
        target_expr = []
        for target in targets:
            target = target.strip()
            if target in expr_series.index:
                target_expr.append(expr_series[target])
            elif target.upper() in expr_series.index:
                target_expr.append(expr_series[target.upper()])
        
        if target_expr:
            mean_target_expr = np.mean(target_expr)
            expr_percentile = (expr_series < mean_target_expr).mean() * 100
        else:
            mean_target_expr = 0
            expr_percentile = 50
        
        # Calculate metabolic vulnerability score
        # Based on flux changes when targets are inhibited
        vulnerability_score = _calculate_vulnerability(model, model_gbm, targets, wt_fluxes, gbm_fluxes)
        
        # Predict sensitivity score (0-100)
        # Higher score = more sensitive
        sensitivity_score = _calculate_sensitivity_score(
            mean_target_expr,
            expr_percentile,
            vulnerability_score,
            clinical_status
        )
        
        # Get IC50 if available
        if has_ic50 and ic50_data is not None:
            drug_ic50 = ic50_data[ic50_data['drug_name'] == drug_name]['IC50_uM']
            ic50_value = drug_ic50.values[0] if len(drug_ic50) > 0 else None
        else:
            ic50_value = None
        
        results.append({
            'drug_name': drug_name,
            'targets': ';'.join(targets),
            'mechanism': mechanism,
            'clinical_status': clinical_status,
            'target_expression': round(mean_target_expr, 3) if mean_target_expr else 0,
            'expression_percentile': round(expr_percentile, 2),
            'vulnerability_score': round(vulnerability_score, 3),
            'sensitivity_score': round(sensitivity_score, 2),
            'predicted_ic50': round(ic50_value, 3) if ic50_value else None,
            'sensitivity_class': classify_sensitivity(sensitivity_score)
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('sensitivity_score', ascending=False)
    
    # Save results
    ensure_dir(os.path.join(output_dir, 'drug_sensitivity.csv'))
    output_path = os.path.join(output_dir, 'drug_sensitivity.csv')
    results_df.to_csv(output_path, index=False)
    
    if verbose:
        print(f"\nResults saved to: {output_path}")
        print(f"\nTop 10 Predicted Drugs:")
        print("-" * 80)
        print(results_df[['drug_name', 'targets', 'sensitivity_score', 'sensitivity_class']].head(10).to_string())
        print()
    
    # Generate visualization
    try:
        fig = plot_drug_sensitivity(results_df)
        plot_path = os.path.join(output_dir, 'drug_sensitivity_results.png')
        save_plot(fig, plot_path)
        if verbose:
            print(f"Plot saved to: {plot_path}")
    except Exception as e:
        if verbose:
            print(f"Warning: Could not generate plot: {e}")
    
    # Save summary
    summary = {
        'total_drugs': len(results_df),
        'highly_sensitive': len(results_df[results_df['sensitivity_class'] == 'HIGH']),
        'moderately_sensitive': len(results_df[results_df['sensitivity_class'] == 'MODERATE']),
        'resistant': len(results_df[results_df['sensitivity_class'] == 'RESISTANT']),
        'top_drug': results_df.iloc[0]['drug_name'] if len(results_df) > 0 else None,
        'top_score': float(results_df.iloc[0]['sensitivity_score']) if len(results_df) > 0 else None
    }
    
    summary_path = os.path.join(output_dir, 'drug_sensitivity_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return results_df


def _calculate_vulnerability(
    model: Model,
    model_gbm: Model,
    targets: List[str],
    wt_fluxes: Dict,
    gbm_fluxes: Dict
) -> float:
    """
    Calculate metabolic vulnerability score based on target inhibition.
    
    Args:
        model: Wild-type model
        model_gbm: GBM-specific model
        targets: List of target genes
        wt_fluxes: Wild-type fluxes
        gbm_fluxes: GBM fluxes
    
    Returns:
        Vulnerability score (0-1)
    """
    # Simple vulnerability based on flux changes
    # More sophisticated: simulate target knockout
    
    target_reactions = []
    for target in targets:
        # Find reactions associated with target
        for rxn in model.reactions:
            if target in rxn.gene_reaction_rule:
                target_reactions.append(rxn.id)
    
    if not target_reactions:
        return 0.5  # Default score
    
    # Calculate flux difference for target reactions
    flux_diffs = []
    for rxn_id in target_reactions:
        if rxn_id in gbm_fluxes and rxn_id in wt_fluxes:
            diff = abs(gbm_fluxes[rxn_id] - wt_fluxes[rxn_id])
            flux_diffs.append(diff)
    
    if flux_diffs:
        # Normalize vulnerability score
        mean_diff = np.mean(flux_diffs)
        vulnerability = min(1.0, mean_diff / 10.0)  # Normalize to 0-1
        return vulnerability
    
    return 0.5


def _calculate_sensitivity_score(
    target_expression: float,
    expr_percentile: float,
    vulnerability_score: float,
    clinical_status: str
) -> float:
    """
    Calculate overall drug sensitivity score.
    
    Args:
        target_expression: Mean target gene expression
        expr_percentile: Expression percentile (0-100)
        vulnerability_score: Metabolic vulnerability (0-1)
        clinical_status: Clinical development status
    
    Returns:
        Sensitivity score (0-100)
    """
    # Base score from expression and vulnerability
    base_score = (
        0.4 * (100 - expr_percentile) +  # High expression → high sensitivity
        0.6 * (vulnerability_score * 100)  # High vulnerability → high sensitivity
    )
    
    # Boost for clinically validated drugs
    clinical_boost = {
        'FDA-approved': 1.2,
        'Clinical Trial': 1.1,
        'Preclinical': 1.0,
        'Investigational': 1.05
    }
    
    boost = clinical_boost.get(clinical_status, 1.0)
    final_score = min(100, base_score * boost)
    
    return final_score


def classify_sensitivity(score: float) -> str:
    """
    Classify drug sensitivity based on score.
    
    Args:
        score: Sensitivity score (0-100)
    
    Returns:
        Sensitivity classification
    """
    if score >= 70:
        return "HIGH"
    elif score >= 40:
        return "MODERATE"
    else:
        return "RESISTANT"


def _get_builtin_drug_targets() -> pd.DataFrame:
    """
    Get built-in drug-target database for GBM.
    
    Returns:
        DataFrame with drug-target information
    """
    drugs = [
        # Metabolic inhibitors
        {'drug_name': 'Temozolomide', 'targets': 'MGMT', 'mechanism': 'DNA alkylating agent', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'CB-839 (Telaglenastat)', 'targets': 'GLS', 'mechanism': 'Glutaminase inhibitor', 'clinical_status': 'Clinical Trial'},
        {'drug_name': 'LDHA-IN-1', 'targets': 'LDHA', 'mechanism': 'Lactate dehydrogenase A inhibitor', 'clinical_status': 'Preclinical'},
        {'drug_name': 'BPTES', 'targets': 'GLS', 'mechanism': 'Glutaminase allosteric inhibitor', 'clinical_status': 'Preclinical'},
        {'drug_name': '2-Deoxyglucose', 'targets': 'HK2;PFKP', 'mechanism': 'Glycolysis inhibitor', 'clinical_status': 'Clinical Trial'},
        {'drug_name': 'Lonidamine', 'targets': 'HK2', 'mechanism': 'Hexokinase inhibitor', 'clinical_status': 'Clinical Trial'},
        
        # Nucleotide synthesis inhibitors
        {'drug_name': 'Hydroxyurea', 'targets': 'RRM2', 'mechanism': 'Ribonucleotide reductase inhibitor', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'Triapine', 'targets': 'RRM2', 'mechanism': 'Ribonucleotide reductase inhibitor', 'clinical_status': 'Clinical Trial'},
        {'drug_name': 'Methotrexate', 'targets': 'DHFR;MTHFD2', 'mechanism': 'Antifolate', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'Pemetrexed', 'targets': 'DHFR;TYMS', 'mechanism': 'Antifolate', 'clinical_status': 'FDA-approved'},
        
        # Lipid metabolism inhibitors
        {'drug_name': 'TVB-2640', 'targets': 'FASN', 'mechanism': 'Fatty acid synthase inhibitor', 'clinical_status': 'Clinical Trial'},
        {'drug_name': 'Orlistat', 'targets': 'FASN', 'mechanism': 'Lipase inhibitor', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'Bemiparin', 'targets': 'ACLY', 'mechanism': 'ATP citrate lyase inhibitor', 'clinical_status': 'Clinical Trial'},
        
        # Serine/Glycine pathway
        {'drug_name': 'PH-99C', 'targets': 'PHGDH', 'mechanism': 'PHGDH inhibitor', 'clinical_status': 'Preclinical'},
        {'drug_name': 'NCT-503', 'targets': 'PHGDH', 'mechanism': 'PHGDH allosteric inhibitor', 'clinical_status': 'Preclinical'},
        
        # Mitochondrial inhibitors
        {'drug_name': 'Metformin', 'targets': 'NDUFV1', 'mechanism': 'Complex I inhibitor', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'Atovaquone', 'targets': 'UQCRB', 'mechanism': 'Complex III inhibitor', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'IACS-010759', 'targets': 'NDUFS1', 'mechanism': 'Complex I inhibitor', 'clinical_status': 'Clinical Trial'},
        
        # IDH inhibitors
        {'drug_name': 'Ivosidenib', 'targets': 'IDH1', 'mechanism': 'IDH1 inhibitor', 'clinical_status': 'FDA-approved'},
        {'drug_name': 'Vorasidenib', 'targets': 'IDH1;IDH2', 'mechanism': 'IDH1/2 inhibitor', 'clinical_status': 'Clinical Trial'},
        
        # MTHFD2 inhibitors
        {'drug_name': 'DS18561882', 'targets': 'MTHFD2', 'mechanism': 'MTHFD2 inhibitor', 'clinical_status': 'Preclinical'},
        {'drug_name': 'AGF347', 'targets': 'MTHFD2', 'mechanism': 'MTHFD2 inhibitor', 'clinical_status': 'Preclinical'},
    ]
    
    return pd.DataFrame(drugs)


def plot_drug_sensitivity(results_df: pd.DataFrame) -> 'go.Figure':
    """
    Create drug sensitivity visualization.
    
    Args:
        results_df: Drug sensitivity results
    
    Returns:
        Plotly Figure object
    """
    if go is None or make_subplots is None:
        raise ImportError("plotly is required for visualization")
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Drug Sensitivity Scores", "Sensitivity Distribution"),
        horizontal_spacing=0.15
    )
    
    # Sort by sensitivity score
    df_sorted = results_df.sort_values('sensitivity_score', ascending=True)
    
    # Color by sensitivity class
    colors = {
        'HIGH': '#E94F37',
        'MODERATE': '#F18F01',
        'RESISTANT': '#2E86AB'
    }
    
    bar_colors = [colors.get(cls, '#CCCCCC') for cls in df_sorted['sensitivity_class']]
    
    # Panel A: Bar chart of sensitivity scores
    fig.add_trace(
        go.Bar(
            x=df_sorted['sensitivity_score'],
            y=df_sorted['drug_name'],
            orientation='h',
            marker_color=bar_colors,
            name='Drugs',
            hovertemplate='<b>%{y}</b><br>Score: %{x:.1f}<br>Class: %{customdata}<extra></extra>',
            customdata=df_sorted['sensitivity_class']
        ),
        row=1, col=1
    )
    
    # Panel B: Distribution histogram
    for cls in ['HIGH', 'MODERATE', 'RESISTANT']:
        cls_data = results_df[results_df['sensitivity_class'] == cls]
        fig.add_trace(
            go.Histogram(
                x=cls_data['sensitivity_score'],
                name=f'{cls} ({len(cls_data)})',
                marker_color=colors[cls],
                opacity=0.7
            ),
            row=1, col=2
        )
    
    # Update layout
    fig.update_layout(
        height=600,
        width=1200,
        title_text="Drug Sensitivity Prediction Results",
        title_font_size=18,
        showlegend=True,
        template="plotly_white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    fig.update_xaxes(title_text="Sensitivity Score (0-100)", row=1, col=1)
    fig.update_xaxes(title_text="Sensitivity Score", row=1, col=2)
    fig.update_yaxes(title_text="Drug", row=1, col=1)
    
    return fig


def _create_mock_drug_results(output_dir: str) -> pd.DataFrame:
    """
    Create mock drug sensitivity results for demonstration.
    
    Args:
        output_dir: Directory for output files
    
    Returns:
        Mock results DataFrame
    """
    np.random.seed(42)
    
    drugs = _get_builtin_drug_targets()
    
    results = []
    for idx, drug_row in drugs.iterrows():
        sensitivity_score = np.random.uniform(30, 95)
        results.append({
            'drug_name': drug_row['drug_name'],
            'targets': drug_row['targets'],
            'mechanism': drug_row.get('mechanism', 'Unknown'),
            'clinical_status': drug_row.get('clinical_status', 'Preclinical'),
            'target_expression': np.random.uniform(5, 15),
            'expression_percentile': np.random.uniform(20, 90),
            'vulnerability_score': np.random.uniform(0.3, 0.9),
            'sensitivity_score': round(sensitivity_score, 2),
            'predicted_ic50': np.random.uniform(0.1, 10),
            'sensitivity_class': classify_sensitivity(sensitivity_score)
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('sensitivity_score', ascending=False)
    
    # Save mock results
    ensure_dir(os.path.join(output_dir, 'drug_sensitivity_demo.csv'))
    output_path = os.path.join(output_dir, 'drug_sensitivity_demo.csv')
    results_df.to_csv(output_path, index=False)
    print(f"Demo results saved to: {output_path}")
    
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Drug Sensitivity Prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python drug_sensitivity.py --expr data/gbm.csv
  python drug_sensitivity.py --expr data/gbm.csv --drugs data/drug_targets.csv
  python drug_sensitivity.py --expr data/gbm.csv --model Human-GEM.xml --output results/drug_results
        """
    )
    
    parser.add_argument('--expr', required=True, help='Path to expression data CSV')
    parser.add_argument('--model', help='Path to SBML model (optional)')
    parser.add_argument('--drugs', help='Path to drug-target database CSV (optional)')
    parser.add_argument('--ic50', help='Path to IC50 data CSV (optional)')
    parser.add_argument('--high', type=float, default=75, help='High expression percentile (default: 75)')
    parser.add_argument('--low', type=float, default=25, help='Low expression percentile (default: 25)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_drug_sensitivity(
        expr_path=args.expr,
        model_path=args.model,
        drug_targets_path=args.drugs,
        ic50_path=args.ic50,
        high_percentile=args.high,
        low_percentile=args.low,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
