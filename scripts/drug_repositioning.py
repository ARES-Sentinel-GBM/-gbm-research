# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Drug Repositioning Screen

Identifies existing FDA-approved drugs with potential anti-GBM activity
through metabolic target analysis.
"""
import os
import sys
import argparse
import json
from datetime import datetime
from typing import Dict, List, Optional

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


def run_drug_repositioning(
    expr_path: str,
    model_path: Optional[str] = None,
    high_percentile: float = 75,
    low_percentile: float = 25,
    output_dir: str = "results",
    verbose: bool = True
) -> pd.DataFrame:
    """
    Screen FDA-approved drugs for potential GBM repositioning.
    
    Args:
        expr_path: Path to gene expression CSV
        model_path: Path to SBML model (optional)
        high_percentile: Percentile for high expression threshold
        low_percentile: Percentile for low expression threshold
        output_dir: Directory for output files
        verbose: Print progress information
    
    Returns:
        DataFrame with drug repositioning predictions
    """
    print("=" * 60)
    print("ARES-GBM Drug Repositioning Screen")
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
    
    # Load FDA-approved drug database
    if verbose:
        print("Loading FDA-approved drug database...")
    drug_db = _get_fda_approved_drugs()
    
    if verbose:
        print(f"  Total drugs: {len(drug_db)}")
        print()
    
    # Load or create metabolic model
    if model_path and os.path.exists(model_path):
        if verbose:
            print(f"Loading model from: {model_path}")
        model = read_sbml_model(model_path)
        model = prepare_model(model)
        
        if verbose:
            print("Building context-specific GBM model...")
        model_gbm = apply_imat(
            model,
            expr_series,
            high_percentile=high_percentile,
            low_percentile=low_percentile,
            verbose=False
        )
        
        wt_obj, _ = optimize_model(model)
        gbm_obj, _ = optimize_model(model_gbm)
    else:
        if verbose:
            print("No model provided. Using demo mode with mock data.")
        model = None
        model_gbm = None
        wt_obj = 1.0
        gbm_obj = 1.0
    
    # Score each drug
    results = []
    
    if verbose:
        print(f"\nScoring {len(drug_db)} drugs...")
    
    for idx, drug_row in drug_db.iterrows():
        drug_name = drug_row['drug_name']
        original_indication = drug_row['original_indication']
        targets = drug_row['targets'].split(';') if isinstance(drug_row['targets'], str) else [drug_row['targets']]
        mechanism = drug_row.get('mechanism', 'Unknown')
        safety_score = drug_row.get('safety_score', 0.5)
        bbb_penetration = drug_row.get('bbb_penetration', 'Unknown')
        
        # Calculate target expression score
        target_expr_score = _calculate_target_expression(expr_series, targets)
        
        # Calculate metabolic vulnerability
        if model is not None:
            vulnerability = _calculate_metabolic_vulnerability(model_gbm, targets)
        else:
            vulnerability = _mock_vulnerability(targets)
        
        # Calculate repositioning score
        # Higher = better candidate
        repositioning_score = _calculate_repositioning_score(
            target_expr_score,
            vulnerability,
            safety_score,
            bbb_penetration
        )
        
        # Priority classification
        priority = classify_priority(repositioning_score, bbb_penetration)
        
        results.append({
            'drug_name': drug_name,
            'original_indication': original_indication,
            'targets': ';'.join(targets),
            'mechanism': mechanism,
            'target_expression_score': round(target_expr_score, 3),
            'metabolic_vulnerability': round(vulnerability, 3),
            'safety_score': round(safety_score, 3),
            'bbb_penetration': bbb_penetration,
            'repositioning_score': round(repositioning_score, 2),
            'priority': priority,
            'gbm_relevance': _get_gbm_relevance(drug_name, targets)
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('repositioning_score', ascending=False)
    
    # Save results
    ensure_dir(os.path.join(output_dir, 'drug_repositioning.csv'))
    output_path = os.path.join(output_dir, 'drug_repositioning.csv')
    results_df.to_csv(output_path, index=False)
    
    # Summary statistics (calculate before verbose block)
    n_high = len(results_df[results_df['priority'] == 'HIGH'])
    n_medium = len(results_df[results_df['priority'] == 'MEDIUM'])
    n_low = len(results_df[results_df['priority'] == 'LOW'])
    
    if verbose:
        print(f"\nResults saved to: {output_path}")
        print(f"\nPriority Distribution:")
        print(f"  HIGH: {n_high}")
        print(f"  MEDIUM: {n_medium}")
        print(f"  LOW: {n_low}")
        
        print(f"\nTop 10 Repositioning Candidates:")
        print("-" * 100)
        print(results_df[['drug_name', 'original_indication', 'repositioning_score', 'priority']].head(10).to_string())
        print()
    
    # Generate visualization
    try:
        fig = plot_drug_repositioning(results_df)
        plot_path = os.path.join(output_dir, 'drug_repositioning_results.png')
        save_plot(fig, plot_path)
        if verbose:
            print(f"Plot saved to: {plot_path}")
    except Exception as e:
        if verbose:
            print(f"Warning: Could not generate plot: {e}")
    
    # Save summary
    summary = {
        'total_drugs': len(results_df),
        'high_priority': int(n_high),
        'medium_priority': int(n_medium),
        'low_priority': int(n_low),
        'top_candidate': results_df.iloc[0]['drug_name'] if len(results_df) > 0 else None,
        'top_score': float(results_df.iloc[0]['repositioning_score']) if len(results_df) > 0 else None
    }
    
    summary_path = os.path.join(output_dir, 'drug_repositioning_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    return results_df


def _get_fda_approved_drugs() -> pd.DataFrame:
    """
    Get database of FDA-approved drugs for repositioning screen.
    
    Returns:
        DataFrame with drug information
    """
    drugs = [
        # Metabolic drugs
        {'drug_name': 'Metformin', 'original_indication': 'Type 2 Diabetes', 'targets': 'NDUFV1;AMPK', 'mechanism': 'Complex I inhibitor', 'safety_score': 0.95, 'bbb_penetration': 'High'},
        {'drug_name': 'Atovaquone', 'original_indication': 'Malaria', 'targets': 'UQCRB', 'mechanism': 'Complex III inhibitor', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Dichloroacetate', 'original_indication': 'Lactic Acidosis', 'targets': 'PDK1;PDK2', 'mechanism': 'PDK inhibitor', 'safety_score': 0.80, 'bbb_penetration': 'High'},
        {'drug_name': 'Orlistat', 'original_indication': 'Obesity', 'targets': 'FASN', 'mechanism': 'Fatty acid synthase inhibitor', 'safety_score': 0.85, 'bbb_penetration': 'Low'},
        {'drug_name': 'Leflunomide', 'original_indication': 'Rheumatoid Arthritis', 'targets': 'DHODH', 'mechanism': 'Pyrimidine synthesis inhibitor', 'safety_score': 0.75, 'bbb_penetration': 'Medium'},
        
        # Cardiovascular drugs
        {'drug_name': 'Aspirin', 'original_indication': 'Cardiovascular Disease', 'targets': 'PTGS1;PTGS2', 'mechanism': 'COX inhibitor', 'safety_score': 0.90, 'bbb_penetration': 'High'},
        {'drug_name': 'Simvastatin', 'original_indication': 'Hypercholesterolemia', 'targets': 'HMGCR', 'mechanism': 'HMG-CoA reductase inhibitor', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Captopril', 'original_indication': 'Hypertension', 'targets': 'ACE', 'mechanism': 'ACE inhibitor', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        {'drug_name': 'Digoxin', 'original_indication': 'Heart Failure', 'targets': 'ATP1A1', 'mechanism': 'Na+/K+ ATPase inhibitor', 'safety_score': 0.65, 'bbb_penetration': 'Low'},
        {'drug_name': 'Propranolol', 'original_indication': 'Hypertension', 'targets': 'ADRB1;ADRB2', 'mechanism': 'Beta blocker', 'safety_score': 0.85, 'bbb_penetration': 'High'},
        
        # CNS drugs
        {'drug_name': 'Sertraline', 'original_indication': 'Depression', 'targets': 'SLC6A4', 'mechanism': 'SSRI', 'safety_score': 0.85, 'bbb_penetration': 'High'},
        {'drug_name': 'Fluoxetine', 'original_indication': 'Depression', 'targets': 'SLC6A4', 'mechanism': 'SSRI', 'safety_score': 0.85, 'bbb_penetration': 'High'},
        {'drug_name': 'Amitriptyline', 'original_indication': 'Depression', 'targets': 'SLC6A2;SLC6A4', 'mechanism': 'TCA', 'safety_score': 0.75, 'bbb_penetration': 'High'},
        {'drug_name': 'Haloperidol', 'original_indication': 'Schizophrenia', 'targets': 'DRD2', 'mechanism': 'Antipsychotic', 'safety_score': 0.70, 'bbb_penetration': 'High'},
        {'drug_name': 'Valproic Acid', 'original_indication': 'Epilepsy', 'targets': 'HDAC', 'mechanism': 'HDAC inhibitor', 'safety_score': 0.80, 'bbb_penetration': 'High'},
        
        # Antibiotics
        {'drug_name': 'Doxycycline', 'original_indication': 'Bacterial Infection', 'targets': 'MMP2;MMP9', 'mechanism': 'Tetracycline antibiotic', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Azithromycin', 'original_indication': 'Bacterial Infection', 'targets': 'RPL4', 'mechanism': 'Macrolide antibiotic', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Ciprofloxacin', 'original_indication': 'Bacterial Infection', 'targets': 'TOP2A', 'mechanism': 'Fluoroquinolone', 'safety_score': 0.80, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Rifampicin', 'original_indication': 'Tuberculosis', 'targets': 'RNAP', 'mechanism': 'RNA polymerase inhibitor', 'safety_score': 0.75, 'bbb_penetration': 'High'},
        {'drug_name': 'Nitroxoline', 'original_indication': 'UTI', 'targets': 'CTSB', 'mechanism': 'Cathepsin B inhibitor', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
        
        # Anti-inflammatory
        {'drug_name': 'Ibuprofen', 'original_indication': 'Pain/Inflammation', 'targets': 'PTGS1;PTGS2', 'mechanism': 'NSAID', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Celecoxib', 'original_indication': 'Arthritis', 'targets': 'PTGS2', 'mechanism': 'COX-2 inhibitor', 'safety_score': 0.75, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Dexamethasone', 'original_indication': 'Inflammation', 'targets': 'NR3C1', 'mechanism': 'Corticosteroid', 'safety_score': 0.80, 'bbb_penetration': 'High'},
        {'drug_name': 'Prednisone', 'original_indication': 'Inflammation', 'targets': 'NR3C1', 'mechanism': 'Corticosteroid', 'safety_score': 0.80, 'bbb_penetration': 'Medium'},
        
        # Antiparasitic
        {'drug_name': 'Ivermectin', 'original_indication': 'Parasitic Infection', 'targets': 'PAX3', 'mechanism': 'Avermectin', 'safety_score': 0.80, 'bbb_penetration': 'High'},
        {'drug_name': 'Niclosamide', 'original_indication': 'Tapeworm', 'targets': 'STAT3', 'mechanism': 'Salicylanilide', 'safety_score': 0.85, 'bbb_penetration': 'Low'},
        {'drug_name': 'Pyrvinium', 'original_indication': 'Pinworm', 'targets': 'WNT', 'mechanism': 'Wnt pathway inhibitor', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        
        # Others
        {'drug_name': 'Disulfiram', 'original_indication': 'Alcohol Dependence', 'targets': 'ALDH2', 'mechanism': 'ALDH inhibitor', 'safety_score': 0.75, 'bbb_penetration': 'High'},
        {'drug_name': 'Thiabendazole', 'original_indication': 'Fungal Infection', 'targets': 'HIF1A', 'mechanism': 'HIF inhibitor', 'safety_score': 0.75, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Auranofin', 'original_indication': 'Rheumatoid Arthritis', 'targets': 'TXNRD1', 'mechanism': 'Thioredoxin reductase inhibitor', 'safety_score': 0.70, 'bbb_penetration': 'Medium'},
        {'drug_name': 'Ebselen', 'original_indication': 'Investigational', 'targets': 'GPX1', 'mechanism': 'Glutathione peroxidase mimetic', 'safety_score': 0.75, 'bbb_penetration': 'High'},
        {'drug_name': 'Pimozide', 'original_indication': 'Schizophrenia', 'targets': 'DRD2', 'mechanism': 'Antipsychotic', 'safety_score': 0.70, 'bbb_penetration': 'High'},
        {'drug_name': 'Penfluridol', 'original_indication': 'Schizophrenia', 'targets': 'DRD2', 'mechanism': 'Antipsychotic', 'safety_score': 0.70, 'bbb_penetration': 'High'},
        {'drug_name': 'Trifluoperazine', 'original_indication': 'Schizophrenia', 'targets': 'DRD2', 'mechanism': 'Antipsychotic', 'safety_score': 0.70, 'bbb_penetration': 'High'},
        {'drug_name': 'Thioridazine', 'original_indication': 'Schizophrenia', 'targets': 'DRD2', 'mechanism': 'Antipsychotic', 'safety_score': 0.65, 'bbb_penetration': 'High'},
        {'drug_name': 'Clotrimazole', 'original_indication': 'Fungal Infection', 'targets': 'HSP90', 'mechanism': 'Azole antifungal', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        {'drug_name': 'Miconazole', 'original_indication': 'Fungal Infection', 'targets': 'CYP51', 'mechanism': 'Azole antifungal', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        {'drug_name': 'Econazole', 'original_indication': 'Fungal Infection', 'targets': 'CYP51', 'mechanism': 'Azole antifungal', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        {'drug_name': 'Tioconazole', 'original_indication': 'Fungal Infection', 'targets': 'CYP51', 'mechanism': 'Azole antifungal', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        {'drug_name': 'Bifonazole', 'original_indication': 'Fungal Infection', 'targets': 'CYP51', 'mechanism': 'Azole antifungal', 'safety_score': 0.80, 'bbb_penetration': 'Low'},
        {'drug_name': 'Statin Combo', 'original_indication': 'Hypercholesterolemia', 'targets': 'HMGCR', 'mechanism': 'HMG-CoA reductase inhibitor', 'safety_score': 0.85, 'bbb_penetration': 'Medium'},
    ]
    
    return pd.DataFrame(drugs)


def _calculate_target_expression(expr_series: pd.Series, targets: List[str]) -> float:
    """
    Calculate target expression score.
    
    Args:
        expr_series: Gene expression series
        targets: List of target genes
    
    Returns:
        Expression score (0-1)
    """
    expr_values = []
    for target in targets:
        if target in expr_series.index:
            expr_values.append(expr_series[target])
        elif target.upper() in expr_series.index:
            expr_values.append(expr_series[target.upper()])
    
    if expr_values:
        mean_expr = np.mean(expr_values)
        # Normalize to percentile
        percentile = (expr_series < mean_expr).mean()
        return percentile
    
    return 0.5


def _calculate_metabolic_vulnerability(model: Model, targets: List[str]) -> float:
    """
    Calculate metabolic vulnerability score.
    
    Args:
        model: Metabolic model
        targets: List of target genes
    
    Returns:
        Vulnerability score (0-1)
    """
    try:
        with model:
            cobra.manipulation.knock_out_model_genes(model, targets)
            obj, _ = optimize_model(model)
            
            if obj is None or obj < 0.1:
                return 1.0  # High vulnerability
            elif obj < 0.5:
                return 0.7
            elif obj < 0.8:
                return 0.4
            else:
                return 0.2
    except:
        return 0.5


def _mock_vulnerability(targets: List[str]) -> float:
    """
    Mock vulnerability for demo mode.
    
    Args:
        targets: List of target genes
    
    Returns:
        Mock vulnerability score
    """
    np.random.seed(hash(''.join(sorted(targets))) % 2**32)
    
    # Known metabolic targets have higher vulnerability
    metabolic_targets = ['NDUFV1', 'UQCRB', 'FASN', 'DHODH', 'HMGCR', 'PDK1']
    
    for target in targets:
        if target in metabolic_targets:
            return np.random.uniform(0.6, 0.95)
    
    return np.random.uniform(0.2, 0.6)


def _calculate_repositioning_score(
    target_expr: float,
    vulnerability: float,
    safety: float,
    bbb: str
) -> float:
    """
    Calculate overall drug repositioning score.
    
    Args:
        target_expr: Target expression score
        vulnerability: Metabolic vulnerability score
        safety: Safety score
        bbb: BBB penetration (High/Medium/Low)
    
    Returns:
        Repositioning score (0-100)
    """
    # BBB penetration boost
    bbb_boost = {'High': 1.3, 'Medium': 1.0, 'Low': 0.7}
    
    base_score = (
        0.3 * target_expr * 100 +
        0.4 * vulnerability * 100 +
        0.3 * safety * 100
    )
    
    return base_score * bbb_boost.get(bbb, 1.0)


def classify_priority(score: float, bbb: str) -> str:
    """
    Classify drug priority.
    
    Args:
        score: Repositioning score
        bbb: BBB penetration
    
    Returns:
        Priority classification
    """
    if score >= 70 and bbb == 'High':
        return 'HIGH'
    elif score >= 50:
        return 'MEDIUM'
    else:
        return 'LOW'


def _get_gbm_relevance(drug_name: str, targets: List[str]) -> str:
    """
    Get GBM relevance annotation.
    
    Args:
        drug_name: Drug name
        targets: Target genes
    
    Returns:
        Relevance description
    """
    # Known anti-GBM activity
    known_gbm_drugs = ['Metformin', 'Atovaquone', 'Dichloroacetate', 'Valproic Acid', 'Disulfiram']
    
    if drug_name in known_gbm_drugs:
        return 'Literature supported'
    
    # Metabolic targets
    metabolic_targets = ['NDUFV1', 'UQCRB', 'FASN', 'DHODH', 'HMGCR', 'PDK1', 'LDHA', 'HK2']
    
    for target in targets:
        if target in metabolic_targets:
            return 'Metabolic target'
    
    return 'Novel prediction'


def plot_drug_repositioning(results_df: pd.DataFrame) -> 'go.Figure':
    """
    Create drug repositioning visualization.
    
    Args:
        results_df: Drug repositioning results
    
    Returns:
        Plotly Figure object
    """
    from plotly.subplots import make_subplots
    
    if go is None or make_subplots is None:
        raise ImportError("plotly is required for visualization")
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Repositioning Scores", "Priority Distribution"),
        horizontal_spacing=0.15
    )
    
    # Sort by score
    df_sorted = results_df.sort_values('repositioning_score', ascending=True)
    
    # Color by priority
    priority_colors = {
        'HIGH': '#E94F37',
        'MEDIUM': '#F18F01',
        'LOW': '#2E86AB'
    }
    
    bar_colors = [priority_colors.get(p, '#CCCCCC') for p in df_sorted['priority']]
    
    # Panel A: Bar chart
    fig.add_trace(
        go.Bar(
            x=df_sorted['repositioning_score'],
            y=df_sorted['drug_name'],
            orientation='h',
            marker_color=bar_colors,
            name='Drugs',
            hovertemplate='<b>%{y}</b><br>Score: %{x:.1f}<br>Priority: %{customdata}<extra></extra>',
            customdata=df_sorted['priority']
        ),
        row=1, col=1
    )
    
    # Panel B: Pie chart
    priority_counts = results_df['priority'].value_counts()
    
    fig.add_trace(
        go.Pie(
            labels=priority_counts.index,
            values=priority_counts.values,
            marker_colors=[priority_colors.get(p, '#CCCCCC') for p in priority_counts.index],
            name='Priority'
        ),
        row=1, col=2
    )
    
    # Update layout
    fig.update_layout(
        height=600,
        width=1200,
        title_text="Drug Repositioning Screen Results",
        title_font_size=18,
        showlegend=False,
        template="plotly_white"
    )
    
    fig.update_xaxes(title_text="Repositioning Score (0-100)", row=1, col=1)
    fig.update_yaxes(title_text="Drug", row=1, col=1)
    
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='ARES-GBM Drug Repositioning Screen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python drug_repositioning.py --expr data/expression.csv
  python drug_repositioning.py --expr data/expression.csv --model Human-GEM.xml
  python drug_repositioning.py --expr data/expression.csv --output results/repositioning
        """
    )
    
    parser.add_argument('--expr', required=True, help='Path to expression data CSV')
    parser.add_argument('--model', help='Path to SBML model (optional)')
    parser.add_argument('--high', type=float, default=75, help='High expression percentile (default: 75)')
    parser.add_argument('--low', type=float, default=25, help='Low expression percentile (default: 25)')
    parser.add_argument('--output', default='results', help='Output directory (default: results)')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    run_drug_repositioning(
        expr_path=args.expr,
        model_path=args.model,
        high_percentile=args.high,
        low_percentile=args.low,
        output_dir=args.output,
        verbose=not args.quiet
    )


if __name__ == '__main__':
    main()
