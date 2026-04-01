# -*- coding: utf-8 -*-
"""
ARES-GBM Research Tools — Utilities

Helper functions for data loading, model preparation, and visualization.
"""
import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from typing import Optional, Dict, List, Tuple


def load_expression_matrix(filepath: str, index_col: int = 0) -> pd.DataFrame:
    """
    Load a gene expression matrix from CSV.
    
    Args:
        filepath: Path to CSV file
        index_col: Column to use as index (default: first column)
    
    Returns:
        DataFrame with genes as rows and samples as columns
    """
    df = pd.read_csv(filepath, index_col=index_col)
    df = df.apply(pd.to_numeric, errors='coerce')
    return df


def load_survival_data(filepath: str) -> pd.DataFrame:
    """
    Load clinical survival data with columns: sample, survival_days, event.
    
    Args:
        filepath: Path to CSV file
    
    Returns:
        DataFrame indexed by sample ID
    """
    df = pd.read_csv(filepath, index_col=0)
    required_cols = ['survival_days', 'event']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    return df


def prepare_model(model, solver: str = "glpk"):
    """
    Prepare a cobrapy model for analysis.
    
    Args:
        model: cobrapy.Model object
        solver: Solver to use (glpk, glpk_exact, scipy, cplex, gurobi)
    
    Returns:
        Prepared model with solver set
    """
    model.solver = solver
    return model


def apply_imat(
    model,
    expression: pd.Series,
    high_percentile: float = 75,
    low_percentile: float = 25,
    verbose: bool = False
):
    """
    Apply iMAT algorithm to integrate gene expression into metabolic model.
    
    Reactions associated with lowly expressed genes are constrained to low flux.
    
    Args:
        model: cobrapy.Model object
        expression: Gene expression values (gene_id -> value)
        high_percentile: Percentile threshold for high expression
        low_percentile: Percentile threshold for low expression
        verbose: Print progress information
    
    Returns:
        Model with constraints applied
    """
    m = model.copy()
    
    # Calculate thresholds
    expr_values = pd.to_numeric(expression, errors='coerce').dropna()
    high_thresh = np.percentile(expr_values.values, high_percentile)
    low_thresh = np.percentile(expr_values.values, low_percentile)
    
    if verbose:
        print(f"Expression thresholds: high={high_thresh:.2f}, low={low_thresh:.2f}")
    
    # Build expression dictionary
    expr_dict = expr_values.to_dict()
    
    # Apply constraints
    constrained = 0
    for rxn in m.reactions:
        if not rxn.genes:
            continue
        
        # Get expression values for genes in this reaction
        gene_expr = [expr_dict.get(g.id) for g in rxn.genes if g.id in expr_dict]
        
        if not gene_expr:
            continue
        
        # If max expression is below low threshold, constrain the reaction
        if max(gene_expr) <= low_thresh and rxn.lower_bound >= 0:
            rxn.upper_bound = 0.01
            constrained += 1
    
    if verbose:
        print(f"Constrained {constrained} reactions based on low expression")
    
    return m


def get_gene_symbol_map(model) -> Dict[str, str]:
    """
    Create a mapping from gene symbols to ENSEMBL IDs.
    
    Args:
        model: cobrapy.Model object
    
    Returns:
        Dictionary mapping symbol -> ensembl_id
    """
    symbol_map = {}
    for gene in model.genes:
        symbol = gene.annotation.get('hgnc.symbol')
        if symbol:
            symbol_map[symbol] = gene.id
    return symbol_map


def optimize_model(model) -> Tuple[float, pd.Series]:
    """
    Optimize a metabolic model and return objective value and fluxes.
    
    Args:
        model: cobrapy.Model object
    
    Returns:
        Tuple of (objective_value, fluxes_series)
    """
    solution = model.optimize()
    return solution.objective_value, solution.fluxes


def plot_flux_comparison(
    flux_gbm: pd.Series,
    flux_astro: pd.Series,
    top_n: int = 20,
    title: str = "Differential Flux Analysis"
) -> go.Figure:
    """
    Create a bar plot comparing fluxes between two conditions.
    
    Args:
        flux_gbm: Flux values for GBM condition
        flux_astro: Flux values for astrocyte condition
        top_n: Number of top differential reactions to show
        title: Plot title
    
    Returns:
        Plotly Figure object
    """
    df = pd.DataFrame({
        'flux_gbm': flux_gbm,
        'flux_astro': flux_astro
    }).fillna(0)
    
    df['delta'] = df['flux_gbm'] - df['flux_astro']
    df['abs_delta'] = df['delta'].abs()
    df['reaction'] = df.index
    
    # Get top N differential reactions
    top = df.nlargest(top_n, 'abs_delta')
    
    fig = px.bar(
        top,
        x='abs_delta',
        y='reaction',
        orientation='h',
        color='abs_delta',
        color_continuous_scale=['#BDC3C7', '#C0392B'],
        title=title,
        labels={'abs_delta': '|Δ Flux|', 'reaction': 'Reaction'}
    )
    
    fig.update_layout(
        showlegend=False,
        height=400,
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=10, r=10, t=40, b=10)
    )
    
    return fig


def plot_knockout_results(
    results: List[Dict],
    title: str = "Gene Knockout Results"
) -> go.Figure:
    """
    Create a bar plot of gene knockout effects.
    
    Args:
        results: List of dicts with keys: genes, ratio, effect
        title: Plot title
    
    Returns:
        Plotly Figure object
    """
    df = pd.DataFrame(results)
    df = df.sort_values('ratio')
    
    fig = px.bar(
        df,
        x='ratio',
        y='genes',
        orientation='h',
        color='ratio',
        color_continuous_scale=['#C0392B', '#27AE60'],
        range_color=[0, 1],
        title=title,
        labels={'ratio': 'Growth Ratio (KO/WT)', 'genes': 'Gene'}
    )
    
    # Add lethal threshold line
    fig.add_vline(x=0.1, line_dash='dash', line_color='#4A4A4A', 
                  annotation_text='Lethal threshold')
    
    fig.update_layout(
        showlegend=False,
        height=400,
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=10, r=10, t=40, b=10)
    )
    
    return fig


def plot_kaplan_meier(
    survival_high: pd.DataFrame,
    survival_low: pd.DataFrame,
    gene_name: str,
    title: Optional[str] = None
) -> go.Figure:
    """
    Create Kaplan-Meier survival curves.
    
    Args:
        survival_high: DataFrame with survival_days and event for high expression group
        survival_low: DataFrame with survival_days and event for low expression group
        gene_name: Name of the gene being analyzed
        title: Optional plot title
    
    Returns:
        Plotly Figure object
    """
    fig = go.Figure()
    
    for group, color, dash in [
        ('HIGH', '#C0392B', 'solid'),
        ('LOW', '#2980B9', 'dash')
    ]:
        data = survival_high if group == 'HIGH' else survival_low
        data = data.sort_values('survival_days')
        
        # Calculate survival curve
        n = len(data)
        times = [0] + data['survival_days'].tolist()
        survival = [1.0]
        s = 1.0
        
        for i, event in enumerate(data['event'].tolist()):
            if event:
                s *= (1 - 1 / (n - i))
            survival.append(s)
        
        fig.add_trace(go.Scatter(
            x=times,
            y=survival,
            mode='lines',
            name=f'{group} (n={len(data)})',
            line=dict(color=color, width=2.5, dash=dash)
        ))
    
    fig.update_layout(
        title=title or f'Kaplan-Meier — {gene_name}',
        xaxis_title='Days',
        yaxis_title='Survival Probability',
        yaxis_range=[-0.05, 1.05],
        height=400,
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=10, r=10, t=40, b=10)
    )
    
    return fig


def save_plot(fig: go.Figure, filepath: str, width: int = 800, height: int = 450):
    """
    Save a Plotly figure to PNG.
    
    Args:
        fig: Plotly Figure object
        filepath: Output file path
        width: Image width in pixels
        height: Image height in pixels
    """
    fig.write_image(filepath, width=width, height=height, scale=2)
    print(f"Plot saved to: {filepath}")


def ensure_dir(filepath: str) -> str:
    """
    Ensure the directory for a filepath exists.
    
    Args:
        filepath: File path
    
    Returns:
        The input filepath
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    return filepath
