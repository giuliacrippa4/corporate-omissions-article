from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class ResponseBoxplotSpec:
    """
    Plot spec for LM response mechanism coefficient boxplots across MC draws.
    """
    # Core coefficient names (statsmodels naming)
    coef_const: str = "const"
    coef_size: str = "at"
    coef_emis: str = "sc1_disclosed_complete"

    # Label mapping (you can override)
    rename_dict: Optional[Dict[str, str]] = None

    # Include sector fixed effects (gsector_*)
    include_sectors: bool = True
    sector_prefix: str = "gsector_"

    # Whether to sort sector dummies by mean effect (within MC)
    sort_sectors_by_mean: bool = True

    # Scaling: multiply slopes by SD of their variable (to interpret as 1 SD change)
    scale_by_sd: bool = True

    # Which df column to use for SD scaling
    size_col_in_df: str = "at"
    emis_col_in_df: str = "sc1_disclosed_complete"

    # Plot aesthetics
    horizontal: bool = True
    figsize: Optional[Tuple[float, float]] = None
    box_width: float = 0.6
    mean_tick_halfwidth: float = 0.30  # in axis units for the categorical axis

    title: str = "Distribution of LM response-mechanism coefficients across MC runs"
    xlabel: str = "Log-odds coefficient"
    ylabel: str = ""  # y-label is usually empty for horizontal coefficient lists

    # Style controls
    show_zero_line: bool = True
    grid: bool = True


def _stack_params(mc_results: Union[pd.DataFrame, Sequence[dict]], params_col: str) -> pd.DataFrame:
    """
    Convert a column of pd.Series coefficients into a (n_mc x n_params) DataFrame.
    """
    if isinstance(mc_results, pd.DataFrame):
        series_list = mc_results[params_col].tolist()
    else:
        series_list = [row[params_col] for row in mc_results]

    # Each element should be a pd.Series-like with index = coefficient names
    return pd.DataFrame(series_list)


def _get_reference_df(
    mc_results: pd.DataFrame,
    df_ref: Optional[pd.DataFrame],
    df_ref_col: str,
    df_ref_row: int = 0,
) -> pd.DataFrame:
    """
    Choose a stable reference df to compute SD scaling.
    """
    if df_ref is not None:
        return df_ref
    if df_ref_col not in mc_results.columns:
        raise KeyError(
            f"df_ref_col='{df_ref_col}' not found in mc_results columns. "
            f"Provide df_ref explicitly or store a reference df in mc_results['{df_ref_col}']."
        )
    return mc_results.loc[df_ref_row, df_ref_col]


def plot_response_lm_coef_boxplot(
    mc_results: pd.DataFrame,
    *,
    params_col: str = "lm_params",
    df_ref: Optional[pd.DataFrame] = None,
    df_ref_col: str = "df_gam",
    df_ref_row: int = 0,
    spec: Optional[ResponseBoxplotSpec] = None,
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    savepath: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
) -> plt.Axes:
    """
    Plot boxplots of LM response-mechanism coefficients across MC runs.

    Parameters
    ----------
    mc_results : pd.DataFrame
        Must contain a column `params_col` whose entries are pd.Series of coefficients.
        Optionally contains a column `df_ref_col` with a df used to compute SD scaling.
    params_col : str
        Column name for coefficient Series. Default 'lm_params'.
    df_ref : Optional[pd.DataFrame]
        Explicit reference df for SD scaling (recommended for maximum control).
        If None, uses mc_results[df_ref_col].iloc[df_ref_row].
    df_ref_col : str
        Column in mc_results holding a reference df (only used if df_ref is None).
    df_ref_row : int
        Row index to pull reference df from (only used if df_ref is None).
    spec : Optional[ResponseBoxplotSpec]
        Plot spec controlling selection, scaling, labels, orientation.
    include/exclude : Optional[Sequence[str]]
        Optional explicit coefficient whitelist/blacklist after automatic selection.
    savepath : Optional[str]
        If provided, saves the figure to this path (png/pdf/etc).
    dpi : int
        Save DPI.
    show : bool
        If True, calls plt.show().

    Returns
    -------
    matplotlib Axes
    """
    if spec is None:
        spec = ResponseBoxplotSpec()

    # Default label mapping
    if spec.rename_dict is None:
        spec.rename_dict = {
            spec.coef_const: "Constant",
            spec.coef_size: "Size (1 SD)",
            spec.coef_emis: "Emissions (1 SD)",
            "green_sector": "Green Sector",
        }

    # 1) Stack coefficients
    B = _stack_params(mc_results, params_col=params_col)

    # 2) Choose coefficients to plot
    base_cols: List[str] = [c for c in [spec.coef_const, spec.coef_size, spec.coef_emis] if c in B.columns]

    sector_cols: List[str] = []
    if spec.include_sectors:
        sector_cols = [c for c in B.columns if c.startswith(spec.sector_prefix)]
        if spec.sort_sectors_by_mean and sector_cols:
            sector_cols = sorted(sector_cols, key=lambda c: float(B[c].mean()))

    plot_cols = base_cols + sector_cols

    # Apply include/exclude filters if provided
    if include is not None:
        include_set = set(include)
        plot_cols = [c for c in plot_cols if c in include_set]
    if exclude is not None:
        exclude_set = set(exclude)
        plot_cols = [c for c in plot_cols if c not in exclude_set]

    if not plot_cols:
        raise ValueError("No coefficients selected for plotting. Check spec/include/exclude.")

    # 3) Scale by SD if requested (only for the slopes, not for dummies/const)
    B_plot = B[plot_cols].copy()

    if spec.scale_by_sd:
        ref = _get_reference_df(mc_results, df_ref, df_ref_col=df_ref_col, df_ref_row=df_ref_row)

        if spec.coef_size in B_plot.columns and spec.size_col_in_df in ref.columns:
            sd_size = float(pd.to_numeric(ref[spec.size_col_in_df], errors="coerce").std())
            if np.isfinite(sd_size) and sd_size > 0:
                B_plot[spec.coef_size] = B_plot[spec.coef_size] * sd_size

        if spec.coef_emis in B_plot.columns and spec.emis_col_in_df in ref.columns:
            sd_emis = float(pd.to_numeric(ref[spec.emis_col_in_df], errors="coerce").std())
            if np.isfinite(sd_emis) and sd_emis > 0:
                B_plot[spec.coef_emis] = B_plot[spec.coef_emis] * sd_emis

    # 4) Build arrays + labels + means
    data = []
    labels = []
    means = []

    for col in plot_cols:
        vals = B_plot[col].dropna().to_numpy()
        data.append(vals)
        means.append(float(np.mean(vals)) if len(vals) else np.nan)

        lab = spec.rename_dict.get(col, col)
        if col.startswith(spec.sector_prefix):
            lab = lab.replace(spec.sector_prefix, "")
        labels.append(lab)

    # 5) Figure size defaults
    if spec.figsize is None:
        # Dynamic height based on number of coefficients for horizontal plot
        if spec.horizontal:
            spec.figsize = (9, max(3.5, 0.38 * len(plot_cols)))
        else:
            spec.figsize = (10, 5)

    fig, ax = plt.subplots(figsize=spec.figsize)

    # 6) Plot (horizontal recommended)
    if spec.horizontal:
        ax.boxplot(
            data,
            vert=False,
            patch_artist=True,
            widths=spec.box_width,
            boxprops=dict(facecolor="white", edgecolor="black", linewidth=1),
            medianprops=dict(color="black", linewidth=1.2),
            whiskerprops=dict(color="black", linewidth=1),
            capprops=dict(color="black", linewidth=1),
            flierprops=dict(marker="o", color="black", markersize=3, alpha=0.6),
        )

        ypos = np.arange(1, len(plot_cols) + 1)

        # red tick = MC mean (matches your caption)
        for y, m in zip(ypos, means):
            if np.isfinite(m):
                ax.vlines(m, y - spec.mean_tick_halfwidth, y + spec.mean_tick_halfwidth,
                          colors="tab:red", linewidth=2)

        if spec.show_zero_line:
            ax.axvline(0, color="gray", linewidth=1, linestyle="--")

        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=10)

        ax.set_xlabel(spec.xlabel, fontsize=10)
        if spec.ylabel:
            ax.set_ylabel(spec.ylabel, fontsize=10)

        if spec.grid:
            ax.grid(axis="x", linestyle="--", alpha=0.6)

    else:
        ax.boxplot(
            data,
            vert=True,
            patch_artist=True,
            widths=0.4,
            boxprops=dict(facecolor="white", edgecolor="black", linewidth=1),
            medianprops=dict(color="black", linewidth=1.2),
            whiskerprops=dict(color="black", linewidth=1),
            capprops=dict(color="black", linewidth=1),
            flierprops=dict(marker="o", color="black", markersize=4, alpha=0.6),
        )

        xpos = np.arange(1, len(plot_cols) + 1)

        # red tick = MC mean (matches your caption)
        for x, m in zip(xpos, means):
            if np.isfinite(m):
                ax.hlines(m, x - 0.18, x + 0.18, colors="tab:red", linewidth=2)

        if spec.show_zero_line:
            ax.axhline(0, color="gray", linewidth=1, linestyle="--")

        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, fontsize=10, rotation=75, ha="right")

        ax.set_ylabel(spec.xlabel, fontsize=10)  # for vertical, x-label is less meaningful
        if spec.grid:
            ax.grid(axis="y", linestyle="--", alpha=0.7)

    ax.set_title(spec.title, fontsize=12)
    plt.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()

    return ax
