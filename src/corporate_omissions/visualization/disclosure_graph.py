"""Scope-1 disclosure composition plot (by year, sector, firm-size decile)."""
from __future__ import annotations

from typing import Literal, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, MaxNLocator


_DISCLOSED_COLOR = "#3B75AF"
_ESTIMATED_COLOR = "#B6C8DE"
_MISSING_EDGE = "black"

_RC_PARAMS = {
    'font.family':       'serif',
    'font.serif':        ['Times New Roman', 'DejaVu Serif'],
    'font.size':         11,
    'axes.titlesize':    12,
    'axes.labelsize':    11,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'legend.fontsize':   10,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.grid':         True,
    'grid.alpha':        0.25,
    'grid.linewidth':    0.5,
    'figure.dpi':        150,
    'savefig.dpi':       300,
}


def _compute_shares(df: pd.DataFrame, by: str, id_col: str,
                    disclosed_col: str, estimated_col: str) -> pd.DataFrame:
    """Per-group shares (%) of Disclosed / Estimated / Missing firms."""
    g = df.groupby(by, sort=True, observed=True)
    n = g[id_col].nunique().rename("n")
    disclosed = (
        df.loc[df[disclosed_col].notna()]
        .groupby(by, observed=True)[id_col].nunique()
        .rename("disclosed")
    )
    estimated = (
        df.loc[df[estimated_col].notna()]
        .groupby(by, observed=True)[id_col].nunique()
        .rename("estimated")
    )
    out = pd.concat([n, disclosed, estimated], axis=1).fillna(0)
    out["Disclosed"] = 100.0 * out["disclosed"] / out["n"]
    out["Estimated"] = 100.0 * out["estimated"] / out["n"]
    out["Missing"] = (100.0 - out["Disclosed"] - out["Estimated"]).clip(lower=0)
    return out[["n", "Disclosed", "Estimated", "Missing"]]


def _stacked_bar(ax, xvals, shares: pd.DataFrame, *, width: float = 0.8) -> None:
    disclosed = shares["Disclosed"].to_numpy()
    estimated = shares["Estimated"].to_numpy()
    missing = shares["Missing"].to_numpy()

    ax.bar(xvals, disclosed, width=width, label="Disclosed", color=_DISCLOSED_COLOR)
    ax.bar(xvals, estimated, width=width, bottom=disclosed,
           label="Estimated", color=_ESTIMATED_COLOR)
    ax.bar(xvals, missing, width=width, bottom=disclosed + estimated,
           label="Missing", facecolor="white",
           edgecolor=_MISSING_EDGE, linestyle="--", linewidth=0.8)

    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _panel_label(ax, text: str) -> None:
    ax.text(-0.06, 1.04, text, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left")


def plot_scope1_disclosure_pct(
    panel: pd.DataFrame,
    *,
    year_col: str = "year",
    id_col: str = "gvkey",
    sector_col: str = "gsector",
    assets_col: str = "at",
    disclosed_col: str = "sc1_disclosed",
    estimated_col: str = "sc1_estimated",
    start_year: Optional[int] = None,
    sub_year: Optional[int] = None,
    n_size_bins: int = 10,
    sector_order: Literal["disclosure", "sample_size", "alpha"] = "disclosure",
    figsize: tuple[float, float] = (12, 12),
) -> plt.Figure:
    """Three-panel stacked-bar figure: disclosure composition by year, sector, size."""
    required = [year_col, id_col, sector_col, assets_col, disclosed_col, estimated_col]
    missing = [c for c in required if c not in panel.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = panel if start_year is None else panel.loc[panel[year_col] >= start_year]
    df_sub = df if sub_year is None else df.loc[df[year_col] == sub_year]

    # ---- Panel A: by year
    by_year = _compute_shares(df, year_col, id_col, disclosed_col, estimated_col)

    # ---- Panel B: by sector (configurable ordering)
    by_sector = _compute_shares(df_sub, sector_col, id_col, disclosed_col, estimated_col)
    if sector_order == "disclosure":
        by_sector = by_sector.sort_values("Disclosed", ascending=False)
    elif sector_order == "sample_size":
        by_sector = by_sector.sort_values("n", ascending=False)
    elif sector_order == "alpha":
        by_sector = by_sector.sort_index()

    # ---- Panel C: by firm-size decile
    x = df_sub[assets_col].astype(float)
    try:
        bins = pd.qcut(x, q=n_size_bins, labels=False, duplicates="drop")
    except ValueError:
        bins = pd.qcut(x.rank(method="first"), q=n_size_bins,
                       labels=False, duplicates="drop")
    df_sub = df_sub.assign(size_bin=bins.astype("Int64") + 1)
    by_size = _compute_shares(
        df_sub.dropna(subset=["size_bin"]),
        "size_bin", id_col, disclosed_col, estimated_col,
    )
    size_labels = [f"Q{int(b)}" for b in by_size.index]
    size_labels[0] = "Smallest"
    size_labels[-1] = "Largest"

    # ---- Figure
    with plt.rc_context(_RC_PARAMS):
        fig, axes = plt.subplots(3, 1, figsize=figsize)
        axA, axB, axC = axes

        _stacked_bar(axA, by_year.index.to_numpy(), by_year)
        axA.set_title("Scope 1 disclosure composition by year")
        axA.set_ylabel("Share of firms")
        axA.xaxis.set_major_locator(MaxNLocator(integer=True))
        axA.set_xticks(by_year.index.to_numpy())
        axA.set_xticklabels([str(int(y)) for y in by_year.index], rotation=0)

        sector_x = np.arange(len(by_sector))
        _stacked_bar(axB, sector_x, by_sector)
        axB.set_title("Scope 1 disclosure composition by sector"
                      + (f" ({sub_year})" if sub_year is not None else ""))
        axB.set_ylabel("Share of firms")
        axB.set_xticks(sector_x)
        axB.set_xticklabels(by_sector.index, rotation=25, ha="right")

        size_x = np.arange(len(by_size))
        _stacked_bar(axC, size_x, by_size)
        axC.set_title(f"Scope 1 disclosure composition by firm size"
                      f" (q={n_size_bins}"
                      + (f", {sub_year}" if sub_year is not None else "") + ")")
        axC.set_ylabel("Share of firms")
        axC.set_xlabel("Firm size (asset quantiles)")
        axC.set_xticks(size_x)
        axC.set_xticklabels(size_labels)

        for ax, tag in zip(axes, ("A", "B", "C")):
            _panel_label(ax, tag)

        handles, labels = axA.get_legend_handles_labels()
        fig.legend(handles, labels, title="Disclosure type",
                   loc="upper left", bbox_to_anchor=(0.85, 0.98),
                   frameon=False)

        fig.tight_layout(rect=(0, 0, 0.90, 0.97))
    return fig
