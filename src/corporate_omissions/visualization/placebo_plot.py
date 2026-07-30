"""Placebo-distribution figure for the SEC event study.

The randomization inference itself is produced by ``Returns/03a_event_study_placebos.R``,
which writes one row per candidate date to ``returns_event_placebo_raw.csv``. This
module redraws the figure in matplotlib so it matches the styling of the other
paper figures (gamma-by-year, gamma-by-sector, gamma-by-IO), which are produced
from the notebooks rather than from ggplot.

Usage
-----
    from corporate_omissions.visualization.placebo_plot import plot_placebo_distribution

    plot_placebo_distribution(
        csv_path="../data/outputs/tables/returns_event_placebo_raw.csv",
        savepath="../data/outputs/figures/returns_event_placebo_dist.pdf",
    )
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Matches the palette used by the gamma figures in 04b: a muted blue for the
# reference mass, warm ochre for the flagged dates, red for the event itself.
C_PLACEBO = "#a8c4e0"
C_CONFOUND = "#f4a020"
C_EVENT = "#c0392b"
C_ZERO = "#555555"

WINDOW_LABELS = {
    "short": r"$[-1,+1]$",
    "medium": r"$[-5,+5]$",
    "long": r"$[-30,+30]$",
}


def _apply_style() -> None:
    """Match the rcParams block used in the imputation notebooks."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Dejavu Serif"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })


def plot_placebo_distribution(
    csv_path: str | Path,
    savepath: Optional[str | Path] = None,
    window: str = "medium",
    bins: Optional[int] = None,
    figsize: tuple[float, float] = (8, 4),
    show: bool = True,
):
    """Histogram of placebo coefficients with the event coefficient marked.

    Parameters
    ----------
    csv_path
        ``returns_event_placebo_raw.csv`` as written by 03a.
    window
        Which event window to plot: ``"short"``, ``"medium"`` or ``"long"``.
    bins
        Histogram bins; defaults to a Freedman-Diaconis-ish rule capped for
        readability, since the annual scheme yields only ~10 placebos while
        the monthly scheme yields ~131.
    """
    csv_path = Path(csv_path)
    d = pd.read_csv(csv_path)

    coef_col = f"coef_{window}"
    if coef_col not in d.columns:
        raise KeyError(f"{coef_col!r} not in {csv_path.name}; columns are {list(d.columns)}")

    # data.table writes logicals as TRUE/FALSE strings
    for c in ("is_event", "confounded"):
        if d[c].dtype == object:
            d[c] = d[c].astype(str).str.upper().eq("TRUE")

    d = d.dropna(subset=[coef_col])
    ev_row = d[d["is_event"]]
    if len(ev_row) != 1:
        raise ValueError(f"expected exactly one event row, found {len(ev_row)}")
    ev = float(ev_row[coef_col].iloc[0])

    pl = d[~d["is_event"]]
    pl_clean = pl[~pl["confounded"]][coef_col].to_numpy()
    pl_conf = pl[pl["confounded"]][coef_col].to_numpy()
    n_below = int((pl[coef_col] <= ev).sum())
    k = len(pl)
    p_emp = (n_below + 1) / (k + 1)

    _apply_style()
    fig, ax = plt.subplots(figsize=figsize)

    if bins is None:
        bins = max(8, min(30, int(np.sqrt(k) * 2)))
    edges = np.histogram_bin_edges(pl[coef_col].to_numpy(), bins=bins)

    # Stacked so the flagged dates are visible without hiding the total mass.
    ax.hist(
        [pl_clean, pl_conf],
        bins=edges,
        stacked=True,
        color=[C_PLACEBO, C_CONFOUND],
        edgecolor="white",
        linewidth=0.4,
        label=["placebo", "placebo (confounded)"],
        zorder=2,
    )

    ax.axvline(0, color=C_ZERO, lw=0.8, ls="--", zorder=1)
    ax.axvline(ev, color=C_EVENT, lw=1.6, zorder=4)

    ymax = ax.get_ylim()[1]
    ax.annotate(
        "SEC proposal",
        xy=(ev, ymax * 0.97),
        xytext=(6, 0),
        textcoords="offset points",
        ha="left",
        va="top",
        color=C_EVENT,
        fontsize=10,
    )

    ax.set_xlabel(
        f"Coefficient on log recovered emissions, CAR {WINDOW_LABELS.get(window, window)}"
    )
    ax.set_ylabel("Placebo dates")
    ax.set_title(
        "Event coefficient against the placebo distribution\n"
        rf"$K={k}$ placebo dates; {n_below} below the event, $p_{{emp}}={p_emp:.3f}$",
        fontsize=11,
    )
    ax.legend(frameon=False, loc="upper left", fontsize=9)

    fig.tight_layout()
    if savepath is not None:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savepath, bbox_inches="tight", dpi=300)
        print(f"Saved: {savepath.name}  (K={k}, n_below={n_below}, p_emp={p_emp:.3f})")
    if show:
        plt.show()
    return fig, ax


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[3]
    plot_placebo_distribution(
        csv_path=root / "data/outputs/tables/returns_event_placebo_raw.csv",
        savepath=root / "data/outputs/figures/returns_event_placebo_dist.pdf",
        show=False,
    )
