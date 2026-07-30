# corporate_omissions/plots/plot_gam.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class GAMMCSpec:
    """
    Spec for aggregating MC partial dependence outputs from PropensityGAM.

    Assumes mc_results has columns:
      - pdep:   dict[int, np.ndarray]
      - confi:  dict[int, np.ndarray] with shape (n_grid, 2)
      - XX_dict: dict[int, np.ndarray] where each is (n_grid, n_features)

    Term indexing convention matches your PropensityGAM plotting code:
      term 0 -> s(size)
      term 1 -> s(y)
      term 2 -> s(gdp) if present
      then factor terms (sector, year)
    """
    pdep_col: str = "pdep"
    confi_col: str = "confi"
    XX_col: str = "XX_dict"

    size_term: int = 0
    emis_term: int = 1

    # Provide these explicitly to avoid guessing:
    sector_term_indices: Optional[List[int]] = None
    sector_labels: Optional[List[str]] = None

    # Standardize curves (your old behavior)
    standardize_continuous: bool = True

    # MC aggregation: keep "mean" to match your earlier paper phrasing
    # (You can extend later to quantiles if desired)
    agg: str = "mean"  # only "mean" implemented here


def _to_ndarray(x) -> np.ndarray:
    if isinstance(x, np.ndarray):
        return x
    return np.asarray(x)


def _zscore(arr: np.ndarray) -> np.ndarray:
    arr = _to_ndarray(arr).astype(float)
    sd = float(np.std(arr))
    if not np.isfinite(sd) or sd == 0.0:
        return arr
    return (arr - float(np.mean(arr))) / sd


def _stack_term_curves(
    mc_results: pd.DataFrame,
    term: int,
    *,
    pdep_col: str,
    confi_col: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Stack pdep and confi for a single term across MC runs.

    Returns:
      values_mat: (n_mc, n_grid)
      lo_mat:     (n_mc, n_grid)
      hi_mat:     (n_mc, n_grid)
    """
    vals = []
    lo = []
    hi = []

    for d_pdep, d_confi in zip(mc_results[pdep_col], mc_results[confi_col]):
        if d_pdep is None or term not in d_pdep:
            continue
        if d_confi is None or term not in d_confi:
            continue

        v = _to_ndarray(d_pdep[term]).astype(float)
        c = _to_ndarray(d_confi[term]).astype(float)  # (n_grid, 2)

        if c.ndim != 2 or c.shape[1] != 2:
            raise ValueError(f"confi[{term}] must be (n_grid, 2). Got {c.shape}.")

        vals.append(v)
        lo.append(c[:, 0])
        hi.append(c[:, 1])

    if len(vals) == 0:
        raise ValueError(f"No runs contained term={term} in both pdep and confi.")

    values_mat = np.vstack(vals)
    lo_mat = np.vstack(lo)
    hi_mat = np.vstack(hi)

    return values_mat, lo_mat, hi_mat


def _aggregate_mean(values_mat: np.ndarray, lo_mat: np.ndarray, hi_mat: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        np.mean(values_mat, axis=0),
        np.mean(lo_mat, axis=0),
        np.mean(hi_mat, axis=0),
    )


def _get_xgrid_1d(
    mc_results: pd.DataFrame,
    term: int,
    *,
    XX_col: str,
    feature_index_hint: Optional[int] = None,
) -> np.ndarray:
    """
    Extract a 1D x-grid for plotting a smooth term.

    We take the column in XX with the largest range unless feature_index_hint is provided.
    """
    XX_dict = mc_results[XX_col].iloc[0]
    if XX_dict is None or term not in XX_dict:
        raise KeyError(f"XX_dict missing term={term}. Store XX_dict in mc_results or pass correct XX_col.")

    XX = _to_ndarray(XX_dict[term]).astype(float)  # (n_grid, n_features)

    if XX.ndim == 1:
        return XX

    if feature_index_hint is not None:
        x = XX[:, feature_index_hint]
        return _to_ndarray(x)

    # heuristic: choose the most varying column
    ranges = np.ptp(XX, axis=0)
    j = int(np.argmax(ranges))
    return XX[:, j]


def _extract_dummy_effect_from_term(
    XX: np.ndarray,
    pdep: np.ndarray,
    confi: np.ndarray,
) -> Tuple[float, float, float]:
    """
    For a factor/dummy term, infer which rows correspond to 'dummy=1' and return:
      effect_mean, lo, hi

    Strategy:
      - Find the XX column with the largest range (should be the toggled dummy column).
      - Take rows where that column is (approximately) 1 as the 'on' state.
      - Aggregate pdep/confi over those rows (mean, and mean lo/hi).
    """
    XX = _to_ndarray(XX).astype(float)
    pdep = _to_ndarray(pdep).astype(float)
    confi = _to_ndarray(confi).astype(float)

    if confi.ndim != 2 or confi.shape[1] != 2:
        raise ValueError(f"confi for dummy term must be (n_grid, 2). Got {confi.shape}.")

    # choose the toggled column
    if XX.ndim == 1:
        xcol = XX
    else:
        ranges = np.ptp(XX, axis=0)
        j = int(np.argmax(ranges))
        xcol = XX[:, j]

    # select "on" rows (robust to float)
    mask_on = np.isclose(xcol, 1.0)
    if mask_on.sum() == 0:
        # fallback: take max xcol group as "on"
        mask_on = np.isclose(xcol, np.max(xcol))

    eff = float(np.mean(pdep[mask_on]))
    lo = float(np.mean(confi[mask_on, 0]))
    hi = float(np.mean(confi[mask_on, 1]))
    return eff, lo, hi


def process_gam_pdep_mc(
    mc_results: pd.DataFrame,
    *,
    spec: Optional[GAMMCSpec] = None,
) -> Dict[str, object]:
    """
    Process MC partial dependence outputs for:
      - size smooth (term spec.size_term)
      - emissions smooth (term spec.emis_term)
      - sector dummies (spec.sector_term_indices)

    Returns a dict ready for plotting, with keys:
      XX_size, mean_values_s0, mean_lower_s0, mean_upper_s0,
      XX_emis, mean_values_s1, mean_lower_s1, mean_upper_s1,
      sector_labels, mean_values_s2, mean_lower_s2, mean_upper_s2
    """
    if spec is None:
        spec = GAMMCSpec()

    # -------- size smooth --------
    vals_s0, lo_s0, hi_s0 = _stack_term_curves(
        mc_results, spec.size_term, pdep_col=spec.pdep_col, confi_col=spec.confi_col
    )
    mean_s0, mean_lo_s0, mean_hi_s0 = _aggregate_mean(vals_s0, lo_s0, hi_s0)
    XX_size = _get_xgrid_1d(mc_results, spec.size_term, XX_col=spec.XX_col)

    # -------- emissions smooth --------
    vals_s1, lo_s1, hi_s1 = _stack_term_curves(
        mc_results, spec.emis_term, pdep_col=spec.pdep_col, confi_col=spec.confi_col
    )
    mean_s1, mean_lo_s1, mean_hi_s1 = _aggregate_mean(vals_s1, lo_s1, hi_s1)
    XX_emis = _get_xgrid_1d(mc_results, spec.emis_term, XX_col=spec.XX_col)

    # Standardize (your previous behavior)
    if spec.standardize_continuous:
        mean_s0 = _zscore(mean_s0)
        mean_lo_s0 = _zscore(mean_lo_s0)
        mean_hi_s0 = _zscore(mean_hi_s0)

        mean_s1 = _zscore(mean_s1)
        mean_lo_s1 = _zscore(mean_lo_s1)
        mean_hi_s1 = _zscore(mean_hi_s1)

    # -------- sector effects --------
    if spec.sector_term_indices is None:
        raise ValueError(
            "spec.sector_term_indices must be provided. "
            "In your GAM these are the f() terms corresponding to gsector dummies."
        )
    if spec.sector_labels is None:
        raise ValueError("spec.sector_labels must be provided (one label per sector term).")

    if len(spec.sector_labels) != len(spec.sector_term_indices):
        raise ValueError("sector_labels and sector_term_indices must have same length.")

    # For each sector dummy term, compute 'on' effect and bounds, then average across runs
    sector_eff_runs = []
    sector_lo_runs = []
    sector_hi_runs = []

    for term in spec.sector_term_indices:
        eff_list = []
        lo_list = []
        hi_list = []

        for d_pdep, d_confi, d_XX in zip(mc_results[spec.pdep_col], mc_results[spec.confi_col], mc_results[spec.XX_col]):
            if d_pdep is None or term not in d_pdep:
                continue
            if d_confi is None or term not in d_confi:
                continue
            if d_XX is None or term not in d_XX:
                continue

            XX = d_XX[term]
            pdep = d_pdep[term]
            confi = d_confi[term]

            eff, lo, hi = _extract_dummy_effect_from_term(XX, pdep, confi)
            eff_list.append(eff)
            lo_list.append(lo)
            hi_list.append(hi)

        if len(eff_list) == 0:
            raise ValueError(f"No runs contained sector term={term} in pdep/confi/XX_dict.")

        sector_eff_runs.append(np.mean(eff_list))
        sector_lo_runs.append(np.mean(lo_list))
        sector_hi_runs.append(np.mean(hi_list))

    out = dict(
        XX_size=_to_ndarray(XX_size),
        mean_values_s0=_to_ndarray(mean_s0),
        mean_lower_s0=_to_ndarray(mean_lo_s0),
        mean_upper_s0=_to_ndarray(mean_hi_s0),
        XX_emis=_to_ndarray(XX_emis),
        mean_values_s1=_to_ndarray(mean_s1),
        mean_lower_s1=_to_ndarray(mean_lo_s1),
        mean_upper_s1=_to_ndarray(mean_hi_s1),
        sector_labels=list(spec.sector_labels),
        mean_values_s2=np.asarray(sector_eff_runs, dtype=float),
        mean_lower_s2=np.asarray(sector_lo_runs, dtype=float),
        mean_upper_s2=np.asarray(sector_hi_runs, dtype=float),
        spec=spec,
    )
    return out


def plot_gam_response_mc_3panel(
    processed: Dict[str, object],
    *,
    figsize: Tuple[float, float] = (14, 4),
    savepath: Optional[str] = None,
    show: bool = True,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Produce the 3-panel plot (Size, Emissions, Sector) from processed MC summaries.
    """
    XX_size = processed["XX_size"]
    XX_emis = processed["XX_emis"]

    fig, axs = plt.subplots(1, 3, figsize=figsize)

    # Size
    axs[0].plot(XX_size, processed["mean_values_s0"], color="tab:red")
    axs[0].fill_between(XX_size, processed["mean_lower_s0"], processed["mean_upper_s0"],
                        color="lightcoral", alpha=0.3)
    axs[0].set_title("Size", fontsize=10)
    axs[0].set_ylabel("Standardized Value")
    axs[0].set_xlabel("Million Total Assets")

    # Emissions
    axs[1].plot(XX_emis, processed["mean_values_s1"], color="tab:red")
    axs[1].fill_between(XX_emis, processed["mean_lower_s1"], processed["mean_upper_s1"],
                        color="lightcoral", alpha=0.3)
    axs[1].set_title("Emissions", fontsize=10)
    axs[1].set_xlabel("CO$_2$ Billion Metric Tonnes")

    # Sector
    labels = processed["sector_labels"]
    x = np.arange(len(labels))
    y = processed["mean_values_s2"]
    lo = processed["mean_lower_s2"]
    hi = processed["mean_upper_s2"]

    axs[2].scatter(x, y, color="tab:red", s=20, zorder=3)
    for i in range(len(x)):
        axs[2].vlines(i, lo[i], hi[i], color="black", linewidth=1, zorder=2)

    axs[2].set_xticks(x)
    axs[2].set_xticklabels(labels, rotation=45, ha="right")
    axs[2].set_title("Sector", fontsize=10)

    for ax in axs:
        ax.grid(linestyle="--", alpha=0.5)

    plt.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, bbox_inches="tight", dpi=300)
    if show:
        plt.show()

    return fig, axs
