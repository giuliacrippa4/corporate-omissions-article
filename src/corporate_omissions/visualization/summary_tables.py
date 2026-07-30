
# src/paper/tables.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
import pandas as pd

# Canonical treatment definition — single source of truth shared with the
# disclosed/estimated split and co2_bool (see data/preprocessing.py).
from corporate_omissions.data.preprocessing import classify_emissions_status, emissions_tier

EmissionsStatus = Literal["Observed", "Missing"]


@dataclass(frozen=True)
class SummaryUnits:
    """
    Unit conversions for presentation.
    """
    market_cap_div: float = 1e6      # report market cap in USD millions
    emissions_div: float = 1e9       # report emissions in (metric tons) billions => "Gt" scale; adjust label in paper
    round_ndigits: int = 2


def add_emissions_status(
    panel: pd.DataFrame,
    *,
    disclosure_col: str = "sc1_disclosure",
    out_col: str = "emissions_status",
    definition: str = "narrow",
) -> pd.DataFrame:
    """
    Return a copy of panel with a new column `emissions_status` in {'Observed','Missing'},
    classified under the given disclosed-`definition` ('narrow' | 'broad').
    """
    if disclosure_col not in panel.columns:
        raise ValueError(
            f"[add_emissions_status] Missing '{disclosure_col}'. "
            f"Available columns: {list(panel.columns)}"
        )

    out = panel.copy()
    out[out_col] = out[disclosure_col].apply(
        lambda v: classify_emissions_status(v, definition)
    )
    return out


def summarize_by_year(
    panel: pd.DataFrame,
    *,
    year_col: str = "year",
    mcap_col: str = "mkvalt",
    emissions_col: str = "sc1",
    intensity_col: str = "int_sc1",
    units: SummaryUnits = SummaryUnits(),
) -> pd.DataFrame:
    """
    Year-level summary (all firms, no Observed/Missing split).

    Returns columns:
      - n_obs (count of rows)
      - market_cap (sum, scaled)
      - scope1_emissions (sum, scaled)
      - co2_intensity (mean)
    """
    required = [year_col, mcap_col, emissions_col, intensity_col]
    missing = [c for c in required if c not in panel.columns]
    if missing:
        raise ValueError(f"[summarize_by_year] Missing required columns: {missing}")

    df = (
        panel.groupby(year_col)
        .agg(
            n_obs=(year_col, "size"),
            market_cap=(mcap_col, "sum"),
            scope1_emissions=(emissions_col, "sum"),
            co2_intensity=(intensity_col, "mean"),
        )
        .reset_index()
    )

    df["market_cap"] = (df["market_cap"] / units.market_cap_div).round(units.round_ndigits)
    df["scope1_emissions"] = (df["scope1_emissions"] / units.emissions_div).round(units.round_ndigits)
    df["co2_intensity"] = df["co2_intensity"].round(units.round_ndigits)
    return df


def summarize_by_year_and_status(
    panel: pd.DataFrame,
    *,
    year_col: str = "year",
    status_col: str = "emissions_status",
    mcap_col: str = "mkvalt",
    emissions_col: str = "sc1",
    intensity_col: str = "int_sc1",
    units: SummaryUnits = SummaryUnits(),
    wide: bool = True,
) -> pd.DataFrame:
    """
    Summary split by year × emissions_status.

    If wide=True (default), returns a wide table with a MultiIndex columns structure:
      metric × {Missing, Observed}

    If wide=False, returns the long (tidy) version:
      year, emissions_status, n_obs, market_cap, scope1_emissions, co2_intensity
    """
    required = [year_col, status_col, mcap_col, emissions_col, intensity_col]
    missing = [c for c in required if c not in panel.columns]
    if missing:
        raise ValueError(f"[summarize_by_year_and_status] Missing required columns: {missing}")

    df = (
        panel.groupby([year_col, status_col])
        .agg(
            n_obs=(status_col, "size"),
            market_cap=(mcap_col, "sum"),
            scope1_emissions=(emissions_col, "sum"),
            co2_intensity=(intensity_col, "mean"),
        )
        .reset_index()
    )

    df["market_cap"] = (df["market_cap"] / units.market_cap_div).round(units.round_ndigits)
    df["scope1_emissions"] = (df["scope1_emissions"] / units.emissions_div).round(units.round_ndigits)
    df["co2_intensity"] = df["co2_intensity"].round(units.round_ndigits)

    if not wide:
        return df

    wide_df = (
        df.set_index([year_col, status_col])
          .sort_index()
          .unstack(status_col)
    )
    return wide_df

# ---------------------------------------------------------------------
# Canonical paper tables
# ---------------------------------------------------------------------

def make_table1_summary(
    panel: pd.DataFrame,
    *,
    units: SummaryUnits = SummaryUnits(),
    definition: str = "narrow",
) -> pd.DataFrame:
    """
    Table 1 — Summary statistics by year (full sample).

    This is the canonical definition of Table 1 in the paper.
    """
    panel = add_emissions_status(panel, definition=definition)
    table1 = summarize_by_year(panel, units=units)

    # enforce stable column order
    cols = ["year", "n_obs", "market_cap", "scope1_emissions", "co2_intensity"]
    table1 = table1[cols]

    return table1

def make_table2_disclosure_types(
    panel: pd.DataFrame,
    *,
    year: int,
    disclosure_col: str = "sc1_disclosure",
    definition: str = "narrow",
) -> pd.DataFrame:
    """
    Table 2 — Breakdown of disclosure types in Trucost data (single year).

    Returns a DataFrame with columns:
      - disclosure_source
      - count
      - share_pct
      - tier       (exact / derived_deterministic / fuel / estimated / none)
      - treatment  (Observed / Missing, under the given `definition`)

    `tier` is definition-agnostic; `treatment` depends on `definition`
    ('narrow' counts only 'exact' as Observed; 'broad' also counts
    'derived_deterministic').
    """
    if disclosure_col not in panel.columns:
        raise ValueError(
            f"[make_table2_disclosure_types] Missing '{disclosure_col}' column. "
            f"Available columns: {list(panel.columns)}"
        )

    df = panel.loc[panel["year"] == year, disclosure_col].copy()

    counts = (
        df.value_counts(dropna=False)
          .rename_axis("disclosure_source")
          .reset_index(name="count")
    )

    total = counts["count"].sum()

    # Shares in percent
    counts["share_pct"] = 100 * counts["count"] / total

    # Provenance tier (definition-agnostic) and treatment (under `definition`).
    counts["tier"] = counts["disclosure_source"].apply(emissions_tier)
    counts["treatment"] = counts["disclosure_source"].apply(
        lambda v: classify_emissions_status(v, definition)
    )

    # Replace NaN label explicitly
    counts["disclosure_source"] = (
        counts["disclosure_source"]
        .astype(str)
        .replace("nan", "NaN")
    )

    # Round shares nicely (LaTeX formatting can handle <0.05 later)
    counts["share_pct"] = counts["share_pct"].round(2)

    # Stable column order
    counts = counts[
        ["disclosure_source", "count", "share_pct", "tier", "treatment"]
    ]

    return counts



def make_table3_summary(
    panel: pd.DataFrame,
    *,
    units: SummaryUnits = SummaryUnits(),
    wide: bool = True,
    definition: str = "narrow",
) -> pd.DataFrame:
    """
    Table 3 — Summary statistics by year and disclosure status.

    If wide=True (default): Observed vs Missing side-by-side (paper version).
    If wide=False: long/tidy version (appendix or robustness).
    """
    panel = add_emissions_status(panel, definition=definition)
    table2 = summarize_by_year_and_status(panel, units=units, wide=wide)

    if wide:
        # enforce metric order at top level
        metric_order = ["n_obs", "market_cap", "scope1_emissions", "co2_intensity"]
        table2 = table2.reindex(metric_order, axis=1, level=0)

    return table2

def to_latex_table(
    df: pd.DataFrame,
    *,
    caption: str,
    label: str,
    float_format: str = "%.2f",
) -> str:
    """
    Standard LaTeX export for paper tables (booktabs).
    """
    return df.to_latex(
        index=False if not isinstance(df.index, pd.MultiIndex) else True,
        float_format=float_format,
        escape=False,
        multicolumn=True,
        multicolumn_format="c",
        caption=caption,
        label=label,
    )