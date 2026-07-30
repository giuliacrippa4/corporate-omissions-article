
# src/data/panel.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PanelPreprocessConfig:
    """
    Canonical post-merge panel preprocessing config.

    Keep this small and stable: these choices define the estimation sample.
    """
    subsample_sec: Optional[str] = None
    require_positive_assets: bool = True
    drop_negative_disclosed: bool = True
    disclosed_definition: str = "narrow"  # "narrow" (exact only) | "broad" (+ deterministic transforms)


def _require_cols(df: pd.DataFrame, cols: list[str], *, name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"[{name}] Missing required columns: {missing}\n"
            f"[{name}] Available columns: {list(df.columns)}"
        )


def emissions_tier(sc1_disclosure: object) -> str:
    """
    Assign a Trucost `sc1_disclosure` string to a provenance tier by keyword,
    robust to source-suffix variants (e.g. ".../CDP" vs ".../Environmental/CSR").

    Tiers (first match wins):
      - 'exact'                 : firm-reported exact value ('Exact Value from ...').
      - 'fuel'                  : firm-provided fuel/activity data converted via
                                  emission factors ('... from fuel use provided ...').
      - 'derived_deterministic' : firm-provided data deterministically transformed
                                  (derived/summed/split '... from data provided ...').
      - 'estimated'             : Trucost estimate, prior-year carry-forward, or
                                  chart/graph approximation.
      - 'none'                  : no coverage (NaN).
      - 'unclassified'          : recognized by none of the above (flagged by
                                  `assert_disclosure_coverage`).

    Note the 'fuel' check precedes 'from data provided': fuel-use strings read
    "... from fuel use provided ...", so they must not fall into the deterministic
    tier.
    """
    if pd.isna(sc1_disclosure):
        return "none"
    s = str(sc1_disclosure).lower()
    if "exact" in s:
        return "exact"
    if "fuel use" in s:
        return "fuel"
    if "from data provided" in s:
        return "derived_deterministic"
    if ("estimate" in s) or ("previous year" in s) or ("approximated" in s):
        return "estimated"
    return "unclassified"


# Which provenance tiers count as *disclosed* (co2_bool == 1) under each definition.
DISCLOSED_TIERS = {
    "narrow": {"exact"},                          # main spec: firm-reported exact values only
    "broad": {"exact", "derived_deterministic"},  # robustness: + deterministic transforms of provided data
}


def assert_disclosure_coverage(disclosure: pd.Series) -> None:
    """
    Raise if any non-null `sc1_disclosure` string is unrecognized by
    `emissions_tier`. Guards the keyword rules against silently sending a new
    Trucost label to 'Missing'.
    """
    unknown = disclosure[disclosure.map(emissions_tier) == "unclassified"]
    if len(unknown):
        vals = sorted(unknown.dropna().astype(str).unique())
        raise ValueError(
            f"[emissions_tier] {len(vals)} unrecognized sc1_disclosure label(s): "
            f"{vals}. Add a keyword rule in emissions_tier() before proceeding."
        )


def classify_emissions_status(sc1_disclosure: object, definition: str = "narrow") -> str:
    """
    Map a Trucost `sc1_disclosure` string to 'Observed' vs 'Missing' under the
    chosen disclosed-definition. Single source of truth for the disclosed/
    estimated split, `co2_bool`, and the paper's disclosure tables.

    definition:
      - 'narrow' (default, main spec): only firm-reported exact values are Observed.
      - 'broad'  (robustness): also count deterministic transforms of firm-provided
        data (derived/summed/split '... from data provided ...') as Observed.
    The 'fuel' and 'estimated' tiers, and no coverage, are Missing under both.
    """
    if definition not in DISCLOSED_TIERS:
        raise ValueError(
            f"definition must be one of {sorted(DISCLOSED_TIERS)}; got {definition!r}"
        )
    return "Observed" if emissions_tier(sc1_disclosure) in DISCLOSED_TIERS[definition] else "Missing"


def add_disclosure_splits(
    df: pd.DataFrame,
    *,
    disclosure_col: str = "sc1_disclosure",
    emissions_col: str = "sc1",
    out_disclosed: str = "sc1_disclosed",
    out_estimated: str = "sc1_estimated",
    definition: str = "narrow",
) -> pd.DataFrame:
    """
    Split Scope 1 emissions into disclosed vs estimated using the canonical
    `classify_emissions_status` treatment definition (see `definition`).

    Rows classified 'Observed' put their `emissions_col` value in `out_disclosed`;
    'Missing' rows put it in `out_estimated`. Raises if any `disclosure_col` label
    is unrecognized by the tier keywords.
    """
    _require_cols(df, [disclosure_col, emissions_col], name="add_disclosure_splits")
    out = df.copy()

    assert_disclosure_coverage(out[disclosure_col])
    is_disclosed = (
        out[disclosure_col].apply(lambda v: classify_emissions_status(v, definition)) == "Observed"
    )

    out[out_disclosed] = np.where(is_disclosed, out[emissions_col], np.nan)
    out[out_estimated] = np.where(is_disclosed, np.nan, out[emissions_col])
    return out


def add_response_indicator(
    df: pd.DataFrame,
    *,
    disclosed_col: str = "sc1_disclosed",
    out_col: str = "co2_bool",
) -> pd.DataFrame:
    """
    Add response indicator: 1 if disclosed emissions observed, 0 otherwise.
    """
    _require_cols(df, [disclosed_col], name="add_response_indicator")
    out = df.copy()
    out[out_col] = (~out[disclosed_col].isna()).astype(int)
    return out


def preprocess_panel(
    panel: pd.DataFrame,
    *,
    config: PanelPreprocessConfig = PanelPreprocessConfig(),
    sector_col: str = "gsector",
    assets_col: str = "at",
) -> pd.DataFrame:
    """
    Canonical post-merge preprocessing used throughout the repo.

    Steps:
      1) optional sector subsample
      2) split disclosed vs estimated emissions
      3) add response indicator co2_bool
      4) apply sample restrictions:
           - assets > 0 (if enabled)
           - drop negative disclosed emissions (if enabled)
    """
    _require_cols(panel, [sector_col, "sc1_disclosure", "sc1", assets_col], name="preprocess_panel")

    df = panel.copy()

    if config.subsample_sec is not None:
        df = df.loc[df[sector_col] == config.subsample_sec].copy()

    df = add_disclosure_splits(df, definition=config.disclosed_definition)
    df = add_response_indicator(df)

    if config.require_positive_assets:
        df = df.loc[df[assets_col] > 0].copy()

    if config.drop_negative_disclosed:
        df = df.loc[(df["sc1_disclosed"] > 0) | (df["sc1_disclosed"].isna())].copy()

    return df.reset_index(drop=True)
