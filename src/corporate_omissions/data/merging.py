from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import pandas as pd


GSIC_TO_SECTOR: Dict[float, str] = {
    10.0: "Energy",
    15.0: "Materials",
    20.0: "Industrials",
    25.0: "Consumer Discr",
    30.0: "Consumer Staples",
    35.0: "Health Care",
    40.0: "Financials",
    45.0: "IT",
    50.0: "Communication",
    55.0: "Utilities",
    60.0: "Real Estate",
}


@dataclass(frozen=True)
class MergeDiagnostics:
    n_fund: int
    n_trucost: int
    n_after_merge: int
    n_dupes_dropped: int
    share_rows_with_trucost: float
    share_rows_with_gdp: float
    year_min: Optional[int]
    year_max: Optional[int]


def _filter_years(df: pd.DataFrame, year_start: int, year_end: int, year_col: str = "year") -> pd.DataFrame:
    if year_col not in df.columns:
        raise ValueError(f"Expected '{year_col}' column, got columns={list(df.columns)}")
    out = df.copy()
    out = out[(out[year_col] >= year_start) & (out[year_col] <= year_end)]
    return out.reset_index(drop=True)


def _infer_merge_keys(fund: pd.DataFrame, trucost: pd.DataFrame) -> Tuple[list[str], str]:
    """
    Prefer (gvkey, year). If gvkey is missing in either, fall back to (ticker, year)
    if ticker exists in both. Otherwise fail loudly.
    """
    if {"gvkey", "year"}.issubset(fund.columns) and {"gvkey", "year"}.issubset(trucost.columns):
        return ["gvkey", "year"], "gvkey_year"
    if {"ticker", "year"}.issubset(fund.columns) and {"ticker", "year"}.issubset(trucost.columns):
        return ["ticker", "year"], "ticker_year"
    raise ValueError(
        "No compatible merge keys found. Need either (gvkey, year) in both or (ticker, year) in both.\n"
        f"fund columns: {list(fund.columns)}\n"
        f"trucost columns: {list(trucost.columns)}"
    )


def merge_panel(
    fund: pd.DataFrame,
    trucost: pd.DataFrame,
    gdp: pd.DataFrame,
    *,
    year_start: int,
    year_end: int,
    gsic_map: Dict[float, str] = GSIC_TO_SECTOR,
    fund_left_join_trucost: bool = True,
    verbose: bool = True,
) -> tuple[pd.DataFrame, MergeDiagnostics]:
    """
    Merge fundamentals + Trucost + GDP, replicating legacy behavior:

    - filter both panels to [year_start, year_end]
    - merge fundamentals (left) with trucost (right) using best available keys
      (prefer gvkey/year else ticker/year)
    - drop duplicates on (year, gvkey) if gvkey exists, else on merge keys
    - map gsector codes to names and drop missing sector
    - merge GDP on year

    Returns:
      merged_df, diagnostics
    """
    fund_f = _filter_years(fund, year_start, year_end, year_col="year")
    trucost_f = _filter_years(trucost, year_start, year_end, year_col="year")

    keys, key_mode = _infer_merge_keys(fund_f, trucost_f)

    how = "left" if fund_left_join_trucost else "inner"
    merged = fund_f.merge(trucost_f, how=how)

    # Determine dedupe keys (legacy: year+gvkey). If gvkey missing, dedupe on merge keys.
    if {"year", "gvkey"}.issubset(merged.columns):
        dedupe_keys = ["year", "gvkey"]
    else:
        dedupe_keys = keys

    n_before = len(merged)
    merged = merged.drop_duplicates(subset=dedupe_keys).reset_index(drop=True)
    n_after = len(merged)
    n_dupes_dropped = n_before - n_after

    # Map gsector if present
    if "gsector" in merged.columns:
        merged["gsector"] = pd.to_numeric(merged["gsector"], errors="coerce")
        merged["gsector"] = merged["gsector"].map(gsic_map)
        merged = merged[~merged["gsector"].isna()].reset_index(drop=True)

    # Merge GDP on year
    if "year" not in gdp.columns:
        raise ValueError(f"GDP must contain 'year'. Got columns={list(gdp.columns)}")

    # allow either 'gdp' (new canonical) or 'GDPC1' (legacy)
    gdp_col = "gdp" if "gdp" in gdp.columns else ("GDPC1" if "GDPC1" in gdp.columns else None)
    if gdp_col is None:
        raise ValueError("GDP must contain either 'gdp' or 'GDPC1' column.")

    merged = merged.merge(gdp[["year", gdp_col]], on="year", how="left")
    if gdp_col != "gdp":
        merged = merged.rename(columns={gdp_col: "gdp"})

    # Diagnostics
    # "has trucost" measure: share rows where sc1 present (preferred), else any non-null trucost column
    if "sc1" in merged.columns:
        share_trucost = float(merged["sc1"].notna().mean())
    elif "scope1_emissions" in merged.columns:
        share_trucost = float(merged["scope1_emissions"].notna().mean())
    else:
        # fallback: if any column from trucost was merged, approximate with non-null on last key
        share_trucost = float(merged[keys[0]].notna().mean())

    share_gdp = float(merged["gdp"].notna().mean())

    year_min = int(merged["year"].min()) if len(merged) else None
    year_max = int(merged["year"].max()) if len(merged) else None

    diag = MergeDiagnostics(
        n_fund=len(fund_f),
        n_trucost=len(trucost_f),
        n_after_merge=len(merged),
        n_dupes_dropped=n_dupes_dropped,
        share_rows_with_trucost=share_trucost,
        share_rows_with_gdp=share_gdp,
        year_min=year_min,
        year_max=year_max,
    )

    if verbose:
        print("=== Merge diagnostics ===")
        print(f"keys: {keys} (mode={key_mode})")
        print(f"fund rows (filtered):   {diag.n_fund:,}")
        print(f"trucost rows (filtered): {diag.n_trucost:,}")
        print(f"dedup keys: {dedupe_keys}")
        print(f"after merge + filters:  {diag.n_after_merge:,}")
        print(f"duplicates dropped:     {diag.n_dupes_dropped:,}")
        print(f"% rows with Trucost:    {diag.share_rows_with_trucost:.3f}")
        print(f"% rows with GDP:        {diag.share_rows_with_gdp:.3f}")
        print(f"year range:             {diag.year_min}–{diag.year_max}")

    return merged, diag
