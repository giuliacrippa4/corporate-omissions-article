from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from corporate_omissions.utils.paths import raw_dir


# -----------------------------
# Column maps (legacy -> canonical)
# -----------------------------

TRUCOST_RENAME = {
    "Fiscal Year": "year",
    "GV Key": "gvkey",
    "Ticker": "ticker",
    "Scope 1 Carbon Disclosure": "sc1_disclosure",
    "Absolute: Greenhouse Gases Scope 1": "sc1",
    "Intensity: GHG Scope 1": "int_sc1",
}

FUND_RENAME = {
    "fyear": "year",
    "tic": "ticker",
    # keep gvkey if present already; otherwise you merge via ticker (not ideal, but legacy)
}

GDP_RENAME = {
    "observation_date": "observation_date",
    "GDPC1": "gdp",
}

def _read_csv(path: Path, *, name: str, index_col: Optional[int] = None) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"[{name}] File not found: {path}")
    return pd.read_csv(path, index_col=index_col)


def _require(df: pd.DataFrame, cols: list[str], *, name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"[{name}] Missing required columns: {missing}\n"
            f"[{name}] Available columns: {list(df.columns)}"
        )


# -----------------------------
# Loaders
# -----------------------------

def load_trucost_raw(path: Optional[Path] = None) -> pd.DataFrame:
    """
    Loads the raw Trucost CSV. Your legacy export used index_col=0.
    """
    path = path or (raw_dir() / "trucost.csv")
    return _read_csv(path, name="trucost", index_col=0)


def load_fundamentals_raw(path: Optional[Path] = None) -> pd.DataFrame:
    """
    Loads the raw fundamentals CSV. Your legacy export used index_col=0.
    """
    path = path or (raw_dir() / "fundamentals.csv")
    return _read_csv(path, name="fundamentals", index_col=0)


def load_gdp_raw(path: Optional[Path] = None) -> pd.DataFrame:
    """
    Loads the raw GDP CSV (no index col in legacy).
    """
    path = path or (raw_dir() / "gdp.csv")
    return _read_csv(path, name="gdp", index_col=None)


def preprocess_trucost(trucost_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Legacy-consistent Trucost preprocessing:
      - select relevant columns
      - rename to canonical names
      - enforce year int
      - sort (year, gvkey)
    """
    _require(trucost_raw, list(TRUCOST_RENAME.keys()), name="trucost_raw")

    df = trucost_raw[list(TRUCOST_RENAME.keys())].copy()
    df = df.rename(columns=TRUCOST_RENAME)

    df["year"] = pd.to_numeric(df["year"], errors="raise").astype(int)
    df = df.sort_values(["year", "gvkey"]).reset_index(drop=True)

    return df


def preprocess_fundamentals(
    fund_raw: pd.DataFrame,
    *,
    filter_fic: bool = True,
    filter_exchg: bool = False,
) -> pd.DataFrame:
    """
    Legacy-consistent fundamentals preprocessing:
      - keep indfmt == 'INDL'
      - drop missing fyear
      - create year=int(fyear)
      - rename tic->ticker
      - drop metadata columns you previously removed

    Optional sample restrictions (applied only if the relevant column exists):
      - filter_fic: keep only fic == 'USA' (US-incorporated firms)
      - filter_exchg: keep only exchg in {11, 12, 14} (NYSE / AMEX / NASDAQ)
    """
    _require(fund_raw, ["indfmt", "fyear"], name="fundamentals_raw")

    fund = fund_raw.loc[fund_raw["indfmt"] == "INDL"].copy()
    fund = fund.dropna(subset=["fyear"]).reset_index()

    if filter_fic and "fic" in fund.columns:
        fund = fund.loc[fund["fic"] == "USA"].copy()

    if filter_exchg and "exchg" in fund.columns:
        fund = fund.loc[fund["exchg"].isin([11, 12, 14])].copy()

    fund["year"] = fund["fyear"].astype(int)

    if "tic" in fund.columns:
        fund = fund.rename(columns={"tic": "ticker"})

    # drop metadata columns (only if present)
    drop_cols = ["datadate", "indfmt", "consol", "popsrc", "datafmt", "curcd", "fyear", "fic"]
    fund = fund.drop(columns=[c for c in drop_cols if c in fund.columns], axis=1)

    return fund


def preprocess_gdp(gdp_raw: pd.DataFrame, *, log_gdp: bool = True) -> pd.DataFrame:
    """
    Legacy-consistent GDP preprocessing:
      - year from observation_date[:4]
      - log transform GDPC1 (renamed to gdp)
      - drop duplicate years (keep last)
    """
    _require(gdp_raw, ["observation_date", "GDPC1"], name="gdp_raw")

    gdp = gdp_raw.rename(columns=GDP_RENAME).copy()
    gdp["year"] = gdp["observation_date"].astype(str).str[:4].astype(int)

    gdp["gdp"] = pd.to_numeric(gdp["gdp"], errors="coerce")
    if log_gdp:
        gdp["gdp"] = gdp["gdp"].apply(lambda x: None if pd.isna(x) else __import__("math").log(x))

    gdp = gdp.drop_duplicates(subset=["year"], keep="last").reset_index(drop=True)

    # keep only what you need downstream
    return gdp[["year", "gdp"]]