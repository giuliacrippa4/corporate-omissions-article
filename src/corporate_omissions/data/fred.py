"""Pull energy price series from FRED and save to data/raw/."""

from pathlib import Path
from typing import Optional

import pandas as pd
from fredapi import Fred

from corporate_omissions.utils.paths import raw_dir


ENERGY_SERIES = {
    "wti_crude":       "DCOILWTICO",
    "henry_hub":       "DHHNGSP",
    "electricity_ind": "APU000072610",
    "ppi_petroleum":   "PCU324110324110",
}

DEFAULT_START = "2012-01-01"
DEFAULT_END = "2023-12-31"


def pull_energy_prices(
    api_key: str,
    start: str = DEFAULT_START,
    end: str = DEFAULT_END,
    series: Optional[dict] = None,
) -> pd.DataFrame:
    """
    Pull energy price series from FRED and collapse to annual means.

    Parameters
    ----------
    api_key : str
        FRED API key.
    start, end : str
        Date range (YYYY-MM-DD).
    series : dict, optional
        Mapping of {name: FRED_ID}. Defaults to ENERGY_SERIES.

    Returns
    -------
    pd.DataFrame with index=year, columns=series names.
    """
    fred = Fred(api_key=api_key)
    series = series or ENERGY_SERIES

    frames = []
    for name, sid in series.items():
        try:
            s = fred.get_series(sid, observation_start=start, observation_end=end)
            annual = s.dropna().resample("YE").mean()
            annual.index = annual.index.year
            annual.index.name = "year"
            frames.append(annual.rename(name))
            print(f"  pulled {name:25s} ({sid})")
        except Exception as e:
            print(f"  FAILED {name}: {e}")

    start_year = int(start[:4])
    end_year = int(end[:4])
    return pd.concat(frames, axis=1).loc[start_year:end_year]


def load_energy_prices(path: Optional[Path] = None) -> pd.DataFrame:
    """Load previously saved energy prices CSV from data/raw/."""
    path = path or (raw_dir() / "fred_energy.csv")
    if not path.exists():
        raise FileNotFoundError(
            f"Energy prices not found at {path}. "
            "Run pull_energy_prices() first and save with save_energy_prices()."
        )
    df = pd.read_csv(path, index_col="year")
    return df


def save_energy_prices(energy: pd.DataFrame, path: Optional[Path] = None) -> Path:
    """Save energy prices to data/raw/fred_energy.csv."""
    path = path or (raw_dir() / "fred_energy.csv")
    energy.to_csv(path)
    print(f"Saved energy prices to {path}")
    return path
