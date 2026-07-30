# preprocessing/empirical_preprocessing.py
"""
Empirical data preprocessing utilities for the Corporate-Omissions pipeline.

This module provides a small, pipeline-friendly preprocessor that:
- validates required columns,
- filters by year,
- winsorizes selected variables *within year*,
- (optionally) creates a green_sector indicator,
- stores convenient splits (observed vs not observed),
- (optionally) plots emissions vs assets.

Designed to be imported and called from your pipeline orchestration code.

Example
-------
from preprocessing.empirical_preprocessing import EmpiricalPreprocessor

pp = EmpiricalPreprocessor(required_columns=[...])
out = pp.fix_df_empirical(df, year_subsample=2010, plot_bool=False)

df_all = out["all"]
df_obs = out["raw_obs"]
df_not_obs = out["raw_not_obs"]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

# Plotting is optional; import lazily in the plotting function to keep deps light.


def _require_columns(df: pd.DataFrame, required: Sequence[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


@dataclass
class EmpiricalPreprocessor:
    """
    Preprocessor for empirical (panel) datasets.

    Notes
    -----
    - Assumes `year` exists and is numeric/int.
    - Assumes `gvkey` exists and identifies firms.
    - Assumes `co2_bool` is 0/1 for missing/observed emissions status.
    - If you want winsorization of additional vars (e.g., scope1, scope2,
      intensities), pass them via `vars_to_winsorize`.
    """

    required_columns: List[str]
    default_vars_to_winsorize: List[str] = field(default_factory=lambda: ["at"])
    winsor_p_low: float = 0.01
    winsor_p_high: float = 0.99

    # stored outputs (mirroring your previous pattern)
    all: Optional[pd.DataFrame] = None
    raw: Optional[pd.DataFrame] = None
    raw_obs: Optional[pd.DataFrame] = None
    raw_not_obs: Optional[pd.DataFrame] = None

    def fix_df_empirical(
        self,
        df: pd.DataFrame,
        year_subsample: Optional[int] = None,
        plot_bool: bool = True,
        create_green_sector: bool = True,
        required_columns: Optional[List[str]] = None,
        vars_to_winsorize: Optional[List[str]] = None,
    ) -> Dict[str, pd.DataFrame]:
        """
        Prepares the dataset for modeling by filtering, handling outliers, and
        optionally creating a green sector indicator.

        Parameters
        ----------
        df:
            Input dataframe.
        year_subsample:
            Minimum year to keep (inclusive). If None, keeps all years.
        plot_bool:
            Whether to plot emissions vs assets (lightweight diagnostics).
        create_green_sector:
            Whether to create a green_sector indicator from gsector membership.
        required_columns:
            Overrides `self.required_columns` if provided.
        vars_to_winsorize:
            Overrides `self.default_vars_to_winsorize` if provided.

        Returns
        -------
        dict
            Dictionary containing:
            - "all": processed full panel
            - "raw": alias to processed full panel
            - "raw_obs": co2_bool == 1
            - "raw_not_obs": co2_bool == 0
        """
        req_cols = required_columns if required_columns is not None else self.required_columns
        _require_columns(df, req_cols)

        # Keep only required columns
        dfp = df.loc[:, req_cols].copy()

        # Basic sanity
        if "gvkey" not in dfp.columns:
            raise ValueError("`gvkey` must be included in required_columns.")
        if "year" not in dfp.columns:
            raise ValueError("`year` must be included in required_columns.")
        if "co2_bool" not in dfp.columns:
            raise ValueError("`co2_bool` must be included in required_columns.")

        # Set index (firm id)
        dfp = dfp.set_index("gvkey", drop=True)

        # Debug length
        print(f"[EmpiricalPreprocessor] Data size (pre-filter): {len(dfp):,} rows")

        # Year filter
        if year_subsample is not None:
            dfp = dfp.loc[dfp["year"] >= year_subsample].copy()
            print(f"[EmpiricalPreprocessor] Data size (year >= {year_subsample}): {len(dfp):,} rows")

        # Drop outliers (rows outside [p_low, p_high] quantiles)
        wvars = vars_to_winsorize if vars_to_winsorize is not None else self.default_vars_to_winsorize
        dfp = self._drop_outliers(
            dfp,
            vars_to_filter=wvars,
            p_low=self.winsor_p_low,
            p_high=self.winsor_p_high,
        )
        print(f"[EmpiricalPreprocessor] Data size (after outlier removal): {len(dfp):,} rows")

        # Optional green sector feature
        if create_green_sector:
            if "gsector" not in dfp.columns:
                raise ValueError("create_green_sector=True requires `gsector` in required_columns.")
            brown_sectors = self._get_brown_sectors()
            dfp["green_sector"] = np.where(dfp["gsector"].isin(brown_sectors), 0, 1)

        # Store (mirroring your previous API)
        self.all = dfp.copy()
        self.raw = dfp.copy()

        self.raw_not_obs = self.raw.loc[self.raw["co2_bool"] == 0].copy()
        self.raw_obs = self.raw.loc[self.raw["co2_bool"] == 1].copy()

        # Optional plot
        if plot_bool:
            self._plot_emissions_vs_assets(self.raw)

        return {
            "all": self.all,
            "raw": self.raw,
            "raw_obs": self.raw_obs,
            "raw_not_obs": self.raw_not_obs,
        }

    @staticmethod
    def _drop_outliers(
        df: pd.DataFrame,
        vars_to_filter: Sequence[str],
        p_low: float = 0.01,
        p_high: float = 0.99,
    ) -> pd.DataFrame:
        """
        Drop rows with values outside [p_low, p_high] quantiles for specified variables.

        Notes
        -----
        - For each variable, removes rows outside the quantile range.
        - Rows with NaN in the filtered variable are kept.
        """
        out = df.copy()

        for v in vars_to_filter:
            if v not in out.columns:
                raise ValueError(f"Filter variable `{v}` not found in dataframe.")
            if not np.issubdtype(out[v].dtype, np.number):
                raise ValueError(f"Filter variable `{v}` must be numeric.")
            s = out[v]
            if s.notna().sum() == 0:
                continue
            ql = s.quantile(p_low)
            qh = s.quantile(p_high)
            out = out[(s >= ql) & (s <= qh) | s.isna()]

        return out

    @staticmethod
    def _plot_emissions_vs_assets(df: pd.DataFrame) -> None:
        """
        Lightweight diagnostic plot (assets vs emissions), if columns exist.

        Tries common column names:
        - assets: `at`
        - emissions: `scope1` or `scope1_emissions` or `co2` (first found)
        """
        # Lazy import to keep module import light in headless pipeline runs.
        import matplotlib.pyplot as plt

        if "at" not in df.columns:
            print("[EmpiricalPreprocessor] Plot skipped: missing `at`.")
            return

        emis_col = None
        for c in ["scope1", "scope1_emissions", "co2", "co2_emissions"]:
            if c in df.columns:
                emis_col = c
                break

        if emis_col is None:
            print("[EmpiricalPreprocessor] Plot skipped: no emissions column found.")
            return

        x = df["at"].astype(float)
        y = df[emis_col].astype(float)

        # Avoid crashing on all-missing
        valid = x.notna() & y.notna()
        if valid.sum() == 0:
            print("[EmpiricalPreprocessor] Plot skipped: no non-missing (at, emissions) pairs.")
            return

        plt.figure()
        plt.scatter(x[valid], y[valid], s=8, alpha=0.35)
        plt.xlabel("Total Assets (at)")
        plt.ylabel(emis_col)
        plt.title("Emissions vs Assets (raw, post-winsor)")
        plt.tight_layout()
        plt.show()

    @staticmethod
    def _get_brown_sectors() -> List[int]:
        """
        Returns your 'brown sectors' list.

        If you already have UsefulFunctions.get_brown_sectors(), this method will
        try to use it. Otherwise, replace the fallback list with your canonical one.
        """
        # Prefer repo canonical source if it exists.
        try:
            # Adjust import path to match your repo structure, if needed.
            from utils.useful_functions import UsefulFunctions  # type: ignore

            return list(UsefulFunctions.get_brown_sectors())
        except Exception:
            # Fallback (REPLACE with your canonical definition if you want this file standalone)
            # Keeping empty would silently label everything "green", which is dangerous.
            raise ImportError(
                "Could not import UsefulFunctions.get_brown_sectors(). "
                "Either (i) ensure utils/useful_functions.py exists with that function, "
                "or (ii) replace `_get_brown_sectors` fallback with your sector list."
            )
