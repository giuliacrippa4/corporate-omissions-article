from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Any, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from pygam import LogisticGAM, s, f
from pygam.terms import Intercept

from corporate_omissions.utils import SamplingUtils


# `iv_spec` lists what to EXCLUDE from the propensity (i.e. the excluded shifter).
#   - "year" is a special token: it drops the year fixed effects.
#   - any other entry names a DataFrame column to exclude. The GAM only ever
#     builds its design from size, emissions, gdp, wti_crude and the fixed
#     effects, so any other named shifter (e.g. a shift-share "z_wti") is
#     excluded by construction.
YEAR_TOKEN = "year"


@dataclass
class PropensityGAM:
    """
    Semi-parametric disclosure propensity via LogisticGAM + exponential tilting.

    Propensity regressors:
      - s(size_col) + s(y_col)
      - optional: s(gdp) if iv_spec includes gdp
      - optional: sector FE (gsector_* or green_sector) if add_sector_fe=True
      - optional: year FE (year_*) if iv_spec includes year

    Tilt restriction:
      tilt_abs = exp( - intercept_hat - s_y_hat(y_i) )
      where s_y_hat is the partial dependence of the y smooth at each observation.
    """
    gam_model: Optional[LogisticGAM] = None
    df_full: Optional[pd.DataFrame] = None

    # Optional diagnostic objects (only populated if add_plot=True)
    pdep: Optional[Dict[int, Any]] = None
    confi: Optional[Dict[int, Any]] = None
    XX_dict: Optional[Dict[int, Any]] = None

    # saved pieces useful for MC
    intercept_hat: Optional[float] = None
    s_y: Optional[np.ndarray] = None
    fit_columns: Optional[pd.Index] = None
    spec_: Optional[dict] = None

    def _build_matrices(
        self,
        df_fit: pd.DataFrame,
        df_full: pd.DataFrame,
        *,
        y_col: str,
        size_col: str,
        use_year: bool,
        use_gdp: bool,
        use_wti: bool,
        add_sector_fe: bool,
        sector_col: str = "gsector",
    ) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
        """
        Build aligned X matrices for fit and full samples.
        Returns X_fit, y_fit, X_full (all with identical columns).
        """
        use_green_sector = "green_sector" in df_fit.columns

        base_cols = [size_col, y_col]
        if use_gdp:
            base_cols.append("gdp")
        if use_wti:
            base_cols.append("wti_crude")

        def _make_X(df: pd.DataFrame, cols: Sequence[str], dummy_cols: Sequence[str]) -> pd.DataFrame:
            X = df[list(cols)].copy()

            if use_green_sector:
                X["green_sector"] = df["green_sector"].astype(int)
            else:
                if dummy_cols:
                    d = pd.get_dummies(df[list(dummy_cols)], drop_first=False, dtype="int")
                    X = pd.concat([X.reset_index(drop=True), d.reset_index(drop=True)], axis=1)

            X = X.replace([np.inf, -np.inf], np.nan)
            return X

        dummy_cols = []
        if (not use_green_sector) and add_sector_fe:
            if sector_col not in df_fit.columns:
                raise KeyError(f"add_sector_fe=True requires '{sector_col}' column.")
            dummy_cols.append(sector_col)
        if use_year:
            if "year" not in df_fit.columns:
                raise KeyError("iv_spec uses year, but 'year' column is missing.")
            dummy_cols.append("year")

        X_fit = _make_X(df_fit, base_cols, dummy_cols)
        cols_fit = X_fit.columns
        X_full = _make_X(df_full, base_cols, dummy_cols).reindex(columns=cols_fit, fill_value=0)

        # drop rows in fit sample with any missing X
        mask = X_fit.notna().all(axis=1)
        X_fit = X_fit.loc[mask].reset_index(drop=True)
        y_fit = df_fit.loc[mask, "co2_bool"].astype(int).reset_index(drop=True)

        # keep for reproducibility
        self.fit_columns = cols_fit

        return X_fit, y_fit, X_full

    def fit_predict(
        self,
        df_unbalanced: pd.DataFrame,
        *,
        y_col: str = "sc1_disclosed",
        size_col: str = "at",
        want_balanced: bool = False,
        add_plot: bool = False,
        clip_quantile: float = 0.01,
        iv_spec: Optional[List[str]] = None,
        add_sector_fe: bool = True,
        sector_col: str = "gsector",
    ) -> pd.DataFrame:
        if "co2_bool" not in df_unbalanced.columns:
            raise KeyError("Target column 'co2_bool' missing from input DataFrame.")
        if y_col not in df_unbalanced.columns:
            raise KeyError(f"y_col='{y_col}' missing from input DataFrame.")
        if size_col not in df_unbalanced.columns:
            raise KeyError(f"size_col='{size_col}' missing from input DataFrame.")

        # ---------- IV spec (list of instruments to EXCLUDE from propensity) ----------
        if iv_spec is None:
            iv_spec = []
        excluded = set(iv_spec)

        # Every non-"year" entry must name an actual column (catches typos).
        unknown = {c for c in excluded - {YEAR_TOKEN} if c not in df_unbalanced.columns}
        if unknown:
            raise ValueError(
                f"iv_spec entries {sorted(unknown)} are neither the '{YEAR_TOKEN}' token "
                f"nor columns of the DataFrame."
            )

        use_year = YEAR_TOKEN not in excluded
        # Include the legacy macro controls only if not excluded and present.
        use_gdp = ("gdp" not in excluded) and ("gdp" in df_unbalanced.columns)
        use_wti = ("wti_crude" not in excluded) and ("wti_crude" in df_unbalanced.columns)

        df_full = df_unbalanced.copy()
        df_fit = SamplingUtils.balance_sample(df_full).copy().reset_index(drop=True) if want_balanced else df_full.copy().reset_index(drop=True)

        X_fit, y_fit, X_full = self._build_matrices(
            df_fit, df_full,
            y_col=y_col, size_col=size_col,
            use_year=use_year, use_gdp=use_gdp, use_wti=use_wti,
            add_sector_fe=add_sector_fe,
            sector_col=sector_col,
        )

        # -----------------------
        # Build GAM terms using column indices
        # -----------------------
        cols = list(X_fit.columns)

        def idx(name: str) -> int:
            return cols.index(name)

        # Core smooth terms
        terms = s(idx(size_col)) + s(idx(y_col))

        # Optional continuous smooth terms
        if use_gdp:
            terms += s(idx("gdp"))
        if use_wti:
            terms += s(idx("wti_crude"))

        # Factor terms
        use_green_sector = "green_sector" in X_fit.columns
        if use_green_sector:
            terms += f(idx("green_sector"))
        else:
            factor_cols = []
            if add_sector_fe:
                factor_cols += [c for c in cols if c.startswith(f"{sector_col}_")]
            if use_year:
                factor_cols += [c for c in cols if c.startswith("year_")]

            for c in factor_cols:
                terms += f(idx(c))

        # Fit
        gam = LogisticGAM(terms).fit(X_fit.values, y_fit.values)
        self.gam_model = gam

        # Predict on full sample
        propensity = gam.predict_proba(X_full.values)

        out = df_full.copy()
        out["propensity"] = propensity

        # Optional PD plots
        if add_plot:
            self.pdep, self.confi, self.XX_dict = {}, {}, {}
            for i, term in enumerate(gam.terms):
                if isinstance(term, Intercept):
                    continue
                XX = gam.generate_X_grid(term=i)
                pdep, confi = gam.partial_dependence(term=i, width=0.95)
                self.pdep[i] = pdep
                self.confi[i] = confi
                self.XX_dict[i] = XX

                #here

        # -----------------------
        # Tilt restricted to intercept + y smooth only
        # term indices: 0=s(size), 1=s(y), then optional s(gdp)/s(wti), then factors
        # -----------------------
        intercept_hat = float(gam.coef_[0])
        y_term_index = 1
        s_y = gam.partial_dependence(term=y_term_index, X=X_full.values)

        tilt_abs = np.exp(-intercept_hat - s_y)
        tilt_abs = pd.Series(tilt_abs, index=out.index)

        lo = tilt_abs.quantile(clip_quantile)
        hi = tilt_abs.quantile(1 - clip_quantile)
        tilt_abs = tilt_abs.clip(lower=lo, upper=hi)

        denom = tilt_abs[out["co2_bool"].astype(int) == 1].mean()
        tilt = tilt_abs / (denom if denom and np.isfinite(denom) else 1.0)

        out["tilt_semiparametric"] = tilt_abs
        out["tilt"] = tilt
        out["sc1_ours"] = out[y_col] * out["tilt"]
        out["sc1_adj"] = np.where(out["co2_bool"].astype(int) == 0, out[y_col] * out["tilt"], out[y_col])

        self.df_full = out
        self.intercept_hat = intercept_hat
        self.s_y = s_y
        self.spec_ = dict(
            model="GAM",
            y_col=y_col,
            size_col=size_col,
            iv_spec=iv_spec,
            want_balanced=want_balanced,
            add_sector_fe=add_sector_fe,
            clip_quantile=clip_quantile,
        )
        return out
