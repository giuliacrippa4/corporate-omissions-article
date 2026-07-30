from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm

from corporate_omissions.utils import SamplingUtils


# `iv_spec` lists what to EXCLUDE from the propensity (i.e. the excluded shifter).
#   - "year" is a special token: it drops the year fixed effects.
#   - any other entry names a DataFrame column to drop from the covariates,
#     e.g. "wti_crude", "gdp", or a shift-share instrument such as "z_wti".
# `gdp` / `wti_crude` are additionally auto-included when *not* excluded, for
# backwards compatibility with earlier specs.
YEAR_TOKEN = "year"
AUTO_ADD_IVS = ("gdp", "wti_crude")


@dataclass
class PropensityModel:
    estimated_params: Optional[pd.Series] = None
    CI: Optional[pd.DataFrame] = None
    df_full: Optional[pd.DataFrame] = None

    def _encode_features(
        self,
        df: pd.DataFrame,
        *,
        x1: Sequence[str],
        size_col: str,
        y_col: str,
        use_year_fe: bool,
        add_sector_fe: bool,
        use_green_sector: bool,
        interactions: bool,
        year_emission_interaction: bool = False,
        sector_col: str = "gsector",
        fit_columns: Optional[pd.Index] = None,
    ) -> Tuple[pd.DataFrame, pd.Index]:
        """
        Build design matrix X. If fit_columns provided, align to those columns (fill missing with 0).
        Returns (X, columns_used).
        """
        # basic checks
        for c in x1:
            if c not in df.columns:
                raise KeyError(f"Column '{c}' required in x1 but missing from DataFrame.")

        if use_green_sector and "green_sector" not in df.columns:
            raise KeyError("use_green_sector=True but 'green_sector' missing from DataFrame.")

        if use_green_sector:
            cols = list(x1) + ["green_sector"]
            X = df[cols].copy()

            if use_year_fe:
                if "year" not in df.columns:
                    raise KeyError("use_year_fe=True but 'year' missing from DataFrame.")
                year_d = pd.get_dummies(df["year"], prefix="year", drop_first=True, dtype="int")
                X = pd.concat([X.reset_index(drop=True), year_d.reset_index(drop=True)], axis=1)

        else:
            dummy_cols = []
            if add_sector_fe:
                if sector_col not in df.columns:
                    raise KeyError(f"add_sector_fe=True but '{sector_col}' missing from DataFrame.")
                dummy_cols.append(sector_col)
            if use_year_fe or year_emission_interaction:
                if "year" not in df.columns:
                    raise KeyError("year column is missing from DataFrame.")
                if "year" not in dummy_cols:
                    dummy_cols.append("year")

            df_enc = (
                pd.get_dummies(df, columns=dummy_cols, drop_first=True, dtype="int")
                if dummy_cols
                else df.copy()
            )

            fe_cols = []
            if add_sector_fe:
                fe_cols += [c for c in df_enc.columns if c.startswith(f"{sector_col}_")]
            if use_year_fe:
                fe_cols += [c for c in df_enc.columns if c.startswith("year_")]

            X = df_enc[list(x1) + fe_cols].copy()

            if interactions and add_sector_fe:
                gcols = [c for c in X.columns if c.startswith(f"{sector_col}_")]
                for col in gcols:
                    X[f"{size_col}_x_{col}"] = X[size_col] * X[col]
                    X[f"{y_col}_x_{col}"] = X[y_col] * X[col]

            if year_emission_interaction:
                ycols = [c for c in df_enc.columns if c.startswith("year_")]
                for col in ycols:
                    X[f"{y_col}_x_{col}"] = df_enc[y_col] * df_enc[col]

        X = sm.add_constant(X, has_constant="add")

        # align columns for prediction matrix
        if fit_columns is not None:
            X = X.reindex(columns=fit_columns, fill_value=0)
            return X, fit_columns

        return X, X.columns

    def fit_predict(
        self,
        df_unbalanced: pd.DataFrame,
        *,
        y_col: str = "sc1_disclosed",
        size_col: str = "at",
        want_balanced: bool = False,
        interactions: bool = False,
        year_emission_interaction: bool = False,
        print_summary: bool = False,
        clip_quantile: float = 0.01,
        base_x: Optional[Sequence[str]] = None,
        iv_spec: Optional[List[str]] = None,
        add_sector_fe: bool = True,
        sector_col: str = "gsector",
        logit_fit_method: str = "newton",
        add_plot: bool = False,
        log_transform: bool = False,
    ) -> pd.DataFrame:
        """
        Fit propensity model (possibly on balanced subsample), predict on full sample,
        and compute tilt-adjusted disclosures.

        Tilt restriction: depends ONLY on y_col via (const, beta_y).
        """
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

        use_year_fe = YEAR_TOKEN not in excluded
        if use_year_fe and "year" not in df_unbalanced.columns:
            raise KeyError("year FE requested but 'year' column is missing.")

        # ---------- estimation sample ----------
        if want_balanced:
            df_fit = SamplingUtils.balance_sample(df_unbalanced).copy().reset_index(drop=True)
        else:
            df_fit = df_unbalanced.copy()

        # ---------- optional log transform of size and y ----------
        # Tilt formula must use whatever scale the propensity model fits on.
        size_col_eff = size_col
        y_col_eff = y_col
        df_full_work = df_unbalanced.copy()
        if log_transform:
            log_size = f"log_{size_col}"
            log_y = f"log_{y_col}"
            for d in (df_fit, df_full_work):
                d[log_size] = np.log(d[size_col].clip(lower=1e-6).astype(float))
                d[log_y] = np.log1p(d[y_col].clip(lower=0).astype(float))
            size_col_eff = log_size
            y_col_eff = log_y

        # ---------- covariates ----------
        if base_x is None:
            x1 = [size_col_eff, y_col_eff]
        else:
            x1 = list(base_x)

        if y_col_eff not in x1:
            x1.append(y_col_eff)
        if size_col_eff not in x1:
            x1.insert(0, size_col_eff)

        # Auto-include the legacy macro controls unless excluded (and if present).
        for col in AUTO_ADD_IVS:
            if col not in excluded and col in df_fit.columns and col not in x1:
                x1.append(col)

        # Drop every excluded column from the propensity covariates.
        x1 = [c for c in x1 if c not in excluded]

        use_green_sector = "green_sector" in df_fit.columns

        # ---------- build X matrices with alignment ----------
        X_fit, fit_cols = self._encode_features(
            df_fit,
            x1=x1,
            size_col=size_col_eff,
            y_col=y_col_eff,
            use_year_fe=use_year_fe,
            add_sector_fe=add_sector_fe,
            use_green_sector=use_green_sector,
            interactions=interactions,
            year_emission_interaction=year_emission_interaction,
            sector_col=sector_col,
            fit_columns=None,
        )
        X_full, _ = self._encode_features(
            df_full_work,
            x1=x1,
            size_col=size_col_eff,
            y_col=y_col_eff,
            use_year_fe=use_year_fe,
            add_sector_fe=add_sector_fe,
            use_green_sector=use_green_sector,
            interactions=interactions,
            year_emission_interaction=year_emission_interaction,
            sector_col=sector_col,
            fit_columns=fit_cols,
        )

        # ---------- hygiene: exog cannot contain nan/inf ----------
        X_fit = X_fit.replace([np.inf, -np.inf], np.nan)
        X_full = X_full.replace([np.inf, -np.inf], np.nan)

        y_fit = df_fit["co2_bool"].astype(int)

        # drop rows with any missing regressor in estimation sample
        mask = X_fit.notna().all(axis=1)
        X_fit = X_fit.loc[mask]
        y_fit = y_fit.loc[mask]

        # drop zero-variance columns (often happens after balancing + FE)
        keep = (X_fit.var() > 0)
        if "const" in X_fit.columns:
            keep["const"] = True
        X_fit = X_fit.loc[:, keep]
        X_full = X_full.reindex(columns=X_fit.columns, fill_value=0)
        
        # ---------- fit logit ----------
        res = sm.Logit(y_fit, X_fit).fit(method=logit_fit_method, disp=False)

        if print_summary:
            print(res.summary())

        self.estimated_params = res.params
        self.CI = res.conf_int(alpha=0.1)

        # ---------- predict propensity ----------
        out = df_full_work.copy()
        out["propensity"] = res.predict(X_full)

        # ---------- tilt: restricted to y_col only ----------
        if y_col_eff not in self.estimated_params.index:
            raise KeyError(f"Model did not estimate coefficient for y_col='{y_col_eff}'. Ensure y_col is in x1/base_x.")

        if "const" in self.estimated_params.index:
            alpha_hat = float(self.estimated_params["const"])
        elif "Intercept" in self.estimated_params.index:
            alpha_hat = float(self.estimated_params["Intercept"])
        else:
            alpha_hat = 0.0  # fallback (shouldn't trigger if const preserved)

        beta_y_hat = float(self.estimated_params[y_col_eff])

        tilt_abs = np.exp(-alpha_hat - beta_y_hat * out[y_col_eff].astype(float))
        tilt_abs = pd.Series(tilt_abs, index=out.index)

        lo = tilt_abs.quantile(clip_quantile)
        hi = tilt_abs.quantile(1 - clip_quantile)
        tilt_abs = tilt_abs.clip(lower=lo, upper=hi)

        denom = tilt_abs[out["co2_bool"].astype(int) == 1].mean()
        tilt = tilt_abs / (denom if denom and np.isfinite(denom) else 1.0)

        out["tilt_absolute"] = tilt_abs
        out["tilt"] = tilt

        out["sc1_ours"] = out[y_col] * out["tilt"]
        out["sc1_adj"] = np.where(
            out["co2_bool"].astype(int) == 0,
            out[y_col] * out["tilt"],
            out[y_col],
        )

        self.df_full = out
        return out
