# corporate_omissions/models/ols_lognormal_imputer.py
"""
OLS-based lognormal imputation for unobserved Scope 1 disclosures.

Concept
-------
Fit an OLS model on observed firms:
    log(y) = X beta + eps,   eps ~ Normal(0, sigma^2)
Then for unobserved rows compute:
    mu = X_unobs beta
    draw(s): y_draw = exp(mu + eps_draw)
Also provide closed-form summaries on original scale:
    median: exp(mu)
    mean:   exp(mu + 0.5*sigma^2)
    q_tau:  exp(mu + z_tau * sigma)

Design goals
------------
- No mutation of caller's DataFrame (copy internally).
- Stable dummy encoding across observed/unobserved via single design matrix.
- Explicit outputs (dict) + stored attributes for convenience.
- Do NOT overwrite the observed column by default; write imputation outputs to new columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm

try:
    from scipy.stats import norm  # for z-scores for quantiles
except Exception:  # pragma: no cover
    norm = None


@dataclass
class OLSLogNormalImputer:
    """
    OLS lognormal imputer for y_col.

    After calling `fit_simulate`, these attributes are populated:
      - model: fitted statsmodels OLSResults
      - estimated_params_ols: pd.Series
      - estimated_params_ols_CI: pd.DataFrame
      - sigma2: float (residual variance)
      - sigma: float (residual std dev)
      - design_columns_: list[str]
    """

    model: Optional[sm.regression.linear_model.RegressionResultsWrapper] = None
    estimated_params_ols: Optional[pd.Series] = None
    estimated_params_ols_CI: Optional[pd.DataFrame] = None
    sigma2: Optional[float] = None
    sigma: Optional[float] = None
    design_columns_: Optional[List[str]] = None

    def fit_simulate(
        self,
        df: pd.DataFrame,
        x_columns: List[str],
        seed: int = 123,
        print_summary: bool = False,
        co2_col: str = "co2_bool",
        y_col: str = "sc1_disclosed",
        # --- output columns ---
        out_prefix: Optional[str] = None,
        imputed_draw_col: Optional[str] = None,
        imputed_mean_col: Optional[str] = None,
        imputed_median_col: Optional[str] = None,
        imputed_sigma_col: Optional[str] = None,
        imputed_mu_col: Optional[str] = None,
        imputed_source_col: Optional[str] = None,
        imputed_method_col: Optional[str] = None,
        imputed_seed_col: Optional[str] = None,
        complete_col: Optional[str] = None,
        # --- distribution summaries ---
        quantiles: Sequence[float] = (0.05, 0.50, 0.95),
        # --- simulation controls ---
        n_draws: int = 1,
        draws_wide: bool = False,
        draws_col_prefix: Optional[str] = None,
        overwrite_y_col: bool = False,
        keep_y_nonpositive_unobs_as_nan: bool = True,
    ) -> Dict[str, object]:
        """
        Fit OLS on observed data and create imputation outputs for unobserved rows.

        Key behavior
        ------------
        - By default, does NOT overwrite `y_col` anywhere.
        - Writes new columns:
            * one (or many) stochastic draw(s)
            * mean / median / quantiles on original scale
            * provenance metadata
            * an explicit `complete` column for convenience

        Parameters
        ----------
        df : pd.DataFrame
            Must contain `co2_col`, `y_col`, and all columns in `x_columns`.
        x_columns : List[str]
            Covariates to include. If includes 'year' and/or 'gsector', they are dummy-encoded.
            If includes 'at', we use log(at) instead of at.
        seed : int
            Random seed for simulation noise.
        print_summary : bool
            Print statsmodels summary.
        co2_col : str
            Column indicating observed status (1 observed, 0 unobserved).
        y_col : str
            Dependent variable (must be positive for observed rows used in log-fit).
        out_prefix : Optional[str]
            Prefix for generated columns if explicit names not provided.
            Default: f"{y_col}_imp".
        quantiles : Sequence[float]
            Quantiles to compute for unobserved y on original scale (requires scipy).
        n_draws : int
            Number of Monte Carlo draws to generate for each unobserved row.
        draws_wide : bool
            If True and n_draws>1, output columns y_imp_draw_1 ... y_imp_draw_K (wide).
            If False, will store only the first draw in `imputed_draw_col`.
        draws_col_prefix : Optional[str]
            Prefix for draw columns when draws_wide=True. Default: f"{out_prefix}_draw".
        overwrite_y_col : bool
            If True, overwrite `y_col` for unobserved rows in df_full with the first draw.
            (Not recommended; default False.)
        keep_y_nonpositive_unobs_as_nan : bool
            If True, force y_col to NaN on unobserved rows in df_full (to avoid mixing).
            Default True.

        Returns
        -------
        dict with keys:
            - "df_full": DataFrame with added imputation columns (and optional complete col)
            - "df_obs": observed subset used for fitting
            - "df_unobs": unobserved subset with imputation columns
            - "model": fitted statsmodels model
            - "params": fitted coefficients
            - "ci": coefficient confidence intervals
            - "sigma2": residual variance
            - "sigma": residual std dev
            - "design_columns": list of design matrix columns used (after dummies + const)
        """
        if n_draws < 1:
            raise ValueError("n_draws must be >= 1")

        rng = np.random.default_rng(seed)

        # --------- column name defaults ----------
        if out_prefix is None:
            out_prefix = f"{y_col}_imp"

        if imputed_draw_col is None:
            imputed_draw_col = f"{out_prefix}_draw"
        if imputed_mean_col is None:
            imputed_mean_col = f"{out_prefix}_mean"
        if imputed_median_col is None:
            imputed_median_col = f"{out_prefix}_median"
        if imputed_sigma_col is None:
            imputed_sigma_col = f"{out_prefix}_sigma"
        if imputed_mu_col is None:
            imputed_mu_col = f"{out_prefix}_mu"
        if imputed_source_col is None:
            imputed_source_col = f"{out_prefix}_source"
        if imputed_method_col is None:
            imputed_method_col = f"{out_prefix}_method"
        if imputed_seed_col is None:
            imputed_seed_col = f"{out_prefix}_seed"
        if complete_col is None:
            complete_col = f"{y_col}_complete"

        if draws_col_prefix is None:
            draws_col_prefix = f"{out_prefix}_draw"

        # --------- validation ----------
        required = [co2_col, y_col] + list(x_columns)
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        d = df.copy()

        # --------- transforms ----------
        x_cols_effective = list(x_columns)
        if "at" in x_cols_effective:
            # log(at) requires positive
            if (d["at"] <= 0).any():
                d = d.loc[d["at"] > 0].copy()
            d["log_at"] = np.log(d["at"].astype(float))
            x_cols_effective = ["log_at" if c == "at" else c for c in x_cols_effective]

        # observed mask (fit rows)
        obs_mask = d[co2_col] == 1
        if not obs_mask.any():
            raise ValueError("No observed rows (co2_bool == 1) available to fit the OLS model.")

        # only positive observed y are usable for log-fit
        y_obs = d.loc[obs_mask, y_col]
        y_obs_pos = y_obs > 0
        keep_obs = obs_mask.copy()
        keep_obs.loc[obs_mask] = y_obs_pos.values
        # drop observed rows with non-positive y (log undefined)
        d = d.loc[~obs_mask | keep_obs].copy()
        obs_mask = d[co2_col] == 1

        if obs_mask.sum() == 0:
            raise ValueError("All observed rows had non-positive y; cannot fit log-OLS.")

        # define unobserved mask
        unobs_mask = d[co2_col] == 0
        if not unobs_mask.any():
            raise ValueError("No unobserved rows (co2_bool == 0) available to simulate.")

        # log(y) for observed
        d["_log_y"] = np.nan
        d.loc[obs_mask, "_log_y"] = np.log(d.loc[obs_mask, y_col].astype(float))

        # --------- design matrix with stable dummies ----------
        cat_cols = [c for c in ["year", "gsector", "ggroup", "gind", "gsubind"] if c in x_cols_effective]
        X_raw = d.loc[:, x_cols_effective].copy()

        if cat_cols:
            X_raw = pd.get_dummies(
                X_raw,
                columns=cat_cols,
                drop_first=True,
                dtype="int",
            )

        X_raw = sm.add_constant(X_raw, has_constant="add")

        X_obs = X_raw.loc[obs_mask]
        y = d.loc[obs_mask, "_log_y"].astype(float)

        # drop zero-variance columns in fit sample (e.g., dummies for categories
        # that appear only in unobserved rows → all-zero in X_obs → NaN coefs)
        keep = X_obs.var() > 0
        if "const" in X_obs.columns:
            keep["const"] = True
        X_obs = X_obs.loc[:, keep]

        X_unobs = X_raw.loc[unobs_mask]
        X_unobs = X_unobs.reindex(columns=X_obs.columns, fill_value=0)

        # --------- fit ----------
        model = sm.OLS(y, X_obs).fit()
        if print_summary:
            print(model.summary())

        resid = y - model.predict(X_obs)
        sigma2 = float(np.var(resid, ddof=0))
        sigma = float(np.sqrt(sigma2))

        # --------- predict mu for unobserved ----------
        mu_unobs = model.predict(X_unobs).astype(float)

        # --------- closed-form summaries on original scale ----------
        # median = exp(mu)
        y_median = np.exp(mu_unobs)
        # mean = exp(mu + 0.5*sigma2)
        y_mean = np.exp(mu_unobs + 0.5 * sigma2)

        # quantiles (requires scipy)
        q_map = {}
        if quantiles:
            if norm is None:
                raise ImportError(
                    "scipy is required to compute quantiles. "
                    "Install scipy or set quantiles=()."
                )
            for q in quantiles:
                if not (0.0 < float(q) < 1.0):
                    raise ValueError(f"Quantiles must be in (0,1). Got {q}.")
                z = float(norm.ppf(q))
                q_map[q] = np.exp(mu_unobs + z * sigma)

        # --------- simulate draw(s) ----------
        # eps ~ N(0, sigma^2)
        if n_draws == 1:
            eps = rng.normal(loc=0.0, scale=sigma, size=len(mu_unobs))
            y_draw_1 = np.exp(mu_unobs + eps)
            draws = None
        else:
            # shape (n_unobs, n_draws)
            eps_mat = rng.normal(loc=0.0, scale=sigma, size=(len(mu_unobs), n_draws))
            mu_arr = mu_unobs.to_numpy()
            y_draws_mat = np.exp(mu_arr[:, None] + eps_mat)
            y_draw_1 = y_draws_mat[:, 0]
            draws = y_draws_mat

        # --------- assemble outputs ----------
        df_obs = d.loc[obs_mask].copy()
        df_unobs = d.loc[unobs_mask].copy()

        # provenance + parameters
        df_unobs[imputed_mu_col] = mu_unobs.values
        df_unobs[imputed_sigma_col] = sigma
        df_unobs[imputed_median_col] = y_median.values
        df_unobs[imputed_mean_col] = y_mean.values
        df_unobs[imputed_draw_col] = y_draw_1

        # quantile columns
        for q, arr in q_map.items():
            qname = str(q).replace(".", "p")
            df_unobs[f"{out_prefix}_q{qname}"] = arr.values

        # multiple draws (wide)
        if n_draws > 1 and draws_wide and draws is not None:
            for k in range(n_draws):
                df_unobs[f"{draws_col_prefix}_{k+1}"] = draws[:, k]

        # provenance metadata
        df_unobs[imputed_source_col] = "imputed"
        df_obs[imputed_source_col] = "observed"

        df_unobs[imputed_method_col] = "ols_lognormal"
        df_obs[imputed_method_col] = "observed"

        df_unobs[imputed_seed_col] = int(seed)
        df_obs[imputed_seed_col] = pd.NA  # no seed for observed

        # concatenate
        df_full = pd.concat([df_obs, df_unobs], axis=0).sort_index()

        # ensure we don't accidentally “complete” y_col by leaving junk values in unobserved
        if keep_y_nonpositive_unobs_as_nan:
            df_full.loc[df_full[co2_col] == 0, y_col] = np.nan

        # explicit complete column (safe downstream default)
        df_full[complete_col] = df_full[y_col]
        df_full.loc[df_full[co2_col] == 0, complete_col] = df_full.loc[
            df_full[co2_col] == 0, imputed_draw_col
        ]

        # optional overwrite (explicitly requested)
        if overwrite_y_col:
            df_full.loc[df_full[co2_col] == 0, y_col] = df_full.loc[
                df_full[co2_col] == 0, imputed_draw_col
            ]

        # --------- store attributes ----------
        self.model = model
        self.estimated_params_ols = model.params
        self.estimated_params_ols_CI = model.conf_int()
        self.sigma2 = sigma2
        self.sigma = sigma
        self.design_columns_ = list(X_obs.columns)

        return {
            "df_full": df_full,
            "df_obs": df_obs,
            "df_unobs": df_unobs,
            "model": model,
            "params": model.params,
            "ci": model.conf_int(),
            "sigma2": sigma2,
            "sigma": sigma,
            "design_columns": list(X_obs.columns),
            "out_prefix": out_prefix,
            "columns": {
                "draw": imputed_draw_col,
                "mean": imputed_mean_col,
                "median": imputed_median_col,
                "mu": imputed_mu_col,
                "sigma": imputed_sigma_col,
                "source": imputed_source_col,
                "method": imputed_method_col,
                "seed": imputed_seed_col,
                "complete": complete_col,
            },
        }
