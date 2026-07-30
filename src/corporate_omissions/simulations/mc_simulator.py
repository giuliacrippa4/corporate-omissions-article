import numpy as np
import pandas as pd

from corporate_omissions.models.ols_lognormal_imputer import OLSLogNormalImputer
from corporate_omissions.models import PropensityEstimator

def run_mc(
    df_raw: pd.DataFrame,
    x_columns=("at", "year", "gsector", "gdp"),
    y_complete="sc1_disclosed_complete",
    n_mc=300,
    mc_seeds=None,
    want_balanced=True,
    balance_seed=0,   # keep constant across MC to avoid extra randomness
    lm_iv_spec=None,
    gam_iv_spec=None,
    store_pdep=False,
    simulated_cols_lm=None,
    simulated_cols_gam=None,
    year_emission_interaction=False,
    logit_fit_method="newton",
    lm_base_x=None,
    interaction_cutoff=2016,
    interaction_col_name="emissions_x_post",
    log_transform=False,
    sector_col="gsector",
    clip_quantile=0.01,
    run_gam=True,
):
    if mc_seeds is None:
        # deterministic list of seeds
        mc_seeds = list(range(1, n_mc + 1))

    records = []
    pdep_y_runs = []   # optional: store GAM y-smooth PD for each run (if you want)
    
    lm_sum = None
    gam_sum = None
    lm_base = None          # first df_lm  – carries constant columns
    gam_base = None
    _sim_lm = None           # resolved column list (set on first iter)
    _sim_gam = None
    n_sims = 0

    imp = OLSLogNormalImputer()

    for s in mc_seeds:
        # 1) impute (ONLY stochastic step)
        out_imp = imp.fit_simulate(
            df_raw,
            x_columns=list(x_columns),
            seed=s,
            n_draws=10,
            draws_wide=True,
            quantiles=(0.05, 0.50, 0.95),
            print_summary=False,
        )
        df_full = out_imp["df_full"]

        # recompute interaction columns if base_x references them
        if lm_base_x is not None and "year" in df_full.columns:
            post_dummy = (df_full["year"] >= interaction_cutoff).astype(int)
            df_full[interaction_col_name] = df_full[y_complete] * post_dummy

        # 2) LM propensity (deterministic given df_full + deterministic balancing)
        est_lm = PropensityEstimator(method="lm")
        lm_kwargs = dict(
            y_col=y_complete,
            size_col="at",
            iv_spec=lm_iv_spec,
            want_balanced=want_balanced,
            year_emission_interaction=year_emission_interaction,
            logit_fit_method=logit_fit_method,
            log_transform=log_transform,
            sector_col=sector_col,
            clip_quantile=clip_quantile,
        )
        if lm_base_x is not None:
            lm_kwargs["base_x"] = list(lm_base_x)
        df_lm = est_lm.fit_predict(
            df_full,
            **lm_kwargs,
        )

        # 3) GAM propensity (deterministic given df_full + deterministic balancing)
        if run_gam:
            est_gam = PropensityEstimator(method="gam")
            df_gam = est_gam.fit_predict(
                df_full,
                y_col=y_complete,
                size_col="at",
                iv_spec=gam_iv_spec,
                want_balanced=want_balanced,
                add_plot=store_pdep,
                sector_col=sector_col,
                clip_quantile=clip_quantile,
            )

        if simulated_cols_lm is not None:
            if lm_sum is None:
                _sim_lm = [c for c in simulated_cols_lm if c in df_lm.columns]
                lm_base = df_lm.copy()
                lm_sum = df_lm[_sim_lm].copy()
            else:
                lm_sum += df_lm[_sim_lm]

        if run_gam and simulated_cols_gam is not None:
            if gam_sum is None:
                _sim_gam = [c for c in simulated_cols_gam if c in df_gam.columns]
                gam_base = df_gam.copy()
                gam_sum = df_gam[_sim_gam].copy()
            else:
                gam_sum += df_gam[_sim_gam]

        n_sims += 1


        # 4) store scalars (full DataFrames are accumulated via lm_sum/gam_sum)
        rec = {
            "seed": s,
            "mean_complete": df_full[y_complete].mean(),
            "mean_adj_lm": df_lm["sc1_adj"].mean(),
            # sd of the regressor THIS draw's propensity was fit on. Needed to express
            # gamma per standard deviation: the coefficient multiplies y_complete, so it
            # must be scaled by sd(y_complete) for that draw -- not by the disclosed-only
            # sd (which is ~2x larger, since non-disclosers are smaller emitters) and not
            # by the sd of the MC-averaged panel (which is smoother than any single draw).
            # Scale per draw, then aggregate; sd varies materially across draws because
            # emissions are heavy-tailed in levels.
            "sd_y_draw": df_full[y_complete].std(),
        }

        if run_gam:
            rec["mean_adj_gam"] = df_gam["sc1_adj"].mean()

        # optional: store params / CIs when available
        rec["lm_params"] = est_lm.estimated_params
        rec["lm_ci"] = est_lm.CI

        if run_gam:
            rec['pdep'] = est_gam.model_.pdep if store_pdep else None
            rec['confi'] = est_gam.model_.confi if store_pdep else None
            rec['XX_dict'] = est_gam.model_.XX_dict if store_pdep else None

        records.append(rec)
    
     # --- build result ---
    result = {"records": pd.DataFrame(records)}

    if simulated_cols_lm is not None and lm_base is not None:
        df_lm_avg = lm_base
        df_lm_avg[_sim_lm] = lm_sum / n_sims
        result["df_lm_avg"] = df_lm_avg

    if run_gam and simulated_cols_gam is not None and gam_base is not None:
        df_gam_avg = gam_base
        df_gam_avg[_sim_gam] = gam_sum / n_sims
        result["df_gam_avg"] = df_gam_avg

    # return pd.DataFrame(records)
    return result
