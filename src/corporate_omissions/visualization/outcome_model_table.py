from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


def _stars(p: float) -> str:
    if p < 0.01:
        return r"$^{***}$"
    if p < 0.05:
        return r"$^{**}$"
    if p < 0.1:
        return r"$^{*}$"
    return ""


def _fmt_coef(x: float) -> str:
    return f"{x:.4f}"


def _fmt_se(x: float) -> str:
    return f"({x:.3f})"


def write_text(path: str | Path, text: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def build_outcome_model_table_from_imputer(
    model,  # statsmodels results (OLSResults or robustcov results)
    *,
    caption: str,
    label: str,
    baseline_sector: str = "Comm. Services",
    baseline_year: int = 2012,
    gdp_term: str = "GDPC1",          # or "gdp" / "_gdp_z" depending on what you used in x_columns
    assets_term: str = "log_at",      # your imputer uses log_at when x_columns contains "at"
    sector_prefix: str = "gsector_",  # because you used pd.get_dummies(columns=["gsector"])
    year_prefix: str = "year_",       # because you used pd.get_dummies(columns=["year"])
    drop_fe_ses: bool = True,         # match your paper style: FE coef only, no SE lines
    year_min: Optional[int] = None,
    year_max: Optional[int] = None,
) -> str:
    """
    Builds the paper-style LaTeX table from the *imputer's* outcome model.

    Assumes your design matrix comes from:
      pd.get_dummies(..., columns=["year","gsector"], drop_first=True)
    so dummy names look like:
      year_2013, year_2014, ...
      gsector_Energy, gsector_Industrials, ...
    and intercept is "const".
    """
    params = pd.Series(model.params, index=model.model.exog_names)
    ses = pd.Series(model.bse, index=model.model.exog_names)
    pvals = pd.Series(model.pvalues, index=model.model.exog_names)

    nobs = int(model.nobs)
    r2 = float(getattr(model, "rsquared", np.nan))

    def line(term: str, name: str) -> str:
        if term not in params.index:
            return ""
        b = float(params[term])
        se = float(ses[term])
        p = float(pvals[term])
        return f"{name} & {_fmt_coef(b)}{_stars(p)} \\\\\n      & {_fmt_se(se)} \\\\\n"

    # main block
    body = ""
    body += r"\textit{Macroeconomic and firm characteristics} \\" + "\n"
    body += line(gdp_term, "GDPC1")  # prints label GDPC1 even if term is e.g. _gdp_z; adjust if needed
    body += line(assets_term, r"(log) Assets")
    body += line("const", "Constant")

    # sector FE
    sector_terms = [t for t in params.index if t.startswith(sector_prefix)]
    sector_terms = sorted(sector_terms, key=lambda t: t.replace(sector_prefix, ""))

    body += "\n" + r"\midrule" + "\n"
    body += rf"\textit{{Sector fixed effects (baseline: {baseline_sector})}} \\" + "\n"
    for t in sector_terms:
        sec = t.replace(sector_prefix, "")
        b = float(params[t]); p = float(pvals[t])
        body += f"{sec} & {_fmt_coef(b)}{_stars(p)} \\\\\n"
        if not drop_fe_ses:
            body += f"      & {_fmt_se(float(ses[t]))} \\\\\n"

    # year FE
    year_terms = [t for t in params.index if t.startswith(year_prefix)]
    def _year_num(t: str) -> int:
        return int(t.replace(year_prefix, ""))
    year_terms = sorted(year_terms, key=_year_num)

    if year_min is not None:
        year_terms = [t for t in year_terms if _year_num(t) >= year_min]
    if year_max is not None:
        year_terms = [t for t in year_terms if _year_num(t) <= year_max]

    body += "\n" + r"\midrule" + "\n"
    body += rf"\textit{{Year fixed effects (baseline: {baseline_year})}} \\" + "\n"
    for t in year_terms:
        yr = t.replace(year_prefix, "")
        b = float(params[t]); p = float(pvals[t])
        body += f"{yr} & {_fmt_coef(b)}{_stars(p)} \\\\\n"
        if not drop_fe_ses:
            body += f"      & {_fmt_se(float(ses[t]))} \\\\\n"

    footer = ""
    footer += "\n" + r"\midrule" + "\n"
    footer += f"Observations & {nobs:,} \\\\\n"
    if np.isfinite(r2):
        footer += f"$R^2$ & {r2:.3f} \\\\\n"

    latex = rf"""
\begin{{table}}[htbp]
\footnotesize
\centering
\caption{{{caption}}}
\label{{{label}}}
\begin{{tabular}}{{l c}}
\toprule
 & Log Scope~1 Emissions \\
\midrule
{body}{footer}\bottomrule
\end{{tabular}}

\begin{{flushleft}}
\footnotesize
\end{{flushleft}}
\end{{table}}
""".strip("\n")

    return latex
