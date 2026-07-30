# Corporate Omissions

This repository contains the replication code for the paper:

> **Corporate Omissions: Scope 1 Emissions Disclosure under Missing Not at Random Reporting**\
> Giulia Crippa, Florian Berg, Roberto Rigobon\
> [SSRN Working Paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5503418)

The project covers emissions measurement, disclosure bias analysis, and adjustment via exponential tilting and statistical reweighting. It is organized as a modular Python package (`src/corporate_omissions`) with a notebook-driven analysis pipeline.

---

## Repository Structure

```
corporate-omissions/
├── src/
│   └── corporate_omissions/       # Core Python package
│       ├── config.py              # Global configuration and column definitions
│       ├── data/                  # Data loading, preprocessing, merging
│       ├── models/                # Propensity estimators (LM, GAM) and imputation
│       ├── preprocessing/         # Filtering, winsorization, outlier removal
│       ├── simulations/           # Monte Carlo simulation framework
│       ├── visualization/         # Plotting and table generation
│       └── utils/                 # Path handling, constants, shared helpers
│
├── notebooks/                     # Analysis notebooks (pipeline stages)
├── Returns/                       # Stock returns analysis
├── tests/                         # Unit tests
├── data/
│   ├── raw/                       # Raw inputs (Trucost, fundamentals, GDP, energy)
│   ├── processed/                 # Intermediate datasets
│   └── outputs/                   # Figures and tables
│
├── environment.yml                # Minimal conda environment
├── environment.lock.yml           # Fully pinned environment (exact reproduction)
├── pyproject.toml                 # Packaging and dependencies
├── .gitignore
└── README.md
```

---

## Environment Setup

After cloning the repo:

```bash
conda env create -f environment.yml
conda activate corporate-omissions
```

To reproduce the exact environment used during development:

```bash
conda env create -f environment.lock.yml
conda activate corporate-omissions
```

---

## Data

Raw data files are included in the repository:

- `data/raw/trucost.csv` — Corporate Scope 1 emissions disclosures
- `data/raw/fundamentals.csv` — Company financials (assets, sector, revenue)
- `data/raw/gdp.csv` — Macroeconomic GDP data
- `data/raw/fred_energy.csv` — Energy prices (WTI crude, petroleum PPI)

All intermediate datasets, figures, and tables produced by the pipeline are written to:

- `data/processed/`
- `data/outputs/figures/`
- `data/outputs/tables/`

---

## Analysis Pipeline

The notebooks are numbered sequentially and should be run in order:

| Notebook | Description |
|----------|-------------|
| `00_pipeline_overview` | Environment validation and path setup |
| `01_load_raw` | Load raw datasets |
| `02_merge_panel` | Merge data sources into a unified panel |
| `03_preprocess_data` | Cleaning, winsorization, outlier removal, filtering |
| `04_imputation` | Outcome model, propensity (LM + GAM), Monte Carlo with exponential tilting |
| `04a_imputation_val` | Validation of imputation (rank correlation, prediction error by quantile) |
| `05_results_analysis` | Emissions bar charts and corporate carbon damages |
| `06_transitioning_companies` | Prediction accuracy for firms transitioning estimated → disclosed |
| `06a_oos_transition_validation` | Out-of-sample validation on transition firms |
| `07_build_returns_data` | Build the firm-month returns panel (WRDS/CRSP/Compustat) |
| `07a_returns_validation` | Carbon-premium regressions on never-disclosing firms |

### Running robustness checks (appendix)

`04_imputation.ipynb` is **parameterized via a `SPEC` dictionary** at the top. To reproduce
the baseline or any of the robustness specifications, edit the `SPEC` cell and re-run.
Pre-configured variants (shown commented out in the notebook):

| Spec name | Description |
|-----------|-------------|
| `baseline` | GDP instrument, 300 MC draws, `gind` (GICS industry) granularity |
| `wti` | WTI crude instrument, 300 MC draws, `gsector` granularity |
| `500mc` | 500 MC draws (convergence check), LM-only |
| `tilt_test` | Adds an emissions × post-2022 interaction to test whether regulatory pressure worsens MNAR bias |

All outputs are tagged with the spec name (e.g. `df_lm_avg_baseline.parquet`,
`response_box_plot_lm_baseline.pdf`) so different specs never overwrite each other.

Reusable logic lives in `src/corporate_omissions/`; notebooks serve as thin orchestration and presentation layers.

---

## Author

Giulia Crippa

