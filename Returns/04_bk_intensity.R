# =============================================================================
# 04_bk_intensity.R
# Aswani-style carbon-premium test: same specification as the BK replication,
# but using EMISSIONS INTENSITY (emissions / denom) instead of raw levels.
#
# Aswani et al. (2024) argue the BK raw-emissions premium correlates
# mechanically with firm size. If no premium appears under intensity, the
# raw-level result was a size artifact.
#
# Two panels:
#   Panel A: intensity = sc1 / revt (Aswani's preferred specification)
#   Panel B: intensity = sc1 / at   (robustness with assets)
#
# Five specifications per panel:
#   (1) BK sample, observed intensity  (sc1_combined / denom)
#   (2) BK sample, MNAR intensity      (sc1_adj / denom)
#   (3) BK sample, MAR intensity       (sc1_disclosed_complete / denom)
#   (4) Full universe, MNAR intensity
#   (5) Full universe, MAR intensity
#
# Output:
#   returns_bk_intensity.tex
# =============================================================================

library(arrow)
source("00_setup.R")


# =============================================================================
# LOAD EMISSIONS PARQUET FOR REVT / AT DENOMINATORS AND MAR SERIES
# =============================================================================
cat("Loading emissions panel for revt and at...\n")
emis <- as.data.table(
  read_parquet("/Users/giulia/Documents/GitHub/corporate-omissions/data/processed/df_lm_avg_baseline_wti_gind.parquet")
)
emis[, gvkey := as.integer(gvkey)]

emis_fy <- unique(emis[, .(gvkey, year,
                           sc1_disclosed,
                           sc1_estimated,
                           sc1_adj,
                           sc1_disclosed_complete,
                           at, revt)])
emis_fy[, sc1_combined := fifelse(!is.na(sc1_disclosed), sc1_disclosed, sc1_estimated)]

safe_div <- function(num, den) {
  ifelse(is.na(num) | is.na(den) | den <= 0 | num <= 0, NA_real_, num / den)
}

emis_fy[, int_obs_rev  := safe_div(sc1_combined,           revt)]
emis_fy[, int_mnar_rev := safe_div(sc1_adj,                revt)]
emis_fy[, int_mar_rev  := safe_div(sc1_disclosed_complete, revt)]
emis_fy[, int_obs_at   := safe_div(sc1_combined,           at)]
emis_fy[, int_mnar_at  := safe_div(sc1_adj,                at)]
emis_fy[, int_mar_at   := safe_div(sc1_disclosed_complete, at)]

for (v in c("int_obs_rev","int_mnar_rev","int_mar_rev",
            "int_obs_at", "int_mnar_at", "int_mar_at")) {
  emis_fy[, (paste0("log_", v)) := log(get(v))]
}

intensity_cols <- c(
  "log_int_obs_rev","log_int_mnar_rev","log_int_mar_rev",
  "log_int_obs_at", "log_int_mnar_at", "log_int_mar_at")

# Merge intensities into the returns panel at gvkey-year level
if (!"gvkey" %in% names(dfm))
  stop("dfm is missing gvkey; check df_lm_baseline_rets.csv schema.")

dfm <- merge(
  dfm,
  emis_fy[, c("gvkey", "year", intensity_cols), with = FALSE],
  by = c("gvkey", "year"), all.x = TRUE, sort = FALSE)


# =============================================================================
# STANDARDISE INTENSITIES (full-sample and BK-subsample, both within month)
# =============================================================================
for (v in intensity_cols)
  dfm[!is.na(get(v)), paste0(v, "_z_full") := as.numeric(scale(get(v))), by = yyyymm]

bk <- dfm[group %in% c("disclosed", "estimated")]
for (v in intensity_cols)
  bk[!is.na(get(v)), paste0(v, "_z_sub") := as.numeric(scale(get(v))), by = yyyymm]

cat("Data ready:", nrow(dfm), "firm-months\n")


# =============================================================================
# REGRESSIONS — 5 specs per denominator
# =============================================================================
fit_and_extract <- function(data, var) {
  f <- as.formula(paste("ret ~", ctrl, "+", var, "| yyyymm + gsector"))
  fit <- feols(f, cluster = ~permno, data = data[!is.na(get(var))])
  ct  <- coeftable(fit)[var, ]
  list(coef = as.numeric(ct["Estimate"]),
       se   = as.numeric(ct["Std. Error"]),
       pval = as.numeric(ct["Pr(>|t|)"]),
       n    = fit$nobs,
       r2   = as.numeric(r2(fit, "r2")))
}

run_panel <- function(denom_label) {
  cat(sprintf("\n--- Panel: intensity = sc1 / %s ---\n", denom_label))
  list(
    bk_obs  = fit_and_extract(bk,  paste0("log_int_obs_",  denom_label, "_z_sub")),
    bk_mnar = fit_and_extract(bk,  paste0("log_int_mnar_", denom_label, "_z_sub")),
    bk_mar  = fit_and_extract(bk,  paste0("log_int_mar_",  denom_label, "_z_sub")),
    fu_mnar = fit_and_extract(dfm, paste0("log_int_mnar_", denom_label, "_z_full")),
    fu_mar  = fit_and_extract(dfm, paste0("log_int_mar_",  denom_label, "_z_full"))
  )
}
r_rev <- run_panel("rev")
r_at  <- run_panel("at")


# =============================================================================
# TABLE
# =============================================================================
row_coef <- function(label, r) paste0(label, " & ",
  fmt(r$bk_obs$coef),  stars(r$bk_obs$pval),  " & ",
  fmt(r$bk_mnar$coef), stars(r$bk_mnar$pval), " & ",
  fmt(r$bk_mar$coef),  stars(r$bk_mar$pval),  " & ",
  fmt(r$fu_mnar$coef), stars(r$fu_mnar$pval), " & ",
  fmt(r$fu_mar$coef),  stars(r$fu_mar$pval),  " \\\\")
row_se <- function(r) paste0(" & (",
  fmt(r$bk_obs$se),  ") & (", fmt(r$bk_mnar$se), ") & (",
  fmt(r$bk_mar$se),  ") & (", fmt(r$fu_mnar$se), ") & (",
  fmt(r$fu_mar$se),  ") \\\\")
row_n <- function(r) paste0("Observations & ",
  fmtN(r$bk_obs$n), " & ", fmtN(r$bk_mnar$n), " & ",
  fmtN(r$bk_mar$n), " & ", fmtN(r$fu_mnar$n), " & ", fmtN(r$fu_mar$n), " \\\\")
row_r2 <- function(r) paste0("$R^2$ & ",
  fmt(r$bk_obs$r2), " & ", fmt(r$bk_mnar$r2), " & ",
  fmt(r$bk_mar$r2), " & ", fmt(r$fu_mnar$r2), " & ", fmt(r$fu_mar$r2), " \\\\")

tex <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Premium under Emissions Intensity: Aswani-style Test.}",
  "Monthly panel regressions of stock returns on log-intensity of Scope~1 emissions,",
  "standardized within month. Following \\citet{aswani2024}, intensity removes",
  "the mechanical size correlation driving the raw-emissions specification of",
  "\\citet{bolton2021}. Panel A uses revenues as denominator;",
  "Panel B uses assets as robustness. Columns~(1)--(3) use the BK",
  "disclosed-plus-estimated subsample; Columns~(4)--(5) use the full Compustat",
  "universe. Standard errors clustered at the firm level.}",
  "\\label{tab:returns_bk_intensity}",
  "\\begin{tabular}{lccccc}", "\\toprule",
  " & (1) & (2) & (3) & (4) & (5) \\\\",
  " & BK Sample & BK Sample & BK Sample & Full Universe & Full Universe \\\\",
  " & Observed  & MNAR      & MAR       & MNAR          & MAR           \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Panel A: Intensity = Scope~1 / Revenues}} \\\\",
  row_coef("Log Emissions Intensity", r_rev),
  row_se(r_rev),
  row_n(r_rev),
  row_r2(r_rev),
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Panel B: Intensity = Scope~1 / Assets}} \\\\",
  row_coef("Log Emissions Intensity", r_at),
  row_se(r_at),
  row_n(r_at),
  row_r2(r_at),
  "\\midrule",
  "Controls    & Yes & Yes & Yes & Yes & Yes \\\\",
  "Time FE     & Yes & Yes & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes & Yes & Yes \\\\",
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")

writeLines(tex, file.path(outdir_tab, "returns_bk_intensity.tex"))
cat("Saved: returns_bk_intensity.tex\n")


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================
print_panel <- function(r, label) {
  cat(sprintf("\n=== %s ===\n", label))
  for (nm in names(r)) {
    x <- r[[nm]]
    s <- ifelse(x$pval < .01, "***", ifelse(x$pval < .05, "**",
         ifelse(x$pval < .10, "*", "")))
    cat(sprintf("  %-10s  coef=%+0.4f  SE=%.4f  p=%.3f  N=%7d  %s\n",
                nm, x$coef, x$se, x$pval, x$n, s))
  }
}
print_panel(r_rev, "Panel A: Scope 1 / Revenues")
print_panel(r_at,  "Panel B: Scope 1 / Assets")

cat("\n=== 04_bk_intensity.R COMPLETE ===\n")
