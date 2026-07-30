# =============================================================================
# 01_diagnostic.R
# Does MNAR correction change the return-emissions relationship, and does it
# reveal a premium that observed data miss?
#
# Tables:
#   returns_main.tex       — Table 1: Panels A/B/C (observed, MNAR, MAR)
#                            across disclosed / combined / full-sample columns
#   returns_timeseries.tex — Table 2: year-by-year horse race
#   returns_estimated.tex  — Appendix: vendor-estimated-only sample
# =============================================================================

source("00_setup.R")


# =============================================================================
# TABLE 1: MAIN DIAGNOSTIC
# Three columns (disclosed / combined / full) x three panels (obs / MNAR / MAR).
# Emissions re-standardised within the relevant subsample to avoid
# full-sample scaling artifacts.
# =============================================================================
cat("\n=== TABLE 1: MAIN DIAGNOSTIC ===\n")

dfm[!is.na(sc1_adj),
    sc1_adj_z_sub := as.numeric(scale(sc1_adj)), by = yyyymm]
dfm[!is.na(sc1_disclosed_filled),
    sc1_filled_z_sub := as.numeric(scale(sc1_disclosed_filled)), by = yyyymm]

# Disclosed sample
fit_d_obs  <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_disclosed_z  | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[group == "disclosed"])
fit_d_mnar <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_sub    | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[group == "disclosed"])
fit_d_mar  <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_filled_z_sub | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[group == "disclosed"])

# Combined sample (disclosed + estimated)
fit_c_obs  <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_combined_z   | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[group %in% c("disclosed", "estimated")])
fit_c_mnar <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_sub    | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[group %in% c("disclosed", "estimated")])
fit_c_mar  <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_filled_z_sub | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[group %in% c("disclosed", "estimated")])

# Full sample
dfm[!is.na(sc1_adj),
    sc1_adj_z_full := as.numeric(scale(sc1_adj)), by = yyyymm]
dfm[!is.na(sc1_disclosed_filled),
    sc1_filled_z_full := as.numeric(scale(sc1_disclosed_filled)), by = yyyymm]

fit_f_mnar <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_full    | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[!is.na(sc1_adj_z_full)])
fit_f_mar  <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_filled_z_full | yyyymm + gsector")),
                    cluster = ~permno, data = dfm[!is.na(sc1_filled_z_full)])

d_obs  <- extr(fit_d_obs,  "sc1_disclosed_z")
d_mnar <- extr(fit_d_mnar, "sc1_adj_z_sub")
d_mar  <- extr(fit_d_mar,  "sc1_filled_z_sub")
c_obs  <- extr(fit_c_obs,  "sc1_combined_z")
c_mnar <- extr(fit_c_mnar, "sc1_adj_z_sub")
c_mar  <- extr(fit_c_mar,  "sc1_filled_z_sub")
f_mnar <- extr(fit_f_mnar, "sc1_adj_z_full")
f_mar  <- extr(fit_f_mar,  "sc1_filled_z_full")

cat(sprintf("  Full sample MNAR: %s (%s) p=%.3f\n", fmt(f_mnar$coef), fmt(f_mnar$se), f_mnar$pval))
cat(sprintf("  Full sample MAR:  %s (%s) p=%.3f\n", fmt(f_mar$coef),  fmt(f_mar$se),  f_mar$pval))

tex1 <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Emissions and Stock Returns.}",
  "Panel regressions of monthly stock returns on standardized Scope~1 CO$_2$ emissions.",
  "All specifications include firm characteristics, time fixed effects, and industry fixed effects,",
  "following \\citet{bolton2021}.",
  "Columns~(1) and~(2) compare observed and corrected emissions on subsamples with available data.",
  "Column~(3) uses corrected emissions on the full Compustat universe.",
  "Panel~B uses MNAR-corrected emissions (\\texttt{sc1\\_adj}).",
  "Panel~C uses MAR-imputed emissions (\\texttt{sc1\\_disclosed\\_filled}: outcome model only).",
  "The comparison between Panels~B and~C isolates the contribution of the selection correction.",
  "Standard errors are clustered at the firm level.}",
  "\\label{tab:returns_main}",
  "\\begin{tabular}{lccc}", "\\toprule",
  " & (1) & (2) & (3) \\\\",
  " & Disclosed & Combined & Full Sample \\\\", "\\midrule",
  "\\textit{Panel A: Observed / vendor emissions} & & & \\\\[2pt]",
  paste0("Scope~1 Emissions & ",
         fmt(d_obs$coef), stars(d_obs$pval), " & ",
         fmt(c_obs$coef), stars(c_obs$pval), " & -- \\\\"),
  paste0(" & (", fmt(d_obs$se), ") & (", fmt(c_obs$se), ") & \\\\"),
  "",
  "\\textit{Panel B: MNAR-corrected emissions (this paper)} & & & \\\\[2pt]",
  paste0("Scope~1 Emissions & ",
         fmt(d_mnar$coef), stars(d_mnar$pval), " & ",
         fmt(c_mnar$coef), stars(c_mnar$pval), " & ",
         fmt(f_mnar$coef), stars(f_mnar$pval), " \\\\"),
  paste0(" & (", fmt(d_mnar$se), ") & (", fmt(c_mnar$se), ") & (", fmt(f_mnar$se), ") \\\\"),
  "",
  "\\textit{Panel C: MAR-imputed emissions (outcome model only, no tilt)} & & & \\\\[2pt]",
  paste0("Scope~1 Emissions & ",
         fmt(d_mar$coef), stars(d_mar$pval), " & ",
         fmt(c_mar$coef), stars(c_mar$pval), " & ",
         fmt(f_mar$coef), stars(f_mar$pval), " \\\\"),
  paste0(" & (", fmt(d_mar$se), ") & (", fmt(c_mar$se), ") & (", fmt(f_mar$se), ") \\\\"),
  "\\midrule",
  "Controls & Yes & Yes & Yes \\\\",
  "Time FE & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes \\\\",
  paste0("Obs (Panel A) & ", fmtN(d_obs$n), " & ", fmtN(c_obs$n), " & -- \\\\"),
  paste0("Obs (Panels B--C) & ", fmtN(d_mnar$n), " & ", fmtN(c_mnar$n), " & ",
         fmtN(f_mnar$n), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Returns are in percentage points.",
  "Emissions are standardized within month on the relevant subsample.",
  "Column~(3) has no Panel~A row because never-disclosing firms have no",
  "reported or vendor-estimated value by construction.",
  "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex1, file.path(outdir_tab, "returns_main.tex"))
cat("returns_main.tex written\n")


# =============================================================================
# TABLE 1b: EMISSIONS INTENSITY WITH A FIRM-REPORTED COLUMN
# Mirrors Table 1, but on log Scope-1 INTENSITY (emissions / revenue), the
# specification of Aswani et al. (2024). The point of this table is the
# firm-reported ("Disclosed") column: Aswani's claim is that any premium is a
# property of VENDOR-ESTIMATED values, so the firm-reported-vs-vendor contrast
# (Panel A, col 1 vs col 2) is the test. Revenue denominator only (Aswani's
# preferred); the assets robustness stays in 04_bk_intensity.R.
#
# NOTE: intensity construction is duplicated from 04_bk_intensity.R. If either
# is edited, keep them in sync (or factor into 00_setup.R).
# =============================================================================
cat("\n=== TABLE 1b: INTENSITY WITH FIRM-REPORTED COLUMN ===\n")

library(arrow)

emis_i <- as.data.table(read_parquet(
  "/Users/giulia/Documents/GitHub/corporate-omissions/data/processed/df_lm_avg_baseline_wti_gind.parquet"))
emis_i[, gvkey := as.integer(gvkey)]
emis_i <- unique(emis_i[, .(gvkey, year, sc1_disclosed, sc1_estimated,
                            sc1_adj, sc1_disclosed_complete, revt)])
emis_i[, sc1_combined := fifelse(!is.na(sc1_disclosed), sc1_disclosed, sc1_estimated)]

safe_div <- function(num, den)
  ifelse(is.na(num) | is.na(den) | den <= 0 | num <= 0, NA_real_, num / den)

emis_i[, log_int_disc := log(safe_div(sc1_disclosed,          revt))]  # firm-reported ONLY
emis_i[, log_int_comb := log(safe_div(sc1_combined,           revt))]  # + vendor estimates
emis_i[, log_int_mnar := log(safe_div(sc1_adj,                revt))]
emis_i[, log_int_mar  := log(safe_div(sc1_disclosed_complete, revt))]

int_cols <- c("log_int_disc", "log_int_comb", "log_int_mnar", "log_int_mar")
dfm <- merge(dfm, emis_i[, c("gvkey", "year", int_cols), with = FALSE],
             by = c("gvkey", "year"), all.x = TRUE, sort = FALSE)

# Standardise each intensity within month, on the subsample it is estimated on.
disc <- dfm[group == "disclosed"]
bk   <- dfm[group %in% c("disclosed", "estimated")]
for (v in int_cols) {
  disc[!is.na(get(v)), paste0(v, "_zd") := as.numeric(scale(get(v))), by = yyyymm]
  bk[  !is.na(get(v)), paste0(v, "_zb") := as.numeric(scale(get(v))), by = yyyymm]
  dfm[ !is.na(get(v)), paste0(v, "_zf") := as.numeric(scale(get(v))), by = yyyymm]
}

fit_i <- function(data, var) {
  fit <- feols(as.formula(paste("ret ~", ctrl, "+", var, "| yyyymm + gsector")),
               cluster = ~permno, data = data[!is.na(get(var))])
  extr(fit, var)
}

# Panel A: observed (firm-reported | vendor-inclusive | -- )
iA_d <- fit_i(disc, "log_int_disc_zd")
iA_c <- fit_i(bk,   "log_int_comb_zb")
# Panel B: MNAR
iB_d <- fit_i(disc, "log_int_mnar_zd")
iB_c <- fit_i(bk,   "log_int_mnar_zb")
iB_f <- fit_i(dfm,  "log_int_mnar_zf")
# Panel C: MAR
iC_d <- fit_i(disc, "log_int_mar_zd")
iC_c <- fit_i(bk,   "log_int_mar_zb")
iC_f <- fit_i(dfm,  "log_int_mar_zf")

cat(sprintf("  Obs firm-reported: %s (%s) p=%.3f\n", fmt(iA_d$coef), fmt(iA_d$se), iA_d$pval))
cat(sprintf("  Obs vendor-incl.:  %s (%s) p=%.3f\n", fmt(iA_c$coef), fmt(iA_c$se), iA_c$pval))

texI <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Emissions Intensity and Stock Returns: Firm-Reported versus Vendor.}",
  "Panel regressions of monthly returns on standardized log Scope~1 emissions",
  "\\emph{intensity} (emissions over revenue), the specification of",
  "\\citet{aswani2024carbon}. Column~(1) uses firm-reported emissions only;",
  "column~(2) adds vendor-estimated firms; column~(3) is the full Compustat",
  "universe. Panel~A uses observed values, Panels~B and~C the MNAR-corrected and",
  "MAR-imputed measures. Standard errors clustered at the firm level.}",
  "\\label{tab:returns_intensity}",
  "\\begin{tabular}{lccc}", "\\toprule",
  " & (1) & (2) & (3) \\\\",
  " & Disclosed & Combined & Full Sample \\\\", "\\midrule",
  "\\textit{Panel A: Observed / vendor intensity} & & & \\\\[2pt]",
  paste0("Log Emissions Intensity & ",
         fmt(iA_d$coef), stars(iA_d$pval), " & ",
         fmt(iA_c$coef), stars(iA_c$pval), " & -- \\\\"),
  paste0(" & (", fmt(iA_d$se), ") & (", fmt(iA_c$se), ") & \\\\"),
  "",
  "\\textit{Panel B: MNAR-corrected intensity (this paper)} & & & \\\\[2pt]",
  paste0("Log Emissions Intensity & ",
         fmt(iB_d$coef), stars(iB_d$pval), " & ",
         fmt(iB_c$coef), stars(iB_c$pval), " & ",
         fmt(iB_f$coef), stars(iB_f$pval), " \\\\"),
  paste0(" & (", fmt(iB_d$se), ") & (", fmt(iB_c$se), ") & (", fmt(iB_f$se), ") \\\\"),
  "",
  "\\textit{Panel C: MAR-imputed intensity (outcome model only, no tilt)} & & & \\\\[2pt]",
  paste0("Log Emissions Intensity & ",
         fmt(iC_d$coef), stars(iC_d$pval), " & ",
         fmt(iC_c$coef), stars(iC_c$pval), " & ",
         fmt(iC_f$coef), stars(iC_f$pval), " \\\\"),
  paste0(" & (", fmt(iC_d$se), ") & (", fmt(iC_c$se), ") & (", fmt(iC_f$se), ") \\\\"),
  "\\midrule",
  "Controls & Yes & Yes & Yes \\\\",
  "Time FE & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes \\\\",
  paste0("Obs (Panel A) & ", fmtN(iA_d$n), " & ", fmtN(iA_c$n), " & -- \\\\"),
  paste0("Obs (Panels B--C) & ", fmtN(iB_d$n), " & ", fmtN(iB_c$n), " & ", fmtN(iB_f$n), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Intensity is Scope~1 over revenue, logged and standardized",
  "within month on the relevant subsample. A firm-reported premium in column~(1)",
  "that is absent while the vendor-inclusive premium in column~(2) is present",
  "would indicate the premium is a property of vendor estimates, per",
  "\\citet{aswani2024carbon}. Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(texI, file.path(outdir_tab, "returns_intensity.tex"))
cat("returns_intensity.tex written\n")


# =============================================================================
# TABLE 2: YEAR-BY-YEAR HORSE RACE
# 2019 (peak ESG salience) is the key year: MNAR reveals a premium that
# observed data miss.
# =============================================================================
cat("\n=== TABLE 2: YEAR-BY-YEAR HORSE RACE ===\n")

years <- sort(unique(dfm$year))

ts_res <- rbindlist(lapply(years, function(y) {
  d <- dfm[year == y]

  fit_obs <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl, "+ sc1_combined_z | yyyymm + gsector")),
          cluster = ~permno,
          data = d[group %in% c("disclosed", "estimated") & !is.na(sc1_combined_z)]),
    error = function(e) NULL)

  d[!is.na(sc1_adj), sc1_adj_z_yr := as.numeric(scale(sc1_adj)), by = yyyymm]
  fit_mnar <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_yr | yyyymm + gsector")),
          cluster = ~permno, data = d[!is.na(sc1_adj_z_yr)]),
    error = function(e) NULL)

  d[!is.na(sc1_disclosed_filled),
    sc1_filled_z_yr := as.numeric(scale(sc1_disclosed_filled)), by = yyyymm]
  fit_mar <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl, "+ sc1_filled_z_yr | yyyymm + gsector")),
          cluster = ~permno, data = d[!is.na(sc1_filled_z_yr)]),
    error = function(e) NULL)

  get_res <- function(fit, var) {
    if (is.null(fit) || !var %in% rownames(coeftable(fit)))
      return(list(coef = NA_real_, se = NA_real_, pval = NA_real_, n = NA_integer_))
    co <- coeftable(fit)[var, ]
    list(coef = round(co["Estimate"], 3),
         se   = round(co["Std. Error"], 3),
         pval = co["Pr(>|t|)"],
         n    = fit$nobs)
  }
  obs  <- get_res(fit_obs,  "sc1_combined_z")
  mnar <- get_res(fit_mnar, "sc1_adj_z_yr")
  mar  <- get_res(fit_mar,  "sc1_filled_z_yr")

  data.table(year    = y,
             obs_c   = obs$coef,  obs_p  = obs$pval,  obs_se  = obs$se,
             mnar_c  = mnar$coef, mnar_p = mnar$pval, mnar_se = mnar$se,
             mar_c   = mar$coef,  mar_p  = mar$pval,  mar_se  = mar$se,
             n_obs   = obs$n,     n_full = mnar$n)
}))

print(ts_res[, .(year, obs_c, obs_p, mnar_c, mnar_p, mar_c, mar_p)])

tex2 <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Emissions and Returns: Year-by-Year.}",
  "Annual panel regressions of monthly returns on Scope~1 CO$_2$ emissions.",
  "Observed uses reported/vendor-combined emissions on the disclosed+estimated subsample.",
  "MNAR and MAR use \\texttt{sc1\\_adj} and \\texttt{sc1\\_disclosed\\_filled} respectively",
  "on the full sample. Emissions re-standardized within year and subsample.",
  "Controls, time and industry FEs. Standard errors clustered at firm level.}",
  "\\label{tab:returns_ts}",
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{Observed} & \\multicolumn{2}{c}{MNAR (this paper)} & \\multicolumn{2}{c}{MAR benchmark} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "Year & Coef. & (SE) & Coef. & (SE) & Coef. & (SE) \\\\", "\\midrule")

na_str <- function(x, p, s) {
  if (is.na(x)) return("-- & ")
  paste0(fmt(x), stars(p), " & (", fmt(s), ")")
}
for (i in seq_len(nrow(ts_res))) {
  r <- ts_res[i]
  tex2 <- c(tex2,
    paste0(r$year, " & ",
           na_str(r$obs_c,  r$obs_p,  r$obs_se),  " & ",
           na_str(r$mnar_c, r$mnar_p, r$mnar_se), " & ",
           na_str(r$mar_c,  r$mar_p,  r$mar_se),  " \\\\"))
}
tex2 <- c(tex2, "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} 2019 coincides with peak ESG institutional momentum.",
  "The divergence between Observed and MNAR columns in 2019 is direct evidence",
  "that the selection correction reveals a premium invisible to conventional methods.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex2, file.path(outdir_tab, "returns_timeseries.tex"))
cat("returns_timeseries.tex written\n")


# =============================================================================
# APPENDIX: VENDOR-ESTIMATED SAMPLE
# Large SE on MNAR reflects extrapolation outside the support of disclosing firms.
# =============================================================================
cat("\n=== APPENDIX: VENDOR-ESTIMATED SAMPLE ===\n")

dfest <- dfm[group == "estimated"]
dfest[!is.na(sc1_estimated), sc1_est_z_sub := as.numeric(scale(sc1_estimated)), by = yyyymm]
dfest[!is.na(sc1_adj),       sc1_adj_z_sub := as.numeric(scale(sc1_adj)),       by = yyyymm]

fit_ea <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_est_z_sub | yyyymm + gsector")),
                cluster = ~permno, data = dfest[!is.na(sc1_est_z_sub)])
fit_em <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_sub | yyyymm + gsector")),
                cluster = ~permno, data = dfest[!is.na(sc1_adj_z_sub)])
ea <- extr(fit_ea, "sc1_est_z_sub")
em <- extr(fit_em, "sc1_adj_z_sub")

tex_app <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Emissions and Returns: Vendor-Estimated Sample (Appendix).}",
  "Vendor-estimated emissions are Trucost model-based values for firms without",
  "firm-reported emissions. The large standard error on MNAR-corrected emissions",
  "reflects extrapolation outside the support of the disclosing sample.",
  "Standard errors clustered at firm level.}",
  "\\label{tab:returns_estimated}",
  "\\begin{tabular}{lcc}", "\\toprule",
  " & Vendor Emissions & MNAR-Corrected \\\\", "\\midrule",
  paste0("Scope~1 Emissions & ", fmt(ea$coef), stars(ea$pval), " & ",
         fmt(em$coef), stars(em$pval), " \\\\"),
  paste0(" & (", fmt(ea$se), ") & (", fmt(em$se), ") \\\\"),
  "\\midrule",
  "Controls & Yes & Yes \\\\", "Time FE & Yes & Yes \\\\",
  "Industry FE & Yes & Yes \\\\",
  paste0("Observations & ", fmtN(ea$n), " & ", fmtN(em$n), " \\\\"),
  paste0("$R^2$ & ", fmt(ea$r2), " & ", fmt(em$r2), " \\\\"),
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_app, file.path(outdir_tab, "returns_estimated.tex"))
cat("returns_estimated.tex written\n")

cat("\n=== 01_diagnostic.R COMPLETE ===\n")
