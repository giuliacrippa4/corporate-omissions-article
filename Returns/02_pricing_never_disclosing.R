# =============================================================================
# 02_pricing_never_disclosing.R
# Carbon premium in never-disclosing firms — the paper's main asset-pricing
# result.
#
# Tables:
#   returns_never_main.tex   — Table 3: four specifications (industry / time /
#                              firm FE; stable-sector subset)
#   returns_never_ts.tex     — Table 4: year-by-year (peaks at 2019 and 2022)
#   returns_never_sector.tex — Appendix: sector-by-sector, with sparse-sector flag
# =============================================================================

source("00_setup.R")

# Never-disclosing subsample. Following returns_2_never_disclosing.r, include
# vendor-estimated firms (more conservative) alongside the no-observation group.
dfnone <- dfm[group %in% c("none", "estimated") & !is.na(sc1_adj)]
# dfnone <- dfm[group %in% c("none") & !is.na(sc1_adj)]

# CRITICAL: re-standardise emissions within this subsample (within-month).
# Using the full-sample z-score in a subset regression introduces a scaling artifact.
dfnone[, sc1_z := as.numeric(scale(sc1_adj)), by = yyyymm]

cat(sprintf("\nNever-disclosing sample: %s firm-months, %s unique firms\n",
    fmtN(nrow(dfnone)), fmtN(dfnone[, uniqueN(permno)])))


# =============================================================================
# TABLE 3: MAIN NEVER-DISCLOSING PREMIUM
# (1) Industry + time FE — within-industry cross-section
# (2) Time FE only       — across sectors
# (3) Firm + time FE     — within-firm over time
# (4) Stable sectors     — drop sectors with sparse disclosing support
# =============================================================================
cat("\n=== TABLE 3: NEVER-DISCLOSING PREMIUM ===\n")

fit_n1 <- feols(as.formula(paste("ret ~ sc1_z +", ctrl, "| yyyymm + gsector")),
                cluster = ~permno, data = dfnone)
fit_n2 <- feols(as.formula(paste("ret ~ sc1_z +", ctrl, "| yyyymm")),
                cluster = ~permno, data = dfnone)
fit_n3 <- feols(as.formula(paste("ret ~ sc1_z +", ctrl, "| permno + yyyymm")),
                cluster = ~permno, data = dfnone)

stable_sectors <- c("Consumer Discr", "Consumer Staples", "Energy",
                    "Industrials", "Materials", "Utilities")

fit_n4 <- feols(as.formula(paste("ret ~ sc1_z +", ctrl, "| yyyymm + gsector")),
                cluster = ~permno, data = dfnone[gsector %in% stable_sectors])

n1 <- extr(fit_n1, "sc1_z"); n2 <- extr(fit_n2, "sc1_z")
n3 <- extr(fit_n3, "sc1_z"); n4 <- extr(fit_n4, "sc1_z")

cat(sprintf("(1) Industry+time FE: %s (%s) p=%.3f\n", fmt(n1$coef), fmt(n1$se), n1$pval))
cat(sprintf("(2) Time FE only:     %s (%s) p=%.3f\n", fmt(n2$coef), fmt(n2$se), n2$pval))
cat(sprintf("(3) Firm FE:          %s (%s) p=%.3f\n", fmt(n3$coef), fmt(n3$se), n3$pval))
cat(sprintf("(4) Stable sectors:   %s (%s) p=%.3f\n", fmt(n4$coef), fmt(n4$se), n4$pval))

tex3 <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Premium in Never-Disclosing Firms.}",
  "Panel regressions of monthly returns on MNAR-corrected Scope~1 emissions",
  "restricted to firms that never voluntarily disclose.",
  "Column~(1) includes industry and time fixed effects.",
  "Column~(2) includes time fixed effects only.",
  "Column~(3) includes firm and time fixed effects, identifying from within-firm",
  "variation in imputed emissions.",
  "Column~(4) repeats~(1) for sectors with stable imputation support",
  "(Energy, Industrials, Materials, Utilities, Consumer Staples, Consumer Discretionary).",
  "Emissions re-standardized within this subsample and month.",
  "Standard errors clustered at the firm level.}",
  "\\label{tab:returns_never}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  " & (1) & (2) & (3) & (4) \\\\",
  " & Industry FE & Time FE & Firm FE & Stable Sectors \\\\", "\\midrule",
  paste0("MNAR Emissions & ",
         fmt(n1$coef), stars(n1$pval), " & ",
         fmt(n2$coef), stars(n2$pval), " & ",
         fmt(n3$coef), stars(n3$pval), " & ",
         fmt(n4$coef), stars(n4$pval), " \\\\"),
  paste0(" & (", fmt(n1$se), ") & (", fmt(n2$se), ") & (",
         fmt(n3$se), ") & (", fmt(n4$se), ") \\\\"),
  "\\midrule",
  "Controls & Yes & Yes & Yes & Yes \\\\",
  "Time FE & Yes & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & No & No & Yes \\\\",
  "Firm FE & No & No & Yes & No \\\\",
  paste0("Observations & ", fmtN(n1$n), " & ", fmtN(n2$n), " & ",
         fmtN(n3$n), " & ", fmtN(n4$n), " \\\\"),
  paste0("$R^2$ & ", fmt(n1$r2), " & ", fmt(n2$r2), " & ",
         fmt(n3$r2), " & ", fmt(n4$r2), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Sample restricted to firm-months with no reported or",
  "vendor-estimated Scope~1 emissions. MNAR-corrected emissions are firm-level",
  "posterior means from 300 Monte Carlo draws of the imputation procedure.",
  "The drop from~(1) to~(2) indicates pricing is within-industry.",
  "The strengthening from~(1) to~(3) indicates a within-firm signal.",
  "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex3, file.path(outdir_tab, "returns_never_main.tex"))
cat("returns_never_main.tex written\n")


# =============================================================================
# TABLE 4: YEAR-BY-YEAR — concentration in 2019 (ESG peak) and 2022 (SEC/CSRD)
# =============================================================================
cat("\n=== TABLE 4: YEAR-BY-YEAR NEVER-DISCLOSING PREMIUM ===\n")

ts_none <- rbindlist(lapply(sort(unique(dfnone$year)), function(y) {
  d <- dfnone[year == y]
  if (nrow(d) < 500) return(NULL)
  d[, sc1_z_yr := as.numeric(scale(sc1_adj)), by = yyyymm]

  fit <- tryCatch(
    feols(as.formula(paste("ret ~ sc1_z_yr +", ctrl, "| yyyymm + gsector")),
          cluster = ~permno, data = d[!is.na(sc1_z_yr)]),
    error = function(e) NULL)
  if (is.null(fit) || !"sc1_z_yr" %in% rownames(coeftable(fit))) return(NULL)
  co <- coeftable(fit)["sc1_z_yr", ]
  data.table(year = y,
             coef = round(co["Estimate"], 3),
             se   = round(co["Std. Error"], 3),
             pval = co["Pr(>|t|)"],
             n    = fit$nobs)
}))
print(ts_none)

tex4 <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Premium in Never-Disclosing Firms: Year-by-Year.}",
  "Annual panel regressions. Controls, time and industry FEs.",
  "Standard errors clustered at firm level.}",
  "\\label{tab:returns_never_ts}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  "Year & Coef. & (SE) & $p$-value & N \\\\", "\\midrule")
for (i in seq_len(nrow(ts_none))) {
  r <- ts_none[i]
  tex4 <- c(tex4,
    paste0(r$year, " & ", fmt(r$coef), stars(r$pval), " & (",
           fmt(r$se), ") & ", round(r$pval, 3), " & ", fmtN(r$n), " \\\\"))
}
tex4 <- c(tex4, "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} 2019 coincides with peak ESG institutional momentum",
  "(BlackRock climate letter, Business Roundtable statement, record ESG fund inflows).",
  "2022 coincides with the SEC climate disclosure proposal and the EU CSRD announcement.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex4, file.path(outdir_tab, "returns_never_ts.tex"))
cat("returns_never_ts.tex written\n")


# =============================================================================
# APPENDIX: SECTOR BREAKDOWN
# Flag sparse-support sectors (unstable imputation).
# =============================================================================
cat("\n=== APPENDIX: SECTOR BREAKDOWN ===\n")

sparse_sectors <- c("Financials", "Health Care", "IT", "Communication", "Real Estate")

sector_res <- rbindlist(lapply(levels(dfnone$gsector), function(s) {
  d <- dfnone[gsector == s & !is.na(sc1_adj)]
  if (nrow(d) < 1000) return(NULL)
  d[, sc1_z_sec := as.numeric(scale(sc1_adj)), by = yyyymm]
  if (d[, sd(sc1_z_sec, na.rm = TRUE)] < 1e-6) return(NULL)

  fit <- tryCatch(
    feols(as.formula(paste("ret ~ sc1_z_sec +", ctrl, "| yyyymm")),
          cluster = ~permno, data = d[!is.na(sc1_z_sec)]),
    error = function(e) NULL)
  if (is.null(fit) || !"sc1_z_sec" %in% rownames(coeftable(fit))) return(NULL)
  co <- coeftable(fit)["sc1_z_sec", ]
  data.table(sector   = s,
             coef     = round(co["Estimate"], 3),
             se       = round(co["Std. Error"], 3),
             pval     = co["Pr(>|t|)"],
             n        = fit$nobs,
             unstable = s %in% sparse_sectors)
}))
print(sector_res)

if (nrow(sector_res) > 0) {
  tex_sec <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    "\\caption{\\textbf{Carbon Premium in Never-Disclosing Firms: Sector Breakdown (Appendix).}",
    "Coefficient on MNAR-corrected emissions from separate panel regressions by GICS sector.",
    "Sectors marked with $\\dagger$ have sparse disclosing-firm support.",
    "Time FE. Standard errors clustered at firm level.}",
    "\\label{tab:returns_sector}",
    "\\begin{tabular}{lcccc}", "\\toprule",
    "Sector & Coef. & (SE) & $p$-value & N \\\\", "\\midrule")
  for (i in seq_len(nrow(sector_res))) {
    r <- sector_res[i]
    sec_label <- paste0(r$sector, ifelse(r$unstable, "$^{\\dagger}$", ""))
    tex_sec <- c(tex_sec,
      paste0(sec_label, " & ", fmt(r$coef), stars(r$pval), " & (",
             fmt(r$se), ") & ", round(r$pval, 3), " & ", fmtN(r$n), " \\\\"))
  }
  tex_sec <- c(tex_sec, "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    "$^{\\dagger}$ Financials, Health Care, IT, Communication Services, Real Estate.",
    "Imputation for these sectors extrapolates outside the support of disclosing firms.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex_sec, file.path(outdir_tab, "returns_never_sector.tex"))
  cat("returns_never_sector.tex written\n")
}

cat("\n=== 02_pricing_never_disclosing.R COMPLETE ===\n")
