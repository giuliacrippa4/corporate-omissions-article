# =============================================================================
# 00_setup.R
# Shared data loader, variable construction, and helper functions.
# Source at the top of every script in this folder.
# =============================================================================

library(data.table)
library(zoo)
library(fixest)

setwd("/Users/giulia/Documents/GitHub/corporate-omissions/Returns")
outdir_tab <- "../data/outputs/tables"
outdir_fig <- "../data/outputs/figures"
dir.create(outdir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# LOAD MONTHLY RETURNS PANEL
# -----------------------------------------------------------------------------
df <- fread("/Users/giulia/Documents/GitHub/corporate-omissions/data/processed/df_lm_baseline_wti_gind_rets_ff_dec_exchg.csv")

df[, mdate   := as.Date(mdate)]
df[, yyyymm  := as.yearmon(yyyymm)]
setkey(df, permno, mdate)

# Emissions variables:
#   sc1_disclosed         — firm-reported (NA for non-disclosers)
#   sc1_estimated         — Trucost vendor estimate (NA if no vendor coverage)
#   sc1_combined          — disclosed if present else estimated
#   sc1_adj               — disclosed (where available) + MNAR-imputed for non-disclosers
#   sc1_disclosed_filled  — disclosed + MAR outcome-model for non-disclosers
#   sc1_ours              — MNAR-imputed for ALL firms (replaces even disclosed; do not use in main specs)
df[, sc1_combined := fifelse(!is.na(sc1_disclosed), sc1_disclosed, sc1_estimated)]
df[, log_me       := log(pmax(me, 1e-6))]
# df[, ret := ret/100]

market_vars <- c("beta12", "mom_12m_excl1", "vol_12m", "log_me")
acct_vars   <- c("log_ppe", "salesgr", "epsgr", "roe", "bm",
                 "hhi", "invest_a", "leverage")
sc1_vars    <- c("sc1_disclosed", "sc1_estimated",
                 "sc1_adj", "sc1_disclosed_filled",
                 "sc1_ours", "sc1_combined")

keep_cols <- intersect(
  c("permno", "gvkey", "fyear", "mdate", "yyyymm", "gsector", "exchg", "ret",
    market_vars, acct_vars, sc1_vars),
  names(df))

dfm <- df[, ..keep_cols][!is.na(ret) & !is.na(yyyymm)]
dfm <- dfm[ret <= 100]   # drop extreme return outliers, following B&K
dfm[, gsector := as.factor(gsector)]
dfm[, year    := as.integer(format(mdate, "%Y"))]

# -----------------------------------------------------------------------------
# WINSORIZATION HELPERS
# -----------------------------------------------------------------------------
winsor_05 <- function(x, p = 0.005) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2])
}
winsor_25 <- function(x, p = 0.025) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE); pmin(pmax(x, q[1]), q[2])
}

# -----------------------------------------------------------------------------
# LAG CONTROLS (monthly market vars, annual accounting vars)
# -----------------------------------------------------------------------------
for (v in market_vars)
  dfm[, paste0(v, "_L1") := shift(get(v), 1L, "lag"), by = permno]

yc <- dfm[, c("permno", "year", acct_vars), with = FALSE][
  , lapply(.SD, function(z) z[1L]), by = .(permno, year)]
for (v in acct_vars)
  yc[, paste0(v, "_L1") := shift(get(v), 1L, "lag"), by = permno]
dfm <- merge(dfm,
             yc[, c("permno", "year", paste0(acct_vars, "_L1")), with = FALSE],
             by = c("permno", "year"), all.x = TRUE, sort = FALSE)

# -----------------------------------------------------------------------------
# STANDARDISE CONTROLS (winsorise first, then within-month z-score)
# -----------------------------------------------------------------------------
for (v in c("salesgr", "epsgr", "mom_12m_excl1", "vol_12m")) {
  vL1 <- paste0(v, "_L1"); vz <- paste0(v, "_L1_z")
  dfm[, (vL1) := winsor_05(get(vL1)), by = yyyymm]
  dfm[, (vz)  := as.numeric(scale(get(vL1))), by = yyyymm]
}
for (v in c("bm", "roe", "invest_a", "leverage")) {
  vL1 <- paste0(v, "_L1"); vz <- paste0(v, "_L1_z")
  dfm[, (vL1) := winsor_25(get(vL1)), by = yyyymm]
  dfm[, (vz)  := as.numeric(scale(get(vL1))), by = yyyymm]
}
for (v in c("hhi", "log_me", "log_ppe", "beta12")) {
  vz <- paste0(v, "_L1_z")
  dfm[, (vz) := as.numeric(scale(get(paste0(v, "_L1")))), by = yyyymm]
}

# -----------------------------------------------------------------------------
# STANDARDISE EMISSIONS WITHIN MONTH (full-sample z).
# For subsample regressions re-standardise within the subsample — using the
# full-sample z-score in a subset introduces a scaling artifact.
# -----------------------------------------------------------------------------
for (v in sc1_vars)
  if (v %in% names(dfm))
    dfm[, paste0(v, "_z") := as.numeric(scale(get(v))), by = yyyymm]

# -----------------------------------------------------------------------------
# DISCLOSURE GROUP
# -----------------------------------------------------------------------------
dfm[, group := fifelse(!is.na(sc1_disclosed_z), "disclosed",
               fifelse(!is.na(sc1_estimated_z), "estimated", "none"))]

ctrl <- paste(
  "beta12_L1_z + mom_12m_excl1_L1_z + vol_12m_L1_z + log_me_L1_z +",
  "invest_a_L1_z + hhi_L1_z + leverage_L1_z + log_ppe_L1_z +",
  "roe_L1_z + bm_L1_z + salesgr_L1_z + epsgr_L1_z"
)

# -----------------------------------------------------------------------------
# OUTPUT FORMATTING HELPERS
# -----------------------------------------------------------------------------
stars <- function(p)
  ifelse(is.na(p), "",
  ifelse(p < .01, "$^{***}$",
  ifelse(p < .05, "$^{**}$",
  ifelse(p < .1,  "$^{*}$", ""))))
fmt  <- function(x) formatC(x, format = "f", digits = 3)
fmtN <- function(x) formatC(x, format = "d", big.mark = ",")
extr <- function(fit, var) {
  co <- coeftable(fit)[var, ]
  list(coef = round(co["Estimate"], 3),
       se   = round(co["Std. Error"], 3),
       pval = co["Pr(>|t|)"],
       n    = fit$nobs,
       r2   = round(r2(fit, "r2"), 3))
}
# =============================================================================
# SAMPLE COMPOSITION DIAGNOSTICS (for Section 4.1 of paper)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. UNIQUE FIRM COUNTS BY DISCLOSURE HISTORY
#    Answers: "What % of unique firms ever disclose? Are 'never disclosers' a
#    persistent population or transitional?"
# -----------------------------------------------------------------------------
firm_history <- dfm[, .(
  n_months         = .N,
  n_disclosed      = sum(group == "disclosed"),
  n_estimated      = sum(group == "estimated"),
  n_none           = sum(group == "none"),
  ever_disclosed   = any(group == "disclosed"),
  ever_estimated   = any(group == "estimated"),
  always_none      = all(group == "none")
), by = permno]

firm_history[, history_type := fcase(
  ever_disclosed,                       "ever_disclosed",
  !ever_disclosed & ever_estimated,     "trucost_only",
  always_none,                          "never_in_trucost"
)]

cat("\n--- Unique firm counts by disclosure history ---\n")
print(firm_history[, .(
  n_firms         = .N,
  share_of_firms  = round(.N / nrow(firm_history), 3),
  total_months    = sum(n_months),
  share_of_months = round(sum(n_months) / nrow(dfm), 3)
), by = history_type][order(-n_firms)])

cat("\nTotal unique firms:", nrow(firm_history), "\n")
cat("Firms that never voluntarily disclose (trucost_only + never_in_trucost):",
    firm_history[history_type != "ever_disclosed", .N],
    "(", round(firm_history[history_type != "ever_disclosed", .N] / nrow(firm_history) * 100, 1), "% )\n")

# -----------------------------------------------------------------------------
# 2. FIRMS PER MONTH / YEAR OVER TIME
#    Answers: "Is the panel stable, expanding, or contracting?"
#    Detects coverage breaks (like the 2009 jump in your earlier annual table).
# -----------------------------------------------------------------------------
firms_by_month <- dfm[, .(n_firms = uniqueN(permno)), by = yyyymm][order(yyyymm)]
firms_by_year  <- dfm[, .(n_firms = uniqueN(permno),
                          n_obs   = .N), by = year][order(year)]

cat("\n--- Unique firms per year ---\n")
print(firms_by_year)

cat("\nMean firms per month:", round(mean(firms_by_month$n_firms), 0), "\n")
cat("Range:", range(firms_by_month$n_firms), "\n")

# -----------------------------------------------------------------------------
# 3. DISCLOSURE COMPOSITION BY YEAR
#    Answers: "Are disclosure rates rising over time? At what rate?"
#    For Figure 1 (top panel) consistency check.
# -----------------------------------------------------------------------------
disclosure_by_year <- dfm[, .(
  n_obs      = .N,
  disclosed  = sum(group == "disclosed"),
  estimated  = sum(group == "estimated"),
  none       = sum(group == "none")
), by = year][order(year)]

disclosure_by_year[, `:=`(
  pct_disclosed = round(disclosed / n_obs * 100, 1),
  pct_estimated = round(estimated / n_obs * 100, 1),
  pct_none      = round(none / n_obs * 100, 1)
)]

cat("\n--- Disclosure composition by year ---\n")
print(disclosure_by_year)

# -----------------------------------------------------------------------------
# 4. BOLTON-KACPERCZYK COMPARABILITY
#    Answers: "How does our sample compare to BK's working sample?"
#    BK ≈ disclosed + estimated (Trucost-covered firms only).
#    Ours adds 'none' (Compustat firms with no Trucost coverage).
# -----------------------------------------------------------------------------
bk_equivalent     <- dfm[group %in% c("disclosed", "estimated"), .N]
ours_extra        <- dfm[group == "none", .N]
total             <- nrow(dfm)

cat("\n--- Sample size relative to Bolton-Kacperczyk universe ---\n")
cat("BK-equivalent (Trucost-covered, disclosed + estimated):", fmtN(bk_equivalent), "firm-months\n")
cat("Our additional firms (Compustat, no Trucost):          ", fmtN(ours_extra),    "firm-months\n")
cat("Total sample:                                          ", fmtN(total),         "firm-months\n")
cat("Ratio (ours / BK-equivalent):                          ", round(total / bk_equivalent, 2), "x\n")

bk_firms   <- dfm[group %in% c("disclosed", "estimated"), uniqueN(permno)]
ours_firms <- dfm[, uniqueN(permno)]
cat("Unique firms BK-equivalent:", bk_firms, "\n")
cat("Unique firms ours:         ", ours_firms, "\n")
cat("Ratio (firm count):        ", round(ours_firms / bk_firms, 2), "x\n")

# -----------------------------------------------------------------------------
# 5. SECTOR COMPOSITION OF NEVER-DISCLOSING POPULATION
#    Answers: "Are the firms we add to BK's sample concentrated in any sector?
#    Are they mostly small caps?"
# -----------------------------------------------------------------------------
sector_by_group <- dfm[, .(n = .N), by = .(gsector, group)]
sector_by_group_wide <- dcast(sector_by_group, gsector ~ group, value.var = "n", fill = 0)

cat("\n--- Sector x group composition (firm-months) ---\n")
print(sector_by_group_wide)

# Size profile
size_by_group <- dfm[, .(
  n_obs       = .N,
  median_me   = round(median(log_me, na.rm = TRUE), 1),
  mean_me     = round(mean(log_me, na.rm = TRUE), 1),
  p10_me      = round(quantile(log_me, 0.10, na.rm = TRUE), 1),
  p90_me      = round(quantile(log_me, 0.90, na.rm = TRUE), 1)
), by = group]

cat("\n--- Market cap (me) profile by group ---\n")
print(size_by_group)

cat("00_setup.R: data loaded,", nrow(dfm), "firm-months.\n")
print(dfm[, .N, by = group])
