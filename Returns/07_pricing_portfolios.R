# =============================================================================
# 15_pricing_portfolios.R
# Portfolio-sort and Fama-MacBeth pricing tests -- the asset-pricing complements
# to the panel regressions of Sections 5.1-5.2.
#
#   T1  Univariate quintile sorts on emissions: value- and equal-weighted
#       high-minus-low spreads, with raw mean, FF5 alpha, and FF5+MOM alpha.
#       Expectation from the panel: a null (the level is not priced).
#   T2  Disclosure long-short (disclosers minus non-disclosers) -- the 0.242
#       transparency premium in tradable form. The test is whether it survives
#       value-weighting and factor adjustment.
#   T3  Fama-MacBeth cross-sectional regressions (Newey-West), the FM counterpart
#       to the panel FE specification.
#
# Factor data: Ken French monthly FF5 (Mkt-RF, SMB, HML, RMW, CMA, RF) and the
# momentum factor, in this folder. Long-short SPREADS are zero-cost, so they are
# regressed on factors directly (no risk-free adjustment needed); alphas are the
# intercepts. Newey-West t-stats, 6 lags.
#
# Run from Returns/:  Rscript 15_pricing_portfolios.R
# =============================================================================

source("00_setup.R")   # dfm, ctrl, fmt/fmtN/stars; leaves full `df` in scope
library(sandwich); library(lmtest)

NW_LAGS <- 6

# --- load Ken French factors -------------------------------------------------
read_ff <- function(path, cols) {
  raw <- readLines(path)
  dat <- raw[grepl("^[0-9]{6},", raw)]          # monthly rows only (6-digit YYYYMM)
  m <- fread(text = paste(dat, collapse = "\n"), header = FALSE)
  setnames(m, c("ym", cols))
  for (c in cols) m[get(c) <= -99, (c) := NA]   # French missing codes
  m[, ym := as.integer(ym)][]
}
ff5 <- read_ff("F-F_Research_Data_5_Factors_2x3.csv",
               c("mktrf", "smb", "hml", "rmw", "cma", "rf"))
mom <- read_ff("F-F_Momentum_Factor.csv", c("mom"))
fac <- merge(ff5, mom, by = "ym")
cat(sprintf("Factors loaded: %d monthly obs, %d--%d\n",
            nrow(fac), min(fac$ym), max(fac$ym)))

# --- estimation panel: bring on market equity and a YYYYMM key ---------------
px <- unique(df[, .(permno, mdate, me)])
dd <- merge(dfm, px, by = c("permno", "mdate"), all.x = TRUE)
setorder(dd, permno, mdate)
dd[, me_lag := shift(me), by = permno]
dd[, disclosed := as.integer(group == "disclosed")]
dd[, ym := as.integer(format(as.Date(yyyymm), "%Y%m"))]
cat(sprintf("Panel: %s firm-months, %d months, %d firms\n",
            fmtN(nrow(dd)), dd[, uniqueN(ym)], dd[, uniqueN(permno)]))

# safe value-weighted mean (drops NA / non-positive weights)
vwm <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok])
}

# --- HAC helpers -------------------------------------------------------------
nw_mean <- function(x, lags = NW_LAGS) {
  x <- x[is.finite(x)]; m <- lm(x ~ 1)
  ct <- coeftest(m, vcov = NeweyWest(m, lag = lags, prewhite = FALSE))
  c(est = unname(ct[1, 1]), t = unname(ct[1, 3]), p = unname(ct[1, 4]), n = length(x))
}
# alpha of a zero-cost spread on a set of factor columns (data.table F)
nw_alpha <- function(y, F, lags = NW_LAGS) {
  ok <- is.finite(y) & complete.cases(F); y <- y[ok]; F <- as.data.frame(F[ok])
  m <- lm(y ~ ., data = cbind(y = y, F))
  ct <- coeftest(m, vcov = NeweyWest(m, lag = lags, prewhite = FALSE))
  c(alpha = unname(ct[1, 1]), t = unname(ct[1, 3]), p = unname(ct[1, 4]))
}
star <- function(p) stars(p)

FF5  <- c("mktrf", "smb", "hml", "rmw", "cma")
FF6  <- c(FF5, "mom")

# spread series -> (mean, FF5 alpha, FF6 alpha), merging factors by ym
alpha_set <- function(spread, ym) {
  d <- merge(data.table(ym = ym, y = spread), fac, by = "ym")[is.finite(y)]
  list(mean = nw_mean(d$y),
       ff5  = nw_alpha(d$y, d[, ..FF5]),
       ff6  = nw_alpha(d$y, d[, ..FF6]))
}
report <- function(a, lab) cat(sprintf(
  "  %-30s mean %+.3f (t=%.2f)%s | FF5 %+.3f (t=%.2f)%s | FF5+MOM %+.3f (t=%.2f)%s\n",
  lab, a$mean["est"], a$mean["t"], star(a$mean["p"]),
  a$ff5["alpha"], a$ff5["t"], star(a$ff5["p"]),
  a$ff6["alpha"], a$ff6["t"], star(a$ff6["p"])))

# =============================================================================
# T1: quintile sorts on emissions
# =============================================================================
quintile_spread <- function(data, sortvar, label) {
  d <- data[is.finite(get(sortvar)) & is.finite(ret)]
  d[, q := as.integer(cut(get(sortvar),
        breaks = quantile(get(sortvar), 0:5/5, na.rm = TRUE),
        include.lowest = TRUE)), by = ym]
  d <- d[!is.na(q)]
  port <- d[, .(vw = vwm(ret, me_lag), ew = mean(ret, na.rm = TRUE)), by = .(ym, q)]
  w_vw <- dcast(port, ym ~ q, value.var = "vw"); w_ew <- dcast(port, ym ~ q, value.var = "ew")
  setnames(w_vw, as.character(1:5), paste0("q", 1:5))
  setnames(w_ew, as.character(1:5), paste0("q", 1:5))
  cat(sprintf("\n[T1] %s  (H = Q5 top emitter, L = Q1)\n", label))
  a_vw <- alpha_set(w_vw$q5 - w_vw$q1, w_vw$ym); report(a_vw, "VW  H-L")
  a_ew <- alpha_set(w_ew$q5 - w_ew$q1, w_ew$ym); report(a_ew, "EW  H-L")
  list(vw = a_vw, ew = a_ew)
}
t1_full <- quintile_spread(dd, "sc1_adj", "Recovered emissions, full universe")
t1_disc <- quintile_spread(dd[group == "disclosed"], "sc1_disclosed",
                           "Disclosed emissions, disclosers only")

# =============================================================================
# T2: the disclosure long-short (transparency premium as a spread)
# =============================================================================
pd <- dd[is.finite(ret), .(
  vw_d  = vwm(ret[disclosed == 1], me_lag[disclosed == 1]),
  ew_d  = mean(ret[disclosed == 1], na.rm = TRUE),
  vw_nd = vwm(ret[disclosed == 0], me_lag[disclosed == 0]),
  ew_nd = mean(ret[disclosed == 0], na.rm = TRUE),
  n_d = sum(disclosed == 1), n_nd = sum(disclosed == 0)), by = ym][order(ym)]
pd <- pd[n_d >= 10 & n_nd >= 10]
cat("\n[T2] Disclosure long-short: disclosers minus non-disclosers\n")
t2_vw <- alpha_set(pd$vw_d - pd$vw_nd, pd$ym); report(t2_vw, "VW  L-S")
t2_ew <- alpha_set(pd$ew_d - pd$ew_nd, pd$ym); report(t2_ew, "EW  L-S")

# =============================================================================
# T3: Fama-MacBeth cross-sectional regressions
# =============================================================================
ctrl_vars <- strsplit(gsub("[+]", " ", ctrl), " +")[[1]]; ctrl_vars <- ctrl_vars[nzchar(ctrl_vars)]
fama_macbeth <- function(rhs, label) {
  vars <- c("ret", rhs); d <- dd[complete.cases(dd[, ..vars])]
  f <- as.formula(paste("ret ~", paste(rhs, collapse = " + ")))
  coefs <- d[, as.list(coef(lm(f, data = .SD))), by = ym]
  cat(sprintf("\n[T3] Fama-MacBeth: %s   (%d months)\n", label, nrow(coefs)))
  out <- rbindlist(lapply(setdiff(names(coefs), "ym"), function(v) {
    s <- nw_mean(coefs[[v]]); data.table(term = v, coef = s["est"], t = s["t"], p = s["p"]) }))
  for (tm in intersect(c("disclosed", "sc1_adj_z"), out$term)) {
    r <- out[term == tm]; cat(sprintf("   %-12s %+.4f (t=%.2f)%s\n", tm, r$coef, r$t, star(r$p))) }
  out
}
fm1 <- fama_macbeth("disclosed", "disclosed only")
fm2 <- fama_macbeth(c("disclosed", "sc1_adj_z", ctrl_vars), "disclosed + emissions + controls")

# =============================================================================
# TABLE
# =============================================================================
r3 <- function(a, lab) c(
  sprintf("%s & %+.3f%s & %+.3f%s & %+.3f%s \\\\", lab,
          a$mean["est"], star(a$mean["p"]), a$ff5["alpha"], star(a$ff5["p"]),
          a$ff6["alpha"], star(a$ff6["p"])),
  sprintf(" & (%.2f) & (%.2f) & (%.2f) \\\\", a$mean["t"], a$ff5["t"], a$ff6["t"]))
fm_cell <- function(out, tm) { r <- out[term == tm]
  if (nrow(r) == 0) "--" else paste0(fmt(round(r$coef, 3)), star(r$p)) }

tex <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Portfolio-sort and Fama-MacBeth pricing tests.}",
  "Panel A reports monthly long-short spreads: emissions high-minus-low",
  "(top-quintile minus bottom-quintile emitter) and the disclosure long-short",
  "(disclosers minus non-disclosers), value- and equal-weighted, with the raw",
  "mean and alphas from the Fama-French five-factor model and that model augmented",
  "with momentum. Panel B reports Fama-MacBeth coefficients on the full sample.",
  "Returns are monthly percentage points; Newey-West $t$-statistics (6 lags) in",
  "parentheses.}",
  "\\label{tab:returns_pricing_portfolios}",
  "\\begin{tabular}{lccc}", "\\toprule",
  "\\multicolumn{4}{l}{\\emph{Panel A: long-short portfolios}} \\\\[2pt]",
  " & Mean & FF5 $\\alpha$ & FF5$+$MOM $\\alpha$ \\\\", "\\midrule",
  r3(t1_full$vw, "Emissions Q5$-$Q1, full (VW)"),
  r3(t1_full$ew, "Emissions Q5$-$Q1, full (EW)"),
  r3(t1_disc$vw, "Emissions Q5$-$Q1, disclosers (VW)"),
  r3(t1_disc$ew, "Emissions Q5$-$Q1, disclosers (EW)"),
  "\\midrule",
  r3(t2_vw, "Disclosed $-$ Non-disclosed (VW)"),
  r3(t2_ew, "Disclosed $-$ Non-disclosed (EW)"),
  "\\midrule",
  "\\multicolumn{4}{l}{\\emph{Panel B: Fama-MacBeth, full sample}} \\\\[2pt]",
  " & Disclosed & Emissions ($z$) & \\\\", "\\midrule",
  paste0("Disclosure only & ", fm_cell(fm1, "disclosed"), " & -- & \\\\"),
  paste0("$+$ emissions, controls & ", fm_cell(fm2, "disclosed"), " & ",
         fm_cell(fm2, "sc1_adj_z"), " & \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Alphas are intercepts from regressing the monthly spread on the",
  "stated factors; zero-cost spreads need no risk-free adjustment. Factors are Ken",
  "French's monthly series. Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex, file.path(outdir_tab, "returns_pricing_portfolios.tex"))
cat("\nreturns_pricing_portfolios.tex written\n")
cat("\n===== 15_pricing_portfolios.R complete =====\n")
