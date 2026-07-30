# =============================================================================
# 08_premium_audience.R
# Is the disclosure premium concentrated where the audience is?
#
# The mechanism behind gamma > 0 is: an audience reads silence, so silence is
# priced, so high emitters stay quiet. Everything testing that audience so far
# lives on BEHAVIOUR (the within-sector IO gradient in gamma, notebook 04b).
# This file puts the same boundary condition on PRICES. If the audience story is
# right, the transparency premium should be larger in high-institutional-
# ownership firms; if the IO gradient in gamma is really a size/capacity
# artefact, the premium should not care about IO at all.
#
# Two objects, one mechanism. That is the point of the file.
#
#   PANEL A  Continuous interaction  ret ~ ... + Disclosed x Audience
#            Audience = within sector-month standing in institutional ownership.
#            Column (2) orthogonalises IO to size FIRST -- IO and log(ME) are
#            strongly correlated and log(ME) is already a control, so the raw
#            interaction is partly a size interaction. (2) is the specification
#            to read; (1) is there to show what the raw measure gives.
#   PANEL B  Same thing in within-sector IO terciles, matching the tercile cut
#            used for gamma in 04b so the two objects are cut the same way.
#   PANEL C  Pricing-side version: the disclosure long-short run separately in
#            the high- and low-IO halves, FF5 alphas and their difference.
#            HALVES, not terciles: the unconditional FF5 alpha is small and only
#            marginally significant, and a three-way split of it has no power.
#
# CAVEAT to carry into any writeup: institutional ownership is chosen, not
# assigned. Nothing here is causal -- this is the same descriptive boundary
# condition as the gamma-by-IO cut, moved onto returns.
#
# Table: returns_premium_audience.tex
#
# Run from Returns/:  Rscript 08_premium_audience.R
# =============================================================================

source("00_setup.R")   # dfm, ctrl, fmt/fmtN/stars/extr; leaves full `df` in scope
library(arrow); library(sandwich); library(lmtest)

NW_LAGS   <- 6
MIN_LEG   <- 10        # min firms per portfolio leg per month (as in 07_)
IO_PATH   <- "../data/processed/io_panel.parquet"

# =============================================================================
# INSTITUTIONAL OWNERSHIP: merge, lag, and build the audience measures
# =============================================================================
cat("\n=== MERGING INSTITUTIONAL OWNERSHIP ===\n")

iop <- as.data.table(read_parquet(IO_PATH))
iop[, gvkey := as.integer(gvkey)]          # zero-padded strings in the parquet
setorder(iop, gvkey, fyear)

# Lag IO by one fiscal year INSIDE the ownership panel, so what gets attached to
# a firm-year is last year's ownership. Guard on the gap: shift() is positional
# and would silently reach across a hole in a firm's fyear sequence.
iop[, `:=`(io_lag = shift(io), fy_gap = fyear - shift(fyear)), by = gvkey]
iop[fy_gap != 1, io_lag := NA_real_]

dd <- merge(dfm, iop[, .(gvkey, fyear, io_L1 = io_lag)],
            by = c("gvkey", "fyear"), all.x = TRUE, sort = FALSE)
setorder(dd, permno, mdate)

cat(sprintf("Firm-months with lagged IO: %s of %s (%.1f%%), %d firms\n",
            fmtN(dd[!is.na(io_L1), .N]), fmtN(nrow(dd)),
            100 * dd[, mean(!is.na(io_L1))], dd[!is.na(io_L1), uniqueN(permno)]))

dd <- dd[!is.na(io_L1)]
dd[, disc_ind := as.integer(group == "disclosed")]

# --- audience measures, all formed WITHIN sector x month ---------------------
# Within-sector is not cosmetic: sectors differ in both ownership and disclosure
# norms, so a raw IO cut would partly reproduce the sector gamma pattern. This
# also mirrors how the tercile cut is formed for gamma in notebook 04b.
zs <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
# residual of IO on size, NA-preserving
resid_on <- function(y, x) {
  out <- rep(NA_real_, length(y))
  ok  <- is.finite(y) & is.finite(x)
  if (sum(ok) < 20) return(out)
  out[ok] <- residuals(lm(y[ok] ~ x[ok]))
  out
}

dd[, io_raw_z  := zs(io_L1),                     by = .(gsector, yyyymm)]
dd[, io_orth   := resid_on(io_L1, log_me_L1),    by = .(gsector, yyyymm)]
dd[, io_orth_z := zs(io_orth),                   by = .(gsector, yyyymm)]

cat(sprintf("corr(IO, log ME) = %.3f  |  corr(IO_orth, log ME) = %.3f\n",
            dd[, cor(io_L1, log_me_L1, use = "complete.obs")],
            dd[, cor(io_orth, log_me_L1, use = "complete.obs")]))

# terciles and halves on the size-orthogonalised measure
terc <- function(x) {
  if (sum(!is.na(x)) < 30) return(rep(NA_character_, length(x)))
  b <- quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (anyDuplicated(b)) return(rep(NA_character_, length(x)))
  as.character(cut(x, breaks = b, labels = c("Low", "Mid", "High"), include.lowest = TRUE))
}
dd[, io_terc := terc(io_orth_z), by = .(gsector, yyyymm)]
dd[, io_half := fifelse(io_orth_z > median(io_orth_z, na.rm = TRUE), "High", "Low"),
   by = .(gsector, yyyymm)]

cat("\nDisclosure rate by within-sector IO tercile:\n")
print(dd[!is.na(io_terc), .(n = .N, io_mean = round(mean(io_L1), 3),
                            disc_rate = round(mean(disc_ind), 3)),
         by = io_terc][order(factor(io_terc, levels = c("Low", "Mid", "High")))])


# =============================================================================
# PANEL A: continuous Disclosed x Audience interaction
# =============================================================================
cat("\n=== PANEL A: CONTINUOUS INTERACTION ===\n")

dA <- dd[!is.na(io_raw_z) & !is.na(io_orth_z)]
dA[, `:=`(disc_x_raw  = disc_ind * io_raw_z,
          disc_x_orth = disc_ind * io_orth_z,
          disc_x_size = disc_ind * log_me_L1_z)]
# sc1_adj has no missings on this panel (MNAR fills every non-discloser), so (3)
# below is estimated on exactly the same rows as (2).
dA[, sc1_adj_z_a := zs(sc1_adj), by = yyyymm]

fit_a1 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind + io_raw_z + disc_x_raw",
                                 "| yyyymm + gsector")),
                cluster = ~permno, data = dA)
fit_a2 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind + io_orth_z + disc_x_orth",
                                 "| yyyymm + gsector")),
                cluster = ~permno, data = dA)
# (3) adds MNAR emissions: does the interaction survive controlling for the level
#     of (recovered) emissions, or is "audience" standing in for brownness?
fit_a3 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind + io_orth_z + disc_x_orth",
                                 "+ sc1_adj_z_a | yyyymm + gsector")),
                cluster = ~permno, data = dA)
# (4) the sharpest guard: let the premium vary with SIZE too. If Disclosed x
#     Audience is really Disclosed x Size in disguise, it dies here.
fit_a4 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind + io_orth_z + disc_x_orth",
                                 "+ disc_x_size | yyyymm + gsector")),
                cluster = ~permno, data = dA)

a1d <- extr(fit_a1, "disc_ind"); a1x <- extr(fit_a1, "disc_x_raw")
a2d <- extr(fit_a2, "disc_ind"); a2x <- extr(fit_a2, "disc_x_orth")
a3d <- extr(fit_a3, "disc_ind"); a3x <- extr(fit_a3, "disc_x_orth")
a4d <- extr(fit_a4, "disc_ind"); a4x <- extr(fit_a4, "disc_x_orth")
a4s <- extr(fit_a4, "disc_x_size")

cat(sprintf("  (1) raw IO       Disc %+.3f  Disc x Aud %+.3f (p %.3f)\n", a1d$coef, a1x$coef, a1x$pval))
cat(sprintf("  (2) size-orth IO Disc %+.3f  Disc x Aud %+.3f (p %.3f)\n", a2d$coef, a2x$coef, a2x$pval))
cat(sprintf("  (3) + emissions  Disc %+.3f  Disc x Aud %+.3f (p %.3f)\n", a3d$coef, a3x$coef, a3x$pval))
cat(sprintf("  (4) + Disc x Size Disc %+.3f  Disc x Aud %+.3f (p %.3f)  Disc x Size %+.3f (p %.3f)\n",
            a4d$coef, a4x$coef, a4x$pval, a4s$coef, a4s$pval))


# =============================================================================
# PANEL B: premium by within-sector IO tercile (cut to match gamma in 04b)
# =============================================================================
cat("\n=== PANEL B: PREMIUM BY IO TERCILE ===\n")

dB <- dd[!is.na(io_terc)]
for (t in c("Low", "Mid", "High")) {
  dB[, paste0("disc_", t) := disc_ind * (io_terc == t)]
  dB[, paste0("iot_",  t) := as.integer(io_terc == t)]
}
# tercile FE via iot_Mid / iot_High (Low is the omitted level)
fit_b <- feols(as.formula(paste("ret ~", ctrl,
                 "+ iot_Mid + iot_High + disc_Low + disc_Mid + disc_High",
                 "| yyyymm + gsector")),
               cluster = ~permno, data = dB)

bt <- lapply(c("Low", "Mid", "High"), function(t) extr(fit_b, paste0("disc_", t)))
names(bt) <- c("Low", "Mid", "High")

# High - Low by delta method off the clustered vcov
V     <- vcov(fit_b)
d_hl  <- coef(fit_b)["disc_High"] - coef(fit_b)["disc_Low"]
se_hl <- sqrt(max(V["disc_High", "disc_High"] + V["disc_Low", "disc_Low"]
                  - 2 * V["disc_High", "disc_Low"], 0))
t_hl  <- d_hl / se_hl
p_hl  <- 2 * pnorm(-abs(t_hl))

for (t in c("Low", "Mid", "High"))
  cat(sprintf("  Premium, %-4s IO: %+.3f (SE %.3f, p %.3f)\n",
              t, bt[[t]]$coef, bt[[t]]$se, bt[[t]]$pval))
cat(sprintf("  High - Low      : %+.3f (SE %.3f, p %.3f)\n", d_hl, se_hl, p_hl))


# =============================================================================
# PANEL C: the disclosure long-short inside each IO half (FF5 alphas)
# =============================================================================
cat("\n=== PANEL C: DISCLOSURE LONG-SHORT BY IO HALF ===\n")

read_ff <- function(path, cols) {
  raw <- readLines(path)
  dat <- raw[grepl("^[0-9]{6},", raw)]
  m <- fread(text = paste(dat, collapse = "\n"), header = FALSE)
  setnames(m, c("ym", cols))
  for (c in cols) m[get(c) <= -99, (c) := NA]
  m[, ym := as.integer(ym)][]
}
ff5 <- read_ff("F-F_Research_Data_5_Factors_2x3.csv",
               c("mktrf", "smb", "hml", "rmw", "cma", "rf"))
mom <- read_ff("F-F_Momentum_Factor.csv", c("mom"))
fac <- merge(ff5, mom, by = "ym")
FF5 <- c("mktrf", "smb", "hml", "rmw", "cma")
FF6 <- c(FF5, "mom")

nw_mean <- function(x, lags = NW_LAGS) {
  x <- x[is.finite(x)]; m <- lm(x ~ 1)
  ct <- coeftest(m, vcov = NeweyWest(m, lag = lags, prewhite = FALSE))
  c(est = unname(ct[1, 1]), t = unname(ct[1, 3]), p = unname(ct[1, 4]), n = length(x))
}
nw_alpha <- function(y, F, lags = NW_LAGS) {
  ok <- is.finite(y) & complete.cases(F); y <- y[ok]; F <- as.data.frame(F[ok])
  m  <- lm(y ~ ., data = cbind(y = y, F))
  ct <- coeftest(m, vcov = NeweyWest(m, lag = lags, prewhite = FALSE))
  c(alpha = unname(ct[1, 1]), t = unname(ct[1, 3]), p = unname(ct[1, 4]))
}
alpha_set <- function(spread, ym) {
  d <- merge(data.table(ym = ym, y = spread), fac, by = "ym")[is.finite(y)]
  list(mean = nw_mean(d$y), ff5 = nw_alpha(d$y, d[, ..FF5]), ff6 = nw_alpha(d$y, d[, ..FF6]))
}
report <- function(a, lab) cat(sprintf(
  "  %-26s mean %+.3f (t=%.2f)%s | FF5 %+.3f (t=%.2f)%s | FF5+MOM %+.3f (t=%.2f)%s\n",
  lab, a$mean["est"], a$mean["t"], stars(a$mean["p"]),
  a$ff5["alpha"], a$ff5["t"], stars(a$ff5["p"]),
  a$ff6["alpha"], a$ff6["t"], stars(a$ff6["p"])))

# bring market equity back for value weights (dropped from dfm in 00_setup)
px <- unique(df[, .(permno, mdate, me)])
dC <- merge(dd[!is.na(io_half)], px, by = c("permno", "mdate"), all.x = TRUE)
setorder(dC, permno, mdate)
dC[, me_lag := shift(me), by = permno]
dC[, ym := as.integer(format(as.Date(yyyymm), "%Y%m"))]

vwm <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok])
}
ls_series <- function(d) {
  p <- d[is.finite(ret), .(
    vw_d  = vwm(ret[disc_ind == 1], me_lag[disc_ind == 1]),
    ew_d  = mean(ret[disc_ind == 1], na.rm = TRUE),
    vw_nd = vwm(ret[disc_ind == 0], me_lag[disc_ind == 0]),
    ew_nd = mean(ret[disc_ind == 0], na.rm = TRUE),
    n_d   = sum(disc_ind == 1), n_nd = sum(disc_ind == 0)), by = ym][order(ym)]
  p <- p[n_d >= MIN_LEG & n_nd >= MIN_LEG]
  p[, `:=`(vw = vw_d - vw_nd, ew = ew_d - ew_nd)][]
}

s_all  <- ls_series(dC)
s_high <- ls_series(dC[io_half == "High"])
s_low  <- ls_series(dC[io_half == "Low"])

c_all  <- alpha_set(s_all$vw,  s_all$ym);  report(c_all,  "All firms (VW)")
c_high <- alpha_set(s_high$vw, s_high$ym); report(c_high, "High-IO half (VW)")
c_low  <- alpha_set(s_low$vw,  s_low$ym);  report(c_low,  "Low-IO half (VW)")

# difference-of-spreads, only on months where both halves survive the leg filter
s_diff <- merge(s_high[, .(ym, hi = vw)], s_low[, .(ym, lo = vw)], by = "ym")
c_diff <- alpha_set(s_diff$hi - s_diff$lo, s_diff$ym)
report(c_diff, "High - Low (VW)")
cat(sprintf("  months: all %d | high %d | low %d | both %d\n",
            nrow(s_all), nrow(s_high), nrow(s_low), nrow(s_diff)))


# =============================================================================
# TABLE
# =============================================================================
cell <- function(e) paste0(fmt(e$coef), stars(e$pval))
se_c <- function(e) paste0("(", fmt(e$se), ")")
# Panel C fills only 4 of the table's 5 columns, so each row carries a trailing
# empty cell before the row break.
r3   <- function(a, lab) c(
  sprintf("%s & %+.3f%s & %+.3f%s & %+.3f%s & \\\\", lab,
          a$mean["est"], stars(a$mean["p"]), a$ff5["alpha"], stars(a$ff5["p"]),
          a$ff6["alpha"], stars(a$ff6["p"])),
  sprintf(" & (%.2f) & (%.2f) & (%.2f) & \\\\", a$mean["t"], a$ff5["t"], a$ff6["t"]))

tex <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Is the Disclosure Premium Concentrated Where the Audience Is?}",
  "Audience is proxied by lagged institutional ownership, measured relative to",
  "other firms in the same sector and month. Panel~A interacts the disclosure",
  "indicator with the continuous audience measure: column~(1) uses raw ownership,",
  "columns~(2)--(4) use ownership orthogonalised to firm size within sector-month,",
  "column~(3) adds MNAR-recovered emissions, and column~(4) additionally lets the",
  "premium vary with size. Panel~B repeats the test in within-sector ownership",
  "terciles, the same cut used for $\\hat\\gamma$. Panel~C reports the",
  "value-weighted disclosure long-short estimated separately in the high- and",
  "low-ownership halves. Panels~A--B cluster standard errors at the firm level;",
  "Panel~C reports Newey-West $t$-statistics (6 lags).}",
  "\\label{tab:returns_premium_audience}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  "\\multicolumn{5}{l}{\\emph{Panel A: continuous interaction}} \\\\[2pt]",
  " & (1) & (2) & (3) & (4) \\\\",
  " & Raw IO & Size-orth.\\ IO & $+$ Emissions & $+$ Disc.\\ $\\times$ Size \\\\",
  "\\midrule",
  paste0("Disclosed & ", cell(a1d), " & ", cell(a2d), " & ", cell(a3d), " & ", cell(a4d), " \\\\"),
  paste0(" & ", se_c(a1d), " & ", se_c(a2d), " & ", se_c(a3d), " & ", se_c(a4d), " \\\\"),
  "",
  paste0("Disclosed $\\times$ Audience & ", cell(a1x), " & ", cell(a2x), " & ",
         cell(a3x), " & ", cell(a4x), " \\\\"),
  paste0(" & ", se_c(a1x), " & ", se_c(a2x), " & ", se_c(a3x), " & ", se_c(a4x), " \\\\"),
  "",
  paste0("Disclosed $\\times$ Size & -- & -- & -- & ", cell(a4s), " \\\\"),
  paste0(" & & & & ", se_c(a4s), " \\\\"),
  "\\midrule",
  "Controls & Yes & Yes & Yes & Yes \\\\",
  "Time FE & Yes & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes & Yes \\\\",
  paste0("Observations & ", fmtN(a1d$n), " & ", fmtN(a2d$n), " & ",
         fmtN(a3d$n), " & ", fmtN(a4d$n), " \\\\"),
  "\\midrule",
  "\\multicolumn{5}{l}{\\emph{Panel B: premium by within-sector ownership tercile}} \\\\[2pt]",
  " & Low & Mid & High & High$-$Low \\\\", "\\midrule",
  paste0("Disclosure premium & ", cell(bt$Low), " & ", cell(bt$Mid), " & ", cell(bt$High),
         " & ", fmt(d_hl), stars(p_hl), " \\\\"),
  paste0(" & ", se_c(bt$Low), " & ", se_c(bt$Mid), " & ", se_c(bt$High),
         " & (", fmt(se_hl), ") \\\\"),
  paste0("Observations & \\multicolumn{4}{c}{", fmtN(bt$Low$n), "} \\\\"),
  "\\midrule",
  "\\multicolumn{5}{l}{\\emph{Panel C: disclosure long-short, value-weighted}} \\\\[2pt]",
  " & Mean & FF5 $\\alpha$ & FF5$+$MOM $\\alpha$ & \\\\", "\\midrule",
  r3(c_all,  "All firms"),
  r3(c_high, "High-IO half"),
  r3(c_low,  "Low-IO half"),
  r3(c_diff, "High $-$ Low"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Institutional ownership is lagged one fiscal year, so it is",
  "predetermined relative to the returns it is interacted with. Ownership is not",
  "randomly assigned: this is the same descriptive boundary condition that the",
  "$\\hat\\gamma$-by-ownership cut applies to disclosure behaviour, applied here to",
  "prices, and it is not a causal estimate. A positive Disclosed $\\times$ Audience",
  "coefficient means transparency is worth more where more institutional capital is",
  "watching -- the pricing counterpart of steeper selection in the same firms.",
  "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")

writeLines(tex, file.path(outdir_tab, "returns_premium_audience.tex"))
cat("\nreturns_premium_audience.tex written\n")
cat("\n=== 08_premium_audience.R COMPLETE ===\n")
