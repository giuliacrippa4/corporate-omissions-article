# =============================================================================
# 05_mechanisms.R
# Mechanism tests for the never-disclosing carbon premium:
#
#   BLOCK 1 — BK replication (level, raw emissions) — appendix benchmark
#   BLOCK 2 — Disclosure premium (Verrecchia/Milgrom test)
#   BLOCK 3 — ESG-salience period split + pooled interaction + disclosing placebo
#   BLOCK 4 — Disclosure surprise at first disclosure + event-year window
#
# The transition DiD formerly lived here; it now lives in 06_did_analysis.R
# together with its event-time, placebo and cohort diagnostics. Nothing in this
# file defines a high-emitter cutoff.
#
# Tables:
#   returns_bk.tex               — BK replication (null, consistent with Aswani)
#   returns_disc_premium.tex     — disclosure premium
#   returns_never_periods.tex    — period split (main-text mechanism)
#   returns_surprise_main.tex    — disclosure surprise
#   returns_surprise_window.tex  — surprise at t-1, t=0, t+1
#
# NOTE: returns_gamma_bias.tex was DROPPED in returns_complete.R
# (the MNAR-MAR gap conflates coverage and severity; use the Python gamma MC
# figure instead). Not reproduced here.
# =============================================================================

source("00_setup.R")

# Sample definition for the never-disclosing analyses (Block 3), matching the
# group-based filter in 02_pricing_never_disclosing.R and 03_event_study.R:
#   "none"           — neither disclosed nor vendor-estimated in that firm-month
#   "none_estimated" — not disclosed; vendor estimate allowed
NEVER_GROUPS <- "none_estimated"


# =============================================================================
# BLOCK 1: BK REPLICATION (appendix)
# =============================================================================
cat("\n=== BLOCK 1: BK REPLICATION ===\n")

bk <- dfm[group %in% c("disclosed", "estimated")]
bk[!is.na(sc1_combined),         sc1_combined_z_sub := as.numeric(scale(sc1_combined)),         by = yyyymm]
bk[!is.na(sc1_adj),              sc1_adj_z_sub      := as.numeric(scale(sc1_adj)),              by = yyyymm]
bk[!is.na(sc1_disclosed_filled), sc1_filled_z_sub   := as.numeric(scale(sc1_disclosed_filled)), by = yyyymm]
dfm[!is.na(sc1_adj),              sc1_adj_z_full    := as.numeric(scale(sc1_adj)),              by = yyyymm]
dfm[!is.na(sc1_disclosed_filled), sc1_filled_z_full := as.numeric(scale(sc1_disclosed_filled)), by = yyyymm]

fit_bk_obs    <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_combined_z_sub | yyyymm + gsector")),
                       cluster = ~permno, data = bk[!is.na(sc1_combined_z_sub)])
fit_bk_mnar   <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_sub      | yyyymm + gsector")),
                       cluster = ~permno, data = bk[!is.na(sc1_adj_z_sub)])
fit_bk_mar    <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_filled_z_sub   | yyyymm + gsector")),
                       cluster = ~permno, data = bk[!is.na(sc1_filled_z_sub)])
fit_full_mnar <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_full    | yyyymm + gsector")),
                       cluster = ~permno, data = dfm[!is.na(sc1_adj_z_full)])
fit_full_mar  <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_filled_z_full | yyyymm + gsector")),
                       cluster = ~permno, data = dfm[!is.na(sc1_filled_z_full)])

bk_obs  <- extr(fit_bk_obs,    "sc1_combined_z_sub")
bk_mnar <- extr(fit_bk_mnar,   "sc1_adj_z_sub")
bk_mar  <- extr(fit_bk_mar,    "sc1_filled_z_sub")
fu_mnar <- extr(fit_full_mnar, "sc1_adj_z_full")
fu_mar  <- extr(fit_full_mar,  "sc1_filled_z_full")

tex_bk <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Premium: Bolton \\& Kacperczyk Replication and MNAR Correction.}",
  "Panel regressions of monthly returns on standardized Scope~1 emissions.",
  "Column~(1) replicates BK on the disclosed-plus-estimated subsample.",
  "Columns~(2)--(3) replace observed emissions with MNAR and MAR imputations.",
  "Columns~(4)--(5) extend to the full Compustat universe.",
  "Standard errors clustered at the firm level.}",
  "\\label{tab:returns_bk}",
  "\\begin{tabular}{lccccc}", "\\toprule",
  " & (1) & (2) & (3) & (4) & (5) \\\\",
  " & BK Sample & BK Sample & BK Sample & Full Universe & Full Universe \\\\",
  " & Observed  & MNAR      & MAR       & MNAR          & MAR           \\\\",
  "\\midrule",
  paste0("Scope~1 Emissions & ",
         fmt(bk_obs$coef),  stars(bk_obs$pval),  " & ",
         fmt(bk_mnar$coef), stars(bk_mnar$pval), " & ",
         fmt(bk_mar$coef),  stars(bk_mar$pval),  " & ",
         fmt(fu_mnar$coef), stars(fu_mnar$pval), " & ",
         fmt(fu_mar$coef),  stars(fu_mar$pval),  " \\\\"),
  paste0(" & (", fmt(bk_obs$se), ") & (", fmt(bk_mnar$se), ") & (",
         fmt(bk_mar$se), ") & (", fmt(fu_mnar$se), ") & (", fmt(fu_mar$se), ") \\\\"),
  "\\midrule",
  "Controls    & Yes & Yes & Yes & Yes & Yes \\\\",
  "Time FE     & Yes & Yes & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes & Yes & Yes \\\\",
  paste0("Observations & ", fmtN(bk_obs$n), " & ", fmtN(bk_mnar$n), " & ",
         fmtN(bk_mar$n), " & ", fmtN(fu_mnar$n), " & ", fmtN(fu_mar$n), " \\\\"),
  paste0("$R^2$ & ", fmt(bk_obs$r2), " & ", fmt(bk_mnar$r2), " & ",
         fmt(bk_mar$r2), " & ", fmt(fu_mnar$r2), " & ", fmt(fu_mar$r2), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} No specification produces a significant carbon premium on",
  "raw emissions, consistent with \\citet{aswani2024}. MNAR correction does not",
  "change this conclusion. Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_bk, file.path(outdir_tab, "returns_bk.tex"))
cat("returns_bk.tex written\n")


# =============================================================================
# BLOCK 2: DISCLOSURE PREMIUM (main text)
# Tests Verrecchia/Milgrom: does disclosure carry a return premium beyond
# what MNAR-corrected emissions explain?
# =============================================================================
cat("\n=== BLOCK 2: DISCLOSURE PREMIUM ===\n")

dfm[, disc_ind := as.integer(group == "disclosed")]

fit_dp1 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind | yyyymm + gsector")),
                 cluster = ~permno, data = dfm)
fit_dp2 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind + sc1_adj_z_full | yyyymm + gsector")),
                 cluster = ~permno, data = dfm[!is.na(sc1_adj_z_full)])
fit_dp3 <- feols(as.formula(paste("ret ~", ctrl, "+ disc_ind + sc1_filled_z_full | yyyymm + gsector")),
                 cluster = ~permno, data = dfm[!is.na(sc1_filled_z_full)])
fit_dp4 <- feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_full | yyyymm + gsector")),
                 cluster = ~permno, data = dfm[!is.na(sc1_adj_z_full)])

dp1  <- extr(fit_dp1, "disc_ind")
dp2d <- extr(fit_dp2, "disc_ind");           dp2e <- extr(fit_dp2, "sc1_adj_z_full")
dp3d <- extr(fit_dp3, "disc_ind");           dp3e <- extr(fit_dp3, "sc1_filled_z_full")
dp4e <- extr(fit_dp4, "sc1_adj_z_full")

tex_dp <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Disclosure Premium: Does MNAR Correction Absorb the Return to Disclosure?}",
  "Panel regressions on the full Compustat sample.",
  "Column~(1) estimates the return premium to voluntary disclosure.",
  "Column~(2) adds MNAR-corrected emissions.",
  "Column~(3) repeats with MAR-imputed emissions.",
  "Column~(4) reports MNAR emissions without the disclosure dummy.",
  "Standard errors clustered at the firm level.}",
  "\\label{tab:returns_disc_premium}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  " & (1) & (2) & (3) & (4) \\\\",
  " & Disc Only & +MNAR & +MAR & MNAR Only \\\\", "\\midrule",
  paste0("Disclosed & ",
         fmt(dp1$coef),  stars(dp1$pval),  " & ",
         fmt(dp2d$coef), stars(dp2d$pval), " & ",
         fmt(dp3d$coef), stars(dp3d$pval), " & -- \\\\"),
  paste0(" & (", fmt(dp1$se), ") & (", fmt(dp2d$se), ") & (", fmt(dp3d$se), ") & \\\\"),
  "",
  paste0("Scope~1 Emissions & -- & ",
         fmt(dp2e$coef), stars(dp2e$pval), " & ",
         fmt(dp3e$coef), stars(dp3e$pval), " & ",
         fmt(dp4e$coef), stars(dp4e$pval), " \\\\"),
  paste0(" & & (", fmt(dp2e$se), ") & (", fmt(dp3e$se), ") & (", fmt(dp4e$se), ") \\\\"),
  "\\midrule",
  "Controls    & Yes & Yes & Yes & Yes \\\\",
  "Time FE     & Yes & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes & Yes \\\\",
  paste0("Observations & ", fmtN(dp1$n), " & ", fmtN(dp2d$n), " & ",
         fmtN(dp3d$n), " & ", fmtN(dp4e$n), " \\\\"),
  paste0("$R^2$ & ", fmt(dp1$r2), " & ", fmt(dp2d$r2), " & ",
         fmt(dp3d$r2), " & ", fmt(dp4e$r2), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Stability of the Disclosed coefficient from~(1) to~(2) indicates",
  "the market prices the act of transparency, not the emissions level,",
  "consistent with \\citet{verrecchia1983} and \\citet{milgrom1981}.",
  "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_dp, file.path(outdir_tab, "returns_disc_premium.tex"))
cat("returns_disc_premium.tex written\n")


# =============================================================================
# BLOCK 3: ESG-SALIENCE PERIOD SPLIT (main text)
# Never-disclosing premium by period + pooled interaction (Wald test) +
# disclosing-firm placebo.
# =============================================================================
cat("\n=== BLOCK 3: ESG-SALIENCE PERIOD SPLIT ===\n")

keep <- switch(NEVER_GROUPS,
               "none"           = "none",
               "none_estimated" = c("none", "estimated"),
               stop("NEVER_GROUPS must be 'none' or 'none_estimated'"))
never <- dfm[group %in% keep]
never[!is.na(sc1_adj), sc1_adj_z_nev := as.numeric(scale(sc1_adj)), by = yyyymm]

assign_period <- function(dt) {
  dt[, period := fcase(
    year %in% 2013:2018, "baseline",
    year %in% 2019:2021, "esg_peak",
    year %in% 2022:2023, "regulatory",
    default = NA_character_)]
  dt[, period := factor(period, levels = c("baseline", "esg_peak", "regulatory"))]
  dt
}
assign_period(never)

periods       <- c("baseline", "esg_peak", "regulatory")
period_labels <- c("2013--2018 (Baseline)", "2019--2021 (ESG Peak)", "2022--2023 (Regulatory)")

period_res <- rbindlist(lapply(seq_along(periods), function(i) {
  p <- periods[i]
  d <- never[period == p & !is.na(sc1_adj_z_nev)]
  d[, sc1_adj_z_per := as.numeric(scale(sc1_adj)), by = yyyymm]
  fit <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl, "+ sc1_adj_z_per | yyyymm + gsector")),
          cluster = ~permno, data = d),
    error = function(e) NULL)
  if (is.null(fit) || !"sc1_adj_z_per" %in% rownames(coeftable(fit)))
    return(data.table(period=p, label=period_labels[i],
                      coef=NA_real_, se=NA_real_, pval=NA_real_, n=NA_integer_, r2=NA_real_))
  co <- coeftable(fit)["sc1_adj_z_per", ]
  data.table(period=p, label=period_labels[i],
             coef=round(co["Estimate"],3), se=round(co["Std. Error"],3),
             pval=co["Pr(>|t|)"], n=fit$nobs, r2=round(r2(fit,"r2"),3))
}))
print(period_res[, .(label, coef, se, pval, n)])

never_pool <- never[!is.na(sc1_adj_z_nev) & !is.na(period)]
never_pool[, sc1_x_esg := sc1_adj_z_nev * (period == "esg_peak")]
never_pool[, sc1_x_reg := sc1_adj_z_nev * (period == "regulatory")]
fit_pool <- feols(
  as.formula(paste("ret ~", ctrl,
    "+ sc1_adj_z_nev + sc1_x_esg + sc1_x_reg | yyyymm + gsector")),
  cluster = ~permno, data = never_pool)
wald_res  <- wald(fit_pool, c("sc1_x_esg", "sc1_x_reg"))
pool_base <- extr(fit_pool, "sc1_adj_z_nev")
pool_esg  <- coeftable(fit_pool)["sc1_x_esg", ]
pool_reg  <- coeftable(fit_pool)["sc1_x_reg", ]
cat(sprintf("Wald p-value (periods jointly zero): %.3f\n", wald_res$p))

# Disclosing placebo
disclosing <- assign_period(dfm[group == "disclosed"])
placebo_res <- rbindlist(lapply(seq_along(periods), function(i) {
  p <- periods[i]
  d <- disclosing[period == p & !is.na(sc1_disclosed)]
  d[, sc1_disc_z_per := as.numeric(scale(sc1_disclosed)), by = yyyymm]
  fit <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl, "+ sc1_disc_z_per | yyyymm + gsector")),
          cluster = ~permno, data = d),
    error = function(e) NULL)
  if (is.null(fit) || !"sc1_disc_z_per" %in% rownames(coeftable(fit)))
    return(data.table(period=p, label=period_labels[i],
                      coef=NA_real_, se=NA_real_, pval=NA_real_, n=NA_integer_))
  co <- coeftable(fit)["sc1_disc_z_per", ]
  data.table(period=p, label=period_labels[i],
             coef=round(co["Estimate"],3), se=round(co["Std. Error"],3),
             pval=co["Pr(>|t|)"], n=fit$nobs)
}))

p1 <- period_res[period=="baseline"];   p2 <- period_res[period=="esg_peak"];   p3 <- period_res[period=="regulatory"]
pl1 <- placebo_res[period=="baseline"]; pl2 <- placebo_res[period=="esg_peak"]; pl3 <- placebo_res[period=="regulatory"]

tex_per <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Carbon Premium in Never-Disclosing Firms: ESG-Salience Periods.}",
  "Columns~(1)--(3): never-disclosing firms by period.",
  "Columns~(4)--(6): voluntarily disclosing firms (placebo).",
  "Emissions re-standardized within period and month.",
  "Controls, time and industry FEs. Standard errors clustered at firm level.}",
  "\\label{tab:returns_never_periods}",
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{3}{c}{Never-Disclosing} & \\multicolumn{3}{c}{Disclosing (Placebo)} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  " & (1) & (2) & (3) & (4) & (5) & (6) \\\\",
  " & 2013--18 & 2019--21 & 2022--23 & 2013--18 & 2019--21 & 2022--23 \\\\",
  "\\midrule",
  paste0("MNAR Emissions & ",
         fmt(p1$coef), stars(p1$pval), " & ",
         fmt(p2$coef), stars(p2$pval), " & ",
         fmt(p3$coef), stars(p3$pval), " & ",
         fmt(pl1$coef), stars(pl1$pval), " & ",
         fmt(pl2$coef), stars(pl2$pval), " & ",
         fmt(pl3$coef), stars(pl3$pval), " \\\\"),
  paste0(" & (", fmt(p1$se), ") & (", fmt(p2$se), ") & (", fmt(p3$se), ") & (",
         fmt(pl1$se), ") & (", fmt(pl2$se), ") & (", fmt(pl3$se), ") \\\\"),
  "\\midrule",
  "Controls    & Yes & Yes & Yes & Yes & Yes & Yes \\\\",
  "Time FE     & Yes & Yes & Yes & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & Yes & Yes & Yes & Yes \\\\",
  paste0("Observations & ", fmtN(p1$n), " & ", fmtN(p2$n), " & ", fmtN(p3$n),
         " & ", fmtN(pl1$n), " & ", fmtN(pl2$n), " & ", fmtN(pl3$n), " \\\\"),
  paste0("$R^2$ & ", fmt(p1$r2), " & ", fmt(p2$r2), " & ", fmt(p3$r2),
         " & -- & -- & -- \\\\"),
  "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Pooled interaction test (never-disclosing)}} \\\\",
  paste0("Baseline & \\multicolumn{2}{c}{",
         fmt(pool_base$coef), stars(pool_base$pval),
         " (", fmt(pool_base$se), ")} & & & & \\\\"),
  paste0("$\\Delta$ ESG peak & \\multicolumn{2}{c}{",
         fmt(pool_esg["Estimate"]), stars(pool_esg["Pr(>|t|)"]),
         " (", fmt(pool_esg["Std. Error"]), ")} & & & & \\\\"),
  paste0("$\\Delta$ Regulatory & \\multicolumn{2}{c}{",
         fmt(pool_reg["Estimate"]), stars(pool_reg["Pr(>|t|)"]),
         " (", fmt(pool_reg["Std. Error"]), ")} & & & & \\\\"),
  paste0("Wald $p$-value & \\multicolumn{2}{c}{",
         formatC(wald_res$p, format="f", digits=3), "} & & & & \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} The Wald test evaluates whether the ESG-peak and regulatory-period",
  "interaction coefficients are jointly zero. Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_per, file.path(outdir_tab, "returns_never_periods.tex"))
cat("returns_never_periods.tex written\n")


# =============================================================================
# BLOCK 4: DISCLOSURE SURPRISE (appendix)
# Surprise = log(disclosed_t0) - log(MNAR_predicted_t-1).
# Negative coefficient = more than expected is bad news.
#
# =============================================================================
cat("\n=== BLOCK 4: DISCLOSURE SURPRISE ===\n")

# Transition firms + event time (shared setup; see also 06_did_analysis.R)
trans <- dfm[!is.na(sc1_disclosed), .(first_disc = min(year)), by = permno]
trans <- trans[permno %in% unique(dfm[is.na(sc1_disclosed), permno])]
dfev  <- merge(dfm, trans, by = "permno", all.x = FALSE)
dfev[, etime := year - first_disc]

pred_tm1 <- dfev[etime == -1 & !is.na(sc1_adj),
                 .(sc1_adj_tm1 = mean(sc1_adj, na.rm = TRUE)), by = permno]
disc_t0  <- dfev[etime ==  0 & !is.na(sc1_disclosed),
                 .(sc1_disc_t0 = mean(sc1_disclosed, na.rm = TRUE)), by = permno]
trans2 <- merge(trans, pred_tm1, by = "permno", all.x = TRUE)
trans2 <- merge(trans2, disc_t0, by = "permno", all.x = TRUE)

trans2[sc1_disc_t0 > 0 & sc1_adj_tm1 > 0,
       surprise_raw := log(sc1_disc_t0) - log(sc1_adj_tm1)]
q_surp <- quantile(trans2$surprise_raw, probs = c(.01, .99), na.rm = TRUE)
trans2[, surprise := pmin(pmax(surprise_raw, q_surp[1]), q_surp[2])]
trans2[!is.na(surprise), surprise_z := as.numeric(scale(surprise))]

dfev2 <- merge(dfev, trans2[, .(permno, surprise_z, sc1_adj_tm1, sc1_disc_t0)],
               by = "permno", all.x = TRUE)
d_t0  <- dfev2[etime == 0 & !is.na(surprise_z)]
cat("Disclosure-year monthly obs:", nrow(d_t0),
    "| firms:", d_t0[, uniqueN(permno)], "\n")

fit_s1 <- feols(as.formula(paste("ret ~ surprise_z +", ctrl, "| yyyymm + gsector")),
                cluster = ~gsector, data = d_t0)
d_t0[, sc1_tm1_z := as.numeric(scale(log(sc1_adj_tm1 + 1)))]
fit_s2 <- feols(as.formula(paste("ret ~ surprise_z + sc1_tm1_z +", ctrl, "| yyyymm + gsector")),
                cluster = ~gsector, data = d_t0[!is.na(sc1_tm1_z)])
fit_s3 <- feols(as.formula(paste("ret ~ surprise_z +", ctrl, "| yyyymm")),
                cluster = ~gsector, data = d_t0)

get_surp <- function(fit) {
  if (!"surprise_z" %in% rownames(coeftable(fit)))
    return(list(coef=NA, se=NA, pval=NA, n=NA, r2=NA))
  extr(fit, "surprise_z")
}
s1 <- get_surp(fit_s1); s2 <- get_surp(fit_s2); s3 <- get_surp(fit_s3)

tex_surp <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Disclosure Surprise and Stock Returns at First Disclosure.}",
  "Monthly panel regressions on standardized surprise",
  "$= \\log(\\text{Disclosed}_{t_0}) - \\log(\\text{MNAR}_{t_{-1}})$.",
  "Negative coefficient indicates information content.",
  "Standard errors clustered at the sector level.}",
  "\\label{tab:returns_surprise}",
  "\\begin{tabular}{lccc}", "\\toprule",
  " & (1) & (2) & (3) \\\\",
  " & Base & + MNAR Level & Time FE Only \\\\", "\\midrule",
  paste0("Surprise & ",
         fmt(s1$coef), stars(s1$pval), " & ",
         fmt(s2$coef), stars(s2$pval), " & ",
         fmt(s3$coef), stars(s3$pval), " \\\\"),
  paste0(" & (", fmt(s1$se), ") & (", fmt(s2$se), ") & (", fmt(s3$se), ") \\\\"),
  "\\midrule",
  "Controls    & Yes & Yes & Yes \\\\",
  "Time FE     & Yes & Yes & Yes \\\\",
  "Industry FE & Yes & Yes & No  \\\\",
  paste0("Observations & ", fmtN(s1$n), " & ", fmtN(s2$n), " & ", fmtN(s3$n), " \\\\"),
  paste0("$R^2$ & ", fmt(s1$r2), " & ", fmt(s2$r2), " & ", fmt(s3$r2), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Sample restricted to disclosure year ($t=0$).",
  "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_surp, file.path(outdir_tab, "returns_surprise_main.tex"))
cat("returns_surprise_main.tex written\n")

# Surprise window: t-1, t=0, t+1
surp_window <- rbindlist(lapply(list(-1, 0, 1), function(et) {
  d <- dfev2[etime == et & !is.na(surprise_z)]
  if (nrow(d) < 50) return(NULL)
  fit <- tryCatch(
    feols(as.formula(paste("ret ~ surprise_z +", ctrl, "| yyyymm + gsector")),
          cluster = ~gsector, data = d),
    error = function(e) NULL)
  if (is.null(fit) || !"surprise_z" %in% rownames(coeftable(fit))) return(NULL)
  co <- coeftable(fit)["surprise_z", ]
  data.table(event_time=et, coef=round(co["Estimate"],3),
             se=round(co["Std. Error"],3), pval=co["Pr(>|t|)"], n=fit$nobs)
}))

tex_sw <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Disclosure Surprise: Return Effect by Event Year.}",
  "Standard errors clustered at sector level.}",
  "\\label{tab:returns_surprise_window}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  "Event Year & Coef. & (SE) & $p$-value & N \\\\", "\\midrule")
for (i in seq_len(nrow(surp_window))) {
  r <- surp_window[i]
  lab <- paste0("$t", ifelse(r$event_time >= 0, "+", ""), r$event_time, "$")
  if (r$event_time == 0)  labf <- "$t=0$ (disclosure year)"
  if (r$event_time == -1) lab <- "$t-1$ (pre-disclosure)"
  if (r$event_time ==  1) lab <- "$t+1$ (post-disclosure)"
  tex_sw <- c(tex_sw,
    paste0(lab, " & ", fmt(r$coef), stars(r$pval), " & (", fmt(r$se), ") & ",
           formatC(r$pval,format="f",digits=3), " & ", fmtN(r$n), " \\\\"))
}
tex_sw <- c(tex_sw, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_sw, file.path(outdir_tab, "returns_surprise_window.tex"))
cat("returns_surprise_window.tex written\n")


cat("\n=== 05_mechanisms.R COMPLETE ===\n")
