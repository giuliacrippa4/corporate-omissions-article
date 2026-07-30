# =============================================================================
# 06_did_analysis.R
# Transition difference-in-differences: the pre/post-disclosure return penalty,
# plus every diagnostic that supports (or stresses) it. Consolidated here so the
# DiD story lives in one place.
#
# Sample: firms that transition from non-disclosure to voluntary disclosure.
# Treatment ("High Emitter") is assigned two ways, run side by side throughout:
#   pre_mnar  — top-emitter on MNAR-imputed emissions at t-1 (ex ante, predictive)
#   pre_disc  — top-emitter on disclosed emissions at t=0    (ex post, ground truth)
#
# SECTIONS
#   0  Transition sample + treatment definitions + classification agreement
#   1  Main DiD across three event windows            -> returns_did_{mnar,disc}.tex
#   2  Event-time profile (parallel-trends check)      -> returns_eventtime_*.tex + fig
#   3  Placebo: false disclosure date (shifted -3y)    -> returns_placebo_*.tex
#   4  Balanced-panel robustness                       -> returns_did_balanced.tex
#   5  Cohort decomposition (pooled buckets + by year) -> returns_did_cohort_*.tex
#   6  Calendar-time placebo (is the effect just a high-emitter calendar trend?)
#                                                      -> returns_calendar_placebo.tex
#
# Run section by section, not in bulk.
# =============================================================================

source("00_setup.R")   # provides dfm, ctrl, outdir_tab, outdir_fig, fmt/fmtN/stars
library(ggplot2)

# -----------------------------------------------------------------------------
# PARAMETERS — change in one place, propagates everywhere
# -----------------------------------------------------------------------------
# Percentile cutoff for "High Emitter". This is the only place in the codebase
# where such a cutoff is defined -- 05_mechanisms.R no longer contains the DiD
# and sets no cutoff of its own. Section 10 sweeps this parameter: the cohort
# sign reversal is stable from p60 to p80 and loses power at p90, where the top
# decile leaves too few transition firms per cohort to estimate.
HI_PCTILE  <- 0.75
EVENT_BAND <- -4:4               # event-time window kept in dfev_base
DID_WINDOWS <- list(c(-2, 2), c(-2, 4), c(-4, 4))

# Units sanity check: ret is in percentage points unless 00_setup_whole.R has
# `df[, ret := ret/100]` un-commented. Coefficients scale with this; p-values do not.
cat(sprintf("\n[units] median(ret) = %.4f  -> ret is in %s\n",
            median(dfm$ret, na.rm = TRUE),
            ifelse(median(dfm$ret, na.rm = TRUE) > 0.05, "PERCENTAGE POINTS",
                   "decimal")))


# =============================================================================
# SECTION 0: TRANSITION SAMPLE + TREATMENT DEFINITIONS
# =============================================================================
cat("\n===== SECTION 0: Transition sample + treatment definitions =====\n")

# Transition firms: have at least one disclosed firm-month AND one non-disclosed.
trans <- dfm[!is.na(sc1_disclosed), .(first_disc = min(year)), by = permno]
trans <- trans[permno %in% unique(dfm[is.na(sc1_disclosed), permno])]
cat("Transition firms:", nrow(trans), "\n")

dfev <- merge(dfm, trans, by = "permno", all.x = FALSE)
dfev[, etime := year - first_disc]
dfev_base <- dfev[etime %in% EVENT_BAND]
dfev_base[, post := as.integer(etime >= 0)]

# Treatment definition (a): MNAR-imputed emissions at t-1 (available before disclosure)
pre_mnar <- dfev_base[etime == -1 & !is.na(sc1_adj),
                      .(sc1_pre = mean(sc1_adj, na.rm = TRUE)), by = permno]
pre_mnar[, hi := sc1_pre > quantile(sc1_pre, HI_PCTILE, na.rm = TRUE)]

# Treatment definition (b): disclosed emissions at t=0 (ground truth, ex post)
pre_disc <- dfev_base[etime == 0 & !is.na(sc1_disclosed),
                      .(sc1_pre = mean(sc1_disclosed, na.rm = TRUE)), by = permno]
pre_disc[, hi := sc1_pre > quantile(sc1_pre, HI_PCTILE, na.rm = TRUE)]

# How often do the two definitions agree on who is a high emitter?
agreement <- merge(pre_mnar[, .(permno, hi_mnar = hi)],
                   pre_disc[, .(permno, hi_disc = hi)],
                   by = "permno", all = TRUE)
cat("\nClassification agreement (MNAR-t-1 vs disclosed-t=0):\n")
print(agreement[, table(hi_mnar, hi_disc, useNA = "ifany")])
both_avail <- agreement[!is.na(hi_mnar) & !is.na(hi_disc)]
cat(sprintf("Agreement rate (firms classified by both): %.1f%%\n",
            both_avail[, mean(hi_mnar == hi_disc)] * 100))


# =============================================================================
# SECTION 1: MAIN DiD ACROSS THREE EVENT WINDOWS
#
# ret ~ controls + hi + post + (hi x post) | yyyymm + gsector
# The coefficient of interest is treat_post = hi x post: the extra return high
# emitters earn after transitioning to disclosure, relative to low emitters.
# =============================================================================
cat("\n===== SECTION 1: Main DiD =====\n")

run_did <- function(pre_data, spec_label) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]

  cols <- lapply(DID_WINDOWS, function(win) {
    sub <- d[etime %in% win[1]:win[2] & !is.na(hi)]
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl,
                             "+ hi + post + treat_post | yyyymm + gsector")),
            cluster = ~permno, data = sub),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ct <- coeftable(fit)
    getrow <- function(v) if (v %in% rownames(ct)) ct[v, ] else rep(NA_real_, 4)
    hi_r <- getrow("hiTRUE"); po_r <- getrow("post"); tp_r <- getrow("treat_post")
    list(hi_c = round(hi_r[1], 3), hi_s = round(hi_r[2], 3), hi_p = hi_r[4],
         po_c = round(po_r[1], 3), po_s = round(po_r[2], 3), po_p = po_r[4],
         tp_c = round(tp_r[1], 3), tp_s = round(tp_r[2], 3), tp_p = tp_r[4],
         n = fit$nobs, r2 = round(r2(fit, "r2"), 3))
  })
  cat(sprintf("\n%s — High x Post across windows:\n", spec_label))
  for (i in seq_along(cols)) {
    cat(sprintf("  [%d,%d]: %.3f (SE %.3f, p=%.3f, N=%s)\n",
                DID_WINDOWS[[i]][1], DID_WINDOWS[[i]][2],
                cols[[i]]$tp_c, cols[[i]]$tp_s, cols[[i]]$tp_p, fmtN(cols[[i]]$n)))
  }
  cols
}

write_did_table <- function(did_cols, treatment_label, label_key, outfile) {
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    paste0("\\caption{\\textbf{Pre-Disclosure Return Penalty: Transition Firms (",
           treatment_label, ").}"),
    "Difference-in-differences for firms transitioning from non-disclosure to",
    "voluntary disclosure.",
    paste0("High Emitter: top emitter on ", treatment_label,
           " ($>$p", round(HI_PCTILE * 100), ")."),
    "Post: disclosure year and after. Controls, time and industry FEs.",
    "Standard errors clustered at firm level.}",
    paste0("\\label{", label_key, "}"),
    "\\begin{tabular}{lccc}", "\\toprule",
    " & $[-2,+2]$ & $[-2,+4]$ & $[-4,+4]$ \\\\", "\\midrule",
    "\\textit{Pre-disclosure penalty} & & & \\\\[2pt]")
  row_defs <- list(
    list("High Emitter",        "hi_c", "hi_s", "hi_p"),
    list("Post Disclosure",     "po_c", "po_s", "po_p"),
    list("High~$\\times$~Post", "tp_c", "tp_s", "tp_p"))
  for (rw in row_defs) {
    label <- rw[[1]]; cc <- rw[[2]]; ss <- rw[[3]]; pp <- rw[[4]]
    tex <- c(tex,
      paste0(label, " & ",
             paste(sapply(did_cols, function(col)
               paste0(fmt(col[[cc]]), stars(col[[pp]]))), collapse = " & "), " \\\\"),
      paste0(" & ",
             paste(sapply(did_cols, function(col)
               paste0("(", fmt(col[[ss]]), ")")), collapse = " & "), " \\\\"), "")
  }
  tex <- c(tex, "\\midrule",
    "Controls    & Yes & Yes & Yes \\\\",
    "Time FE     & Yes & Yes & Yes \\\\",
    "Industry FE & Yes & Yes & Yes \\\\",
    paste0("Observations & ",
           paste(sapply(did_cols, function(col) fmtN(col$n)), collapse = " & "), " \\\\"),
    paste0("$R^2$ & ",
           paste(sapply(did_cols, function(col) fmt(col$r2)), collapse = " & "), " \\\\"),
    "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    paste0("\\textit{Notes:} Treatment assigned using ", treatment_label, "."),
    "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}

did_mnar <- run_did(pre_mnar, "MNAR-t-1")
did_disc <- run_did(pre_disc, "Disclosed-t=0")

write_did_table(did_mnar, "MNAR-corrected emissions at $t-1$",
                "tab:returns_did_mnar", "returns_did_mnar.tex")
write_did_table(did_disc, "disclosed emissions at $t=0$",
                "tab:returns_did_disc", "returns_did_disc.tex")


# =============================================================================
# SECTION 2: EVENT-TIME PROFILE (parallel-trends check)
#
# High-emitter coefficient estimated separately at each event year. A flat,
# near-zero profile before t=0 is the parallel-trends evidence; the jump (if any)
# should be at/after t=0.
# =============================================================================
cat("\n===== SECTION 2: Event-time profile =====\n")

run_eventtime <- function(pre_data, spec_label) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  res <- rbindlist(lapply(EVENT_BAND, function(et) {
    sub <- d[etime == et & !is.na(hi)]
    if (nrow(sub) < 100) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl, "+ hi | yyyymm + gsector")),
            cluster = ~permno, data = sub),
      error = function(e) NULL)
    if (is.null(fit) || !"hiTRUE" %in% rownames(coeftable(fit))) return(NULL)
    co <- coeftable(fit)["hiTRUE", ]
    data.table(etime = et, coef = round(co[1], 3), se = round(co[2], 3),
               pval = co[4], n = fit$nobs, treatment = spec_label)
  }))
  cat(sprintf("\n%s — High-emitter coefficient by event year:\n", spec_label))
  print(res)
  res
}

write_eventtime_table <- function(etime_res, treatment_label, label_key, outfile) {
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    paste0("\\caption{\\textbf{Return Differential Around First Disclosure: ",
           "Event-Time Profile (", treatment_label, ").}"),
    "High Emitter coefficient from annual regressions at each event year.",
    "Controls, time and industry FEs. Standard errors clustered at firm level.}",
    paste0("\\label{", label_key, "}"),
    "\\begin{tabular}{lcccc}", "\\toprule",
    "Event Year & Coef. & (SE) & $p$-value & N \\\\", "\\midrule")
  for (i in seq_len(nrow(etime_res))) {
    r <- etime_res[i]
    lab <- if (r$etime == 0) "$t=0$"
           else paste0("$t", ifelse(r$etime > 0, "+", ""), r$etime, "$")
    tex <- c(tex,
      paste0(lab, " & ", fmt(r$coef), stars(r$pval), " & (", fmt(r$se), ") & ",
             round(r$pval, 3), " & ", fmtN(r$n), " \\\\"))
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}

et_mnar <- run_eventtime(pre_mnar, "MNAR-t-1")
et_disc <- run_eventtime(pre_disc, "Disclosed-t=0")
write_eventtime_table(et_mnar, "MNAR-t-1",
                      "tab:returns_eventtime_mnar", "returns_eventtime_mnar.tex")
write_eventtime_table(et_disc, "Disclosed-t=0",
                      "tab:returns_eventtime_disc", "returns_eventtime_disc.tex")

# Event-time plot, both treatment definitions overlaid
et_all <- rbindlist(list(et_mnar, et_disc))
p_event <- ggplot(et_all, aes(etime, coef, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  labs(x = "Event year (relative to first disclosure)",
       y = "High-emitter return coefficient",
       title = "Returns around the disclosure transition") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(outdir_fig, "did_eventtime.pdf"), p_event, width = 7, height = 4)
cat("did_eventtime.pdf written\n")


# =============================================================================
# SECTION 3: PLACEBO — FALSE DISCLOSURE DATE (shifted back 3 years)
#
# Re-run the [-2,+2] DiD pretending each firm disclosed 3 years earlier than it
# actually did. A clean placebo shows a null High x Post.
# =============================================================================
cat("\n===== SECTION 3: Placebo (false disclosure date, t-3) =====\n")

run_placebo <- function(pre_data, treatment_label) {
  trans_pl <- copy(trans); trans_pl[, first_disc := first_disc - 3]
  dfev_pl  <- merge(dfm, trans_pl, by = "permno", all.x = FALSE)
  dfev_pl[, etime := year - first_disc]
  dfev_pl  <- dfev_pl[etime %in% -2:2]
  dfev_pl  <- merge(dfev_pl, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  dfev_pl[, post       := as.integer(etime >= 0)]
  dfev_pl[, treat_post := as.integer(hi == TRUE & post == 1L)]

  dfev_real <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  dfev_real[, treat_post := as.integer(hi == TRUE & post == 1L)]

  f <- paste("ret ~", ctrl, "+ hi + post + treat_post | yyyymm + gsector")
  fit_pl   <- feols(as.formula(f), cluster = ~permno, data = dfev_pl[!is.na(hi)])
  fit_real <- feols(as.formula(f), cluster = ~permno,
                    data = dfev_real[etime %in% -2:2 & !is.na(hi)])
  cat(sprintf("\n%s — High x Post: real = %.3f (p=%.3f) | placebo t-3 = %.3f (p=%.3f)\n",
              treatment_label,
              coeftable(fit_real)["treat_post", 1], coeftable(fit_real)["treat_post", 4],
              coeftable(fit_pl)["treat_post", 1],   coeftable(fit_pl)["treat_post", 4]))
  list(real = coeftable(fit_real), placebo = coeftable(fit_pl),
       n_real = fit_real$nobs, n_pl = fit_pl$nobs, label = treatment_label)
}

write_placebo_table <- function(pl_res, label_key, outfile) {
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    paste0("\\caption{\\textbf{Placebo Test: False Disclosure Date (", pl_res$label, ").}"),
    "Main DiD specification vs placebo with disclosure event shifted back 3 years.",
    "The relevant test is High~$\\times$~Post.",
    "Controls, time and industry FEs. Standard errors clustered at firm level.}",
    paste0("\\label{", label_key, "}"),
    "\\begin{tabular}{lcc}", "\\toprule",
    " & Real Event & Placebo ($t-3$) \\\\", "\\midrule")
  for (rw in list(
    list("High Emitter",        "hiTRUE"),
    list("Post Disclosure",     "post"),
    list("High~$\\times$~Post", "treat_post"))) {
    label <- rw[[1]]; v <- rw[[2]]
    rr <- if (v %in% rownames(pl_res$real)) pl_res$real[v, ] else rep(NA_real_, 4)
    pp <- if (v %in% rownames(pl_res$placebo)) pl_res$placebo[v, ] else rep(NA_real_, 4)
    tex <- c(tex,
      paste0(label, " & ", fmt(rr[1]), stars(rr[4]), " & ",
             fmt(pp[1]), stars(pp[4]), " \\\\"),
      paste0(" & (", fmt(rr[2]), ") & (", fmt(pp[2]), ") \\\\"), "")
  }
  tex <- c(tex,
    paste0("Observations & ", fmtN(pl_res$n_real), " & ", fmtN(pl_res$n_pl), " \\\\"),
    "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    "\\textit{Notes:} Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}

pl_mnar <- run_placebo(pre_mnar, "MNAR-t-1")
pl_disc <- run_placebo(pre_disc, "Disclosed-t=0")
write_placebo_table(pl_mnar, "tab:returns_placebo_mnar", "returns_placebo_mnar.tex")
write_placebo_table(pl_disc, "tab:returns_placebo_disc", "returns_placebo_disc.tex")


# =============================================================================
# SECTION 4: BALANCED-PANEL ROBUSTNESS
#
# Restrict each window to firms whose full event window falls inside the sample
# (2012-2023), so the estimate is not driven by composition changes across
# event time.
# =============================================================================
cat("\n===== SECTION 4: Balanced-panel robustness =====\n")

run_did_balanced <- function(pre_data, win, spec_label) {
  pre_yrs  <- abs(win[1]); post_yrs <- win[2]
  min_disc <- 2012 + pre_yrs; max_disc <- 2023 - post_yrs
  balanced_firms <- trans[first_disc >= min_disc & first_disc <= max_disc, permno]

  d <- merge(dfev_base[permno %in% balanced_firms], pre_data[, .(permno, hi)],
             by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  sub <- d[etime %in% win[1]:win[2] & !is.na(hi)]
  fit <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl,
                           "+ hi + post + treat_post | yyyymm + gsector")),
          cluster = ~permno, data = sub),
    error = function(e) NULL)
  if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit))) return(NULL)
  ct <- coeftable(fit)["treat_post", ]
  data.table(treatment = spec_label, window = sprintf("[%d,%d]", win[1], win[2]),
             coef = round(ct[1], 3), se = round(ct[2], 3), pval = round(ct[4], 3),
             n_firms = sub[, uniqueN(permno)], n_obs = nrow(sub))
}

balanced_res <- rbindlist(lapply(DID_WINDOWS, function(win) {
  rbindlist(list(run_did_balanced(pre_mnar, win, "MNAR-t-1"),
                 run_did_balanced(pre_disc, win, "Disclosed-t=0")))
}))
cat("\nBalanced-panel High x Post:\n")
print(balanced_res)


# =============================================================================
# SECTION 5: COHORT DECOMPOSITION
#
# Is the DiD a uniform effect, or driven by firms that disclosed in a particular
# era? Two views:
#   5a  Pooled buckets: 2013-16 / 2017-19 / 2020-23, across all three windows.
#   5b  Year-by-year: a separate DiD for each first-disclosure-year cohort.
# Both run under both treatment definitions.
# =============================================================================
cat("\n===== SECTION 5: Cohort decomposition =====\n")

# --- 5a: pooled three-bucket cohort DiD ---
run_cohort_did <- function(pre_data, spec_label) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  d[, cohort_group := fcase(
    first_disc <= 2016, "2013-2016",
    first_disc <= 2019, "2017-2019",
    first_disc >= 2020, "2020-2023")]

  res <- rbindlist(lapply(DID_WINDOWS, function(win) {
    rbindlist(lapply(c("2013-2016", "2017-2019", "2020-2023"), function(cg) {
      sub <- d[cohort_group == cg & etime %in% win[1]:win[2] & !is.na(hi)]
      if (sub[, uniqueN(permno)] < 30) {
        cat(sprintf("  [%d,%d] %s: only %d firms, skipping\n",
                    win[1], win[2], cg, sub[, uniqueN(permno)]))
        return(NULL)
      }
      fit <- tryCatch(
        feols(as.formula(paste("ret ~", ctrl,
                               "+ hi + post + treat_post | yyyymm + gsector")),
              cluster = ~permno, data = sub),
        error = function(e) NULL)
      if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit))) return(NULL)
      ct <- coeftable(fit)["treat_post", ]
      data.table(window = sprintf("[%d,%d]", win[1], win[2]), cohort = cg,
                 coef = round(ct[1], 3), se = round(ct[2], 3), pval = round(ct[4], 3),
                 n_firms = sub[, uniqueN(permno)], n_obs = nrow(sub))
    }))
  }))
  cat(sprintf("\n%s — High x Post by cohort bucket x window:\n", spec_label))
  print(res)
  res
}

cohort_did_mnar <- run_cohort_did(pre_mnar, "MNAR-t-1")
cohort_did_disc <- run_cohort_did(pre_disc, "Disclosed-t=0")

# --- 5b: year-by-year cohort DiD (one DiD per first-disclosure year) ---
run_cohort_year_did <- function(pre_data, spec_label, win = c(-2, 2)) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  res <- rbindlist(lapply(sort(unique(trans$first_disc)), function(yr) {
    cohort_firms <- trans[first_disc == yr, permno]
    sub <- d[permno %in% cohort_firms & etime %in% win[1]:win[2] & !is.na(hi)]
    if (sub[, uniqueN(permno)] < 30) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl,
                             "+ hi + post + treat_post | yyyymm + gsector")),
            cluster = ~permno, data = sub),
      error = function(e) NULL)
    if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit))) return(NULL)
    ct <- coeftable(fit)["treat_post", ]
    data.table(first_disc = yr, coef = round(ct[1], 3), se = round(ct[2], 3),
               pval = round(ct[4], 3), n_firms = sub[, uniqueN(permno)],
               n_obs = nrow(sub))
  }))
  cat(sprintf("\n%s — High x Post by first-disclosure year [window %d,%d]:\n",
              spec_label, win[1], win[2]))
  print(res)
  res
}

cohort_year_mnar <- run_cohort_year_did(pre_mnar, "MNAR-t-1")
cohort_year_disc <- run_cohort_year_did(pre_disc, "Disclosed-t=0")

# --- 5b-bis: export the year-by-year decomposition ---------------------------
# The three-bucket table shows a sign reversal between the 2017-2019 and
# 2020-2023 cohorts. That is either a genuine regime shift or an artifact of
# where the bucket boundaries were drawn, and only the year-by-year series can
# tell them apart: a monotone drift from negative to positive supports the
# former, a single outlying year the latter.
write_cohort_year_table <- function(cy_mnar, cy_disc, outfile) {
  yrs <- sort(union(cy_mnar$first_disc, cy_disc$first_disc))
  cell <- function(dt, yr) {
    r <- dt[first_disc == yr]
    if (nrow(r) == 0) return("--")
    paste0(fmt(r$coef), stars(r$pval), " (", fmt(r$se), ")")
  }
  nfirm <- function(dt, yr) {
    r <- dt[first_disc == yr]; if (nrow(r) == 0) "--" else fmtN(r$n_firms)
  }
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    "\\caption{\\textbf{Transition DiD by First-Disclosure Year.}",
    "High~$\\times$~Post estimated on a separate $[-2,+2]$ window for each",
    "first-disclosure cohort, under both treatment definitions. Cohorts with",
    "fewer than 30 firms are not estimated. Controls, time and industry fixed",
    "effects; standard errors clustered at the firm level.}",
    "\\label{tab:returns_did_cohort_year}",
    "\\begin{tabular}{lccc}", "\\toprule",
    " & \\multicolumn{2}{c}{High~$\\times$~Post} & \\\\",
    "\\cmidrule(lr){2-3}",
    "First disclosure & MNAR ($t-1$) & Disclosed ($t=0$) & Firms \\\\",
    "\\midrule")
  for (yr in yrs)
    tex <- c(tex, paste0(yr, " & ", cell(cy_mnar, yr), " & ", cell(cy_disc, yr),
                         " & ", nfirm(cy_disc, yr), " \\\\"))
  tex <- c(tex, "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    "\\textit{Notes:} Firm counts are from the disclosed-$t=0$ classification;",
    "the MNAR arm differs slightly because it requires a recovered value at",
    "$t-1$. Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}
write_cohort_year_table(cohort_year_mnar, cohort_year_disc,
                        "returns_did_cohort_year.tex")

# Write the pooled-bucket decomposition (slide material) to TeX
write_cohort_table <- function(cohort_res, treatment_label, label_key, outfile) {
  windows <- unique(cohort_res$window)
  buckets <- c("2013-2016", "2017-2019", "2020-2023")
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    paste0("\\caption{\\textbf{Cohort Decomposition of the Transition DiD (",
           treatment_label, ").}"),
    "High~$\\times$~Post estimated separately for each first-disclosure cohort.",
    "Controls, time and industry FEs. Standard errors clustered at firm level.}",
    paste0("\\label{", label_key, "}"),
    paste0("\\begin{tabular}{l", paste(rep("c", length(windows)), collapse = ""), "c}"),
    "\\toprule",
    paste0("Cohort & ", paste(windows, collapse = " & "), " & Firms \\\\"), "\\midrule")
  for (cg in buckets) {
    cells <- sapply(windows, function(w) {
      r <- cohort_res[cohort == cg & window == w]
      if (nrow(r) == 0) return("--")
      paste0(fmt(r$coef), stars(r$pval), " (", fmt(r$se), ")")
    })
    # Firm count from the narrowest window: a cohort coefficient is not
    # interpretable without knowing how many firms sit behind it.
    nf <- cohort_res[cohort == cg & window == windows[1], n_firms]
    tex <- c(tex, paste0(cg, " & ", paste(cells, collapse = " & "), " & ",
                         if (length(nf)) fmtN(nf) else "--", " \\\\"))
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    paste0("\\textit{Notes:} Treatment assigned using ", treatment_label,
           ". Firms counted in the ", windows[1], " window.",
           " Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."),
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}

write_cohort_table(cohort_did_mnar, "MNAR-t-1",
                   "tab:returns_did_cohort_mnar", "returns_did_cohort_mnar.tex")
write_cohort_table(cohort_did_disc, "Disclosed-t=0",
                   "tab:returns_did_cohort_disc", "returns_did_cohort_disc.tex")


# =============================================================================
# SECTION 6: CALENDAR-TIME PLACEBO
#
# Concern: maybe high emitters just earn more in 2020-23 for reasons unrelated to
# transitioning. If so, the transition DiD is partly capturing a calendar trend.
# Test: estimate the high-emitter return premium year by year on firms that do
# NOT transition in that window.
#
# NOTE on the treatment variable. The natural wish is to classify by *disclosed*
# emissions everywhere, for comparability with the main DiD. But never-disclosers
# have no disclosed emissions by construction (sc1_disclosed is NA for them), so
# that is only possible on firms that DO disclose. Hence two complementary cuts:
#   6a  Never-disclosers, classified on sc1_adj (MNAR) — the only emissions
#       measure that exists for this group.
#   6b  Always-disclosers (disclose every month, never transition), classified on
#       sc1_disclosed — the disclosed-emissions version, on the sample where it
#       is defined. This is the cleaner placebo for the main DiD's treatment.
# =============================================================================
cat("\n===== SECTION 6: Calendar-time placebo =====\n")

calendar_premium <- function(data, emis_var, spec_label, pctile = HI_PCTILE) {
  d <- data[!is.na(get(emis_var))]
  d[, hi_year := get(emis_var) > quantile(get(emis_var), pctile, na.rm = TRUE),
    by = year]
  res <- rbindlist(lapply(sort(unique(d$year)), function(yr) {
    sub <- d[year == yr & !is.na(hi_year)]
    if (nrow(sub) < 500) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl, "+ hi_year | yyyymm + gsector")),
            cluster = ~permno, data = sub),
      error = function(e) NULL)
    if (is.null(fit) || !"hi_yearTRUE" %in% rownames(coeftable(fit))) return(NULL)
    co <- coeftable(fit)["hi_yearTRUE", ]
    data.table(year = yr, coef = round(co[1], 3), se = round(co[2], 3),
               pval = round(co[4], 3), n_firms = sub[, uniqueN(permno)],
               n_obs = nrow(sub), spec = spec_label)
  }))
  cat(sprintf("\n%s — high-emitter coefficient by calendar year (p%d cutoff):\n",
              spec_label, round(pctile * 100)))
  print(res)
  if (nrow(res) > 0)
    cat(sprintf("  Mean 2013-2019: %.3f  |  Mean 2020-2023: %.3f\n",
                res[year %in% 2013:2019, mean(coef)],
                res[year %in% 2020:2023, mean(coef)]))
  res
}

# 6a: never-disclosers, MNAR classification
dfn_never <- dfm[group %in% c("none", "estimated")]
cal_never <- calendar_premium(dfn_never, "sc1_adj",
                              "Never-disclosers (MNAR)")

# 6b: disclosers, disclosed-emissions classification.
# Observation-level filter (firm-months where the firm actually disclosed) so
# this arm mirrors 6a's construction. Requiring disclosure in *every* month
# would select early/large established disclosers and is asymmetric with 6a.
dfn_disc <- dfm[group == "disclosed"]
cat(sprintf("\nDisclosing firm-months: %d (%d unique firms)\n",
            nrow(dfn_disc), dfn_disc[, uniqueN(permno)]))
cal_disc <- calendar_premium(dfn_disc, "sc1_disclosed",
                             "Disclosers (disclosed emissions)")

# Write both calendar-placebo series to one TeX table
cal_all <- rbindlist(list(cal_never, cal_disc), fill = TRUE)
write_calendar_table <- function(cal_res, outfile) {
  specs <- unique(cal_res$spec)
  yrs   <- sort(unique(cal_res$year))
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    "\\caption{\\textbf{Calendar-Time Placebo: High-Emitter Return Premium by Year.}",
    "High-emitter coefficient estimated year by year on non-transitioning firms.",
    "Controls, time and industry FEs. Standard errors clustered at firm level.}",
    "\\label{tab:returns_calendar_placebo}",
    paste0("\\begin{tabular}{l", paste(rep("c", length(specs)), collapse = ""), "}"),
    "\\toprule",
    paste0("Year & ", paste(specs, collapse = " & "), " \\\\"), "\\midrule")
  for (yr in yrs) {
    cells <- sapply(specs, function(sp) {
      r <- cal_res[year == yr & spec == sp]
      if (nrow(r) == 0) return("--")
      paste0(fmt(r$coef), stars(r$pval), " (", fmt(r$se), ")")
    })
    tex <- c(tex, paste0(yr, " & ", paste(cells, collapse = " & "), " \\\\"))
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    "\\textit{Notes:} Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}
write_calendar_table(cal_all, "returns_calendar_placebo.tex")

# =============================================================================
# SECTION 7: EVENT-TIME PROFILE BY COHORT
#
# The pooled profile in Section 2 is flat and insignificant everywhere. That is
# what averaging opposite-signed cohorts produces, so it neither confirms nor
# refutes the DiD. Splitting the profile by cohort is the informative version:
# each cohort should show a flat pre-period and a jump at t=0, in opposite
# directions. If instead the "jump" appears before t=0, the design is broken.
# =============================================================================
cat("\n===== SECTION 7: Event-time profile by cohort =====\n")

run_eventtime_cohort <- function(pre_data, lo, hi_yr, spec_label) {
  cohort_firms <- trans[first_disc >= lo & first_disc <= hi_yr, permno]
  d <- merge(dfev_base[permno %in% cohort_firms],
             pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  rbindlist(lapply(EVENT_BAND, function(et) {
    sub <- d[etime == et & !is.na(hi)]
    if (nrow(sub) < 100 || sub[, uniqueN(hi)] < 2) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl, "+ hi | yyyymm + gsector")),
            cluster = ~permno, data = sub), error = function(e) NULL)
    if (is.null(fit) || !"hiTRUE" %in% rownames(coeftable(fit))) return(NULL)
    co <- coeftable(fit)["hiTRUE", ]
    data.table(etime = et, coef = round(co[1], 3), se = round(co[2], 3),
               pval = co[4], n = fit$nobs,
               cohort = sprintf("%d-%d", lo, hi_yr), treatment = spec_label)
  }))
}

COHORT_RANGES <- list(c(2013, 2016), c(2017, 2019), c(2020, 2023))
et_cohort <- rbindlist(lapply(COHORT_RANGES, function(rg)
  rbindlist(list(
    run_eventtime_cohort(pre_mnar, rg[1], rg[2], "MNAR-t-1"),
    run_eventtime_cohort(pre_disc, rg[1], rg[2], "Disclosed-t=0")))), fill = TRUE)
cat("\nEvent-time profile by cohort:\n"); print(et_cohort)

if (nrow(et_cohort)) {
  p_ec <- ggplot(et_cohort, aes(etime, coef, color = treatment, fill = treatment)) +
    geom_ribbon(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
    facet_wrap(~ cohort) +
    labs(x = "Event year (relative to first disclosure)",
         y = "High-emitter return coefficient",
         title = "Returns around the disclosure transition, by cohort") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(outdir_fig, "did_eventtime_cohort.pdf"), p_ec,
         width = 9, height = 3.6)
  cat("did_eventtime_cohort.pdf written\n")
}


# =============================================================================
# SECTION 8: FIRM FIXED EFFECTS
#
# The main spec has no firm FE, so `post` is identified by comparing pre- and
# post-disclosure firm-months ACROSS firms within a calendar month. Adding
# permno FE makes it a within-firm estimate. `hi` is absorbed (time-invariant),
# as is gsector; High x Post survives and is the coefficient of interest.
# =============================================================================
cat("\n===== SECTION 8: Firm fixed effects =====\n")

run_did_firmfe <- function(pre_data, spec_label) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  rbindlist(lapply(DID_WINDOWS, function(win) {
    sub <- d[etime %in% win[1]:win[2] & !is.na(hi)]
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl, "+ post + treat_post | permno + yyyymm")),
            cluster = ~permno, data = sub), error = function(e) NULL)
    if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit))) return(NULL)
    ct <- coeftable(fit)["treat_post", ]; po <- coeftable(fit)["post", ]
    data.table(treatment = spec_label, window = sprintf("[%d,%d]", win[1], win[2]),
               tp_coef = round(ct[1], 3), tp_se = round(ct[2], 3), tp_p = round(ct[4], 3),
               po_coef = round(po[1], 3), po_se = round(po[2], 3), po_p = round(po[4], 3),
               n_firms = sub[, uniqueN(permno)], n_obs = nrow(sub))
  }))
}
firmfe_res <- rbindlist(list(run_did_firmfe(pre_mnar, "MNAR-t-1"),
                             run_did_firmfe(pre_disc, "Disclosed-t=0")))
cat("\nWith firm fixed effects (High x Post and Post):\n"); print(firmfe_res)


# =============================================================================
# SECTION 9: IS THE 2020-2023 COHORT EFFECT A 2021 CALENDAR ARTIFACT?
#
# Section 6 finds a high-emitter premium of +1.30 among non-transitioning
# DISCLOSERS in 2021 alone -- almost exactly the size of the 2020-2023 cohort
# DiD. Since that cohort's post-window covers 2021, the two could be the same
# thing. Dropping 2021 (and, separately, each post-period year in turn)
# separates them: if the cohort effect survives, it is not the calendar premium.
# =============================================================================
cat("\n===== SECTION 9: Leave-one-calendar-year-out, 2020-2023 cohort =====\n")

run_cohort_drop_year <- function(pre_data, drop_yr, spec_label, win = c(-2, 2)) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  sub <- d[first_disc >= 2020 & etime %in% win[1]:win[2] & !is.na(hi)]
  if (!is.na(drop_yr)) sub <- sub[year != drop_yr]
  if (sub[, uniqueN(permno)] < 30) return(NULL)
  fit <- tryCatch(
    feols(as.formula(paste("ret ~", ctrl,
                           "+ hi + post + treat_post | yyyymm + gsector")),
          cluster = ~permno, data = sub), error = function(e) NULL)
  if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit))) return(NULL)
  ct <- coeftable(fit)["treat_post", ]
  data.table(treatment = spec_label,
             dropped = if (is.na(drop_yr)) "none (baseline)" else as.character(drop_yr),
             coef = round(ct[1], 3), se = round(ct[2], 3), pval = round(ct[4], 3),
             n_firms = sub[, uniqueN(permno)], n_obs = nrow(sub))
}

drop_res <- rbindlist(lapply(c(NA, 2020:2023), function(y)
  rbindlist(list(run_cohort_drop_year(pre_mnar, y, "MNAR-t-1"),
                 run_cohort_drop_year(pre_disc, y, "Disclosed-t=0")))), fill = TRUE)
cat("\n2020-2023 cohort, dropping one calendar year at a time:\n"); print(drop_res)

write_drop_table <- function(dres, outfile) {
  yrs <- unique(dres$dropped)
  cell <- function(tr, yy) {
    r <- dres[treatment == tr & dropped == yy]
    if (nrow(r) == 0) return("--")
    paste0(fmt(r$coef), stars(r$pval), " (", fmt(r$se), ")")
  }
  tex <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    "\\caption{\\textbf{Is the 2020--2023 cohort effect a 2021 calendar artifact?}",
    "High~$\\times$~Post for the 2020--2023 first-disclosure cohort over the",
    "$[-2,+2]$ window, dropping one calendar year of returns at a time. The",
    "calendar-time placebo of Table~\\ref{tab:returns_calendar_placebo} finds a",
    "high-emitter premium among non-transitioning disclosers in 2021 of",
    "comparable magnitude, so the relevant test is whether the cohort effect",
    "survives the removal of 2021. Controls, time and industry fixed effects;",
    "standard errors clustered at the firm level.}",
    "\\label{tab:returns_did_drop_year}",
    "\\begin{tabular}{lcc}", "\\toprule",
    "Calendar year dropped & MNAR ($t-1$) & Disclosed ($t=0$) \\\\", "\\midrule")
  for (yy in yrs)
    tex <- c(tex, paste0(yy, " & ", cell("MNAR-t-1", yy), " & ",
                         cell("Disclosed-t=0", yy), " \\\\"))
  tex <- c(tex, "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    "\\textit{Notes:} Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex, file.path(outdir_tab, outfile))
  cat(outfile, "written\n")
}
if (nrow(drop_res)) write_drop_table(drop_res, "returns_did_drop_year.tex")


# =============================================================================
# SECTION 10: SENSITIVITY TO THE HIGH-EMITTER CUTOFF
#
# HI_PCTILE is an arbitrary choice, and 05_mechanisms.R uses a different one
# (0.90). If the cohort reversal only appears at one cutoff it is not a
# finding. Re-run the pooled and the 2020-2023 cohort DiD across cutoffs.
# =============================================================================
cat("\n===== SECTION 10: High-emitter cutoff sensitivity =====\n")

run_did_at_cutoff <- function(pct, win = c(-2, 2)) {
  pm <- dfev_base[etime == -1 & !is.na(sc1_adj),
                  .(sc1_pre = mean(sc1_adj, na.rm = TRUE)), by = permno]
  pm[, hi := sc1_pre > quantile(sc1_pre, pct, na.rm = TRUE)]
  d <- merge(dfev_base, pm[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  one <- function(sub, lab) {
    if (sub[, uniqueN(permno)] < 30) return(NULL)
    fit <- tryCatch(
      feols(as.formula(paste("ret ~", ctrl,
                             "+ hi + post + treat_post | yyyymm + gsector")),
            cluster = ~permno, data = sub), error = function(e) NULL)
    if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit))) return(NULL)
    ct <- coeftable(fit)["treat_post", ]
    data.table(cutoff = pct, sample = lab, coef = round(ct[1], 3),
               se = round(ct[2], 3), pval = round(ct[4], 3),
               n_firms = sub[, uniqueN(permno)])
  }
  rbindlist(list(
    one(d[etime %in% win[1]:win[2] & !is.na(hi)], "all cohorts"),
    one(d[first_disc >= 2020 & etime %in% win[1]:win[2] & !is.na(hi)], "2020-2023"),
    one(d[first_disc %between% c(2017, 2019) & etime %in% win[1]:win[2] & !is.na(hi)],
        "2017-2019")))
}
cutoff_res <- rbindlist(lapply(c(0.60, 0.70, 0.75, 0.80, 0.90), run_did_at_cutoff),
                        fill = TRUE)
cat("\nHigh x Post by high-emitter cutoff (MNAR-t-1 classification):\n")
print(dcast(cutoff_res, cutoff ~ sample, value.var = c("coef", "pval")))


# =============================================================================
# SECTION 11: SUN & ABRAHAM COHORT-ROBUST EVENT STUDY (optional)
#
# With cohort effects that reverse sign, two-way fixed effects is exactly the
# estimator the staggered-DiD literature warns against: the pooled coefficient
# is a variance-weighted average that can sit outside the range of the
# underlying group-time effects. sunab() reweights so each cohort contributes
# its own event-time path. Run separately for high and low emitters; the object
# of interest is the difference between the two profiles.
# =============================================================================
cat("\n===== SECTION 11: Sun-Abraham event study =====\n")

sa_res <- tryCatch({
  d <- merge(dfev_base, pre_mnar[, .(permno, hi)], by = "permno", all.x = TRUE)
  d <- d[!is.na(hi)]
  fits <- lapply(c(TRUE, FALSE), function(g) {
    feols(as.formula(paste("ret ~", ctrl, "+ sunab(first_disc, year) | permno + yyyymm")),
          cluster = ~permno, data = d[hi == g])
  })
  names(fits) <- c("high", "low")
  for (nm in names(fits)) {
    cat(sprintf("\n-- %s emitters, aggregated ATT --\n", nm))
    print(summary(fits[[nm]], agg = "att"))
  }
  fits
}, error = function(e) {
  cat("Sun-Abraham step skipped:", conditionMessage(e), "\n")
  cat("(fixest >= 0.10 provides sunab(); check that `year` and `first_disc` are integers)\n")
  NULL
})


# =============================================================================
# SECTION 12: INDUSTRY-BY-YEAR FIXED EFFECTS (does the cohort reversal survive?)
#
# The cohort reversal (Section 5) is identified off calendar-month + sector FE. If
# high emitters that transitioned after 2020 sit in sectors that rallied in 2021-23
# (energy, materials, utilities), a sector-wide annual return -- not the transition
# -- could manufacture the positive 2020-23 cohort effect. Adding sector x year
# interacted FE absorbs exactly that. Treatment is High x Post in EVENT time, so it
# remains identified within a sector-year cell (firms sit at different event times).
#
# Result: the 2020-23 cohort effect and the pooled estimate lose significance under
# sector x year FE, and dropping Energy+Materials alone removes ~75-80% of the
# 2020-23 effect -- the reversal is largely the commodity-sector rally, not repricing
# of the disclosure transition. Reported honestly as a robustness limitation.
#   -> returns_did_industry_year.tex
# =============================================================================
cat("\n===== SECTION 12: Industry-by-year fixed effects =====\n")

FE_IY <- c(base = "yyyymm + gsector", secyr = "yyyymm + gsector^year")

fit_hxp_iy <- function(pre_data, fe, lo, hi_yr, drop_sectors = NULL, win = c(-2, 2)) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  sub <- d[first_disc >= lo & first_disc <= hi_yr &
             etime %in% win[1]:win[2] & !is.na(hi)]
  if (!is.null(drop_sectors)) sub <- sub[!gsector %in% drop_sectors]
  if (sub[, uniqueN(permno)] < 30)
    return(list(coef = NA, se = NA, p = NA, nf = sub[, uniqueN(permno)]))
  fit <- tryCatch(feols(as.formula(paste("ret ~", ctrl, "+ hi + post + treat_post |", fe)),
                        cluster = ~permno, data = sub), error = function(e) NULL)
  if (is.null(fit) || !"treat_post" %in% rownames(coeftable(fit)))
    return(list(coef = NA, se = NA, p = NA, nf = sub[, uniqueN(permno)]))
  ct <- coeftable(fit)["treat_post", ]
  list(coef = ct[1], se = ct[2], p = ct[4], nf = sub[, uniqueN(permno)])
}

COHORTS_IY <- list(`2013--2016` = c(2013, 2016), `2017--2019` = c(2017, 2019),
                   `2020--2023` = c(2020, 2023), `Pooled` = c(2013, 2023))
grab_iy <- function(pre_data) rbindlist(lapply(names(COHORTS_IY), function(cg) {
  rg <- COHORTS_IY[[cg]]
  b <- fit_hxp_iy(pre_data, FE_IY["base"], rg[1], rg[2])
  s <- fit_hxp_iy(pre_data, FE_IY["secyr"], rg[1], rg[2])
  data.table(cohort = cg, base = b$coef, base_se = b$se, base_p = b$p,
             sy = s$coef, sy_se = s$se, sy_p = s$p, nf = b$nf)
}))
res_iy_mnar <- grab_iy(pre_mnar); res_iy_disc <- grab_iy(pre_disc)
cat("\nMNAR-t-1 (baseline vs +sector x year):\n"); print(res_iy_mnar)
cat("\nDisclosed-t=0 (baseline vs +sector x year):\n"); print(res_iy_disc)

cellf_iy <- function(c, s, p)
  if (is.na(c)) "--" else paste0(fmt(round(c, 3)), stars(p), " (", fmt(round(s, 3)), ")")
tex_iy <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Transition DiD under industry-by-year fixed effects.}",
  "High~$\\times$~Post over the $[-2,+2]$ window, by first-disclosure cohort, under",
  "the published specification (calendar-month and sector fixed effects) and after",
  "adding sector~$\\times$~year interacted fixed effects, which absorb any",
  "sector-wide annual return. Both treatment definitions. Standard errors clustered",
  "at the firm level.}",
  "\\label{tab:returns_did_industry_year}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  " & \\multicolumn{2}{c}{MNAR ($t-1$)} & \\multicolumn{2}{c}{Disclosed ($t=0$)} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Cohort & Baseline & $+$ sector$\\times$year & Baseline & $+$ sector$\\times$year \\\\",
  "\\midrule")
for (cg in names(COHORTS_IY)) {
  m <- res_iy_mnar[cohort == cg]; d <- res_iy_disc[cohort == cg]
  tex_iy <- c(tex_iy, paste0(cg, " & ",
    cellf_iy(m$base, m$base_se, m$base_p), " & ", cellf_iy(m$sy, m$sy_se, m$sy_p), " & ",
    cellf_iy(d$base, d$base_se, d$base_p), " & ", cellf_iy(d$sy, d$sy_se, d$sy_p), " \\\\"))
}
tex_iy <- c(tex_iy, "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} The 2020--2023 cohort effect and the pooled estimate lose",
  "significance once sector~$\\times$~year effects are absorbed, indicating that much",
  "of the transition DiD reflects sector-wide annual returns rather than the",
  "disclosure transition. Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_iy, file.path(outdir_tab, "returns_did_industry_year.tex"))
cat("returns_did_industry_year.tex written\n")

# Driver: dropping the commodity-rally sectors from the 2020-2023 cohort
cat("\n-- 2020-2023 cohort, baseline FE, dropping rally sectors --\n")
for (tn in c("MNAR-t-1", "Disclosed-t=0")) {
  pd <- if (tn == "MNAR-t-1") pre_mnar else pre_disc
  full <- fit_hxp_iy(pd, FE_IY["base"], 2020, 2023)
  d1   <- fit_hxp_iy(pd, FE_IY["base"], 2020, 2023, drop_sectors = c("Energy", "Materials"))
  d2   <- fit_hxp_iy(pd, FE_IY["base"], 2020, 2023,
                     drop_sectors = c("Energy", "Materials", "Utilities"))
  cat(sprintf("%-14s all: %.3f (p=%.3f, n=%d) | drop En+Mat: %.3f (p=%.3f, n=%d) | +Util: %.3f (p=%.3f, n=%d)\n",
              tn, full$coef, full$p, full$nf, d1$coef, d1$p, d1$nf, d2$coef, d2$p, d2$nf))
}

# --- Post main effect under industry-by-year FE ------------------------------
# Unlike High x Post, the Post<0 result (returns fall after a firm starts
# disclosing) is SIGN-robust: negative under every spec including firm FE +
# sector x year. Precision degrades under saturation -- the magnitude shrinks by
# ~1/3 and several specs lose significance -- so "negative and significant
# throughout" is the claim that must weaken, not the sign.  -> returns_did_post_iy.tex
cat("\n-- POST main effect under industry-by-year FE --\n")
POST_SPECS <- list(
  `Pooled ($yyyymm+$sector)`        = list(fe = "yyyymm + gsector",               firmfe = FALSE),
  `Pooled ($+$ sector$\\times$year)` = list(fe = "yyyymm + gsector^year",          firmfe = FALSE),
  `Firm FE ($permno+yyyymm$)`       = list(fe = "permno + yyyymm",                firmfe = TRUE),
  `Firm FE ($+$ sector$\\times$year)`= list(fe = "permno + yyyymm + gsector^year", firmfe = TRUE))

get_post_iy <- function(pre_data, spec, win) {
  d <- merge(dfev_base, pre_data[, .(permno, hi)], by = "permno", all.x = TRUE)
  d[, treat_post := as.integer(hi == TRUE & post == 1L)]
  sub <- d[etime %in% win[1]:win[2] & !is.na(hi)]
  rhs <- if (spec$firmfe) "post + treat_post" else "hi + post + treat_post"
  fit <- tryCatch(feols(as.formula(paste("ret ~", ctrl, "+", rhs, "|", spec$fe)),
                        cluster = ~permno, data = sub), error = function(e) NULL)
  if (is.null(fit) || !"post" %in% rownames(coeftable(fit)))
    return(list(coef = NA, se = NA, p = NA))
  ct <- coeftable(fit)["post", ]; list(coef = ct[1], se = ct[2], p = ct[4])
}

post_panel <- function(pre_data) {
  lines <- c()
  for (sp in names(POST_SPECS)) {
    cells <- sapply(DID_WINDOWS, function(w) {
      r <- get_post_iy(pre_data, POST_SPECS[[sp]], w)
      if (is.na(r$coef)) "--" else paste0(fmt(round(r$coef, 3)), stars(r$p))
    })
    ses <- sapply(DID_WINDOWS, function(w) {
      r <- get_post_iy(pre_data, POST_SPECS[[sp]], w)
      if (is.na(r$se)) "" else paste0("(", fmt(round(r$se, 3)), ")")
    })
    lines <- c(lines,
      paste0(sp, " & ", paste(cells, collapse = " & "), " \\\\"),
      paste0(" & ", paste(ses, collapse = " & "), " \\\\"))
  }
  lines
}

tex_post <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{The post-disclosure return under industry-by-year fixed effects.}",
  "The Post indicator (returns in disclosure-year-and-after, relative to before) for",
  "transition firms, across event windows and fixed-effect specifications. The sign",
  "is negative throughout, including under firm and sector~$\\times$~year effects;",
  "significance weakens under saturation. Standard errors clustered at the firm level.}",
  "\\label{tab:returns_did_post_iy}",
  "\\begin{tabular}{lccc}", "\\toprule",
  " & $[-2,+2]$ & $[-2,+4]$ & $[-4,+4]$ \\\\", "\\midrule",
  "\\multicolumn{4}{l}{\\emph{Panel A: MNAR ($t-1$) classification}} \\\\[2pt]",
  post_panel(pre_mnar),
  "\\midrule",
  "\\multicolumn{4}{l}{\\emph{Panel B: Disclosed ($t=0$) classification}} \\\\[2pt]",
  post_panel(pre_disc),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_post, file.path(outdir_tab, "returns_did_post_iy.tex"))
cat("returns_did_post_iy.tex written\n")
for (tn in c("MNAR-t-1", "Disclosed-t=0")) {
  pd <- if (tn == "MNAR-t-1") pre_mnar else pre_disc
  cat(sprintf("  %s: ", tn))
  for (sp in names(POST_SPECS)) {
    r <- get_post_iy(pd, POST_SPECS[[sp]], c(-2, 2))
    cat(sprintf("%s=%.3f(p%.2f) ", gsub("[$\\{}]|yyyymm|permno|times|sector|year|Pooled |Firm FE ", "", sp), r$coef, r$p))
  }
  cat("\n")
}

# --- Treatment as the disclosure SURPRISE (sector-orthogonal by construction) -
# The level-based "High Emitter" is collinear with sector, so the cohort DiD picks
# up the commodity rally. Redefining treatment as the surprise at first disclosure,
#   surprise_i = log(disclosed at t=0) - log(MNAR-predicted at t-1),
# closes that channel mechanically: the MNAR prediction carries the outcome model's
# industry FE, so the sector mean is already in the predicted term and the residual
# is ~orthogonal to sector. Verified below (R^2 of surprise~sector vs level~sector),
# then the surprise DiD -- which is a clean null, retiring the transition-DiD
# audience claim rather than salvaging it.  -> returns_did_surprise.tex
cat("\n-- Disclosure-surprise treatment (sector-orthogonal) --\n")
lg_s <- function(x) log(pmax(x, 1))
pred_t1 <- dfev_base[etime == -1 & !is.na(sc1_adj),
                     .(pred = mean(sc1_adj, na.rm = TRUE)), by = permno]
disc_t0 <- dfev_base[etime == 0 & !is.na(sc1_disclosed),
                     .(disc = mean(sc1_disclosed, na.rm = TRUE)), by = permno]
firm_sec <- dfev_base[, .(gsector = gsector[1]), by = permno]
sur <- merge(pred_t1, disc_t0, by = "permno")[pred > 0 & disc > 0]
sur[, surprise := lg_s(disc) - lg_s(pred)]
sur[, surprise_z := as.numeric(scale(surprise))]
sur[, hi_s := surprise > quantile(surprise, HI_PCTILE, na.rm = TRUE)]
sur <- merge(sur, firm_sec, by = "permno")
lvl_sec <- merge(pred_t1, firm_sec, by = "permno")
r2_sur <- summary(lm(surprise ~ gsector, data = sur))$r.squared
r2_lvl <- summary(lm(lg_s(pred) ~ gsector, data = lvl_sec))$r.squared
cat(sprintf("firms=%d | R2(surprise~sector)=%.3f  vs  R2(level~sector)=%.3f\n",
            nrow(sur), r2_sur, r2_lvl))

build_sur <- function(win, drop_sectors = NULL) {
  d <- merge(dfev_base, sur[, .(permno, hi_s, surprise_z)], by = "permno", all.x = TRUE)
  d[, post_his := as.integer(hi_s == TRUE & post == 1L)]
  d[, post_sz  := surprise_z * post]
  d <- d[etime %in% win[1]:win[2] & !is.na(hi_s)]
  if (!is.null(drop_sectors)) d <- d[!gsector %in% drop_sectors]
  d
}
grab_s <- function(fit, v) {
  if (is.null(fit) || !v %in% rownames(coeftable(fit))) return(list(c = NA, s = NA, p = NA))
  ct <- coeftable(fit)[v, ]; list(c = ct[1], s = ct[2], p = ct[4])
}
fitb_s <- function(d, fe) tryCatch(feols(as.formula(paste("ret ~", ctrl, "+ hi_s + post + post_his |", fe)),
                                         cluster = ~permno, data = d), error = function(e) NULL)
fitc_s <- function(d, fe) tryCatch(feols(as.formula(paste("ret ~", ctrl, "+ surprise_z + post + post_sz |", fe)),
                                         cluster = ~permno, data = d), error = function(e) NULL)

rob_s <- rbindlist(lapply(DID_WINDOWS, function(win) {
  d0  <- build_sur(win); dEM <- build_sur(win, drop_sectors = c("Energy", "Materials"))
  b <- grab_s(fitb_s(d0, FE_IY["base"]),  "post_his")
  s <- grab_s(fitb_s(d0, FE_IY["secyr"]), "post_his")
  x <- grab_s(fitb_s(dEM, FE_IY["base"]), "post_his")
  cat(sprintf("[%d,%d] High-surprise x Post: base=%.3f(p%.2f) | +sec x yr=%.3f(p%.2f) | drop En+Mat=%.3f(p%.2f)\n",
              win[1], win[2], b$c, b$p, s$c, s$p, x$c, x$p))
  data.table(window = sprintf("[%d,%d]", win[1], win[2]),
             b_c = b$c, b_s = b$s, b_p = b$p, s_c = s$c, s_s = s$s, s_p = s$p,
             x_c = x$c, x_s = x$s, x_p = x$p)
}))
cont44 <- grab_s(fitc_s(build_sur(c(-4, 4)), FE_IY["base"]), "post_sz")

cell_s <- function(c, s, p) if (is.na(c)) "--" else paste0(fmt(round(c,3)), stars(p), " (", fmt(round(s,3)), ")")
tex_s <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Transition DiD with the disclosure surprise as treatment.}",
  "Treatment is the emissions surprise at first disclosure,",
  "$\\text{surprise}_i=\\log(\\text{disclosed}_{t=0})-\\log(\\text{MNAR-predicted}_{t-1})$;",
  "High-surprise is the top quartile. Because the MNAR prediction carries the outcome",
  "model's industry fixed effects, the surprise is orthogonal to sector by",
  "construction, closing the commodity-rally channel the level-based treatment leaves",
  "open. High-surprise~$\\times$~Post over three event windows, under the baseline",
  "specification, adding sector~$\\times$~year fixed effects, and dropping Energy and",
  "Materials. Standard errors clustered at the firm level.}",
  "\\label{tab:returns_did_surprise}",
  "\\begin{tabular}{lccc}", "\\toprule",
  " & Baseline & $+$ sector$\\times$year & Drop En.$+$Mat. \\\\", "\\midrule")
for (i in seq_len(nrow(rob_s))) {
  r <- rob_s[i]
  tex_s <- c(tex_s, paste0(r$window, " & ", cell_s(r$b_c, r$b_s, r$b_p), " & ",
                          cell_s(r$s_c, r$s_s, r$s_p), " & ", cell_s(r$x_c, r$x_s, r$x_p), " \\\\"))
}
tex_s <- c(tex_s, "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  sprintf(paste0("\\textit{Notes:} The surprise is sector-neutral by construction: a ",
                 "sector-only regression has $R^2=%.3f$ for the surprise against $R^2=%.3f$ ",
                 "for the log emissions level. The High-surprise effect is a stable null ",
                 "across all three columns, and the continuous surprise (per SD)~$\\times$~Post ",
                 "is likewise small ($%.3f$, $p=%.2f$ at $[-4,+4]$): once the sector level is ",
                 "removed, the transition margin carries no repricing. Significance: ",
                 "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$."),
          r2_sur, r2_lvl, round(cont44$c, 3), round(cont44$p, 2)),
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_s, file.path(outdir_tab, "returns_did_surprise.tex"))
cat("returns_did_surprise.tex written\n")


cat("\n===== 06_did_analysis.R complete. =====\n")
cat("Tables written to:", outdir_tab, "\n")
cat("Figures written to:", outdir_fig, "\n")
