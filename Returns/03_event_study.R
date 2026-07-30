# =============================================================================
# 03_event_study.R
# Event study around the SEC climate disclosure rule proposal (2022-03-21).
#
# Merges the four event-study scripts into one pipeline:
#   - Main event vs placebo (2019-03-21, no climate-policy news)
#   - Heteroskedasticity-robust SEs with sector FEs
#   - CAR decile plot around the main event
#   - MAR vs MNAR horse race (same sample)
#   - Trucost vs MNAR horse race (same sample)
#
# Sample: firms NOT DISCLOSING AS OF THE PRE-EVENT YEAR (2021 for the main
# event, 2018 for the placebo) — i.e. firms for which the market had no
# firm-reported emissions figure when the announcement landed. This is a
# single-year condition, not "never discloses"; a subset of these firms carry a
# Trucost vendor estimate, which the horse race below exploits.
#
# Tables written:
#   returns_event_study.tex                        — MAIN TEXT: event + placebo, robust SEs
#   returns_event_study_mar_vs_mnar_main.tex       — MAR vs MNAR, main
#   returns_event_study_mar_vs_mnar_plac.tex       — MAR vs MNAR, placebo
#   returns_event_study_trucost_vs_mnar_main.tex   — Trucost vs MNAR, main
#   returns_event_study_trucost_vs_mnar_plac.tex   — Trucost vs MNAR, placebo
#
# Figure written:
#   returns_event_study_car.pdf — CAR [-5,+5] by MNAR decile around main event
#
# Data: df_lm_avg_baseline.parquet (firm-year MNAR panel) + WRDS daily CRSP.
# =============================================================================

library(data.table)
library(arrow)
library(fixest)
library(RPostgres)
library(DBI)
library(ggplot2)

setwd("/Users/giulia/Documents/GitHub/corporate-omissions/Returns")
outdir_tab <- "../data/outputs/tables"
outdir_fig <- "../data/outputs/figures"
dir.create(outdir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# TABLE PRESENTATION
# Human-readable names for every variable that reaches a table, set once and
# inherited by all etable() calls below. Without this the tables print raw
# column names (log_sc1_adj, car_short, gsector), which are not publishable.
# =============================================================================
setFixest_dict(c(
  car_short    = "CAR $[-1,+1]$",
  car_medium   = "CAR $[-5,+5]$",
  car_long     = "CAR $[-30,+30]$",
  log_sc1_adj  = "Log recovered emissions (MNAR)",
  log_sc1_mar  = "Log imputed emissions (MAR)",
  log_sc1_truc = "Log vendor emissions (Trucost)",
  log_at       = "Log assets",
  gsector      = "Sector"
))
# depvar row is redundant once `headers` carries the window; drop Within-R2.
setFixest_etable(
  depvar      = FALSE,
  fitstat     = ~ n + r2,
  digits      = 4,
  digits.stats = 3,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10)
)


# =============================================================================
# CONFIG
# =============================================================================
EVENT_DATE_MAIN    <- as.Date("2022-03-21")   # SEC proposal
PRE_YEAR_MAIN      <- 2021
EVENT_DATE_PLACEBO <- as.Date("2019-03-21")   # same calendar month, no climate news
PRE_YEAR_PLACEBO   <- 2018

ESTIMATION_WINDOW <- c(-250L, -30L)           # market-model estimation
EVENT_WINDOWS <- list(
  short  = c(-1L,  1L),
  medium = c(-5L,  5L),
  long   = c(-30L, 30L)
)
MIN_EST_OBS <- 100L

# Sample definition: firms NOT DISCLOSING in pre_year. Note this is a
# single-year condition, not "never discloses over the sample" — a firm that
# disclosed in an earlier year and stopped is included, which is the correct
# population for this test (what matters is whether the market had a reported
# figure at the announcement).
#   "none"           — neither disclosed nor vendor-estimated in pre_year
#   "none_estimated" — not disclosed; vendor estimate allowed (default; required
#                      for the Trucost horse race, which needs sc1_estimated)
SAMPLE_GROUPS <- "none_estimated"

WRDS_USER <- "gcrippa4"


# =============================================================================
# LOAD EMISSIONS PANEL ONCE, OPEN WRDS
# =============================================================================
cat("Loading MNAR panel...\n")
emis <- as.data.table(
  read_parquet("/Users/giulia/Documents/GitHub/corporate-omissions/data/processed/df_lm_avg_baseline_wti_gind_exchg.parquet")
)
setnames(emis, "fyear", "year")

cat("Connecting to WRDS...\n")
wrds <- dbConnect(
  Postgres(),
  host    = "wrds-pgdata.wharton.upenn.edu",
  port    = 9737,
  dbname  = "wrds",
  sslmode = "require",
  user    = WRDS_USER
  # password via ~/.pgpass
)


# =============================================================================
# HELPERS
# =============================================================================
build_nondisclosing_sample <- function(pre_year) {
  d <- emis[year == pre_year]
  d[, group := fcase(
    !is.na(sc1_disclosed),                        "disclosed",
    is.na(sc1_disclosed) & !is.na(sc1_estimated), "estimated",
    default = "none")]
  keep <- switch(SAMPLE_GROUPS,
                 "none"           = "none",
                 "none_estimated" = c("none", "estimated"),
                 stop("SAMPLE_GROUPS must be 'none' or 'none_estimated'"))
  d <- d[group %in% keep,
         .(gvkey, gsector, at_pre = `at`, sc1_adj_pre = sc1_adj)
        ][!is.na(sc1_adj_pre) & sc1_adj_pre > 0]
  d[, log_sc1_adj := log(sc1_adj_pre)]
  d[, log_at      := log(pmax(at_pre, 1e-6))]
  d[, gvkey       := as.integer(gvkey)]
  d
}

link_permno <- function(gvkeys, event_date) {
  gv_sql <- paste(sprintf("'%06d'", gvkeys), collapse = ",")
  ccm <- as.data.table(dbGetQuery(wrds, sprintf("
    SELECT gvkey, lpermno AS permno, linktype, linkprim, linkdt, linkenddt
    FROM crsp.ccmxpf_linktable
    WHERE gvkey IN (%s)
      AND linktype IN ('LU','LC','LN','LS','LX','LD')
      AND linkprim IN ('P','C')
  ", gv_sql)))
  ccm[is.na(linkenddt), linkenddt := as.Date("2099-12-31")]
  ccm <- ccm[event_date %between% list(linkdt, linkenddt)]
  ccm[, gvkey := as.integer(gvkey)]
  unique(ccm[, .(gvkey, permno)])
}

pull_daily_returns <- function(permnos, date_lo, date_hi) {
  pn_sql <- paste(sprintf("'%d'", permnos), collapse = ",")
  dsf <- as.data.table(dbGetQuery(wrds, sprintf("
    SELECT permno, date, ret FROM crsp.dsf
    WHERE permno IN (%s) AND date BETWEEN '%s' AND '%s'
  ", pn_sql, date_lo, date_hi)))
  dsi <- as.data.table(dbGetQuery(wrds, sprintf("
    SELECT date, vwretd AS mktret FROM crsp.dsi
    WHERE date BETWEEN '%s' AND '%s'
  ", date_lo, date_hi)))
  r <- merge(dsf, dsi, by = "date")
  r[, `:=`(date   = as.Date(date),
           ret    = as.numeric(ret),
           mktret = as.numeric(mktret))]
  r[!is.na(ret) & !is.na(mktret)]
}

compute_cars <- function(ret, event_date) {
  setorder(ret, permno, date)
  trading_days <- sort(unique(ret$date))
  t0_idx  <- which(trading_days == max(trading_days[trading_days <= event_date]))
  day_idx <- data.table(date  = trading_days,
                        t_rel = seq_along(trading_days) - t0_idx)
  r <- merge(ret, day_idx, by = "date")

  est_mm <- function(g) {
    d <- g[t_rel %between% ESTIMATION_WINDOW]
    if (nrow(d) < MIN_EST_OBS) return(NULL)
    fit <- lm(ret ~ mktret, data = d)
    co  <- coef(fit)
    data.table(alpha = co[1], beta = co[2],
               sigma_e = sd(residuals(fit)), n_est = nrow(d))
  }
  mm <- r[, est_mm(.SD), by = permno]
  r  <- merge(r, mm, by = "permno")
  r[, ar := ret - (alpha + beta * mktret)]

  cars <- lapply(names(EVENT_WINDOWS), function(w) {
    lo <- EVENT_WINDOWS[[w]][1]; hi <- EVENT_WINDOWS[[w]][2]
    d  <- r[t_rel %between% c(lo, hi),
            .(car = sum(ar, na.rm = TRUE), n_days = .N), by = permno]
    setnames(d, c("car", "n_days"),
             c(paste0("car_", w), paste0("n_", w)))
    d
  })
  list(ar_panel = r,
       cars     = Reduce(function(a, b) merge(a, b, by = "permno", all = TRUE), cars))
}

run_event <- function(event_date, pre_year, label) {
  cat(sprintf("\n--- Event: %s (%s) ---\n", label, event_date))
  never <- build_nondisclosing_sample(pre_year)
  cat(sprintf("  [1] Non-disclosers in %d with MNAR:   %d\n", pre_year, nrow(never)))
  if (nrow(never) == 0) return(list(sample = never, ar = NULL))

  lk <- link_permno(never$gvkey, event_date)
  never <- merge(never, lk, by = "gvkey")
  cat(sprintf("  [2] Linked to permno:                %d\n", nrow(never)))

  lo <- event_date + ESTIMATION_WINDOW[1] - 10
  hi <- event_date + max(sapply(EVENT_WINDOWS, `[`, 2)) + 10
  ret <- pull_daily_returns(unique(never$permno), lo, hi)
  cat(sprintf("  [3] CRSP rows (%s..%s):               %d\n", lo, hi, nrow(ret)))

  car_out <- compute_cars(ret, event_date)
  cat(sprintf("  [4] Permnos with CARs:               %d\n", nrow(car_out$cars)))

  d <- merge(never, car_out$cars, by = "permno")
  # Require complete short window (primary test). Long windows may have gaps
  # (e.g. around COVID); don't let them kill the sample.
  expect_short <- diff(EVENT_WINDOWS$short) + 1L
  d <- d[n_short >= expect_short - 1]
  cat(sprintf("  [5] With complete [-1,+1] window:    %d\n", nrow(d)))

  d[, gsector := as.factor(gsector)]
  list(sample = d, ar = car_out$ar_panel)
}

augment_mar <- function(d, pre_year) {
  mar <- emis[year == pre_year & !is.na(sc1_disclosed_complete) &
              sc1_disclosed_complete > 0,
              .(gvkey = as.integer(gvkey),
                sc1_mar_pre = sc1_disclosed_complete)]
  d <- merge(d, mar, by = "gvkey", all.x = TRUE)
  d[, log_sc1_mar := log(sc1_mar_pre)]
  d
}

augment_trucost <- function(d, pre_year) {
  tru <- emis[year == pre_year & !is.na(sc1_estimated) & sc1_estimated > 0,
              .(gvkey = as.integer(gvkey),
                sc1_truc_pre = sc1_estimated)]
  d <- merge(d, tru, by = "gvkey", all.x = TRUE)
  d[, log_sc1_truc := log(sc1_truc_pre)]
  d
}


# =============================================================================
# RUN PIPELINE: MAIN + PLACEBO
# =============================================================================
main_res    <- run_event(EVENT_DATE_MAIN,    PRE_YEAR_MAIN,    "SEC proposal")
placebo_res <- run_event(EVENT_DATE_PLACEBO, PRE_YEAR_PLACEBO, "Placebo (no climate news)")

dbDisconnect(wrds)

sample_main    <- main_res$sample
sample_placebo <- placebo_res$sample
sample_main    <- augment_mar(augment_trucost(sample_main,    PRE_YEAR_MAIN),    PRE_YEAR_MAIN)
sample_placebo <- augment_mar(augment_trucost(sample_placebo, PRE_YEAR_PLACEBO), PRE_YEAR_PLACEBO)


# =============================================================================
# TABLE: MAIN-EVENT CARs
# Heteroskedasticity-robust SEs. Note the two vcov choices we do NOT use:
# fixest's default here is IID, which is untenable for CARs; and clustering on
# gsector would rest on 11 clusters, far below where cluster-robust inference
# is reliable. White SEs are the defensible option.
# =============================================================================
fit_mnar_all <- function(d) list(
  short  = feols(car_short  ~ log_sc1_adj + log_at | gsector, data = d, vcov = "hetero"),
  medium = feols(car_medium ~ log_sc1_adj + log_at | gsector, data = d, vcov = "hetero"),
  long   = feols(car_long   ~ log_sc1_adj + log_at | gsector, data = d, vcov = "hetero")
)
fits_main    <- fit_mnar_all(sample_main)
fits_placebo <- fit_mnar_all(sample_placebo)

# Hand-built rather than via etable(): this is the main-text exhibit, and it
# needs the event and the placebo side by side under spanning headers, which
# etable's `headers` argument does not render reliably.
ex <- function(fit, var) {
  ct <- coeftable(fit)[var, ]
  list(b  = as.numeric(ct["Estimate"]),
       se = as.numeric(ct["Std. Error"]),
       p  = as.numeric(ct["Pr(>|t|)"]),
       n  = fit$nobs,
       r2 = as.numeric(fitstat(fit, "r2", simplify = TRUE)))
}
star <- function(p) if (is.na(p)) "" else
  if (p < .01) "$^{***}$" else if (p < .05) "$^{**}$" else
  if (p < .10) "$^{*}$" else ""
f4 <- function(x) sprintf("%.4f", x)
fN <- function(x) formatC(x, format = "d", big.mark = ",")

COLS <- list(fits_main$short, fits_main$medium, fits_main$long,
             fits_placebo$short, fits_placebo$medium, fits_placebo$long)

coef_row <- function(label, var) {
  e <- lapply(COLS, ex, var = var)
  c(paste0(label, " & ",
           paste(sapply(e, function(z) paste0(f4(z$b), star(z$p))), collapse = " & "),
           " \\\\"),
    paste0(" & ",
           paste(sapply(e, function(z) paste0("(", f4(z$se), ")")), collapse = " & "),
           " \\\\"))
}
e1 <- lapply(COLS, ex, var = "log_sc1_adj")

tex_main <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Market reaction to the SEC climate-disclosure proposal.}",
  "Cumulative abnormal returns around the proposal of 21 March 2022, regressed",
  "on log recovered Scope~1 emissions measured in the pre-event year. The sample",
  "is firms that had not voluntarily disclosed as of that year, for whom the",
  "emissions measure is imputed rather than reported; abnormal returns come from",
  "a market model estimated over days $[-250,-30]$. The placebo columns repeat",
  "the identical construction on the same calendar date in 2019, a date with no",
  "climate-policy news, rebuilding the sample from firms not disclosing as of",
  "2018. Sector fixed effects throughout; heteroskedasticity-robust standard",
  "errors in parentheses.}",
  "\\label{tab:returns_event_study}",
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{3}{c}{SEC proposal (2022-03-21)} & \\multicolumn{3}{c}{Placebo (2019-03-21)} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  " & $[-1,+1]$ & $[-5,+5]$ & $[-30,+30]$ & $[-1,+1]$ & $[-5,+5]$ & $[-30,+30]$ \\\\",
  "\\midrule",
  coef_row("Log recovered emissions (MNAR)", "log_sc1_adj"),
  coef_row("Log assets", "log_at"),
  "\\midrule",
  paste0("Sector fixed effects & ", paste(rep("Yes", 6), collapse = " & "), " \\\\"),
  paste0("Observations & ",
         paste(sapply(e1, function(z) fN(z$n)), collapse = " & "), " \\\\"),
  paste0("$R^2$ & ",
         paste(sapply(e1, function(z) sprintf("%.3f", z$r2)), collapse = " & "), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Coefficients are per one-log-unit (roughly $2.7$-fold)",
  "increase in recovered emissions; CARs are in decimal returns, so $-0.0051$",
  "is $-0.51$ percentage points. A multi-date version of the placebo, comparing",
  "the event coefficient with the distribution over the twenty-first of every",
  "month in the sample, is reported in Table~\\ref{tab:returns_event_placebo_dist}.",
  "Significance: $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$.",
  "\\end{flushleft}", "\\end{table}")
writeLines(tex_main, file.path(outdir_tab, "returns_event_study.tex"))
cat("Saved: returns_event_study.tex (event + placebo, main text)\n")


# =============================================================================
# TABLE: CAR decile plot (main event, [-5,+5])
# =============================================================================
sample_main[, mnar_decile := as.integer(cut(log_sc1_adj,
  breaks = quantile(log_sc1_adj, probs = 0:10/10, na.rm = TRUE),
  include.lowest = TRUE))]

decile_summary <- sample_main[, .(
  car_medium_mean = mean(car_medium, na.rm = TRUE),
  car_medium_se   = sd(car_medium, na.rm = TRUE) / sqrt(.N),
  n               = .N
), by = mnar_decile][order(mnar_decile)]
print(decile_summary)

p <- ggplot(decile_summary, aes(x = mnar_decile, y = car_medium_mean)) +
  geom_col(fill = "#1B4F8A", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = car_medium_mean - 1.96 * car_medium_se,
        ymax = car_medium_mean + 1.96 * car_medium_se),
    width = 0.25, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "MNAR emissions decile (2021, firms not disclosing as of 2021)",
       y = "CAR [-5, +5] around 2022-03-21",
       title = "CAR around SEC climate proposal by MNAR-implied emissions",
       caption = "Source: CRSP daily, market-model AR. Firms not disclosing as of 2021.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(outdir_fig, "returns_event_study_car.pdf"),
       p, width = 8, height = 4.5)
cat("Saved: returns_event_study_car.pdf\n")


# NOTE: the single-date placebo no longer gets its own table -- it is columns
# (4)-(6) of returns_event_study.tex above. The multi-date placebo distribution,
# which supersedes it as evidence, is produced by 03a_event_study_placebos.R.


# =============================================================================
# MAR vs MNAR HORSE RACE (same firms)
# =============================================================================
fit_mar_vs_mnar <- function(d) {
  d <- d[!is.na(log_sc1_mar) & !is.na(log_sc1_adj)]
  list(
    mnar_short  = feols(car_short  ~ log_sc1_adj + log_at | gsector, data = d, vcov = "hetero"),
    mnar_medium = feols(car_medium ~ log_sc1_adj + log_at | gsector, data = d, vcov = "hetero"),
    mnar_long   = feols(car_long   ~ log_sc1_adj + log_at | gsector, data = d, vcov = "hetero"),
    mar_short   = feols(car_short  ~ log_sc1_mar + log_at | gsector, data = d, vcov = "hetero"),
    mar_medium  = feols(car_medium ~ log_sc1_mar + log_at | gsector, data = d, vcov = "hetero"),
    mar_long    = feols(car_long   ~ log_sc1_mar + log_at | gsector, data = d, vcov = "hetero")
  )
}
mv_main <- fit_mar_vs_mnar(sample_main)
mv_plac <- fit_mar_vs_mnar(sample_placebo)

etable(mv_main$mnar_short, mv_main$mar_short,
       mv_main$mnar_medium, mv_main$mar_medium,
       mv_main$mnar_long,   mv_main$mar_long,
       headers = c("MNAR [-1,+1]", "MAR [-1,+1]",
                   "MNAR [-5,+5]", "MAR [-5,+5]",
                   "MNAR [-30,+30]", "MAR [-30,+30]"),
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
       tex = TRUE,
       file = file.path(outdir_tab, "returns_event_study_mar_vs_mnar_main.tex"),
       replace = TRUE,
       title = paste0(
         "MAR vs MNAR event-study CARs around the SEC climate disclosure ",
         "proposal (2022-03-21). Same sample; MNAR = $\\log(\\text{sc1\\_adj})$; ",
         "MAR = $\\log(\\text{sc1\\_disclosed\\_complete})$."),
       label = "tab:returns_event_study_mar_vs_mnar_main")

etable(mv_plac$mnar_short, mv_plac$mar_short,
       mv_plac$mnar_medium, mv_plac$mar_medium,
       mv_plac$mnar_long,   mv_plac$mar_long,
       headers = c("MNAR [-1,+1]", "MAR [-1,+1]",
                   "MNAR [-5,+5]", "MAR [-5,+5]",
                   "MNAR [-30,+30]", "MAR [-30,+30]"),
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
       tex = TRUE,
       file = file.path(outdir_tab, "returns_event_study_mar_vs_mnar_plac.tex"),
       replace = TRUE,
       title = "MAR vs MNAR placebo event study (2019-03-21, no climate-policy news).",
       label = "tab:returns_event_study_mar_vs_mnar_plac")
cat("Saved: returns_event_study_mar_vs_mnar_main.tex\n")
cat("Saved: returns_event_study_mar_vs_mnar_plac.tex\n")


# =============================================================================
# TRUCOST vs MNAR HORSE RACE (same firms)
# =============================================================================
fit_truc_vs_mnar <- function(d, window) {
  car_col <- paste0("car_", window)
  d <- d[!is.na(log_sc1_adj) & !is.na(log_sc1_truc) & !is.na(get(car_col))]
  f_mnar <- as.formula(sprintf("%s ~ log_sc1_adj + log_at | gsector", car_col))
  f_truc <- as.formula(sprintf("%s ~ log_sc1_truc + log_at | gsector", car_col))
  f_both <- as.formula(sprintf("%s ~ log_sc1_adj + log_sc1_truc + log_at | gsector", car_col))
  list(
    mnar = feols(f_mnar, data = d, vcov = "hetero"),
    truc = feols(f_truc, data = d, vcov = "hetero"),
    both = feols(f_both, data = d, vcov = "hetero")
  )
}
tv_main <- setNames(lapply(c("short","medium","long"),
                           function(w) fit_truc_vs_mnar(sample_main, w)),
                    c("short","medium","long"))
tv_plac <- setNames(lapply(c("short","medium","long"),
                           function(w) fit_truc_vs_mnar(sample_placebo, w)),
                    c("short","medium","long"))

cat(sprintf("\nCorr log_sc1_adj vs log_sc1_truc (main):    %.3f\n",
            cor(sample_main$log_sc1_adj, sample_main$log_sc1_truc, use="complete.obs")))
cat(sprintf("Corr log_sc1_adj vs log_sc1_truc (placebo): %.3f\n",
            cor(sample_placebo$log_sc1_adj, sample_placebo$log_sc1_truc, use="complete.obs")))

write_truc_table <- function(fits, outfile, title_str, label_str) {
  etable(
    fits$short$mnar,  fits$short$truc,  fits$short$both,
    fits$medium$mnar, fits$medium$truc, fits$medium$both,
    fits$long$mnar,   fits$long$truc,   fits$long$both,
    headers = c("MNAR [-1,+1]",   "Truc [-1,+1]",   "Both [-1,+1]",
                "MNAR [-5,+5]",   "Truc [-5,+5]",   "Both [-5,+5]",
                "MNAR [-30,+30]", "Truc [-30,+30]", "Both [-30,+30]"),
    signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    tex = TRUE, file = file.path(outdir_tab, outfile), replace = TRUE,
    title = title_str, label = label_str)
}
write_truc_table(tv_main,
  "returns_event_study_trucost_vs_mnar_main.tex",
  paste0("Horse race: Trucost vs MNAR around the SEC climate disclosure proposal ",
         "(2022-03-21). MNAR = $\\log(\\text{sc1\\_adj})$; ",
         "Trucost = $\\log(\\text{sc1\\_estimated})$."),
  "tab:returns_event_study_trucost_vs_mnar_main")
write_truc_table(tv_plac,
  "returns_event_study_trucost_vs_mnar_plac.tex",
  "Horse race: Trucost vs MNAR on placebo date (2019-03-21).",
  "tab:returns_event_study_trucost_vs_mnar_plac")
cat("Saved: returns_event_study_trucost_vs_mnar_main.tex\n")
cat("Saved: returns_event_study_trucost_vs_mnar_plac.tex\n")


# =============================================================================
# CONSOLE SUMMARY
# =============================================================================
cat("\n===== EVENT-STUDY SUMMARY =====\n")
for (lbl in c("MAIN", "PLACEBO")) {
  d  <- if (lbl == "MAIN") sample_main      else sample_placebo
  dt <- if (lbl == "MAIN") EVENT_DATE_MAIN  else EVENT_DATE_PLACEBO
  cat(sprintf("\n%s (%s), N=%d\n", lbl, dt, nrow(d)))
  for (w in names(EVENT_WINDOWS)) {
    m <- mean(d[[paste0("car_", w)]], na.rm = TRUE)
    s <- sd(d[[paste0("car_", w)]],   na.rm = TRUE) / sqrt(nrow(d))
    cat(sprintf("  Mean CAR %-10s = %+.4f (SE %.4f, t=%.2f)\n",
                paste0("[", paste(EVENT_WINDOWS[[w]], collapse = ","), "]"),
                m, s, m / s))
  }
}

cat("\n=== 03_event_study.R COMPLETE ===\n")
