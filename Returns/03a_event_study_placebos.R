# =============================================================================
# 03a_event_study_placebos.R
# Multi-date placebo / randomization inference for the SEC event study.
#
# WHY THIS EXISTS
# 03_event_study.R compares the SEC proposal (2022-03-21) against a SINGLE
# placebo date (2019-03-21). One placebo cannot say whether the event-date
# coefficient is unusual -- it only says one other date was quiet. Here we
# re-estimate the IDENTICAL regression on many non-event dates chosen by a
# pre-specified rule and ask where the event-date coefficient falls in the
# resulting distribution. The rule picks the dates, not the researcher, which
# is what makes this immune to the date-shopping objection.
#
# INFERENCE
# The reported p-value is an empirical (randomization) p-value: the share of
# placebo coefficients at least as negative as the event coefficient. Note the
# resolution limit -- with K placebos the smallest achievable p is 1/(K+1), so
# the "annual" scheme (10 placebos) cannot go below p = 0.09. Use "monthly"
# if you need to reject at 5%.
#
# SCHEMES
#   "annual"  — 21 March of every other year in the sample. Holds calendar
#               seasonality exactly fixed; only ~10 dates.
#   "monthly" — the 21st of every month. ~130 dates, finer resolution, but
#               does not hold seasonality fixed and adjacent [-30,+30] windows
#               overlap (short/medium windows do not).
#
# EFFICIENCY
# Daily CRSP returns are pulled ONCE for the union of permnos over the full
# span, then CARs are recomputed per candidate date in memory. Pulling per
# date would mean one WRDS round-trip per placebo.
#
# Outputs:
#   returns_event_placebo_dist.tex — event vs placebo distribution, per window
#   returns_event_placebo_dates.tex — per-date coefficients (annual scheme)
#   returns_event_placebo_dist.pdf  — placebo histogram with event marked
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
# CONFIG  (kept identical to 03_event_study.R where it matters)
# =============================================================================
EVENT_DATE  <- as.Date("2022-03-21")   # SEC proposal
EVENT_PRE_Y <- 2021

PLACEBO_SCHEME <- "annual"             # "annual" | "monthly"
PLACEBO_DAY    <- 21L                  # day-of-month used by both schemes

ESTIMATION_WINDOW <- c(-250L, -30L)
EVENT_WINDOWS <- list(
  short  = c(-1L,  1L),
  medium = c(-5L,  5L),
  long   = c(-30L, 30L)
)
MIN_EST_OBS   <- 100L
SAMPLE_GROUPS <- "none_estimated"      # see 03_event_study.R

# Dates we do NOT treat as clean placebos. Reported, but flagged: a reader
# should not have to discover these confounds themselves.
CONFOUNDED <- list(
  "2020" = "COVID-19 crash; extreme cross-sectional dispersion",
  "2021" = "US re-entry to Paris (Feb 2021); climate policy actively in the news"
)

WRDS_USER <- "gcrippa4"


# =============================================================================
# LOAD PANEL, OPEN WRDS
# =============================================================================
cat("Loading MNAR panel...\n")
emis <- as.data.table(
  read_parquet("/Users/giulia/Documents/GitHub/corporate-omissions/data/processed/df_lm_avg_baseline_wti_gind_exchg.parquet")
)
setnames(emis, "fyear", "year")

YEARS_AVAIL <- sort(unique(emis$year))

# Candidate dates: an event needs a pre-year present in the panel.
build_candidate_dates <- function(scheme) {
  yrs <- YEARS_AVAIL[YEARS_AVAIL + 1L <= max(YEARS_AVAIL)] + 1L   # event year = pre_year + 1
  if (scheme == "annual") {
    d <- as.Date(sprintf("%d-03-%02d", yrs, PLACEBO_DAY))
  } else if (scheme == "monthly") {
    d <- as.Date(unlist(lapply(yrs, function(y)
      sprintf("%d-%02d-%02d", y, 1:12, PLACEBO_DAY))), origin = "1970-01-01")
  } else stop("PLACEBO_SCHEME must be 'annual' or 'monthly'")
  d <- d[!is.na(d)]
  sort(unique(c(d, EVENT_DATE)))
}
cand_dates <- build_candidate_dates(PLACEBO_SCHEME)
cat(sprintf("Candidate dates (%s scheme): %d\n", PLACEBO_SCHEME, length(cand_dates)))

cat("Connecting to WRDS...\n")
wrds <- dbConnect(
  Postgres(),
  host    = "wrds-pgdata.wharton.upenn.edu",
  port    = 9737,
  dbname  = "wrds",
  sslmode = "require",
  user    = WRDS_USER
)


# =============================================================================
# HELPERS  (mirrors 03_event_study.R; kept in sync by hand)
# =============================================================================
build_nondisclosing_sample <- function(pre_year) {
  d <- emis[year == pre_year]
  if (nrow(d) == 0) return(NULL)
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
  if (nrow(d) == 0) return(NULL)
  d[, log_sc1_adj := log(sc1_adj_pre)]
  d[, log_at      := log(pmax(at_pre, 1e-6))]
  d[, gvkey       := as.integer(gvkey)]
  d
}

# Pull the CCM link table ONCE for every gvkey we might need, then filter by
# date in memory rather than re-querying per candidate date.
all_gvkeys <- sort(unique(as.integer(emis$gvkey)))
cat(sprintf("Pulling CCM link table for %d gvkeys...\n", length(all_gvkeys)))
gv_sql <- paste(sprintf("'%06d'", all_gvkeys), collapse = ",")
CCM <- as.data.table(dbGetQuery(wrds, sprintf("
  SELECT gvkey, lpermno AS permno, linktype, linkprim, linkdt, linkenddt
  FROM crsp.ccmxpf_linktable
  WHERE gvkey IN (%s)
    AND linktype IN ('LU','LC','LN','LS','LX','LD')
    AND linkprim IN ('P','C')
", gv_sql)))
CCM[is.na(linkenddt), linkenddt := as.Date("2099-12-31")]
CCM[, gvkey := as.integer(gvkey)]

link_permno_cached <- function(gvkeys, event_date) {
  x <- CCM[gvkey %in% gvkeys & event_date >= linkdt & event_date <= linkenddt]
  unique(x[, .(gvkey, permno)])
}

# One CRSP pull covering every candidate window.
span_lo <- min(cand_dates) + ESTIMATION_WINDOW[1] - 10
span_hi <- max(cand_dates) + max(sapply(EVENT_WINDOWS, `[`, 2)) + 10
all_permnos <- sort(unique(CCM$permno))
cat(sprintf("Pulling CRSP daily %s..%s for %d permnos (one query)...\n",
            span_lo, span_hi, length(all_permnos)))
pn_sql <- paste(sprintf("'%d'", all_permnos), collapse = ",")
DSF <- as.data.table(dbGetQuery(wrds, sprintf("
  SELECT permno, date, ret FROM crsp.dsf
  WHERE permno IN (%s) AND date BETWEEN '%s' AND '%s'
", pn_sql, span_lo, span_hi)))
DSI <- as.data.table(dbGetQuery(wrds, sprintf("
  SELECT date, vwretd AS mktret FROM crsp.dsi
  WHERE date BETWEEN '%s' AND '%s'
", span_lo, span_hi)))
RET <- merge(DSF, DSI, by = "date")
RET[, `:=`(date = as.Date(date), ret = as.numeric(ret), mktret = as.numeric(mktret))]
RET <- RET[!is.na(ret) & !is.na(mktret)]
setkey(RET, permno, date)
dbDisconnect(wrds)
cat(sprintf("CRSP panel: %s rows\n", format(nrow(RET), big.mark = ",")))

TRADING_DAYS <- sort(unique(RET$date))

# CARs for one date, computed from the cached panel.
cars_for_date <- function(permnos, event_date) {
  lo_d <- event_date + ESTIMATION_WINDOW[1] - 10
  hi_d <- event_date + max(sapply(EVENT_WINDOWS, `[`, 2)) + 10
  r <- RET[permno %in% permnos & date %between% list(lo_d, hi_d)]
  if (nrow(r) == 0) return(NULL)

  td <- TRADING_DAYS[TRADING_DAYS %between% list(lo_d, hi_d)]
  t0 <- max(td[td <= event_date])
  if (!length(t0) || is.na(t0)) return(NULL)
  idx <- data.table(date = td, t_rel = seq_along(td) - which(td == t0))
  r <- merge(r, idx, by = "date")

  est_mm <- function(g) {
    d <- g[t_rel %between% ESTIMATION_WINDOW]
    if (nrow(d) < MIN_EST_OBS) return(NULL)
    fit <- lm(ret ~ mktret, data = d)
    co <- coef(fit)
    data.table(alpha = co[1], beta = co[2])
  }
  mm <- r[, est_mm(.SD), by = permno]
  if (nrow(mm) == 0) return(NULL)
  r <- merge(r, mm, by = "permno")
  r[, ar := ret - (alpha + beta * mktret)]

  out <- lapply(names(EVENT_WINDOWS), function(w) {
    lo <- EVENT_WINDOWS[[w]][1]; hi <- EVENT_WINDOWS[[w]][2]
    d <- r[t_rel %between% c(lo, hi),
           .(car = sum(ar, na.rm = TRUE), n_days = .N), by = permno]
    setnames(d, c("car", "n_days"), c(paste0("car_", w), paste0("n_", w)))
    d
  })
  Reduce(function(a, b) merge(a, b, by = "permno", all = TRUE), out)
}

# Full pipeline for one candidate date -> one row of coefficients.
run_one_date <- function(event_date) {
  pre_year <- as.integer(format(event_date, "%Y")) - 1L
  samp <- build_nondisclosing_sample(pre_year)
  if (is.null(samp)) return(NULL)

  lk <- link_permno_cached(samp$gvkey, event_date)
  if (nrow(lk) == 0) return(NULL)
  samp <- merge(samp, lk, by = "gvkey")

  cars <- cars_for_date(unique(samp$permno), event_date)
  if (is.null(cars)) return(NULL)
  d <- merge(samp, cars, by = "permno")
  d <- d[n_short >= (diff(EVENT_WINDOWS$short) + 1L) - 1L]
  if (nrow(d) < 50 || uniqueN(d$gsector) < 2) return(NULL)
  d[, gsector := as.factor(gsector)]

  res <- data.table(event_date = event_date, pre_year = pre_year, n = nrow(d))
  for (w in names(EVENT_WINDOWS)) {
    cc <- paste0("car_", w)
    dd <- d[!is.na(get(cc))]
    fit <- tryCatch(
      feols(as.formula(sprintf("%s ~ log_sc1_adj + log_at | gsector", cc)),
            data = dd, vcov = "hetero"),
      error = function(e) NULL)
    if (is.null(fit) || !"log_sc1_adj" %in% rownames(coeftable(fit))) {
      res[, (paste0(c("coef_","se_","p_"), w)) := NA_real_]
    } else {
      ct <- coeftable(fit)["log_sc1_adj", ]
      res[, (paste0("coef_", w)) := as.numeric(ct["Estimate"])]
      res[, (paste0("se_",   w)) := as.numeric(ct["Std. Error"])]
      res[, (paste0("p_",    w)) := as.numeric(ct["Pr(>|t|)"])]
    }
  }
  res
}


# =============================================================================
# RUN ALL DATES
# =============================================================================
cat("\n=== Estimating across candidate dates ===\n")
rows <- vector("list", length(cand_dates))
for (i in seq_along(cand_dates)) {
  dt <- cand_dates[i]
  r  <- run_one_date(dt)
  rows[[i]] <- r
  cat(sprintf("  %s  %s\n", format(dt),
              if (is.null(r)) "skipped (insufficient data)"
              else sprintf("N=%5d  b[-1,+1]=%+.4f  b[-5,+5]=%+.4f",
                           r$n, r$coef_short, r$coef_medium)))
}
res <- rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
res[, is_event := event_date == EVENT_DATE]
res[, yr := format(event_date, "%Y")]
res[, confounded := yr %in% names(CONFOUNDED) & !is_event]

stopifnot(sum(res$is_event) == 1)


# =============================================================================
# RANDOMIZATION INFERENCE
# =============================================================================
# One-sided: the hypothesis is that the announcement moved high-emission firms
# DOWN, so we ask how often a placebo date is at least as negative.
summarise_window <- function(w) {
  cf <- paste0("coef_", w)
  ev <- res[is_event == TRUE][[cf]]
  pl <- res[is_event == FALSE & !is.na(get(cf))][[cf]]
  pl_clean <- res[is_event == FALSE & confounded == FALSE & !is.na(get(cf))][[cf]]
  data.table(
    window     = w,
    event_coef = ev,
    k_placebo  = length(pl),
    pl_mean    = mean(pl),
    pl_sd      = sd(pl),
    pl_min     = min(pl),
    pl_max     = max(pl),
    n_below    = sum(pl <= ev),
    p_emp      = (sum(pl <= ev) + 1) / (length(pl) + 1),
    p_emp_clean = (sum(pl_clean <= ev) + 1) / (length(pl_clean) + 1),
    z_vs_pl    = (ev - mean(pl)) / sd(pl)
  )
}
ri <- rbindlist(lapply(names(EVENT_WINDOWS), summarise_window))

cat("\n=== RANDOMIZATION INFERENCE ===\n")
print(ri[, .(window, event_coef = round(event_coef, 4), k_placebo,
             pl_mean = round(pl_mean, 4), pl_sd = round(pl_sd, 4),
             n_below, p_emp = round(p_emp, 3),
             p_emp_clean = round(p_emp_clean, 3), z = round(z_vs_pl, 2))])
cat(sprintf("\nResolution floor with %d placebos: smallest achievable p = %.3f\n",
            ri$k_placebo[1], 1 / (ri$k_placebo[1] + 1)))
if (any(res$confounded))
  cat("Confounded dates (reported, excluded from p_emp_clean): ",
      paste(sort(unique(res[confounded == TRUE]$yr)), collapse = ", "), "\n")


# =============================================================================
# TABLE 1: randomization summary
# =============================================================================
fmt3 <- function(x) sprintf("%.4f", x)
wlab <- c(short = "$[-1,+1]$", medium = "$[-5,+5]$", long = "$[-30,+30]$")

tex1 <- c(
  "\\begin{table}[htbp]", "\\footnotesize\\centering",
  "\\caption{\\textbf{Multi-date placebo: randomization inference.}",
  "The event-date coefficient on log recovered emissions is compared with the",
  "distribution of coefficients from the identical regression estimated on",
  sprintf("non-event dates selected by a pre-specified rule (%s scheme: %s).",
          PLACEBO_SCHEME,
          if (PLACEBO_SCHEME == "annual") "21 March of every other sample year"
          else "the 21st of every month"),
  "$p_{\\text{emp}}$ is the share of placebo coefficients at least as negative as",
  "the event coefficient, $(\\#\\{b_{p} \\le b_{e}\\}+1)/(K+1)$.",
  "$p_{\\text{clean}}$ excludes dates flagged as confounded (see notes).",
  "Standard errors are heteroskedasticity-robust; sector fixed effects and a",
  "log-assets control throughout.}",
  "\\label{tab:returns_event_placebo_dist}",
  "\\begin{tabular}{lccccccc}", "\\toprule",
  " & Event & \\multicolumn{4}{c}{Placebo distribution} & \\multicolumn{2}{c}{Randomization} \\\\",
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-8}",
  "Window & coef. & mean & s.d. & min & max & $p_{\\text{emp}}$ & $p_{\\text{clean}}$ \\\\",
  "\\midrule")
for (w in names(EVENT_WINDOWS)) {
  r <- ri[window == w]
  tex1 <- c(tex1, paste0(
    wlab[[w]], " & ", fmt3(r$event_coef), " & ", fmt3(r$pl_mean), " & ",
    fmt3(r$pl_sd), " & ", fmt3(r$pl_min), " & ", fmt3(r$pl_max), " & ",
    sprintf("%.3f", r$p_emp), " & ", sprintf("%.3f", r$p_emp_clean), " \\\\"))
}
tex1 <- c(tex1, "\\midrule",
  paste0("Placebo dates ($K$) & \\multicolumn{7}{c}{", ri$k_placebo[1], "} \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{flushleft}\\footnotesize",
  "\\textit{Notes:} Sample for each date is firms not disclosing as of the",
  "preceding year, so both the sample and the emissions measure are rebuilt",
  paste0("date by date. With $K$ placebos the smallest attainable empirical ",
         "$p$-value is $1/(K+1) = ", sprintf("%.3f", 1/(ri$k_placebo[1]+1)), "$."),
  "Dates flagged as confounded and excluded from $p_{\\text{clean}}$:",
  paste0(paste(sprintf("%s (%s)", names(CONFOUNDED), unlist(CONFOUNDED)),
               collapse = "; "), "."),
  "\\end{flushleft}", "\\end{table}")
writeLines(tex1, file.path(outdir_tab, "returns_event_placebo_dist.tex"))
cat("Saved: returns_event_placebo_dist.tex\n")


# =============================================================================
# TABLE 2: per-date coefficients (only worth printing for the annual scheme)
# =============================================================================
if (PLACEBO_SCHEME == "annual") {
  tex2 <- c(
    "\\begin{table}[htbp]", "\\footnotesize\\centering",
    "\\caption{\\textbf{Per-date placebo coefficients.}",
    "Coefficient on log recovered emissions from the event-study regression,",
    "estimated on 21 March of each sample year. The SEC proposal date is in",
    "bold. Sector fixed effects, log-assets control,",
    "heteroskedasticity-robust standard errors.}",
    "\\label{tab:returns_event_placebo_dates}",
    "\\begin{tabular}{lrcccc}", "\\toprule",
    "Date & $N$ & $[-1,+1]$ & $[-5,+5]$ & $[-30,+30]$ & Note \\\\", "\\midrule")
  setorder(res, event_date)
  for (i in seq_len(nrow(res))) {
    r <- res[i]
    lab <- format(r$event_date)
    nb  <- if (r$is_event) "\\textbf{SEC proposal}"
           else if (r$confounded) CONFOUNDED[[r$yr]] else ""
    bold <- function(x) if (r$is_event) paste0("\\textbf{", x, "}") else x
    tex2 <- c(tex2, paste0(
      bold(lab), " & ", format(r$n, big.mark = ","), " & ",
      bold(fmt3(r$coef_short)), " & ", bold(fmt3(r$coef_medium)), " & ",
      bold(fmt3(r$coef_long)), " & \\footnotesize ", nb, " \\\\"))
  }
  tex2 <- c(tex2, "\\bottomrule", "\\end{tabular}",
    "\\begin{flushleft}\\footnotesize",
    "\\textit{Notes:} Each row rebuilds the sample from firms not disclosing as",
    "of the preceding year, so $N$ varies with panel coverage.",
    "\\end{flushleft}", "\\end{table}")
  writeLines(tex2, file.path(outdir_tab, "returns_event_placebo_dates.tex"))
  cat("Saved: returns_event_placebo_dates.tex\n")
}


# =============================================================================
# FIGURE: placebo distribution with the event coefficient marked
# =============================================================================
pl <- res[is_event == FALSE, .(event_date, yr, confounded, coef_medium)]
ev <- res[is_event == TRUE, coef_medium]

p <- ggplot(pl, aes(x = coef_medium)) +
  geom_histogram(aes(fill = confounded), bins = max(8, nrow(pl) %/% 3),
                 colour = "white", alpha = 0.85) +
  scale_fill_manual(values = c(`FALSE` = "#9EB9D4", `TRUE` = "#E0B080"),
                    labels = c("placebo", "placebo (confounded)"),
                    name = NULL) +
  geom_vline(xintercept = ev, colour = "#B2182B", linewidth = 0.9) +
  annotate("text", x = ev, y = Inf, label = "  SEC proposal", hjust = 0, vjust = 1.6,
           colour = "#B2182B", size = 3.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(x = "Coefficient on log recovered emissions, CAR [-5,+5]",
       y = "Placebo dates",
       title = "Event coefficient against the placebo distribution",
       caption = sprintf("%s scheme; %d placebo dates. Sample rebuilt per date.",
                         PLACEBO_SCHEME, nrow(pl))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(outdir_fig, "returns_event_placebo_dist.pdf"), p,
       width = 7.5, height = 4.2)
cat("Saved: returns_event_placebo_dist.pdf\n")

fwrite(res, file.path(outdir_tab, "returns_event_placebo_raw.csv"))
cat("Saved: returns_event_placebo_raw.csv (per-date raw output)\n")

cat("\n=== 03a_event_study_placebos.R COMPLETE ===\n")
