# ================================
# Paper Pack — CATE Summaries (v1, patched)
# Requires outputs from v7 toolkit
# ================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  ok_patch <- requireNamespace("patchwork", quietly = TRUE)
})

# -------- Paths / params --------
RESULTS_Q <- "nft_trait_v7_results_with_q.csv"
OUT_DIR   <- "nft_cf_cate_reports_v7"     # where v7 wrote per-combo files
PAPER_DIR <- file.path(OUT_DIR, "paper_assets")
dir.create(PAPER_DIR, showWarnings = FALSE, recursive = TRUE)

K_SET  <- 1:4
ROUND1 <- function(x) ifelse(is.finite(as.numeric(x)), round(as.numeric(x), 1), NA_real_)
ROUND2 <- function(x) ifelse(is.finite(as.numeric(x)), round(as.numeric(x), 2), NA_real_)

shorten_combo <- function(s) {
  s <- gsub("Background_", "BG:", s)
  s <- gsub("Earring_",   "Ear:", s)
  s <- gsub("Clothes_",   "Cl:",  s)
  s <- gsub("Eyes_",      "Eyes:", s)
  s <- gsub("Fur_",       "Fur:", s)
  s <- gsub("Hat_",       "Hat:", s)
  s <- gsub("Mouth_",     "Mouth:", s)
  s <- gsub("\\.", "·", s)
  s
}

safe_combo_name <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

# -------- Load flagship list (top-N per K from v7). Pick top-1 per K --------
flag_path <- file.path(OUT_DIR, "selected_flagship_treatments.csv")
if (!file.exists(flag_path)) stop("Missing ", flag_path, " — run the v7 CATE Reporting Toolkit first.")
flag <- data.table::fread(flag_path)

# restrict to K_SET and take the first-ranked (rank_k) per K
pick <- flag[k %in% K_SET][order(k, rank_k), .SD[1], by = k]
if (!nrow(pick)) stop("No flagship combos found for K in ", paste(K_SET, collapse = ","))

# -------- Load ATE% and kept% from results_with_q; build join stub --------
resq <- fread(RESULTS_Q)
if (!"traits_combination" %in% names(resq)) stop("traits_combination missing in results CSV")
resq[, join_stub := safe_combo_name(traits_combination)]

# be defensive about numeric columns coming in as character-with-commas
num_cols <- c("ATE_pct","CI_lower_pct","CI_upper_pct","kept_pct","n_kept","q_value_bh")
for (cc in intersect(num_cols, names(resq))) {
  if (is.character(resq[[cc]])) resq[[cc]] <- as.numeric(gsub(",", "", resq[[cc]]))
}

# -------- Helpers to read per-combo artifacts safely --------
ftry <- function(path) if (file.exists(path)) fread(path) else data.table()

path_gs   <- function(k, stub) file.path(OUT_DIR, sprintf("global_stats_k%d_%s.csv",    k, stub))
path_gt   <- function(k, stub) file.path(OUT_DIR, sprintf("gates_k%d_%s.csv",           k, stub))
path_seas <- function(k, stub) file.path(OUT_DIR, sprintf("seasonality_k%d_%s.csv",     k, stub))
path_cal  <- function(k, stub) file.path(OUT_DIR, sprintf("calibration_k%d_%s.csv",     k, stub))
path_blp  <- function(k, stub) file.path(OUT_DIR, sprintf("blp_k%d_%s.csv",             k, stub))
path_rob  <- function(k, stub) file.path(OUT_DIR, sprintf("robustness_k%d_%s.csv",      k, stub))
path_vi   <- function(k, stub) file.path(OUT_DIR, sprintf("varimp_top25_k%d_%s.csv",    k, stub))

# =========================
# TABLE A — CATE OVERVIEW
# =========================
build_row_A <- function(one_k) {
  k     <- one_k$k
  combo <- one_k$traits_combination
  stub  <- one_k$treatment_file_stub
  lab   <- shorten_combo(combo)
  
  gs <- ftry(path_gs(k, stub))
  gt <- ftry(path_gt(k, stub))
  
  # ATE% & kept% from results_with_q (join on sanitized stub + K)
  rq <- resq[k == one_k$k & join_stub == one_k$treatment_file_stub]
  # if multiple rows match (rare stub collision), prefer exact text match; else smallest q
  if (nrow(rq) > 1L) {
    exact <- rq[traits_combination == combo]
    rq <- if (nrow(exact)) exact else rq[order(q_value_bh, -abs(ATE_pct), -kept_pct)][1]
  }
  
  ATEp <- if (nrow(rq)) rq$ATE_pct[1] else NA_real_
  lo   <- if (nrow(rq)) rq$CI_lower_pct[1] else NA_real_
  hi   <- if (nrow(rq)) rq$CI_upper_pct[1] else NA_real_
  kept <- if (nrow(rq)) rq$kept_pct[1] else NA_real_
  # fallback: global_stats n is #kept
  nkept<- if (nrow(rq) && is.finite(rq$n_kept[1])) rq$n_kept[1] else if (nrow(gs)) gs$n[1] else NA_real_
  
  # GATES: top vs bottom decile
  if (nrow(gt)) {
    if (!is.numeric(gt$gates_dec)) gt[, gates_dec := as.numeric(as.character(gates_dec))]
    gtop <- gt[gates_dec == max(gates_dec, na.rm = TRUE)]
    gbot <- gt[gates_dec == min(gates_dec, na.rm = TRUE)]
    bot_mean <- if (nrow(gbot)) gbot$mean_pct[1] else NA_real_
    top_mean <- if (nrow(gtop)) gtop$mean_pct[1] else NA_real_
    lift <- if (all(is.finite(c(bot_mean, top_mean)))) top_mean - bot_mean else NA_real_
  } else {
    bot_mean <- top_mean <- lift <- NA_real_
  }
  
  data.table(
    K = k,
    Combination = lab,
    n_kept = ROUND1(nkept),
    kept_pct = ROUND1(kept),
    ATE_pct = ROUND1(ATEp),
    CI95_lo = ROUND1(lo),
    CI95_hi = ROUND1(hi),
    Median_pct = ROUND1(gs$median_cate_pct[1]),
    IQR_lo = ROUND1(gs$iqr_cate_pct_lo[1]),
    IQR_hi = ROUND1(gs$iqr_cate_pct_hi[1]),
    Share_pos_pct = ROUND1(100 * gs$share_positive[1]),
    GATES_bot_mean = ROUND1(bot_mean),
    GATES_top_mean = ROUND1(top_mean),
    GATES_lift = ROUND1(lift),
    Combo_full = combo,    # traceability
    Stub = stub
  )
}

tabA <- rbindlist(lapply(split(pick, pick$k), function(row) build_row_A(row[1,])), fill = TRUE)
setorder(tabA, K)
fwrite(tabA, file.path(PAPER_DIR, "paper_TableA_CATE_overview.csv"))

# =========================
# TABLE B — CHECKS & ROBUSTNESS
# =========================
build_row_B <- function(one_k) {
  k     <- one_k$k
  combo <- one_k$traits_combination
  stub  <- one_k$treatment_file_stub
  lab   <- shorten_combo(combo)
  
  cal <- ftry(path_cal(k, stub))
  blp <- ftry(path_blp(k, stub))
  rob <- ftry(path_rob(k, stub))
  
  data.table(
    K = k,
    Combination = lab,
    Cal_intercept = ROUND2(cal$intercept[1]),
    Cal_slope     = ROUND2(cal$slope[1]),
    BLP_month_beta= ROUND2(blp$coef[1]),
    BLP_month_se  = ROUND2(blp$se[1]),
    Base_ATE_all      = ROUND2(rob$baseline_est_all[1]),
    Robust_ATE_all    = ROUND2(rob$robust_est_all[1]),
    Delta_all         = ROUND2(abs(rob$robust_est_all[1] - rob$baseline_est_all[1])),
    Base_ATE_treated  = ROUND2(rob$baseline_est_treat[1]),
    Robust_ATE_treated= ROUND2(rob$robust_est_treat[1]),
    Delta_treated     = ROUND2(abs(rob$robust_est_treat[1] - rob$baseline_est_treat[1])),
    Combo_full = combo,
    Stub = stub
  )
}

tabB <- rbindlist(lapply(split(pick, pick$k), function(row) build_row_B(row[1,])), fill = TRUE)
setorder(tabB, K)
fwrite(tabB, file.path(PAPER_DIR, "paper_TableB_checks_robustness.csv"))

# =========================================
# FIGURE — 4-up GATES (one panel per K)
# =========================================
g_plot <- function(k, stub, title_txt){
  gt <- ftry(path_gt(k, stub))
  if (!nrow(gt)) return(ggplot() + theme_void() + ggtitle(paste0(title_txt, " (no GATES)")))
  if (!is.numeric(gt$gates_dec)) gt[, gates_dec := as.numeric(as.character(gates_dec))]
  ggplot(gt, aes(x=gates_dec, y=mean_pct)) +
    geom_point() +
    geom_errorbar(aes(ymin=lo_pct, ymax=hi_pct), width=0.2) +
    labs(x="CATE decile (low → high)", y="Mean CATE% (95% CI)", title=title_txt) +
    theme_minimal(base_size = 11)
}

plots <- lapply(seq_len(nrow(pick)), function(i){
  k <- pick$k[i]; stub <- pick$treatment_file_stub[i]; combo <- pick$traits_combination[i]
  ttl <- sprintf("K=%d: %s", k, shorten_combo(combo))
  g_plot(k, stub, ttl)
})

if (ok_patch) {
  fig <- patchwork::wrap_plots(plots, ncol = 2)
  ggsave(file.path(PAPER_DIR, "paper_Figure_GATES_4up.png"), fig, width = 10, height = 7, dpi = 300)
} else {
  warning("Package 'patchwork' not available; saving panels separately.")
  for (i in seq_along(plots)) {
    k <- pick$k[i]
    ggsave(file.path(PAPER_DIR, sprintf("paper_Figure_GATES_K%d.png", k)),
           plots[[i]], width = 5, height = 4, dpi = 300)
  }
}

# ======================================================
# OPTIONAL — Consensus variable-importance (across 4 combos)
#   - average rank (lower = more important)
#   - frequency in top-5
# ======================================================
vi_all <- rbindlist(lapply(seq_len(nrow(pick)), function(i){
  k <- pick$k[i]; stub <- pick$treatment_file_stub[i]; combo <- pick$traits_combination[i]
  vi <- ftry(path_vi(k, stub))
  if (!nrow(vi)) return(data.table())
  vi[, `:=`(K = k, Combination = shorten_combo(combo))]
  vi
}), fill = TRUE)

if (nrow(vi_all)) {
  # normalize feature names for readability (dots -> underscores)
  if ("feature" %in% names(vi_all)) vi_all[, feature := gsub("\\.", "_", feature)]
  
  vi_all[, top5 := rank <= 5]
  agg <- vi_all[, .(
    avg_rank = mean(rank, na.rm = TRUE),
    times_in_top5 = sum(top5, na.rm = TRUE),
    combos_covered = uniqueN(Combination)
  ), by = feature][order(avg_rank, -times_in_top5)]
  fwrite(agg, file.path(PAPER_DIR, "paper_VI_consensus.csv"))
  
  # tiny bar plot of top 10 features by avg_rank
  top10 <- head(agg[is.finite(avg_rank)], 10)
  if (nrow(top10)) {
    pvi <- ggplot(top10, aes(x = reorder(feature, avg_rank), y = avg_rank)) +
      geom_col() +
      coord_flip() +
      labs(x = "Feature", y = "Average rank (lower is better)",
           title = "Consensus variable importance (top 10 by avg rank)") +
      theme_minimal(base_size = 11)
    ggsave(file.path(PAPER_DIR, "paper_VI_consensus_top10.png"), pvi, width = 7, height = 5, dpi = 300)
  }
}

# -------- Done --------
message("✅ Wrote paper assets to: ", normalizePath(PAPER_DIR))
message(" - Table A: paper_TableA_CATE_overview.csv")
message(" - Table B: paper_TableB_checks_robustness.csv")
message(" - Figure:  paper_Figure_GATES_4up.png (or per-K panels if patchwork missing)")
if (nrow(vi_all)) message(" - Optional: paper_VI_consensus.csv (+ PNG)")
