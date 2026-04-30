# ================================
# Paper CATE Pack (v2)
# Uses outputs from v7 toolkit
# ================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  ok_patch <- requireNamespace("patchwork", quietly = TRUE)
  if (ok_patch) library(patchwork)
})

# -------- Paths / params --------
RESULTS_Q <- "nft_trait_v7_results_with_q.csv"
OUT_DIR   <- "nft_cf_cate_reports_v7"           # where v7 wrote per-combo files
PAPER_DIR <- file.path(OUT_DIR, "paper_assets")  # same place other paper files went
dir.create(PAPER_DIR, showWarnings = FALSE, recursive = TRUE)

K_SET          <- 1:4
TOP_N_PER_K    <- 5            # set to 1 if one combo per K
SHOW_DRILLDOWN <- TRUE         # set FALSE if only the main table + lift plot

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
ftry <- function(path) if (file.exists(path)) fread(path) else data.table()

# -------- Helper: file builders --------
path_gs   <- function(k, stub) file.path(OUT_DIR, sprintf("global_stats_k%d_%s.csv",    k, stub))
path_gt   <- function(k, stub) file.path(OUT_DIR, sprintf("gates_k%d_%s.csv",           k, stub))
path_cal  <- function(k, stub) file.path(OUT_DIR, sprintf("calibration_k%d_%s.csv",     k, stub))
path_vi   <- function(k, stub) file.path(OUT_DIR, sprintf("varimp_top25_k%d_%s.csv",    k, stub))
path_sub  <- function(k, stub, fam) file.path(OUT_DIR, sprintf("subgroup_%s_k%d_%s.csv", fam, k, stub))

# -------- Load flagship list (top-N per K from v7) --------
flag_path <- file.path(OUT_DIR, "selected_flagship_treatments.csv")
if (!file.exists(flag_path)) stop("Missing ", flag_path, " — run the v7 CATE Reporting Toolkit first.")
flag <- fread(flag_path)

# Filter to K_SET and take top-N per K (by rank_k)
pick <- flag[k %in% K_SET][order(k, rank_k)]
pick <- pick[, .SD[seq_len(min(.N, TOP_N_PER_K))], by = k]
if (!nrow(pick)) stop("No flagship combos found for K in ", paste(K_SET, collapse = ","))

# -------- Load ATE% and kept% from results_with_q; build join stub --------
resq <- fread(RESULTS_Q)
stopifnot("traits_combination" %in% names(resq))
resq[, join_stub := safe_combo_name(traits_combination)]
num_cols <- c("ATE_pct","CI_lower_pct","CI_upper_pct","kept_pct","n_kept","q_value_bh")
for (cc in intersect(num_cols, names(resq))) {
  if (is.character(resq[[cc]])) resq[[cc]] <- as.numeric(gsub(",", "", resq[[cc]]))
}

# -------- Main table builder (rows = combos) --------
build_row_combo <- function(one){
  k     <- one$k
  combo <- one$traits_combination
  stub  <- one$treatment_file_stub
  lab   <- shorten_combo(combo)
  
  # ATE and kept
  rq <- resq[k == k & join_stub == stub]
  if (nrow(rq) > 1L) {
    exact <- rq[traits_combination == combo]
    rq <- if (nrow(exact)) exact else rq[order(q_value_bh, -abs(ATE_pct), -kept_pct)][1]
  }
  ATEp <- if (nrow(rq)) rq$ATE_pct[1] else NA_real_
  lo   <- if (nrow(rq)) rq$CI_lower_pct[1] else NA_real_
  hi   <- if (nrow(rq)) rq$CI_upper_pct[1] else NA_real_
  kept <- if (nrow(rq)) rq$kept_pct[1] else NA_real_
  nkept<- if (nrow(rq) && is.finite(rq$n_kept[1])) rq$n_kept[1] else NA_real_
  qval <- if (nrow(rq)) rq$q_value_bh[1] else NA_real_
  
  # Global stats (median/IQR/share+)
  gs <- ftry(path_gs(k, stub))
  med <- if (nrow(gs)) gs$median_cate_pct[1] else NA_real_
  iqr_lo <- if (nrow(gs)) gs$iqr_cate_pct_lo[1] else NA_real_
  iqr_hi <- if (nrow(gs)) gs$iqr_cate_pct_hi[1] else NA_real_
  share_pos <- if (nrow(gs)) 100*gs$share_positive[1] else NA_real_
  
  # GATES top-bottom & lift
  gt <- ftry(path_gt(k, stub))
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
  
  # Calibration
  cal <- ftry(path_cal(k, stub))
  cal_slope <- if (nrow(cal)) cal$slope[1] else NA_real_
  
  # Tiny VI summary (top-3 by rank)
  vi <- ftry(path_vi(k, stub))
  vi3 <- if (nrow(vi)) paste(head(gsub("\\.", "_", vi$feature[order(vi$rank)]), 3), collapse = ", ") else NA_character_
  
  data.table(
    K = k,
    Combination = lab,
    ATE_pct = ROUND1(ATEp),
    CI95_lo = ROUND1(lo),
    CI95_hi = ROUND1(hi),
    n_kept = ROUND1(nkept),
    kept_pct = ROUND1(kept),
    GATES_bot_mean = ROUND1(bot_mean),
    GATES_top_mean = ROUND1(top_mean),
    GATES_lift = ROUND1(lift),
    Median_pct = ROUND1(med),
    IQR_lo = ROUND1(iqr_lo),
    IQR_hi = ROUND1(iqr_hi),
    Share_pos_pct = ROUND1(share_pos),
    Cal_slope = ROUND2(cal_slope),
    q_value_bh = ifelse(is.finite(qval), round(qval, 3), NA_real_),
    Combo_full = combo,
    Stub = stub
  )
}

tab_combo <- rbindlist(lapply(seq_len(nrow(pick)), function(i) build_row_combo(pick[i])), fill = TRUE)
tab_combo[, abs_ATE := abs(ATE_pct)]
setorder(tab_combo, K, -GATES_lift, -abs_ATE)
tab_combo[, abs_ATE := NULL]   # optional cleanup
tab_combo <- tab_combo[order(K, -GATES_lift, -abs(ATE_pct))]

fwrite(tab_combo, file.path(PAPER_DIR, "paper_TableC_CATE_main.csv"))

# -------- Figure A: Lift-by-Combo (bars), grouped by K --------
# Only plot rows with finite lift
lift_dat <- tab_combo[is.finite(GATES_lift)]
if (nrow(lift_dat)) {
  lift_dat[, row_lab := paste0("K", K, " · ", Combination)]
  p_lift <- ggplot(lift_dat, aes(x = reorder(row_lab, GATES_lift), y = GATES_lift, fill = factor(K))) +
    geom_col(width = 0.7, show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~ K, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Heterogeneity lift (Top decile − Bottom decile)  [pp]",
         title = "GATES heterogeneity (Top−Bottom) by combo",
         subtitle = "Rows = selected combinations (Top-N per K). Higher lift = stronger within-sample heterogeneity.") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(PAPER_DIR, "paper_FigureA_CATE_Lift_byCombo.png"), p_lift, width = 10, height = 7, dpi = 300)
}

# ================================
# Figure B: Drill-down subgroup CATEs
#   For the #1 combo per K, plot the two richest non-treatment families
# ================================
if (SHOW_DRILLDOWN) {
  # pick rank 1 per K
  top1 <- pick[order(k, rank_k)][, .SD[1], by = k]
  make_drill_plots <- function(one){
    k <- one$k; stub <- one$treatment_file_stub; lab <- shorten_combo(one$traits_combination)
    
    # discover which subgroup_* files exist for this combo
    sub_files <- list.files(OUT_DIR, pattern = sprintf("^subgroup_.*_k%d_%s\\.csv$", k, stub), full.names = TRUE)
    if (!length(sub_files)) return(NULL)
    
    # read all, rank families by count, keep top 2
    subs <- rbindlist(lapply(sub_files, ftry), fill = TRUE)
    if (!nrow(subs)) return(NULL)
    fam_rank <- subs[, .(N = sum(is.finite(mean_pct))), by = family][order(-N)]
    fam_use <- head(fam_rank$family, 2)
    
    panels <- lapply(fam_use, function(tp){
      dt <- subs[family == tp & is.finite(mean_pct)]
      if (!nrow(dt)) return(NULL)
      dt[, level := factor(level, levels = level[order(mean_pct, decreasing = TRUE)])]
      ggplot(dt, aes(y = level, x = mean_pct)) +
        geom_point() +
        geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct), height = 0.2) +
        labs(y = paste0(tp, " levels"),
             x = "Mean CATE% (95% CI)",
             title = sprintf("K=%d  %s — subgroup CATEs in %s", k, lab, tp)) +
        theme_minimal(base_size = 11)
    })
    panels <- Filter(Negate(is.null), panels)
    if (!length(panels)) return(NULL)
    
    if (ok_patch && length(panels) > 1) {
      fig <- wrap_plots(panels, ncol = 2)
    } else {
      fig <- panels[[1]]
    }
    out <- file.path(PAPER_DIR, sprintf("paper_FigureB_Drilldown_K%d_%s.png", k, stub))
    ggsave(out, fig, width = 10, height = 6, dpi = 300)
    out
  }
  
  drill_paths <- lapply(seq_len(nrow(top1)), function(i) make_drill_plots(top1[i]))
}

# -------- Console summary --------
message("✅ Wrote CATE assets to: ", normalizePath(PAPER_DIR))
message(" - TableC (rows = combos): paper_TableC_CATE_main.csv")
if (nrow(lift_dat)) message(" - FigureA: paper_FigureA_CATE_Lift_byCombo.png")
if (SHOW_DRILLDOWN) message(" - FigureB: per-K drilldowns saved (if subgroup files exist)")
