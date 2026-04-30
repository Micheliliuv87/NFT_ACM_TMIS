# =========================================================
# CATE Reporting Toolkit (v7) — for BAYC Price 1 Sale CF
# Outputs CSVs + PNGs for Methods/Results
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grf)
  library(caret)
  library(arules)   # not required here but harmless
})

set.seed(87)

# ------------------ CONFIG ------------------
RESULTS_Q  <- "nft_trait_v7_results_with_q.csv"
CATE_DIR   <- "nft_cf_diag_v7"             # from your CF run
OUT_DIR    <- "nft_cf_cate_reports_v7"     # new outputs
dir.create(OUT_DIR, showWarnings = FALSE)

# selection for "flagship" deep dives
K_SET          <- 1:4
TOP_N_PER_K    <- 5
KEPT_PCT_MIN   <- 60
REFIT_EACH     <- TRUE   # TRUE => re-fit cf to compute calibration/BLP/VI
RNG_SEED       <- 87
N_BOOT         <- 500     # bootstrap for CI of group means on % scale
ROBUST_MIN_NODE<- 50
ROBUST_TREES   <- 1000

# ------------------ UTILS ------------------
fmt_num  <- function(x) format(x, big.mark = ",", scientific = FALSE)
trait_key <- function(x) sub("_.*$", "", x)
safe_combo_name <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

rename_dummy_cols_clean <- function(X) {
  nm <- colnames(X)
  nm <- gsub("\\.", "_", nm)
  nm <- gsub("__+", "_", nm)
  nm <- gsub("^_+|_+$", "", nm)
  colnames(X) <- make.names(nm, unique = TRUE, allow_ = TRUE)
  X
}

# group mean + bootstrap CI on *percent* scale
boot_mean_ci_pct <- function(v, nboot=N_BOOT, seed=RNG_SEED) {
  set.seed(seed)
  m <- mean(v)
  if (length(v) < 2) return(list(mean=m, lo=NA_real_, hi=NA_real_))
  bs <- replicate(nboot, mean(sample(v, length(v), replace=TRUE)))
  lo <- quantile(bs, 0.025, names=FALSE)
  hi <- quantile(bs, 0.975, names=FALSE)
  list(mean=m, lo=lo, hi=hi)
}

# -------- safe extractors for GRF outputs (version-agnostic) ----------
safe_get_named <- function(x, nm) {
  if (is.null(x)) return(NA_real_)
  # Works for lists and named vectors
  out <- tryCatch({
    if (!is.null(names(x)) && nm %in% names(x)) x[[nm]] else NA_real_
  }, error=function(e) NA_real_)
  if (is.null(out) || (is.atomic(out) && length(out) == 0)) return(NA_real_)
  as.numeric(out)
}

safe_extract_calibration <- function(tc) {
  data.table(
    intercept    = safe_get_named(tc, "calibration_intercept"),
    intercept_se = safe_get_named(tc, "calibration_intercept_se"),
    slope        = safe_get_named(tc, "calibration_slope"),
    slope_se     = safe_get_named(tc, "calibration_slope_se"),
    p_intercept  = safe_get_named(tc, "calibration_intercept_pval"),
    p_slope      = safe_get_named(tc, "calibration_slope_pval")
  )
}

safe_extract_blp <- function(blp_obj, var_name="month_sale", default_index=2L) {
  if (is.null(blp_obj)) return(NULL)
  
  # Coefficients
  if (is.list(blp_obj) && !is.null(blp_obj$coefficients)) {
    coefs <- as.numeric(blp_obj$coefficients)
    cn    <- names(blp_obj$coefficients)
  } else {
    coefs <- as.numeric(blp_obj)   # atomic vector fallback
    cn    <- names(blp_obj)
  }
  
  # Determine which coefficient corresponds to var_name (or fallback index)
  idx <- if (!is.null(cn) && var_name %in% cn) which(cn == var_name) else default_index
  if (is.na(idx) || length(coefs) < idx) return(NULL)
  
  # SEs (if available)
  if (is.list(blp_obj) && !is.null(blp_obj$variance.estimates)) {
    se_vec <- sqrt(diag(blp_obj$variance.estimates))
    se <- if (length(se_vec) >= idx) as.numeric(se_vec[idx]) else NA_real_
  } else {
    se <- NA_real_
  }
  
  list(coef = as.numeric(coefs[idx]), se = se)
}

# ------------------ REBUILD df / X (to join covariates & optionally refit) ------------------
build_df_X <- function() {
  df1 <- read.csv("df_table1.csv")
  df2 <- read.csv("df_table3.csv")
  # Use base::merge to avoid data.table coercion
  df  <- base::merge(df1[,c("token_id","time_1_sale","price_1_sale")], df2, by="token_id", all.x=TRUE)
  rm(df1, df2)
  
  traits <- c("Background","Clothes","Earring","Eyes","Fur","Hat","Mouth")
  
  # in-place recode
  rare_thresh <- 0.01
  for (tr in traits) {
    df[[tr]][is.na(df[[tr]]) | df[[tr]]==""] <- "No"
    ft <- prop.table(table(df[[tr]]))
    rare_levels <- setdiff(names(ft)[ft < rare_thresh], "No")
    df[[tr]] <- ifelse(df[[tr]] %in% rare_levels, "Other", as.character(df[[tr]]))
    df[[tr]] <- factor(df[[tr]])
  }
  
  # outcome + time
  df$row_id   <- seq_len(nrow(df))
  if ("token_id" %in% names(df)) df$token_id <- NULL
  if ("rarity.rank" %in% names(df)) df$rarity.rank <- NULL
  
  # winsorize price
  p05 <- as.numeric(quantile(df$price_1_sale, 0.05, na.rm=TRUE))
  p95 <- as.numeric(quantile(df$price_1_sale, 0.95, na.rm=TRUE))
  df$price_1_sale <- pmin(pmax(df$price_1_sale, p05), p95)
  df$log_price <- log1p(df$price_1_sale)
  
  df$month_sale <- as.numeric(format(as.POSIXct(df$time_1_sale, origin="1970-01-01"), "%m"))
  
  Y <- df$log_price
  
  # one-hot for traits
  dummies <- caret::dummyVars(~ ., data = df[, traits])
  X       <- as.data.frame(predict(dummies, df))
  X       <- rename_dummy_cols_clean(X)
  X$month_sale <- df$month_sale
  
  list(df=df, X=X, Y=Y, traits=traits)
}

# ------------------ LOAD & SELECT FLAGSHIP TREATMENTS ------------------
sel_table <- function() {
  if (!file.exists(RESULTS_Q)) stop("Cannot find results_with_q CSV: ", RESULTS_Q)
  res <- data.table::fread(RESULTS_Q)
  
  # ensure numeric for sort/filter
  for (col in c("kept_pct","ATE_pct","q_value_bh")) {
    if (is.character(res[[col]])) res[[col]] <- as.numeric(gsub(",", "", res[[col]]))
  }
  
  res <- res[k %in% K_SET & kept_pct >= KEPT_PCT_MIN]
  if (!nrow(res)) stop("No rows after filtering by K_SET and kept_pct.")
  
  # sort: q ascending, then |ATE_pct| descending
  res <- res[order(q_value_bh, -abs(ATE_pct))]
  res[, rank_k := seq_len(.N), by = k]
  res <- res[rank_k <= TOP_N_PER_K]
  if (!nrow(res)) stop("Top-N per K selection is empty. Loosen filters?")
  
  res[, treatment_file_stub := safe_combo_name(traits_combination)]
  res[]
}

# ------------------ PER-TREATMENT REPORTS ------------------
report_for_combo <- function(one, df, X, Y, traits) {
  combo <- one$traits_combination
  k     <- one$k
  stub  <- one$treatment_file_stub
  fam_in <- unique(trait_key(unlist(strsplit(combo, "\\+"))))
  
  message(sprintf("▶ CATE report for k=%d | %s", k, combo))
  
  # read per-unit CATE file
  f_units <- file.path(CATE_DIR, sprintf("CATE_units_k%d_%s.csv", k, stub))
  if (!file.exists(f_units)) {
    warning("Missing CATE units file: ", f_units)
    return(invisible(FALSE))
  }
  dt <- fread(f_units)  # row_id, cate_log, cate_pct, CI, etc.
  
  # join month for slicing
  dt <- merge(dt, df[, .(row_id, month_sale)], by="row_id", all.x=TRUE)
  
  # ---------- 1) GLOBAL CATE STATS ----------
  global <- data.table(
    treatment = combo, k = k,
    n = nrow(dt),
    median_cate_pct = median(dt$cate_pct),
    iqr_cate_pct_lo = quantile(dt$cate_pct, 0.25, names = FALSE),
    iqr_cate_pct_hi = quantile(dt$cate_pct, 0.75, names = FALSE),
    share_positive  = mean(dt$cate_pct > 0),
    min_pct = min(dt$cate_pct),
    max_pct = max(dt$cate_pct)
  )
  fwrite(global, file.path(OUT_DIR, sprintf("global_stats_k%d_%s.csv", k, stub)))
  
  # ---------- 2) GATES (deciles by τ̂%) ----------
  brks <- quantile(dt$cate_pct, seq(0,1,0.1), na.rm=TRUE)
  brks <- unique(as.numeric(brks))        # protect against duplicate breaks
  if (length(brks) >= 3) {
    dt[, gates_dec := cut(cate_pct, brks, include.lowest=TRUE, labels=FALSE)]
    gates_stats <- dt[, {
      ci <- boot_mean_ci_pct(cate_pct)
      .(n=.N, mean_pct=ci$mean, lo_pct=ci$lo, hi_pct=ci$hi)
    }, by=gates_dec][order(gates_dec)]
    
    fwrite(gates_stats, file.path(OUT_DIR, sprintf("gates_k%d_%s.csv", k, stub)))
    
    p_gates <- ggplot(gates_stats, aes(x=gates_dec, y=mean_pct)) +
      geom_point() +
      geom_errorbar(aes(ymin=lo_pct, ymax=hi_pct), width=0.2) +
      labs(x="CATE decile (low → high)", y="Mean CATE% (95% CI)",
           title=sprintf("GATES — k=%d: %s", k, combo)) +
      theme_minimal(base_size = 12)
    ggsave(file.path(OUT_DIR, sprintf("gates_k%d_%s.png", k, stub)), p_gates, width=7, height=4, dpi=200)
  } else {
    message("  • Skipping GATES (insufficient unique quantiles for deciles).")
  }
  
  # ---------- 3) SUBGROUP CATEs (families NOT in treatment) ----------
  fam_all <- c("Background","Clothes","Earring","Eyes","Fur","Hat","Mouth")
  fam_ok  <- setdiff(fam_all, fam_in)
  
  cols_ok <- unlist(lapply(fam_ok, function(tp) grep(paste0("^", tp, "_"), names(X), value = TRUE)))
  if (length(cols_ok)) {
    sub_dt  <- cbind(row_id = df$row_id, X[, cols_ok, drop=FALSE])
    sub_dt  <- as.data.table(sub_dt)
    dt <- merge(dt, sub_dt, by="row_id", all.x=TRUE)
    
    fam_candidates <- lapply(fam_ok, function(tp){
      lv <- grep(paste0("^", tp, "_"), names(dt), value=TRUE)
      if (!length(lv)) return(NULL)
      prev <- colMeans(as.matrix(dt[, ..lv]) > 0, na.rm=TRUE)
      common <- names(prev)[prev > 0.05]
      data.table(family=tp, n_common=length(common))
    })
    fam_rank <- rbindlist(Filter(Negate(is.null), fam_candidates), fill=TRUE)
    if (nrow(fam_rank)) {
      fam_rank <- fam_rank[order(-n_common)]
      fam_use  <- head(fam_rank$family, 2)
      
      for (tp in fam_use) {
        lv <- grep(paste0("^", tp, "_"), names(dt), value=TRUE)
        if (!length(lv)) next
        prev <- colMeans(as.matrix(dt[, ..lv]) > 0, na.rm=TRUE)
        lv_common <- names(prev)[prev > 0.05]
        if (!length(lv_common)) next
        
        tabs <- rbindlist(lapply(lv_common, function(v){
          idx <- dt[[v]] > 0
          n_i <- sum(idx)
          if (n_i < 100) return(NULL)
          ci <- boot_mean_ci_pct(dt$cate_pct[idx])
          data.table(
            treatment = combo, k = k, family = tp,
            level = sub(paste0("^", tp, "_"), "", v),
            n = n_i,
            mean_pct = ci$mean, lo_pct = ci$lo, hi_pct = ci$hi
          )
        }), fill=TRUE)
        
        if (is.data.frame(tabs) && nrow(tabs)) {
          tabs[, level := factor(level, levels=level[order(mean_pct, decreasing=TRUE)])]
          p_sub <- ggplot(tabs, aes(y=level, x=mean_pct)) +
            geom_point() +
            geom_errorbarh(aes(xmin=lo_pct, xmax=hi_pct), height=0.2) +
            labs(y=paste0(tp, " levels"), x="CATE% (95% CI)",
                 title=sprintf("Subgroup CATE — %s (k=%d)", tp, k)) +
            theme_minimal(base_size=12)
          ggsave(file.path(OUT_DIR, sprintf("subgroup_%s_k%d_%s.png", tp, k, stub)), p_sub, width=7, height=5, dpi=200)
          
          fwrite(tabs, file.path(OUT_DIR, sprintf("subgroup_%s_k%d_%s.csv", tp, k, stub)))
        }
      }
    } else {
      message("  • No eligible families for subgroup CATEs.")
    }
  } else {
    message("  • No non-treatment family dummies found in X; skipping subgroup plots.")
  }
  
  # ---------- 4) MONTH / SEASONALITY ----------
  mon <- dt[is.finite(month_sale) & month_sale >=1 & month_sale <= 12,
            { ci <- boot_mean_ci_pct(cate_pct); .(n=.N, mean_pct=ci$mean, lo_pct=ci$lo, hi_pct=ci$hi) },
            by=month_sale][order(month_sale)]
  if (nrow(mon)) {
    p_mon <- ggplot(mon, aes(x=month_sale, y=mean_pct)) +
      geom_line() + geom_point() +
      geom_ribbon(aes(ymin=lo_pct, ymax=hi_pct), alpha=0.15) +
      scale_x_continuous(breaks=1:12) +
      labs(x="Month of first sale", y="Mean CATE% (95% CI)",
           title=sprintf("Seasonality — k=%d: %s", k, combo)) +
      theme_minimal(base_size=12)
    ggsave(file.path(OUT_DIR, sprintf("seasonality_k%d_%s.png", k, stub)), p_mon, width=7, height=4, dpi=200)
    fwrite(mon, file.path(OUT_DIR, sprintf("seasonality_k%d_%s.csv", k, stub)))
  } else {
    message("  • Skipping seasonality (no month data after merges).")
  }
  
  # ---------- 5) OPTIONAL: REFIT for calibration/BLP/VI & ROBUSTNESS ----------
  if (REFIT_EACH) {
    set.seed(RNG_SEED)
    
    all_cols <- names(X)
    comb_vec <- strsplit(combo, "\\+")[[1]]
    
    # same-family drop
    same_fam_cols <- unique(unlist(lapply(fam_in, function(tp) {
      grep(paste0("^", tp, "_"), all_cols, value = TRUE)
    })))
    X_cov <- X[, setdiff(all_cols, c(comb_vec, same_fam_cols)), drop = FALSE]
    
    # treatment + ok rows
    W  <- as.numeric(rowSums(X[, comb_vec, drop=FALSE]) == k)
    ok <- stats::complete.cases(Y, W, X_cov)
    if (sum(ok) >= 100 && length(unique(W[ok])) == 2) {
      
      pf <- tryCatch(
        probability_forest(X_cov[ok,], W[ok], num.trees=2000, min.node.size=20, honesty=TRUE),
        error=function(e) regression_forest(X_cov[ok,], W[ok], num.trees=2000, min.node.size=20, honesty=TRUE)
      )
      ehat <- predict(pf)$predictions
      
      keep <- rep(FALSE, length(ok)); alpha_used <- NA_real_
      for (a in c(0.05, 0.02, 0.00)) {
        keep_try <- if (a == 0) ok else (ok & ehat >= a & ehat <= 1 - a)
        if (sum(keep_try) >= 100 && length(unique(W[keep_try])) == 2) {
          keep <- keep_try; alpha_used <- a; break
        }
      }
      if (sum(keep) >= 100) {
        of <- regression_forest(X_cov[keep,], Y[keep], num.trees=2000, min.node.size=20, honesty=TRUE)
        mu <- predict(of)$predictions
        cf <- causal_forest(X_cov[keep,], Y[keep], W[keep], W.hat = ehat[keep], Y.hat = mu,
                            num.trees=2000, min.node.size=20, honesty=TRUE)
        
        # Calibration (version-agnostic)
        tc <- tryCatch(test_calibration(cf), error=function(e) NULL)
        calib_row <- safe_extract_calibration(tc)
        calib_row[, `:=`(treatment = combo, k = k)]
        fwrite(calib_row, file.path(OUT_DIR, sprintf("calibration_k%d_%s.csv", k, stub)))
        
        # BLP (month only; version-agnostic accessor)
        Z <- as.matrix(data.frame(month_sale = X_cov$month_sale))
        blp <- tryCatch(best_linear_projection(cf, Z), warning=function(w) w, error=function(e) NULL)
        if (inherits(blp, "warning")) {
          # keep moving, still try to extract
          blp <- tryCatch(suppressWarnings(best_linear_projection(cf, Z)), error=function(e) NULL)
        }
        blpex <- safe_extract_blp(blp, "month_sale", default_index = 2L)
        if (!is.null(blpex)) {
          fwrite(data.table(treatment=combo, k=k, covariate="month_sale",
                            coef=blpex$coef, se=blpex$se),
                 file.path(OUT_DIR, sprintf("blp_k%d_%s.csv", k, stub)))
        } else {
          fwrite(data.table(treatment=combo, k=k, covariate="month_sale",
                            coef=NA_real_, se=NA_real_),
                 file.path(OUT_DIR, sprintf("blp_k%d_%s.csv", k, stub)))
        }
        
        # Variable importance
        vi <- tryCatch(variable_importance(cf), error=function(e) NULL)
        if (!is.null(vi)) {
          vi_tbl <- data.table(feature = colnames(X_cov), importance = as.numeric(vi))
          data.table::setorder(vi_tbl, -importance)
          vi_tbl[, `:=`(rank = .I, treatment = combo, k = k)]
          fwrite(vi_tbl[1:min(25, .N)], file.path(OUT_DIR, sprintf("varimp_top25_k%d_%s.csv", k, stub)))
        }
        
        # Robustness: smaller forest / larger node (also record treated-target ATEs)
        cf_r <- causal_forest(X_cov[keep,], Y[keep], W[keep],
                              W.hat = ehat[keep], Y.hat = mu,
                              num.trees=ROBUST_TREES, min.node.size=ROBUST_MIN_NODE, honesty=TRUE)
        
        ts_base <- if (!is.na(alpha_used) && alpha_used > 0) "overlap" else "all"
        ate0_all <- average_treatment_effect(cf,   target.sample = ts_base)
        ateR_all <- average_treatment_effect(cf_r, target.sample = ts_base)
        ate0_tr  <- average_treatment_effect(cf,   target.sample = "treated")
        ateR_tr  <- average_treatment_effect(cf_r, target.sample = "treated")
        
        robust_row <- data.table(
          treatment=combo, k=k, alpha_used=alpha_used, target_sample=ts_base,
          baseline_est_all = as.numeric(ate0_all["estimate"]),
          baseline_se_all  = as.numeric(ate0_all["std.err"]),
          robust_est_all   = as.numeric(ateR_all["estimate"]),
          robust_se_all    = as.numeric(ateR_all["std.err"]),
          baseline_est_treat = as.numeric(ate0_tr["estimate"]),
          baseline_se_treat  = as.numeric(ate0_tr["std.err"]),
          robust_est_treat   = as.numeric(ateR_tr["estimate"]),
          robust_se_treat    = as.numeric(ateR_tr["std.err"])
        )
        fwrite(robust_row, file.path(OUT_DIR, sprintf("robustness_k%d_%s.csv", k, stub)))
      }
    }
  }
  
  invisible(TRUE)
}  # <-- CLOSES report_for_combo()

# ------------------ RUN EVERYTHING ------------------
set.seed(RNG_SEED)
built <- build_df_X()
df <- as.data.table(built$df); X <- built$X; Y <- built$Y; traits <- built$traits
rm(built)

flagship <- sel_table()
data.table::fwrite(flagship, file.path(OUT_DIR, "selected_flagship_treatments.csv"))

for (i in seq_len(nrow(flagship))) {
  one <- flagship[i]
  report_for_combo(one, df, X, Y, traits)
}

message("✅ CATE reporting complete. Outputs written to: ", OUT_DIR)

# ------------------ OPTIONAL: Aggregate story tables ------------------
all_globals <- rbindlist(lapply(list.files(OUT_DIR, pattern="^global_stats_.*\\.csv$", full.names=TRUE), fread), fill=TRUE)
if (nrow(all_globals)) fwrite(all_globals, file.path(OUT_DIR, "ALL_global_stats.csv"))

all_gates <- rbindlist(lapply(list.files(OUT_DIR, pattern="^gates_.*\\.csv$", full.names=TRUE), fread), fill=TRUE)
if (nrow(all_gates)) fwrite(all_gates, file.path(OUT_DIR, "ALL_gates.csv"))

all_season <- rbindlist(lapply(list.files(OUT_DIR, pattern="^seasonality_.*\\.csv$", full.names=TRUE), fread), fill=TRUE)
if (nrow(all_season)) fwrite(all_season, file.path(OUT_DIR, "ALL_seasonality.csv"))
