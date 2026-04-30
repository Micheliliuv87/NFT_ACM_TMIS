
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fastDummies)
  library(grf)
  library(hdm)
})

set.seed(87)

# ============================================================
# M2 temporal-nuisance robustness based directly on M2V11 logic
#
# Purpose:
#   Address the Section 4.2 concern about whether
#   standard residualization is sufficient in a volatile NFT
#   secondary market with temporal dependence.
#
# Design principle:
#   Preserve the original M2V11 control architecture as the BASE.
#   Do NOT replace the model family. Do NOT move to daily FE or
#   rolling windows. Only add a modest temporal robustness layer
#   that matches the actual event-time structure of the data.
#
# Conditions:
#   1) base  : all Section 4.2 observations
#   2) samp  : Buyer(N-1) == Seller(N)
#   3) diffp : Buyer(N-1) != Seller(N)
#
# Robustness idea:
#   BASE_M2V11:
#     - exactly the M2V11-style nuisance structure for each condition
#   TEMPORAL_AUG_M2V11:
#     - BASE_M2V11 plus:
#         * lag_return_available indicator
#         * exact log holding period log_days_since_prev
#
# Saved outputs:
#   1) m2_temporal_sample_diagnostics.csv
#   2) cf_m2_M2V11_temporal_results.csv
#   3) cf_m2_M2V11_temporal_summary.csv
#   4) pds_m2_M2V11_temporal_results.csv
#   5) ols_m2_M2V11_month_residual_diag.csv
# ============================================================

RUN_CONDITIONS <- c("base", "samp", "diffp")
NUM_TREES      <- 2000
MIN_NODE_SIZE  <- 50
ALPHA_SEQ      <- c(0.10, 0.05, 0.02)
MIN_KEEP_N     <- 200
TRAIT_FAMILIES <- c("Background", "Clothes", "Earring", "Eyes", "Fur", "Hat", "Mouth")

DROP_COLLINEAR_DEFAULT <- c(
  "log_sellern_rolling_std_value_last10",
  "log_buyern_1_rolling_std_value_last10",
  "log_total_offers",
  "log_duration_offer_days",
  "log_buyern_rolling_std_value_last10",
  "log_sellern_1_rolling_std_value_last10"
)

winsorize_vec <- function(x, p = c(0.05, 0.95)) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(x)
  qs <- quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  x[x < qs[1]] <- qs[1]
  x[x > qs[2]] <- qs[2]
  x
}

winsorize_log1p_cols <- function(dt, cols, p = c(0.05, 0.95)) {
  cols <- intersect(cols, names(dt))
  if (!length(cols)) return(dt)
  for (cc in cols) {
    v <- winsorize_vec(suppressWarnings(as.numeric(dt[[cc]])), p = p)
    dt[[cc]] <- log1p(v)
    data.table::setnames(dt, cc, paste0("log_", cc))
  }
  dt
}

drop_constant_cols <- function(DT) {
  stopifnot(data.table::is.data.table(DT))
  keep <- vapply(DT, function(z) {
    u <- unique(z)
    sum(!is.na(u)) > 1
  }, logical(1))
  cols <- names(DT)[keep]
  if (!length(cols)) return(DT[, .(dummy_col = 0)])
  DT[, ..cols]
}

as_numeric_design <- function(df_like) {
  df <- as.data.frame(df_like)
  ok <- vapply(df, function(x) is.numeric(x) || is.integer(x) || is.logical(x), logical(1))
  if (!all(ok)) df <- df[ok]
  if (!ncol(df)) df$dummy_col <- 0
  stats::model.matrix(~ . - 1, data = df)
}

to_binary_with_meta <- function(v, var_name = NULL) {
  v <- suppressWarnings(as.numeric(v))
  v[!is.finite(v)] <- NA_real_
  if (all(is.na(v))) {
    return(list(w = rep(NA_integer_, length(v)), rule = NA_character_,
                cutoff_on_scale = NA_real_, cutoff_raw = NA_real_, share_high = NA_real_))
  }
  is_log <- !is.null(var_name) && startsWith(var_name, "log_")
  med <- stats::median(v, na.rm = TRUE)
  w   <- as.integer(v >= med)
  rule <- "median"; cut <- med
  if (length(unique(w[!is.na(w)])) < 2) {
    q75 <- suppressWarnings(as.numeric(quantile(v, 0.75, na.rm = TRUE)))
    if (is.finite(q75)) {
      w <- as.integer(v >= q75); rule <- "q75"; cut <- q75
    }
    if (length(unique(w[!is.na(w)])) < 2) {
      w <- as.integer(v > 0); rule <- ">0"; cut <- 0
      if (length(unique(w[!is.na(w)])) < 2) {
        return(list(w = rep(NA_integer_, length(v)), rule = NA_character_,
                    cutoff_on_scale = NA_real_, cutoff_raw = NA_real_, share_high = NA_real_))
      }
    }
  }
  cutoff_raw <- if (is.na(cut)) NA_real_ else if (is_log) expm1(cut) else cut
  list(w = w, rule = rule, cutoff_on_scale = cut, cutoff_raw = cutoff_raw,
       share_high = mean(w == 1, na.rm = TRUE))
}

compute_overlap_keep <- function(X_cov, W, y,
                                 alpha_seq = c(0.10, 0.05, 0.02),
                                 min_total = 200,
                                 min_per_arm = 5) {
  base_idx <- which(is.finite(y) & !is.na(W))
  Wbin <- as.integer(W > 0)
  out <- list(keep = integer(0), pHat = NULL, alpha = NA_real_,
              n_total = length(base_idx), n_kept = 0, kept_pct = 0)
  if (length(unique(Wbin[base_idx])) < 2 || length(base_idx) < min_total) return(out)
  
  pHat_base <- tryCatch({
    pf <- grf::probability_forest(X_cov[base_idx, , drop = FALSE], Wbin[base_idx],
                                  num.trees = 1000, honesty = TRUE)
    predict(pf)$predictions
  }, error = function(e) {
    rf <- grf::regression_forest(X_cov[base_idx, , drop = FALSE], Wbin[base_idx],
                                 num.trees = 1000, honesty = TRUE)
    predict(rf)$predictions
  })
  
  pHat_base <- pmin(pmax(pHat_base, 1e-6), 1 - 1e-6)
  
  for (a in alpha_seq) {
    ok <- if (a <= 0) seq_along(pHat_base) else which(pHat_base >= a & pHat_base <= (1 - a))
    keep <- base_idx[ok]
    tab <- table(Wbin[keep])
    if (length(keep) >= min_total && length(tab) == 2 && all(tab >= min_per_arm)) {
      pHat_all <- rep(NA_real_, length(W))
      pHat_all[base_idx] <- pHat_base
      out$keep <- keep
      out$pHat <- pHat_all
      out$alpha <- a
      out$n_kept <- length(keep)
      out$kept_pct <- 100 * length(keep) / length(base_idx)
      return(out)
    }
  }
  out
}

safe_calibration_v2 <- function(cf) {
  out <- list(intercept = NA_real_, slope = NA_real_, p_slope = NA_real_)
  cal <- tryCatch(grf::test_calibration(cf), error = function(e) NULL)
  if (is.null(cal)) return(out)
  if (inherits(cal, "lm")) {
    co <- tryCatch(coef(summary(cal)), error = function(e) NULL)
    if (!is.null(co)) {
      if ("(Intercept)" %in% rownames(co)) out$intercept <- unname(co["(Intercept)", "Estimate"])
      slope_row <- rownames(co)[grepl("pred", rownames(co), ignore.case = TRUE)]
      if (length(slope_row)) {
        out$slope <- unname(co[slope_row[1], "Estimate"])
        pcol <- intersect(colnames(co), c("Pr(>|t|)", "Pr(>|z|)", "p.value"))
        if (length(pcol)) out$p_slope <- unname(co[slope_row[1], pcol[1]])
      }
    }
  }
  if (is.list(cal) && !is.null(cal$regression) && inherits(cal$regression, "lm")) {
    co <- tryCatch(coef(summary(cal$regression)), error = function(e) NULL)
    if (!is.null(co)) {
      if ("(Intercept)" %in% rownames(co)) out$intercept <- unname(co["(Intercept)", "Estimate"])
      slope_row <- rownames(co)[grepl("pred", rownames(co), ignore.case = TRUE)]
      if (length(slope_row)) {
        out$slope <- unname(co[slope_row[1], "Estimate"])
        pcol <- intersect(colnames(co), c("Pr(>|t|)", "Pr(>|z|)", "p.value"))
        if (length(pcol)) out$p_slope <- unname(co[slope_row[1], pcol[1]])
      }
    }
  }
  out
}

compute_r2tau_cv <- function(X, y, K = 5, seed = 87) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X)
  if (!is.finite(n) || n < max(50, K + 5)) return(c(r2_tau = NA_real_, r2_tau_oos = NA_real_))
  set.seed(seed)
  foldid <- sample(rep(1:K, length.out = n))
  yhat <- rep(NA_real_, n)
  for (k in 1:K) {
    tr <- which(foldid != k); te <- which(foldid == k)
    if (length(tr) < 20 || length(te) < 5) next
    if (requireNamespace("glmnet", quietly = TRUE)) {
      fit <- try(glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr], alpha = 0, nfolds = 5), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        yhat[te] <- as.numeric(predict(fit, newx = X[te, , drop = FALSE], s = "lambda.min"))
      }
    }
    if (all(is.na(yhat[te]))) {
      fit2 <- try(hdm::rlasso(X[tr, , drop = FALSE], y[tr]), silent = TRUE)
      if (!inherits(fit2, "try-error")) {
        yhat[te] <- as.numeric(predict(fit2, newx = X[te, , drop = FALSE]))
      }
    }
  }
  if (all(is.na(yhat))) return(c(r2_tau = NA_real_, r2_tau_oos = NA_real_))
  vY <- stats::var(y, na.rm = TRUE)
  vYhat <- stats::var(yhat, na.rm = TRUE)
  mse <- mean((y - yhat)^2, na.rm = TRUE)
  c(
    r2_tau = if (is.finite(vY) && vY > 0) vYhat / vY else NA_real_,
    r2_tau_oos = if (is.finite(vY) && vY > 0) 1 - (mse / vY) else NA_real_
  )
}

parse_ate_obj <- function(obj) {
  if (is.null(obj)) return(list(est = NA_real_, se = NA_real_))
  if (is.list(obj)) {
    est <- NA_real_; se <- NA_real_
    if (!is.null(obj$estimate)) est <- as.numeric(obj$estimate)
    if (!is.null(obj$std.err))  se  <- as.numeric(obj$std.err)
    if (!is.null(obj$std.errs)) se  <- as.numeric(obj$std.errs)[1]
    return(list(est = est, se = se))
  }
  if (is.numeric(obj) && length(obj) >= 2) return(list(est = as.numeric(obj[1]), se = as.numeric(obj[2])))
  list(est = NA_real_, se = NA_real_)
}

standardize_trait_name <- function(x, trait_name) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- paste0("No_", trait_name)
  x
}

pool_rare_levels <- function(v, trait_name, min_prop = 0.01) {
  tab <- prop.table(table(v, useNA = "no"))
  rare <- setdiff(names(tab)[tab < min_prop], paste0("No_", trait_name))
  if (length(rare)) v[v %in% rare] <- paste0("Other_", trait_name)
  v
}

monthly_residual_diag <- function(month_str, resid_vec) {
  out <- data.frame(
    n_obs = length(resid_vec),
    n_months = NA_integer_,
    dw_like = NA_real_,
    ar1_coef = NA_real_,
    ar1_p = NA_real_,
    lb_p_lag1 = NA_real_,
    lb_p_lag2 = NA_real_,
    lb_p_lag3 = NA_real_
  )
  
  if (length(resid_vec) < 20) return(out)
  
  dd <- data.frame(
    month = as.character(month_str),
    resid = as.numeric(resid_vec),
    stringsAsFactors = FALSE
  )
  
  dd <- dd[is.finite(dd$resid) & !is.na(dd$month), , drop = FALSE]
  if (nrow(dd) == 0) return(out)
  
  agg <- stats::aggregate(resid ~ month, data = dd, FUN = function(z) mean(z, na.rm = TRUE))
  agg <- agg[order(agg$month), , drop = FALSE]
  
  out$n_months <- nrow(agg)
  if (nrow(agg) < 6) return(out)
  
  r <- agg$resid
  r_lag <- c(NA, head(r, -1))
  ok <- is.finite(r) & is.finite(r_lag)
  
  if (sum(ok) >= 5) {
    fit_ar1 <- lm(r[ok] ~ r_lag[ok])
    sm <- summary(fit_ar1)$coefficients
    out$ar1_coef <- if ("r_lag[ok]" %in% rownames(sm)) sm["r_lag[ok]", "Estimate"] else NA_real_
    out$ar1_p    <- if ("r_lag[ok]" %in% rownames(sm)) sm["r_lag[ok]", "Pr(>|t|)"] else NA_real_
    num <- sum(diff(r)^2, na.rm = TRUE)
    den <- sum(r^2, na.rm = TRUE)
    out$dw_like <- if (is.finite(den) && den > 0) num / den else NA_real_
  }
  
  maxlag <- max(1, min(3, floor(length(r) / 3)))
  out$lb_p_lag1 <- tryCatch(Box.test(r, lag = min(1, maxlag), type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  out$lb_p_lag2 <- tryCatch(Box.test(r, lag = min(2, maxlag), type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  out$lb_p_lag3 <- tryCatch(Box.test(r, lag = min(3, maxlag), type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  
  out
}

build_m2v11_condition_dataset <- function(condition = c("base", "samp", "diffp")) {
  condition <- match.arg(condition)
  
  df1 <- fread(
    "df_table1.csv",
    select = c("token_id", "time_n_sale", "time_n-1_sale", "time_n-2_sale",
               "price_n_sale", "price_n-1_sale", "price_n-2_sale", "price_1_sale",
               "buyer_n_sale", "seller_n_sale", "buyer_n-1_sale", "seller_n-1_sale")
  )
  df1 <- df1[!is.na(`time_n-1_sale`) & !is.na(`price_n-1_sale`) & (`price_n-1_sale` > 0)]
  
  df1[, time_n_sale_dt := as.POSIXct(time_n_sale, origin = "1970-01-01", tz = "UTC")]
  df1[, `time_n-1_sale_dt` := as.POSIXct(`time_n-1_sale`, origin = "1970-01-01", tz = "UTC")]
  df1[, days_since_prev := as.numeric(difftime(time_n_sale_dt, `time_n-1_sale_dt`, units = "days"))]
  df1 <- df1[!is.na(days_since_prev)]
  df1[, sold_after_30d := as.integer(days_since_prev >= 30)]
  
  if (condition == "samp") {
    df1 <- df1[`buyer_n-1_sale` == `seller_n_sale`]
  } else if (condition == "diffp") {
    df1 <- df1[`buyer_n-1_sale` != `seller_n_sale`]
  }
  
  df1[, lag_return_available := as.integer(!is.na(`price_n-2_sale`) & (`price_n-2_sale` > 0) & (`price_n-1_sale` > 0))]
  
  df_buyern    <- fread("df_table4.csv")
  df_sellern   <- fread("df_table5.csv")
  df_sellern_1 <- fread("df_table7.csv")
  df_traits    <- fread("df_table3.csv")
  df_offer     <- fread("Panel_for_Model2.csv",
                        select = c("token_id", "time_of_sale", "total_offers", "unique_makers_count", "median_offer_price", "duration_offer_days"))
  df_buyern_1 <- if (condition != "samp") fread("df_table6.csv") else NULL
  
  rename_side_features <- function(df, prefix) {
    old <- c("transaction_count", "active_period", "total_value", "total_gasUsed", "avg_gasPrice", "avg_gasLimit",
             "rolling_avg_value_last10", "rolling_std_value_last10")
    new <- c(
      paste0(prefix, "_tscount"), paste0(prefix, "_act_period"), paste0(prefix, "_total_value"),
      paste0(prefix, "_total_gasUsed"), paste0(prefix, "_avg_gasPrice"), paste0(prefix, "_avg_gasLimit"),
      paste0(prefix, "_rolling_avg_value_last10"), paste0(prefix, "_rolling_std_value_last10")
    )
    keep <- intersect(old, names(df))
    map <- setNames(new[match(keep, old)], keep)
    setnames(df, old = names(map), new = unname(map))
    df
  }
  
  df_buyern    <- rename_side_features(df_buyern,   "buyern")
  df_sellern   <- rename_side_features(df_sellern,  "sellern")
  df_sellern_1 <- rename_side_features(df_sellern_1,"sellern_1")
  if (!is.null(df_buyern_1)) df_buyern_1 <- rename_side_features(df_buyern_1, "buyern_1")
  
  dt <- copy(df1)
  dt <- merge(dt, df_buyern, by.x = "buyer_n_sale", by.y = names(df_buyern)[1], all.x = TRUE)
  dt <- merge(dt, df_sellern, by.x = "seller_n_sale", by.y = names(df_sellern)[1], all.x = TRUE)
  if (!is.null(df_buyern_1)) dt <- merge(dt, df_buyern_1, by.x = "buyer_n-1_sale", by.y = names(df_buyern_1)[1], all.x = TRUE)
  dt <- merge(dt, df_sellern_1, by.x = "seller_n-1_sale", by.y = names(df_sellern_1)[1], all.x = TRUE)
  
  df_offer_pre <- df_offer[time_of_sale == "time_n-1_sale"]
  setorder(df_offer_pre, token_id)
  df_offer_pre <- df_offer_pre[, .SD[.N], by = .(token_id)]
  offer_keep_cols <- setdiff(names(df_offer_pre), "time_of_sale")
  dt <- merge(dt, df_offer_pre[, ..offer_keep_cols], by = "token_id", all.x = TRUE)
  
  dt <- merge(dt, df_traits, by = "token_id", all.x = TRUE)
  data.table::setDT(dt)
  data.table::setalloccol(dt)
  
  drop_cols <- c("buyern_tscount", "buyern_act_period", "sellern_tscount", "sellern_act_period",
                 "sellern_1_tscount", "sellern_1_act_period")
  if (condition != "samp") drop_cols <- c(drop_cols, "buyern_1_tscount", "buyern_1_act_period")
  dt[, (intersect(drop_cols, names(dt))) := NULL]
  
  dt[, row_id := .I]
  dt[, time_n_sale_str := format(as.POSIXct(time_n_sale, origin = "1970-01-01", tz = "UTC"), "%Y-%m")]
  dt[, `time_n-1_sale_str` := format(as.POSIXct(`time_n-1_sale`, origin = "1970-01-01", tz = "UTC"), "%Y-%m")]
  dt[, c("time_n_sale", "time_n-1_sale", "time_n-2_sale", "time_n_sale_dt", "time_n-1_sale_dt") := NULL]
  
  month_n_vec  <- dt$time_n_sale_str
  month_n1_vec <- dt$`time_n-1_sale_str`
  
  WEI_TO_ETH <- 1 / 1e18
  gas_cost_eth_from_gwei <- function(units, price_wei) {
    u <- suppressWarnings(as.numeric(units))
    p <- suppressWarnings(as.numeric(price_wei))
    ifelse(is.finite(u) & is.finite(p), u * p * WEI_TO_ETH, NA_real_)
  }
  
  for (pref in c("buyern", if (condition != "samp") "buyern_1", "sellern", "sellern_1")) {
    used_col <- paste0(pref, "_total_gasUsed")
    price_col <- paste0(pref, "_avg_gasPrice")
    limit_col <- paste0(pref, "_avg_gasLimit")
    if (all(c(used_col, price_col) %in% names(dt))) dt[[paste0(pref, "_avg_total_gaspaid")]] <- gas_cost_eth_from_gwei(dt[[used_col]], dt[[price_col]])
    if (all(c(limit_col, price_col) %in% names(dt))) dt[[paste0(pref, "_avg_gasLimit_cost")]] <- gas_cost_eth_from_gwei(dt[[limit_col]], dt[[price_col]])
  }
  
  raw_gas_drop <- c(
    "buyern_total_gasUsed", "buyern_avg_gasPrice", "buyern_avg_gasLimit",
    "sellern_total_gasUsed", "sellern_avg_gasPrice", "sellern_avg_gasLimit",
    "sellern_1_total_gasUsed", "sellern_1_avg_gasPrice", "sellern_1_avg_gasLimit"
  )
  if (condition != "samp") raw_gas_drop <- c(raw_gas_drop, "buyern_1_total_gasUsed", "buyern_1_avg_gasPrice", "buyern_1_avg_gasLimit")
  dt[, (intersect(raw_gas_drop, names(dt))) := NULL]
  
  for (tr in TRAIT_FAMILIES) {
    if (!tr %in% names(dt)) dt[[tr]] <- NA_character_
    dt[[tr]] <- standardize_trait_name(dt[[tr]], tr)
    dt[[paste0(tr, "_recoded")]] <- pool_rare_levels(dt[[tr]], tr, min_prop = 0.01)
  }
  
  cat_cols <- c("time_n_sale_str", "time_n-1_sale_str", paste0(TRAIT_FAMILIES, "_recoded"))
  month_n_vec  <- dt$time_n_sale_str
  month_n1_vec <- dt$`time_n-1_sale_str`
  dt <- fastDummies::dummy_cols(
    dt,
    select_columns = intersect(cat_cols, names(dt)),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
  data.table::setDT(dt)
  data.table::setalloccol(dt)
  
  dt[, (intersect(TRAIT_FAMILIES, names(dt))) := NULL]
  
  cols_to_w <- c(
    "price_n_sale", "price_n-1_sale", "price_n-2_sale", "price_1_sale",
    "buyern_total_value", "buyern_rolling_avg_value_last10", "buyern_rolling_std_value_last10",
    "sellern_total_value", "sellern_rolling_avg_value_last10", "sellern_rolling_std_value_last10",
    "sellern_1_total_value", "sellern_1_rolling_avg_value_last10", "sellern_1_rolling_std_value_last10",
    "median_offer_price", "duration_offer_days", "total_offers", "unique_makers_count",
    "buyern_avg_total_gaspaid", "sellern_avg_total_gaspaid", "sellern_1_avg_total_gaspaid",
    "buyern_avg_gasLimit_cost", "sellern_avg_gasLimit_cost", "sellern_1_avg_gasLimit_cost"
  )
  if (condition != "samp") {
    cols_to_w <- c(cols_to_w,
                   "buyern_1_total_value", "buyern_1_rolling_avg_value_last10", "buyern_1_rolling_std_value_last10",
                   "buyern_1_avg_total_gaspaid", "buyern_1_avg_gasLimit_cost")
  }
  dt <- winsorize_log1p_cols(dt, cols = intersect(cols_to_w, names(dt)), p = c(0.05, 0.95))
  data.table::setDT(dt)
  data.table::setalloccol(dt)
  
  dt[, log_price_change := `log_price_n_sale` - `log_price_n-1_sale`]
  if ("log_price_n-2_sale" %in% names(dt)) {
    dt[, log_price_change_n_1 := `log_price_n-1_sale` - `log_price_n-2_sale`]
  } else {
    dt[, log_price_change_n_1 := 0]
  }
  dt[, log_days_since_prev := log1p(winsorize_vec(days_since_prev, c(0.05, 0.95)))]
  
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  for (nm in num_cols) {
    idx_na <- which(is.na(dt[[nm]]))
    if (length(idx_na)) data.table::set(dt, i = idx_na, j = nm, value = 0)
  }
  
  drop_from_X <- c(
    "token_id", "row_id", "buyer_n_sale", "seller_n_sale", "buyer_n-1_sale", "seller_n-1_sale",
    "log_price_change", "log_price_n_sale", "log_price_n-1_sale", "log_price_n-2_sale", "days_since_prev",
    "log_sellern_rolling_std_value_last10", "log_total_offers", "log_duration_offer_days",
    "log_buyern_rolling_std_value_last10", "log_sellern_1_rolling_std_value_last10"
  )
  if (condition != "samp") drop_from_X <- c(drop_from_X, "log_buyern_1_rolling_std_value_last10")
  X_all_dt <- dt[, setdiff(names(dt), intersect(drop_from_X, names(dt))), with = FALSE]
  X_all_dt <- drop_constant_cols(X_all_dt)
  
  list(
    dt = dt,
    X_all_dt = X_all_dt,
    y = dt$log_price_change,
    month_n = month_n_vec,
    month_n1 = month_n1_vec,
    sample_diag = data.table(
      condition = condition,
      n_obs = nrow(dt),
      n_lag_available = sum(dt$lag_return_available == 1, na.rm = TRUE),
      pct_lag_available = 100 * mean(dt$lag_return_available == 1, na.rm = TRUE),
      n_months_n = uniqueN(month_n_vec),
      n_months_n1 = uniqueN(month_n1_vec),
      median_days_since_prev = median(dt$days_since_prev, na.rm = TRUE),
      p75_days_since_prev = quantile(dt$days_since_prev, 0.75, na.rm = TRUE),
      p90_days_since_prev = quantile(dt$days_since_prev, 0.90, na.rm = TRUE)
    )
  )
}

make_condition_specs <- function(condition_obj) {
  temporal_extra <- c("lag_return_available", "log_days_since_prev")
  base_cols <- setdiff(colnames(condition_obj$X_all_dt), temporal_extra)
  aug_cols <- unique(c(base_cols, intersect(temporal_extra, names(condition_obj$dt))))
  list(
    BASE_M2V11 = base_cols,
    TEMPORAL_AUG_M2V11 = aug_cols
  )
}

treatment_var_map <- function(dt, condition) {
  map <- c(
    seller_N_total_value_high = "log_sellern_total_value",
    seller_N_1_total_value_high = "log_sellern_1_total_value",
    buyer_N_total_value_high = "log_buyern_total_value",
    seller_N_gaspaid_high = "log_sellern_avg_total_gaspaid",
    seller_N_1_gaspaid_high = "log_sellern_1_avg_total_gaspaid",
    buyer_N_gaspaid_high = "log_buyern_avg_total_gaspaid",
    seller_N_gaslimit_cost_high = "log_sellern_avg_gasLimit_cost",
    seller_N_1_gaslimit_cost_high = "log_sellern_1_avg_gasLimit_cost",
    buyer_N_gaslimit_cost_high = "log_buyern_avg_gasLimit_cost"
  )
  if (condition != "samp") {
    map <- c(map,
             buyer_N_1_total_value_high = "log_buyern_1_total_value",
             buyer_N_1_gaspaid_high = "log_buyern_1_avg_total_gaspaid",
             buyer_N_1_gaslimit_cost_high = "log_buyern_1_avg_gasLimit_cost")
  }
  map[map %in% names(dt)]
}

leakage_drop_for <- function(tname, condition) {
  map <- c(
    seller_N_total_value_high = "log_sellern_total_value",
    seller_N_1_total_value_high = "log_sellern_1_total_value",
    buyer_N_total_value_high = "log_buyern_total_value",
    seller_N_gaspaid_high = "log_sellern_avg_total_gaspaid",
    seller_N_1_gaspaid_high = "log_sellern_1_avg_total_gaspaid",
    buyer_N_gaspaid_high = "log_buyern_avg_total_gaspaid",
    seller_N_gaslimit_cost_high = "log_sellern_avg_gasLimit_cost",
    seller_N_1_gaslimit_cost_high = "log_sellern_1_avg_gasLimit_cost",
    buyer_N_gaslimit_cost_high = "log_buyern_avg_gasLimit_cost"
  )
  if (condition != "samp") {
    map <- c(map,
             buyer_N_1_total_value_high = "log_buyern_1_total_value",
             buyer_N_1_gaspaid_high = "log_buyern_1_avg_total_gaspaid",
             buyer_N_1_gaslimit_cost_high = "log_buyern_1_avg_gasLimit_cost")
  }
  unname(map[tname])
}

run_cf_for_condition <- function(condition_name, condition_obj) {
  dt <- condition_obj$dt
  y  <- condition_obj$y
  specs <- make_condition_specs(condition_obj)
  tmap  <- treatment_var_map(dt, condition_name)
  res_rows <- list()
  for (spec_name in names(specs)) {
    X_cols_full <- intersect(specs[[spec_name]], names(dt))
    X_all_dt <- drop_constant_cols(as.data.table(dt[, ..X_cols_full]))
    for (tname in names(tmap)) {
      tvar <- tmap[[tname]]
      bin <- to_binary_with_meta(dt[[tvar]], var_name = tvar)
      W <- bin$w
      if (all(is.na(W)) || length(unique(W[!is.na(W)])) < 2) next
      leak_drop <- leakage_drop_for(tname, condition_name)
      keep_names <- setdiff(names(X_all_dt), leak_drop)
      X_cov_dt <- drop_constant_cols(X_all_dt[, ..keep_names])
      X_cov <- as_numeric_design(X_cov_dt)
      ov <- compute_overlap_keep(X_cov, W, y, alpha_seq = ALPHA_SEQ, min_total = MIN_KEEP_N, min_per_arm = 5)
      if (!length(ov$keep)) next
      idx <- ov$keep
      Xk <- X_cov[idx, , drop = FALSE]
      Yk <- y[idx]
      Wk <- as.integer(W[idx])
      cf <- tryCatch(causal_forest(Xk, Yk, Wk, num.trees = NUM_TREES, honesty = TRUE, min.node.size = MIN_NODE_SIZE),
                     error = function(e) NULL)
      if (is.null(cf)) next
      ate_obj <- tryCatch(average_treatment_effect(cf, target.sample = ifelse(ov$alpha > 0, "overlap", "all")),
                          error = function(e) NULL)
      ate <- parse_ate_obj(ate_obj)
      ci_lo <- if (is.finite(ate$est) && is.finite(ate$se)) ate$est - 1.96 * ate$se else NA_real_
      ci_hi <- if (is.finite(ate$est) && is.finite(ate$se)) ate$est + 1.96 * ate$se else NA_real_
      pval <- if (is.finite(ate$est) && is.finite(ate$se) && ate$se > 0) 2 * pnorm(abs(ate$est / ate$se), lower.tail = FALSE) else NA_real_
      cate <- tryCatch(predict(cf)$predictions, error = function(e) rep(NA_real_, nrow(Xk)))
      r2s <- compute_r2tau_cv(Xk, cate, K = 5, seed = 87)
      cal <- safe_calibration_v2(cf)
      res_rows[[paste(condition_name, spec_name, tname, sep = "|")]] <- data.table(
        condition = condition_name, nuisance_spec = spec_name, treatment = tname,
        treatment_var = tvar, split_rule = bin$rule, cutoff_raw = bin$cutoff_raw, share_high = bin$share_high,
        ATE_log = ate$est, SE_log = ate$se,
        ATE_pct = 100 * (exp(ate$est) - 1),
        CI95_lo_pct = 100 * (exp(ci_lo) - 1),
        CI95_hi_pct = 100 * (exp(ci_hi) - 1),
        p_value = pval, kept_pct = ov$kept_pct, alpha = ov$alpha, n_kept = ov$n_kept,
        calib_intercept = cal$intercept, calib_slope = cal$slope, calib_p_slope = cal$p_slope,
        r2_tau = r2s[["r2_tau"]], r2_tau_oos = r2s[["r2_tau_oos"]]
      )
    }
  }
  if (!length(res_rows)) return(NULL)
  rbindlist(res_rows, use.names = TRUE, fill = TRUE)
}

run_pds_for_condition <- function(condition_name, condition_obj) {
  dt <- condition_obj$dt
  y  <- condition_obj$y
  specs <- make_condition_specs(condition_obj)
  tmap  <- treatment_var_map(dt, condition_name)
  focal_treats <- intersect(c("seller_N_total_value_high", "seller_N_1_total_value_high", "buyer_N_total_value_high", "buyer_N_1_total_value_high"), names(tmap))
  res_rows <- list()
  for (spec_name in names(specs)) {
    X_cols_full <- intersect(specs[[spec_name]], names(dt))
    X_all_dt <- drop_constant_cols(as.data.table(dt[, ..X_cols_full]))
    for (tname in focal_treats) {
      tvar <- tmap[[tname]]
      bin <- to_binary_with_meta(dt[[tvar]], var_name = tvar)
      W <- bin$w
      if (all(is.na(W)) || length(unique(W[!is.na(W)])) < 2) next
      leak_drop <- leakage_drop_for(tname, condition_name)
      keep_names <- setdiff(names(X_all_dt), leak_drop)
      X_cov_dt <- drop_constant_cols(X_all_dt[, ..keep_names])
      X_cov <- as_numeric_design(X_cov_dt)
      ov <- compute_overlap_keep(X_cov, W, y, alpha_seq = ALPHA_SEQ, min_total = MIN_KEEP_N, min_per_arm = 5)
      if (!length(ov$keep)) next
      idx <- ov$keep
      Xk <- X_cov[idx, , drop = FALSE]
      Yk <- y[idx]
      Wk <- as.integer(W[idx])
      fit <- tryCatch(hdm::rlassoEffect(x = Xk, y = Yk, d = Wk, method = "double selection"), error = function(e) NULL)
      if (is.null(fit)) next
      est <- tryCatch(as.numeric(coef(fit)), error = function(e) NA_real_)
      se  <- tryCatch(as.numeric(sqrt(diag(vcov(fit))))[1], error = function(e) NA_real_)
      if (!is.finite(se) || se <= 0) {
        ci <- c(NA_real_, NA_real_); pval <- NA_real_
      } else {
        ci <- est + c(-1.96, 1.96) * se
        pval <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
      }
      res_rows[[paste(condition_name, spec_name, tname, sep = "|")]] <- data.table(
        condition = condition_name, nuisance_spec = spec_name, treatment = tname, treatment_var = tvar,
        split_rule = bin$rule, cutoff_raw = bin$cutoff_raw, share_high = bin$share_high,
        ATE_log = est, SE_log = se,
        ATE_pct = 100 * (exp(est) - 1),
        CI95_lo_pct = 100 * (exp(ci[1]) - 1),
        CI95_hi_pct = 100 * (exp(ci[2]) - 1),
        p_value = pval, kept_pct = ov$kept_pct, alpha = ov$alpha, n_kept = ov$n_kept
      )
    }
  }
  if (!length(res_rows)) return(NULL)
  rbindlist(res_rows, use.names = TRUE, fill = TRUE)
}

run_month_resid_diag <- function(condition_name, condition_obj) {
  dt <- condition_obj$dt
  y  <- condition_obj$y
  specs <- make_condition_specs(condition_obj)
  rows <- list()
  for (spec_name in names(specs)) {
    X_cols_full <- intersect(specs[[spec_name]], names(dt))
    X_cols_full <- setdiff(X_cols_full, DROP_COLLINEAR_DEFAULT)
    Xdt <- drop_constant_cols(as.data.table(dt[, ..X_cols_full]))
    X <- as_numeric_design(Xdt)
    fit <- tryCatch(lm(y ~ ., data = data.frame(y = y, X)), error = function(e) NULL)
    if (is.null(fit)) next
    stopifnot(length(condition_obj$month_n) == length(residuals(fit)))
    d <- monthly_residual_diag(condition_obj$month_n, residuals(fit))
    rows[[paste(condition_name, spec_name, sep = "|")]] <- data.table(
      condition = condition_name, nuisance_spec = spec_name,
      n_obs = d$n_obs, n_months = d$n_months, dw_like = d$dw_like,
      ar1_coef = d$ar1_coef, ar1_p = d$ar1_p,
      lb_p_lag1 = d$lb_p_lag1, lb_p_lag2 = d$lb_p_lag2, lb_p_lag3 = d$lb_p_lag3
    )
  }
  if (!length(rows)) return(NULL)
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

sample_diags <- list()
all_cf <- list()
all_pds <- list()
all_diag <- list()

for (cond in RUN_CONDITIONS) {
  message("==== Running condition: ", cond, " ====")
  obj <- build_m2v11_condition_dataset(cond)
  sample_diags[[cond]] <- obj$sample_diag
  cf_res <- run_cf_for_condition(cond, obj)
  if (!is.null(cf_res)) all_cf[[cond]] <- cf_res
  pds_res <- run_pds_for_condition(cond, obj)
  if (!is.null(pds_res)) all_pds[[cond]] <- pds_res
  diag_res <- run_month_resid_diag(cond, obj)
  if (!is.null(diag_res)) all_diag[[cond]] <- diag_res
}

sample_diag_out <- if (length(sample_diags)) rbindlist(sample_diags, use.names = TRUE, fill = TRUE) else data.table()
cf_out   <- if (length(all_cf))   rbindlist(all_cf,   use.names = TRUE, fill = TRUE) else data.table()
pds_out  <- if (length(all_pds))  rbindlist(all_pds,  use.names = TRUE, fill = TRUE) else data.table()
diag_out <- if (length(all_diag)) rbindlist(all_diag, use.names = TRUE, fill = TRUE) else data.table()

if (nrow(sample_diag_out)) fwrite(sample_diag_out, "m2_temporal_sample_diagnostics.csv")

if (nrow(cf_out)) {
  fwrite(cf_out, "cf_m2_M2V11_temporal_results.csv")
  cf_wide <- dcast(
    cf_out,
    condition + treatment ~ nuisance_spec,
    value.var = c("ATE_pct", "CI95_lo_pct", "CI95_hi_pct", "p_value", "kept_pct", "r2_tau", "calib_slope"),
    sep = "__"
  )
  if (all(c("ATE_pct__BASE_M2V11", "ATE_pct__TEMPORAL_AUG_M2V11") %in% names(cf_wide))) {
    cf_wide[, same_sign := sign(ATE_pct__BASE_M2V11) == sign(ATE_pct__TEMPORAL_AUG_M2V11)]
    cf_wide[, abs_delta_ate := abs(ATE_pct__TEMPORAL_AUG_M2V11 - ATE_pct__BASE_M2V11)]
    cf_wide[, pct_change_vs_base := ifelse(abs(ATE_pct__BASE_M2V11) > 1e-8,
                                           100 * (ATE_pct__TEMPORAL_AUG_M2V11 - ATE_pct__BASE_M2V11) / abs(ATE_pct__BASE_M2V11),
                                           NA_real_)]
  }
  fwrite(cf_wide, "cf_m2_M2V11_temporal_summary.csv")
}

if (nrow(pds_out)) fwrite(pds_out, "pds_m2_M2V11_temporal_results.csv")
if (nrow(diag_out)) fwrite(diag_out, "ols_m2_M2V11_month_residual_diag.csv")

cat("\nSaved files:\n")
cat("  - m2_temporal_sample_diagnostics.csv\n")
cat("  - cf_m2_M2V11_temporal_results.csv\n")
cat("  - cf_m2_M2V11_temporal_summary.csv\n")
cat("  - pds_m2_M2V11_temporal_results.csv\n")
cat("  - ols_m2_M2V11_month_residual_diag.csv\n")
