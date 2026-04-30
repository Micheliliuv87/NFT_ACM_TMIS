suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(fastDummies)
  library(hdm)   # PDS: rlassoEffect()
  library(grf)   # ONLY for overlap trimming (probability_forest); no GRF effects
})

set.seed(87)

# ======================= Helpers =======================
winsorize_vec <- function(x, p = c(0.05, 0.95)) {
  x <- as.numeric(x)
  qs <- quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  x[x < qs[1]] <- qs[1]
  x[x > qs[2]] <- qs[2]
  x
}

winsorize_log1p_cols <- function(dt, cols, p = c(0.05, 0.95)) {
  cols <- intersect(cols, names(dt))
  if (!length(cols)) return(dt)
  for (cc in cols) {
    v <- winsorize_vec(suppressWarnings(as.numeric(dt[[cc]])), p)
    dt[[cc]] <- log1p(v)
    data.table::setnames(dt, cc, paste0("log_", cc))
  }
  dt
}

# Binary split with metadata; handles log_* names for raw cutoffs
to_binary_with_meta <- function(v, var_name = NULL) {
  v <- suppressWarnings(as.numeric(v))
  v[!is.finite(v)] <- NA_real_
  if (all(is.na(v))) {
    return(list(w = rep(NA_integer_, length(v)),
                rule = NA_character_, cutoff_on_scale = NA_real_,
                cutoff_raw = NA_real_, share_high = NA_real_))
  }
  is_log <- !is.null(var_name) && startsWith(var_name, "log_")
  # pick cut: median → q75 → >0
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
        return(list(w = rep(NA_integer_, length(v)),
                    rule = NA_character_, cutoff_on_scale = NA_real_,
                    cutoff_raw = NA_real_, share_high = NA_real_))
      }
    }
  }
  cutoff_raw <- if (is.na(cut)) NA_real_ else if (is_log) expm1(cut) else cut
  list(
    w = w,
    rule = rule,
    cutoff_on_scale = cut,
    cutoff_raw = cutoff_raw,
    share_high = mean(w == 1, na.rm = TRUE)
  )
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
  if (!all(ok)) {
    bad <- names(df)[!ok]
    message("[design] Dropping non-numeric columns: ", paste(bad, collapse = ", "))
    df <- df[ok]
  }
  if (!ncol(df)) df$dummy_col <- 0
  stats::model.matrix(~ . - 1, data = df)
}

smd_vec <- function(X, W) {
  X <- as.matrix(X); W <- as.integer(W)
  treat <- X[W==1,,drop=FALSE]; ctrl <- X[W==0,,drop=FALSE]
  if (nrow(treat) < 2 || nrow(ctrl) < 2) return(rep(NA_real_, ncol(X)))
  m1 <- colMeans(treat); m0 <- colMeans(ctrl)
  v1 <- apply(treat, 2, var); v0 <- apply(ctrl, 2, var)
  sd_pool <- sqrt(pmax(1e-12, 0.5*(v1+v0)))
  (m1 - m0) / sd_pool
}

quant_str <- function(x) {
  qs <- quantile(x, probs = c(0, .05, .25, .50, .75, .95, 1), na.rm = TRUE)
  paste0("min=", sprintf("%.3f", qs[1]),
         ", p5=", sprintf("%.3f", qs[2]),
         ", p25=", sprintf("%.3f", qs[3]),
         ", p50=", sprintf("%.3f", qs[4]),
         ", p75=", sprintf("%.3f", qs[5]),
         ", p95=", sprintf("%.3f", qs[6]),
         ", max=", sprintf("%.3f", qs[7]))
}

# ---- Overlap with probability_forest; stricter trimming by default ----
compute_overlap_keep <- function(X_cov, W, y, alpha_seq = c(0.10, 0.05, 0.02),
                                 min_total = 200, min_per_arm = 5) {
  base_idx <- which(is.finite(y) & !is.na(W))
  Wbin <- as.integer(W > 0)
  out <- list(keep = integer(0), pHat = NULL, alpha = NA_real_,
              n_total = length(base_idx), n_kept = 0, kept_pct = 0)
  if (length(unique(Wbin[base_idx])) < 2 || length(base_idx) < min_total) return(out)
  
  pHat_base <- tryCatch({
    pf <- grf::probability_forest(X_cov[base_idx, , drop=FALSE], Wbin[base_idx],
                                  num.trees = 1000, honesty = TRUE)
    predict(pf)$predictions
  }, error = function(e) {
    rf <- grf::regression_forest(X_cov[base_idx, , drop=FALSE], Wbin[base_idx],
                                 num.trees = 1000, honesty = TRUE)
    predict(rf)$predictions
  })
  pHat_base <- pmin(pmax(pHat_base, 1e-6), 1-1e-6)
  
  for (a in alpha_seq) {
    ok <- if (a <= 0) seq_along(pHat_base) else which(pHat_base >= a & pHat_base <= (1 - a))
    keep <- base_idx[ok]
    tab  <- table(Wbin[keep])
    if (length(keep) >= min_total && length(tab) == 2 && all(tab >= min_per_arm)) {
      pHat_all <- rep(NA_real_, length(W)); pHat_all[base_idx] <- pHat_base
      out$keep <- keep; out$pHat <- pHat_all; out$alpha <- a
      out$n_kept <- length(keep); out$kept_pct <- 100*length(keep)/length(base_idx)
      return(out)
    }
  }
  out
}

# ======================= Load & engineer =======================
df1 <- fread(
  "df_table1.csv",
  select = c("token_id","time_n_sale","time_n-1_sale","time_n-2_sale",
             "price_n_sale","price_n-1_sale","price_n-2_sale", "price_1_sale", 
             "buyer_n_sale","seller_n_sale","buyer_n-1_sale","seller_n-1_sale")
)
df1 <- df1[!is.na(`time_n-1_sale`) & !is.na(`price_n-1_sale`) & (`price_n-1_sale` > 0)]

df1[, time_n_sale_dt     := as.POSIXct(time_n_sale,     origin = "1970-01-01", tz = "UTC")]
df1[, `time_n-1_sale_dt` := as.POSIXct(`time_n-1_sale`, origin = "1970-01-01", tz = "UTC")]
df1[, days_since_prev    := as.numeric(difftime(time_n_sale_dt, `time_n-1_sale_dt`, units = "days"))]
df1 <- df1[!is.na(days_since_prev)]
df1[, sold_after_30d := as.integer(days_since_prev >= 30)]

# Sides, traits, offers
df_buyern    <- fread("df_table4.csv")
df_buyern_1  <- fread("df_table6.csv")
df_sellern   <- fread("df_table5.csv")
df_sellern_1 <- fread("df_table7.csv")
df_traits    <- fread("df_table3.csv")
df_offer     <- fread("Panel_for_Model2.csv",
                      select = c("token_id","time_of_sale","total_offers","unique_makers_count",
                                 "median_offer_price", "duration_offer_days"))

# Standardize side cols
data.table::setnames(df_buyern,
                     c("transaction_count","active_period","total_value","total_gasUsed","avg_gasPrice","avg_gasLimit",
                       "rolling_avg_value_last10","rolling_std_value_last10"),
                     c("buyern_tscount","buyern_act_period","buyern_total_value","buyern_total_gasUsed","buyern_avg_gasPrice",
                       "buyern_avg_gasLimit","buyern_rolling_avg_value_last10","buyern_rolling_std_value_last10"))
data.table::setnames(df_buyern_1,
                     c("transaction_count","active_period","total_value","total_gasUsed","avg_gasPrice","avg_gasLimit",
                       "rolling_avg_value_last10","rolling_std_value_last10"),
                     c("buyern_1_tscount","buyern_1_act_period","buyern_1_total_value","buyern_1_total_gasUsed","buyern_1_avg_gasPrice",
                       "buyern_1_avg_gasLimit","buyern_1_rolling_avg_value_last10","buyern_1_rolling_std_value_last10"))
data.table::setnames(df_sellern,
                     c("transaction_count","active_period","total_value","total_gasUsed","avg_gasPrice","avg_gasLimit",
                       "rolling_avg_value_last10","rolling_std_value_last10"),
                     c("sellern_tscount","sellern_act_period","sellern_total_value","sellern_total_gasUsed","sellern_avg_gasPrice",
                       "sellern_avg_gasLimit","sellern_rolling_avg_value_last10","sellern_rolling_std_value_last10"))
data.table::setnames(df_sellern_1,
                     c("transaction_count","active_period","total_value","total_gasUsed","avg_gasPrice","avg_gasLimit",
                       "rolling_avg_value_last10","rolling_std_value_last10"),
                     c("sellern_1_tscount","sellern_1_act_period","sellern_1_total_value","sellern_1_total_gasUsed","sellern_1_avg_gasPrice",
                       "sellern_1_avg_gasLimit","sellern_1_rolling_avg_value_last10","sellern_1_rolling_std_value_last10"))

# Merge
dt <- merge(df1, df_buyern,    by.x = "buyer_n_sale",    by.y = "buyer_n_address",     all.x = TRUE)
dt <- merge(dt,  df_sellern,   by.x = "seller_n_sale",   by.y = "seller_n_address",    all.x = TRUE)
dt <- merge(dt,  df_buyern_1,  by.x = "buyer_n-1_sale",  by.y = "buyer_n-1_address",   all.x = TRUE)
dt <- merge(dt,  df_sellern_1, by.x = "seller_n-1_sale", by.y = "seller_n-1_address",  all.x = TRUE)

df_offer_pre <- df_offer[time_of_sale == "time_n-1_sale"]
setorder(df_offer_pre, token_id)
df_offer_pre <- df_offer_pre[, .SD[.N], by = token_id]
dt <- merge(dt, df_offer_pre[, !"time_of_sale"], by = "token_id", all.x = TRUE)
dt <- merge(dt, df_traits, by = "token_id", all.x = TRUE)

# Row id
dt[, row_id := .I]

# Time FE strings (month)
dt[, time_n_sale_str     := format(as.POSIXct(time_n_sale,     origin="1970-01-01", tz="UTC"), "%Y-%m")]
dt[, `time_n-1_sale_str` := format(as.POSIXct(`time_n-1_sale`, origin="1970-01-01", tz="UTC"), "%Y-%m")]

# Remove raw stamps so they never leak into X
cols_never_in_X <- c("time_n_sale", "time_n-1_sale", "time_n-2_sale",
                     "time_n_sale_dt", "time_n-1_sale_dt","rarity.rank","time_of_sale")
dt[, (intersect(cols_never_in_X, names(dt))) := NULL]

# ---------------- Gas features (ETH) then drop raw gas parts ----------------
WEI_TO_ETH <- 1/1e18
gas_cost_eth_from_wei <- function(units, price_wei) {
  u <- suppressWarnings(as.numeric(units))
  p <- suppressWarnings(as.numeric(price_wei))
  ifelse(is.finite(u) & is.finite(p), u * p * WEI_TO_ETH, NA_real_)
}
dt[, buyern_avg_total_gaspaid    := gas_cost_eth_from_wei(buyern_total_gasUsed,    buyern_avg_gasPrice)]
dt[, buyern_1_avg_total_gaspaid  := gas_cost_eth_from_wei(buyern_1_total_gasUsed,  buyern_1_avg_gasPrice)]
dt[, sellern_avg_total_gaspaid   := gas_cost_eth_from_wei(sellern_total_gasUsed,   sellern_avg_gasPrice)]
dt[, sellern_1_avg_total_gaspaid := gas_cost_eth_from_wei(sellern_1_total_gasUsed, sellern_1_avg_gasPrice)]
dt[, buyern_avg_gasLimit_cost    := gas_cost_eth_from_wei(buyern_avg_gasLimit,     buyern_avg_gasPrice)]
dt[, buyern_1_avg_gasLimit_cost  := gas_cost_eth_from_wei(buyern_1_avg_gasLimit,   buyern_1_avg_gasPrice)]
dt[, sellern_avg_gasLimit_cost   := gas_cost_eth_from_wei(sellern_avg_gasLimit,    sellern_avg_gasPrice)]
dt[, sellern_1_avg_gasLimit_cost := gas_cost_eth_from_wei(sellern_1_avg_gasLimit,  sellern_1_avg_gasPrice)]

# Drop raw gas internals (avoid leakage & re-logging)
dt[, c("buyern_total_gasUsed","buyern_avg_gasPrice",
       "buyern_1_total_gasUsed","buyern_1_avg_gasPrice",
       "sellern_total_gasUsed","sellern_avg_gasPrice",
       "sellern_1_total_gasUsed","sellern_1_avg_gasPrice",
       "buyern_avg_gasLimit","buyern_1_avg_gasLimit",
       "sellern_avg_gasLimit","sellern_1_avg_gasLimit") := NULL]

# ---------------- Trait recoding ----------------
traits <- c("Background","Clothes","Earring","Eyes","Fur","Hat","Mouth")
threshold <- 0.01
for (tr in traits) {
  if (!tr %in% names(dt)) dt[[tr]] <- NA_character_
  dt[[tr]][is.na(dt[[tr]]) | dt[[tr]] == ""] <- paste0("No_", tr)
}
for (tr in traits) {
  ft   <- prop.table(table(dt[[tr]]))
  rare <- setdiff(names(ft)[ft < threshold], paste0("No_", tr))
  dt[[paste0(tr, "_recoded")]] <- ifelse(dt[[tr]] %in% rare, paste0("Other_", tr), dt[[tr]])
}

# One-hot months + recoded traits
cat_cols <- c("time_n_sale_str","time_n-1_sale_str", paste0(traits, "_recoded"))
dt <- fastDummies::dummy_cols(
  dt,
  select_columns = cat_cols,
  remove_first_dummy = TRUE,
  remove_selected_columns = TRUE
)
dt[, c("Background","Clothes","Earring","Eyes","Fur","Hat","Mouth") := NULL]

# ---------------- Winsorize + log1p ----------------
cols_to_w <- c(
  "price_n_sale","price_n-1_sale","price_n-2_sale", "price_1_sale",
  "buyern_total_value","buyern_rolling_avg_value_last10","buyern_rolling_std_value_last10",
  "sellern_total_value","sellern_rolling_avg_value_last10","sellern_rolling_std_value_last10",
  "buyern_1_total_value","buyern_1_rolling_avg_value_last10","buyern_1_rolling_std_value_last10",
  "sellern_1_total_value","sellern_1_rolling_avg_value_last10","sellern_1_rolling_std_value_last10",
  "total_offers","unique_makers_count","median_offer_price","duration_offer_days",
  # Keep gas in X:
  "buyern_avg_total_gaspaid","buyern_1_avg_total_gaspaid",
  "sellern_avg_total_gaspaid","sellern_1_avg_total_gaspaid",
  "buyern_avg_gasLimit_cost","buyern_1_avg_gasLimit_cost",
  "sellern_avg_gasLimit_cost","sellern_1_avg_gasLimit_cost"
)
dt <- winsorize_log1p_cols(dt, cols_to_w)

# ---------------- Outcome ----------------
dt[, log_price_change := `log_price_n_sale` - `log_price_n-1_sale`]
if ("log_price_n-2_sale" %in% names(dt)) {
  dt[, `log_price_change_n-1` := `log_price_n-1_sale` - `log_price_n-2_sale`]
} else {
  dt[, `log_price_change_n-1` := 0]
}

# ---------------- Fill remaining NA numerics with 0 ----------------
num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
for (nm in num_cols) data.table::set(dt, i = which(is.na(dt[[nm]])), j = nm, value = 0)

# ---------------- Design matrix (X) ----------------
# Exclude: identifiers, outcome pieces, ALL *tscount*, and days_since_prev from X
drop_from_X <- intersect(
  c("token_id","row_id","buyer_n_sale","seller_n_sale","buyer_n-1_sale","seller_n-1_sale",
    "log_price_change","log_price_n_sale","log_price_n-1_sale","log_price_n-2_sale",
    "days_since_prev",
    "buyern_tscount","buyern_1_tscount","sellern_tscount","sellern_1_tscount", 
    "log_sellern_rolling_std_value_last10", "log_buyern_1_rolling_std_value_last10",
    "log_total_offers", "log_duration_offer_days", "log_buyern_rolling_std_value_last10", 
    "log_sellern_1_rolling_std_value_last10"),
  names(dt)
)
X_all_dt <- dt[, setdiff(names(dt), drop_from_X), with = FALSE]
X_all_dt <- drop_constant_cols(X_all_dt)
X_all <- as_numeric_design(X_all_dt)
y <- dt$log_price_change

### >>> NEW: Multicollinearity diagnostics (GLOBAL on full X pool) -------------
# Config
.top_corr_threshold <- 0.90
.max_cols_for_corr  <- 1500  # sample if too wide for a quick pass

# Condition number of standardized full X pool
X_std <- scale(as.matrix(X_all_dt))
X_std[!is.finite(X_std)] <- 0
sv <- svd(X_std, nu = 0, nv = 0)
cond_number <- if (length(sv$d)) max(sv$d, na.rm=TRUE) / pmax(min(sv$d, na.rm=TRUE), 1e-12) else NA_real_
data.table::fwrite(
  data.frame(scope = "global_full_X", condition_number = cond_number),
  "PDS_multicollinearity_global_condition_number.csv"
)

# Top |corr| pairs among controls (sample columns if very wide)
cols_check <- colnames(X_all_dt)
if (length(cols_check) > .max_cols_for_corr) { set.seed(87); cols_check <- sample(cols_check, .max_cols_for_corr) }
Xc <- as.matrix(X_all_dt[, ..cols_check])
cr <- suppressWarnings(cor(Xc, use = "pairwise.complete.obs"))
cr[upper.tri(cr, diag = TRUE)] <- NA
hits <- which(abs(cr) >= .top_corr_threshold, arr.ind = TRUE)
if (nrow(hits)) {
  out_pairs <- data.frame(
    var1 = rownames(cr)[hits[, "row"]],
    var2 = colnames(cr)[hits[, "col"]],
    r    = cr[hits],
    stringsAsFactors = FALSE
  )
  out_pairs <- out_pairs[order(-abs(out_pairs$r)), ]
  data.table::fwrite(out_pairs, "PDS_multicollinearity_global_top_pairs.csv")
} else {
  data.table::fwrite(
    data.frame(note = sprintf("No |corr| >= %.2f in sampled set", .top_corr_threshold)),
    "PDS_multicollinearity_global_top_pairs.csv"
  )
}
### <<< NEW -------------------------------------------------------------------

# ================= Treatments (ONLY total value, logged) =================
treatments_info <- list(
  seller_N_total_value_high   = "log_sellern_total_value",
  seller_N_1_total_value_high = "log_sellern_1_total_value",
  buyer_N_total_value_high    = "log_buyern_total_value",
  buyer_N_1_total_value_high  = "log_buyern_1_total_value"
)

bin_treatments <- list()
treat_threshold_meta <- list()

for (tname in names(treatments_info)) {
  var_name <- treatments_info[[tname]]
  if (!var_name %in% names(dt)) {
    message(sprintf("[treat] %s → missing column %s; skipping.", tname, var_name))
    next
  }
  res <- to_binary_with_meta(dt[[var_name]], var_name = var_name)
  w <- res$w
  if (all(is.na(w)) || length(unique(w[!is.na(w)])) < 2) {
    message(sprintf("[treat] %s → no variation after binarization; skipping.", tname))
    next
  }
  bin_treatments[[tname]] <- w
  treat_threshold_meta[[tname]] <- data.frame(
    treatment        = tname,
    split_rule       = res$rule,
    cutoff_on_scale  = res$cutoff_on_scale,
    cutoff_raw       = res$cutoff_raw,
    share_high       = res$share_high,
    stringsAsFactors = FALSE
  )
  cat(sprintf("[treat] %-28s rule=%-6s cutoff(on-scale)=% .3f cutoff_raw=% .3f share_high=%.1f%%\n",
              tname, res$rule, res$cutoff_on_scale, res$cutoff_raw, 100*res$share_high))
}
if (!length(bin_treatments)) stop("No usable treatments after binarization. Check inputs/variation.")
treat_threshold_df <- dplyr::bind_rows(treat_threshold_meta)
data.table::fwrite(treat_threshold_df, "PDS_treatment_thresholds.csv")

# ================= Leakage map (drop definers from X) =================
leakage_drop_for <- function(tname) {
  drops <- switch(tname,
                  seller_N_total_value_high    = "log_sellern_total_value",
                  seller_N_1_total_value_high  = "log_sellern_1_total_value",
                  buyer_N_total_value_high     = "log_buyern_total_value",
                  buyer_N_1_total_value_high   = "log_buyern_1_total_value",
                  character(0)
  )
  intersect(drops, colnames(X_all_dt))
}

rebuild_X_without <- function(drop_names) {
  keep_names <- setdiff(colnames(X_all_dt), drop_names)
  X_cov_dt <- X_all_dt[, ..keep_names]
  X_cov_dt <- drop_constant_cols(X_cov_dt)
  as_numeric_design(X_cov_dt)
}

# ================= Containers =================
diag_overlap_rows  <- list()
diag_balance_rows  <- list()
pds_rows           <- list()

# ================= Runner (PDS only) =================
run_pds_for_treatment <- function(W, tname, y) {
  drop_names <- leakage_drop_for(tname)
  if (length(drop_names)) {
    message(sprintf("  [info %s] Dropping leakage cols: %s", tname, paste(drop_names, collapse = ", ")))
  }
  X_cov <- rebuild_X_without(drop_names)
  if (!ncol(X_cov)) {
    message(sprintf("  [skip %s] No covariates left after leakage/constant-drop.", tname))
    return(NULL)
  }
  
  # Overlap trimming
  ov <- compute_overlap_keep(X_cov, W, y, alpha_seq = c(0.10, 0.05, 0.02))
  keep <- ov$keep; pHat_all <- ov$pHat; alpha_used <- ov$alpha
  if (!length(keep)) {
    message(sprintf("  [skip %s] insufficient overlap/sample.", tname))
    return(NULL)
  }
  n_used <- length(keep); treated_share <- mean(W[keep])
  
  # Overlap diagnostics
  prop_str_all  <- if (!is.null(pHat_all)) quant_str(pHat_all[is.finite(pHat_all)]) else NA_character_
  prop_str_keep <- if (!is.null(pHat_all)) quant_str(pHat_all[keep]) else NA_character_
  diag_overlap_rows[[tname]] <<- data.frame(
    treatment   = tname,
    alpha_used  = alpha_used,
    n_total     = ov$n_total,
    n_kept      = ov$n_kept,
    kept_pct    = ov$kept_pct,
    prop_quant_all  = prop_str_all,
    prop_quant_kept = prop_str_keep,
    stringsAsFactors = FALSE
  )
  
  # Balance diagnostics (SMD) pre vs post
  smd_raw  <- smd_vec(X_cov, W)
  smd_keep <- smd_vec(X_cov[keep, , drop=FALSE], W[keep])
  balance_tbl <- data.frame(
    treatment = tname,
    covariate = colnames(X_cov),
    SMD_raw   = as.numeric(smd_raw),
    SMD_kept  = as.numeric(smd_keep),
    flag_gt_0_1  = as.integer(abs(smd_keep) > 0.10),
    flag_gt_0_25 = as.integer(abs(smd_keep) > 0.25),
    stringsAsFactors = FALSE
  )
  diag_balance_rows[[tname]] <<- balance_tbl
  
  ### >>> NEW: Multicollinearity diagnostics on KEPT sample (per treatment) ----
  Xi_kept <- X_cov[keep, , drop = FALSE]
  
  # Condition number on standardized kept X
  Xi_std <- scale(as.matrix(Xi_kept)); Xi_std[!is.finite(Xi_std)] <- 0
  sv_i <- svd(Xi_std, nu = 0, nv = 0)
  cond_i <- if (length(sv_i$d)) max(sv_i$d, na.rm=TRUE) / pmax(min(sv_i$d, na.rm=TRUE), 1e-12) else NA_real_
  data.table::fwrite(
    data.frame(treatment = tname, condition_number_kept = cond_i),
    "PDS_multicollinearity_condition_number_by_treatment.csv",
    append = file.exists("PDS_multicollinearity_condition_number_by_treatment.csv")
  )
  
  # Top |corr| pairs on kept sample (sample columns if very wide)
  cols_i <- colnames(Xi_kept)
  .max_cols_for_corr_i <- 1200
  if (length(cols_i) > .max_cols_for_corr_i) { set.seed(87); cols_i <- sample(cols_i, .max_cols_for_corr_i) }
  Xs <- as.matrix(Xi_kept[, cols_i, drop = FALSE])
  cr_i <- suppressWarnings(cor(Xs, use = "pairwise.complete.obs"))
  cr_i[upper.tri(cr_i, diag = TRUE)] <- NA
  thr_i <- 0.90
  hits_i <- which(abs(cr_i) >= thr_i, arr.ind = TRUE)
  if (nrow(hits_i)) {
    out_i <- data.frame(
      treatment = tname,
      var1 = rownames(cr_i)[hits_i[, "row"]],
      var2 = colnames(cr_i)[hits_i[, "col"]],
      r    = cr_i[hits_i],
      stringsAsFactors = FALSE
    )
    out_i <- out_i[order(-abs(out_i$r)), ]
    data.table::fwrite(out_i,
                       "PDS_multicollinearity_top_pairs_by_treatment.csv",
                       append = file.exists("PDS_multicollinearity_top_pairs_by_treatment.csv")
    )
  }
  ### <<< NEW ------------------------------------------------------------------
  
  # ---------- PDS on kept sample ----------
  Xk <- X_cov[keep, , drop=FALSE]; yk <- y[keep]; Wk <- W[keep]
  pds_fit <- hdm::rlassoEffect(x = Xk, y = yk, d = Wk, method = "double selection")
  sm_pds  <- summary(pds_fit)
  coef_mat <- tryCatch(sm_pds$coefficients, error = function(e) NULL)
  if (is.null(coef_mat)) coef_mat <- tryCatch(as.matrix(sm_pds), error = function(e) NULL)
  
  est <- as.numeric(coef_mat[1, 1]); se <- as.numeric(coef_mat[1, 2])
  tval <- as.numeric(coef_mat[1, 3])
  pval <- suppressWarnings(as.numeric(coef_mat[1, 4]))
  if (!is.finite(pval)) {
    df_pds <- max(nrow(Xk) - ncol(Xk) - 1, 1)
    pval <- ifelse(is.finite(tval), 2*pt(-abs(tval), df = df_pds), NA_real_)
  }
  lo <- est - 1.96 * se
  hi <- est + 1.96 * se
  
  row_pds <- data.frame(
    treatment       = tname,
    model           = "PDS",
    n_used          = n_used,
    treated_share   = treated_share,
    overlap_alpha   = alpha_used,
    ATE_log         = est, SE_log = se,
    CI95_lo_log     = lo,  CI95_hi_log = hi,
    z_or_t          = tval,   p_value = pval,
    ATE_pp          = 100*est, CI95_lo_pp = 100*lo, CI95_hi_pp = 100*hi,
    ATE_pct         = 100*(exp(est)-1),
    CI95_lo_pct     = 100*(exp(lo)-1), CI95_hi_pct = 100*(exp(hi)-1),
    stringsAsFactors = FALSE
  )
  pds_rows[[tname]] <<- row_pds
  
  cat(sprintf("  [PDS] %s  ATE(log)=%+.4f (%+.4f,%+.4f), t=%.2f, p=%.3g  → %+.2f%%\n",
              tname, est, lo, hi, tval, pval, 100*(exp(est)-1)))
  invisible(TRUE)
}

# ================= Run =================
for (nm in names(bin_treatments)) {
  W <- bin_treatments[[nm]]
  cat(sprintf("\n[Run] %s ...\n", nm))
  if (is.null(W) || length(unique(W[!is.na(W)])) < 2) {
    cat(sprintf("  → skipped %s (no variation).\n", nm)); next
  }
  run_pds_for_treatment(W, nm, y)
}

# ================= Save core outputs + dashboard =================
if (length(pds_rows)) {
  out_pds <- dplyr::bind_rows(pds_rows)
  out_pds$q_value_bh <- p.adjust(out_pds$p_value, method = "BH")
  data.table::fwrite(out_pds, "PDS_results_PDS.csv")
  
  if (length(diag_overlap_rows)) {
    data.table::fwrite(dplyr::bind_rows(diag_overlap_rows), "PDS_overlap_diagnostics.csv")
  }
  if (length(diag_balance_rows)) {
    data.table::fwrite(dplyr::bind_rows(diag_balance_rows), "PDS_balance_diagnostics.csv")
  }
  
  # Always save thresholds too
  data.table::fwrite(treat_threshold_df, "PDS_treatment_thresholds.csv")
  
  # ---------- Dashboard (PDS) ----------
  # Balance summary (robust to empty)
  if (length(diag_balance_rows)) {
    bal <- dplyr::bind_rows(diag_balance_rows)
    bal$abs_smd_kept <- abs(bal$SMD_kept)
    smd_summary <- bal %>%
      dplyr::group_by(treatment) %>%
      dplyr::summarise(
        total_covariates = sum(!is.na(abs_smd_kept)),
        count_gt_0_10    = sum(abs_smd_kept > 0.10, na.rm = TRUE),
        count_gt_0_25    = sum(abs_smd_kept > 0.25, na.rm = TRUE),
        median_abs_smd   = if (sum(!is.na(abs_smd_kept)) == 0) NA_real_ else median(abs_smd_kept, na.rm = TRUE),
        p90_abs_smd      = if (sum(!is.na(abs_smd_kept)) == 0) NA_real_ else as.numeric(quantile(abs_smd_kept, 0.90, na.rm = TRUE, names = FALSE)),
        max_abs_smd      = if (sum(!is.na(abs_smd_kept)) == 0) NA_real_ else max(abs_smd_kept, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    smd_summary <- tibble::tibble(
      treatment        = out_pds$treatment %>% unique(),
      total_covariates = NA_integer_,
      count_gt_0_10    = NA_integer_,
      count_gt_0_25    = NA_integer_,
      median_abs_smd   = NA_real_,
      p90_abs_smd      = NA_real_,
      max_abs_smd      = NA_real_
    )
  }
  
  # Overlap kept% summary
  if (length(diag_overlap_rows)) {
    ov_keep <- dplyr::bind_rows(diag_overlap_rows) %>% dplyr::select(treatment, kept_pct)
  } else {
    ov_keep <- tibble::tibble(treatment = out_pds$treatment %>% unique(), kept_pct = NA_real_)
  }
  
  # thresholds for join
  thresholds_for_join <- treat_threshold_df %>% dplyr::select(treatment, split_rule, cutoff_raw)
  
  dash <- out_pds %>%
    dplyr::left_join(ov_keep,             by = "treatment") %>%
    dplyr::left_join(smd_summary,         by = "treatment") %>%
    dplyr::left_join(thresholds_for_join, by = "treatment") %>%
    dplyr::mutate(
      flag_signif = dplyr::case_when(
        q_value_bh < 0.05 ~ "***",
        q_value_bh < 0.10 ~ "**",
        p_value   < 0.05  ~ "*",
        TRUE ~ ""
      ),
      flag_overlap = ifelse(!is.na(kept_pct) & kept_pct < 50, "poor_overlap", ""),
      flag_balance = ifelse(
        (!is.na(median_abs_smd) & median_abs_smd > 0.10) |
          (!is.na(p90_abs_smd)   & p90_abs_smd   > 0.25), "poor_balance", ""
      ),
      flag_tiny    = ifelse(!is.na(ATE_pct) & abs(ATE_pct) < 1, "tiny", ""),
      flags = paste0(
        flag_signif,
        ifelse(flag_overlap != "", paste0(",", flag_overlap), ""),
        ifelse(flag_balance != "", paste0(",", flag_balance), ""),
        ifelse(flag_tiny    != "", paste0(",", flag_tiny), "")
      ) %>% sub("^,", "", .)
    ) %>%
    dplyr::select(
      model, treatment,
      ATE_pct,
      cutoff_raw, split_rule,
      CI95_lo_pct, CI95_hi_pct,
      p_value, q_value_bh,
      kept_pct,
      total_covariates, count_gt_0_10, count_gt_0_25,
      median_abs_smd, p90_abs_smd, max_abs_smd,
      flags
    )
  
  data.table::fwrite(dash, "PDS_dashboard_PDS.csv")
  
  ### >>> NEW: simple combined multicollinearity summary -----------------------
  # Combine global condition number with per-treatment ones (if present)
  if (file.exists("PDS_multicollinearity_condition_number_by_treatment.csv")) {
    cond_by_t <- fread("PDS_multicollinearity_condition_number_by_treatment.csv")
    glob_cond <- fread("PDS_multicollinearity_global_condition_number.csv")
    glob_cond$treatment <- "GLOBAL_X_POOL"
    setnames(glob_cond, "condition_number", "condition_number_kept")
    all_cond <- rbindlist(list(glob_cond[, .(treatment, condition_number_kept)],
                               cond_by_t), use.names = TRUE, fill = TRUE)
    fwrite(all_cond, "PDS_multicollinearity_summary.csv")
  }
  ### <<< NEW ------------------------------------------------------------------
  
  cat("\n[Saved]\n",
      "  - PDS_results_PDS.csv\n",
      "  - PDS_overlap_diagnostics.csv\n",
      "  - PDS_balance_diagnostics.csv\n",
      "  - PDS_treatment_thresholds.csv\n",
      "  - PDS_dashboard_PDS.csv\n",
      "  - PDS_multicollinearity_global_condition_number.csv\n",
      "  - PDS_multicollinearity_global_top_pairs.csv\n",
      "  - PDS_multicollinearity_condition_number_by_treatment.csv\n",
      "  - PDS_multicollinearity_top_pairs_by_treatment.csv\n",
      "  - PDS_multicollinearity_summary.csv\n",
      sep = "")
} else {
  cat("\n[Note] No PDS results produced; check overlap trimming & treatment variation.\n")
}
