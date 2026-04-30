suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(fastDummies)
  library(grf)
  library(hdm)
  # robust-OLS fallback for BLP:
  library(sandwich)
  library(lmtest)
})

set.seed(87)

# ----------------------- Helpers -----------------------
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

to_binary_with_meta <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  v[!is.finite(v)] <- NA_real_
  if (all(is.na(v))) {
    return(list(w = rep(NA_integer_, length(v)),
                rule = NA_character_, cutoff_log = NA_real_,
                cutoff_raw = NA_real_, share_high = NA_real_))
  }
  med <- stats::median(v, na.rm = TRUE)
  w   <- as.integer(v >= med)
  if (length(unique(w[!is.na(w)])) == 2) {
    return(list(w = w, rule = "median",
                cutoff_log = med, cutoff_raw = expm1(med),
                share_high = mean(w == 1, na.rm = TRUE)))
  }
  q75 <- suppressWarnings(as.numeric(quantile(v, 0.75, na.rm = TRUE)))
  if (is.finite(q75)) {
    w <- as.integer(v >= q75)
    if (length(unique(w[!is.na(w)])) == 2) {
      return(list(w = w, rule = "q75",
                  cutoff_log = q75, cutoff_raw = expm1(q75),
                  share_high = mean(w == 1, na.rm = TRUE)))
    }
  }
  w <- as.integer(v > 0)
  if (length(unique(w[!is.na(w)])) == 2) {
    return(list(w = w, rule = ">0",
                cutoff_log = 0, cutoff_raw = 0,
                share_high = mean(w == 1, na.rm = TRUE)))
  }
  list(w = rep(NA_integer_, length(v)),
       rule = NA_character_, cutoff_log = NA_real_,
       cutoff_raw = NA_real_, share_high = NA_real_)
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

# ---- Balance & summary helpers ----
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

# ---- Safe calibration extractor ----
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
        out$slope  <- unname(co[slope_row[1], "Estimate"])
        pcol <- intersect(colnames(co), c("Pr(>|t|)", "Pr(>|z|)", "p.value"))
        if (length(pcol)) out$p_slope <- unname(co[slope_row[1], pcol[1]])
      }
    }
    return(out)
  }
  if (is.list(cal) && !is.null(cal$regression) && inherits(cal$regression, "lm")) {
    co <- tryCatch(coef(summary(cal$regression)), error = function(e) NULL)
    if (!is.null(co)) {
      if ("(Intercept)" %in% rownames(co)) out$intercept <- unname(co["(Intercept)", "Estimate"])
      slope_row <- rownames(co)[grepl("pred", rownames(co), ignore.case = TRUE)]
      if (length(slope_row)) {
        out$slope  <- unname(co[slope_row[1], "Estimate"])
        pcol <- intersect(colnames(co), c("Pr(>|t|)", "Pr(>|z|)", "p.value"))
        if (length(pcol)) out$p_slope <- unname(co[slope_row[1], pcol[1]])
      }
    }
    return(out)
  }
  if (is.list(cal) && !is.null(cal$coefficients)) {
    co <- as.matrix(cal$coefficients); rn <- rownames(co)
    if (!is.null(rn)) {
      if ("(Intercept)" %in% rn) out$intercept <- co["(Intercept)", 1]
      slope_row <- rn[grepl("pred", rn, ignore.case = TRUE)]
      if (length(slope_row)) {
        out$slope <- co[slope_row[1], 1]
        if (ncol(co) >= 4) out$p_slope <- co[slope_row[1], 4]
      }
    }
    return(out)
  }
  if (is.atomic(cal) && is.numeric(cal) && length(cal) == 1L && is.finite(cal)) {
    out$p_slope <- as.numeric(cal)
    return(out)
  }
  out
}

# ---- Overlap helper (parametrized) ----
compute_overlap_keep <- function(X_cov, W, y,
                                 alpha_seq = c(0.10, 0.05, 0.02),
                                 n_min = 200) {
  base_idx <- which(is.finite(y) & !is.na(W))
  Wbin <- as.integer(W > 0)
  out <- list(keep = integer(0), pHat = NULL, alpha = NA_real_,
              n_total = length(base_idx), n_kept = 0, kept_pct = 0)
  if (length(unique(Wbin[base_idx])) < 2 || length(base_idx) < n_min) return(out)
  pHat_base <- tryCatch({
    pf <- probability_forest(X_cov[base_idx, , drop=FALSE], Wbin[base_idx],
                             num.trees = 1000, honesty = TRUE)
    predict(pf)$predictions
  }, error = function(e) {
    rf <- regression_forest(X_cov[base_idx, , drop=FALSE], Wbin[base_idx],
                            num.trees = 1000, honesty = TRUE)
    predict(rf)$predictions
  })
  pHat_base <- pmin(pmax(pHat_base, 1e-6), 1-1e-6)
  for (a in alpha_seq) {
    ok <- if (a <= 0) seq_along(pHat_base) else which(pHat_base >= a & pHat_base <= (1 - a))
    keep <- base_idx[ok]
    if (length(keep) >= n_min && length(unique(Wbin[keep])) == 2) {
      pHat_all <- rep(NA_real_, length(W)); pHat_all[base_idx] <- pHat_base
      out$keep <- keep; out$pHat <- pHat_all; out$alpha <- a
      out$n_kept <- length(keep); out$kept_pct <- 100*length(keep)/length(base_idx)
      return(out)
    }
  }
  out
}

# ---- CV R2_tau ----
compute_r2tau_cv <- function(X, y, K = 5, seed = 87) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X)
  if (!is.finite(n) || n < max(50, K + 5)) {
    return(c(r2_tau = NA_real_, r2_tau_oos = NA_real_))
  }
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
  vY    <- stats::var(y,    na.rm = TRUE)
  vYhat <- stats::var(yhat, na.rm = TRUE)
  mse   <- mean((y - yhat)^2, na.rm = TRUE)
  r2_var <- if (is.finite(vY) && vY > 0) vYhat / vY else NA_real_
  r2_oos <- if (is.finite(vY) && vY > 0) 1 - (mse / vY) else NA_real_
  c(r2_tau = r2_var, r2_tau_oos = r2_oos)
}

# ----------------------- Load & base filter (Condition B) -----------------------
df1 <- fread(
  "df_table1.csv",
  select = c("token_id","time_n_sale","time_n-1_sale","time_n-2_sale",
             "price_n_sale","price_n-1_sale","price_n-2_sale","price_1_sale",
             "buyer_n_sale","seller_n_sale","buyer_n-1_sale","seller_n-1_sale")
)
df1 <- df1[!is.na(`time_n-1_sale`) & !is.na(`price_n-1_sale`) & (`price_n-1_sale` > 0)]
# Condition B: require mismatch
df1 <- df1[`buyer_n-1_sale` != `seller_n_sale`]

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
                                 "median_offer_price","duration_offer_days"))

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

# drop tscount / act_period (sparse & collinear with total_value)
dt[, c("buyern_tscount","buyern_act_period",
       "sellern_tscount","sellern_act_period",
       "buyern_1_tscount","buyern_1_act_period",
       "sellern_1_tscount","sellern_1_act_period") := NULL]

# Row id
dt[, row_id := .I]

# Time FE strings
dt[, time_n_sale_str     := format(as.POSIXct(time_n_sale,     origin="1970-01-01", tz="UTC"), "%Y-%m")]
dt[, `time_n-1_sale_str` := format(as.POSIXct(`time_n-1_sale`, origin="1970-01-01", tz="UTC"), "%Y-%m")]

# Never allow raw time columns in X
cols_never_in_X <- c("time_n_sale","time_n-1_sale","time_n-2_sale",
                     "time_n_sale_dt","time_n-1_sale_dt","rarity.rank","time_of_sale")
dt[, (intersect(cols_never_in_X, names(dt))) := NULL]

# ----------------------- Gas features (ETH) then drop raw gas cols -----------------------
WEI_TO_ETH <- 1/1e18
gas_cost_eth_from_gwei <- function(units, price_wei) {
  u <- suppressWarnings(as.numeric(units))
  p <- suppressWarnings(as.numeric(price_wei))
  ifelse(is.finite(u) & is.finite(p), u * p * WEI_TO_ETH, NA_real_)
}

# Gas paid (gasUsed * gasPrice) in ETH
dt[, buyern_avg_total_gaspaid    := gas_cost_eth_from_gwei(buyern_total_gasUsed,    buyern_avg_gasPrice)]
dt[, buyern_1_avg_total_gaspaid  := gas_cost_eth_from_gwei(buyern_1_total_gasUsed,  buyern_1_avg_gasPrice)]
dt[, sellern_avg_total_gaspaid   := gas_cost_eth_from_gwei(sellern_total_gasUsed,   sellern_avg_gasPrice)]
dt[, sellern_1_avg_total_gaspaid := gas_cost_eth_from_gwei(sellern_1_total_gasUsed, sellern_1_avg_gasPrice)]

# Gas limit cost in ETH
dt[, buyern_avg_gasLimit_cost    := gas_cost_eth_from_gwei(buyern_avg_gasLimit,    buyern_avg_gasPrice)]
dt[, buyern_1_avg_gasLimit_cost  := gas_cost_eth_from_gwei(buyern_1_avg_gasLimit,  buyern_1_avg_gasPrice)]
dt[, sellern_avg_gasLimit_cost   := gas_cost_eth_from_gwei(sellern_avg_gasLimit,   sellern_avg_gasPrice)]
dt[, sellern_1_avg_gasLimit_cost := gas_cost_eth_from_gwei(sellern_1_avg_gasLimit, sellern_1_avg_gasPrice)]

# Drop raw gas inputs
dt[, c("buyern_total_gasUsed","buyern_avg_gasPrice",
       "buyern_1_total_gasUsed","buyern_1_avg_gasPrice",
       "sellern_total_gasUsed","sellern_avg_gasPrice",
       "sellern_1_total_gasUsed","sellern_1_avg_gasPrice",
       "buyern_avg_gasLimit","buyern_1_avg_gasLimit",
       "sellern_avg_gasLimit","sellern_1_avg_gasLimit") := NULL]

# ----------------------- Trait recoding + dummies -----------------------
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

cat_cols <- c("time_n_sale_str","time_n-1_sale_str", paste0(traits, "_recoded"))
dt <- fastDummies::dummy_cols(
  dt,
  select_columns = cat_cols,
  remove_first_dummy = TRUE,
  remove_selected_columns = TRUE
)
dt[, c("Background","Clothes","Earring","Eyes","Fur","Hat","Mouth") := NULL]

# ----------------------- Winsorize + log1p -----------------------
cols_to_w <- c(
  "price_n_sale","price_n-1_sale","price_n-2_sale","price_1_sale",
  "buyern_total_value","buyern_rolling_avg_value_last10","buyern_rolling_std_value_last10",
  "sellern_total_value","sellern_rolling_avg_value_last10","sellern_rolling_std_value_last10",
  "buyern_1_total_value","buyern_1_rolling_avg_value_last10","buyern_1_rolling_std_value_last10",
  "sellern_1_total_value","sellern_1_rolling_avg_value_last10","sellern_1_rolling_std_value_last10",
  "median_offer_price","duration_offer_days","total_offers","unique_makers_count",
  "buyern_avg_total_gaspaid","buyern_1_avg_total_gaspaid",
  "sellern_avg_total_gaspaid","sellern_1_avg_total_gaspaid",
  "buyern_avg_gasLimit_cost","buyern_1_avg_gasLimit_cost",
  "sellern_avg_gasLimit_cost","sellern_1_avg_gasLimit_cost"
)
dt <- winsorize_log1p_cols(dt, cols_to_w)

# Outcomes
dt[, log_price_change := `log_price_n_sale` - `log_price_n-1_sale`]
if ("log_price_n-2_sale" %in% names(dt)) {
  dt[, "log_price_change_n-1" := `log_price_n-1_sale` - `log_price_n-2_sale`]
} else {
  dt[, "log_price_change_n-1" := 0]
}

# Replace remaining NA numerics with 0
num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
for (nm in num_cols) data.table::set(dt, i = which(is.na(dt[[nm]])), j = nm, value = 0)

# ----------------------- Design matrix (X) -----------------------
drop_from_X <- intersect(
  c("token_id","row_id","buyer_n_sale","seller_n_sale","buyer_n-1_sale","seller_n-1_sale",
    "log_price_change","log_price_n_sale","log_price_n-1_sale","log_price_n-2_sale","days_since_prev",
    "log_sellern_rolling_std_value_last10","log_buyern_1_rolling_std_value_last10",
    "log_total_offers","log_duration_offer_days","log_buyern_rolling_std_value_last10",
    "log_sellern_1_rolling_std_value_last10"),
  names(dt)
)
X_all_dt <- dt[, setdiff(names(dt), drop_from_X), with = FALSE]
X_all_dt <- drop_constant_cols(X_all_dt)
X_all <- as_numeric_design(X_all_dt)

y <- dt$log_price_change

# ----------------------- Treatments (binary) -----------------------
treat_defs <- list(
  seller_N_total_value_high    = dt$log_sellern_total_value,
  seller_N_1_total_value_high  = dt$log_sellern_1_total_value,
  buyer_N_total_value_high     = dt$log_buyern_total_value,
  buyer_N_1_total_value_high   = dt$log_buyern_1_total_value,
  buyer_N_gaspaid_high         = dt$log_buyern_avg_total_gaspaid,
  buyer_N_1_gaspaid_high       = dt$log_buyern_1_avg_total_gaspaid,
  seller_N_gaspaid_high        = dt$log_sellern_avg_total_gaspaid,
  seller_N_1_gaspaid_high      = dt$log_sellern_1_avg_total_gaspaid,
  buyer_N_gaslimit_cost_high   = dt$log_buyern_avg_gasLimit_cost,
  buyer_N_1_gaslimit_cost_high = dt$log_buyern_1_avg_gasLimit_cost,
  seller_N_gaslimit_cost_high  = dt$log_sellern_avg_gasLimit_cost,
  seller_N_1_gaslimit_cost_high= dt$log_sellern_1_avg_gasLimit_cost
)

bin_treatments <- list()
treat_threshold_meta <- list()

for (nm in names(treat_defs)) {
  v <- treat_defs[[nm]]; if (is.null(v)) next
  res <- to_binary_with_meta(v)
  w <- res$w
  if (all(is.na(w)) || length(unique(w[!is.na(w)])) < 2) {
    message(sprintf("[treat] %s → no variation after binarization; skipping.", nm))
    next
  }
  bin_treatments[[nm]] <- w
  treat_threshold_meta[[nm]] <- data.frame(
    treatment   = nm,
    split_rule  = res$rule,
    cutoff_log  = res$cutoff_log,
    cutoff_raw  = res$cutoff_raw,
    share_high0 = res$share_high,
    stringsAsFactors = FALSE
  )
  cat(sprintf("[treat] %-32s rule=%-6s cutoff_log=% .3f cutoff_raw=% .3f share_high=%.1f%%\n",
              nm, res$rule, res$cutoff_log, res$cutoff_raw, 100*res$share_high))
}
if (!length(bin_treatments)) stop("No usable treatments after binarization. Check inputs/variation.")
treat_threshold_df <- dplyr::bind_rows(treat_threshold_meta)
data.table::fwrite(treat_threshold_df, "cf_M3_treatment_thresholds.csv")

# ----------------------- Combinations & Strategies — Condition B -----------------------
# Helpers
make_11vs00 <- function(A, B) ifelse(A==1 & B==1, 1L, ifelse(A==0 & B==0, 0L, NA_integer_))
make_multi_11vs00 <- function(...) {
  L <- list(...)
  all1 <- Reduce(`&`, lapply(L, function(v) v == 1L))
  all0 <- Reduce(`&`, lapply(L, function(v) v == 0L))
  ifelse(all1, 1L, ifelse(all0, 0L, NA_integer_))
}
inv01  <- function(x) ifelse(is.na(x), NA_integer_, 1L - as.integer(x))

# Shorthand to safely add only when all inputs exist
.safe_add <- function(name, ...) {
  parts <- list(...)
  if (all(vapply(parts, function(z) !is.null(z) && is.integer(z), logical(1)))) {
    bin_treatments[[name]] <<- make_multi_11vs00(...)
  } else {
    message(sprintf("[combo] Skipped %s (some inputs missing / no variation).", name))
  }
}

# Clear any previously-defined combos/strategies to keep ONLY what we want for Condition B
rm_names <- grep(
  paste(c("_Tvalue_.*_11vs00$", "_Gas(limit|paid)_given_", "^all_total_value_11vs00$",
          "^buyer_total_value_11vs00$", "^seller_total_value_11vs00$", "^combo"),
        collapse="|"),
  names(bin_treatments), value = TRUE
)
if (length(rm_names)) bin_treatments[rm_names] <- NULL

# Singles shorthand (from bin_treatments created above)
bnTVH  <- bin_treatments[["buyer_N_total_value_high"]]
bn1TVH <- bin_treatments[["buyer_N_1_total_value_high"]]
snTVH  <- bin_treatments[["seller_N_total_value_high"]]
sn1TVH <- bin_treatments[["seller_N_1_total_value_high"]]

bnGPH  <- bin_treatments[["buyer_N_gaspaid_high"]]
bn1GPH <- bin_treatments[["buyer_N_1_gaspaid_high"]]
snGPH  <- bin_treatments[["seller_N_gaspaid_high"]]
sn1GPH <- bin_treatments[["seller_N_1_gaspaid_high"]]

bnGLH  <- bin_treatments[["buyer_N_gaslimit_cost_high"]]
bn1GLH <- bin_treatments[["buyer_N_1_gaslimit_cost_high"]]
snGLH  <- bin_treatments[["seller_N_gaslimit_cost_high"]]
sn1GLH <- bin_treatments[["seller_N_1_gaslimit_cost_high"]]

# Needed “low” versions (locals only)
bn1TVL <- inv01(bn1TVH)
snTVL  <- inv01(snTVH)
sn1TVL <- inv01(sn1TVH)

# ----------------------- COMBINATIONS -----------------------
# Buyer N High: only have (Total Value + Gas Paid) and (Total Value + Gas Limit Cost)
.safe_add("buyer_N_Tvalue_Gaspaid_11vs00",  bnTVH, bnGPH)
.safe_add("buyer_N_Tvalue_Gaslimit_11vs00", bnTVH, bnGLH)

# Buyer N-1 High:
#   (1) Total value + Gas Paid
#   (2) Total value + Gas Limit Cost
#   (3) Total value + Gas Paid + Gas Limit Cost (3-way)
.safe_add("buyer_N_1_Tvalue_Gaspaid_11vs00",               bn1TVH, bn1GPH)
.safe_add("buyer_N_1_Tvalue_Gaslimit_11vs00",              bn1TVH, bn1GLH)
.safe_add("buyer_N_1_Tvalue_Gaspaid_Gaslimit_11vs00",      bn1TVH, bn1GPH, bn1GLH)

# Seller N High: Total Value + Gas Limit Cost only
.safe_add("seller_N_Tvalue_Gaslimit_11vs00",               snTVH,  snGLH)

# Seller N-1 High: Total Value + Gas Paid only
.safe_add("seller_N_1_Tvalue_Gaspaid_11vs00",              sn1TVH, sn1GPH)

# Combination All:
#   Buyer N-1 TV High + Seller N TV High + Seller N-1 TV High  (NOTE: excludes Buyer N TV)
.safe_add("combo_all_TV_BN1_SN_SN1_11vs00",                 bn1TVH, snTVH, sn1TVH)

# Seller only:
#   Seller N TV High + Seller N-1 TV High
.safe_add("seller_total_value_11vs00",                       snTVH,  sn1TVH)

# ----------------------- STRATEGIES (prefix 'combo' → no trimming) -----------------------
# (1) BN-1 low TV + SN low TV + SN-1 low TV
.safe_add("comboB_lowTV_all_11vs00",                         bn1TVL, snTVL, sn1TVL)

# (2) BN-1 low TV + BN-1 high GL + BN-1 high GP
.safe_add("comboB_BN1_lowTV_highFees_11vs00",                bn1TVL, bn1GLH, bn1GPH)

# (3) SN low TV + SN high GP + SN high GL
.safe_add("comboB_SN_lowTV_highFees_11vs00",                 snTVL,  snGPH,  snGLH)

# (4) Mix of (2) and (3): BN-1 low TV + BN-1 high GL + BN-1 high GP + SN low TV + SN high GP + SN high GL
.safe_add("comboB_mix_lowTV_highFees_11vs00",                bn1TVL, bn1GLH, bn1GPH, snTVL, snGPH, snGLH)

# ----------------------- Leakage map -----------------------
leakage_drop_for <- function(tname) {
  map <- list(
    # singles
    seller_N_total_value_high      = "log_sellern_total_value",
    seller_N_1_total_value_high    = "log_sellern_1_total_value",
    buyer_N_total_value_high       = "log_buyern_total_value",
    buyer_N_1_total_value_high     = "log_buyern_1_total_value",
    buyer_N_gaspaid_high           = "log_buyern_avg_total_gaspaid",
    buyer_N_1_gaspaid_high         = "log_buyern_1_avg_total_gaspaid",
    seller_N_gaspaid_high          = "log_sellern_avg_total_gaspaid",
    seller_N_1_gaspaid_high        = "log_sellern_1_avg_total_gaspaid",
    buyer_N_gaslimit_cost_high     = "log_buyern_avg_gasLimit_cost",
    buyer_N_1_gaslimit_cost_high   = "log_buyern_1_avg_gasLimit_cost",
    seller_N_gaslimit_cost_high    = "log_sellern_avg_gasLimit_cost",
    seller_N_1_gaslimit_cost_high  = "log_sellern_1_avg_gasLimit_cost",
    
    # --- Condition B combos ---
    buyer_N_Tvalue_Gaspaid_11vs00              = c("log_buyern_total_value","log_buyern_avg_total_gaspaid"),
    buyer_N_Tvalue_Gaslimit_11vs00             = c("log_buyern_total_value","log_buyern_avg_gasLimit_cost"),
    
    buyer_N_1_Tvalue_Gaspaid_11vs00            = c("log_buyern_1_total_value","log_buyern_1_avg_total_gaspaid"),
    buyer_N_1_Tvalue_Gaslimit_11vs00           = c("log_buyern_1_total_value","log_buyern_1_avg_gasLimit_cost"),
    buyer_N_1_Tvalue_Gaspaid_Gaslimit_11vs00   = c("log_buyern_1_total_value",
                                                   "log_buyern_1_avg_total_gaspaid",
                                                   "log_buyern_1_avg_gasLimit_cost"),
    
    seller_N_Tvalue_Gaslimit_11vs00            = c("log_sellern_total_value","log_sellern_avg_gasLimit_cost"),
    seller_N_1_Tvalue_Gaspaid_11vs00           = c("log_sellern_1_total_value","log_sellern_1_avg_total_gaspaid"),
    
    combo_all_TV_BN1_SN_SN1_11vs00             = c("log_buyern_1_total_value",
                                                   "log_sellern_total_value",
                                                   "log_sellern_1_total_value"),
    seller_total_value_11vs00                  = c("log_sellern_total_value","log_sellern_1_total_value"),
    
    # --- Strategies (drop all sources used, even when a leg is 'low') ---
    comboB_lowTV_all_11vs00                    = c("log_buyern_1_total_value",
                                                   "log_sellern_total_value",
                                                   "log_sellern_1_total_value"),
    comboB_BN1_lowTV_highFees_11vs00           = c("log_buyern_1_total_value",
                                                   "log_buyern_1_avg_total_gaspaid",
                                                   "log_buyern_1_avg_gasLimit_cost"),
    comboB_SN_lowTV_highFees_11vs00            = c("log_sellern_total_value",
                                                   "log_sellern_avg_total_gaspaid",
                                                   "log_sellern_avg_gasLimit_cost"),
    comboB_mix_lowTV_highFees_11vs00           = c("log_buyern_1_total_value",
                                                   "log_buyern_1_avg_total_gaspaid",
                                                   "log_buyern_1_avg_gasLimit_cost",
                                                   "log_sellern_total_value",
                                                   "log_sellern_avg_total_gaspaid",
                                                   "log_sellern_avg_gasLimit_cost")
  )
  drops <- map[[tname]]
  if (is.null(drops)) drops <- character(0)
  intersect(drops, colnames(X_all_dt))
}


rebuild_X_without <- function(drop_names) {
  keep_names <- setdiff(colnames(X_all_dt), drop_names)
  X_cov_dt <- X_all_dt[, ..keep_names]
  X_cov_dt <- drop_constant_cols(X_cov_dt)
  as_numeric_design(X_cov_dt)
}

# ----------------------- Diagnostics containers -----------------------
diag_overlap_rows  <- list()
diag_balance_rows  <- list()
diag_robust_rows   <- list()
diag_cate_cov_overview <- list()
diag_cate_cov_bins     <- list()
diag_blp_rows <- list()
diag_r2tau_rows <- list()

# ---- BLP helpers ----
safe_blp_single <- function(cf, x_col_matrix) {
  out <- tryCatch(grf::best_linear_projection(cf, x_col_matrix, compute.se = TRUE), error = function(e) NULL)
  if (!is.null(out)) return(out)
  out <- tryCatch(grf::best_linear_projection(cf, X = x_col_matrix, compute.se = TRUE), error = function(e) NULL)
  if (!is.null(out)) return(out)
  out <- tryCatch(grf::best_linear_projection(cf, x_col_matrix), error = function(e) NULL)
  if (!is.null(out)) return(out)
  out <- tryCatch(grf::best_linear_projection(cf, X = x_col_matrix), error = function(e) NULL)
  if (!is.null(out)) return(out)
  NULL
}
parse_blp_coef_se <- function(blp_obj) {
  if (is.null(blp_obj)) return(NULL)
  if (is.list(blp_obj)) {
    if (!is.null(blp_obj$coefficients) && !is.null(blp_obj$standard.errors)) {
      return(list(beta = as.numeric(blp_obj$coefficients)[2],
                  se   = as.numeric(blp_obj$standard.errors)[2]))
    }
    if (!is.null(blp_obj$beta) && !is.null(blp_obj$se)) {
      return(list(beta = as.numeric(blp_obj$beta)[2],
                  se   = as.numeric(blp_obj$se)[2]))
    }
    if (!is.null(blp_obj$gamma_hat) && !is.null(blp_obj$gamma_se)) {
      return(list(beta = as.numeric(blp_obj$gamma_hat)[2],
                  se   = as.numeric(blp_obj$gamma_se)[2]))
    }
  }
  if (is.matrix(blp_obj) || is.data.frame(blp_obj)) {
    est <- suppressWarnings(as.numeric(blp_obj[2, 1]))
    se  <- suppressWarnings(as.numeric(blp_obj[2, 2]))
    if (is.finite(est) && is.finite(se)) return(list(beta = est, se = se))
  }
  if (is.atomic(blp_obj) && is.numeric(blp_obj)) return(NULL)
  NULL
}

# ----------------------- Runner -----------------------
run_models_for_treatment <- function(W, tname, y) {
  drop_names <- leakage_drop_for(tname)
  if (length(drop_names)) {
    message(sprintf("  [info %s] Dropping leakage cols: %s", tname, paste(drop_names, collapse = ", ")))
  }
  X_cov <- rebuild_X_without(drop_names)
  if (!ncol(X_cov)) {
    message(sprintf("  [skip %s] No covariates left after leakage/constant-drop.", tname))
    return(NULL)
  }
  
  # Overlap trimming (lighter for strategies)
  is_strategy <- grepl("^combo", tname)
  ov <- compute_overlap_keep(
    X_cov, W, y,
    alpha_seq = if (is_strategy) c(0.00) else c(0.05, 0.02, 0.00),
    n_min     = if (is_strategy) 50 else 200
  )
  keep <- ov$keep; pHat_all <- ov$pHat; alpha_used <- ov$alpha
  
  # Fallback if still empty but both groups present
  if (!length(keep)) {
    base_idx <- which(is.finite(y) & !is.na(W))
    if (length(unique(W[base_idx])) == 2 && length(base_idx) >= 50) {
      message(sprintf("  [fallback %s] using untrimmed sample (n=%d).", tname, length(base_idx)))
      keep       <- base_idx
      pHat_all   <- NULL
      alpha_used <- NA_real_
    } else {
      message(sprintf("  [skip %s] insufficient overlap/sample.", tname))
      return(NULL)
    }
  }
  
  n_used <- length(keep); treated_share <- mean(W[keep])
  
  # Overlap diagnostics
  n_total_base <- if (is.null(ov$n_total)) length(which(is.finite(y) & !is.na(W))) else ov$n_total
  kept_pct_calc <- if (n_total_base > 0) 100 * length(keep) / n_total_base else NA_real_
  diag_overlap_rows[[tname]] <<- data.frame(
    treatment       = tname,
    alpha_used      = alpha_used,
    n_total         = n_total_base,
    n_kept          = length(keep),
    kept_pct        = kept_pct_calc,
    prop_quant_all  = if (!is.null(pHat_all)) quant_str(pHat_all[is.finite(pHat_all)]) else NA_character_,
    prop_quant_kept = if (!is.null(pHat_all)) quant_str(pHat_all[keep])               else NA_character_,
    stringsAsFactors = FALSE
  )
  
  # GRF
  X_kept <- X_cov[keep, , drop=FALSE]
  of <- regression_forest(X_kept, y[keep], num.trees = 2000, honesty = TRUE)
  mu <- predict(of)$predictions
  cf <- causal_forest(
    X_kept, y[keep], W[keep],
    Y.hat = mu,
    W.hat = if (!is.null(pHat_all)) pHat_all[keep] else NULL,
    num.trees = 2000, min.node.size = 50, honesty = TRUE
  )
  
  target_primary <- if (!is.na(alpha_used) && alpha_used > 0) "overlap" else "all"
  ate <- average_treatment_effect(cf, target.sample = target_primary)
  est <- as.numeric(ate[["estimate"]]); se <- as.numeric(ate[["std.err"]])
  z   <- ifelse(se > 0, est / se, NA_real_)
  p   <- ifelse(is.finite(z), 2*pnorm(-abs(z)), NA_real_)
  lo  <- est - 1.96*se; hi <- est + 1.96*se
  
  cal_res <- safe_calibration_v2(cf)
  calib_int   <- cal_res$intercept
  calib_slope <- cal_res$slope
  calib_p     <- cal_res$p_slope
  
  ate_all_log <- NA_real_; ate_overlap_log <- NA_real_; ate_treated_log <- NA_real_
  if (!is.null(pHat_all)) {
    e_kept <- pmin(pmax(pHat_all[keep], 1e-6), 1-1e-6)
    extreme <- (min(e_kept) < 0.05) || (max(e_kept) > 0.95)
    if (!extreme) {
      ate_all_log     <- as.numeric(average_treatment_effect(cf, target.sample="all")["estimate"])
      ate_treated_log <- as.numeric(average_treatment_effect(cf, target.sample="treated")["estimate"])
    }
    ate_overlap_log <- as.numeric(average_treatment_effect(cf, target.sample="overlap")["estimate"])
  }
  
  # AIPW on kept (optional)
  aipw_est <- NA_real_; aipw_se <- NA_real_
  if (!is.null(pHat_all)) {
    y_kept <- y[keep]; W_kept <- W[keep]; e <- pmin(pmax(pHat_all[keep], 1e-6), 1-1e-6)
    idx1 <- which(W_kept == 1L); idx0 <- which(W_kept == 0L)
    rf1 <- regression_forest(X_kept[idx1, , drop=FALSE], y_kept[idx1], num.trees=1000, honesty=TRUE)
    rf0 <- regression_forest(X_kept[idx0, , drop=FALSE], y_kept[idx0], num.trees=1000, honesty=TRUE)
    m1  <- predict(rf1, X_kept)$predictions; m0 <- predict(rf0, X_kept)$predictions
    resY <- y_kept - (W_kept*m1 + (1-W_kept)*m0)
    dr   <- (W_kept - e) * resY / (e*(1-e)) + (m1 - m0)
    aipw_est <- mean(dr); aipw_se <- sd(dr)/sqrt(length(dr))
  }
  
  # VI (best-effort)
  vi  <- try(variable_importance(cf), silent = TRUE)
  ord <- if (inherits(vi, "try-error") || all(!is.finite(vi))) integer(0) else order(vi, decreasing = TRUE)
  top_imp <- if (!length(ord)) "(n/a)" else {
    k <- min(8, length(ord)); paste0(colnames(X_kept)[ord][seq_len(k)], ":", sprintf("%.4f", vi[ord][seq_len(k)]), collapse = ", ")
  }
  
  # CATEs (kept)
  pred <- predict(cf, estimate.variance = TRUE)
  cate_log <- as.numeric(pred$predictions)
  cate_se  <- sqrt(pmax(0, pred$variance.estimates))
  cate_lo  <- cate_log - 1.96 * cate_se
  cate_hi  <- cate_log + 1.96 * cate_se
  cate_pct <- 100 * (exp(cate_log) - 1)
  cate_df <- data.frame(
    treatment   = tname,
    row_id      = dt$row_id[keep],
    token_id    = dt$token_id[keep],
    cate_log    = cate_log,
    cate_se_log = cate_se,
    cate_lo_log = cate_lo,
    cate_hi_log = cate_hi,
    cate_pp     = 100 * cate_log,
    cate_pct    = cate_pct,
    stringsAsFactors = FALSE
  )
  
  # Per-treatment R2_tau
  r2pair <- compute_r2tau_cv(X_kept, cate_log, K = 5, seed = 87)
  diag_r2tau_rows[[tname]] <<- data.frame(
    treatment   = tname,
    r2_tau      = as.numeric(r2pair["r2_tau"]),
    r2_tau_oos  = as.numeric(r2pair["r2_tau_oos"]),
    n_kept      = length(keep),
    stringsAsFactors = FALSE
  )
  
  # ATE rows
  row_grf <- data.frame(
    treatment       = tname,
    model           = "GRF",
    n_used          = length(keep),
    treated_share   = treated_share,
    overlap_alpha   = alpha_used,
    target_primary  = target_primary,
    ATE_log         = est, SE_log = se,
    CI95_lo_log     = lo,  CI95_hi_log = hi,
    z_or_t          = z,   p_value     = p,
    ATE_pp          = 100*est, CI95_lo_pp = 100*lo, CI95_hi_pp = 100*hi,
    ATE_pct         = 100*(exp(est)-1),
    CI95_lo_pct     = 100*(exp(lo)-1), CI95_hi_pct = 100*(exp(hi)-1),
    top_importances = top_imp,
    calib_intercept = calib_int,
    calib_slope     = calib_slope,
    calib_p_slope   = calib_p,
    AIPW_log        = aipw_est,
    AIPW_se_log     = aipw_se,
    stringsAsFactors = FALSE
  )
  
  # PDS
  Xk <- X_kept; yk <- y[keep]; Wk <- W[keep]
  pds_fit <- hdm::rlassoEffect(x = Xk, y = yk, d = Wk, method = "double selection")
  sm_pds  <- summary(pds_fit)
  coef_mat <- tryCatch(sm_pds$coefficients, error = function(e) NULL)
  if (is.null(coef_mat)) coef_mat <- tryCatch(as.matrix(sm_pds), error = function(e) NULL)
  est_pds <- as.numeric(coef_mat[1, 1]); se_pds <- as.numeric(coef_mat[1, 2])
  t_pds   <- as.numeric(coef_mat[1, 3])
  p_pds   <- suppressWarnings(as.numeric(coef_mat[1, 4]))
  if (!is.finite(p_pds)) {
    df_pds <- max(nrow(Xk) - ncol(Xk) - 1, 1); p_pds <- ifelse(is.finite(t_pds), 2*pt(-abs(t_pds), df=df_pds), NA_real_)
  }
  ci_pds <- tryCatch(confint(pds_fit, level = 0.95), error = function(e) NULL)
  if (!is.null(ci_pds)) {
    if (is.matrix(ci_pds)) { lo_pds <- as.numeric(ci_pds[1, 1]); hi_pds <- as.numeric(ci_pds[1, 2]) }
    else { lo_pds <- as.numeric(ci_pds[1]); hi_pds <- as.numeric(ci_pds[2]) }
  } else { lo_pds <- est_pds - 1.96 * se_pds; hi_pds <- est_pds + 1.96 * se_pds }
  row_pds <- data.frame(
    treatment       = tname,
    model           = "PDS",
    n_used          = length(keep),
    treated_share   = treated_share,
    overlap_alpha   = alpha_used,
    target_primary  = NA_character_,
    ATE_log         = est_pds, SE_log = se_pds,
    CI95_lo_log     = lo_pds,  CI95_hi_log = hi_pds,
    z_or_t          = t_pds,   p_value     = p_pds,
    ATE_pp          = 100*est_pds, CI95_lo_pp = 100*lo_pds, CI95_hi_pp = 100*hi_pds,
    ATE_pct         = 100*(exp(est_pds)-1),
    CI95_lo_pct     = 100*(exp(lo_pds)-1), CI95_hi_pct = 100*(exp(hi_pds)-1),
    top_importances = NA_character_,
    calib_intercept = NA_real_,
    calib_slope     = NA_real_,
    calib_p_slope   = NA_real_,
    AIPW_log        = NA_real_,
    AIPW_se_log     = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # ---- CATE vs covariate (kept) ----
  cov_names <- colnames(X_kept)
  for (cv in cov_names) {
    v <- as.numeric(X_kept[, cv]); idx <- which(is.finite(v) & is.finite(cate_pct))
    if (length(idx) < 30) next
    vv <- v[idx]; cc <- cate_pct[idx]; ux <- sort(unique(vv))
    if (length(ux) <= 2 && all(ux %in% c(0,1))) {
      mean0 <- mean(cc[vv==0], na.rm=TRUE); mean1 <- mean(cc[vv==1], na.rm=TRUE)
      n0 <- sum(vv==0); n1 <- sum(vv==1)
      diag_cate_cov_overview[[paste(tname, cv, sep="|")]] <<- data.frame(
        treatment=tname, covariate=cv, type="binary", n=n0+n1,
        n0=n0, n1=n1, mean0=mean0, mean1=mean1, diff1minus0=mean1-mean0,
        stringsAsFactors=FALSE
      )
    } else {
      cor_val <- suppressWarnings(cor(vv, cc, use="pairwise.complete.obs"))
      slope <- tryCatch(unname(coef(lm(cc ~ vv))[2]), error=function(e) NA_real_)
      diag_cate_cov_overview[[paste(tname, cv, sep="|")]] <<- data.frame(
        treatment=tname, covariate=cv, type="numeric", n=length(vv),
        cor=cor_val, slope_catepct_per_unit=slope, stringsAsFactors=FALSE
      )
      q <- unique(quantile(vv, probs = seq(0,1,length.out=6), na.rm=TRUE))
      if (length(q) > 2) {
        b <- cut(vv, breaks=q, include.lowest=TRUE)
        tmp <- data.table(bin=b, cate_pct=cc)
        agg <- tmp[, .(n=.N, mean_cate_pct=mean(cate_pct), sd_cate_pct=sd(cate_pct)), by=bin]
        agg[, `:=`(treatment=tname, covariate=cv)]
        diag_cate_cov_bins[[paste(tname, cv, sep="|")]] <<- as.data.frame(agg)
      }
    }
  }
  
  # ---- BLP per covariate (kept) ----
  blp_rows <- vector("list", length = ncol(X_kept)); names(blp_rows) <- colnames(X_kept)
  for (j in seq_len(ncol(X_kept))) {
    cv <- colnames(X_kept)[j]
    xj <- as.numeric(X_kept[, j])
    ok <- is.finite(xj) & is.finite(cate_log)
    if (sum(ok) < 30 || sd(xj[ok]) < 1e-8) next
    blp_obj <- safe_blp_single(cf, matrix(xj[ok], ncol = 1, dimnames = list(NULL, cv)))
    parsed  <- tryCatch(parse_blp_coef_se(blp_obj), error = function(e) NULL)
    if (!is.null(parsed) && is.finite(parsed$beta) && is.finite(parsed$se) && parsed$se > 0) {
      beta <- parsed$beta; se <- parsed$se; method <- "grf_blp_single"
    } else {
      df_tmp <- data.frame(cate_log = cate_log[ok], xj = xj[ok])
      fit <- lm(cate_log ~ xj, data = df_tmp)
      V   <- tryCatch(sandwich::vcovHC(fit, type = "HC2"), error = function(e) vcov(fit))
      beta <- unname(coef(fit)["xj"])
      se   <- sqrt(pmax(0, diag(V)[names(coef(fit))=="xj"]))
      method <- "ols_HC2"
      if (!is.finite(beta) || !is.finite(se) || se <= 0) next
    }
    z   <- beta / se
    pvl <- 2*pnorm(-abs(z))
    lo  <- beta - 1.96*se
    hi  <- beta + 1.96*se
    beta_pct <- 100*(exp(beta)-1)
    lo_pct   <- 100*(exp(lo)-1)
    hi_pct   <- 100*(exp(hi)-1)
    blp_rows[[j]] <- data.frame(
      treatment   = tname,
      covariate   = cv,
      method      = method,
      n_obs       = sum(ok),
      coef_log    = beta,
      se_log      = se,
      z_or_t      = z,
      p_value     = pvl,
      ci_lo_log   = lo,
      ci_hi_log   = hi,
      coef_pct    = beta_pct,
      ci_lo_pct   = lo_pct,
      ci_hi_pct   = hi_pct,
      stringsAsFactors = FALSE
    )
  }
  blp_tbl <- dplyr::bind_rows(blp_rows[sapply(blp_rows, is.data.frame)])
  if (nrow(blp_tbl)) diag_blp_rows[[tname]] <<- blp_tbl
  
  cat(sprintf("  [GRF] %s ATE(log)=%+.4f (%+.4f,%+.4f), z=%.2f, p=%.3g  → %+.2f%%\n",
              toupper(target_primary), est, lo, hi, z, p, 100*(exp(est)-1)))
  if (is.finite(ate_overlap_log)) cat(sprintf("  ATE(overlap)=%+.4f\n", ate_overlap_log))
  if (is.finite(ate_all_log))     cat(sprintf("  ATE(all)=%+.4f\n", ate_all_log))
  if (is.finite(ate_treated_log)) cat(sprintf("  ATT=%+.4f\n", ate_treated_log))
  if (is.finite(aipw_est))        cat(sprintf("  [AIPW] ATE(log)=%+.4f (se=%.4f)\n", aipw_est, aipw_se))
  
  list(rows = dplyr::bind_rows(row_grf, row_pds), cf = cf, cates = cate_df)
}

# ----------------------- Run all treatments -----------------------
all_rows  <- list()
forests   <- list()
all_cates <- list()

for (nm in names(bin_treatments)) {
  W <- bin_treatments[[nm]]
  cat(sprintf("\n[Run] %s ...\n", nm))
  if (is.null(W) || length(unique(W[!is.na(W)])) < 2) {
    cat(sprintf("  → skipped %s (no variation).\n", nm)); next
  }
  ans <- run_models_for_treatment(W, nm, y)
  if (is.null(ans)) next
  all_rows[[nm]]  <- ans$rows
  forests[[nm]]   <- ans$cf
  all_cates[[nm]] <- ans$cates
}

if (length(all_rows)) {
  out <- dplyr::bind_rows(all_rows)
  data.table::fwrite(out, "cf_pds_M3_results_combined.csv")
  out_grf <- out %>% dplyr::filter(model == "GRF")
  out_pds <- out %>% dplyr::filter(model == "PDS")
  if (nrow(out_grf)) out_grf$q_value_bh <- p.adjust(out_grf$p_value, method = "BH")
  if (nrow(out_pds)) out_pds$q_value_bh <- p.adjust(out_pds$p_value, method = "BH")
  data.table::fwrite(out_grf, "cf_M3_results_GRF.csv")
  data.table::fwrite(out_pds, "cf_M3_results_PDS.csv")
  
  # CATE rows
  if (length(all_cates)) data.table::fwrite(dplyr::bind_rows(all_cates), "cf_M3_results_CATEs.csv")
  
  # Diagnostics tables
  if (length(diag_overlap_rows)) data.table::fwrite(dplyr::bind_rows(diag_overlap_rows), "cf_M3_overlap_diagnostics.csv")
  if (length(diag_balance_rows)) data.table::fwrite(dplyr::bind_rows(diag_balance_rows), "cf_M3_balance_diagnostics.csv")
  if (length(diag_robust_rows))  data.table::fwrite(dplyr::bind_rows(diag_robust_rows),  "cf_M3_robustness.csv")
  
  # ---- GATES from in-memory CATEs ----
  safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
  boot_mean_ci_pct <- function(v, nboot = 500, seed = 87) {
    v <- v[is.finite(v)]; set.seed(seed)
    if (length(v) < 2) return(list(mean = mean(v), lo = NA_real_, hi = NA_real_))
    bs <- replicate(nboot, mean(sample(v, length(v), replace = TRUE)))
    list(mean = mean(v), lo = as.numeric(quantile(bs, 0.025)), hi = as.numeric(quantile(bs, 0.975)))
  }
  safe_stub <- function(s) gsub("[^A-Za-z0-9]+", "_", s)
  
  gates_rows <- list()
  if (length(all_cates)) {
    dir.create("gates_by_treatment", showWarnings = FALSE)
    for (tnm in names(all_cates)) {
      cdf <- all_cates[[tnm]]
      if (!is.data.frame(cdf) || !"cate_pct" %in% names(cdf)) next
      v <- cdf$cate_pct
      brks <- unique(as.numeric(quantile(v, seq(0, 1, 0.1), na.rm = TRUE)))
      if (length(brks) < 3) next
      dec <- cut(v, breaks = brks, include.lowest = TRUE, labels = FALSE)
      dtab <- data.table(decile = dec, cate_pct = v)[, {
        ci <- boot_mean_ci_pct(cate_pct)
        .(n = .N, mean_pct = ci$mean, lo_pct = ci$lo, hi_pct = ci$hi)
      }, by = decile][order(decile)]
      fwrite(dtab, file.path("gates_by_treatment",
                             sprintf("gates_deciles_%s.csv", safe_stub(tnm))))
      d1  <- dtab[1]; d10 <- dtab[.N]
      gates_rows[[tnm]] <- data.table(
        treatment        = tnm,
        GATES_bot_mean   = d1$mean_pct,  GATES_bot_lo = d1$lo_pct,  GATES_bot_hi = d1$hi_pct,
        GATES_top_mean   = d10$mean_pct, GATES_top_lo = d10$lo_pct, GATES_top_hi = d10$hi_pct,
        GATES_lift       = d10$mean_pct - d1$mean_pct
      )
    }
  }
  gates_tbl <- if (length(gates_rows)) data.table::rbindlist(gates_rows, fill = TRUE) else data.table()
  
  # ---------- Forest-ready treatment summary (GRF only) ----------
  parse_treatment_label <- function(tnm) {
    # singles
    m <- regexpr("^(seller|buyer)_(N|N_1)_(total_value|gaspaid|gaslimit_cost)_high$", tnm, perl = TRUE)
    if (m[1] >= 0) {
      side <- if (grepl("^seller", tnm)) "Seller" else "Buyer"
      when <- if (grepl("_N_1_", tnm)) "N-1" else "N"
      metric_raw <- sub("^.*_(total_value|gaspaid|gaslimit_cost)_high$", "\\1", tnm)
      metric <- switch(metric_raw,
                       "total_value"   = "Total Value",
                       "gaspaid"       = "Gas Paid",
                       "gaslimit_cost" = "Gas Limit Cost")
      subgroup <- sprintf("%s %s High", side, when)
      subgroup_order <- switch(subgroup,
                               "Seller N High"   = 1L, "Seller N-1 High" = 2L,
                               "Buyer N High"    = 3L, "Buyer N-1 High"  = 4L, 99L)
      metric_order <- switch(metric, "Total Value"=1L, "Gas Paid"=2L, "Gas Limit Cost"=3L, 99L)
      return(list(subgroup=subgroup, metric=metric, pretty=sprintf("%s — %s", subgroup, metric),
                  subgroup_order=subgroup_order, metric_order=metric_order))
    }
    
    # 2-way 11vs00
    if (tnm == "seller_N_Tvalue_Gaslimit_11vs00")
      return(list(subgroup="Combination — Seller N High",   metric="Total Value + Gas Limit",
                  pretty="Combination — Seller N High — Total Value + Gas Limit", subgroup_order=5L, metric_order=1L))
    if (tnm == "seller_N_Tvalue_Gaspaid_11vs00")
      return(list(subgroup="Combination — Seller N High",   metric="Total Value + Gas Paid",
                  pretty="Combination — Seller N High — Total Value + Gas Paid", subgroup_order=5L, metric_order=2L))
    if (tnm == "seller_N_1_Tvalue_Gaslimit_11vs00")
      return(list(subgroup="Combination — Seller N-1 High", metric="Total Value + Gas Limit",
                  pretty="Combination — Seller N-1 High — Total Value + Gas Limit", subgroup_order=6L, metric_order=1L))
    if (tnm == "seller_N_1_Tvalue_Gaspaid_11vs00")
      return(list(subgroup="Combination — Seller N-1 High", metric="Total Value + Gas Paid",
                  pretty="Combination — Seller N-1 High — Total Value + Gas Paid", subgroup_order=6L, metric_order=2L))
    if (tnm == "buyer_N_Tvalue_Gaslimit_11vs00")
      return(list(subgroup="Combination — Buyer N High",    metric="Total Value + Gas Limit",
                  pretty="Combination — Buyer N High — Total Value + Gas Limit", subgroup_order=7L, metric_order=1L))
    if (tnm == "buyer_N_Tvalue_Gaspaid_11vs00")
      return(list(subgroup="Combination — Buyer N High",    metric="Total Value + Gas Paid",
                  pretty="Combination — Buyer N High — Total Value + Gas Paid", subgroup_order=7L, metric_order=2L))
    if (tnm == "buyer_N_1_Tvalue_Gaslimit_11vs00")
      return(list(subgroup="Combination — Buyer N-1 High",  metric="Total Value + Gas Limit",
                  pretty="Combination — Buyer N-1 High — Total Value + Gas Limit", subgroup_order=8L, metric_order=1L))
    if (tnm == "buyer_N_1_Tvalue_Gaspaid_11vs00")
      return(list(subgroup="Combination — Buyer N-1 High",  metric="Total Value + Gas Paid",
                  pretty="Combination — Buyer N-1 High — Total Value + Gas Paid", subgroup_order=8L, metric_order=2L))
    
    # conditionals
    if (tnm == "seller_N_Gaslimit_given_TvalueHigh")
      return(list(subgroup="Combination — Seller N High",   metric="Gas Limit | Total Value",
                  pretty="Combination — Seller N High — Gas Limit | Total Value", subgroup_order=5L, metric_order=3L))
    if (tnm == "seller_N_Gaspaid_given_TvalueHigh")
      return(list(subgroup="Combination — Seller N High",   metric="Gas Paid | Total Value",
                  pretty="Combination — Seller N High — Gas Paid | Total Value", subgroup_order=5L, metric_order=4L))
    if (tnm == "seller_N_1_Gaslimit_given_TvalueHigh")
      return(list(subgroup="Combination — Seller N-1 High", metric="Gas Limit | Total Value",
                  pretty="Combination — Seller N-1 High — Gas Limit | Total Value", subgroup_order=6L, metric_order=3L))
    if (tnm == "seller_N_1_Gaspaid_given_TvalueHigh")
      return(list(subgroup="Combination — Seller N-1 High", metric="Gas Paid | Total Value",
                  pretty="Combination — Seller N-1 High — Gas Paid | Total Value", subgroup_order=6L, metric_order=4L))
    if (tnm == "buyer_N_Gaslimit_given_TvalueHigh")
      return(list(subgroup="Combination — Buyer N High",    metric="Gas Limit | Total Value",
                  pretty="Combination — Buyer N High — Gas Limit | Total Value", subgroup_order=7L, metric_order=3L))
    if (tnm == "buyer_N_Gaspaid_given_TvalueHigh")
      return(list(subgroup="Combination — Buyer N High",    metric="Gas Paid | Total Value",
                  pretty="Combination — Buyer N High — Gas Paid | Total Value", subgroup_order=7L, metric_order=4L))
    if (tnm == "buyer_N_1_Gaslimit_given_TvalueHigh")
      return(list(subgroup="Combination — Buyer N-1 High",  metric="Gas Limit | Total Value",
                  pretty="Combination — Buyer N-1 High — Gas Limit | Total Value", subgroup_order=8L, metric_order=3L))
    if (tnm == "buyer_N_1_Gaspaid_given_TvalueHigh")
      return(list(subgroup="Combination — Buyer N-1 High",  metric="Gas Paid | Total Value",
                  pretty="Combination — Buyer N-1 High — Gas Paid | Total Value", subgroup_order=8L, metric_order=4L))
    # -------------- Add these cases inside parse_treatment_label(tnm) --------------
    if (tnm == "buyer_N_1_Tvalue_Gaspaid_Gaslimit_11vs00")
      return(list(subgroup="Combination — Buyer N-1 High", metric="Total Value + Gas Paid + Gas Limit",
                  pretty="Combination — Buyer N-1 High — TV + Gas Paid + Gas Limit",
                  subgroup_order=4L, metric_order=5L))
    
    if (tnm == "combo_all_TV_BN1_SN_SN1_11vs00")
      return(list(subgroup="Combination — All", metric="BN-1 TV + S:N TV + S:N-1 TV",
                  pretty="Combination — All — BN-1 TV + S:N TV + S:N-1 TV",
                  subgroup_order=8L, metric_order=1L))
    
    if (tnm == "comboB_lowTV_all_11vs00")
      return(list(subgroup="Strategy", metric="BN-1 TV ↓ + S:N TV ↓ + S:N-1 TV ↓",
                  pretty="Strategy — BN-1 TV ↓ + S:N TV ↓ + S:N-1 TV ↓",
                  subgroup_order=20L, metric_order=10L))
    
    if (tnm == "comboB_BN1_lowTV_highFees_11vs00")
      return(list(subgroup="Strategy", metric="BN-1 TV ↓ + BN-1 GL ↑ + BN-1 GP ↑",
                  pretty="Strategy — BN-1 TV ↓ + BN-1 GL ↑ + BN-1 GP ↑",
                  subgroup_order=20L, metric_order=11L))
    
    if (tnm == "comboB_SN_lowTV_highFees_11vs00")
      return(list(subgroup="Strategy", metric="S:N TV ↓ + S:N GL ↑ + S:N GP ↑",
                  pretty="Strategy — S:N TV ↓ + S:N GL ↑ + S:N GP ↑",
                  subgroup_order=20L, metric_order=12L))
    
    if (tnm == "comboB_mix_lowTV_highFees_11vs00")
      return(list(subgroup="Strategy", metric="BN-1 (TV ↓, GL ↑, GP ↑) + S:N (TV ↓, GL ↑, GP ↑)",
                  pretty="Strategy — BN-1 (TV ↓, GL ↑, GP ↑) + S:N (TV ↓, GL ↑, GP ↑)",
                  subgroup_order=20L, metric_order=13L))
    
    # multi-party TVs
    if (tnm == "all_total_value_11vs00")
      return(list(subgroup="Combination", metric="All Sides: Total Value",
                  pretty="Combination — All Sides: Total Value", subgroup_order=9L, metric_order=1L))
    if (tnm == "buyer_total_value_11vs00")
      return(list(subgroup="Combination — Buyer", metric="B:N + B:N-1: Total Value",
                  pretty="Combination — Buyer — B:N + B:N-1: Total Value", subgroup_order=7L, metric_order=1L))
    if (tnm == "seller_total_value_11vs00")
      return(list(subgroup="Combination — Seller", metric="S:N + S:N-1: Total Value",
                  pretty="Combination — Seller — S:N + S:N-1: Total Value", subgroup_order=5L, metric_order=1L))
    
    # strategies
    if (tnm == "comboB1_strategy_buyer_push_11vs00")
      return(list(subgroup="Strategy", metric="B:N TV↑ + B:N-1 TV↑ + B:N GL↑ + S:N TV↓ + S:N-1 TV↓ + S:N GP↑",
                  pretty="Strategy — Buyer Push (BN/BN-1 TV↑, GL↑; Sellers TV↓, GP↑)", subgroup_order=20L, metric_order=99L))
    if (tnm == "comboB2_strategy_cross_fees_11vs00")
      return(list(subgroup="Strategy", metric="B:N GP↓ + B:N-1 GP↓ + S:N GP↑ + S:N-1 GP↑ + B:N TV↑ + S:N TV↓",
                  pretty="Strategy — Cross Fees (Buyers GP↓; Sellers GP↑; BN TV↑; SN TV↓)", subgroup_order=20L, metric_order=99L))
    if (tnm == "comboB3_strategy_liquidity_cycle_11vs00")
      return(list(subgroup="Strategy", metric="B:N GL↑ + B:N-1 GL↑ + S:N GL↑ + B:N TV↑ + S:N TV↓",
                  pretty="Strategy — Liquidity Cycle (Buyers GL↑; Sellers GL↑; TV diverge)", subgroup_order=20L, metric_order=99L))
    if (tnm == "comboB4_strategy_sellers_strong_11vs00")
      return(list(subgroup="Strategy", metric="S:N TV↑ + S:N-1 TV↑ + B:N TV↓ + B:N-1 TV↓",
                  pretty="Strategy — Sellers Strong (S TVs↑; B TVs↓)", subgroup_order=20L, metric_order=99L))
    
    list(subgroup = tnm, metric = tnm, pretty = tnm, subgroup_order = 99L, metric_order = 99L)
  }
  
  ov_keep <- if (length(diag_overlap_rows)) {
    as.data.table(dplyr::bind_rows(diag_overlap_rows))[, .(treatment, n_kept, kept_pct)]
  } else data.table(treatment = character(), n_kept = integer(), kept_pct = numeric())
  
  # Join pieces & write forest-ready CSV
  if (nrow(out_grf <- out_grf) == 0) out_grf <- out %>% dplyr::filter(model == "GRF")
  lab <- lapply(out_grf$treatment, parse_treatment_label)
  lab_df <- dplyr::bind_rows(lapply(lab, as.data.frame))
  out_grf2 <- cbind(out_grf, lab_df)
  
  # attach GATES if available
  if (nrow(gates_tbl)) {
    out_grf2 <- out_grf2 %>%
      dplyr::left_join(gates_tbl, by = c("treatment" = "treatment"))
  }
  
  treat_summary <- out_grf2 %>%
    dplyr::left_join(ov_keep, by = "treatment") %>%
    dplyr::mutate(
      ATE      = ATE_log,    CI_lower = CI95_lo_log, CI_upper = CI95_hi_log, StdErr = SE_log,
      ATE_pct      = 100 * (exp(ATE_log)     - 1),
      CI_lower_pct = 100 * (exp(CI95_lo_log) - 1),
      CI_upper_pct = 100 * (exp(CI95_hi_log) - 1),
      treated_prop = treated_share,
      treated_kept = treated_share * n_kept,
      control_kept = (1 - treated_share) * n_kept
    ) %>%
    dplyr::arrange(subgroup_order, metric_order, dplyr::desc(abs(ATE_pct))) %>%
    as.data.table()
  
  if (!"q_value_bh" %in% names(treat_summary)) {
    treat_summary[, q_value_bh := p.adjust(p_value, method = "BH")]
  }
  
  fwrite(treat_summary[, .(
    treatment,
    display_label = pretty, subgroup, metric,
    n_kept, kept_pct, treated_prop, treated_kept, control_kept,
    ATE, CI_lower, CI_upper, StdErr,
    ATE_pct, CI_lower_pct, CI_upper_pct,
    p_value, q_value_bh,
    GATES_lift, GATES_bot_mean, GATES_top_mean
  )], "cf_M3_treatment_results_with_q.csv")
  
  # ---------------- Dashboards ----------------
  thresholds_for_join <- treat_threshold_df %>%
    dplyr::select(treatment, split_rule, cutoff_raw)
  
  make_dash <- function(res_tbl, label) {
    if (!nrow(res_tbl)) return(NULL)
    res_tbl$q_value_bh <- p.adjust(res_tbl$p_value, method = "BH")
    
    bal_join <- if (length(diag_balance_rows)) {
      bal_all <- dplyr::bind_rows(diag_balance_rows)
      bal_all$abs_smd_kept <- abs(bal_all$SMD_kept)
      bal_all %>% dplyr::group_by(treatment) %>% dplyr::summarise(
        total_covariates = sum(!is.na(abs_smd_kept)),
        count_gt_0_10    = sum(abs_smd_kept > 0.10, na.rm = TRUE),
        count_gt_0_25    = sum(abs_smd_kept > 0.25, na.rm = TRUE),
        median_abs_smd   = median(abs_smd_kept, na.rm = TRUE),
        p90_abs_smd      = quantile(abs_smd_kept, 0.90, na.rm = TRUE),
        max_abs_smd      = max(abs_smd_kept, na.rm = TRUE),
        .groups = "drop"
      )
    } else {
      # ← ensure columns exist so mutate() doesn't error
      tibble::tibble(
        treatment        = unique(res_tbl$treatment),
        total_covariates = NA_integer_,
        count_gt_0_10    = NA_integer_,
        count_gt_0_25    = NA_integer_,
        median_abs_smd   = NA_real_,
        p90_abs_smd      = NA_real_,
        max_abs_smd      = NA_real_
      )
    }
    
    gates_min <- if (exists("gates_tbl") && nrow(gates_tbl)) {
      gates_tbl[, .(treatment, GATES_lift, GATES_bot_mean, GATES_top_mean)]
    } else data.table(treatment = character())
    
    dash <- res_tbl %>%
      dplyr::left_join(dplyr::select(treat_summary, treatment, subgroup, metric), by="treatment") %>%
      dplyr::left_join(dplyr::select(dplyr::bind_rows(diag_overlap_rows), treatment, kept_pct), by="treatment") %>%
      dplyr::left_join(bal_join, by="treatment") %>%
      dplyr::left_join(thresholds_for_join, by="treatment") %>%
      dplyr::left_join(gates_min, by="treatment") %>%
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
        flag_tiny = ifelse(!is.na(ATE_pct) & abs(ATE_pct) < 1, "tiny", ""),
        flags = paste0(flag_signif,
                       ifelse(flag_overlap != "", paste0(",", flag_overlap), ""),
                       ifelse(flag_balance != "", paste0(",", flag_balance), ""),
                       ifelse(flag_tiny    != "", paste0(",", flag_tiny), "")) %>% sub("^,", "", .)
      ) %>%
      dplyr::select(
        model, treatment, subgroup, metric,
        ATE_pct, CI95_lo_pct, CI95_hi_pct,
        GATES_lift, GATES_bot_mean, GATES_top_mean,
        cutoff_raw, split_rule, p_value, q_value_bh, kept_pct,
        total_covariates, count_gt_0_10, count_gt_0_25, median_abs_smd, p90_abs_smd, max_abs_smd, flags
      )
    
    out_file <- paste0("cf_M3_dashboard_", label, ".csv")
    data.table::fwrite(dash, out_file)
    out_file
  }
  
  file_dash_grf <- make_dash(out_grf, "GRF")
  file_dash_pds <- make_dash(out_pds, "PDS")
  
  if (!is.null(file_dash_grf) && !is.null(file_dash_pds)) {
    comp <- dplyr::full_join(
      out_grf %>% dplyr::select(treatment, ATE_pct_grf = ATE_pct, p_grf = p_value),
      out_pds %>% dplyr::select(treatment, ATE_pct_pds = ATE_pct, p_pds = p_value),
      by = "treatment"
    )
    data.table::fwrite(comp, "cf_M3_dashboard_compare_GRF_vs_PDS.csv")
  }
  
  cat("\n[Saved]\n",
      "  - cf_pds_M3_results_combined.csv\n",
      "  - cf_M3_results_GRF.csv\n",
      "  - cf_M3_results_PDS.csv\n",
      "  - cf_M3_treatment_thresholds.csv\n",
      "  - cf_M3_results_CATEs.csv\n",
      "  - cf_M3_overlap_diagnostics.csv\n",
      "  - cf_M3_balance_diagnostics.csv\n",
      "  - cf_M3_treatment_results_with_q.csv\n",
      "  - cf_M3_dashboard_GRF.csv\n",
      "  - cf_M3_dashboard_PDS.csv\n",
      "  - cf_M3_dashboard_compare_GRF_vs_PDS.csv\n",
      sep = "")
} else {
  cat("\n[Note] No results produced; check overlap trimming & treatment variation.\n")
}

# ----------------------- Save R2_tau (per treatment) -----------------------
if (length(diag_r2tau_rows)) {
  r2tau_tab <- dplyr::bind_rows(diag_r2tau_rows)
  data.table::fwrite(r2tau_tab, "cf_M3_R2tau_by_treatment.csv")
  cat("  - cf_M3_R2tau_by_treatment.csv\n")
}

# ----------------------- Save BLP outputs (+ join R2_tau) -----------------------
if (length(diag_blp_rows)) {
  blp_out <- dplyr::bind_rows(diag_blp_rows)
  if (exists("r2tau_tab")) {
    blp_out <- dplyr::left_join(
      blp_out,
      r2tau_tab %>% dplyr::select(treatment, r2_tau, r2_tau_oos),
      by = "treatment"
    )
  }
  data.table::fwrite(blp_out, "cf_M3_BLP_per_covariate.csv")
  cat("  - cf_M3_BLP_per_covariate.csv\n")
}
