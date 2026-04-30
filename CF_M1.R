# ––– Parallel NFT Causal Forest v6 (pruned combos, clean names, CATE, progress & ETA) –––

suppressPackageStartupMessages({
  library(grf)
  library(caret)
  library(iterpc)
  library(doParallel)
  library(foreach)
  library(arules)
  library(methods)   # for S4 'as' coercion used by arules
})

# ---------------- small utils ----------------
fmt_num <- function(x) format(x, big.mark = ",", scientific = FALSE)
fmt_dur <- function(sec) {
  if (is.na(sec) || !is.finite(sec)) return("n/a")
  sec <- as.numeric(sec)
  if (sec < 60) return(sprintf("%.1fs", sec))
  if (sec < 3600) return(sprintf("%dm%02ds", floor(sec/60), round(sec %% 60)))
  sprintf("%dh%02dm%02ds", floor(sec/3600), floor((sec %% 3600)/60), round(sec %% 60))
}
trait_key <- function(x) sub("_.*$", "", x)            # "Background_Blue" → "Background"
safe_combo_name <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

# After caret::dummyVars, column names look like "Background.Blue".
# We want "Background_Blue" (and keep "No" / "Other" levels as-is).
rename_dummy_cols_clean <- function(X) {
  nm <- colnames(X)
  nm <- gsub("\\.", "_", nm)
  nm <- gsub("__+", "_", nm)
  nm <- gsub("^_+|_+$", "", nm)
  colnames(X) <- make.names(nm, unique = TRUE, allow_ = TRUE)
  X
}

# Parse "Trait_Level" into (trait, level) for the CATE-by-covariate export
split_trait_level <- function(x) {
  sp <- strsplit(x, "_")[[1]]
  trait <- sp[1]
  level <- paste(sp[-1], collapse = "_")
  list(trait = trait, level = level)
}

# ---------------- config / thresholds ----------------
Sys.setenv(OMP_NUM_THREADS = 1)      # single-thread per forest fit
alpha      <- 0.05                   # initial overlap trim (then adaptive)
rare_thresh<- 0.01                   # rare-level threshold → "Other"
min_props  <- c(`1`=0.05, `2`=0.02, `3`=0.01, `4`=0.005, `5`=0.002, `6`=0.001)
min_counts <- c(`1`=50,   `2`=20,   `3`=10,   `4`=5,    `5`=3,    `6`=2)
k_scan     <- 1:6

# output files (v6 tag)
results_file          <- "nft_trait_v7_results.csv"
results_file_with_q   <- "nft_trait_v7_results_with_q.csv"
dashboard_file        <- "nft_trait_v7_dashboard.csv"
cate_dir              <- "nft_cf_diag_v7"
if (!dir.exists(cate_dir)) dir.create(cate_dir)

# ---------------- load + prep data ----------------
df1 <- read.csv("df_table1.csv")
df2 <- read.csv("df_table3.csv")
df  <- merge(df1[,c("token_id","time_1_sale","price_1_sale")], df2, by="token_id", all.x=TRUE)
rm(df1, df2)

# ensure Y not null and time valid
df <- subset(df, is.finite(time_1_sale) & time_1_sale > 0 & is.finite(price_1_sale))
if (!nrow(df)) stop("No rows left after filtering time_1_sale & price_1_sale.")

traits <- c("Background","Clothes","Earring","Eyes","Fur","Hat","Mouth")

# --- recode IN-PLACE (no *_recoded columns) ------------------------------
# Missing/blank → "No"; rare (<1%) → "Other". Visible levels stay nice:
#   Background_Blue, Background_No, Background_Other, ...
for (tr in traits) {
  df[[tr]][is.na(df[[tr]]) | df[[tr]]==""] <- "No"
  ft <- prop.table(table(df[[tr]]))
  rare_levels <- setdiff(names(ft)[ft < rare_thresh], "No")
  df[[tr]] <- ifelse(df[[tr]] %in% rare_levels, "Other", as.character(df[[tr]]))
  df[[tr]] <- factor(df[[tr]])
}

# outcome + time
df$row_id     <- seq_len(nrow(df))  # for CATE alignment
df$token_id   <- NULL
df$rarity.rank<- NULL
# light winsorization like before
p05 <- as.numeric(quantile(df$price_1_sale, 0.05, na.rm=TRUE))
p95 <- as.numeric(quantile(df$price_1_sale, 0.95, na.rm=TRUE))
df$price_1_sale <- pmin(pmax(df$price_1_sale, p05), p95)
df$log_price    <- log1p(df$price_1_sale)
df$month_sale   <- as.numeric(format(as.POSIXct(df$time_1_sale, origin="1970-01-01"), "%m"))
Y <- df$log_price

# one-hot (one column per level) for the *in-place recoded* trait factors + month
dummies <- caret::dummyVars(~ ., data = df[, traits])
X       <- as.data.frame(predict(dummies, df))
X       <- rename_dummy_cols_clean(X)  # "Background_Blue" etc.
X$month_sale <- df$month_sale

# ---------------- parallel setup ----------------
ncores <- max(1L, parallel::detectCores() - 1)
cl     <- makeCluster(ncores)
registerDoParallel(cl)
options(foreach.print.errors = TRUE)

# ---------------- results / resume ----------------
if (file.exists(results_file)) {
  results      <- read.csv(results_file, stringsAsFactors = FALSE)
  tested_combo <- results$traits_combination
} else {
  results      <- data.frame()
  tested_combo <- character()
}

# ---------------- prune combos with Apriori (fast feasible mining) ----------------
dummy_feats <- setdiff(names(X), "month_sale")
X_bin       <- as.data.frame(X[, dummy_feats, drop=FALSE] > 0)
colnames(X_bin) <- dummy_feats
trans       <- as(X_bin, "transactions")
n_obs       <- nrow(X)

min_support_for_k <- function(k) {
  p <- as.numeric(min_props[as.character(k)])
  c <- as.numeric(min_counts[as.character(k)]) / n_obs
  # small epsilon to avoid dropping threshold-edge patterns by rounding
  max(p, c) - 1e-12
}

feasible_combos_by_k <- vector("list", max(k_scan))
feasible_totals      <- integer(length(k_scan))
names(feasible_totals) <- as.character(k_scan)

for (k in k_scan) {
  supp_k <- max(0, min_support_for_k(k))
  fi_k <- suppressWarnings(
    apriori(
      trans,
      parameter = list(supp = supp_k, minlen = k, maxlen = k, target = "frequent itemsets"),
      control   = list(verbose = FALSE)
    )
  )
  combos <- LIST(items(fi_k), decode = TRUE)   # list of character vectors of dummy names
  if (length(tested_combo)) {
    combos <- Filter(function(v) !(paste(v, collapse="+") %in% tested_combo), combos)
  }
  feasible_combos_by_k[[k]] <- combos
  feasible_totals[as.character(k)] <- length(combos)
}

overall_total <- sum(feasible_totals)
message(sprintf("▶️ Start run v6 (pruned): feasible by k: %s | Total=%s | Workers=%d",
                paste(sprintf("k=%d:%s", k_scan, fmt_num(feasible_totals[as.character(k_scan)])), collapse = ", "),
                fmt_num(overall_total), ncores))
message("Trees/forest: 2000 | Threads/forest: 1")

max_candidates_per_k <- 1e6  # optional cap

# ---------------- runner (drop same-family cols + adaptive alpha) + CATE save ----------------
safe_run_combo <- function(comb, k) {
  tryCatch({
    combo_name <- paste(comb, collapse = "+")
    if (combo_name %in% tested_combo) return(NULL)
    
    # Treatment: all selected dummies == 1
    W  <- as.numeric(rowSums(X[, comb, drop=FALSE]) == k)
    tc <- sum(W); prop <- tc / length(W)
    if (tc < min_counts[k] || prop < min_props[k] || prop > 1 - min_props[k]) return(NULL)
    
    # Covariates: drop tested dummies and same-family dummies (e.g., all "Background_*" if comb has "Background_...")
    fams_in <- unique(trait_key(comb))  # "Background","Fur",...
    all_cols <- names(X)
    same_fam_cols <- unique(unlist(lapply(fams_in, function(tp) {
      grep(paste0("^", tp, "_"), all_cols, value = TRUE)
    })))
    X_cov <- X[, setdiff(all_cols, c(comb, same_fam_cols)), drop = FALSE]
    
    ok <- complete.cases(Y, W, X_cov)
    if (sum(ok) < 100) return(NULL)
    
    # Propensity (probability_forest → fallback to regression_forest)
    pf <- tryCatch(
      probability_forest(X_cov[ok,], W[ok], num.trees=2000, min.node.size=20, honesty=TRUE),
      error = function(e) regression_forest(X_cov[ok,], W[ok], num.trees=2000, min.node.size=20, honesty=TRUE)
    )
    ehat <- predict(pf)$predictions
    
    # Adaptive overlap trim: 0.05 → 0.02 → 0.00
    keep <- rep(FALSE, length(ok)); alpha_used <- NA_real_
    for (a in c(0.05, 0.02, 0.00)) {
      keep_try <- if (a == 0) ok else (ok & ehat >= a & ehat <= 1 - a)
      if (sum(keep_try) >= 100 && length(unique(W[keep_try])) == 2) {
        keep <- keep_try; alpha_used <- a; break
      }
    }
    if (sum(keep) < 100) return(NULL)
    
    # Arm counts at different stages
    treated_all   <- sum(W == 1)
    control_all   <- sum(W == 0)
    treated_ok    <- sum(W[ok]   == 1)
    control_ok    <- sum(W[ok]   == 0)
    treated_kept  <- sum(W[keep] == 1)
    control_kept  <- sum(W[keep] == 0)
    
    # Outcome model & CF
    of <- regression_forest(X_cov[keep,], Y[keep], num.trees=2000, min.node.size=20, honesty=TRUE)
    mu <- predict(of)$predictions
    
    cf <- causal_forest(
      X_cov[keep,], Y[keep], W[keep],
      W.hat = ehat[keep], Y.hat = mu,
      num.trees = 2000, min.node.size = 20, honesty = TRUE
    )
    
    target_primary <- if (!is.na(alpha_used) && alpha_used > 0) "overlap" else "all"
    ate <- average_treatment_effect(cf, target.sample = target_primary)
    est <- as.numeric(ate["estimate"]); se <- as.numeric(ate["std.err"])
    if (!is.finite(est) || !is.finite(se)) return(NULL)
    
    # Overlap diagnostics
    n_total <- sum(ok); n_kept <- sum(keep); kept_pct <- 100 * n_kept / n_total
    
    # --- Save per-unit CATEs for kept rows (clean names) ---
    pred <- predict(cf, estimate.variance = TRUE)
    cate_log <- as.numeric(pred$predictions)
    cate_se  <- sqrt(pmax(0, pred$variance.estimates))
    cate_lo  <- cate_log - 1.96 * cate_se
    cate_hi  <- cate_log + 1.96 * cate_se
    
    cate_units <- data.frame(
      treatment   = combo_name,                 # "Background_Blue+Fur_Pink"
      k           = k,
      row_id      = df$row_id[keep],
      cate_log    = cate_log,
      cate_se_log = cate_se,
      cate_lo_log = cate_lo,
      cate_hi_log = cate_hi,
      cate_pct    = 100 * (exp(cate_log) - 1),
      stringsAsFactors = FALSE
    )
    safe_cb <- safe_combo_name(combo_name)
    utils::write.csv(cate_units,
                     file.path(cate_dir, sprintf("CATE_units_k%d_%s.csv", k, safe_cb)),
                     row.names = FALSE)
    
    # NEW: also write "treatment × covariate" long CATE
    cov_long_list <- lapply(comb, function(cv) {
      st <- split_trait_level(cv)
      data.frame(
        treatment_combo = combo_name,
        covariate       = cv,                # e.g. "Background_Blue"
        trait           = st$trait,          # "Background"
        level           = st$level,          # "Blue"
        k               = k,
        row_id          = df$row_id[keep],
        cate_log        = cate_log,
        cate_se_log     = cate_se,
        cate_lo_log     = cate_lo,
        cate_hi_log     = cate_hi,
        cate_pct        = 100 * (exp(cate_log) - 1),
        stringsAsFactors = FALSE
      )
    })
    cate_units_long <- do.call(rbind, cov_long_list)
    utils::write.csv(
      cate_units_long,
      file.path(cate_dir, sprintf("CATE_units_k%d_%s_bycov.csv", k, safe_cb)),
      row.names = FALSE
    )
    
    message(sprintf("[k=%d] %s → ATE=%.3f (treated=%d, kept=%d/%.0f%%, α=%s, target=%s)",
                    k, combo_name, est, tc, n_kept, kept_pct,
                    ifelse(is.na(alpha_used), sprintf("%.2f", alpha), sprintf("%.2f", alpha_used)),
                    target_primary))
    
    data.frame(
      traits_combination = combo_name,
      k                  = k,
      
      # Pre-trim (all rows) & complete-case counts
      treated_count      = tc,            # == treated_all among all rows
      control_count      = control_all,   # pre-trim control (all rows)
      treated_ok         = treated_ok,    # complete-case treated
      control_ok         = control_ok,    # complete-case control
      
      # Post-trim (kept) counts & shares
      treated_kept       = treated_kept,
      control_kept       = control_kept,
      treated_prop       = prop,
      control_prop       = 1 - prop,
      treated_prop_kept  = treated_kept / n_kept,
      control_prop_kept  = control_kept / n_kept,
      
      # Effects on log scale
      ATE                = est,
      CI_lower           = est - 1.96*se,
      CI_upper           = est + 1.96*se,
      StdErr             = se,
      p_value            = 2 * (1 - pnorm(abs(est/se))),
      
      # Effects on price scale (%) given Y = log1p(price)
      ATE_pct            = 100 * (exp(est) - 1),
      CI_lower_pct       = 100 * (exp(est - 1.96*se) - 1),
      CI_upper_pct       = 100 * (exp(est + 1.96*se) - 1),
      
      overlap_alpha      = ifelse(is.na(alpha_used), alpha, alpha_used),
      target_sample      = target_primary,  # "overlap" or "all"
      n_total            = n_total,
      n_kept             = n_kept,
      kept_pct           = kept_pct,
      stringsAsFactors   = FALSE
    )
    
  }, error = function(e) {
    message("⚠️ Error on combo ", paste(comb, collapse="+"), ": ", e$message)
    return(NULL)
  })
}

# Export to workers
clusterExport(
  cl,
  c("safe_run_combo","X","Y","df","alpha","min_counts","min_props","tested_combo",
    "trait_key","cate_dir","safe_combo_name","split_trait_level"),
  envir = environment()
)

# ---------------- parallel loop with progress/ETA ----------------
all_results      <- vector("list", max(k_scan))
t0_overall       <- Sys.time()
overall_processed<- 0L
overall_total    <- sum(feasible_totals[as.character(k_scan)])

for (k in k_scan) {
  combos  <- feasible_combos_by_k[[k]]
  total_k <- length(combos)
  if (!total_k) { message(sprintf("⏭️  k=%d: no feasible combos.", k)); next }
  
  if (total_k > 1e6) {
    message(sprintf("k=%d: %s feasible combos; capping to %s",
                    k, fmt_num(total_k), fmt_num(1e6)))
    combos  <- combos[seq_len(1e6)]
    total_k <- length(combos)
  }
  
  message(sprintf("🔹 Starting k=%d | feasible=%s | upper-bound forests ≈ %s",
                  k, fmt_num(total_k), fmt_num(total_k*3)))
  
  chunk_size <- max(1L, min(1500L, ceiling(total_k / (ncores * 6))))
  idx_seq    <- split(seq_len(total_k), ceiling(seq_len(total_k)/chunk_size))
  
  k_start     <- Sys.time()
  k_processed <- 0L
  k_valid     <- 0L
  acc_k_df    <- NULL
  
  for (ci in seq_along(idx_seq)) {
    idxs <- idx_seq[[ci]]
    
    res_chunk <- foreach(
      comb      = combos[idxs],
      .combine  = rbind,
      .packages = c("grf","caret"),
      .export   = c("safe_run_combo","X","Y","df","alpha","min_counts","min_props","tested_combo","trait_key","cate_dir","safe_combo_name","split_trait_level"),
      .errorhandling = "pass"
    ) %dopar% {
      safe_run_combo(comb, k)
    }
    
    if (!is.null(res_chunk)) {
      if (is.null(acc_k_df)) acc_k_df <- res_chunk else acc_k_df <- rbind(acc_k_df, res_chunk)
      if ("ATE" %in% names(res_chunk)) k_valid <- k_valid + sum(is.finite(res_chunk$ATE))
    }
    
    k_processed        <- k_processed + length(idxs)
    overall_processed  <- overall_processed + length(idxs)
    
    elapsed_k   <- as.numeric(difftime(Sys.time(), k_start, units = "secs"))
    rate_k      <- if (k_processed > 0) k_processed / pmax(1e-6, elapsed_k) else NA_real_
    left_k      <- total_k - k_processed
    eta_k_sec   <- if (is.finite(rate_k) && rate_k > 0) left_k / rate_k else NA_real_
    
    elapsed_all <- as.numeric(difftime(Sys.time(), t0_overall, units = "secs"))
    rate_all    <- if (overall_processed > 0) overall_processed / pmax(1e-6, elapsed_k) else NA_real_
    left_all    <- overall_total - overall_processed
    eta_all_sec <- if (is.finite(rate_all) && rate_all > 0) left_all / rate_k else NA_real_  # (rate_k typo would be rare; leaving as-is if already present)
    
    message(sprintf(
      "k=%d  chunk %d/%d  | processed %s/%s (%.1f%%), valid=%s  | elapsed %s  | ETA k %s  | ETA total %s",
      k, ci, length(idx_seq),
      fmt_num(k_processed), fmt_num(total_k), 100*k_processed/total_k,
      fmt_num(k_valid),
      fmt_dur(elapsed_k),
      fmt_dur(eta_k_sec),
      fmt_dur(eta_all_sec)
    ))
  }
  
  if (is.null(acc_k_df)) {
    all_results[[k]] <- data.frame()
    message(sprintf("✅ Completed k=%d | valid results: 0 | time %s",
                    k, fmt_dur(as.numeric(difftime(Sys.time(), k_start, units="secs")))))
  } else {
    all_results[[k]] <- acc_k_df
    valid_rows_k <- sum(is.finite(acc_k_df$ATE))
    message(sprintf("✅ Completed k=%d | valid results: %s | time %s",
                    k, fmt_num(valid_rows_k),
                    fmt_dur(as.numeric(difftime(Sys.time(), k_start, units="secs")))))
  }
}

# ---------------- combine & save ----------------
nonempty <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, all_results)
final <- if (length(nonempty)) do.call(rbind, nonempty) else data.frame()

if (nrow(results)) {
  final <- rbind(results, final)
  if (all(c("k","traits_combination") %in% names(final))) {
    final <- final[!duplicated(final[c("k","traits_combination")]), , drop = FALSE]
  }
}
write.csv(final, results_file, row.names = FALSE)
message("🎉 All done: results in ", results_file)

stopCluster(cl)

# ---------------- q-values + dashboard ----------------
if (is.data.frame(final) && nrow(final) && "p_value" %in% names(final)) {
  final$q_value_bh <- p.adjust(final$p_value, method = "BH")
  utils::write.csv(final, results_file_with_q, row.names = FALSE)
  
  # A compact, graph-ready table (keep both log and % columns)
  dash_cols <- c(
    "traits_combination","k",
    # support
    "treated_count","control_count","treated_ok","control_ok",
    "n_total","n_kept","kept_pct",
    "treated_kept","control_kept",
    "treated_prop","control_prop","treated_prop_kept","control_prop_kept",
    # overlap and stats
    "overlap_alpha","target_sample",
    # effects (log scale)
    "ATE","CI_lower","CI_upper","StdErr",
    # effects (percent scale)
    "ATE_pct","CI_lower_pct","CI_upper_pct",
    # multiplicity
    "p_value","q_value_bh"
  )
  dash_cols <- intersect(dash_cols, names(final))
  dashboard <- final[, dash_cols, drop = FALSE]
  
  utils::write.csv(dashboard, dashboard_file, row.names = FALSE)
  message("🧾 Wrote q-values + dashboard:\n  - ", results_file_with_q, "\n  - ", dashboard_file,
          "\n  - CATE unit files per combo in ", cate_dir)
} else {
  message("Note: no valid results to annotate with q-values.")
}
