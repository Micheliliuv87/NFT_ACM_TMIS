# ---- NFT Forest Plot (with GATES lift, robust-normalized box sizes) ----
# Reads: nft_trait_v7_results_with_q.csv  +  paper_TableC_CATE_main.csv
# Output: nft_forest_simple_K1to4_top5_V3.pdf + .png
# Requires: forestploter, grid, scales (and ragg for PNG)

suppressPackageStartupMessages({
  library(forestploter)
  library(grid)
  library(scales)
})

# ---------------- Config ----------------
infile       <- "nft_trait_v7_results_with_q.csv"
cate_main_in <- file.path("nft_cf_cate_reports_v7", "paper_assets", "paper_TableC_CATE_main.csv")
outfile_pdf  <- "nft_forest_simple_K1to4_top5_V3.pdf"
outfile_png  <- "nft_forest_simple_K1to4_top5_V3.png"

ks_to_plot  <- 1:4
top_n       <- 5
ci_spaces   <- 40         # width of CI lane (spaces)
q_max       <- 0.10       # <-- NEW: BH q-value threshold
wrap_width  <- 70         # wrap width for combo labels
use_shorten <- TRUE       # abbreviate long trait family prefixes
base_font   <- 12         # forest_theme(base_size)
keep_min    <- 50         # require kept_pct >= 60

# Box-size normalization band (robust)
size_floor <- 0.5        # min box size
size_cap   <- 0.9        # max box size
winsor_lo  <- 0.15        # lower quantile for winsorization
winsor_hi  <- 0.85        # upper quantile for winsorization

# ---------------- Helpers ----------------
fmt_int  <- function(x) format(round(as.numeric(x)), big.mark = ",", trim = TRUE)
fmt_pct1 <- function(x) sprintf("%.1f", as.numeric(x))
fmt_q3   <- function(x) sprintf("%.3f", as.numeric(x))
to_pct   <- function(z) 100 * (exp(as.numeric(z)) - 1)

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

wrap_combo_vec <- function(x, width = 30) {
  vapply(x, function(s) {
    if (is.na(s) || s == "") return("")
    s <- gsub("\\+", " + ", s)
    s <- gsub("_",  "_ ",  s)
    s <- gsub("\\.", "·",  s)
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1))
}

safe_combo_name <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

warn_once <- local({
  warned <- FALSE
  function(msg) { if (!warned) { message(msg); warned <<- TRUE } }
})

# ---------------- Load main results ----------------
df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

need <- c("k","traits_combination",
          "treated_count","control_count","treated_prop",
          "ATE","CI_lower","CI_upper","StdErr","q_value_bh",
          "n_kept","kept_pct")
miss <- setdiff(need, names(df))
if (length(miss)) stop("Missing columns in input: ", paste(miss, collapse = ", "))

# Optional post-trim counts
has_treated_kept <- "treated_kept" %in% names(df)
has_control_kept <- "control_kept" %in% names(df)

# Use % effect columns from file; fallback if absent/bad
have_pct_cols <- all(c("ATE_pct","CI_lower_pct","CI_upper_pct") %in% names(df))
if (have_pct_cols) {
  df$ATE_pct      <- as.numeric(df$ATE_pct)
  df$CI_lower_pct <- as.numeric(df$CI_lower_pct)
  df$CI_upper_pct <- as.numeric(df$CI_upper_pct)
  if (any(!is.finite(df$ATE_pct)) || any(!is.finite(df$CI_lower_pct)) || any(!is.finite(df$CI_upper_pct))) {
    df$ATE_pct      <- to_pct(df$ATE)
    df$CI_lower_pct <- to_pct(df$CI_lower)
    df$CI_upper_pct <- to_pct(df$CI_upper)
  }
} else {
  df$ATE_pct      <- to_pct(df$ATE)
  df$CI_lower_pct <- to_pct(df$CI_lower)
  df$CI_upper_pct <- to_pct(df$CI_upper)
}

# Ensure numerics
df$kept_pct <- as.numeric(df$kept_pct)
df$n_kept   <- as.numeric(df$n_kept)

# Post-trim counts: use provided, else approximate
if (!has_treated_kept || !has_control_kept) {
  warn_once("treated_kept / control_kept not found; approximating via treated_prop × n_kept.")
  df$treated_kept <- df$treated_prop * df$n_kept
  df$control_kept <- (1 - df$treated_prop) * df$n_kept
  df$kept_approx  <- TRUE
} else {
  df$treated_kept <- as.numeric(df$treated_kept)
  df$control_kept <- as.numeric(df$control_kept)
  df$kept_approx  <- FALSE
}

# ---------------- Load GATES lift from TableC (if available) ----------------
gates_map <- NULL
if (file.exists(cate_main_in)) {
  cm <- read.csv(cate_main_in, stringsAsFactors = FALSE, check.names = FALSE)
  # Expected columns: Combo_full, GATES_lift, GATES_bot_mean, GATES_top_mean
  want <- c("Combo_full","GATES_lift","GATES_bot_mean","GATES_top_mean")
  miss_c <- setdiff(want, names(cm))
  if (length(miss_c)) {
    warning("paper_TableC_CATE_main.csv missing columns: ", paste(miss_c, collapse = ", "),
            " — GATES column will be blank.")
  } else {
    cm$Stub <- safe_combo_name(cm$Combo_full)
    gates_map <- list(
      by_combo = list(
        lift = setNames(as.numeric(cm$GATES_lift),       cm$Combo_full),
        bot  = setNames(as.numeric(cm$GATES_bot_mean),   cm$Combo_full),
        top  = setNames(as.numeric(cm$GATES_top_mean),   cm$Combo_full)
      ),
      by_stub = list(
        lift = setNames(as.numeric(cm$GATES_lift),       cm$Stub),
        bot  = setNames(as.numeric(cm$GATES_bot_mean),   cm$Stub),
        top  = setNames(as.numeric(cm$GATES_top_mean),   cm$Stub)
      )
    )
  }
} else {
  warning("Cannot find ", cate_main_in, " — GATES column will be blank.")
}

# ---------------- Filter first, then rank ----------------
df_f <- df[is.finite(df$kept_pct) & df$kept_pct >= keep_min &
             is.finite(df$q_value_bh) & df$q_value_bh <= q_max, , drop = FALSE]

# If nothing passes globally, stop early
if (!nrow(df_f)) stop("No rows pass kept_pct >= ", keep_min, "% and q ≤ ", q_max, ".")

df_f$abs_ATE_pct <- abs(df_f$ATE_pct)

# Top-N per K after filter with order: q asc -> |ATE_pct| desc -> kept_pct desc
topN <- do.call(rbind, lapply(ks_to_plot, function(kv) {
  sub <- df_f[df_f$k == kv, , drop = FALSE]
  if (!nrow(sub)) return(sub)
  sub <- sub[order(sub$q_value_bh, -sub$abs_ATE_pct, -sub$kept_pct), , drop = FALSE]
  sub[seq_len(min(nrow(sub), top_n)), , drop = FALSE]
}))

# ---------------- Build rows (no Pre (C/T); add GATES) ----------------
make_rows_for_k <- function(dat_k, k_val) {
  head_row <- data.frame(
    Subgroup                        = paste0("K = ", k_val),
    `Post (C/T)`                    = "",
    `Kept (n %)`                    = "",
    ` `                             = paste(rep(" ", ci_spaces), collapse = " "),
    `ATE% (95% CI)`                 = "",
    `GATES lift (pp)`  = "",
    `q`                             = "",
    est = NA_real_, low = NA_real_, hi = NA_real_, se = NA_real_,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!nrow(dat_k)) return(head_row)
  
  # Labels
  labs <- dat_k$traits_combination
  if (isTRUE(use_shorten)) labs <- shorten_combo(labs)
  labs <- wrap_combo_vec(labs, width = wrap_width)
  combo <- paste0("   ", labs)
  
  # Counts and text columns
  post_ct   <- paste0(fmt_int(dat_k$control_kept),  " / ", fmt_int(dat_k$treated_kept))
  kept_info <- paste0(fmt_int(dat_k$n_kept), " (", fmt_pct1(dat_k$kept_pct), "%)")
  eff_txt   <- paste0(fmt_pct1(dat_k$ATE_pct), "% [",
                      fmt_pct1(dat_k$CI_lower_pct), ", ",
                      fmt_pct1(dat_k$CI_upper_pct), "%]")
  q_txt     <- fmt_q3(dat_k$q_value_bh)
  
  # GATES mapping (try exact combo first, then stub fallback)
  gates_txt <- rep("", nrow(dat_k))
  if (!is.null(gates_map)) {
    combos <- dat_k$traits_combination
    stubs  <- safe_combo_name(combos)
    
    gl <- unname(gates_map$by_combo$lift[combos]); gb <- unname(gates_map$by_combo$bot[combos]); gt <- unname(gates_map$by_combo$top[combos])
    need_stub <- !(is.finite(gl) & is.finite(gb) & is.finite(gt))
    if (any(need_stub)) {
      gl[need_stub] <- unname(gates_map$by_stub$lift[stubs[need_stub]])
      gb[need_stub] <- unname(gates_map$by_stub$bot[stubs[need_stub]])
      gt[need_stub] <- unname(gates_map$by_stub$top[stubs[need_stub]])
    }
    ok <- is.finite(gl) & is.finite(gb) & is.finite(gt)
    gates_txt[ok] <- paste0(fmt_pct1(gl[ok]), " [", fmt_pct1(gb[ok]), "/", fmt_pct1(gt[ok]), "]")
  }
  
  # --- Box size: robust-normalized (not tiny, not huge) ---
  w1 <- 1 / pmax(1e-9, dat_k$StdErr)                                # precision
  w2 <- sqrt(pmax(0, pmin(dat_k$treated_kept, dat_k$control_kept)))  # support
  w  <- 0.6*w1 + 0.4*w2                                              # blend
  
  # optional stabilization if highly skewed:
  # w <- log1p(w)
  
  qlo <- as.numeric(quantile(w, winsor_lo, na.rm = TRUE))
  qhi <- as.numeric(quantile(w, winsor_hi, na.rm = TRUE))
  w_clamped <- pmin(pmax(w, qlo), qhi)
  
  rng <- max(w_clamped, na.rm = TRUE) - min(w_clamped, na.rm = TRUE)
  if (!is.finite(rng) || rng <= 0) {
    se <- rep((size_floor + size_cap)/2, length(w_clamped))  # fallback: constant
  } else {
    se <- size_floor + (w_clamped - min(w_clamped, na.rm = TRUE)) / rng * (size_cap - size_floor)
  }
  
  det <- data.frame(
    Subgroup                        = combo,
    `Post (C/T)`                    = post_ct,
    `Kept (n %)`                    = kept_info,
    ` `                             = paste(rep(" ", ci_spaces), collapse = " "),
    `ATE% (95% CI)`                 = eff_txt,
    `GATES lift (pp)`  = gates_txt,
    `q`                             = q_txt,
    est = dat_k$ATE_pct, low = dat_k$CI_lower_pct, hi = dat_k$CI_upper_pct, se = se,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  rbind(head_row, det)
}

plot_table <- do.call(rbind, lapply(ks_to_plot, function(kv) {
  make_rows_for_k(topN[topN$k == kv, , drop = FALSE], kv)
}))

# Clean text NAs
txt_cols <- c("Subgroup","Post (C/T)","Kept (n %)"," ","ATE% (95% CI)","GATES lift (pp)","q")
for (cc in txt_cols) if (cc %in% names(plot_table)) plot_table[[cc]][is.na(plot_table[[cc]])] <- ""

# ----- Axis: build breaks first, then set xlim to the breaks' range
rng_vals <- range(c(plot_table$low, plot_table$hi), na.rm = TRUE)
if (!all(is.finite(rng_vals))) rng_vals <- c(-50, 50)
breaks <- pretty(rng_vals, n = 5)
xlim   <- range(breaks)
ticks  <- breaks

# ---------------- Theme (right-aligned columns, no zebra stripes)
n_vis_cols <- 7

thm <- forest_theme(
  base_size   = base_font,
  ci_pch      = 15, ci_col = "#4575b4", ci_fill = "#4575b4",
  ci_alpha    = 0.9, ci_lwd = 1.1, ci_Theight = 0.18,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey40"),
  core = list(
    fg_params = list(
      # Right-align everything EXCEPT the GATES column (6th), which we left-align
      hjust = c(1, 1, 1, 1, 1, 1, 1),
      x     = c(1, 1, 1, 1, 1, 1, 1),
      lineheight = rep(1.05, n_vis_cols)
    ),
    bg_params = list(fill = "white", col = NA)
  ),
  colhead = list(
    fg_params = list(
      # Center headers, but left-align the GATES header
      hjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
      x     = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
      fontface = rep(2, n_vis_cols),
      lineheight = rep(1.05, n_vis_cols)
    ),
    bg_params = list(fill = "white", col = NA)
  ),
  summary_col  = "#333333",
  summary_fill = "#333333"
)


# Bold subgroup rows (K = ...)
is_sum <- is.na(plot_table$est)

# Visible columns (blank column is 5th → ci_column = 5)
visible_cols <- c("Subgroup","Post (C/T)","Kept (n %)"," ","ATE% (95% CI)","GATES lift (pp)", "q")

ci_col_idx <- which(colnames(plot_table[, visible_cols, drop = FALSE]) == " ")

# Draw forest
p <- forest(
  plot_table[, visible_cols],
  est        = plot_table$est,
  lower      = plot_table$low,
  upper      = plot_table$hi,
  sizes      = plot_table$se,            # robust-normalized size
  ci_column  = ci_col_idx,                        # plot CI in the blank " " column
  ref_line   = 0,
  xlim       = xlim,
  ticks_at   = ticks,
  is_summary = is_sum,                   # bold K rows
  arrow_lab  = c("Lower price", "Higher price"),
  theme      = thm,
  # footnote   = if (any(df$kept_approx)) {
  #   "Note: Post-trim counts may be approximations based on treated_prop × n_kept. GATES lift = Top−Bottom decile mean CATE (pp)."
  # } else { "GATES lift = Top−Bottom decile mean CATE (pp)." }
)

# --- Header rule (horizontal line under column names) ---
ncols <- length(visible_cols)
p <- add_border(
  p,
  part  = "header",
  row   = 1,
  col   = 1:ncols,
  where = "bottom",
  gp    = gpar(lwd = 1.2, col = "grey40")
)

# -------- Save PDF --------
pdf(outfile_pdf, width = 18, height = 10)
plot(p)
dev.off()
message("Saved PDF: ", normalizePath(outfile_pdf))

# -------- Save PNG (high-res via ragg) --------
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
ragg::agg_png(outfile_png, width = 15, height = 7, units = "in", res = 900)
plot(p); dev.off()
message("Saved PNG: ", normalizePath(outfile_png))
