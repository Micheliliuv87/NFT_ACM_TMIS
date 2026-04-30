# ---- CF M2 V7 Up — Forest Plot (ATE/GATES + BLP badges + R²(τ)) ----
# Reads:
#   1) cf_M2_V7_up_treatment_results_with_q.csv   (ATE/GATES by treatment; from CF pipeline)
#   2) cf_M2_V7_up_BLP_per_covariate.csv          (BLP per covariate; from BLP step)  [optional]
#   3) cf_M2_V7_up_R2tau_by_treatment.csv         (R2_tau per treatment)              [optional fallback]
#
# Output: cf_M2_V7_up_forest.pdf + .png
# Requires: forestploter, grid, scales, data.table, dplyr (and ragg for PNG)

suppressPackageStartupMessages({
  library(forestploter)
  library(grid)
  library(scales)
  library(data.table)
  library(dplyr)
})

# ---------- Encoding & font ----------
try(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")), silent = TRUE)

# Prefer a Libertine family (different OSes ship different names).
# If none is found, fall back to a Unicode-safe font.
font_candidates <- c(      # modern fork, common on many systems
  "Linux Libertine O",     # TeX/OTF name
  # older name
  "Libertine")             # generic
font_family <- "DejaVu Sans"                  # fallback (has Greek τ)

if (requireNamespace("systemfonts", quietly = TRUE)) {
  fams <- unique(systemfonts::system_fonts()$family)
  hit  <- intersect(font_candidates, fams)
  if (length(hit)) font_family <- hit[1]
}
message("Using font family: ", font_family)



# ---------------- Config ----------------
infile_ate  <- "cf_M2_V7_up_treatment_results_with_q.csv"
infile_blp  <- "cf_M2_V7_up_BLP_per_covariate.csv"
infile_r2t  <- "cf_M2_V7_up_R2tau_by_treatment.csv"

outfile_pdf <- "cf_M2_V7_up_forest.pdf"
outfile_png <- "cf_M2_V7_up_forest.png"

ci_spaces   <- 30   # width of CI lane (spaces)
wrap_width  <- 70
base_font   <- 12

# Robust box-size normalization
size_floor <- 0.5
size_cap   <- 0.9
winsor_lo  <- 0.15
winsor_hi  <- 0.85

# BLP significance cutoff shown in the table badge
Q_BLP_CUTOFF <- 0.10

# Pretty mathy label for R2tau
epxres <- expression(R^2 * tau)
# make the expression to display as a string for later use
epxres_str <- "R²(τ)"
R2TAU_LABEL <- epxres_str     # "R²(τ)"
R2TAU_COL   <- paste0(R2TAU_LABEL, " (%)")

# ---------------- Helpers ----------------
fmt_int  <- function(x) format(round(as.numeric(x)), big.mark = ",", trim = TRUE)
fmt_pct1 <- function(x) sprintf("%.1f", as.numeric(x))
fmt_q3   <- function(x) ifelse(is.finite(as.numeric(x)), sprintf("%.3f", as.numeric(x)), "")
wrap_text <- function(x, width = 30) {
  vapply(x, function(s) if (is.na(s) || s == "") "" else paste(strwrap(s, width), collapse = "\n"), "")
}

metric_short <- function(m) {
  m <- as.character(m)
  m[m == "Total Value"]    <- "Total Value"
  m[m == "Gas Paid"]       <- "Gas Paid"
  m[m == "Gas Limit Cost"] <- "Gas Limit Cost"
  m
}

# Subgroup ordering for headers
subgroup_order_map <- function(s) {
  ifelse(s == "Seller N High",   1L,
         ifelse(s == "Seller N-1 High", 2L,
                ifelse(s == "Buyer N High",    3L,
                       ifelse(s == "Buyer N-1 High",  4L, 99L))))
}
metric_order_map <- function(m) {
  ifelse(m == "Total Value", 1L,
         ifelse(m == "Gas Paid", 2L,
                ifelse(m == "Gas Limit Cost", 3L, 99L)))
}

# Vectorized badge formatter: "23/253 (9.1%)"
badge_str <- function(sig, tot) {
  sig <- as.numeric(sig); tot <- as.numeric(tot)
  pct <- 100 * sig / pmax(tot, 1)
  out <- sprintf("%d/%d (%s%%)", sig, tot, fmt_pct1(pct))
  out[!is.finite(sig) | !is.finite(tot) | tot <= 0] <- ""
  out
}

# ---------------- Load ATE/GATES data ----------------
if (!file.exists(infile_ate)) stop("Cannot find: ", infile_ate)
df <- data.table::fread(infile_ate)

need <- c("treatment","display_label","subgroup","metric",
          "n_kept","kept_pct","treated_prop","treated_kept","control_kept",
          "ATE","CI_lower","CI_upper","StdErr",
          "ATE_pct","CI_lower_pct","CI_upper_pct",
          "p_value","q_value_bh",
          "GATES_lift","GATES_bot_mean","GATES_top_mean")
miss <- setdiff(need, names(df))
if (length(miss)) stop("Missing columns in input: ", paste(miss, collapse = ", "))

# Coerce numerics safely
num_cols <- c("n_kept","kept_pct","treated_prop","treated_kept","control_kept",
              "ATE","CI_lower","CI_upper","StdErr",
              "ATE_pct","CI_lower_pct","CI_upper_pct","p_value","q_value_bh",
              "GATES_lift","GATES_bot_mean","GATES_top_mean")
for (cc in intersect(num_cols, names(df))) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))

# Order (report everything)
df$subgroup_ord <- subgroup_order_map(df$subgroup)
df$metric_ord   <- metric_order_map(df$metric)
data.table::setorder(df, subgroup_ord, metric_ord, treatment)

# ---------------- Load BLP (optional) and compute per-treatment badges + R²(τ) ----------------
treats <- unique(df$treatment)
blp_badge_raw <- setNames(rep("", length(treats)), treats)
r2tau_map_pct <- setNames(rep(NA_real_, length(treats)), treats)

if (file.exists(infile_blp)) {
  blp_raw <- data.table::fread(infile_blp)
  # normalize names
  names(blp_raw) <- tolower(gsub("\\s+", "_", names(blp_raw)))
  if (!("treatment" %in% names(blp_raw)) && "treat" %in% names(blp_raw)) {
    setnames(blp_raw, "treat", "treatment")
  }
  if (!("covariate" %in% names(blp_raw))) {
    if ("variable" %in% names(blp_raw)) setnames(blp_raw, "variable", "covariate")
    else blp_raw$covariate <- NA_character_
  }
  blp_raw$treatment <- as.character(blp_raw$treatment)
  blp_raw$covariate <- as.character(blp_raw$covariate)
  
  # Identify q-value column; if absent, derive from p-values
  qcol <- intersect(names(blp_raw), c("q_blp","q","q_value_bh","q_bh"))
  if (!length(qcol)) {
    pcol <- intersect(names(blp_raw), c("p_blp","p_value","p"))
    if (length(pcol)) {
      blp_raw$q_blp <- p.adjust(blp_raw[[pcol[1]]], method = "BH")
      qcol <- "q_blp"
    }
  }
  
  # Per-treatment R²(τ) (if provided inside the BLP file)
  r2cand <- intersect(names(blp_raw), c("r2_tau","r2t","r2_hetero","r2"))
  if (length(r2cand)) {
    r2_by_treat <- blp_raw %>%
      group_by(treatment) %>%
      summarise(r2_tau = suppressWarnings(max(.data[[r2cand[1]]], na.rm = TRUE)),
                .groups = "drop")
    idx <- match(names(r2tau_map_pct), r2_by_treat$treatment)
    pick <- !is.na(idx)
    r2tau_map_pct[pick] <- 100 * r2_by_treat$r2_tau[idx[pick]]
  } else if (file.exists(infile_r2t)) {
    r2t <- data.table::fread(infile_r2t)
    names(r2t) <- tolower(gsub("\\s+", "_", names(r2t)))
    if (!("r2_tau" %in% names(r2t))) {
      cand <- intersect(names(r2t), c("r2tau","r2","r2_hetero"))
      if (length(cand)) setnames(r2t, cand[1], "r2_tau")
    }
    if (all(c("treatment","r2_tau") %in% names(r2t))) {
      r2t$treatment <- as.character(r2t$treatment)
      idx2 <- match(names(r2tau_map_pct), r2t$treatment)
      pick2 <- !is.na(idx2)
      r2tau_map_pct[pick2] <- 100 * r2t$r2_tau[idx2[pick2]]
    }
  }
  
  # Print overall R²(τ) summary for copy/paste (not plotted)
  if (any(is.finite(r2tau_map_pct))) {
    r2vals <- r2tau_map_pct[is.finite(r2tau_map_pct)]
    cat(
      sprintf("\nOverall %s (median [min,max]): %s%% [%s%%, %s%%]\n",
              R2TAU_LABEL,
              fmt_pct1(median(r2vals)), fmt_pct1(min(r2vals)), fmt_pct1(max(r2vals)))
    )
  }
  
  # Badges (RAW only)
  if (length(qcol)) {
    blp_sum <- blp_raw %>%
      filter(!is.na(.data[[qcol[1]]])) %>%
      group_by(treatment) %>%
      summarise(
        blp_sig = sum(.data[[qcol[1]]] < Q_BLP_CUTOFF, na.rm = TRUE),
        blp_tot = n(),
        .groups = "drop"
      )
    blp_badge_raw[blp_sum$treatment] <- badge_str(blp_sum$blp_sig, blp_sum$blp_tot)
  }
}

# ---------------- Build rows (subgroup headers + detail rows) ----------------
make_rows_for_subgroup <- function(dat_sg, sg_name, blp_badge_raw, r2tau_map_pct) {
  head_row <- data.frame(
    Subgroup              = sg_name,
    `Post (C/T)`          = "",
    `Kept (n %)`          = "",
    ` `                   = paste(rep(" ", ci_spaces), collapse = " "),
    `ATE% (95% CI)`       = "",
    `GATES lift (pp)`     = "",
    `BLP q<0.10`          = "",
    `q`                   = "",
    est = NA_real_, low = NA_real_, hi = NA_real_, se = NA_real_,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  head_row[[R2TAU_COL]] <- ""  # add pretty R²(τ) column
  
  if (!nrow(dat_sg)) return(head_row)
  
  # Label: show metric under subgroup (concise)
  lab_metric <- metric_short(dat_sg$metric)
  combo <- paste0("   ", wrap_text(lab_metric, width = wrap_width))
  
  # Text columns
  post_ct   <- paste0(fmt_int(dat_sg$control_kept), " / ", fmt_int(dat_sg$treated_kept))
  kept_info <- paste0(fmt_int(dat_sg$n_kept), " (", fmt_pct1(dat_sg$kept_pct), "%)")
  eff_txt   <- paste0(fmt_pct1(dat_sg$ATE_pct), "% [",
                      fmt_pct1(dat_sg$CI_lower_pct), ", ",
                      fmt_pct1(dat_sg$CI_upper_pct), "%]")
  q_txt     <- fmt_q3(dat_sg$q_value_bh)
  
  # GATES text
  gates_txt <- ifelse(
    is.finite(dat_sg$GATES_lift) & is.finite(dat_sg$GATES_bot_mean) & is.finite(dat_sg$GATES_top_mean),
    paste0(fmt_pct1(dat_sg$GATES_lift), " [",
           fmt_pct1(dat_sg$GATES_bot_mean), "/", fmt_pct1(dat_sg$GATES_top_mean), "]"),
    ""
  )
  
  # --- Box size: precision + support, winsorized, then scaled ---
  w1 <- 1 / pmax(1e-9, dat_sg$StdErr)
  w2 <- sqrt(pmax(0, pmin(dat_sg$treated_kept, dat_sg$control_kept)))
  w  <- 0.6*w1 + 0.4*w2
  qlo <- as.numeric(quantile(w, winsor_lo, na.rm = TRUE))
  qhi <- as.numeric(quantile(w, winsor_hi, na.rm = TRUE))
  w_clamped <- pmin(pmax(w, qlo), qhi)
  rng <- max(w_clamped, na.rm = TRUE) - min(w_clamped, na.rm = TRUE)
  if (!is.finite(rng) || rng <= 0) {
    se <- rep((size_floor + size_cap)/2, length(w_clamped))
  } else {
    se <- size_floor + (w_clamped - min(w_clamped, na.rm = TRUE)) / rng * (size_cap - size_floor)
  }
  
  # Vectorized badges & R²(τ) by row
  tnms <- as.character(dat_sg$treatment)
  blp_raw_txt  <- blp_badge_raw[tnms]; blp_raw_txt[is.na(blp_raw_txt)] <- ""
  r2vals <- r2tau_map_pct[tnms]
  r2tau_txt <- ifelse(is.finite(r2vals), fmt_pct1(r2vals), "")
  
  det <- data.frame(
    Subgroup              = combo,
    `Post (C/T)`          = post_ct,
    `Kept (n %)`          = kept_info,
    ` `                   = paste(rep(" ", ci_spaces), collapse = " "),
    `ATE% (95% CI)`       = eff_txt,
    `GATES lift (pp)`     = gates_txt,
    `BLP q<0.10`          = blp_raw_txt,
    `q`                   = q_txt,
    est = dat_sg$ATE_pct, low = dat_sg$CI_lower_pct, hi = dat_sg$CI_upper_pct, se = se,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  det[[R2TAU_COL]] <- r2tau_txt   # fill pretty R²(τ) column
  
  rbind(head_row, det)
}

plot_table <- do.call(
  rbind,
  lapply(split(df, df$subgroup), function(d) {
    make_rows_for_subgroup(d, unique(d$subgroup)[1], blp_badge_raw, r2tau_map_pct)
  })
)

# Clean text NAs
txt_cols <- c("Subgroup","Post (C/T)","Kept (n %)"," ",
              "ATE% (95% CI)","GATES lift (pp)","BLP q<0.10",
              R2TAU_COL,"q")
for (cc in intersect(txt_cols, names(plot_table))) {
  plot_table[[cc]][is.na(plot_table[[cc]])] <- ""
}

# Axis from CI range on % scale
rng_vals <- range(c(plot_table$low, plot_table$hi), na.rm = TRUE)
if (!all(is.finite(rng_vals))) rng_vals <- c(-50, 50)
breaks <- pretty(rng_vals, n = 5)
xlim   <- range(breaks)
ticks  <- breaks

# ---------------- Theme ----------------
bold_names <- c("Seller N High","Seller N-1 High","Buyer N High","Buyer N-1 High")
is_sum <- plot_table$Subgroup %in% bold_names  # no overall R²(τ) header row

# Visible columns (blank column is 4th → ci_column = 4)
visible_cols <- c("Subgroup","Post (C/T)","Kept (n %)"," ",
                  "ATE% (95% CI)","GATES lift (pp)","BLP q<0.10",
                  R2TAU_COL,"q")
ci_col_idx <- which(colnames(plot_table[, visible_cols, drop = FALSE]) == " ")
n_vis_cols <- length(visible_cols)

thm <- forest_theme(
  base_size   = base_font,
  ci_pch      = 15, ci_col = "#4575b4", ci_fill = "#4575b4",
  ci_alpha    = 0.9, ci_lwd = 1.1, ci_Theight = 0.18,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey40"),
  core = list(
    fg_params = list(
      hjust = rep(1, n_vis_cols),
      x     = rep(1, n_vis_cols),
      lineheight = rep(1.05, n_vis_cols)
    ),
    bg_params = list(fill = "white", col = NA)
  ),
  colhead = list(
    fg_params = list(
      hjust = rep(0.5, n_vis_cols),
      x     = rep(0.6, n_vis_cols),
      fontface = rep(2, n_vis_cols),
      lineheight = rep(1.05, n_vis_cols)
    ),
    bg_params = list(fill = "white", col = NA)
  ),
  summary_col  = "#333333",
  summary_fill = "#333333"
)

# Draw forest
p <- forest(
  plot_table[, visible_cols],
  est        = plot_table$est,
  lower      = plot_table$low,
  upper      = plot_table$hi,
  sizes      = plot_table$se,
  ci_column  = ci_col_idx,
  ref_line   = 0,
  xlim       = xlim,
  ticks_at   = ticks,
  is_summary = is_sum,
  arrow_lab  = c("Lower price", "Higher price"),
  theme      = thm
)

# Header rule
ncols <- length(visible_cols)
p <- add_border(
  p, part = "header", row = 1, col = 1:ncols, where = "bottom",
  gp = gpar(lwd = 1.2, col = "grey40")
)

# -------- Save PDF --------
pdf(outfile_pdf, width = 18, height = 10)
plot(p)
dev.off()
message("Saved PDF: ", normalizePath(outfile_pdf))

# -------- Save PNG (high-res via ragg) --------
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
ragg::agg_png(outfile_png, width = 13, height = 4, units = "in", res = 900)
plot(p); dev.off()
message("Saved PNG: ", normalizePath(outfile_png))
