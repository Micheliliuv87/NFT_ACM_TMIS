# R script: compact LaTeX table (binary) with <= 7 cols total
# INPUT/OUTPUT
infile  <- "m2_cate_binary_top3.csv"
outfile <- "m2_cate_binary_grouped_table.tex"

# ---- Read
df <- read.csv(infile, stringsAsFactors = FALSE)

# ---- Ensure/derive useful columns
if (!"abs_smd" %in% names(df) && "SMD_kept" %in% names(df)) df$abs_smd <- abs(df$SMD_kept)
if (!"n_obs" %in% names(df) && all(c("n0","n1") %in% names(df))) df$n_obs <- df$n0 + df$n1
if (!"q_value_bh" %in% names(df) && "p_value" %in% names(df)) df$q_value_bh <- df$p_value

# ---- Minimal required columns must exist
need_min <- c("treatment","covariate")
stopifnot(all(need_min %in% names(df)))

# ---- Grouping & subcategory
grp <- ifelse(grepl("^buyer_N_1_",  df$treatment), "Buyer N-1 High",
              ifelse(grepl("^buyer_N_",    df$treatment), "Buyer N High",
                     ifelse(grepl("^seller_N_1_", df$treatment), "Seller N-1 High",
                            ifelse(grepl("^seller_N_",   df$treatment), "Seller N High", "Other"))))
subcat <- sub("_high$", "", sub("^((buyer|seller)_N(_1)?_)", "", df$treatment))

# ---- Escape underscores as "\_ "
esc_us <- function(x) gsub("_", "\\\\_ ", x, fixed = TRUE)
df$group       <- grp
df$subcategory <- esc_us(subcat)
df$covariate   <- esc_us(df$covariate)

# ---- Build combined counts column n(0,1) = "(n0, n1)"
fmt0_raw <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(is.na(x), "", formatC(x, format="f", digits=0)) }
df$n01 <- paste0("(", fmt0_raw(df$n0), ", ", fmt0_raw(df$n1), ")")

# ---- Formatters
fmt3 <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(is.na(x), "", formatC(x, format="f", digits=3)) }
# show 0.000 when < 0.001
fmtp <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", ifelse(x < 0.001, "0.000", formatC(x, format="f", digits=3)))
}

# ---- Choose extra columns (5 max): n01 replaces n0/n1; include abs_smd
EXTRA_COLS <- c("n01","difference","abs_smd","q_value_bh","r2_tau")  # total columns = 2 + 5 = 7
stopifnot(length(EXTRA_COLS) <= 5)
missing_extras <- setdiff(EXTRA_COLS, names(df))
if (length(missing_extras)) for (m in missing_extras) df[[m]] <- NA

# Map column -> header and formatter
header_map <- c(
  n01="n (0,1)", difference="Difference",
  abs_smd="$|\\mathrm{SMD}|$", q_value_bh="BLP q-value",
  r2_tau="$R^2_{\\tau}$"
)
fmt_for <- function(col, v) switch(col,
                                   n01=as.character(v), difference=fmt3(v), abs_smd=fmt3(v),
                                   q_value_bh=fmtp(v),  r2_tau=fmt3(v),
                                   fmt3(v)
)

# ---- Ordering: by q-value asc, then R2_tau desc
df$q_value_bh <- suppressWarnings(as.numeric(df$q_value_bh))
df$r2_tau     <- suppressWarnings(as.numeric(df$r2_tau))
group_order   <- c("Buyer N-1 High", "Buyer N High", "Seller N-1 High", "Seller N High", "Other")
df$group      <- factor(df$group, levels = group_order)

qsort  <- ifelse(is.na(df$q_value_bh),
                 if ("p_value" %in% names(df)) suppressWarnings(as.numeric(df$p_value)) else Inf,
                 df$q_value_bh)
r2sort <- -ifelse(is.na(df$r2_tau), -Inf, df$r2_tau)  # negative for descending
ord <- order(df$group, df$subcategory, qsort, r2sort, df$covariate)
df  <- df[ord, , drop = FALSE]

# ---- Build LaTeX rows
Ncols <- 2 + length(EXTRA_COLS)
mkrow <- function(cells) paste(paste(cells, collapse=" & "), "\\\\")
rows <- character(0)

for (g in group_order[group_order %in% df$group]) {
  dfg <- df[df$group == g, , drop = FALSE]
  if (!nrow(dfg)) next
  
  # Group header + spacer
  rows <- c(rows, mkrow(c(g, rep("", Ncols - 1))))
  rows <- c(rows, mkrow(rep("", Ncols)))
  
  for (sc in unique(dfg$subcategory)) {
    dfs <- dfg[dfg$subcategory == sc, , drop = FALSE]
    for (i in seq_len(nrow(dfs))) {
      treat_col <- if (i == 1) dfs$subcategory[i] else ""
      # >>> FIX: index each extra by row i so the formatter returns length-1 <<<
      extras <- vapply(EXTRA_COLS, function(cc) fmt_for(cc, dfs[[cc]][i]), character(1))
      rows <- c(rows, mkrow(c(treat_col, dfs$covariate[i], extras)))
    }
  }
}

# ---- Headers + table spec
extra_headers <- vapply(EXTRA_COLS, function(cc) if (!is.null(header_map[[cc]])) header_map[[cc]] else cc, character(1))
colspec <- paste(rep("c", Ncols), collapse = "")
header_line <- mkrow(c("Treatment", "Covariate", extra_headers))

latex <- c(
  "\\begin{table}",
  "\\caption{Binary covariates by treatment: counts $n(0,1)$, difference in CATE (\\% pts), adjusted BLP q-value (values < 0.001 shown as 0.000), and heterogeneity $R^2_{\\tau}$.}",
  "\\label{tab:m2_cate_binary}",
  paste0("\\begin{tabular}{", colspec, "}"),
  "\\toprule",
  header_line,
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex, outfile)
cat("Wrote:", outfile, "\n")
