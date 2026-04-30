# R script: compact LaTeX table (exactly 6 columns)
# INPUT/OUTPUT
infile  <- "m2_cate_numeric_top3.csv"          # adjust path if needed
outfile <- "m2_cate_numeric_top3_grouped_table.tex"

# -------- Read
df <- read.csv(infile, stringsAsFactors = FALSE)

# Ensure needed columns; add fallbacks
if (!"abs_smd" %in% names(df)) {
  if ("SMD_kept" %in% names(df)) df$abs_smd <- abs(df$SMD_kept) else df$abs_smd <- NA_real_
}
need <- c("treatment","covariate","slope_catepct_per_unit","abs_smd","q_value_bh","r2_tau")
miss <- setdiff(need, names(df))
if (length(miss)) df[, miss] <- NA

# -------- Group + subcategory from treatment
grp <- ifelse(grepl("^buyer_N_1_",  df$treatment), "Buyer N-1 High",
              ifelse(grepl("^buyer_N_",    df$treatment), "Buyer N High",
                     ifelse(grepl("^seller_N_1_", df$treatment), "Seller N-1 High",
                            ifelse(grepl("^seller_N_",   df$treatment), "Seller N High", "Other"))))
subcat <- sub("_high$", "", sub("^((buyer|seller)_N(_1)?_)", "", df$treatment))

# -------- Escape underscores as "\_ "
esc_us <- function(x) gsub("_", "\\\\_ ", x, fixed = TRUE)
df$group       <- grp
df$subcategory <- esc_us(subcat)
df$covariate   <- esc_us(df$covariate)

# -------- Formatting
fmt3 <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 3))
# p-like formatter: if < 0.001 show 0.000 (per request)
fmtp <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "",
         ifelse(x < 0.001, "0.000", formatC(x, format = "f", digits = 3)))
}

# -------- Ordering: by q-value asc, then R2_tau desc
df$q_value_bh <- suppressWarnings(as.numeric(df$q_value_bh))
df$r2_tau     <- suppressWarnings(as.numeric(df$r2_tau))

group_order <- c("Buyer N-1 High", "Buyer N High", "Seller N-1 High", "Seller N High")
df$group <- factor(df$group, levels = group_order)

qsort  <- ifelse(is.na(df$q_value_bh), Inf, df$q_value_bh)
r2sort <- -ifelse(is.na(df$r2_tau), -Inf, df$r2_tau)  # negative for descending
ord <- order(df$group, df$subcategory, qsort, r2sort, df$covariate)
df  <- df[ord, ]

# -------- Build LaTeX rows
# Columns (6 total): Treatment | Covariate | Slope(CATE%/unit) | |SMD| | BLP q-value | R^2_tau
rows <- character(0)
for (g in group_order[group_order %in% df$group]) {
  dfg <- df[df$group == g, , drop = FALSE]
  if (nrow(dfg) == 0) next
  
  # group header + spacer (6 cols -> 5 ampersands)
  rows <- c(rows, sprintf("%s &  &  &  &  & \\\\", g))
  rows <- c(rows, " &  &  &  &  & \\\\")
  
  for (sc in unique(dfg$subcategory)) {
    dfs <- dfg[dfg$subcategory == sc, , drop = FALSE]
    for (i in seq_len(nrow(dfs))) {
      treat_col <- if (i == 1) dfs$subcategory[i] else ""
      rows <- c(rows, sprintf(
        "%s & %s & %s & %s & %s & %s\\\\",
        treat_col,
        dfs$covariate[i],
        fmt3(dfs$slope_catepct_per_unit[i]),
        fmt3(dfs$abs_smd[i]),
        fmtp(dfs$q_value_bh[i]),
        fmt3(dfs$r2_tau[i])
      ))
    }
  }
}

# -------- Assemble LaTeX (6 columns)
latex <- c(
  "\\begin{table}",
  "\\caption{Top covariates per treatment (numeric): slope of CATE (\\% per unit), balance, BLP q-value (BH; values < 0.001 shown as 0.000), and heterogeneity $R^2_{\\tau}$.}",
  "\\label{tab:m2_cate_numeric_top5}",
  "\\begin{tabular}{llcccc}",
  "\\toprule",
  "Treatment & Covariate & Slope $\\mathrm{CATE}_{\\%}$/unit & $|\\mathrm{SMD}|$ & BLP q-value & $R^2_{\\tau}$\\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex, outfile)
cat("Wrote:", outfile, "\n")
