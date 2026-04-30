# R script: build a LaTeX table from a *binary* CSV where columns are:
# treatment, covariate, n0, n1, mean0, mean1, difference, SMD_kept
# Underscores in LaTeX fields are escaped as "\_ " (backslash-underscore-space).
# Output uses ONLY the minimal table structure requested.

infile  <- "m2_cate_binary_top5.csv"        # <-- change to your path
outfile <- "m2_cate_binary_grouped_table.tex"

# ---------- Read & select ----------
df <- read.csv(infile, stringsAsFactors = FALSE)

# keep/rename to expected names (in case of casing differences)
need <- c("treatment","covariate","n0","n1","mean0","mean1","difference","SMD_kept")
stopifnot(all(need %in% names(df)))
df <- df[need]

# ---------- Grouping ----------
grp <- ifelse(grepl("^buyer_N_1_",  df$treatment), "Buyer N-1 High",
              ifelse(grepl("^buyer_N_",    df$treatment), "Buyer N High",
                     ifelse(grepl("^seller_N_1_", df$treatment), "Seller N-1 High",
                            ifelse(grepl("^seller_N_",   df$treatment), "Seller N High", "Other"))))

# subcategory from treatment (e.g., gaslimit_cost, gaspaid, total_value)
subcat <- sub("_high$", "", sub("^((buyer|seller)_N(_1)?_)", "", df$treatment))

# ---------- Escape underscores as "\_ " ----------
esc_us <- function(x) gsub("_", "\\\\_ ", x, fixed = TRUE)

df$group       <- grp
df$subcategory <- esc_us(subcat)
df$covariate   <- esc_us(df$covariate)

# ---------- Numeric formatting ----------
fmt0 <- function(x) formatC(as.numeric(x), format = "f", digits = 0)  # counts
fmt3 <- function(x) formatC(as.numeric(x), format = "f", digits = 3)  # decimals

# ---------- Order ----------
group_order <- c("Buyer N-1 High", "Buyer N High", "Seller N-1 High", "Seller N High", "Other")
df$group <- factor(df$group, levels = group_order)
ord <- order(df$group, df$subcategory, df$covariate)
df  <- df[ord, ]

# ---------- Build LaTeX rows ----------
rows <- character(0)
for (g in group_order[group_order %in% df$group]) {
  dfg <- df[df$group == g, ]
  if (nrow(dfg) == 0) next
  
  # Group header + blank line (8 columns => 7 &'s per row)
  rows <- c(rows, sprintf("%s &  &  &  &  &  &  & \\\\", g))
  rows <- c(rows, " &  &  &  &  &  &  & \\\\")
  
  # Within group, iterate subcategories
  for (sc in unique(dfg$subcategory)) {
    dfs <- dfg[dfg$subcategory == sc, ]
    for (i in seq_len(nrow(dfs))) {
      treat_col <- if (i == 1) dfs$subcategory[i] else ""
      rows <- c(rows, sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s\\\\",
        treat_col,
        dfs$covariate[i],
        fmt0(dfs$n0[i]),
        fmt0(dfs$n1[i]),
        fmt3(dfs$mean0[i]),
        fmt3(dfs$mean1[i]),
        fmt3(dfs$difference[i]),
        fmt3(dfs$SMD_kept[i])
      ))
    }
  }
}

# ---------- Assemble *minimal* LaTeX table ----------
latex <- c(
  "\\begin{table}",
  "\\centering",
  "\\begin{tabular}{cccccccc}",
  "Treatment & Covariate & n0 & n1 & mean0 & mean1 & Difference & SMD kept\\\\",
  rows,
  "\\end{tabular}",
  "\\caption{Caption}",
  "\\label{tab:placeholder}",
  "\\end{table}"
)

writeLines(latex, outfile)
