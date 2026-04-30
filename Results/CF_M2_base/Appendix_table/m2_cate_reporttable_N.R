# R script: build a LaTeX table from m2_cate_numeric_top5.csv with underscores escaped as "\_ "
# Input/Output
infile  <- "m2_cate_numeric_top5.csv"      # change to your path
outfile <- "m2_cate_numeric_top5_grouped_table.tex"

# Read
df <- read.csv(infile, stringsAsFactors = FALSE)

# Keep needed columns
df <- df[, c("treatment", "covariate", "cor", "slope_catepct_per_unit", "SMD_kept")]

# Group label from treatment
grp <- ifelse(grepl("^buyer_N_1_",  df$treatment), "Buyer N-1 High",
              ifelse(grepl("^buyer_N_",    df$treatment), "Buyer N High",
                     ifelse(grepl("^seller_N_1_", df$treatment), "Seller N-1 High",
                            ifelse(grepl("^seller_N_",   df$treatment), "Seller N High", "Other"))))

# Subcategory from treatment (e.g., gaslimit_cost, gaspaid, total_value)
subcat <- sub("_high$", "", sub("^((buyer|seller)_N(_1)?_)", "", df$treatment))

# Escape underscores as "\_ " per requirement
esc_us <- function(x) gsub("_", "\\\\_ ", x, fixed = TRUE)

df$group      <- grp
df$subcategory<- subcat
df$treatment  <- esc_us(df$treatment)       # not used directly but kept for reference
df$covariate  <- esc_us(df$covariate)
df$subcategory<- esc_us(df$subcategory)

# Number formatting (3 decimals)
fmt3 <- function(x) formatC(x, format = "f", digits = 3)

# Desired group order
group_order <- c("Buyer N-1 High", "Buyer N High", "Seller N-1 High", "Seller N High")
df$group <- factor(df$group, levels = group_order)

# Sort within group by subcategory then covariate (stable)
ord <- order(df$group, df$subcategory, df$covariate)
df  <- df[ord, ]

# Build LaTeX rows with hierarchical blocks:
# - Group header row: "<Group> & & & & \\"
# - Blank row: "& & & & \\"
# - For each subcategory: first row shows subcategory in Treatment col; subsequent rows leave it blank
rows <- character(0)
for (g in group_order[group_order %in% df$group]) {
  dfg <- df[df$group == g, ]
  if (nrow(dfg) == 0) next
  
  rows <- c(rows, sprintf("%s &  &  &  & \\\\", g))
  rows <- c(rows, " &  &  &  & \\\\")
  
  for (sc in unique(dfg$subcategory)) {
    dfs <- dfg[dfg$subcategory == sc, ]
    for (i in seq_len(nrow(dfs))) {
      treat_col <- if (i == 1) dfs$subcategory[i] else ""
      rows <- c(rows, sprintf(
        "%s & %s & %s & %s & %s\\\\",
        treat_col,
        dfs$covariate[i],
        fmt3(dfs$cor[i]),
        fmt3(dfs$slope_catepct_per_unit[i]),
        fmt3(dfs$SMD_kept[i])
      ))
    }
  }
}

# Assemble minimal LaTeX table ONLY using the specified structure
latex <- c(
  "\\begin{table}",
  "\\centering",
  "\\begin{tabular}{ccccc}",
  "Treatment & Covariate  & Cor  & Slope $\\text{CATE}_{pct}$ Per Unit& SMD kept\\\\",
  rows,
  "\\end{tabular}",
  "\\caption{Caption}",
  "\\label{tab:placeholder}",
  "\\end{table}"
)

writeLines(latex, outfile)
