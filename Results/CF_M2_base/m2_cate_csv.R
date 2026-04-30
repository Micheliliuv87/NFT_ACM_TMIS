suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# ---------------- params ----------------
top_k <- 3

# ---------------- inputs ----------------
df_1 <- read.csv("cf_M2_V7_up_cate_by_covariate_overview.csv")
df_2 <- read.csv("cf_M2_V7_up_covariate_smd_and_cate.csv")

# BLP per covariate (already has r2_tau joined in your main script)
# expected cols: treatment, covariate, p_value, coef_pct, n_obs, method, r2_tau, r2_tau_oos
blp <- read.csv("cf_M2_V7_up_BLP_per_covariate.csv") %>%
  select(treatment, covariate, p_value, coef_pct, n_obs, method,
         r2_tau, r2_tau_oos)

# q-values per treatment (BH)
blp <- blp %>%
  group_by(treatment) %>%
  mutate(q_value_bh = if (all(is.na(p_value))) NA_real_ else p.adjust(p_value, method = "BH")) %>%
  ungroup()

# ---------------- merge + filter ----------------
# Only bring SMD_kept from df_2
df <- merge(
  df_1,
  df_2[, c("treatment", "covariate", "SMD_kept")],
  by = c("treatment", "covariate"),
  all.x = TRUE
)

# drop time FE + the sold_after_30d flag covariate
df <- df %>%
  filter(!str_detect(covariate, "^`?time_")) %>%   # handles optional backtick
  filter(!str_detect(covariate, "^sold_after_30d"))

# attach BLP & R2_tau to every treatment×covariate row
df <- df %>%
  left_join(blp, by = c("treatment", "covariate")) %>%
  mutate(abs_smd = abs(SMD_kept))

# ---------------- numeric table ----------------
# keep numeric rows + relevant columns
df_numeric <- df %>%
  filter(type == "numeric") %>%
  select(
    treatment, covariate,
    cor, slope_catepct_per_unit,      # descriptive (from CATE-by-cov overview)
    SMD_kept, abs_smd,
    # BLP fields:
    coef_pct, p_value, q_value_bh, n_obs, method,
    r2_tau, r2_tau_oos
  )

# rank within treatment: q-value ↑, R2_tau ↓, abs_smd ↑ (better balance first → smaller), |slope| ↓ as tiebreaker
df_numeric_topk <- df_numeric %>%
  group_by(treatment) %>%
  arrange(
    is.na(q_value_bh), q_value_bh,          # NA last, then smallest q first
    desc(r2_tau),                            # larger heterogeneity explained first
    abs_smd,                                 # better balance (smaller SMD) first
    desc(abs(slope_catepct_per_unit))        # stronger descriptive slope as final tie-breaker
  ) %>%
  slice_head(n = top_k) %>%
  ungroup()

# ---------------- binary table ----------------
# build difference column (mean1 - mean0) and keep relevant columns
df_binary <- df %>%
  filter(type == "binary") %>%
  mutate(difference = diff1minus0) %>%
  select(
    treatment, covariate,
    n0, n1, mean0, mean1, difference,
    SMD_kept, abs_smd,
    # BLP fields:
    coef_pct, p_value, q_value_bh, n_obs, method,
    r2_tau, r2_tau_oos
  )

# rank within treatment: q-value ↑, R2_tau ↓, abs_smd ↑, |difference| ↓
df_binary_topk <- df_binary %>%
  group_by(treatment) %>%
  arrange(
    is.na(q_value_bh), q_value_bh,
    desc(r2_tau),
    abs_smd,
    desc(abs(difference))
  ) %>%
  slice_head(n = top_k) %>%
  ungroup()

# ---------------- outputs ----------------
write.csv(df_numeric_topk, "m2_cate_numeric_top3.csv", row.names = FALSE)
write.csv(df_binary_topk,  "m2_cate_binary_top3.csv",  row.names = FALSE)
cat("Done! Wrote m2_cate_numeric_top5.csv and m2_cate_binary_top5.csv with BLP q-values and R2_tau.\n")
