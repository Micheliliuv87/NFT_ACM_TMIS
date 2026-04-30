# NFT_ACM_TMIS


## Paper Assets

More details on paper will be available after the paper is published.

## Analysis Code Files


`CF_M1`:Causal Forest Model for Bored Ape Yacht Club (BAYC) NFTs' first sale price with traits as treatment + time fixed effect control for sale month. 

`PDS_M1`: Contains autocorrelation checks and multicollinearity checks for the covariates used in CF_M1 and run separet PDS model for CF_M1.

`CF_M2_base`: Causal Forest Model for BAYC NFTs' price change under base condition

`CF_M2_samp`: Causal Forest Model for BAYC NFTs' price change under condition A where $\text{Buyer}_{N-1} = \text{Seller}_{N}$

`CF_M2_diffp`: Causal Forest Model for BAYC NFTs' price change under condition B where $\text{Buyer}_{N-1} \neq \text{Seller}_{N}$

`CF_M2_temporal_robustness`: Robustness check for temporal stability of CF_M2 base, A, and B conditions. 

`CF_M3_samp`: Cuasal Forest Model for BAYC NFTs' price change under condition A with joint treatment combining wallet characteristics. 

`CF_M3_diffp`: Causal Forest Model for BAYC NFTs' price change under condition B with joint treatment combining wallet characteristics.

## Data 

Data is collected from both OpenSea and Etherscan APIs and upto 2025-03-07.

Data used for all models under this repository is listed from `df_table1.csv` to`df_table7.csv` and `Panel_for_Model2.csv`. 

Data collection method is displayed in the previous paper code files repository **[NFT_ICIS25](https://github.com/Micheliliuv87/NFT_ICIS25.git)**.

## Results

Results of all models are stored in the `Results` folder. Under each folder, R code files for creating visual tables and further analysis are also included. (please run separately)

## Code File Structure 
```markdown
NFT_ACM_TMIS/
├── CF_M1.R
├── CF_M2_Temporal_Robustness.R
├── CF_M2_base.R
├── CF_M2_diffp.R
├── CF_M2_samp.R
├── CF_M3_Samp.R
├── CF_M3_diffp.R
├── LICENSE
├── PDS_M1.R
├── Panel_for_Model2.csv
├── README.md
├── Results
│   ├── CF_M1_V2
│   │   ├── Output_V7_Final
│   │   ├── V7viz7.R
│   │   ├── cate_reporting_v7.R
│   │   └── nft_trait_v7_results_with_q.numbers
│   ├── CF_M2_base
│   │   ├── Appendix_table
│   │   ├── V11viz.R
│   │   ├── cf_M2_BLP_per_covariate.csv
│   │   ├── cf_M2_Base.png
│   │   ├── cf_M2_R2tau_by_treatment.csv
│   │   ├── cf_M2_balance_diagnostics.csv
│   │   ├── cf_M2_cate_by_covariate_bins.csv
│   │   ├── cf_M2_cate_by_covariate_overview.csv
│   │   ├── cf_M2_covariate_smd_and_cate.csv
│   │   ├── cf_M2_dashboard_GRF.csv
│   │   ├── cf_M2_dashboard_GRF_treatview.csv
│   │   ├── cf_M2_dashboard_PDS.csv
│   │   ├── cf_M2_dashboard_PDS_treatview.csv
│   │   ├── cf_M2_dashboard_compare_GRF_vs_PDS.csv
│   │   ├── cf_M2_forest.pdf
│   │   ├── cf_M2_overlap_diagnostics.csv
│   │   ├── cf_M2_results_CATEs.csv
│   │   ├── cf_M2_results_GRF.csv
│   │   ├── cf_M2_results_PDS.csv
│   │   ├── cf_M2_treatment_results_with_q.csv
│   │   ├── cf_M2_treatment_thresholds.csv
│   │   ├── cf_pds_up_results_combined.csv
│   │   ├── gates_by_treatment
│   │   ├── m2_cate_binary_grouped_table.tex
│   │   ├── m2_cate_binary_top3.csv
│   │   ├── m2_cate_csv.R
│   │   ├── m2_cate_numeric_top3.csv
│   │   ├── m2_cate_numeric_top3_grouped_table.tex
│   │   ├── m2_cate_reporttable_N.R
│   │   ├── m2_cate_repottable_B.R
│   │   └── paper_assets
│   ├── CF_M2_diffp
│   │   ├── V11viz_diffp.R
│   │   ├── cf_M2_diffp.png
│   │   ├── cf_M2_diffp_BLP_per_covariate.csv
│   │   ├── cf_M2_diffp_R2tau_by_treatment.csv
│   │   ├── cf_M2_diffp_balance_diagnostics.csv
│   │   ├── cf_M2_diffp_cate_by_covariate_bins.csv
│   │   ├── cf_M2_diffp_cate_by_covariate_overview.csv
│   │   ├── cf_M2_diffp_covariate_smd_and_cate.csv
│   │   ├── cf_M2_diffp_dashboard_GRF.csv
│   │   ├── cf_M2_diffp_dashboard_GRF_treatview.csv
│   │   ├── cf_M2_diffp_dashboard_PDS.csv
│   │   ├── cf_M2_diffp_dashboard_PDS_treatview.csv
│   │   ├── cf_M2_diffp_dashboard_compare_GRF_vs_PDS.csv
│   │   ├── cf_M2_diffp_forest.pdf
│   │   ├── cf_M2_diffp_overlap_diagnostics.csv
│   │   ├── cf_M2_diffp_results_CATEs.csv
│   │   ├── cf_M2_diffp_results_GRF.csv
│   │   ├── cf_M2_diffp_results_PDS.csv
│   │   ├── cf_M2_diffp_treatment_results_with_q.csv
│   │   ├── cf_M2_diffp_treatment_thresholds.csv
│   │   ├── cf_pds_diffp_results_combined.csv
│   │   ├── gates_by_treatment
│   │   └── paper_assets
│   ├── CF_M2_samp
│   │   ├── V11viz_samp.R
│   │   ├── cf_M2_samp.png
│   │   ├── cf_M2_samp_BLP_per_covariate.csv
│   │   ├── cf_M2_samp_R2tau_by_treatment.csv
│   │   ├── cf_M2_samp_cate_by_covariate_bins.csv
│   │   ├── cf_M2_samp_cate_by_covariate_overview.csv
│   │   ├── cf_M2_samp_dashboard_GRF.csv
│   │   ├── cf_M2_samp_dashboard_GRF_treatview.csv
│   │   ├── cf_M2_samp_dashboard_PDS.csv
│   │   ├── cf_M2_samp_dashboard_PDS_treatview.csv
│   │   ├── cf_M2_samp_dashboard_compare_GRF_vs_PDS.csv
│   │   ├── cf_M2_samp_forest.pdf
│   │   ├── cf_M2_samp_overlap_diagnostics.csv
│   │   ├── cf_M2_samp_results_CATEs.csv
│   │   ├── cf_M2_samp_results_GRF.csv
│   │   ├── cf_M2_samp_results_PDS.csv
│   │   ├── cf_M2_samp_treatment_results_with_q.csv
│   │   ├── cf_M2_samp_treatment_thresholds.csv
│   │   ├── cf_m2_samp_gates_by_treatment
│   │   ├── cf_pds_samp_M2_results_combined.csv
│   │   └── paper_assets
│   ├── CF_M3_diffp
│   │   ├── M3Viz_diffp.R
│   │   ├── cf_M3_BLP_per_covariate.csv
│   │   ├── cf_M3_R2tau_by_treatment.csv
│   │   ├── cf_M3_dashboard_GRF.csv
│   │   ├── cf_M3_dashboard_PDS.csv
│   │   ├── cf_M3_dashboard_compare_GRF_vs_PDS.csv
│   │   ├── cf_M3_diffp_combo.png
│   │   ├── cf_M3_overlap_diagnostics.csv
│   │   ├── cf_M3_results_CATEs.csv
│   │   ├── cf_M3_results_GRF.csv
│   │   ├── cf_M3_results_PDS.csv
│   │   ├── cf_M3_treatment_results_with_q.csv
│   │   ├── cf_M3_treatment_thresholds.csv
│   │   ├── cf_pds_M3_results_combined.csv
│   │   └── gates_by_treatment
│   ├── CF_M3_samp
│   │   ├── M3viz_samp.R
│   │   ├── cf_M3_BLP_per_covariate.csv
│   │   ├── cf_M3_R2tau_by_treatment.csv
│   │   ├── cf_M3_cate_by_covariate_bins.csv
│   │   ├── cf_M3_cate_by_covariate_overview.csv
│   │   ├── cf_M3_dashboard_GRF.csv
│   │   ├── cf_M3_dashboard_GRF_treatview.csv
│   │   ├── cf_M3_dashboard_PDS.csv
│   │   ├── cf_M3_dashboard_PDS_treatview.csv
│   │   ├── cf_M3_dashboard_compare_GRF_vs_PDS.csv
│   │   ├── cf_M3_overlap_diagnostics.csv
│   │   ├── cf_M3_results_CATEs.csv
│   │   ├── cf_M3_results_GRF.csv
│   │   ├── cf_M3_results_PDS.csv
│   │   ├── cf_M3_samp_combo.png
│   │   ├── cf_M3_treatment_results_with_q.csv
│   │   ├── cf_M3_treatment_thresholds.csv
│   │   ├── cf_pds_M3_results_combined.csv
│   │   ├── gates_by_treatment
│   │   └── paper_assets
│   ├── PDS_M1
│   │   ├── PDS_balance_diagnostics.csv
│   │   ├── PDS_dashboard_PDS.csv
│   │   ├── PDS_multicollinearity_condition_number_by_treatment.csv
│   │   ├── PDS_multicollinearity_global_condition_number.csv
│   │   ├── PDS_multicollinearity_global_top_pairs.csv
│   │   ├── PDS_multicollinearity_summary.csv
│   │   ├── PDS_overlap_diagnostics.csv
│   │   ├── PDS_results_PDS.csv
│   │   └── PDS_treatment_thresholds.csv
│   └── Robustness
│       ├── cf_m2_M2V11_temporal_results.csv
│       ├── cf_m2_M2V11_temporal_summary.csv
│       ├── m2_temporal_sample_diagnostics.csv
│       ├── ols_m2_M2V11_month_residual_diag.csv
│       └── pds_m2_M2V11_temporal_results.csv
├── df_table1.csv
├── df_table2.csv
├── df_table3.csv
├── df_table4.csv
├── df_table5.csv
├── df_table6.csv
└── df_table7.csv
```