# Causal



| Folder name  | File name | Description |
| ----------- | ----------- | ----------- |
| 1_data_import |        |        |
|    | 1_data_into_dir.sh        | Creates the folder structure for the running of the analysis, copies over the summary data (filtered down to 0.01 pvalue), gene expression data, individual level data and covariates       |
|  |  2_covariate_information.R   |  Get information on how many samples with/without hypertension, age and gender and Chi-squared test to see if gender is sig. dif between conditions      |
|    | 3_reduce_pval_of_summary_data.sh   |  Create a new summary stat file filtered down to lower pvalues.      |
|    | 4_mapping_SNPs_to_RS.sh  |  Map the SNP to RSS ID using the UKBB GWAS imputed v3   |
|  2_clumping  |        |        |
|    |  1_seperating_by_probe.sh   |  Combine individual level data together and create folder for each probe and create list of SNPs associated with each probe within the folder      |
|    |  2_clumping_data.R   |   Clump the SNPs for each probe using TwoSampleMR (clump_kb = 10000, clump_r2 = 0.1, clump_p1 = 1, clump_p2 = 1,  pop = "EUR")     |
|    |  3_clumping_results.R  | Identify how many SNPs are associated with each probe before and after clumping    |
| 3_preparing_data   |        |        |
|    |   1_Sorting_individual_level.sh     |   get the clumped SNPs and extract individual level data for these     |
|    |    2_Preparing_for_Causal.R    |        |
|    |        |        |
|    |        |        |
|    |        |        |
|    |        |        |
| 4_clustering_data   |        |        |
|    |        |        |
|    |        |        |
| 5_MRPC   |        |        |
|    |        |        |
