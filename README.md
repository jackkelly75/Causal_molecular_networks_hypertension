# Causal molecular networks

##### List of files in each directory

</br>

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
|    |    2_Preparing_for_Causal.R   |   remove the probes that don't have any SNPs associated with them after clumping   |
|    | 3_Getting_data_size.R | Get the number of probes for each chromsome file (not being used as not concerened with chromosome information)  |
|    |    4_SVA_prop_var_testing.R | Identify the best number of SV using proportion of variance explained  |
|    |     5.1_remove SVA+batch.R   | Control for SVA and batch effects at the same time using limma |
|    | 5.2_remove SVA no batch.R  |  Control for SVA using limma  |
|    | 5_remove_batch_effects.R  |  Remove the effect of sex, age and transplant status using limma, then calculate SVs and remove using limma   |
|    | 6_map_to_ID.R  |  Map gene to ID, if any duplicate then only keep highest MAD gene |
| 4_clustering_data   |        |        |
|    |  1_hierachical_clustering_TOM.R      | Generate topological overlap matrix for the data  |
|    | 2_hierachical_clustering.R  | WGCNA clustering followed by k-means  |
| 5_MRPC   |        |        |
|    |  1_genotype_data.R    |  For each gene in a module, get the SNPs associated with them created previously and make the R object for analysis    |
|    | 2_MRPC.R  |  Run MRPC with dualPC on each module   |
|    |  3_Identify_important_networks.R | Searches through results of MRPC on each modules and gives vectors with name of modules with causal predictors of hypertension and number of genes in network  |
|    | 4_plotting.R  | plotting the important networks with igraph |
|    | 5_speed_test.R  |  Run MRPC and our method on all modules to get running time |
|    | 6_plot_speed_test.R  | Plot the times for MRPC and our method to visualise the running time |
|    |  MRdualPC.R |  R file containing the MRdualPC functions for running   |
| 6_shiny   |        |        |
|    |  1_shiny_plotting.R  |   Create the data for shiny app and generate the interactive networks in Rstudio     |
|    | server.R  |   server file for deploying shiny app  |
|    | ui.R  |  ui file for deploying shiny app  |


##### Shiny app

https://jack-kelly-manchester.shinyapps.io/Causal_molecular_networks_Hypertension/
