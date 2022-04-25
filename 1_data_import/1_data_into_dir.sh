#Making data folders
mkdir causal_hypertension
cd causal_hypertension
mkdir eQTL
mkdir mQTL
mkdir sQTL


#
# eQTL data
#
#individual level data unfiltered is available at 
cd /mnt/bmh01-rds/xiaoguang-cells/eQTL/data/



# copy over the the given data for eqtl
cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL
mkdir pval_01
cd pval_01
mkdir data
cd data
mkdir map_to_rsId
mkdir  supplied_data
cd supplied_data

#copy over the summary data (filtered to 0.01 pvalue)
cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data
cp /mnt/bmh01-rds/xiaoguang-cells/eQTL/output/eqtl_results_p0.01.csv summary_p0.01.csv
#copy over the covariate info
cp /mnt/bmh01-rds/xiaoguang-cells/eQTL/data/covariates.txt ./
cp /mnt/bmh01-rds/xiaoguang-cells/eQTL/hypertension_phenotype_eQTL.csv ./
#copy over the gene expression data
cp /mnt/bmh01-rds/xiaoguang-cells/eQTL/data/gene_expression.txt ./
#copy over the individual level data
cp /mnt/bmh01-rds/xiaoguang-cells/eQTL/data/chr* ./



#set up the files for clumping
cd /mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01
mkdir clumping
cd clumping
mkdir single_probes
mkdir scripts
cd scripts
mkdir clump_0.1
mkdir clump_0.01
                              eQTL
                                |
                              pval_01
                                |
        -----------------------------------
        |                                 |
       data                            clumping 
        |                                 |
   ------------------                --------------               
   |                |                |            |
map_to_rsId  supplied_data       single_probes  scripts
                                                  |          
                                              ------------
                                              |          |
                                         clump_0.1   clump_0.01





#
# copy over the the given data for mqtl
#


#
# copy over the the given data for sqtl
#
