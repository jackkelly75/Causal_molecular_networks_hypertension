#module load apps/gcc/R/4.1.0

#eQTL
eqtl <- read.table("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/covariates.txt", header = 1, row.names = 1, sep = "\t")
hypertension <- read.table("/mnt/iusers01/bk01/m20349jk/causal_hypertension/eQTL/pval_01/data/supplied_data/hypertension_phenotype_eQTL.csv", header = 1, row.names = 1, sep = ",")
eqtl <- rbind(eqtl, t(hypertension))


eqtl_control <- eqtl[, eqtl["Hypertension",] == "Y" ]
eqtl_hyper <- eqtl[, eqtl["Hypertension",] == "N" ]





eqtl_control_age <- as.numeric(eqtl_control[1,])
eqtl_hyper_age <- as.numeric(eqtl_hyper[1,])
t.test(eqtl_control_age, eqtl_hyper_age)
#        Welch Two Sample t-test
#data:  eqtl_control_age and eqtl_hyper_age
#t = 8.788, df = 264.04, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  9.046103 14.270241
#sample estimates:
#mean of x mean of y
# 61.70333  50.04516


##
# Chi squared test
##
table(as.numeric(eqtl_control[2,]))
#1 is male
#  0   1
#106 194
table(as.numeric(eqtl_hyper[2,]))
# 0  1
#64 91
#Results
 		hypertension	Control				Row Totals
Female	106  			64 					 170
Male	194  			91 					 285
#Column Totals	300		155					455  (Grand Total)
#The chi-square statistic is 1.5496. The p-value is .213192. The result is not significant at p < .05.



