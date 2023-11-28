
# load packages
library(dplyr)
library(tidyr)
library(psych)

# set working directory
setwd("/cluster/projects/p471/people/Leo/conduct/")

# read Rdata file which was edited in TSD and copied to the cluster
load("PGS_preProcessed.RData")


# write function to regress out covariates and then perform the function
model <- function(y) {
  m <- lm(y ~ processed_prs$PC1 + processed_prs$PC2 + 
            processed_prs$PC3 + processed_prs$PC4 + processed_prs$PC5 + 
            processed_prs$PC6 + processed_prs$PC7 + processed_prs$PC8 + 
            processed_prs$PC9 + processed_prs$PC10 + processed_prs$genotyping_batch + 
            processed_prs$Plate_id + processed_prs$imputation_batch +
              processed_prs$sex, processed_prs$YOB, 
          na.action = na.exclude)
  return(rstandard(m))
}

# check time needed for adjusting all PGS

# try with small subset

processed_prs_orig <- processed_prs
processed_prs <- processed_prs_orig[1:10000, ]

timestamp()

# only feasible on cluster
processed_prsRes <- processed_prs %>% dplyr::mutate_at(vars(tidyr::matches("PGS_")), 
                                                       list(res = ~model(.))) %>% dplyr::select(IID, FID, 
                                                                                                tidyr::matches("PREG_ID"), BARN_NR, tidyr::matches("F_ID"), 
                                                                                                tidyr::matches("M_ID"), Role, tidyr::matches("PGS_"), 
                                                                                                tidyr::matches("PC"), genotyping_batch, Plate_id, imputation_batch)

# write file as output 
save(processed_prsRes, file = "PGS_processed.RData")

