# intergen_conduct

This is a repository for analyses described in "Examining intergenerational risk factors for conduct problems using polygenic scores in the Norwegian Mother, Father and Child Cohort Study" by Frach et al. (2023). Preprint at https://psyarxiv.com/vgkcu/

## Overview
* 01_data_preparation.Rmd - data prep of phenotypic and genetic data from MoBa performed on TSD virtual environment
* 01b_process_PGS_cluster.R - processing of all computed polygenic scores (adjusting for batch effects) performed on the TSD high performance cluster Colossus
* 02_Multiple_imputation_prep_and_diagnostics.R - pre MI: testing associations between a variety of auxilliary variables (registry and questionnaire data) and the phenotype as well as missingness; post MI: diagnostics plots found in the supplementary materials
* MI.R - script used for imputation of phenotype data at age 8 years performed on the TSD HPC Colossus
* 03_Main_analyses_completeData.Rmd - structural equation models and descriptive statistics for analyses using complete data (ages 8 and 14 years) and sensitivity analyses
* model_SEM.R - example script of structural equation models used for imputed data performed for each type polygenic score separately on Colosssus
* 03_Main_analyses_MIdata.Rmd - loading model outputs of all SEMs with imputed data, parameter estimates and further descriptive statistics

The pipeline to compute polygenic scores with LDpred2 in MoBa (on TSD HPC Colossus) can be found at https://github.com/AndreAllegrini/LDpred2 
