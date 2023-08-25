# Preparation for multiple imputatiin after analyses with complete data


####################################################
# Impute missing Phenotye values for RS_DBD CD scores
####################################################

# clear workspace
rm(list=ls())

# load packages
library(psych)
library(tidyverse)
library(haven)

# Read in processed data from 01_data_prep.R

# load("N:/data/durable/people/Leo/conduct/data/combinedData.RData")
# alldata_old <- alldata

load("N:/durable/people/Leo/conduct/v1/data/combinedData_new.RData")

regPath = "N:/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_MBRN_541_v12.sav"


#---------------------------------------------------
# Prepare auxilliary variables to use for multiple imputation from pheno questionnaires and registry data files
#---------------------------------------------------

# read in relevant variables from registry data

mbrn <- foreign::read.spss(regPath, to.data.frame = TRUE) %>% 
  select(PREG_ID_2306, BARN_NR, KJONN, FAAR, MORS_ALDER, FARS_ALDER, 
         PARITET_5, ROYK_BEG, ROYK_BEG_ANT, VEKT, KMI_SLUTT, 
         SIVST_2, DIABETES_MELLITUS, HYPERTENSJON_ALENE,
         PREEKL, BOHELSEREGION_DAGENS)

mbrn %>% summary()

# high number of missingness for variables DIABETES_MELLITUS, HYPERTENSJON, PREEKL (all introduced in later versions)
# remove these variables (also for maternal BMI but keep this)


mbrn <- mbrn %>% select(-c(DIABETES_MELLITUS, HYPERTENSJON_ALENE, PREEKL))


# other variables which are not available but should be
#HYPERTENSIV_TILSTAND, KOMPLIKASJONER, BOKOMM, BOFYLKE

# data available at each measurement wave
info <- foreign::read.spss("N:/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_status_v12.sav",
                           to.data.frame = T)

info_excl <- foreign::read.spss("N:/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/PDB2306_SV_INFO_v12.sav",
                           to.data.frame = T)


# which variables are in both data sets

names(alldata)[names(alldata) %in% names(mbrn)]


# remove individuals that withdrew consent and remove variables KJONN, FAAR, and age of parents (bc missingness higher than in mbrn)
alldata <- alldata %>% 
  filter(PREG_ID_2306 %in% info_excl$PREG_ID_2306) %>%
  select(-c(KJONN, FAAR, MORS_ALDER, PARITET_5))


# read in ids of trios with pgs to merge with registry data (only keep registry data for those trios with pgs)

mbrn2 <- merge(alldata, mbrn, all.x = T, all.y = F)
info2 <- merge(alldata, info, all.x = T, all.y = F)

summary(mbrn2$VEKT) # birth weight in grams
summary(mbrn2$BOHELSEREGION_DAGENS) # mothers region of residence 

# mom cpd during start of pregnancy - set to 0 if indicated no to smoking during preg
mbrn2$momCPDstPreg <- ifelse(mbrn2$ROYK_BEG == "No", 0, mbrn2$ROYK_BEG_ANT)

hist(mbrn2$momCPDstPreg)
plot(mbrn2$ROYK_BEG)

summary(mbrn2$ROYK_BEG)

table(is.na(mbrn2$momCPDstPreg)) # available for 29461, missing 6665 

mbrn <- mbrn2 %>% select(PREG_ID_2306, BARN_NR, KJONN, FAAR,
                         MORS_ALDER, FARS_ALDER, PARITET_5,
                         momCPDstPreg, VEKT, SIVST_2, BOHELSEREGION_DAGENS, KMI_SLUTT) 

str(mbrn)

mbrn <-  mutate_at(mbrn, c("MORS_ALDER", "FARS_ALDER", "PARITET_5"), function(x) as.numeric(x)) # note, numeric values do not reflect levels, but in correct order so fine for MI

# set missing factor levels of child sex to NA

mbrn$KJONN <- as.factor(ifelse(mbrn$KJONN == "Male", "Male", 
                      ifelse(mbrn$KJONN == "Female", "Female", NA)))

summary(mbrn$KJONN)

mbrn$BOHELSEREGION_DAGENS <- droplevels(mbrn$BOHELSEREGION_DAGENS)

# questionnaire available at each wave
names(info)
table(info2$status_mfr_v12)
table(info2$status_skj1_v12)
table(info2$status_skj3_v12)
table(info2$status_skj4_v12)
table(info2$status_skj5_v12)
table(info2$status_skj6_v12)
table(info2$status_skj5y_v12)
table(info2$status_skj8y_v12)


# individuals participated in questionnaire wave, recode to yes/no
info2 <- info2 %>% 
  select(c(PREG_ID_2306, BARN_NR, status_mfr_v12, status_skj1_v12, status_skj3_v12, status_skj4_v12,
           status_skj5_v12, status_skj6_v12, status_skj5y_v12, status_skj8y_v12)) %>%
  mutate_at(c("status_mfr_v12", "status_skj1_v12", "status_skj3_v12",
                     "status_skj4_v12", "status_skj5_v12", "status_skj6_v12", 
                     "status_skj5y_v12", "status_skj8y_v12"), function(x) ifelse(x == 1, 1, 0))

head(info2)

# data frame including all variables
data <- alldata %>% left_join(mbrn) %>% left_join(info2)

table(data$status_skj8y_v12) # almost all individuals with missing CD values age 8 missed the whole wave/questionnaire


#### test associations of variables with missingness ###
data$missing <- is.na(data$cdSum)
table(data$missing)/length(data$missing) # 53.65% missing pheno data

# pgs (should all be included anyway because they are in the final model)

ttests_PGS_Pval <- data %>% select(contains("_res_")) %>% 
  lapply(., function(x) summary(lm(x ~ data$missing))$coefficients[2, 4]) %>% unlist()

ttests_PGS_R2 <- data %>% select(contains("_res_")) %>% 
  lapply(., function(x) summary(lm(x ~ data$missing))$r.squared) %>% unlist()

sig_pgs_miss <- as.data.frame(ttests_PGS_Pval[which(ttests_PGS_Pval < 0.05)]) # 37 out of 57 PGS significant
goodR2 <- as.data.frame(ttests_PGS_R2[which(ttests_PGS_R2 >= 0.01)]) # none


# auxilliary variables
ttests_aux_Pval <- data %>% select(c(FAAR, MORS_ALDER, FARS_ALDER, PARITET_5,
                                     VEKT, momCPDstPreg, KMI_SLUTT)) %>% 
  lapply(., function(x) summary(lm(x ~ data$missing))$coefficients[2, 4]) %>% unlist()

ttests_aux_R2 <- data %>% select(c(FAAR, MORS_ALDER, FARS_ALDER, PARITET_5,
                               VEKT, momCPDstPreg, KMI_SLUTT)) %>% 
  lapply(., function(x) summary(lm(x ~ data$missing))$r.squared) %>% unlist()

sig_aux_miss <- as.data.frame(ttests_aux_Pval[which(ttests_aux_Pval < 0.05)]) # only FAAR, VEKT, MORS_ALDER, FARS_ALDER, PARITET_5, momCPDstPreg and BMI
goodR2 <- as.data.frame(ttests_aux_R2[which(ttests_aux_R2 >= 0.01)])

table(is.na(data$KMI_SLUTT))/nrow(data) # only 5% have data on maternal BMI

# chi square test for SIVST_2 and BOHELSEREGION_DAGENS and KJONN

chisq.test(data$missing, data$SIVST_2) # significant
effSize <- sqrt(39.177/length(data$SIVST_2)) # very small (0.032)

chisq.test(data$missing, data$BOHELSEREGION_DAGENS) # significant
effSize <- sqrt(6.2934/length(data$SIVST_2)) # very small (0.017)


chisq.test(data$missing, data$KJONN) # not associated with missingness


### association with outcome (CD sum score) ###


# pgs (should all be included anyway because they are in the final model)

lm_PGS_Pval <- data %>% select(contains("_res_")) %>% 
  lapply(., function(x) summary(lm(data$cdSum ~ x))$coefficients[2, 4]) %>% unlist()

lm_PGS_R2 <- data %>% select(contains("_res_")) %>% 
  lapply(., function(x) summary(lm(data$cdSum ~ x))$r.squared) %>% unlist()

sig_pgs_CD <- as.data.frame(lm_PGS_Pval[which(lm_PGS_Pval < 0.05)]) # 27 PGS significant
goodR2 <- as.data.frame(lm_PGS_R2[which(lm_PGS_R2 >= 0.01)]) # none


# auxilliary variables
lm_aux_Pval <- data %>% select(c(FAAR, MORS_ALDER, FARS_ALDER, PARITET_5,
                                     VEKT, momCPDstPreg, KMI_SLUTT)) %>% 
  lapply(., function(x) summary(lm(data$cdSum ~ x))$coefficients[2, 4]) %>% unlist()

lm_aux_R2 <- data %>% select(c(FAAR, MORS_ALDER, FARS_ALDER, PARITET_5,
                                   VEKT, momCPDstPreg, KMI_SLUTT)) %>% 
  lapply(., function(x) summary(lm(data$cdSum ~ x))$r.squared) %>% unlist()

sig_aux_CD <- as.data.frame(lm_aux_Pval[which(lm_aux_Pval < 0.05)]) # only FAAR, MORS_ALDER, momCPDstPreg and VEKT
goodR2 <- as.data.frame(lm_aux_R2[which(lm_aux_R2 >= 0.01)]) # none


### correlations between missingness, outcome (cdSum) and PGS ###


cor_PGS_miss_r <- data %>% select(contains("_res_")) %>% 
  lapply(., function(x) corr.test(x, data$missing)$r) %>% unlist()

goodR <- as.data.frame(cor_PGS_miss_r[which(cor_PGS_miss_r >= 0.1 | cor_PGS_miss_r <= - 0.1)]) # no PGS cor higher than 0.1

cor_PGS_CD_r <- data %>% select(contains("_res_")) %>% 
  lapply(., function(x) corr.test(x, data$cdSum)$r) %>% unlist()

goodR <- as.data.frame(cor_PGS_CD_r[which(cor_PGS_CD_r >= 0.1 | cor_PGS_CD_r <= -0.1)]) # no PGS cor higher than 0.1


# correlations between missingness, outcome (cdSum) and auxilliary variables 

cor_aux_miss_r <- data %>% select(c(FAAR, MORS_ALDER, FARS_ALDER, PARITET_5,
                                    VEKT, momCPDstPreg, KMI_SLUTT)) %>% 
  lapply(., function(x) corr.test(x, data$missing)$r) %>% unlist()

goodR <- as.data.frame(cor_aux_miss_r[which(cor_aux_miss_r >= 0.1 | cor_aux_miss_r <= -0.1)]) # KMI_SLUTT


cor_aux_CD_r <- data %>% select(c(FAAR, MORS_ALDER, FARS_ALDER, PARITET_5,
                                    VEKT, momCPDstPreg, KMI_SLUTT)) %>% 
  lapply(., function(x) corr.test(x, data$cdSum)$r) %>% unlist()

goodR <- as.data.frame(cor_aux_CD_r[which(cor_aux_CD_r >= 0.1 | cor_aux_CD_r <= -0.1)]) # no aux cor higher than 0.1


## CONCLUSION: Use all final PGS, child sex and birth year 




# only auxilliary variables, outcome phenotype (items) and polygenic scores
data_preImp <- data %>% select(PREG_ID_2306, BARN_NR, KJONN, FAAR, KMI_SLUTT,
                               contains("_res_"),
                               NN111:NN118)


### Read data from previous questionnaires ###

# Q1 - Education - For mothers AA1124, fathers AA1126 (code 7 - other education - as missing)
#Household income (Q1, p11) - Mother (AA1315), Father (AA1316)


Q1data <- foreign::read.spss(file='N:/durable/data/MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q1_v12.sav', to.data.frame = TRUE) %>%
  select(PREG_ID_2306, AA1124, AA1126, AA1315, AA1316, AA86, AA87)  

# father income has "do not know" option - code as NA
table(Q1data$AA1316)
Q1data$AA1316 <- ifelse(Q1data$AA1316 == "Do not know", NA, Q1data$AA1316)
Q1data$AA1316 <- Q1data$AA1316 - 1 


library(phenotools)
myphenovars <-  available_variables(source = "moba") %>% 
  variable_search(c("bmi"), where="anywhere") %>% #You can search for specifc variable names, or use words in the variable description to get what you want
  .$var_name

myphenovars <- myphenovars[1]

mydata <- curate_dataset(variables_required= myphenovars,
                         pheno_data_root_dir="N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/",
                         PDB="2306",
                         moba_data_version=12,
                         completion_threshold=0.5,
                         return_items=FALSE,
                         consistent_items=FALSE,
                         out_format="merged_df")

names(mydata)[1] <- "PREG_ID_2306"
mydata <- mydata %>% select(c(PREG_ID_2306, BARN_NR, bmi_derived_m_q1))
  
Q1data <- merge(Q1data, mydata, all.x = T)
Q1data$bmi_derived_m_q1 <- round(Q1data$bmi_derived_m_q1, 4)
summary(Q1data$bmi_derived_m_q1)
table(is.na(Q1data$bmi_derived_m_q1))

## test associations with missingness and CD ##

# correlations between missingness, outcome (cdSum) and auxilliary variables 

data <- merge(data, Q1data, all.x = T)
data_preImp <- merge(data_preImp, Q1data, all.x = T)


# recode to numeric
data <- data %>% mutate_at(vars(AA1124, AA1126, AA1315, AA1316), as.numeric) 

# association with missingness
cor_auxQ1_miss_r <- data %>% select(c(AA1124, AA1126, AA1315, AA1316, bmi_derived_m_q1)) %>% 
  lapply(., function(x) corr.test(x, data$missing)$r) %>% unlist()
cor_auxQ1_miss_p <- data %>% select(c(AA1124, AA1126, AA1315, AA1316, bmi_derived_m_q1)) %>% 
  lapply(., function(x) corr.test(x, data$missing)$p) %>% unlist()

sig_auxQ1_miss <- as.data.frame(cor_auxQ1_miss_p[which(cor_auxQ1_miss_p < 0.05)]) # all significant
goodR <- as.data.frame(cor_auxQ1_miss_r[which(cor_auxQ1_miss_r >= 0.1 | cor_auxQ1_miss_r <= -0.1)])
# only r > |0.1| is AA1124 = maternal education

summary(lm(AA1124 ~ missing, data = data)) # R2 = 0.02


cor_auxQ1_CD_r <- data %>% select(c(AA1124, AA1126, AA1315, AA1316, bmi_derived_m_q1)) %>% 
  lapply(., function(x) corr.test(x, data$cdSum)$r) %>% unlist()
cor_auxQ1_CD_p <- data %>% select(c(AA1124, AA1126, AA1315, AA1316, bmi_derived_m_q1)) %>% 
  lapply(., function(x) corr.test(x, data$cdSum)$p) %>% unlist()

sig_auxQ1_CD <- as.data.frame(cor_auxQ1_CD_p[which(cor_auxQ1_CD_p < 0.05)]) # all significant
goodR <- as.data.frame(cor_auxQ1_CD_r[which(cor_auxQ1_CD_r >= 0.1 | cor_auxQ1_CD_r <= -0.1)]) # no aux cor higher than 0.1



## Q3 ##
# CC1073 testosterone products

testo <- foreign::read.spss(file='N:/durable/data/MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q3_v12.sav', to.data.frame = TRUE) %>%
  select(PREG_ID_2306, CC1073)  

# no mother used testo products during pregnancy (level 4 = 0) and only 3 used it 6 months before preg (and 29 previously)
# 
# DO NOT INCLUDE due to little variation


## Q4 ##
# child mood and temperament, measured by Infant Characteristic Questionnaire (ICQ-6) with the fussy/difficult subscale

icq <-  available_variables(source = "moba") %>% 
  variable_search(c("icq_fussy"), where="anywhere") %>% #You can search for specifc variable names, or use words in the variable description to get what you want
  .$var_name

icq <- curate_dataset(variables_required= icq,
                         pheno_data_root_dir="N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/",
                         PDB="2306",
                         moba_data_version=12,
                         completion_threshold=0.5,
                         return_items=FALSE,
                         consistent_items=FALSE,
                         out_format="merged_df")

names(icq)[1] <- "PREG_ID_2306"
icq <- icq %>% select(c(PREG_ID_2306, BARN_NR, icq_fussy_c_6m))


# check association with missingness and CD sum score

data <- merge(data, icq, all.x = T)
summary(lm(missing ~ icq_fussy_c_6m, data)) # not significant
summary(lm(cdSum ~ icq_fussy_c_6m, data)) # significant but low R2



# Child behaviour checklist - externalizing scale (for ages 18m, 3y and 5y) from phenotools

myphenovars <-  available_variables(source = "moba") %>% 
variable_search(c("CBCL"), where="anywhere") %>% #You can search for specifc variable names, or use words in the variable description to get what you want
  .$var_name

myphenovars <- myphenovars[c(4:6)]

cbcl <- curate_dataset(variables_required= myphenovars,
                      pheno_data_root_dir="N:/data/durable/data/MoBaPhenoData/PDB2306_MoBa_V12/SPSS/",
                      PDB="2306",
                      moba_data_version=12,
                      completion_threshold=0.5,
                      return_items=FALSE,
                      consistent_items=FALSE,
                      out_format="merged_df")


names(cbcl)[1] <- "PREG_ID_2306"
cbcl <- cbcl %>% select(c(PREG_ID_2306, BARN_NR, cbcl_ext_c_18m, cbcl_ext_c_3yr, cbcl_ext_c_5yr))

data <- merge(data, cbcl, all.x = T)


## associations with missingness and CD sum

summary(lm(cbcl_ext_c_18m ~ missing, data)) # significant but low R2
summary(lm(cbcl_ext_c_3yr ~ missing, data)) # significant but low R2
summary(lm(cbcl_ext_c_5yr ~ missing, data)) # significant but low R2

summary(lm(cdSum ~ cbcl_ext_c_18m, data)) # significant and R2 = 0.029
summary(lm(cdSum ~ cbcl_ext_c_3yr, data)) # significant and R2 = 0.063
summary(lm(cdSum ~ cbcl_ext_c_5yr, data)) # significant and R2 = 0.126



## Q5 ##
# behavioural problems (EE847-EE850), living with child's father (EE488), and coping with financial situations (EE584)

Q5data <- foreign::read.spss(file='N:/durable/data/MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q5_18months_v12.sav', to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR, EE488, EE584, EE847:EE850)  

summary(Q5data) # no variation at all in behavioural problems 
table(Q5data$EE850) # only 122 answered question 4 (referral to professional), 93 no, 29 yes

# DO NOT USE behavioural problems

Q5data <- Q5data %>% select(PREG_ID_2306, BARN_NR, EE488, EE584)
str(Q5data)

table(Q5data$EE488)
table(Q5data$EE584)

length(unique(Q5data$PREG_ID_2306))

# duplicated preg IDs
#Q5data <- Q5data %>% arrange(PREG_ID_2306, trio_pheno) %>% distinct(PREG_ID_2306, .keep_all = T)

check <- Q5data[duplicated(Q5data$PREG_ID_2306), ]

Q5data$EE488 <- as.factor(ifelse(Q5data$EE488 == "Yes", "Yes",
                       ifelse(Q5data$EE488 == "No", "No", NA)))

Q5data$EE584 <- as.ordered(dplyr::recode_factor(Q5data$EE584, 'No, never' = 1,
                                     'Yes, sometimes' = 2,
                                     'Yes, sometimes_duplicated_3' = 3,
                                     'Yes, often' = 4
                                     ))

## test associations with missingness and Cd sum

data <- left_join(data, Q5data)

chisq.test(data$EE488, data$missing)
effSize <- sqrt(19.612/length(data$EE488)) # low effect size

chisq.test(data$EE584, data$missing) # significant
effSize <- sqrt(72.084/length(data$EE584)) # small effect size


summary(lm(cdSum ~ EE488, data)) # significant but low R2
summary(lm(cdSum ~ EE584, data)) # significant but low R2



## Q6 ##
# living with child's father GG383 
# not needed as previous time point it was significant but low R2?



## Q 5 years ##
# LL94 behavioural problems

Q5ydata <- foreign::read.spss(file='N:/durable/data/MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q5yrs_v12.sav', to.data.frame = TRUE) %>%
  select(PREG_ID_2306, LL94, BARN_NR)  

table(Q5ydata$LL94) # 851 yes (ever had)
Q5ydata$LL94 <- as.factor(ifelse(Q5ydata$LL94 == "Yes", "Yes",
                       ifelse(Q5ydata$LL94 == "No", "No", NA)))

data <- merge(data, Q5ydata, all.x = T)

## associations with missingness and CD sum score

chisq.test(data$LL94, data$missing) # not significant

summary(lm(cdSum ~ LL94, data)) # significant and R2 = 0.045



## Q 8 years ##
# NN56-NN59 - Behavioural problems
# NN271 - maternal education

Q8ydata <- foreign::read.spss(file='N:/durable/data/MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q8yrs_v12.sav', to.data.frame = TRUE) %>%
  select(PREG_ID_2306, BARN_NR, NN271, NN56:NN59)  

table(Q8ydata$NN56) # no (never had)
table(Q8ydata$NN57) # yes currently
table(Q8ydata$NN58) # yes in the past (already covered by Q5y?)
table(Q8ydata$NN59) # if yes, referral to specialist? avoid this as conditional

Q8ydata$behav_prob[Q8ydata$NN57 == 1 | Q8ydata$NN58 == 1] <- 1
Q8ydata$behav_prob[Q8ydata$NN56 == 1] <- 0
Q8ydata$behav_prob[Q8ydata$behav_prob == ""] <- NA

table(Q8ydata$behav_prob)
summary(Q8ydata$behav_prob)

Q8ydata <- Q8ydata %>% select(PREG_ID_2306, BARN_NR, NN271, behav_prob)

table(Q8ydata$behav_prob)
Q8ydata$NN271 <- as.ordered(Q8ydata$NN271)

## associations with missingness and CD sum

data <- left_join(data, Q8ydata)

summary(lm(behav_prob ~ missing, data)) # not significant
chisq.test(data$behav_prob, data$missing) # neither

summary(lm(cdSum ~ behav_prob, data)) # significant, R2 = 0.16

chisq.test(data$NN271, data$missing) # not significant
summary(lm(cdSum ~ NN271, data)) # significant, very low R2



#---------------------------------------------------
# Additional variables using PhenoTools
#---------------------------------------------------

library(phenotools)
vignette('phenotools')

# Check available ICD codes that correspond with diagnosis of hyperkinetic disorder, conduct disorder, or mixed disorders of conduct and emotions

  
diagn <-  available_variables(source = "npr") %>% 
  variable_search(c("F90", "F91", "F92"), where="anywhere") %>% #You can search for specifc variable names, or use words in the variable description to get what you want
  .$var_name

# read in child ADHD and disruptive behaviour disorders diagnosess:
diagn <- diagn[c(1, 12)]



diagn_child <- curate_npr(diagnoses = diagn,
                          out_format = "merged_df", 
                          npr_filename = "18_34161_NPR.sav")


table(diagn_child$received_dx_F90_npr)

table(diagn_child$received_dx_F90_npr, useNA = "ifany") #4590  of 114114 (4%)

table(diagn_child$received_dx_F91_npr, useNA = "ifany") #596 of 114114 (0.52%)
#table(diagn_child$received_dx_F92_npr, useNA = "ifany") #391 of 114114 (0.45%)

diagn_child_F90 <- diagn_child %>% select(contains("F90") & contains("received"))

# diagn_child$received_dx_F91orF92 <- ifelse(diagn_child$received_dx_F91_npr == "yes" | 
#                                              diagn_child$received_dx_F92_npr == "yes", "yes", 
#                                            ifelse(is.na(diagn_child$received_dx_F91_npr) & 
#                                                     is.na(diagn_child$received_dx_F92_npr), NA, "no"))



# table(diagn_child$received_dx_F91orF92, useNA = "ifany") #945 of 114114 (0.83%)

NPRdata_child <- diagn_child %>% select(preg_id, BARN_NR, received_dx_F90_npr, received_dx_F91_npr) %>%
  rename(PREG_ID_2306 = preg_id, ADHDdiag_child = received_dx_F90_npr, CDdiag_child = received_dx_F91_npr)

# Set preg id and birth order as integer to join with other dataframes:
NPRdata_child$PREG_ID_2306 <- as.integer(NPRdata_child$PREG_ID_2306)
NPRdata_child$BARN_NR <- as.integer(NPRdata_child$BARN_NR)

# Set missing values to zero for MI and convert to factor for MI:

# NPRdata_child$ADHDdiag_child[is.na(NPRdata_child$ADHDdiag_child)] = 0
# NPRdata_child$ADHDdiag_child <- ifelse(NPRdata_child$ADHDdiag_child == "yes", "Yes", "No")
# 
# NPRdata_child$CDdiag_child[is.na(NPRdata_child$CDdiag_child)] = 0
# NPRdata_child$CDdiag_child <- ifelse(NPRdata_child$CDdiag_child == "Yes", "Yes", "No")

NPRdata_child$ADHDdiag_child <- as.factor(NPRdata_child$ADHDdiag_child)
NPRdata_child$CDdiag_child <- as.factor(NPRdata_child$CDdiag_child)

table(NPRdata_child$ADHDdiag_child) # 6001    of 114013
table(NPRdata_child$CDdiag_child) #691  of 114013


#-------------- Parental ADHD and CD diagnoses ----------------#

## Mothers:
diagn_parents <-  available_variables(source = "npr") %>% 
  variable_search(c("F90", "F91", "F10", 
                    "F32", "F33", "F34", "F602"), 
                  where="anywhere") %>% #You can search for specifc variable names, or use words in the variable description to get what you want
  .$var_name


diagn_mother <- curate_npr(diagnoses = diagn_parents,
                           return_items = F, 
                           dx_owner = "mother",
                           npr_filename = "18_34161_NPR.sav")

head(diagn_mother)

table(diagn_mother$received_dx_F90_npr, useNA = "ifany") #1431of 114114 (1.25%)
table(diagn_mother$received_dx_F91_npr, useNA = "ifany") #16 of 114114 (0.01%)

table(diagn_mother$received_dx_F32_npr, useNA = "ifany") #6782  of 114114 (5.94%)
table(diagn_mother$received_dx_F33_npr, useNA = "ifany") #4404  of 114114 (3.86%)
table(diagn_mother$received_dx_F34_npr, useNA = "ifany") #647  of 114114 (0.57%)

table(diagn_mother$received_dx_F602_npr, useNA = "ifany") #10 of 114013

diagn_mother_F90 <- diagn_mother %>% select(contains("F90") & contains("received"))


table(diagn_mother$received_dx_F10_npr, useNA = "ifany") #928        of 114013 (0.67%)
# table(diagn_mother$received_dx_F11_npr, useNA = "ifany") #183   of 114013 (0.16%)
# table(diagn_mother$received_dx_F12_npr, useNA = "ifany") #141  of 114013 (0.12%)
# table(diagn_mother$received_dx_F13_npr, useNA = "ifany") #274  of 114013 (0.24%)
# table(diagn_mother$received_dx_F14_npr, useNA = "ifany") #22  of 114013 (0.01%)
# table(diagn_mother$received_dx_F15_npr, useNA = "ifany") #166  of 114013 (0.15%)
# table(diagn_mother$received_dx_F16_npr, useNA = "ifany") #10  of 114013 
# table(diagn_mother$received_dx_F17_npr, useNA = "ifany") #90  of 114013 (0.08%)
# table(diagn_mother$received_dx_F18_npr, useNA = "ifany") #2  of 114013 
# table(diagn_mother$received_dx_F19_npr, useNA = "ifany") #285  of 114013 (0.25%)


# diagn_mother$received_dx_F10toF19 <- ifelse(diagn_mother$received_dx_F10_npr == "yes" | 
#                                                   diagn_mother$received_dx_F11_npr == "yes" |
#                                                   diagn_mother$received_dx_F12_npr == "yes" |
#                                                   diagn_mother$received_dx_F13_npr == "yes" |
#                                                   diagn_mother$received_dx_F14_npr == "yes" |
#                                                   diagn_mother$received_dx_F15_npr == "yes" |
#                                                   diagn_mother$received_dx_F16_npr == "yes" |
#                                                   diagn_mother$received_dx_F17_npr == "yes" |
#                                                   diagn_mother$received_dx_F18_npr == "yes" |
#                                                   diagn_mother$received_dx_F19_npr == "yes",
#                                                 "Yes", "No")
# 
# table(diagn_mother$received_dx_F10toF19) # 1261 of 114233 (1.1%)

diagn_mother$received_dx_F32toF34 <- ifelse(diagn_mother$received_dx_F32_npr == "yes" | 
                                              diagn_mother$received_dx_F33_npr == "yes" |
                                              diagn_mother$received_dx_F34_npr == "yes",
                                                "Yes", "No")

table(diagn_mother$received_dx_F32toF34) # 11395 of 114013 (9.99%)

NPRdata_mother <- diagn_mother %>% 
  select(preg_id, BARN_NR, received_dx_F90_npr, received_dx_F32toF34, received_dx_F10_npr) %>%
  rename(PREG_ID_2306 = preg_id, ADHDdiag_mother = received_dx_F90_npr, 
         DEPRdiag_mother = received_dx_F32toF34, SUBSTdiag_mother = received_dx_F10_npr)


# Set preg id and birth order as integer to join with other dataframes:
NPRdata_mother$PREG_ID_2306 <- as.integer(NPRdata_mother$PREG_ID_2306)
NPRdata_mother$BARN_NR <- as.integer(NPRdata_mother$BARN_NR)

# Set missing values to zero for MI and convert to factor for MI:

# NPRdata_mother$ADHDdiag_mother[is.na(NPRdata_mother$ADHDdiag_mother)] = 0
# NPRdata_mother$ADHDdiag_mother <- ifelse(NPRdata_mother$ADHDdiag_mother == "yes", "Yes", "No")
# 
# NPRdata_mother$DEPRdiag_mother[is.na(NPRdata_mother$DEPRdiag_mother)] = 0
# NPRdata_mother$DEPRdiag_mother <- ifelse(NPRdata_mother$DEPRdiag_mother == "Yes", "Yes", "No")
# 
# NPRdata_mother$SUBSTdiag_mother[is.na(NPRdata_mother$SUBSTdiag_mother)] = 0
# NPRdata_mother$SUBSTdiag_mother <- ifelse(NPRdata_mother$SUBSTdiag_mother == "Yes", "Yes", "No")


NPRdata_mother$ADHDdiag_mother <- as.factor(NPRdata_mother$ADHDdiag_mother)
NPRdata_mother$SUBSTdiag_mother <- as.factor(NPRdata_mother$SUBSTdiag_mother)
NPRdata_mother$DEPRdiag_mother <- as.factor(NPRdata_mother$DEPRdiag_mother)



## Fathers:

rm(info2, mbrn2, mbrn, Q1data, Q5data, Q5ydata, Q8ydata, diagn_child, diagn_child_F90, testo, cbcl)
rm(alldata, check, icq, info, info_excl)
rm(diagn_mother, diagn_mother_F90, mydata)
data_preImp <- data_preImp[1:10, ]

gc()

diagn_parents1 <- diagn_parents[1:52]
diagn_parents2 <- diagn_parents[53:88]

diagn_father <- curate_npr(diagnoses = diagn_parents1,
                           return_items = F, 
                           dx_owner = "father",
                           npr_filename = "18_34161_NPR.sav")

gc()

diagn_father2 <- curate_npr(diagnoses = diagn_parents2,
                           return_items = F, 
                           dx_owner = "father",
                           npr_filename = "18_34161_NPR.sav")
gc()

diagn_father1 <- diagn_father

diagn_father <- dplyr::left_join(diagn_father1, diagn_father2)


#diagn_father_npr <- diagn_father$npr

table(diagn_father$received_dx_F90_npr, useNA = "ifany") #1146 of 87882 (1.03%)
table(diagn_father$received_dx_F91_npr, useNA = "ifany") #8 of 87882 

table(diagn_father$received_dx_F32_npr, useNA = "ifany") #3484   of 87882 (3.31%)
table(diagn_father$received_dx_F33_npr, useNA = "ifany") #2107   of 87882 (2.01%)
table(diagn_father$received_dx_F34_npr, useNA = "ifany") #396  of 87882 (0.40%)

table(diagn_father$received_dx_F602_npr, useNA = "ifany") #51 of 87882 (0.05%)


table(diagn_father$received_dx_F10_npr, useNA = "ifany") #1373   of 87882 (1.28%)
# table(diagn_father_npr$received_dx_F11_npr, useNA = "ifany") #231   of 87979 (0.26%)
# table(diagn_father_npr$received_dx_F12_npr, useNA = "ifany") #275  of 87979 (0.31%)
# table(diagn_father_npr$received_dx_F13_npr, useNA = "ifany") #261  of 87979 (0.30%)
# table(diagn_father_npr$received_dx_F14_npr, useNA = "ifany") #48  of 87979 (0.05%)
# table(diagn_father_npr$received_dx_F15_npr, useNA = "ifany") #243  of 87979 (0.28%)
# table(diagn_father_npr$received_dx_F16_npr, useNA = "ifany") #18  of 87979 
# table(diagn_father_npr$received_dx_F17_npr, useNA = "ifany") #61  of 87979 (0.07%)
# table(diagn_father_npr$received_dx_F18_npr, useNA = "ifany") #4  of 87979 
# table(diagn_father_npr$received_dx_F19_npr, useNA = "ifany") #345  of 87979 (0.39%)
# 
# 
# diagn_father_npr$received_dx_F10toF19 <- ifelse(diagn_father_npr$received_dx_F10_npr == "yes" | 
#                                                   diagn_father_npr$received_dx_F11_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F12_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F13_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F14_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F15_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F16_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F17_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F18_npr == "yes" |
#                                                   diagn_father_npr$received_dx_F19_npr == "yes",
#                                                 "Yes", "No")
# 
# table(diagn_father_npr$received_dx_F10toF19) # 1636 of 87979 (1.86%)

diagn_father$received_dx_F32toF34 <- ifelse(diagn_father$received_dx_F32_npr == "yes" | 
                                              diagn_father$received_dx_F33_npr == "yes" |
                                              diagn_father$received_dx_F34_npr == "yes",
                                            "yes", ifelse(is.na(diagn_father$received_dx_F32_npr) & 
                                                            is.na(diagn_father$received_dx_F33_npr) &
                                                            is.na(diagn_father$received_dx_F34_npr), NA,
                                                          "No"))


table(diagn_father$received_dx_F32toF34) # 4882  of 87882 (4.72%)

library(tidyverse)
NPRdata_father <- diagn_father %>% 
  select(preg_id, BARN_NR, received_dx_F90_npr, received_dx_F32toF34, received_dx_F10_npr) %>%
  rename(PREG_ID_2306 = preg_id, ADHDdiag_father = received_dx_F90_npr, 
         DEPRdiag_father = received_dx_F32toF34, SUBSTdiag_father = received_dx_F10_npr)


# Set preg id and birth order as integer to join with other dataframes:
NPRdata_father$PREG_ID_2306 <- as.integer(NPRdata_father$PREG_ID_2306)
NPRdata_father$BARN_NR <- as.integer(NPRdata_father$BARN_NR)

# Set missing values to zero for MI and convert to factor for MI:

# NPRdata_father$ADHDdiag_father[is.na(NPRdata_father$ADHDdiag_father)] = 0
# NPRdata_father$ADHDdiag_father <- ifelse(NPRdata_father$ADHDdiag_father == "yes", "Yes", "No")
# 
# NPRdata_father$DEPRdiag_father[is.na(NPRdata_father$DEPRdiag_father)] = 0
# NPRdata_father$DEPRdiag_father <- ifelse(NPRdata_father$DEPRdiag_father == "Yes", "Yes", "No")
# 
# NPRdata_father$SUBSTdiag_father[is.na(NPRdata_father$SUBSTdiag_father)] = 0
# NPRdata_father$SUBSTdiag_father <- ifelse(NPRdata_father$SUBSTdiag_father == "Yes", "Yes", "No")


NPRdata_father$ADHDdiag_father <- as.factor(NPRdata_father$ADHDdiag_father)
NPRdata_father$SUBSTdiag_father <- as.factor(NPRdata_father$SUBSTdiag_father)
NPRdata_father$DEPRdiag_father <- as.factor(NPRdata_father$DEPRdiag_father)


#----------------------------------------------------------#


# merge NPR data
NPRdata <- NPRdata_child %>% left_join(NPRdata_father) %>% left_join(NPRdata_mother)


#---------------------------------------------------
# Merge precomputed pgs and registry data files
#---------------------------------------------------


data <- data %>% left_join(NPRdata)


# test associations between NPRdata and missingness and cd sum score
# 
# child
summary(lm(missing ~ ADHDdiag_child, data)) # significant but low R2
summary(lm(cdSum ~ ADHDdiag_child, data)) # significant, R2 = 5.8%

summary(lm(missing ~ CDdiag_child, data)) # significant but low R2
summary(lm(cdSum ~ CDdiag_child, data)) # significant, R2 = 3.5%


# mother
summary(lm(missing ~ ADHDdiag_mother, data)) # significant, R2 = 0.1%
summary(lm(cdSum ~ ADHDdiag_mother, data)) # significant, R2 = 0.03%

summary(lm(missing ~ DEPRdiag_mother, data)) # significant, R2 = 0.2%
summary(lm(cdSum ~ DEPRdiag_mother, data)) # significant, R2 = 0.04%

summary(lm(missing ~ SUBSTdiag_mother, data)) # significant, low R2 
summary(lm(cdSum ~ SUBSTdiag_mother, data)) # significant, low R2 


# father
summary(lm(missing ~ ADHDdiag_father, data)) # significant, low R2 
summary(lm(cdSum ~ ADHDdiag_father, data)) # significant, R2 = 0.2%

summary(lm(missing ~ DEPRdiag_father, data)) # significant, low R2 
summary(lm(cdSum ~ DEPRdiag_father, data)) # significant, R2 = 0.1%

summary(lm(missing ~ SUBSTdiag_father, data)) # significant, low R2 
summary(lm(cdSum ~ SUBSTdiag_father, data)) # significant, R2 = 0.1%


## only keep variables used for MI ##

MIdata <- data %>% select(c(names(data_preImp), 
                            ADHDdiag_child, CDdiag_child, cbcl_ext_c_18m, cbcl_ext_c_3yr, cbcl_ext_c_5yr, LL94, NN271, behav_prob))

MIdata <- MIdata %>% select(-c(bmi_derived_m_q1, AA1315, AA1316, AA1126, AA86, AA87,
                               NN271, KMI_SLUTT))


#recode items as factors
MIdata <- MIdata %>% mutate(across(NN111:NN118, as.ordered))
class(MIdata$NN111)
MIdata$behav_prob <- as.factor(MIdata$behav_prob)

head(MIdata)
str(MIdata)

table(MIdata$ADHDdiag_child, useNA = "ifany") #1533  of 31337 (3.7%)

summary(lm(missing ~ BARN_NR, data = data))
summary(lm(cdSum ~ BARN_NR, data = data))

MIdata$ADHDdiag_child <- as.factor(MIdata$ADHDdiag_child)
MIdata$CDdiag_child <- as.factor(MIdata$CDdiag_child)

# drop levels of binary factors for impuation
#MIdata <- MIdata %>% mutate(across(c(KJONN, LL94, behav_prob, ADHDdiag_child, CDdiag_child), as.factor))

# save file
save(MIdata, file = "N:/durable/people/Leo/conduct/v1/data/MIdata.Rdata")


# load MIdata 
load("N:/durable/people/Leo/conduct/v1/data/MIdata.Rdata")


# # install packages
# install.packages("withr")
# install.packages("mice")
# install.packages("Rcpp")

# convert all factors to numeric
cordata <- dplyr::mutate_all(MIdata, function(x) as.numeric(x))
str(cordata) # for cor to work
corMat <- cor(cordata[, 44:59], use = "pairwise.complete.obs")
corMat <- ifelse(corMat == 1, NA, corMat)
max(corMat, na.rm = T)

# compute correlations between rsdbd8 and all other vars
# i1 <- sapply(cordata, is.numeric)
# y1 <- "cdSum" 
# x1 <- setdiff(names(cordata)[i1], y1)
# (corRes <- cor(cordata[x1], cordata[[y1]],use="pairwise.complete.obs"))


library(mice)

#percentage of missing data
p <- function(x) {sum(is.na(x))/length(x)*100}
apply(MIdata,2,p)
(missDatPerc <- md.pairs(MIdata))

test <- MIdata[1:50, ]

time_pre <- Sys.time()
MItest <- mice(test[,3:dim(test)[2]], m = 100, seed = 123) #note, exlude preg ID and BARN_NR here
time_post <- Sys.time()

## the imputation for the whole data set was done on the cluster using the script MI.R and process_MI.bash in scripts/PGS_bash 

# MIfull <- mice(MIdata[, 3:dim(test)[2]], m = 100, maxit = 30, seed = 123) #note, exlude preg ID and BARN_NR here
# 
# # create list object with element for each of 100 imputed datasets, then add preg id and birth order to each:
# 
# mice.imp <- NULL
# 
# for(i in 1:100) {
#   mice.imp[[i]] <- complete(MIfull, action=i, include=FALSE)
#   mice.imp[[i]]$PREG_ID_2306 <- MIdata$PREG_ID_2306
#   mice.imp[[i]]$BARN_NR <- MIdata$BARN_NR
# }
# 
# head(mice.imp[[10]])



## load test data to compare results from different imputations
load("N:/durable/people/Leo/conduct/v1/data/MI_test/Test_MIfull.RData")
load("N:/durable/people/Leo/conduct/v1/data/MI_test/Test_MIfull_num.RData")
load("N:/durable/people/Leo/conduct/v1/data/MIdata.Rdata") # non-imputed data

library(tidyverse)

# check correlations between items for non-imputed data
cor_raw2 <- MIdata %>% select(NN111:NN118) %>% 
  mutate_all(., as.numeric) %>% cor(., use = "complete.obs")

round(cor_raw2, digits = 2)

# correlation of ordinal values
polycorr <- MIdata %>% select(NN111:NN118) %>% psych::polychoric()
round(polycorr$rho, digits = 2)

# resulting estimates are much higher


## check correlations between items for both types of imputation (test data set, not complete)

# to perform correlation on pooled data, fisher transformation before and after pooling is recommended
fisher.trans <- function(x) 1/2 * log((1 + x) / (1 - x))
fisher.backtrans <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)

names <- MIdata %>% select(NN111:NN118) %>% names()

# recode items
# 


polycorr_mi <- MIfull %>%
  mice::complete("all") %>%
  map(select, NN111:NN118) %>%
  map(psych::polychoric)

# check only for first 5
round(polycorr_mi$`1`$rho, digits = 2)
round(polycorr_mi$`2`$rho, digits = 2)
round(polycorr_mi$`3`$rho, digits = 2)
round(polycorr_mi$`4`$rho, digits = 2)
round(polycorr_mi$`5`$rho, digits = 2)

# all very similar to unimputed values

MIfull$data$NN111 <- MIfull$data$NN111 %>% as.numeric()
MIfull$data$NN112 <- MIfull$data$NN112 %>% as.numeric()
MIfull$data$NN113 <- MIfull$data$NN113 %>% as.numeric()
MIfull$data$NN114 <- MIfull$data$NN114 %>% as.numeric()
MIfull$data$NN115 <- MIfull$data$NN115 %>% as.numeric()
MIfull$data$NN116 <- MIfull$data$NN116 %>% as.numeric()
MIfull$data$NN117 <- MIfull$data$NN117 %>% as.numeric()
MIfull$data$NN118 <- MIfull$data$NN118 %>% as.numeric()

MIfull$imp$NN111 <- MIfull$imp$NN111 %>% mutate_all(as.numeric)
MIfull$imp$NN112 <- MIfull$imp$NN112 %>% mutate_all(as.numeric)
MIfull$imp$NN113 <- MIfull$imp$NN113 %>% mutate_all(as.numeric)
MIfull$imp$NN114 <- MIfull$imp$NN114 %>% mutate_all(as.numeric)
MIfull$imp$NN115 <- MIfull$imp$NN115 %>% mutate_all(as.numeric)
MIfull$imp$NN116 <- MIfull$imp$NN116 %>% mutate_all(as.numeric)
MIfull$imp$NN117 <- MIfull$imp$NN117 %>% mutate_all(as.numeric)
MIfull$imp$NN118 <- MIfull$imp$NN118 %>% mutate_all(as.numeric)


cor <- MIfull %>%
  mice::complete("all") %>%
  map(select, NN111:NN118) %>%
  map(stats::cor) %>%
  map(fisher.trans)

cor.rect <- Reduce("+", cor) / length(cor) # m is equal to the length of the list
cor.rect <- fisher.backtrans(cor.rect)

round(cor.rect, digits = 2)

# correlations using ordered values are unrealistically high!


# correlate results of numeric and ordered values

# cor_num <- MIfull_num %>%
#   mice::complete("all") %>%
#   map(select, NN111:NN118) %>%
#   map(stats::cor) %>%
#   map(fisher.trans)
# 
# cor.rect_num <- Reduce("+", cor_num) / length(cor_num) # m is equal to the length of the list
# cor.rect_num <- fisher.backtrans(cor.rect_num)
# 
# round(cor.rect_num, digits = 2)
# round(cor_raw2, digits = 2)



## load the data resulting from imputation
# 

load("N:/durable/people/Leo/conduct/v1/data/MIfull.RData") # results from MI
load("N:/durable/people/Leo/conduct/v1/data/MIfull_with_IDs.RData") # results from MI merged with

#save(mice.imp, file = "N:/durable/people/Leo/conduct/v1/data/MIfull_with_IDs.RData")

# repeat correlation for the complete data set (imputed)
# 
cor <- MIfull %>%
mice::complete("all") %>%
  map(select, NN111:NN118) %>%
  map(stats::cor) %>%
  map(fisher.trans)

cor.rect <- Reduce("+", cor) / length(cor) # m is equal to the length of the list
cor.rect <- fisher.backtrans(cor.rect)

(cor.rect <- round(cor.rect, digits = 2))

cor_raw <- MIfull$data %>% select(NN111:NN118) %>% 
   cor(., use = "complete.obs")

(cor_raw <- round(cor_raw, digits = 2))


apaTables::apa.cor.table(MIfull$data[names], filename = "N:/durable/people/Leo/conduct/v1/output/cor_Items_complete.docx")

write.table(cor_raw, file = "N:/durable/people/Leo/conduct/v1/output/cor_Items_raw.txt")
write.table(cor.rect, file = "N:/durable/people/Leo/conduct/v1/output/cor_Items_imputed.txt")

# MIfull %>%
#   mice::complete("all") %>%
#   map(select, NN111:NN118) %>%
#   map(apaTables::apa.cor.table, filename = "cor_Items_imputed.docx") %>% mice::pool()

  

# plot the trace lines against iteration number to check convergence
# modify plot function provided by package mice


# used from stackoverflow: https://stackoverflow.com/questions/34541682/xyplot-bottom-axis-when-last-row-has-fewer-panels-than-columns

lattice::trellis.par.set(clip = list(panel = "off"))
myPanel <- function(...){
  lattice::panel.xyplot(...)
  if(lattice::panel.number() == 28) {
    at = seq(0, 30, by = 5)
    lattice::panel.axis("bottom", at = at, outside = T,
                        labels = T, half = F)
  }
  if(lattice::panel.number() == 27) {
    at = seq(0, 30, by = 5)
    lattice::panel.axis("bottom", at = at, outside = T,
                        labels = T, half = F)
  }
  
  lattice::panel.xyplot(...)
  if(lattice::panel.number() == 26) {
    at = seq(0, 30, by = 5)
    lattice::panel.axis("bottom", at = at, outside = T,
                        labels = T, half = F)
  }
  if(lattice::panel.number() == 25) {
    at = seq(0, 30, by = 5)
    lattice::panel.axis("bottom", at = at, outside = T,
                        labels = T, half = F)
  }
}


#' Plot the trace lines of the MICE algorithm
#'
#' Trace line plots portray the value of an estimate
#' against the iteration number. The estimate can be anything that you can calculate, but
#' typically are chosen as parameter of scientific interest. The \code{plot} method for
#' a \code{mids} object plots the mean and standard deviation of the imputed (not observed)
#' values against the iteration number for each of the $m$ replications. By default,
#' the function plot the development of the mean and standard deviation for each incomplete
#' variable. On convergence, the streams should intermingle and be free of any trend.
#'
#' @param x      An object of class \code{mids}
#' @param y      A formula that specifies which variables, stream and iterations are plotted.
#'               If omitted, all streams, variables and iterations are plotted.
#' @param theme  The trellis theme to applied to the graphs. The default is \code{mice.theme()}.
#' @param layout A vector of length 2 given the number of columns and rows in the plot.
#'               The default is \code{c(2, 3)}.
#' @param type   Parameter \code{type} of \code{\link{panel.xyplot}}.
#' @param col    Parameter \code{col} of \code{\link{panel.xyplot}}.
#' @param lty    Parameter \code{lty} of \code{\link{panel.xyplot}}.
#' @param ...    Extra arguments for \code{\link{xyplot}}.
#' @return An object of class \code{"trellis"}.
#' @author Stef van Buuren 2011
#' @seealso \code{\link{mice}}, \code{\link[=mids-class]{mids}},
#' \code{\link{xyplot}}
#' @method plot mids
#' @examples
#' imp <- mice(nhanes, print = FALSE)
#' plot(imp, bmi + chl ~ .it | .ms, layout = c(2, 1))
#' @export

plot.mids <- function(x, y = NULL, theme = mice.theme(), layout = c(4, 9),
                      type = "l", col = 1:10, lty = 1, ...) {
  strip.combined <- function(which.given, which.panel, factor.levels, ...) {
    if (which.given == 1) {
      lattice::panel.rect(0, 0, 1, 1,
                          col = theme$strip.background$col, border = 1
      )
      lattice::panel.text(
        x = 0, y = 0.5, pos = 4,
        lab = factor.levels[which.panel[which.given]]
      )
    }
    if (which.given == 2) {
      lattice::panel.text(
        x = 1, y = 0.5, pos = 2,
        lab = factor.levels[which.panel[which.given]]
      )
    }
  }
  
  call <- match.call()
  if (!is.mids(x)) {
    stop("argument 'x' must be a 'mids' object", call. = FALSE)
  }
  if (is.null(x$chainMean)) {
    stop("no convergence diagnostics found", call. = FALSE)
  }
  
  mn <- x$chainMean
  sm <- sqrt(x$chainVar)
  
  # select subset of nonmissing entries
  obs <- apply(!(is.nan(mn) | is.na(mn)), 1, all)
  varlist <- names(obs)[obs]
  
  ## create formula if not given in y
  if (missing(y)) {
    formula <- as.formula(paste0(
      paste0(varlist, collapse = "+"),
      "~.it|.ms"
    ))
  } else {
    formula <- NULL
    if (is.null(y)) {
      formula <- as.formula(paste0(
        paste0(varlist, collapse = "+"),
        "~.it|.ms"
      ))
    }
    if (is.character(y)) {
      formula <- if (length(y) == 1) {
        as.formula(paste0(y, "~.it|.ms"))
      } else {
        as.formula(paste0(paste0(y, collapse = "+"), "~.it|.ms"))
      }
    }
    if (is.integer(y) || is.logical(y)) {
      vars <- varlist[y]
      formula <- if (length(vars) == 1) {
        as.formula(paste0(vars, "~.it|.ms"))
      } else {
        as.formula(paste0(paste0(vars, collapse = "+"), "~.it|.ms"))
      }
    }
    if (is.null(formula)) {
      formula <- as.formula(y)
    }
  }
  
  m <- x$m
  it <- x$iteration
  mn <- matrix(aperm(mn[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it)
  sm <- matrix(aperm(sm[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it)
  
  adm <- expand.grid(seq_len(it), seq_len(m), c("mean", "sd"))
  data <- cbind(adm, rbind(mn, sm))
  colnames(data) <- c(".it", ".m", ".ms", varlist)
  ## Dummy to trick R CMD check
  .m <- NULL
  rm(.m)
  
  tp <- xyplot(
    x = formula, data = data, groups = .m,
    type = type, lty = lty, col = col, layout = layout,
    scales = list(
      y = list(relation = "free"),
      x = list(alternating = FALSE)
    ),
    as.table = TRUE,
    xlab = "Iteration",
    ylab = "",
    strip = strip.combined,
    par.strip.text = list(lines = 0.5),
    panel = myPanel,
    ...
  )
  update(tp, par.settings = theme)
}

p1 <- plot.mids(MIfull)

#p1[, c(1:9, 11:15)]

pdf("N:/durable/people/Leo/conduct/v1/output/MI_convergencePlot.pdf")
lattice::trellis.par.set(clip = list(panel = "off"))
p1
dev.off()

png("N:/durable/people/Leo/conduct/v1/output/MI_convergencePlot.png", height = 3000, width = 2000, res = 200, type = "cairo")
lattice::trellis.par.set(clip = list(panel = "off"))
p1
dev.off()


# other plots to check the imputation
strip1 <- mice::stripplot(MIfull, NN111~.imp, pch=20, cex=2)
strip2 <- mice::stripplot(MIfull, NN112~.imp, pch=20, cex=2)
strip3 <- mice::stripplot(MIfull, NN113~.imp, pch=20, cex=2)
strip4 <- mice::stripplot(MIfull, NN114~.imp, pch=20, cex=2)
strip5 <- mice::stripplot(MIfull, NN115~.imp, pch=20, cex=2)
strip6 <- mice::stripplot(MIfull, NN116~.imp, pch=20, cex=2)
strip7 <- mice::stripplot(MIfull, NN117~.imp, pch=20, cex=2)
strip8 <- mice::stripplot(MIfull, NN118~.imp, pch=20, cex=2)


#stripPlot <- mice::stripplot(MIfull)

png(filename = "stripPlot_items.png", height = 4000, width = 2000, res = 200, type = "cairo")
gridExtra::grid.arrange(strip1, strip2, strip3, strip4, strip5, strip6, strip7, strip8, ncol = 2)
dev.off()


# other plots to inspect imputation

dens1 <- mice::densityplot(MIfull)
#[c(1:9, 11:13)]

png(filename = "N:/durable/people/Leo/conduct/v1/output/densityPlot_imp.png", height = 4000, width = 4000, res = 300, type = "cairo")
dens1
dev.off()

# plots for continuous variables used in MI
bw1 <- mice::bwplot(MIfull)[38:45]
bw2 <- mice::bwplot(MIfull)[c(1, 46:49)]

png(filename = "N:/durable/people/Leo/conduct/v1/output/bwPlot_items.png", height = 8000, width = 8000, res = 400, type = "cairo")
bw1
dev.off()

png(filename = "N:/durable/people/Leo/conduct/v1/output/bwPlot_aux.png", height = 8000, width = 8000, res = 400, type = "cairo")
bw2
dev.off()


#---------------------------------------------------------------------------------------

