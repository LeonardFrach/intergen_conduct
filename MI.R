
# 
# install and load package
library("mice")
library("dplyr")

# set working directory
#setwd("/cluster/projects/p471/people/Leo/conduct/")

# read Rdata file which was created in TSD and copied to the cluster
load("MIdata.Rdata")


MIdata <- MIdata %>% mutate(across(NN111:NN118, as.numeric))
class(MIdata$NN111)

# check time before running it
pre <- timestamp()

MIfull <- mice(MIdata[, 3:dim(MIdata)[2]],
                         m = 100, maxit = 30, seed = 123) #note, exlude preg ID and BARN_NR here

# time after running
post <- timestamp()


# create list object with element for each of 100 imputed datasets, then add preg id and birth order to each:

mice.imp <- NULL

for(i in 1:100) {
  mice.imp[[i]] <- complete(MIfull, action=i, include=FALSE)
  mice.imp[[i]]$PREG_ID_2306 <- MIdata$PREG_ID_2306
  mice.imp[[i]]$BARN_NR <- MIdata$BARN_NR
}


# write file as output 
save(MIfull, file = "MIfull.RData")
save(mice.imp, file = "MIfull_with_IDs.RData")
