# empty results data frame
results_sem_trio <- data.frame(matrix(NA, nrow = length(vars)*7, ncol = 10))


# loop to run all univariate models

MIdata$KJONN <- as.numeric(MIdata$KJONN)

for (i in 1:length(vars)) {
  
  sem_trio =  c(
    #measurement model
    "conduct_lf =~ NN111 + NN112 + NN113 + NN114 + NN115 + NN116 + NN117 + NN118",
    
    # genetic relatedness
    paste0(pspc_scores[3*i-2], " ~~ ft*", pspc_scores[3*i-1]," + mt*", pspc_scores[3*i], " + FAAR"),
    
    #regressions (direct effect + non-transmitted)
    paste0("conduct_lf ~ c1*", pspc_scores[3*i-2]," + fnt*", pspc_scores[3*i-1], " + mnt*", pspc_scores[3*i], " + FAAR"),
    
    # genetic transmission (indirect) effects
    "gt_f := ft*c1", 
    "gt_m := mt*c1",
    
    "total_f := gt_f + fnt",
    "total_m := gt_m + mnt"
    
  )
  
  sem_fit_trio <- sem(model = sem_trio, data = MIdata, std.lv = T)
  
  # compare with models where total_f = 
  
  # outputs
  out_sem_trio <- parameterEstimates(sem_fit_trio, standardized = T)
  
  results_sem_trio[7*i-6, ] <- out_sem_trio[out_sem_trio$label == "c1", ]
  results_sem_trio[7*i-5, ] <- out_sem_trio[out_sem_trio$label == "fnt", ]
  results_sem_trio[7*i-4, ] <- out_sem_trio[out_sem_trio$label == "mnt", ]
  results_sem_trio[7*i-3, ] <- out_sem_trio[out_sem_trio$label == "gt_f", ]
  results_sem_trio[7*i-2, ] <- out_sem_trio[out_sem_trio$label == "gt_m", ]
  results_sem_trio[7*i-1, ] <- out_sem_trio[out_sem_trio$label == "total_f", ]
  results_sem_trio[7*i, ] <- out_sem_trio[out_sem_trio$label == "total_m", ]
}