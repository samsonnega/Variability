cortAUC_full <- cortImputed %>% ## cortImputed = the dataset with imputed values replacing the missing values
  group_by(Year, Litter, SquirrelID, Sex, Mum.Treat, complete, TreatmentDuring, JulianDate, age)
%>% 
  ## grouped by all these variables so they would stay in the dataset after summarizing
  summarize(AUC_imputed = auc(sample,Cort3, from = 2, to = 4, type = "spline", na.rm = TRUE)) 
## auc calculates the area under the curve using 'sample' as the key (base = 1, dex = 2, 30 = 3, 60 = 4) and 'Cort3' as the variable