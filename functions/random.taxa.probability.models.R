random.taxa.prob.models <- function(novel.freq.df, 
                                   test.model=FALSE){

require(multcomp)
require(lme4)
require(merTools)
  
taxa.prob.models <- lapply(1:6, function(n){
  
  success.var = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")[n]
  failure.var = c("non.all.instant", "non.all.cumul", "non.back", 
                  "non.instant", "non.cumul", "non.novel")[n]
  
  print(success.var)
  
  temp.df <- novel.freq.df
  temp.df$success = temp.df[, success.var]
  temp.df$failure = temp.df[, failure.var]
  
  taxa.prob.m <- glmer(cbind(success, failure) ~ 1 + (1|taxa),
                         data=temp.df, family=binomial)
  
  # model tests
  if(test.model){
    print("Running post-model tests...")
    sim.res <- simulateResiduals(taxa.prob.m)
    disp.test <- testDispersion(sim.res)
    print(ifelse(disp.test$p.value <=0.05, 
                 paste0("Dispersal not okay :( -- ", round(disp.test$statistic, 2)),
                 "Dispersal okay :)"))
    
    gc()
  }
  
  # group predictions
  pred.df <- as.data.frame(summary(taxa.prob.m)$coefficients)
  pred.df$taxa.rand <- summary(taxa.prob.m)$varcor$taxa[1,1]

  write.csv(pred.df,
            date.wrap(paste0("./outputs/", success.var,
                             " random taxa model summary"),
                      ".csv")) 
  
  return(list(model=taxa.prob.m,
              pred.df = pred.df))
  
})

return(taxa.prob.models)

}