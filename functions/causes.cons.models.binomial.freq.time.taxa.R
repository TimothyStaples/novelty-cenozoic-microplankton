demographic.binom.glmm.time.fixed <- function(success.var,
                                     trial.var,
                                     local.causes.list,
                                     edge.var = NA,
                                     taxa.vect){
  
  all.local$success <- all.local[,success.var]
  all.local$trial <- all.local[,trial.var]
  all.local$failure <- all.local$trial - all.local$success
  
  all.local$cat.bef <- as.factor(all.local$cat.bef)
  all.local$cat.aft <- as.factor(all.local$cat.aft)
  
  all.local <- all.local[complete.cases(all.local[,c("cat.bef","cat.aft","bin.lag")]),]
  
  # generate data for model
  print("Running model...")
  
  cc.model <-  glmer(cbind(success, failure) ~ (cat.bef + cat.aft) * taxa + bin.lag + (1|site/bins),
                     family=binomial, data=all.local,
                     control=glmerControl(optimizer="bobyqa"))
  summary(cc.model)
  
  return(list(model = cc.model,
              coefs = summary(cc.model)$coefficients,
              data = all.local))
  
}