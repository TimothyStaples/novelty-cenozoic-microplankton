causes.cons.models.binom.edge.fixed <- function(success.var,
                                          trial.var,
                                          edge.var = NA,
                                          turnover.df){

  turnover.df$success <- turnover.df[,success.var]
  turnover.df$trial <- turnover.df[,trial.var]
  turnover.df$failure <- turnover.df$trial - turnover.df$success
  
  turnover.df$cat.bef <- as.factor(turnover.df$cat.bef)
  turnover.df$cat.aft <- as.factor(turnover.df$cat.aft)
  
  turnover.df <- turnover.df[complete.cases(turnover.df[,c("cat.bef","cat.aft","bin.lag")]),]
  turnover.df$lag.scaled <- scale(log(turnover.df$bin.lag))
  
  # generate data for model
  print("Running model...")
  if(is.na(edge.var)){
  cc.model <-  glmer(cbind(success, failure) ~ (cat.bef + cat.aft) * taxa + lag.scaled + 
                     (1|site),
                     family=binomial, data=turnover.df,
                     control=glmerControl(optimizer="bobyqa"))
  } else {
    turnover.df$edge.scaled <- scale(turnover.df[, edge.var])
    
    cc.model <-  glmer(cbind(success, failure) ~ (cat.bef + cat.aft) * taxa + 
                       lag.scaled + edge.scaled + (1|site),
                       family=binomial, data=turnover.df,
                       control=glmerControl(optimizer="bobyqa"))
    
  }
  
  return(list(model = cc.model,
              coefs = summary(cc.model)$coefficients,
              data = turnover.df))
  
}