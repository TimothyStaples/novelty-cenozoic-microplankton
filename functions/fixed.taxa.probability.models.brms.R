fixed.taxa.prob.models.brms <- function(novel.freq.df, 
                                   test.model=FALSE,
                                   iter=10000,
                                   warmup=1000){

require(brms)
  
model.env <- environment()  
  
taxa.prob.models <- lapply(1:6, function(n){
  
  success.var = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")[n]
  failure.var = c("non.all.instant", "non.all.cumul", "non.back", 
                  "non.instant", "non.cumul", "non.novel")[n]
  
  print(success.var)
  
  temp.df <- novel.freq.df  
  temp.df$success = temp.df[, success.var]
  temp.df$failure = temp.df[, failure.var]
  temp.df$trials <- temp.df$success + temp.df$failure
  
  if(n==1){
  
  temp.prob.m <- brm(success | trials(trials) ~ -1 + taxa,
                     data = temp.df, family = binomial, cores = 4,
                     control = list(adapt_delta = 0.99),
                     iter=iter, warmup=warmup,
                     prior = c(set_prior("normal(0,10)", class = "b")))
  
  assign("taxa.prob.m", temp.prob.m, envir = model.env)
  
  } else {
    
    temp.prob.m <- update(taxa.prob.m, newdata = temp.df,
                          cores = 4,
                          control = list(adapt_delta = 0.99),
                          iter=iter, warmup=warmup)
    
  }

  plot(temp.prob.m)
  
  # model tests
  if(test.model){
    print("Running post-model tests...")
  }
  
  # group predictions
  pred.df <- summary(temp.prob.m)$fixed
  
  return(list(model=temp.prob.m,
              pred.df = pred.df))
  
})

return(taxa.prob.models)

}