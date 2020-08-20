taxa.turnover.models.brms <- function(novel.list,
                                 test.edge.effects = FALSE,
                                 ext.cutoff = NA,
                                 orig.cutoff = NA,
                                 warmup=1000,
                                 iter=2000){
  
  # calculate bin-to-bin turnover for each taxa
  turnover.list <- lapply(1:length(novel.list), function(x.n){
    print(x.n)
    x <- novel.list[[x.n]]
    taxa <- c("nano", "foram", "radio", "diatom")[x.n]
    
    # for each site
    do.call("rbind", lapply(1:length(x$sites), function(n){
      
      temp <- local.extorig(ssmat = x$ssmats[[n]],
                            novel = x$novel[[n]],
                            site.name = x$sites[n],
                            raw.data = x$data)
      temp$taxa <- taxa
      return(temp)
    }))
  })
  turnover.df <- do.call("rbind", turnover.list)
  
  # old version including overall loss and gain
  # success.vars <- c("ext", "orig", "emig", "immig", "loss", "gain")
  # trial.vars <- c("ext.rich", "orig.rich", "ext.rich", "orig.rich", "ext.rich", "orig.rich")
  # edge.vars <- c("n.from.end", "n.from.start", "n.from.end", "n.from.start", NA, NA)
   
  success.vars <- c("ext", "orig", "emig", "immig")
  trial.vars <- c("ext.rich", "orig.rich", "ext.rich", "orig.rich")
  edge.vars <- c("n.from.end", "n.from.start", "n.from.end", "n.from.start")
  
  
  if(test.edge.effects){
    # for extinction and origination, there are likely edge effects. We need to work out
    # where those quadratic edge effects balance out
    turnover.edge.effect(turnover.df)
  }
  
  demographic.models <- lapply(1:4, function(n){
    
    temp.df <- turnover.df  
    
    print(n)
    success.var = success.vars[n]
    trial.var = trial.vars[n]
    edge.var = edge.vars[n]
    
    if(n %in% c(1,3) & 
       !is.na(ext.cutoff)){
      temp.df <- temp.df[temp.df$n.from.end >= ext.cutoff,]  
    }
    
    if(n %in% c(2,4) & 
      !is.na(orig.cutoff)){
      temp.df <- temp.df[temp.df$n.from.start >= orig.cutoff,]  
      }
    
    temp.df$success <- temp.df[,success.var]
    temp.df$trials <- temp.df[,trial.var]
    temp.df$failure <- temp.df$trial - temp.df$success
    
    temp.df$cat.bef <- as.factor(temp.df$cat.bef)
    temp.df$cat.aft <- as.factor(temp.df$cat.aft)
    
    temp.df <- temp.df[complete.cases(temp.df[,c("cat.bef","cat.aft","bin.lag")]),]
    temp.df$lag.scaled <- scale(log(temp.df$bin.lag))
    
    # generate data for model
    print("Running model...")
    if(is.na(edge.var)){
      
    dem.model <- brm(success | trials(trials) ~ 
                       (cat.bef + cat.aft) + lag.scaled +
                        (1|taxa/site),
                       data = temp.df, family = binomial, cores = 4,
                       control = list(adapt_delta = 0.99),
                     iter=iter, warmup=warmup,
                       prior = c(set_prior("normal(0,10)", class = "b"),
                                 set_prior("normal(0,10)", class = "Intercept"),
                                 set_prior("cauchy(0,5)", class = "sd")))
      
    } else {
      
      temp.df$edge.scaled <- scale(temp.df[, edge.var])
      
      dem.model <- brm(success | trials(trials) ~ 
                       (cat.bef + cat.aft) + lag.scaled + edge.scaled +
                       (1|taxa/site),
                       data = temp.df, family = binomial, cores = 4,
                       control = list(adapt_delta = 0.99),
                       iter=iter, warmup=warmup,
                       prior = c(set_prior("normal(0,10)", class = "b"),
                                 set_prior("normal(0,10)", class = "Intercept"),
                                 set_prior("cauchy(0,5)", class = "sd")))
      
    }
    
    cat.marg <- marginal_effects(dem.model, plot=FALSE, 
                             effects="cat.bef",
                             conditions=data.frame(cat.aft=levels(temp.df$cat.aft),
                                                   trials =1),
                             method = "fitted",
                             re_formula = NA)
    
    return(list(model = dem.model,
                coefs = summary(dem.model)$fixed,
                marg_effects = cat.marg$cat.bef))
    
  })
  
  return(demographic.models)
}
