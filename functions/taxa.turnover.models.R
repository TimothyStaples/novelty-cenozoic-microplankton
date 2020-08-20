taxa.turnover.models <- function(novel.list,
                                 test.edge.effects = FALSE,
                                 ext.cutoff = NA,
                                 orig.cutoff = NA){
  
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
  
  success.vars <- c("ext", "orig", "emig", "immig", "loss", "gain")
  trial.vars <- c("ext.rich", "orig.rich", "ext.rich", "orig.rich", "ext.rich", "orig.rich")
  edge.vars <- c("n.from.end", "n.from.start", NA, NA, NA, NA)
  
  if(test.edge.effects){
    # for extinction and origination, there are likely edge effects. We need to work out
    # where those quadratic edge effects balance out
    turnover.edge.effect(turnover.df)
    
  }
  
  demographic.models <- lapply(1:6, function(n){
    
    sub.turnover <- turnover.df
    
    print(n)
    success.var = success.vars[n]
    trail.var = trial.vars[n]
    edge.var = edge.vars[n]
    
    if(n %in% c(1,3) & 
       !is.na(ext.cutoff)){
      sub.turnover <- turnover.df[turnover.df$n.from.end >= ext.cutoff,]  
    }
    
    if(n %in% c(2,4) & 
      !is.na(orig.cutoff)){
      sub.turnover <- turnover.df[turnover.df$n.from.start >= orig.cutoff,]  
      }
          
    model.list <- causes.cons.models.binom.edge(success.var = success.var,
                                  trial.var = trail.var,
                                  edge.var = edge.var,
                                  turnover.df = sub.turnover)
    
    coefs <- as.data.frame(summary(model.list$model)$coefficients)
    coefs$taxa.rand <- summary(model.list$model)$varcor$taxa[1,1]
    coefs$site.rand <- summary(model.list$model)$varcor$`site:taxa`[1,1]
    
    write.csv(coefs,
              date.wrap(paste0("./outputs/", success.var,
                               "demographic model summary"),
                        ".csv")) 
    
    return(model.list)
    
  })
  
  return(demographic.models)
}