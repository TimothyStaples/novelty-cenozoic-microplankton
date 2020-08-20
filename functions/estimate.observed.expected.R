estimate.observed.expected <- function(prob.model.list,
                                       novel.list,
                                       dist.draws = 1e6){
  
  # extract observed classification probabilities from models
  print("Preparing observed transition data...")
  
  # condense community data into single data-frame.
  comm.dat <- do.call('rbind', lapply(1:4, function(n){
          x <- novel.list[[n]]
          temp <- do.call("rbind", x$novel)
          temp$taxa <- c("nano", "foram", "radio", "diatom")[n]
          return(temp)
          }))
  
  # set site to be a taxa-specific factor
  comm.dat$site <- as.factor(paste0(comm.dat$taxa,":", comm.dat$site))

  # observed transition data-frame
  obs.df <- as.data.frame(table(comm.dat$cat.bef,
                                comm.dat$cat,
                                comm.dat$site))
  colnames(obs.df) <- c("cat.bef", "cat.aft", "site", "obs")
  obs.df$taxa <- substr(obs.df$site, 1, regexpr(":", obs.df$site)-1)
  obs.df$taxa.site <- paste0(obs.df$taxa, ":", obs.df$site)
  
  # get all transition combinations
  trans.list <- expand.grid(levels(obs.df$cat.bef), 
                            levels(obs.df$cat.aft),
                            stringsAsFactors = FALSE)
  
  # model observed probabilities of transition occurring
  print("Running observed transition models...")
  trans.preds <- lapply(1:16, function(n){
    
    target.trans <- trans.list[n,]
    
    # subset just frequency of target transition
    trans.df <- do.call("rbind", lapply(split(obs.df, f=obs.df$site), 
                                        function(x){
                                          
                                          x$cat.bef <- as.character(x$cat.bef)
                                          x$cat.aft <- as.character(x$cat.aft)
                                          
                                          data.frame(site=x$site[1],
                                                     taxa=x$taxa[1],
                                                     success = x$obs[x$cat.bef == target.trans[1,1] &
                                                                       x$cat.aft == target.trans[1,2]],
                                                     failure = sum(x$obs[x$cat.bef != target.trans[1,1] |
                                                                           x$cat.aft != target.trans[1,2]]))
                                          
                                        }))
    
    trans.m <- glmer(cbind(success, failure) ~ 1 + (1|taxa),
                   data = trans.df, family = binomial,
                   glmerControl(optimizer ="bobyqa"))
    
    
    return(list(transition = target.trans,
                model = trans.m,
                fixed = summary(trans.m)$coefficients))
    
  })
  
  # calculate expected probabilities from shift probability models
  print("Calculating expected transition probabilities...")
  base.probs <- do.call("rbind", lapply(prob.model.list$random.prob.models,
                                        function(x){x$pred.df}))
  
  rownames(base.probs) = 1:nrow(base.probs)
  base.probs <- as.data.frame(base.probs)
  base.probs$cat = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")
  base.probs <- base.probs[-(1:2),]
  
  print("Drawing from expected and observed distributions...")
  obs.exp.probs.df <- do.call("rbind", lapply(1:16, function(n){
    
    target.trans <- cbind(trans.preds[[n]]$transition,
                          trans.preds[[n]]$fixed)[,1:4]
    colnames(target.trans)[1:2] = c("cat.bef", "cat.aft")
    
    trans1 <- base.probs[base.probs$cat == target.trans[1,1],]
    trans2 <- base.probs[base.probs$cat == target.trans[1,2],]
    
    final.df <- target.trans
    colnames(final.df)[3:4] = paste0("obs.", colnames(final.df)[3:4])
    
    # draw randomly from probability distribution on logit-scale,
    # then calculate combined probabilities, return distribution
    # of expected probabilities
    exp.probs.bef <- rnorm(1e6, trans1$Estimate, trans1$`Std. Error`)
    exp.probs.aft <- rnorm(1e6, trans2$Estimate, trans2$`Std. Error`)
    
    exp.probs.bef.raw = plogis(exp.probs.bef)
    exp.probs.aft.raw = plogis(exp.probs.aft)
    
    exp.probs.comb = exp.probs.bef.raw * exp.probs.aft.raw
    exp.probs.df = data.frame(exp.mean = mean(exp.probs.comb),
                              exp.upper = quantile(exp.probs.comb, 0.975),
                              exp.lower = quantile(exp.probs.comb, 0.025))
    
    obs.probs = rnorm(1e6, final.df$obs.Estimate, final.df$`obs.Std. Error`)
    
    ratio.probs =  plogis(obs.probs) / exp.probs.comb
    
    ratio.probs.df = data.frame(ratio.mean = mean(ratio.probs),
                                ratio.upper = quantile(ratio.probs, 0.975),
                                ratio.lower = quantile(ratio.probs, 0.025))
    
    return(cbind(final.df, exp.probs.df, ratio.probs.df))
    
  }))
  
  obs.exp.probs.df$non.zero <- obs.exp.probs.df$ratio.lower > 1 & obs.exp.probs.df$ratio.upper > 1 |
    obs.exp.probs.df$ratio.lower < 1 & obs.exp.probs.df$ratio.upper < 1
  
  obs.exp.probs.df <- obs.exp.probs.df[order(obs.exp.probs.df$non.zero,
                                             obs.exp.probs.df$ratio.mean),]
  
  return(obs.exp.probs.df)
  
}