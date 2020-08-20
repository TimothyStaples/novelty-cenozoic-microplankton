estimate.observed.expected.brms <- function(prob.model.list,
                                            novel.list,
                                            iter=10000,
                                            warmup=1000){
  
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
  trans.list <- expand.grid(cat.bef = levels(obs.df$cat.bef), 
                            cat.aft = levels(obs.df$cat.aft),
                            stringsAsFactors = FALSE)
  
  print("Calculating expected transition probabilities...")
  base.prob.models <- lapply(prob.model.list$random.prob.models,
                             function(x){x$model})[3:6]
  names(base.prob.models) = c("back", "instant", "cumul", "novel")
  
  exp.trans.list <- lapply(1:16, function(n){
    
    # what transiation are we after?
    target.trans = trans.list[n,]
    
    # get posterior probability for each state
    exp.post1 <- posterior_samples(base.prob.models[target.trans[1,1]][[1]],
                                   pars = "b_Intercept")[,1]
    exp.post2 <- posterior_samples(base.prob.models[target.trans[1,2]][[1]],
                                   pars = "b_Intercept")[,1]
    
    # calculate expected transition probability by multiplying 
    # back-transformed posterior of state 1 and state 2
    # (probability of state 1 occurring, followed by state 2)
    exp.post <- plogis(exp.post1) * plogis(exp.post2)
    exp.post <- log(exp.post / (1 - exp.post))
    
  })
  
  # model observed probabilities of transition occurring
  print("Running observed transition models...")
  model.env <- environment()  
  trans.preds <- lapply(1:16, function(n){
    
    target.trans <- trans.list[n,]
    
    print(paste0(target.trans, collapse = " -> "))
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
    trans.df$trials <- trans.df$success + trans.df$failure
    
    exp.trans <- exp.trans.list[[n]]
    # int.prior <- c(mean(exp.trans), sd(exp.trans))
    # int.prior <- paste0("normal(",int.prior[1],",",int.prior[2],")")
    # model probability of transition
    if(n==1){

    trans.m <- brm(success | trials(trials) ~ 0 + Intercept + (1|taxa),
                   data = trans.df, family = binomial, cores = 4,
                   control = list(adapt_delta = 0.99),
                   sample_prior = TRUE,
                   iter=iter, warmup=warmup,
                   prior = c(set_prior("normal(0,10)", class = "b"),
                             set_prior("cauchy(0,5)", class = "sd")),
                   refresh=0)
    
    assign("trans.m.template", trans.m, envir = model.env)
    
    } else {
      
    trans.m <- update(trans.m.template,
                      newdata = trans.df,
                      cores = 4,
                      control = list(adapt_delta = 0.99),
                      iter=iter, warmup=warmup,
                      prior = c(set_prior("normal(0,10)", class = "b"),
                                  set_prior("cauchy(0,5)", class = "sd")),
                      refresh=0)
      
    }

    return(list(transition = target.trans,
                data = trans.df,
                model = trans.m,
                fixed = summary(trans.m)$fixed))
    
  })

  
  # now compare observed and expected posteriors
  trans.comp.df <- do.call("rbind", lapply(1:16, function(n){
  
  obs.post <- posterior_samples(trans.preds[[n]]$model, pars = "b_Intercept")[,1]
  exp.post <- posterior_samples(exp.trans.list[[n]])[,1]
 
  diff.post <- plogis(obs.post) / plogis(exp.post)
  
  diff.df <- data.frame(obs.mean = plogis(mean(obs.post)),
             obs.upper = plogis(quantile(obs.post, 0.975)),
             obs.lower = plogis(quantile(obs.post, 0.025)),
             
             exp.mean = plogis(mean(exp.post)),
             exp.upper = plogis(quantile(exp.post, 0.975)),
             exp.lower = plogis(quantile(exp.post, 0.025)),
             
             lower.prop = sum(diff.post <= 1) / length(diff.post),
             
             ratio.mean = mean(diff.post),
             ratio.lower = quantile(diff.post, 0.025),
             ratio.upper = quantile(diff.post, 0.975)) 
  
  rownames(diff.df) <- paste0(trans.list[n,], collapse = " -> ")
  return(diff.df)
  }))
  
  trans.comp.df$non.zero <- trans.comp.df$lower.prop <= 0.05 |
                            trans.comp.df$lower.prop >= 0.95

  trans.comp.df <- cbind(trans.list, trans.comp.df)
  
  trans.comp.df <- trans.comp.df[order(trans.comp.df$non.zero,
                                       trans.comp.df$ratio.mean),]

  return(trans.comp.df)
  
}