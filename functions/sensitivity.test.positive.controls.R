novel.sensitivity.test.controls <- function(observed.data,
                                   probability.models,
                                   turnover.models,
                                   nsims){
  
# calculate turnover probabilities from demographic models
# observed.data <- neptune.novel.cut
# probability.models <- probability.models$random.prob.models
# turnover.models <- demographic.df

print("Estimating taxa turnover probabilities")
cc.model.preds <- lapply(1:4, function(n){
    
    print(n)
    temp.model <- turnover.models[[n]]$model
    
    covariate.names <- colnames(coefficients(temp.model)[[1]])
    covariate.names <- covariate.names[!grepl("Intercept|cat", covariate.names)]
    pred.df <- expand.grid(cat.bef=levels(temp.model@frame$cat.bef),
                           cat.aft=levels(temp.model@frame$cat.aft))
    
    covariate.df <- do.call("cbind", lapply(covariate.names, function(x){
      temp <- data.frame(y = rep(0, nrow(pred.df)))
      colnames(temp) = x
      return(temp)
    }))
    
    pred.df <- cbind(pred.df, covariate.df)
    pred.df$site = "1"
    pred.df$bins = "1"
    pred.df$taxa = "1"
    
    mer.preds <- cbind(pred.df,
                       man.fit = predict(temp.model, newdata=pred.df, re.form=NA),
                       predictInterval(temp.model,
                                       newdata = pred.df,
                                              n.sims = 999,
                                              which="fixed",
                                              level=0.68,
                                              include.resid.var=FALSE,
                                              type="linear.prediction"))
    
  })

# extract values out of prediction data-frame, and converting 68% confidence intervals
# into standard error (difference from mean prediction)
ext.prob.df <- cc.model.preds[[1]]
ext.prob.df <- ext.prob.df[ext.prob.df$cat.bef == "back",]
ext.prob.df[,"se"] = abs(ext.prob.df[,"fit"] - ext.prob.df[,"upr"])

orig.prob.df <- cc.model.preds[[2]]
orig.prob.df <- orig.prob.df[orig.prob.df$cat.bef == "back",]
orig.prob.df[,"se"] = abs(orig.prob.df[,"fit"] - orig.prob.df[,"upr"])

emig.prob.df <- cc.model.preds[[3]]
emig.prob.df <- emig.prob.df[emig.prob.df$cat.bef == "back",]
emig.prob.df[,"se"] = abs(emig.prob.df[,"fit"] - emig.prob.df[,"upr"])

immig.prob.df <- cc.model.preds[[4]]
immig.prob.df <- immig.prob.df[immig.prob.df$cat.bef == "back",]
immig.prob.df[,"se"] = abs(immig.prob.df[,"fit"] - immig.prob.df[,"upr"])

# get observed distributions (all drawn as if from poisson distribution) 
observed.data <- do.call("c", lapply(observed.data, function(x){x$ssmats}))
obs.ts.length <- sapply(observed.data, nrow)
obs.gamma.diversity <- sapply(observed.data, ncol)
obs.alpha.diversity <- unlist(sapply(observed.data, rowSums))

# simulate matrix dimensions
sim.ts.length <- rpois(nsims, mean(obs.ts.length))
sim.gamma <- rpois(nsims, mean(obs.gamma.diversity))

# get novel classification probabilities
class.probs <- do.call("rbind", lapply(probability.models[4:6],
                                       function(x){x$pred.df}))
rownames(class.probs) = c("instant", "cumul", "novel")

# identify which time-points will be treated as other classifications
novel.locs <- lapply(sim.ts.length, function(x){
  
  instant = 1 ; cumul = 1 ; novel = 1
  
  # make sure we don't end up overlapping classifications
  while(sum(instant %in% cumul, instant %in% novel, instant %in% novel) != 0){
  instant = sample(2:(x-1), ceiling(plogis(rnorm(1, class.probs["instant", 1], class.probs["instant", 2]))*x))
  cumul = sample(2:(x-1), ceiling(plogis(rnorm(1, class.probs["cumul", 1], class.probs["cumul", 2]))*x))
  novel = sample(2:(x-1), ceiling(plogis(rnorm(1, class.probs["novel", 1], class.probs["novel", 2]))*x))
}
  
  return(list(instant = instant, cumul = cumul, novel = novel))
})

# for each simulation
print("Creating simulated matrices")
sim.mat.list <- lapply(1:nsims, function(n){
  
print(paste0("Creating sim ", n, "..."))
  
# seed matrix with initial species drawn from alpha diversity distribution,
# avoiding zeros
init.div = 0
while(init.div == 0){
init.div = rpois(1, mean(obs.alpha.diversity))
}

# set up empty matrix
sim.mat <- matrix(0, nrow=sim.ts.length[n], ncol=sim.gamma[n])

# seed initial community
sim.mat[1,sample(1:ncol(sim.mat), init.div)] = 1

# now for each time-point, we pick existing species to see if they disappear, and
# then in the succeding time-point see if any extra species appear.

# This method treats our background model estimates as Gaussian on the link 
# (log odds) scale, draws a random probability within those bounds, and then
# tests whether each event is successful for each possible taxa, and iterates
# along the time-series.

for(i in 1:(nrow(sim.mat)-1)){
  
# pull out the probabilities for this time-point, based on whether it's been flagged
# as novelty
class = "back"
if(i %in% novel.locs[[n]]$instant){class = "instant"}
if(i %in% novel.locs[[n]]$cumul){class = "cumul"}
if(i %in% novel.locs[[n]]$novel){class = "novel"}

ext.prob = unlist(ext.prob.df[ext.prob.df$cat.aft == class, c("fit", "se")])
orig.prob = unlist(orig.prob.df[orig.prob.df$cat.aft == class, c("fit", "se")])
emig.prob = unlist(emig.prob.df[emig.prob.df$cat.aft == class, c("fit", "se")])
immig.prob = unlist(immig.prob.df[immig.prob.df$cat.aft == class, c("fit", "se")])

# for each row, find col # of taxa that have not appeared yet
if(i>1){not.appeared <- colSums(sim.mat[1:i,]) == 0} else {not.appeared = sim.mat[i,] == 0}
  
# set extintions to 0
if(i==1){extinct = rep(FALSE, ncol(sim.mat))}
  
# ORIGINATIONS ##
  
# origination can only occur for taxa that haven't appeared before  
sim.orig <- rbinom(sum(not.appeared), 1, plogis(rnorm(1, orig.prob[1], orig.prob[2])))
sim.mat[i+1, sim.orig==1] = 1

# EXTINCTIONS ##
sim.ext <- rbinom(sum(sim.mat[i,][!extinct]), 1, plogis(rnorm(1, ext.prob[1], ext.prob[2])))
sim.mat[i+1, sim.ext==1] = 0

# add on extinctions so they can be excluded later
extinct[sim.mat[i,]==1][sim.ext==1]=TRUE

# EMIGRATION ##

# non-extinct present taxa can disappear
sim.emig <- rbinom(sum(sim.mat[i,][!extinct]), 1, plogis(rnorm(1, emig.prob[1], emig.prob[2])))
sim.mat[i+1, sim.mat[i,]==1 & !extinct][sim.emig==1] = 0

# non-emigrations stay present
sim.mat[i+1, sim.mat[i,]==1 & !extinct][sim.emig==0] = 1

# IMMIGRATION

# immigration only possible for i > 1
if(i>1){

# candidates for immigration must be:
# - not extinct
# - have previously appeared in the time-series
# - must not be present currently
immig.cand <- sim.mat[i,] == 0 & !extinct & !not.appeared
  
sim.immig <- rbinom(sum(immig.cand), 1, plogis(rnorm(1, immig.prob[1], immig.prob[2])))
  
sim.mat[i+1,immig.cand] = 1 
  
}

  
}

rownames(sim.mat) <- nrow(sim.mat):1

return(sim.mat[rowSums(sim.mat) > 0 , colSums(sim.mat) > 0])

})

# test for novel communities
sim.novel.list <- lapply(1:nsims, function(n){

print(n)
ssmat <- sim.mat.list[[n]]

if(is.null(nrow(ssmat)) | is.null(ncol(ssmat))){return(NULL)}

temp <- identify.novel.gam(ssmat, alpha=0.05, metric="jaccard",
                     site=1, plot=FALSE)

temp$sim <- n

return(temp)

})
sim.novel.df <- do.call("rbind", sim.novel.list)

return(sim.novel.df)

}