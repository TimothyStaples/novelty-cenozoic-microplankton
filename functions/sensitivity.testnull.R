novel.sensitivity.test.null <- function(observed.data,
                                   probability.models,
                                   turnover.models,
                                   nsims){
  
# calculate turnover probabilities from demographic models
print("Estimating taxa turnover probabilities")
cc.model.preds <- lapply(1:4, function(n){
    
    print(n)
    temp.model <- turnover.models[[n]]$model
    
    pred.df <- expand.grid(cat.bef=levels(temp.model@frame$cat.bef),
                           cat.aft=levels(temp.model@frame$cat.aft))
    pred.df$lag.scaled = 0
    
    if(n %in% 1:2){
      pred.df$edge.scaled = 0 
      pred.mat <- model.matrix(~ (cat.bef + cat.aft) + lag.scaled + edge.scaled, data=pred.df)
    } else {
      pred.mat <- model.matrix(~ (cat.bef + cat.aft) + lag.scaled, data=pred.df)  
    }
    
    glht.preds <- glht(temp.model, linfct=pred.mat)
    glht.preds <- as.data.frame(summary(glht.preds)$test[c("coefficients", 
                                                           "sigma")])
    colnames(glht.preds) <- c("fit", "se.fit")
    #cor(cbind(plogis(glht.preds[,"fit"]), mer.preds$r.fit))
    
    model.preds <-cbind(pred.df, glht.preds)
    model.preds$comb <- as.factor(paste0(model.preds$cat.bef, model.preds$cat.aft))
    model.preds$model <- c("ext", "orig", "emig", "immig")[n]
    
    model.preds$estimate <- plogis(model.preds$fit)
    model.preds$upper <- plogis(model.preds$fit + 1.96 * model.preds$se.fit)
    model.preds$lower <- plogis(model.preds$fit - 1.96 * model.preds$se.fit)
    
    return(model.preds)
  })
  
# extract values out of prediction data-frame, and converting 68% confidence intervals
# into standard error (difference from mean prediction)
ext.prob.df <- cc.model.preds[[1]]
ext.prob.df <- ext.prob.df[ext.prob.df$cat.bef == "back",]
colnames(ext.prob.df) <- gsub("__", "", colnames(ext.prob.df))

orig.prob.df <- cc.model.preds[[2]]
orig.prob.df <- orig.prob.df[orig.prob.df$cat.bef == "back",]
colnames(orig.prob.df) <- gsub("__", "", colnames(orig.prob.df))

emig.prob.df <- cc.model.preds[[3]]
emig.prob.df <- emig.prob.df[emig.prob.df$cat.bef == "back",]
colnames(emig.prob.df) <- gsub("__", "", colnames(emig.prob.df))

immig.prob.df <- cc.model.preds[[4]]
immig.prob.df <- immig.prob.df[immig.prob.df$cat.bef == "back",]
colnames(immig.prob.df) <- gsub("__", "", colnames(immig.prob.df))

# get novel classification probabilities
class.probs <- do.call("rbind", lapply(probability.models[4:6], function(x){x$pred.df}))
rownames(class.probs) = c("instant", "cumul", "novel")

# get observed distributions (all drawn as if from poisson distribution) 
observed.data <- do.call("c", lapply(observed.data, function(x){x$ssmats}))
obs.ts.length <- sapply(observed.data, nrow)
obs.gamma.diversity <- sapply(observed.data, ncol)
obs.alpha.diversity <- unlist(sapply(observed.data, rowSums))

# simulate matrix dimensions
sim.ts.length <- rpois(nsims, mean(obs.ts.length))
sim.gamma <- rpois(nsims, mean(obs.gamma.diversity))

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

# brms automatically back-transform to probability scale, we want to draw from normal
# distribution where the standard errors are symmetrical

# simulate probability of selecting what type of community classification appears
class="back"
ext.prob <- unlist(ext.prob.df[ext.prob.df$cat.aft == class, c("fit", "se.fit")])
orig.prob <- unlist(orig.prob.df[orig.prob.df$cat.aft == class, c("fit", "se.fit")])
emig.prob <- unlist(emig.prob.df[emig.prob.df$cat.aft == class, c("fit", "se.fit")])
immig.prob <- unlist(immig.prob.df[immig.prob.df$cat.aft == class, c("fit", "se.fit")])

for(i in 1:(nrow(sim.mat)-1)){

# for each row, find col # of taxa that have not appeared yet
if(i>1){ not.appeared <- colSums(sim.mat[1:i,]) == 0
} else { not.appeared = sim.mat[i,] == 0 }
  
# set extinctions to 0
if(i==1){extinct = rep(FALSE, ncol(sim.mat))}
  
# ORIGINATIONS ##
  
# origination can only occur for taxa that haven't appeared before  
if(sum(not.appeared)>0){
sim.orig <- rbinom(sum(not.appeared), 1, plogis(rnorm(1, orig.prob[1], orig.prob[2])))
sim.mat[i+1, sim.orig==1] = 1
}

# EXTINCTIONS ##
# extinctions can happen to any non-extinct taxa, except in the time-point they originate
sim.ext <- rbinom(sum(sim.mat[i,][!extinct]), 1, plogis(rnorm(1, ext.prob[1], ext.prob[2])))
sim.mat[i+1, sim.ext==1] = 0

# add on extinctions to our records so they cannot re-immigrate later
extinct[sim.mat[i,]==1][sim.ext==1]=TRUE

# EMIGRATION ##

# non-extinct present taxa can disappear
if(sum(!extinct) > 0){
sim.emig <- rbinom(sum(sim.mat[i,][!extinct]), 1, plogis(rnorm(1, emig.prob[1], emig.prob[2])))
sim.mat[i+1, sim.mat[i,]==1 & !extinct][sim.emig==1] = 0
}

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
  
if(sum(immig.cand) > 0){
sim.immig <- rbinom(sum(immig.cand), 1, plogis(rnorm(1, immig.prob[1], immig.prob[2])))
sim.mat[i+1,immig.cand] = 1 
}
  
}

  
}

rownames(sim.mat) <- nrow(sim.mat):1

return(sim.mat[rowSums(sim.mat) > 0 , colSums(sim.mat) > 0])

})

# test for novel communities
sim.novel.list <- lapply(1:nsims, function(n){

print(paste0("Detecting novelty in sim ", n, "..."))
ssmat <- sim.mat.list[[n]]

if(is.null(nrow(ssmat)) | is.null(ncol(ssmat))){return(NULL)}

temp <- identify.novel.gam(ssmat, alpha=0.05, metric="jaccard",
                     site=1, plot=FALSE)

temp$sim <- n

return(temp)

})
sim.novel.df <- do.call("rbind", sim.novel.list)

cumul.col <- rgb(0.373,0.651,0.765)

pdf(date.wrap("./plots/simulated novel plot", ".pdf"),
    height=4, width=4.25, useDingbats = FALSE)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(3,3,1,1), ps=8, mgp=c(3,0.5,0), tcl=-0.25, las=1)

sim.novel.df$color <- c("grey80", cumul.col, "red", "orange")[as.factor(sim.novel.df$cat)]
sim.cats <- split(sim.novel.df, f=as.factor(sim.novel.df$cat))

hulls <- lapply(sim.cats, function(x){
  x <- x[complete.cases(x), c("raw.min.dist","seq.dist")]
  return(x[chull(x),])})

sapply(1:4,
       function(n){
         
         plot(x=NULL, y=NULL, xlim=c(0,1.02), ylim=c(0,1.02), xaxt="n", yaxt="n", xlab="", ylab="")
         
         if(n %in% 3:4){axis(side=1, mgp=c(3,0,0))}
         if(n %in% c(1,3)){axis(side=2)}
         
         if(n==3){
         par(xpd=NA)
         mtext(side=1, line=1, at = par("usr")[2], text="Dissimilarity to preceding time-point")
         mtext(side=2, line=1.75, at = par("usr")[4], text="Dissimilarity to most similar past time-point", las=0)
         par(xpd=FALSE)
         }
         
        # if(n != 1){polygon(x=hulls[[1]][,2], y=hulls[[1]][,1], 
        #                    lwd=0.5, col=rgb(0.7,0.7,0.7,0.5))}
        # if(n != 2){polygon(x=hulls[[2]][,2], y=hulls[[2]][,1], 
        #                     lwd=0.5, col=rgb(1,0,0,0.5))
         
         # remove overlapping points to reduce plot size
          un.cats <- do.call("rbind", sim.cats[(1:4)[-n]])
          un.cats <- un.cats[!duplicated(paste0(round(un.cats$raw.min.dist, 2), ":",
                                                round(un.cats$seq.dist, 2))),]
          
         points(un.cats$raw.min.dist ~ un.cats$seq.dist, pch=16, col="grey90")
         
         sub.cats <- sim.cats[[n]]
         sub.cats <- sub.cats[!duplicated(paste0(round(sub.cats$raw.min.dist, 2), ":",
                                               round(sub.cats$seq.dist, 2))),]
         
         with(sub.cats, points(raw.min.dist ~ seq.dist, pch=21, bg=as.character(color)),
              lwd=0.5, cex=0.75)
         
         perc <- round(nrow(sim.cats[[n]]) / nrow(sim.novel.df), 3) * 100
           
         text(x=relative.axis.point(0.02, "x"),
              y=relative.axis.point(0.95, "y"),
              labels=paste0("(",LETTERS[n],") ", c("Background",
                                                   "Cumulative novelty",
                                                   "Instantaneous novelty",
                                                   "Novel community")[n],
                            ": ", perc, "%"),
              adj=0, font=2)
       })
dev.off()

return(sim.novel.df)

}