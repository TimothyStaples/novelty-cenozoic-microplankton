test.richness.sensitivity <- function(all.novel.list){
  
require(lme4)
require(merTools)
  cumul.col <- rgb(0.373,0.651,0.765)

# get richness of each community in each site
all.novel.comb <- do.call("rbind", lapply(1:4, function(n){
  temp <- do.call("rbind", all.novel.list[[n]]$novel)
  temp$taxa <- c("nano", "foram", "radio", "diatom")[n]
  return(temp)
}))

all.div.com <- do.call("rbind", lapply(1:4, function(n){
    
    novel.list <- all.novel.list[[n]]
    taxa <- c("nano", "foram", "radio", "diatom")[n]
    
    print(taxa[n])
    
    all.div <- do.call("rbind", lapply(novel.list$site, function(x){
      
      temp <- data.frame(site = x,
                         bins = rownames(novel.list$ssmat[[x]]),
                         sp.rich = rowSums(novel.list$ssmat[[x]]),
                         t(estimateR(novel.list$ssmat[[x]])))
      
      temp$bin.lag <- c(NA, diff(as.numeric(as.character(temp$bins))))
      temp$taxa <- taxa
      
      return(temp)
    }))
  
  
}))

all.div.com$site.bin <- paste0(all.div.com$site, ":", all.div.com$bins)
all.novel.comb$site.bin <- paste0(all.novel.comb$site, ":", all.novel.comb$bins)

all.div <- merge(all.div.com, 
                 all.novel.comb[, !colnames(all.novel.comb) %in% 
                                  c("site", "bins", "taxa", "bin.lag")], 
                 all.x=TRUE, all.y=FALSE,
                 by.x="site.bin", by.y="site.bin", sort=FALSE)

all.div$cat <- as.factor(all.div$cat)
all.div$cat.bef <- as.factor(all.div$cat.bef)

# run models 
richness.prob.models <- lapply(c("instant", "cumul", "novel"), function(x){
  print(x)
  temp.data <- all.div  
  temp.data$success <- temp.data$cat==x  
  temp.data$log.rich <- log(temp.data$sp.rich)
  
  temp.model <- glmer(success ~ log.rich + I(log.rich^2) + (1|taxa/site), 
                      data=temp.data,
                      family=binomial)
  
  return(temp.model)
})

rich.prob.preds <- lapply(richness.prob.models,
                          function(temp){
                            
                            pred.df <- data.frame(log.rich = seq(min(temp@frame$log.rich),
                                                                 max(temp@frame$log.rich),
                                                                 len=200))
                            
                            # pred.df <- do.call("rbind", replicate(4, pred.df, simplify=FALSE))
                            #  pred.df$taxa <- factor(rep(c("nano", "foram", "radio", "diatom"),
                            #                             each=200),
                            #                        levels=levels(as.factor(temp$data$taxa)))
                            pred.df$taxa <- temp@frame$taxa[1]
                            pred.df$site <- temp@frame$site[1]
                            
                            pred.df <- cbind(pred.df,
                                             sum.fit = predict(temp,
                                                               newdata=pred.df,
                                                               re.form=NA),
                                             predictInterval(temp, 
                                                             newdata = pred.df,
                                                             n.sims = 999,
                                                             which="fixed",
                                                             level=0.95,
                                                             include.resid.var=FALSE,
                                                             type="linear.prediction"))
                            
                            return(pred.df)
                          })

pdf("./plots/prob by richness.pdf", height=6, width=3, useDingbats=FALSE)
par(mfrow=c(3,1), mar=c(0,0,0,0), oma=c(2.5,3.5,1,1), 
    tcl=-0.25, mgp=c(3,0.5,0), las=1, ps=10)

lapply(1:3, function(n){
  
  data <- richness.prob.models[[n]]@frame
  preds <- rich.prob.preds[[n]]
  
  summary(exp(data$log.rich))
  
  plot(x=NULL, y=NULL, xlim=c(1, 100), ylim=c(-0.005,0.17), yaxs="i",
       xlab="", ylab="", xaxt="n", log="x", yaxt="n")
  
  if(n==3){
    mtext(side=1, line=1.25, text="Species richness", cex=0.8)
    axis(side=1, mgp=c(3,0.2,0))
    axis(side=1, at=c(seq(1,10,1), seq(10,70,10)), tcl=-0.125, labels=NA)
  } else {
    axis(side=1, mgp=c(3,0.2,0), labels=NA)
    axis(side=1, at=c(seq(1,10,1), seq(10,70,10)), tcl=-0.125, labels=NA)
  }
  
  mtext(side=2, line=2.5, text="Probability", las=0, cex=0.8)
  axis(side=2, mgp=c(3,0.5,0))
  
  per.ts.prop <- tapply(data$success,
                        exp(data$log.rich),
                        function(x){sum(x) / length(x)})
  
  sample.size <- tapply(data$log.rich,
                        exp(data$log.rich),
                        length)
  
  summary(0.25 + log(sample.size+1) / max(log(sample.size))*1.5)
  
  col.mat <- col2rgb(c("red",cumul.col,"orange"))/255
  points(per.ts.prop ~ as.numeric(names(per.ts.prop)),
         cex = 0.4 + log(sample.size+1) / max(log(sample.size))*1.5,
         pch=16, col=rgb(col.mat[1,n],
                         col.mat[2,n],
                         col.mat[3,n],0.25), lwd=0.1)
  
  lines(x=exp(preds$log.rich),
        y=plogis(preds$upr),
        col="black", lwd=1, lty="31")
  
  lines(x=exp(preds$log.rich),
        y=plogis(preds$lwr),
        col="black", lwd=1, lty="31")
  
  lines(x=exp(preds$log.rich),
        y=plogis(preds$fit),
        col="black", lwd=2)
  
  text(x=0.9,
       y=relative.axis.point(0.93, "y"),
       labels=c("(A) Instantaneous Novelty",
                "(B) Cumulative Novelty",
                "(C) Novel community")[n], adj=0, font=2)
  
  if(n==2){
    legend(x=35, y=0.24,
           legend=c("1", "10", "100", "500"),
           pch=21, pt.bg="grey90", yjust=0.5, bty="n",
           pt.cex=0.4 + log(c(1,10,100,500)+1) / max(log(sample.size))*1.5,
           x.intersp=0.75, y.intersp=0.8, title="Sample\nsize")
  }
})
dev.off()

}