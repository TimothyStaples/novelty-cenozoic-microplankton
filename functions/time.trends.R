time.trends <- function(novel.list, max.k=-1){
  
  require(mgcv)
  
  div.df <- do.call("rbind", lapply(1:length(novel.list), function(n1){
    # for each taxa
    taxa.samp <- novel.list[[n1]]$raw_samp
    
    # in each longhurst province
    do.call("rbind", lapply(1:length(taxa.samp), function(n){
      
      temp.samp <- taxa.samp[[n]]
      
      # ignore longhurst provinces with no sampling
      if(is.null(dim(temp.samp))){return(NULL)}
      
      temp.samp$site <- names(taxa.samp)[n]
      temp.samp$taxa <- c("nano", "foram", "radio", "diatom")[n1]
      return(temp.samp)
    }))
    
  }))
  div.df$taxa <- as.factor(div.df$taxa)
  div.df$site <- as.factor(div.df$site)
  
  # Diversity model
  print("Running diversity model")
  div.model <- gam(divSIB ~ s(mid.point, bs="cr", by=taxa, k=max.k) + taxa + s(site, bs="re"), family=nb, data=div.df)
  
  pred.df <- data.frame(mid.point=rep(seq(min(div.df$mid.point), max(div.df$mid.point), len=200), 4),
                        site="a",
                        taxa=rep(factor(levels(div.df$taxa), levels=levels(div.df$taxa)), each=200))
  div.pred <- cbind(pred.df, predict(div.model, newdata=pred.df, se.fit=TRUE))
  div.pred$upper <- exp(div.pred$fit + 1.96 * div.pred$se.fit)
  div.pred$lower <- exp(div.pred$fit - 1.96 * div.pred$se.fit)
  div.pred$fit <- exp(div.pred$fit)
  div.pred.list <- split(div.pred, f=div.pred$taxa)
  
  # Sampling completeness model
  print("Running sample completeness model")
  
  # sample completness model (by interval by longhurst province)
  # continuous proportion data, modelled as beta distribution
  samp.model <- gam(samp3t ~ s(mid.point, bs="cr", by=taxa, k=max.k) + taxa + s(site, bs="re"), 
                    family=betar, data=div.df)
  
  pred.df <- data.frame(mid.point=rep(seq(min(div.df$mid.point), max(div.df$mid.point), len=200), 4),
                        site="a",
                        taxa=rep(factor(levels(div.df$taxa), levels=levels(div.df$taxa)), each=200))
  samp.pred <- cbind(pred.df, predict(samp.model, newdata=pred.df, se.fit=TRUE))
  samp.pred$upper <- plogis(samp.pred$fit + 1.96 * samp.pred$se.fit)
  samp.pred$lower <- plogis(samp.pred$fit - 1.96 * samp.pred$se.fit)
  samp.pred$fit <- plogis(samp.pred$fit)
  samp.pred.list <- split(samp.pred, f=samp.pred$taxa)
  
  # Novelty through time model
  print("Running novelty through time model")
  zachos <- read.csv("./raw.datafiles/zachos2001.csv")
  zachos$genus <- gsub(" ", "", as.character(zachos$genus))
  
  bins <- seq(0,66, 0.1)
  
  zachos.100k <- cut(zachos$age.mya, breaks = bins)
  
  zachos.agg <- data.frame(bin = levels(zachos.100k),
                           d18O = as.vector(tapply(zachos$d18O,
                                                   zachos.100k,
                                                   mean, na.rm=TRUE)))
  
  age.bins <- as.character(zachos.agg$bin)
  zachos.agg$age <- round(as.numeric(substr(age.bins, regexpr(",", age.bins)+1,
                                            nchar(age.bins)-1)) - 0.5*(bins[2]-bins[1]),2)
  
  zachos.agg <- zachos.agg[order(zachos.agg$age, decreasing=TRUE),]
  zachos.agg <- zachos.agg[!is.na(zachos.agg$d18O),]
  zachos.agg$diff.d18O <- c(NA, diff(zachos.agg$d18O))
  zachos.agg$abs.diff.d18O <- abs(zachos.agg$diff.d18O)
  zachos.agg$bin.lag <- round(abs(c(NA, diff(zachos.agg$age))),2)
  
  # model novel probability over time ####
  
  novel.df <- do.call("rbind", lapply(1:length(novel.list), function(n){
    
    x <- novel.list[[n]]
    
    x.df <- do.call("rbind", x$novel)
    x.df$taxa <- c("nano", "foram", "radio", "diatom")[n]
    
    return(x.df)    
  }))
  novel.df$bin.num <- as.numeric(as.character(novel.df$bins))
  
  novel.prop <- data.frame(novel.prop = tapply(novel.df$novel,
                                               cut(novel.df$bin.num, breaks=bins),
                                               mean))
  novel.prop$age <- rowMeans(cbind(bins[-1],
                                   bins[-length(bins)]))
  
  zachos.novel <- merge(novel.prop, zachos.agg,
                        all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  zachos.novel.long <- merge(novel.df, zachos.agg,
                             by.x = "bins", by.y="age",
                             all.x=TRUE, all.y=FALSE, sort=FALSE)
  zachos.novel.long$scaled.n <- as.vector(scale(zachos.novel.long$bin.n))
  zachos.novel.long$scaled.lag <- as.vector(scale(zachos.novel.long$bin.lag.x))
  zachos.novel.long$taxa <- as.factor(zachos.novel.long$taxa)
  zachos.novel.long$site = as.factor(zachos.novel.long$site)
  
  novel.trends <- gam(novel ~ s(bin.num, bs="cr", k=max.k, by=taxa) + scaled.n + scaled.lag + taxa + s(site, bs="re"),
                      family=binomial, data = zachos.novel.long)
  
  zachos.agg$scaled.lag <- as.vector(scale(zachos.agg$age))
  
  temp.trends <- gam(abs.diff.d18O ~ s(age, bs="cr", k=max.k), data=zachos.agg)
  
  novel.pred.df <- lapply(1:4, function(n){
    
    temp.pred <-  data.frame(bin.num = seq(min(div.df$mid.point),
                                           max(div.df$mid.point),0.1),
                             scaled.n=0,
                             scaled.lag=0,
                             taxa = levels(zachos.novel.long$taxa)[n],
                             site="a")
    
    cbind(temp.pred,
          as.data.frame(predict(novel.trends, newdata=temp.pred, se.fit=TRUE)))
  })
  
  zachos.pred.df <- data.frame(age = seq(min(div.df$mid.point),
                                         max(div.df$mid.point),0.1),
                               scaled.lag=0)
  zachos.pred.df <- cbind(zachos.pred.df,
                          as.data.frame(predict(temp.trends, newdata=zachos.pred.df, se.fit=TRUE)))
  
  # PLOT ####
  
  pdf("./plots/time_trends.pdf", height=6, width=5, useDingbats=FALSE)
  
 split.screen(rbind(c(0.14,0.99,0.75,0.975),
                     c(0.14,0.99,0.525,0.75),
                     c(0.14,0.99,0.3, 0.525),
                     c(0.14, 0.99, 0.075, 0.3)))
  
  screen(1)
  par(mar=c(0,0,0,0), tcl=-0.25, mgp=c(3,0.5,0), ps=8, las=1)
  plot(x=NULL, y=NULL, xlim=c(68,0.5), ylim=c(0,0.075),
       axes=FALSE, xlab="", ylab="", xaxs="i")
  
  axis(side=1, at=seq(0,65,5), labels=NA)
  axis(side=1, at=seq(0,65,1), labels=NA, tcl=-0.125)
  
  axis(side=2, at=seq(0,0.08, 0.02))
  axis(side=2, at=seq(0,0.08, 0.01), tcl=-0.125, labels=NA)
  mtext(side=2, line=2, text="Probability of\nnovel community", las=0)
  
  lapply(1:length(novel.pred.df), function(n){
    
    temp.col = c("red", "blue", "darkgreen", "orange")[n]
    temp.light.col <- col2rgb(temp.col) / 255
    temp.alpha <-  ifelse(n==4, 0.3, 0.15)
    
    with(novel.pred.df[[n]], 
         polygon(y=c(plogis(fit + 1.96*se.fit), rev(plogis(fit - 1.96*se.fit))),
                 x = c(bin.num, rev(bin.num)),
                 lwd=2, col= rgb(temp.light.col[1], temp.light.col[2], temp.light.col[3], temp.alpha),
                 type="l", border=NA))
    
    with(novel.pred.df[[n]], 
         lines(plogis(fit) ~ bin.num, lwd=3, 
               col=temp.col, 
               type="l"))
    
    text(x=rev(novel.pred.df[[n]]$bin.num)[1],
         y=rev(plogis(novel.pred.df[[n]]$fit))[1],
         pos=2, offset=0.25,
         labels = c("D", "F", "N", "R")[n],
         col=temp.col, font=2)
    
  })
  box()
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.935, "y"),
       font=2, labels="(A)", adj=0)
  
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), tcl=-0.25, mgp=c(3,0.5,0), ps=8, las=1)
  with(zachos.pred.df,
       plot(fit ~ age, lwd=2, col="black",
            type="l", xlim=c(68,0.5), ylim=c(0,0.55),
            axes=FALSE, xlab="", ylab="", xaxs="i"))
  
  axis(side=1, at=seq(0,65,5), labels=NA)
  axis(side=1, at=seq(0,65,1), labels=NA, tcl=-0.125)
  
  axis(side=2, at=seq(0,0.5,0.1))
  axis(side=2, at=seq(0,0.5,0.05), tcl=-0.125, labels=NA)
  mtext(side=2, line=2.5, text="Temperature proxy change", las=0)
  mtext(side=2, line=1.55, text=expression("|"*Delta*" "*delta^18*"O| / 100K years"), las=0)
  
  box()
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.915, "y"),
       font=2, labels="(B)", adj=0)
  close.screen(2)
  
  screen(3)
  par(mar=c(0,0,0,0), tcl=-0.25, mgp=c(3,0.5,0), ps=8, las=1)
  plot(x=NULL, y=NULL, xlim=c(68,-0.5), ylim=c(0,33), axes=FALSE, xlab="", ylab="",
       xaxs="i", yaxs="i")
  
  axis(side=1, at=seq(0,65,5), labels=NA)
  axis(side=1, at=seq(0,65,1), labels=NA, tcl=-0.125)
  
  
  axis(side=2, mgp=c(3,0.4,0))
  mtext(side=2, line=1.5, text="Sample bin richness", las=0)
  
  lapply(1:length(div.pred.list), function(n){
    
    temp.col = c("red", "blue", "darkgreen", "orange")[n]
    temp.light.col <- col2rgb(temp.col) / 255
    temp.alpha <-  ifelse(n==4, 0.3, 0.15)
    
    with(div.pred.list[[n]], 
         polygon(y=c(upper, rev(lower)),
                 x = c(mid.point, rev(mid.point)),
                 lwd=2, col= rgb(temp.light.col[1], temp.light.col[2], temp.light.col[3], temp.alpha),
                 type="l", border=NA))
    
    with(div.pred.list[[n]], 
         lines(fit ~ mid.point, lwd=3, 
               col=temp.col, 
               type="l"))
    
    text(x=rev(div.pred.list[[n]]$mid.point)[1],
         y=rev(div.pred.list[[n]]$fit)[1],
         pos=2, offset=0.25,
         labels = c("D", "F", "N", "R")[n],
         col=temp.col, font=2)
    
  })
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.915, "y"),
       font=2, labels="(C)", adj=0)
  
  
  box()
  close.screen(3)
  
  screen(4)
  par(mar=c(0,0,0,0), tcl=-0.25, mgp=c(3,0.5,0), ps=8, las=1)
  plot(x=NULL, y=NULL, xlim=c(68,-0.5), ylim=c(0.6,1), axes=FALSE, xlab="", ylab="",
       xaxs="i", yaxs="i")
  
  axis(side=1, at=seq(0,65,5), mgp=c(3,0,0))
  axis(side=1, at=seq(0,65,1), labels=NA, tcl=-0.125)
  mtext(side=1, line=0.85, text="Millions of years BP")
  
  
  axis(side=2, mgp=c(3,0.4,0))
  mtext(side=2, line=1.5, text="Sample bin completeness", las=0)
  
  lapply(1:length(samp.pred.list), function(n){
    
    temp.col = c("red", "blue", "darkgreen", "orange")[n]
    temp.light.col <- col2rgb(temp.col) / 255
    temp.alpha <-  ifelse(n==4, 0.3, 0.15)
    
    with(samp.pred.list[[n]], 
         polygon(y=c(upper, rev(lower)),
                 x = c(mid.point, rev(mid.point)),
                 lwd=2, col= rgb(temp.light.col[1], temp.light.col[2], temp.light.col[3], temp.alpha),
                 type="l", border=NA))
    
    with(samp.pred.list[[n]], 
         lines(fit ~ mid.point, lwd=3, 
               col=temp.col, 
               type="l"))
    
    text(x=rev(samp.pred.list[[n]]$mid.point)[1],
         y=rev(samp.pred.list[[n]]$fit)[1],
         pos=2, offset=0.25,
         labels = c("D", "F", "N", "R")[n],
         col=temp.col, font=2)
    
  })
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.915, "y"),
       font=2, labels="(D)", adj=0)
  
  
  box()
  close.screen(4)
  close.screen(all.screens=TRUE)
  
  dev.off()
  
}