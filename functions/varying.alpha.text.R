varying.alpha.test <- function(probs,
                               threshold.prob.list,
                               circle.1 = 0.01,
                               circle.2 = 0.05,
                               circle.3 = 0.1,
                               ylims = c(0.005,max(probs)+0.005)){
  
  require(gamm4)
  
  combined.list <- do.call("rbind", 
                           lapply(1:length(threshold.prob.list), function(n){
                             
                             sub <- threshold.prob.list[[n]]$data
                             sub$threshold <- thresholds[n]
                             return(sub)
                           }))
  
  I.prob <- gamm4(cbind(instant, non.instant) ~ s(threshold, bs="cr", k=10),
                         random = ~(1|taxa/site),
                  family=binomial, 
                  data=combined.list[combined.list$all.instant > 0,])
  C.prob <- gamm4(cbind(cumul, non.cumul) ~ s(threshold, bs="cr", k=10),
                  random = ~(1|taxa/site),
                  family=binomial, 
                  data=combined.list[combined.list$all.cumul > 0,])
  N.prob <- gamm4(cbind(novel, non.novel) ~ s(threshold, bs="cr", k=10),
                  random = ~(1|taxa/site),
                  family=binomial, 
                  data=combined.list[combined.list$all.cumul > 0,])
  
  large.overlap <- gamm4(cbind(novel, all.instant) ~ s(threshold, bs="cr", k=10),
                         random = ~(1|taxa/site),
                         family=binomial, 
                         data=combined.list[combined.list$all.instant > 0,])
  new.overlap <- gamm4(cbind(novel, all.cumul) ~ s(threshold, bs="cr", k=10),
                       random = ~(1|taxa/site),
                       family=binomial, 
                       data=combined.list[combined.list$all.cumul > 0,])
  
  pred.df <- data.frame(threshold = unique(sort(c(seq(min(thresholds), max(thresholds), 0.01),
                                                  seq(min(thresholds), max(thresholds), len=200)))))
  
  pred.list <- lapply(list(I.prob, C.prob, N.prob, large.overlap, new.overlap),
                      function(model){
                        
                        new.preds <- cbind(pred.df, predict(model$gam, newdata=pred.df, se.fit=TRUE))
                        new.preds$upper <- plogis(new.preds$fit + 1.96 * new.preds$se.fit)
                        new.preds$lower <- plogis(new.preds$fit - 1.96 * new.preds$se.fit)
                        new.preds$fit <- plogis(new.preds$fit)
                        
                        return(new.preds)
                        
                      })
  
  cumul.col <- c(95,166,195)/255
  novel.col <- col2rgb("orange") /255
  
  pdf(date.wrap("./plots/novel overlap gams", ".pdf"),
      height=4.5, width=5.5, useDingbats = FALSE)
  
  # split.screen(rbind(c(0.1,0.9,0.175,0.975),
  #                    c(0.3,0.6,0.175,0.975), # venn circles
  #                    c(0.1,0.3,0.175,0.975), # number prob plot
  #                    c(0.6,0.75,0.175,0.975), # I overlap prob
  #                    c(0.75,0.9,0.175,0.975))) # C overlap prb
  
  split.screen(rbind(c(0.15,0.95,0.1,0.925),
                     c(0.15,0.95,0.5,0.725), # venn circles
                     
                     c(0.15,0.95,0.725,0.925), # number prob plot
                     c(0.15,0.95,0.3,0.5), # I overlap prob
                     c(0.15,0.95,0.1,0.3))) # C overlap prb
  
  screen(1)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, ylim=c(0,1), xlim=ylims, axes=FALSE,
       xlab="", ylab="", xaxs="i", yaxs="i")
  abline(v=c(circle.1, circle.2, circle.3), lty="32", col="grey70")
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, ylim=c(0,1), xlim=ylims, 
       axes=FALSE, xaxt="n", yaxt="n", 
       xlab="", ylab="", xaxs="i", yaxs="i")
  
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       border=NA, col="white")
  
  segments(y0=par("usr")[3], y1=relative.axis.point(0.35, "y"),
           x0=circle.1, x1=0.015, lty="32", col="grey70")
  
  segments(y0=relative.axis.point(0.65, "y"), y1=par("usr")[4],
           x0=0.015, x1=circle.1, lty="32", col="grey70")
  
  novel.venn.circles(y.pos = 0.5, 
                     x.pos = c(0.013,0.0165), 
                     radius = c(0.004,0.003))
  
  circle.1.inst <- pred.list[[1]][pred.list[[1]]$threshold == circle.1,]
  circle.1.cumul <- pred.list[[2]][pred.list[[2]]$threshold == circle.1,]
  circle.1.novel <- pred.list[[3]][pred.list[[3]]$threshold == circle.1,]
  
  segments(y0=par("usr")[3], y1=relative.axis.point(0.35, "y"),
           x0=circle.2, x1=circle.2, lty="32", col="grey70")
  
  segments(y0=relative.axis.point(0.65, "y"), y1=par("usr")[4],
           x0=circle.2, x1=circle.2, lty="32", col="grey70")
  novel.venn.circles(y.pos = 0.5, 
                     x.pos = c(0.0475,0.0525), 
                     radius = c(0.006,0.005))
  
  segments(y0=par("usr")[3], y1=relative.axis.point(0.35, "y"),
           x0=circle.3, x1=0.08, lty="32", col="grey70")
  segments(y0=relative.axis.point(0.65, "y"), y1=par("usr")[4],
           x0=0.08, x1=circle.3, lty="32", col="grey70")
  
  novel.venn.circles(y.pos = 0.5, 
                     x.pos = c(0.085,0.0925), 
                     radius = c(0.009,0.0075))
  
  close.screen(2)
  
  screen(3)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, ylim=c(0,0.06), xlim=ylims, axes=FALSE,
       xlab="", ylab="", xaxs="i", yaxs="i")
  
  axis(side=2, at=seq(0,0.06,0.02))
  axis(side=2, at=seq(0,0.06, 0.005), tcl=-0.125, labels=NA)
  mtext(side=2, text="Probability", las=0, line=2, las=0)
  axis(side=3, mgp=c(3,0.2,0))
  axis(side=3, at=seq(0,0.1, 0.01), tcl=-0.125, labels=NA)
  mtext(side=3, text=expression(alpha*" outlier detection threshold"), line=1)
  
  polygon(y=c(pred.list[[1]]$upper, rev(pred.list[[1]]$lower)),
          x=c(pred.list[[1]]$threshold, rev(pred.list[[1]]$threshold)),
          border=NA, col=rgb(1,0,0,0.25))
  lines(pred.list[[1]]$fit ~ pred.list[[1]]$threshold, col="red")
  
  polygon(y=c(pred.list[[2]]$upper, rev(pred.list[[2]]$lower)),
          x=c(pred.list[[2]]$threshold, rev(pred.list[[2]]$threshold)),
          border=NA, col=rgb(cumul.col[1],cumul.col[2],cumul.col[3],0.25))
  
  lines(pred.list[[2]]$fit ~ pred.list[[2]]$threshold, 
        col=rgb(cumul.col[1],cumul.col[2],cumul.col[3],1))
  
  polygon(y=c(pred.list[[3]]$upper, rev(pred.list[[3]]$lower)),
          x=c(pred.list[[3]]$threshold, rev(pred.list[[3]]$threshold)),
          border=NA, col=rgb(novel.col[1], novel.col[2], novel.col[3], 0.4))
  
  text(x=max(pred.list[[1]]$threshold),
       y=c(rev(pred.list[[1]]$fit)[1],
           rev(pred.list[[2]]$fit)[1],
           rev(pred.list[[3]]$fit)[1]),
       labels=c("I","C","N"), col=c("red", rgb(cumul.col[1], 
                                               cumul.col[2], 
                                               cumul.col[3]), "orange"),
       font=2, pos=4, offset=0.15)
  
  lines(pred.list[[3]]$fit ~ pred.list[[3]]$threshold, col="orange")
  
  text(x = relative.axis.point(0.025, "x"),
       y = relative.axis.point(0.9, "y"),
       labels = "(A)", font=2, adj=0.5)
  
  box()
  close.screen(3)
  
  screen(4)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, ylim=c(0.15,0.30), xlim=ylims, axes=FALSE,
       xlab="", ylab="", xaxs="i", yaxs="i")
  axis(side=2, at=seq(0.05,0.5,0.05), las=1)
  axis(side=2, at=seq(0,0.4, 0.01), tcl=-0.125, labels=NA)
  
  mtext(side=2, text="Probability of", las=0, line=2.5)
  mtext(side=2, text=bquote(bold("N")*phantom(" given ")*phantom(bold("I"))), las=0, line=1.75, col="orange")
  mtext(side=2, text=bquote(phantom(bold("N"))*" given "*phantom(bold("I"))), las=0, line=1.75)
  mtext(side=2, text=bquote(phantom(bold("N"))*phantom(" given ")*bold("I")), las=0, line=1.75, col="red")
  
  axis(side=1, labels=NA)
  axis(side=1, at=seq(0,0.1, 0.01), tcl=-0.125, labels=NA)
  
  polygon(y=c(pred.list[[4]]$upper, rev(pred.list[[4]]$lower)),
          x=c(pred.list[[4]]$threshold, rev(pred.list[[4]]$threshold)),
          border=NA, col=rgb(1,0,0,0.25))
  lines(pred.list[[4]]$fit ~ pred.list[[4]]$threshold, col="red")
  
  line.seg<-which(pred.list[[4]]$threshold==0.05)
  segments(y0=pred.list[[4]]$upper[line.seg],
           x0=pred.list[[4]]$threshold[1],
           y1=pred.list[[4]]$upper[line.seg],
           x1=max(pred.list[[4]]$threshold), col=rgb(1,0,0,0.4), lty="11")
  segments(y0=pred.list[[4]]$lower[line.seg],
           x0=pred.list[[4]]$threshold[1],
           y1=pred.list[[4]]$lower[line.seg],
           x1=max(pred.list[[4]]$threshold), col=rgb(1,0,0,0.4), lty="11")
  box()
  
  text(x = relative.axis.point(0.025, "x"),
       y = relative.axis.point(0.9, "y"),
       labels = "(B)", font=2, adj=0.5)
  
  close.screen(4)
  
  screen(5)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, ylim=c(0.22,0.39), xlim=ylims, axes=FALSE,
       xlab="", ylab="", xaxs="i", yaxs="i")
  
  axis(side=2, at=seq(0.05,0.5,0.05))
  axis(side=2, at=seq(0,0.4, 0.01), tcl=-0.125, labels=NA)
  
  mtext(side=2, text="Probability of", las=0, line=2.5)
  mtext(side=2, text=bquote(bold("N")*phantom(" given ")*phantom(bold("C"))), las=0, line=1.75, col="orange")
  mtext(side=2, text=bquote(phantom(bold("N"))*" given "*phantom(bold("C"))), las=0, line=1.75)
  mtext(side=2, text=bquote(phantom(bold("N"))*phantom(" given ")*bold("C")), las=0, line=1.75, 
        col=rgb(cumul.col[1], 
                cumul.col[2], 
                cumul.col[3]))
  
  axis(side=1, mgp=c(3,0,0))
  axis(side=1, at=seq(0,0.1, 0.01), tcl=-0.125, labels=NA)
  mtext(side=1, text=expression(alpha*" outlier detection threshold"), line=0.75)
  
  polygon(y=c(pred.list[[5]]$upper, rev(pred.list[[5]]$lower)),
          x=c(pred.list[[5]]$threshold, rev(pred.list[[5]]$threshold)),
          border=NA, col=rgb(cumul.col[1],cumul.col[2],cumul.col[3],0.25))
  
  lines(pred.list[[5]]$fit ~ pred.list[[5]]$threshold, 
        col=rgb(cumul.col[1],cumul.col[2],cumul.col[3]))
  box()
  
  text(x = relative.axis.point(0.025, "x"),
       y = relative.axis.point(0.9, "y"),
       labels = "(C)", font=2, adj=0.5)
  
  close.screen(5)
  
  dev.off()
  
}
