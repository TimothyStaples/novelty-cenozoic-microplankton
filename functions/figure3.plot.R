figure3.plot <- function(turnover.models, 
                         plot.name,
                         left.xlims,
                         right.xlims,
                         ylims){
  
  print("Predicting from models...")
  
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
  
  light.cols <- rbind(c(205,205,205),
                      c(207,228,237),
                      c(255,179,179),
                      c(255,228,179)) / 255
  light.cols <- rgb(light.cols[,1], light.cols[,2], light.cols[,3])
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  library(shape)
  library(plotrix)
  pdf(date.wrap(paste0("./plots/cc combined - temp vs perm (",
                       plot.name,
                       ")"),
                ".pdf"), 
      height=4.75, width=11, useDingbats = FALSE)
  
  framemat <- rbind(c(0.075,0.45,0.1,0.95),
                    c(0.525,0.9,0.1,0.95),
                    c(0.9,1,0.25,0.75))
  split.screen(framemat)
  
  screen(1)
  
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim=left.xlims, ylim=ylims,
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  
  abline(a=0, b=1, lwd=2, col="grey60", lty="31")

  text(y = 0.18, x = 0.07, adj=0.5, font=1, cex=1, col="grey60",
       labels="Background\nshifts")
  text(x = 0.12, y = 0.51, adj=0, font=1, cex=1, col="grey60",
       labels="Shifts to\ninstantaneous\nnovelty")
  text(x = 0.05, y = 0.4, adj=1, font=1, cex=1, col="grey60",
       labels="Shifts to\ncumulative\nnovelty")
  text(x = 0.15, y = 0.45, adj=0.5, font=1, cex=1, col="grey60",
       labels="Shifts to\nnovel community")
  
  custom.circle(x=0.04, y=0.15, r=1.55, col=rgb(0.5,0.5,0.5,0.15),
                screen=framemat[2,], border=NA)
  
  custom.circle(x=0.07, y=0.4, r=1.55, col=rgb(0.5,0.5,0.5,0.15),
                screen=framemat[2,], border=NA)
  
  custom.circle(x=0.12, y=0.4, r=1.55, col=rgb(0.5,0.5,0.5,0.15),
                screen=framemat[2,], border=NA)
  
  custom.circle(x=0.15, y=0.4, r=1.55, col=rgb(0.5,0.5,0.5,0.15),
                screen=framemat[2,], border=NA)
  
  axis(side=1, mgp=c(3,0.2,0))
  axis(side=2)
  mtext(side=1, line=1.25, text="Probability of local extinction")
  mtext(side=2, line=2, text="Probability of emigration", las=0)
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.965, "y"),
       adj=0, labels="(A) Local extinction and emigration", font=2)
  
  segments(y0=cc.model.preds[[3]]$estimate,
           y1=cc.model.preds[[3]]$estimate,
           x0=cc.model.preds[[1]]$upper,
           x1=cc.model.preds[[1]]$lower,
           col="black", lwd=1)
  
  segments(x0=cc.model.preds[[1]]$estimate,
           x1=cc.model.preds[[1]]$estimate,
           y0=cc.model.preds[[3]]$upper,
           y1=cc.model.preds[[3]]$lower,
           col="black", lwd=1)
  
  sapply(1:dim(cc.model.preds[[1]])[1], function(x){
    print(x)
    
    arrow.shape(x=cc.model.preds[[1]]$estimate[x],
                y=cc.model.preds[[3]]$estimate[x],
                r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[1]]$cat.aft[x]],
                shape="circle",
                border=NA, plot=TRUE)
    
    arrow.shape(x=cc.model.preds[[1]]$estimate[x],
                y=cc.model.preds[[3]]$estimate[x],
                r=0.2, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[1]]$cat.bef[x]],
                plot=TRUE, border="black", add.arrow=TRUE)
    
    arrow.shape(x=cc.model.preds[[1]]$estimate[x],
                y=cc.model.preds[[3]]$estimate[x],
                r=0.2, screen=framemat[1,], rads=c(0, 2),
                shape="circle",
                plot=TRUE, border="black", add.arrow=FALSE)
    
  })
  
  # # novel "expectations"
  # # add diff from background
  # x.back <- cc.model.preds[[1]]$estimate[cc.model.preds[[1]]$cat.aft=="back"]
  # y.back <- cc.model.preds[[3]]$estimate[cc.model.preds[[3]]$cat.aft=="back"]
  # 
  # cc.x <- x.back + 
  #         (cc.model.preds[[1]]$fit[cc.model.preds[[1]]$cat.aft=="instant"] - x.back) +
  #         (cc.model.preds[[1]]$fit[cc.model.preds[[1]]$cat.aft=="cumul"] - x.back)
  # 
  # cc.y <- y.back +
  #         (cc.model.preds[[3]]$fit[cc.model.preds[[3]]$cat.aft=="instant"] - y.back) +
  #         (cc.model.preds[[3]]$fit[cc.model.preds[[3]]$cat.aft=="cumul"] - y.back)
  # 
  # left.cols <- light.cols
  # 
  # sapply(1:length(cc.x), function(x){
  #   print(x)
  #   
  #   arrow.shape(x=cc.x[x],
  #               y=cc.y[x],
  #               r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
  #               col=left.cols[4],
  #               shape="circle",
  #               border=NA, plot=TRUE)
  #   
  #   arrow.shape(x=cc.x[x],
  #               y=cc.y[x],
  #               r=0.2, screen=framemat[1,],
  #               rads=c(0.5, 1.5), lwd=0.5,
  #               shape="circle",
  #               col=left.cols[x],
  #               plot=TRUE, border="grey70", 
  #               add.arrow=TRUE, lty="solid")
  #   
  #   arrow.shape(x=cc.x[x],
  #               y=cc.y[x],
  #               r=0.2, screen=framemat[1,], rads=c(0, 2),
  #               shape="circle",
  #               plot=TRUE, border="grey70", 
  #               add.arrow=FALSE,
  #               lty="solid")
  #   
  # })
  # 
  
  box()
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=right.xlims, ylim=ylims,
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  abline(a=0, b=1, lwd=2, col="grey60", lty="31")
  
  text(y = 0.38, x = 0.045, adj=0.5, font=1, cex=1, col="grey60",
       labels="Background\nshifts")
  text(x = 0.115, y = 0.61, adj=0, font=1, cex=1, col="grey60",
       labels="Shifts to\ninstantaneous\nnovelty")
  text(x = 0.175, y = 0.375, adj=0.5, font=1, cex=1, col="grey60",
       labels="Shifts to\ncumulative\nnovelty")
  text(x = 0.35, y = 0.41, adj=0.5, font=1, cex=1, col="grey60",
       labels="Shifts to\nnovel community")
  
  draw.ellipse(x=0.045, y=0.23, a=0.03, b=0.12, 
               border=NA, col=rgb(0.5,0.5,0.5,0.15),
               angle=0)
  
  draw.ellipse(x=0.075, y=0.605, a=0.035, b=0.15, 
               border=NA, col=rgb(0.5,0.5,0.5,0.15),
               angle=0)
  
  custom.circle(x=0.175, y=0.205, r=1.55, col=rgb(0.5,0.5,0.5,0.15),
                screen=framemat[2,], border=NA)
  
  draw.ellipse(x=0.35, y=0.25, a=0.09, b=0.13, 
               border=NA, col=rgb(0.5,0.5,0.5,0.15),
               angle=360-15)
  
  axis(side=1, mgp=c(3,0.2,0))
  axis(side=2)
  
  mtext(side=1, line=1.25, text="Probability of local origination")
  mtext(side=2, line=2, text="Probability of immigration", las=0)
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.965, "y"),
       adj=0, labels="(B) Local origination and immigration", font=2)
  
  # UNCOMMENT THESE LINES TO GENERATE 95% CIS ON PLOT
  segments(y0=cc.model.preds[[4]]$estimate,
           y1=cc.model.preds[[4]]$estimate,
           x0=cc.model.preds[[2]]$upper,
           x1=cc.model.preds[[2]]$lower,
           col="black", lwd=1)
  
  segments(x0=cc.model.preds[[2]]$estimate,
           x1=cc.model.preds[[2]]$estimate,
           y0=cc.model.preds[[4]]$upper,
           y1=cc.model.preds[[4]]$lower,
           col="black", lwd=1)
  
  sapply(1:dim(cc.model.preds[[3]])[1], function(x){
    print(x)
    
    arrow.shape(x=cc.model.preds[[2]]$estimate[x],
                y=cc.model.preds[[4]]$estimate[x],
                r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[2]]$cat.aft[x]],
                shape="circle",
                border=NA, plot=TRUE)
    
    arrow.shape(x=cc.model.preds[[2]]$estimate[x],
                y=cc.model.preds[[4]]$estimate[x],
                r=0.2, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[2]]$cat.bef[x]],
                plot=TRUE, border="black", add.arrow=TRUE)
    
    arrow.shape(x=cc.model.preds[[2]]$estimate[x],
                y=cc.model.preds[[4]]$estimate[x],
                r=0.2, screen=framemat[1,], rads=c(0, 2),
                shape="circle",
                plot=TRUE, border="black", add.arrow=FALSE)
    
  })
  box()
  
  # # novel "expectations"
  # # add diff from background
  # x.back <- cc.model.preds[[2]]$estimate[cc.model.preds[[2]]$cat.aft=="back"]
  # y.back <- cc.model.preds[[4]]$estimate[cc.model.preds[[4]]$cat.aft=="back"]
  # 
  # cc.x <- x.back + 
  #   (cc.model.preds[[2]]$estimate[cc.model.preds[[2]]$cat.aft=="instant"] - x.back) +
  #   (cc.model.preds[[2]]$estimate[cc.model.preds[[2]]$cat.aft=="cumul"] - x.back)
  # 
  # cc.y <- y.back +
  #   (cc.model.preds[[4]]$estimate[cc.model.preds[[4]]$cat.aft=="instant"] - y.back) +
  #   (cc.model.preds[[4]]$estimate[cc.model.preds[[4]]$cat.aft=="cumul"] - y.back)
  # 
  # left.cols <- light.cols
  # 
  # sapply(1:length(cc.x), function(x){
  #   print(x)
  #   
  #   arrow.shape(x=cc.x[x],
  #               y=cc.y[x],
  #               r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
  #               col=left.cols[4],
  #               shape="circle",
  #               border=NA, plot=TRUE)
  #   
  #   arrow.shape(x=cc.x[x],
  #               y=cc.y[x],
  #               r=0.2, screen=framemat[1,],
  #               rads=c(0.5, 1.5), lwd=0.5,
  #               shape="circle",
  #               col=left.cols[x],
  #               plot=TRUE, border="grey70", 
  #               add.arrow=TRUE, lty="solid")
  #   
  #   arrow.shape(x=cc.x[x],
  #               y=cc.y[x],
  #               r=0.2, screen=framemat[1,], rads=c(0, 2),
  #               shape="circle",
  #               plot=TRUE, border="grey70", 
  #               add.arrow=FALSE,
  #               lty="solid")
  #   
  # })
  # 
  
  close.screen(2)
  
  screen(3)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  
  par(xpd=NA)
  
  taxa.pos <- rev(seq(0.375,0.675,len=4))
  
  # text(x=0.575, y=0.75, adj=0.5,
  #      labels=bquote(bold(underline("Taxa"))), font=2, cex=0.8)
  # sapply(taxa.pos, function(y){
  #   arrow.shape(x=0.125,
  #               y=y,
  #               r=0.2, screen=framemat[3,], rads=c(1.5,0.5),
  #               col="white",
  #               shape=c("circle", "square", "diamond", "triangle")[which(taxa.pos==y)],
  #               border="black", plot=TRUE)
  #   
  #   arrow.shape(x=0.125,
  #               y=y,
  #               r=0.2, screen=framemat[3,], rads=c(0.5, 1.5),
  #               col="white",
  #               shape=c("circle", "square", "diamond", "triangle")[which(taxa.pos==y)],
  #               border="black", plot=TRUE, add.arrow=TRUE, lwd=0.5)
  #   
  #   arrow.shape(x=0.125,
  #               y=y,
  #               r=0.2, screen=framemat[3,], rads=c(0,2),
  #               shape=c("circle", "square", "diamond", "triangle")[which(taxa.pos==y)],
  #               border="black", plot=TRUE)
  #   
  #   
  # })
  # 
  # text(x=0.13, pos=4, y=taxa.pos,
  #      labels=c("Nannoplankton", "Foraminifera", "Radiolarians", "Diatoms"), 
  #      cex=0.8, offset=0.75)
  
  par(lheight=0.85)
  text(x=0.15, y=1.05, labels="Preceding\ncommunity", adj=0, cex=0.8,
       lheight=0.5)
  
  text(x=0.85, y=0.925, labels="Succeeding\ncommunity", adj=1, cex=0.8,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.55, y=1.2,
              r=0.25, screen=framemat[3,], rads=c(1.5,0.5),
              col="white", shape="circle",
              border="black", plot=TRUE)
  
  arrow.shape(x=0.45, y=1.2,
              r=0.25, screen=framemat[3,], rads=c(0.5,1.5),
              col="white", shape="circle",
              border="black", plot=TRUE, add.arrow = TRUE)
  
  segments(x0=c(0.125, 0.075, 0.875, 0.925),
           x1=c(0.075, 0.075, 0.925, 0.925),
           y0=c(1.05, 1.05, 0.925, 0.925),
           y1=c(1.05, 1.15, 0.925, 1.15))
  
  Arrows(x0=c(0.075, 0.925),
         x1=c(0.3, 0.7),
         y0=c(1.15, 1.15),
         y1=c(1.2, 1.2),
         arr.length = 0.1, arr.width = 0.1,
         arr.type = "triangle")
  
  rect.pos <- rev(seq(-0.275,0.05, len=4))
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[1] + 0.03, ybottom=rect.pos[1] - 0.03, col="grey35")
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[2] + 0.03, ybottom=rect.pos[2] - 0.03, col="red")
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[3] + 0.03, ybottom=rect.pos[3] - 0.03, col=cumul.col)
  
  rect(xleft=0.05, xright=0.2,
       ytop=rect.pos[4] + 0.03, ybottom=rect.pos[4] - 0.03, col="orange")
  
  par(lheight=0.7)
  text(x=0.125, y=rect.pos - c(0, 0.015, 0.015, 0.015), pos=4, offset=0.75, cex=0.8,
       labels=c("Background", "Instantaneous\nnovelty", 
                "Cumulative\nnovelty", "Novel\ncommunity"))
  par(lheight=1)
  
  text(x=0.575, y=0.215, adj=0.5,
       labels=bquote(bold(underline("Community"))), font=2, cex=0.8)
  text(x=0.575, y=0.15, adj=0.5,
       labels=bquote(bold(underline("classification"))), font=2, cex=0.8)
  
  par(xpd=FALSE)
  
  close.screen(3)
  
  dev.off()
  
}