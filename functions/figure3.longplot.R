figure3.longplot <- function(turnover.models, plot.name){
  
  print("Predicting from models...")
  cc.model.preds <- lapply(1:6, function(n){
    
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
  pdf(date.wrap(paste0("./plots/cc combined long (", plot.name, ")"), ".pdf"), 
      height=4.75, width=2.25, useDingbats = FALSE)
  
  framemat <- rbind(c(0.16,0.95,0.5,0.86),
                    c(0.16,0.95,0.065,0.425),
                    c(0.16,0.95, 0.875, 0.98))
  split.screen(framemat)
  
  #                             taxon extinction and emigration ####
  screen(1)
  par(mar=c(0,0,0,0), ps=6, tcl=-0.2, mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim=c(0,0.2), ylim=c(0,0.7),
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  
  abline(a=0, b=1, lwd=0.75, lty="31")
  
  text(x = 0.2725, y = 0.415, adj=0, font=1, cex=1, col="grey60",
       labels="Transitions to\nnovel\ncommunity")
  text(y = 0.135, x = 0.08, adj=0, font=1, cex=1, col="grey60",
       labels="Background\ntransitions")
  text(x = 0.12, y = 0.535, adj=0.5, font=1, cex=1, col="grey60",
       labels="Transitions to\ncumulative or\ninstantaneous\nnovelty")
  
  custom.circle(x=0.045, y=0.135, r=0.6, col=rgb(0.5,0.5,0.5,0.35),
                screen=framemat[1,], border=NA)
  
  custom.circle(x=0.12, y=0.35, r=1.2, col=rgb(0.5,0.5,0.5,0.35),
                screen=framemat[1,], border=NA)
  
  draw.ellipse(x=0.21, y=0.4, a=0.06, b=0.09, border=NA, 
               col=rgb(0.5,0.5,0.5,0.35), angle=360-05)
  
  axis(side=1, mgp=c(3,-0.1,0))
  axis(side=2, mgp=c(3,0.3,0))
  mtext(side=1, line=0.5, text="Probability of local extinction")
  mtext(side=2, line=1.1, text="Probability of emigration", las=0)
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.955, "y"),
       adj=0, labels="(A) Local extinction and emigration", font=2)
  
  # UNCOMMENT THESE LINES TO GENERATE 95% CIS ON PLOT
  segments(y0=plogis(cc.model.preds[[3]]$fit),
           y1=plogis(cc.model.preds[[3]]$fit),
           x0=plogis(cc.model.preds[[1]]$fit + 1.96*cc.model.preds[[1]]$se.fit),
           x1=plogis(cc.model.preds[[1]]$fit - 1.96*cc.model.preds[[1]]$se.fit),
           col="black", lwd=0.75)
  
  segments(x0=plogis(cc.model.preds[[1]]$fit),
           x1=plogis(cc.model.preds[[1]]$fit),
           y0=plogis(cc.model.preds[[3]]$fit + 1.96*cc.model.preds[[3]]$se.fit),
           y1=plogis(cc.model.preds[[3]]$fit - 1.96*cc.model.preds[[3]]$se.fit),
           col="black", lwd=0.75)
  
  sapply(1:dim(cc.model.preds[[1]])[1], function(x){
    print(x)
    
    arrow.shape(x=plogis(cc.model.preds[[1]]$fit[x]),
                y=plogis(cc.model.preds[[3]]$fit[x]),
                r=0.125, screen=framemat[1,], rads=c(1.5, 0.5),
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[1]]$cat.aft[x]],
                shape="circle",
                border=NA, plot=TRUE)
    
    arrow.shape(x=plogis(cc.model.preds[[1]]$fit[x]),
                y=plogis(cc.model.preds[[3]]$fit[x]),
                r=0.125, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[1]]$cat.bef[x]],
                plot=TRUE, border="black", add.arrow=TRUE)
    
    arrow.shape(x=plogis(cc.model.preds[[1]]$fit[x]),
                y=plogis(cc.model.preds[[3]]$fit[x]),
                r=0.125, screen=framemat[1,], rads=c(0, 2),
                shape="circle", lwd=0.5,
                plot=TRUE, border="black", add.arrow=FALSE)
    
  })
  
  box()
  close.screen(1)
  
  #                             taxon origination and immigration ####
  screen(2)
  par(mar=c(0,0,0,0), ps=6, tcl=-0.2, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(0,0.35), ylim=c(0,0.82),
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  abline(a=0, b=1, lwd=0.75, lty="31")
  
  # rect(xleft=0.015, xright=0.09, ybottom=0.015, ytop=0.081,
  #      col=rgb(1,1,1,0.8), border=NA)
  
  rect(xleft=0.38, xright=0.4125, ybottom=0.375, ytop=0.44,
       col=rgb(1,1,1,0.8), border=NA)
  
  text(y = 0.38, x = 0.045, adj=0.5, font=1, cex=1, col="grey60",
       labels="Background\ntransitions")
  text(x = 0.117, y = 0.6075, adj=0, font=1, cex=1, col="grey60",
       labels="Transitions to\ninstantaneous\nnovelty")
  text(x = 0.175, y = 0.38, adj=0.5, font=1, cex=1, col="grey60",
       labels="Transitions to\ncumulative\nnovelty")
  text(x = 0.35, y = 0.41, adj=0.5, font=1, cex=1, col="grey60",
       labels="Transitions to\nnovel community")
  
  draw.ellipse(x=0.045, y=0.23, a=0.03, b=0.12, 
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=0)
  
  draw.ellipse(x=0.075, y=0.61, a=0.04, b=0.15, 
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=0)
  
  custom.circle(x=0.175, y=0.205, r=1.3, col=rgb(0.5,0.5,0.5,0.35),
                screen=framemat[2,], border=NA)
  
  draw.ellipse(x=0.35, y=0.25, a=0.09, b=0.13, 
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=360-15)
  
  axis(side=1, mgp=c(3,-0.1,0))
  axis(side=2, mgp=c(3,0.3,0))
  
  mtext(side=1, line=0.5, text="Probability of local origination")
  mtext(side=2, line=1.1, text="Probability of immigration", las=0)
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.955, "y"),
       adj=0, labels="(B) Local orgination and immigration", font=2)
  
  # UNCOMMENT THESE LINES TO GENERATE 95% CIS ON PLOT
  segments(y0=plogis(cc.model.preds[[4]]$fit),
           y1=plogis(cc.model.preds[[4]]$fit),
           x0=plogis(cc.model.preds[[2]]$fit + 1.96*cc.model.preds[[2]]$se.fit),
           x1=plogis(cc.model.preds[[2]]$fit - 1.96*cc.model.preds[[2]]$se.fit),
           col="black", lwd=0.75)
  
  segments(x0=plogis(cc.model.preds[[2]]$fit),
           x1=plogis(cc.model.preds[[2]]$fit),
           y0=plogis(cc.model.preds[[4]]$fit + 1.96*cc.model.preds[[4]]$se.fit),
           y1=plogis(cc.model.preds[[4]]$fit - 1.96*cc.model.preds[[4]]$se.fit),
           col="black", lwd=0.75)
  
  sapply(1:dim(cc.model.preds[[3]])[1], function(x){
    print(x)
    
    arrow.shape(x=plogis(cc.model.preds[[2]]$fit[x]),
                y=plogis(cc.model.preds[[4]]$fit[x]),
                r=0.125, screen=framemat[2,], rads=c(1.5, 0.5),
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[2]]$cat.aft[x]],
                shape="circle",
                border=NA, plot=TRUE)
    
    arrow.shape(x=plogis(cc.model.preds[[2]]$fit[x]),
                y=plogis(cc.model.preds[[4]]$fit[x]),
                r=0.125, screen=framemat[2,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[2]]$cat.bef[x]],
                plot=TRUE, border="black", add.arrow=TRUE)
    
    arrow.shape(x=plogis(cc.model.preds[[2]]$fit[x]),
                y=plogis(cc.model.preds[[4]]$fit[x]),
                r=0.125, screen=framemat[2,], rads=c(0, 2),
                shape="circle", lwd=0.5,
                plot=TRUE, border="black", add.arrow=FALSE)
    
  })
  box()
  close.screen(2)
  
  #                             legend ####
  screen(3)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  # axis(side=1)
  # axis(side=2)
  par(xpd=NA)
  
  par(lheight=0.85)
  text(x=-0.015, y=0.55, labels="Preceding\nstate", adj=0, cex=0.6,
       lheight=0.5)
  
  text(x=0.39, y=0.2, labels="Succeeding\nstate", adj=1, cex=0.6,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.205, y=0.9,
              r=0.2, screen=framemat[3,], rads=c(1.5,0.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE)
  
  arrow.shape(x=0.155, y=0.9,
              r=0.2, screen=framemat[3,], rads=c(0.5,1.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE, add.arrow = TRUE)
  
  segments(x0=c(-0.03,-0.05, 0.4,  0.42),
           x1=c(-0.05, -0.05, 0.42, 0.42),
           y0=c(0.55, 0.55, 0.15,  0.15),
           y1=c(0.55, 0.725, 0.15,  0.45), lwd=0.5)
  
  Arrows(x0=c(-0.05, 0.42),
         x1=c(0.075, 0.285),
         y0=c(0.725, 0.45),
         y1=c(0.85, 0.775),
         arr.length = 0.1, arr.width = 0.1,
         arr.type = "triangle", lwd=0.5)
  
  rect.pos <- rev(seq(0.135,0.785, len=4))
  rect.height <- 0.09
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[1] + rect.height, ybottom=rect.pos[1] - rect.height, 
       col="grey35", lwd=0.5)
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[2] + rect.height, ybottom=rect.pos[2] - rect.height, 
       col="red", lwd=0.5)
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[3] + rect.height, ybottom=rect.pos[3] - rect.height, 
       col=cumul.col, lwd=0.5)
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[4] + rect.height, ybottom=rect.pos[4] - rect.height, 
       col="orange", lwd=0.5)
  
  par(lheight=0.7)
  text(x=0.5, y=rect.pos, pos=4, offset=0.25, cex=0.6,
       labels=c("Background", "Instantaneous novelty",
                "Cumulative novelty", "Novel community"))
  par(lheight=1)

  par(xpd=FALSE)
  close.screen(3)
  close.screen(all.screens=TRUE)
  dev.off()

  # LOSS/GAIN ####
  close.screen(all.screens=TRUE)
  
  pdf(date.wrap(paste0("./plots/overall loss-gain (", plot.name, ")"), ".pdf"),
      height=2.65, width=2.25, useDingbats = FALSE)

  framemat <- rbind(c(0.16,0.95,0.115,0.775),
                    c(0.16,0.95, 0.775, 0.98))
  split.screen(framemat)

  screen(1)
  par(mar=c(0,0,0,0), ps=6, tcl=-0.2, mgp=c(3,0.5,0), las=1)

  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1),
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")

  abline(a=0, b=1, lwd=0.75, lty="31")

  draw.ellipse(x=0.215, y=0.305, a=0.085, b=0.15,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=0)
  draw.ellipse(x=0.53, y=0.425, a=0.115, b=0.17,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=0)
  draw.ellipse(x=0.59, y=0.73, a=0.105, b=0.14,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=0)
  draw.ellipse(x=0.72, y=0.66, a=0.105, b=0.14,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle=0)

  text(y = 0.5, x = 0.215, adj=0.5, font=1, cex=1, col="grey60",
       labels="Background\ntransitions")
  text(x = 0.475, y = 0.73, adj=1, font=1, cex=1, col="grey60",
       labels="Transitions to\ninstantaneous\nnovelty")
  text(x = 0.53, y = 0.19, adj=0.5, font=1, cex=1, col="grey60",
       labels="Transitions to\ncumulative\nnovelty")
  text(x = 0.825, y = 0.66, adj=0, font=1, cex=1, col="grey60",
       labels="Transitions to\nnovel\ncommunity")

  axis(side=1, mgp=c(3,-0.1,0))
  axis(side=2, mgp=c(3,0.3,0))
  mtext(side=1, line=0.5, text="Probability of taxon loss")
  mtext(side=2, line=1.1, text="Probability of taxon gain", las=0)
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.965, "y"),
       adj=0, labels="(A) Taxon loss and gain", font=2)

  # UNCOMMENT THESE LINES TO GENERATE 95% CIS ON PLOT
  segments(y0=plogis(cc.model.preds[[6]]$fit),
           y1=plogis(cc.model.preds[[6]]$fit),
           x0=plogis(cc.model.preds[[5]]$fit + 1.96*cc.model.preds[[5]]$se.fit),
           x1=plogis(cc.model.preds[[5]]$fit - 1.96*cc.model.preds[[5]]$se.fit),
           col="black", lwd=0.75)

  segments(x0=plogis(cc.model.preds[[5]]$fit),
           x1=plogis(cc.model.preds[[5]]$fit),
           y0=plogis(cc.model.preds[[6]]$fit + 1.96*cc.model.preds[[6]]$se.fit),
           y1=plogis(cc.model.preds[[6]]$fit - 1.96*cc.model.preds[[6]]$se.fit),
           col="black", lwd=0.75)

  sapply(1:dim(cc.model.preds[[1]])[1], function(x){
    print(x)

    arrow.shape(x=plogis(cc.model.preds[[5]]$fit[x]),
                y=plogis(cc.model.preds[[6]]$fit[x]),
                r=0.125, screen=framemat[1,], rads=c(1.5, 0.5),
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[5]]$cat.aft[x]],
                shape="circle",
                border=NA, plot=TRUE)

    arrow.shape(x=plogis(cc.model.preds[[5]]$fit[x]),
                y=plogis(cc.model.preds[[6]]$fit[x]),
                r=0.125, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col=c("grey35", cumul.col, "red", "orange")[cc.model.preds[[5]]$cat.bef[x]],
                plot=TRUE, border="black", add.arrow=TRUE)

    arrow.shape(x=plogis(cc.model.preds[[5]]$fit[x]),
                y=plogis(cc.model.preds[[6]]$fit[x]),
                r=0.125, screen=framemat[1,], rads=c(0, 2),
                shape="circle", lwd=0.5,
                plot=TRUE, border="black", add.arrow=FALSE)

  })

  box()

  close.screen(1)

  #                             legend ####
  screen(2)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  # axis(side=1)
  # axis(side=2)
  par(xpd=NA)
  
  par(lheight=0.85)
  text(x=-0.015, y=0.55, labels="Preceding\nstate", adj=0, cex=0.6,
       lheight=0.5)
  
  text(x=0.39, y=0.2, labels="Succeeding\nstate", adj=1, cex=0.6,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.205, y=0.9,
              r=0.2, screen=framemat[2,], rads=c(1.5,0.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE)
  
  arrow.shape(x=0.155, y=0.9,
              r=0.2, screen=framemat[2,], rads=c(0.5,1.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE, add.arrow = TRUE)
  
  segments(x0=c(-0.03,-0.05, 0.4,  0.42),
           x1=c(-0.05, -0.05, 0.42, 0.42),
           y0=c(0.55, 0.55, 0.15,  0.15),
           y1=c(0.55, 0.725, 0.15,  0.45), lwd=0.5)
  
  Arrows(x0=c(-0.05, 0.42),
         x1=c(0.075, 0.285),
         y0=c(0.725, 0.45),
         y1=c(0.85, 0.775),
         arr.length = 0.1, arr.width = 0.1,
         arr.type = "triangle", lwd=0.5)
  
  rect.pos <- rev(seq(0.135,0.785, len=4))
  rect.height <- 0.09
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[1] + rect.height, ybottom=rect.pos[1] - rect.height, 
       col="grey35", lwd=0.5)
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[2] + rect.height, ybottom=rect.pos[2] - rect.height, 
       col="red", lwd=0.5)
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[3] + rect.height, ybottom=rect.pos[3] - rect.height, 
       col=cumul.col, lwd=0.5)
  
  rect(xleft=0.45, xright=0.515,
       ytop=rect.pos[4] + rect.height, ybottom=rect.pos[4] - rect.height, 
       col="orange", lwd=0.5)
  
  par(lheight=0.7)
  text(x=0.5, y=rect.pos, pos=4, offset=0.25, cex=0.6,
       labels=c("Background", "Instantaneous novelty",
                "Cumulative novelty", "Novel community"))
  par(lheight=1)
  
  par(xpd=FALSE)
  close.screen(2)
  close.screen(all.screens=TRUE)
  dev.off()
  
}  