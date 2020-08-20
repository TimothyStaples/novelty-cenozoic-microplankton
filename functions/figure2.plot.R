figure2.plot <- function(trans.df, plot.name, ylims){
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  light.cols <- rbind(c(205,205,205),
                      c(207,228,237),
                      c(255,179,179),
                      c(255,228,179)) / 255
  light.cols <- rgb(light.cols[,1], light.cols[,2], light.cols[,3])
  
  library(shape)
  library(plotrix)
  
  pdf(date.wrap(paste0("./plots/transition null model (", 
                       plot.name, 
                       ")"), ".pdf"), 
      height=2, width=2.25, useDingbats = FALSE, colormodel = "cmyk")
  
  framemat<-rbind(c(0.19,0.99,0.175,0.99),
                  c(0.8,1,0.25,0.75))
  
  split.screen(framemat)
  
  par(mar=c(0,0,0,0), ps=6, tcl=-0.2, las=1, mgp=c(3,0,0))
  plot(x=NULL, y=NULL, xlim=log(c(2e-4,1)), ylim=ylims,
       axes=FALSE, xlab="", ylab="", yaxs="i")
  
  abline(h=0, lwd=0.75, lty="31")
  
  custom.circle(x=log(1), y=log(1), r=0.5, col=rgb(0.5,0.5,0.5,0.35),
                screen=framemat[1,], border=NA)
  
  draw.ellipse(x=log(0.05), y=log(0.775), a=0.85, b=0.35,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle= 0)
  
  draw.ellipse(x=log(0.0002), y=log(2.8), a=2, b=1,
               border=NA, col=rgb(0.5,0.5,0.5,0.35),
               angle = 0)
  
  mtext(side=1, line=0.4, cex=1,
        text="Expected probability of transition")
  mtext(side=1, line=0.8, cex=0.8,
        text="(calculated from occurrence probabilities)")
  
  par(lheight=0.85)
  mtext(side=2, line=1.1, las=0,
        text="ln(Ratio of observed to expected probability)")
  
  text(x=log(1), y=log(1.05), pos=3, adj=0.5, offset=0.8,
       labels="Background\ntransitions", col="grey60", font=1)
  
  text(x=log(0.05), y=log(0.75),  adj=0,
       labels="Transitions\nto and\nfrom\nnovelty", col="grey60", font=1)
  
  text(x=log(0.0002), y=log(2.8),  adj=0,
       labels="Transitions\nbetween\nnovelty\ncategories", col="grey60", font=1)
  
  axis(side=1,
       at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
       labels=c(0.0001, 0.001,0.01,0.1,1,10,100),
       mgp=c(3,-0.15,0))
  
  axis(side=1, at=log(c(seq(0.0001,0.001,0.0001),
                        seq(0.001,0.01,0.001),
                        seq(0.01,0.1,0.01),
                        seq(0.1,1,0.1),
                        seq(1,10,1),
                        seq(10,100,10))), tcl=-0.1, labels=NA)
  
  y.locs <- c(seq(1, 15, 0.5),
              1/seq(1, 15, 0.5))
  
  y.locs <- c(seq(1,20,1), seq(10,100,10))
  y.locs <- c(y.locs, 1/y.locs)
  
  axis(side=2, at=log(y.locs), tcl=-0.125, labels=NA)
  
  axis(side = 2, 
       at = log(c(0.5,1,2,3,5,10,20,30,50,100)),
       labels = c(0.5,1,2,3,5,10,20,30,50,100), 
       mgp=c(3,0.4,0))
  
  sapply(1:dim(trans.df)[1], function(x){
    
    if(trans.df$non.zero[x]){
    
      aft.col = c("grey35", cumul.col, "red", "orange")[as.factor(trans.df$cat.aft)[x]]
      bef.col = c("grey35", cumul.col, "red", "orange")[as.factor(trans.df$cat.bef)[x]]
      border.col = "black"
      
    } else {
      
      aft.col = light.cols[as.factor(trans.df$cat.aft)[x]]
      bef.col = light.cols[as.factor(trans.df$cat.bef)[x]]
      border.col = "grey50"
      
    }
    
    segments(y0=log(trans.df$ratio.mean[x]),
             y1=log(trans.df$ratio.mean[x]),
             x0=log(trans.df$exp.upper[x]),
             x1=log(trans.df$exp.lower[x]),
             col=ifelse(trans.df$non.zero[x],
                        "black", "grey50"), lwd=0.75)
    
    segments(y0=log(trans.df$ratio.lower[x]),
             y1=log(trans.df$ratio.upper[x]),
             x0=log(trans.df$exp.mean[x]),
             x1=log(trans.df$exp.mean[x]),
             col=ifelse(trans.df$non.zero[x],
                        "black", "grey50"), lwd=0.75)
    
    arrow.shape(x=log(trans.df$exp.mean[x]),
                y=log(trans.df$ratio.mean[x]),
                r=0.125, screen=framemat[1,], rads=c(1.5, 0.5),
                col = aft.col,
                shape="circle",
                border=NA)
    
    arrow.shape(x=log(trans.df$exp.mean[x]),
                y=log(trans.df$ratio.mean[x]),
                r=0.125, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape="circle",
                col = bef.col, 
                border=border.col, add.arrow=TRUE)
    
    arrow.shape(x=log(trans.df$exp.mean[x]),
                y=log(trans.df$ratio.mean[x]),
                r=0.125, screen=framemat[1,], rads=c(0, 2),
                shape="circle", lwd=0.5,
                plot=TRUE, border=border.col, add.arrow=FALSE)

    
  })
  box()
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  
  par(xpd=NA)
  
  taxa.pos <- rev(seq(0.375,0.675,len=4))
  
  par(lheight=0.85)
  text(x=0.1, y=1.3, labels="Preceding\nstate", adj=0, cex=0.45,
       lheight=0.5)

  text(x=0.85, y=1.135, labels="Succeeding\nstate", adj=1, cex=0.45,
       lheight=0.5)
  par(lheight=1)

  arrow.shape(x=0.55, y=1.45,
              r=0.125, screen=framemat[2,], rads=c(1.5,0.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE)

  arrow.shape(x=0.45, y=1.45,
              r=0.125, screen=framemat[2,], rads=c(0.5,1.5),
              col="white", shape="circle", lwd=0.5,
              border="black", plot=TRUE, add.arrow = TRUE)

  segments(x0=c(0.0725, 0.025, 0.875, 0.925),
           x1=c(0.025, 0.025, 0.925, 0.925),
           y0=c(1.29, 1.29, 1.125, 1.125),
           y1=c(1.29, 1.39, 1.125, 1.35), lwd=0.5)

  Arrows(x0=c(0.025, 0.925),
         x1=c(0.275, 0.725),
         y0=c(1.39, 1.35),
         y1=c(1.44, 1.44), lwd=0.5,
         arr.length = 0.05, arr.width = 0.05,
         arr.type = "triangle")
  # 
  rect.pos <- rev(seq(1.175,1.475, len=4))

  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[1] + 0.04, ybottom=rect.pos[1] - 0.04, col="grey35")

  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[2] + 0.04, ybottom=rect.pos[2] - 0.04, col="red")

  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[3] + 0.04, ybottom=rect.pos[3] - 0.04, col=cumul.col)

  rect(xleft=-1.95, xright=-1.75, lwd=0.5,
       ytop=rect.pos[4] + 0.04, ybottom=rect.pos[4] - 0.04, col="orange")

  par(lheight=0.7)
  text(x=-2.05, y=rect.pos - 0.01, pos=4, offset=0.75, cex=0.45,
       labels=c("Background", "Instantaneous novelty",
                "Cumulative novelty", "Novel community"))
  par(lheight=1)

  par(xpd=FALSE)
  close.screen(2)
  
  dev.off()
}