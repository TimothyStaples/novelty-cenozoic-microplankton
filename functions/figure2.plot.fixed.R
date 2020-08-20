figure2.plot.fixed <- function(trans.df, plot.name, xlims, ylims){
  
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
                       ")"),
                ".pdf"), 
      height=5.5, width=7.5, useDingbats = FALSE)
  
  framemat<-rbind(c(0.10,0.47,0.525,0.95),
                  c(0.47,0.84,0.525,0.95),
                  c(0.1,0.47,0.125,0.525),
                  c(0.47,0.84,0.125,0.525),
                  c(0.85,1,0.35,0.75))
  
  # framemat<-rbind(c(0.10,0.8,0.135,0.95),
  #                 c(0.8,1,0.25,0.75))
  # 
  split.screen(framemat)
  
  sapply(1:4, function(n){
    
  screen(n)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, las=1, mgp=c(3,0.2,0))  
  
  sub.trans <- trans.df[trans.df$taxa == levels(trans.df$taxa)[n],]
  
  plot(x=NULL, y=NULL, xlim = xlims, ylim = ylims,
       axes=FALSE, xlab="", ylab="", yaxs="i")
  
  abline(h=0, lwd=2,col="grey60", lty="31")
  
  # custom.circle(x=log(1), y=log(1), r=0.5, col=rgb(0.5,0.5,0.5,0.15),
  #               screen=framemat[n,], border=NA)
  # 
  # draw.ellipse(x=log(0.05), y=log(0.775), a=0.85, b=0.35,
  #              border=NA, col=rgb(0.5,0.5,0.5,0.15),
  #              angle= 0)
  # 
  # draw.ellipse(x=log(0.0002), y=log(2.8), a=2, b=1,
  #              border=NA, col=rgb(0.5,0.5,0.5,0.15),
  #              angle = 0)
  # 
  if(n == 3){
  mtext(side=1, line=2.25, at = par("usr")[2],
        text="Expected probability of transition\n(calculated from occurrence probabilities)")
  }

  if(n == 1){
  mtext(side=2, line=1.85, las=0, at = par("usr")[3],
        text="Observed probability / Expected probability")
  }

  # text(x=log(1), y=log(1.05), pos=3, adj=0.5, offset=0.8,
  #      labels="Background\nshifts", col="grey60", font=1)
  # 
  # text(x=log(0.05), y=log(0.75),  adj=0,
  #      labels="Shifts\nto and\nfrom\nnovelty", col="grey60", font=1)
  # 
  # text(x=log(0.0002), y=log(2.8),  adj=0,
  #      labels="Shifts\nbetween\nnovelty\ncategories", col="grey60", font=1)
  
  if(n %in% 3:4){
  axis(side=1, 
       at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
       labels=c(0.0001, 0.001,0.01,0.1,1,10,100))
  } else {
    axis(side=1, 
         at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
         labels=NA)
  }
  
  axis(side=1, at=log(c(seq(0.0001,0.001,0.0001),
                        seq(0.001,0.01,0.001),
                        seq(0.01,0.1,0.01),
                        seq(0.1,1,0.1),
                        seq(1,10,1),
                        seq(10,100,10))), tcl=-0.125, labels=NA)
  
  y.locs <- c(seq(1, 10, 1),
              1/seq(1, 10, 0.5))
  
  axis(side=2, at=log(y.locs), tcl=-0.125, labels=NA)
  
  if(n %in% c(1,3)){
  axis(side = 2, 
       at = log(c(0.5,1,2,3,5,10,20)),
       labels = c(0.5,1,2,3,5,10,20), 
       mgp=c(3,0.5,0))
  } else {
    axis(side = 2, 
         at = log(c(0.5,1,2,3,5,10,20)),
         labels = NA, 
         mgp=c(3,0.5,0))
  }
  
  sub.trans <- sub.trans[order(sub.trans$non.zero),]
  
  sapply(1:dim(sub.trans)[1], function(x){
    print(x)
    
    segments(y0=log(sub.trans$ratio.mean[x]),
             y1=log(sub.trans$ratio.mean[x]),
             x0=log(sub.trans$exp.upper[x]),
             x1=log(sub.trans$exp.lower[x]),
             col=ifelse(sub.trans$non.zero[x],
                        "black", "grey70"))
    
    segments(y0=log(sub.trans$ratio.upper[x]),
             y1=log(sub.trans$ratio.lower[x]),
             x0=log(sub.trans$exp.mean[x]),
             x1=log(sub.trans$exp.mean[x]),
             col=ifelse(sub.trans$non.zero[x],
                        "black", "grey70"))
    
    if(sub.trans$non.zero[x]){
    
      aft.col = c("grey35", cumul.col, "red", "orange")[as.factor(sub.trans$cat.aft)[x]]
      bef.col = c("grey35", cumul.col, "red", "orange")[as.factor(sub.trans$cat.bef)[x]]
      shape = c("triangle", "square", "circle", "diamond")[sub.trans$taxa[x]]
      border.col = "black"
      
    } else {

      aft.col = light.cols[as.factor(sub.trans$cat.aft)[x]]
      bef.col = light.cols[as.factor(sub.trans$cat.bef)[x]]
      shape = c("triangle", "square", "circle", "diamond")[sub.trans$taxa[x]]
      border.col = "grey70"

    }

    arrow.shape(x=log(sub.trans$exp.mean[x]),
                y=log(sub.trans$ratio.mean[x]),
                r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
                col = aft.col,
                shape=shape,
                border=NA)
    
    arrow.shape(x=log(sub.trans$exp.mean[x]),
                y=log(sub.trans$ratio.mean[x]),
                r=0.2, screen=framemat[1,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape=shape,
                col = bef.col, 
                border=border.col, add.arrow=TRUE)
    
    arrow.shape(x=log(sub.trans$exp.mean[x]),
                y=log(sub.trans$ratio.mean[x]),
                r=0.2, screen=framemat[1,], rads=c(0, 2),
                shape=shape, border=border.col, add.arrow=FALSE)

    
  })
  
  text(x=relative.axis.point(0.02, "x"),
       y = relative.axis.point(0.935, "y"),
       labels = paste0("(",LETTERS[n],") ",
                       c("Diatoms", "Foraminifera", "Nannoplankton", "Radiolarians")[n]),
                       font=2, adj=0)
  box()
  close.screen(n)
  })
  
  screen(5)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  
  par(xpd=NA)
  
  taxa.pos <- rev(seq(0.375,0.675,len=4))
  
  text(x=0.575, y=0.75, adj=0.5,
       labels=bquote(bold(underline("Taxa"))), font=2, cex=0.8)
  sapply(taxa.pos, function(y){
    
    arrow.shape(x=0.125,
                y=y,
                r=0.2, screen=framemat[5,], rads=c(1.5,0.5),
                col="white",
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE)
    
    arrow.shape(x=0.125,
                y=y,
                r=0.2, screen=framemat[5,], rads=c(0.5, 1.5),
                col="white",
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE, add.arrow=TRUE, lwd=0.5)
    
    arrow.shape(x=0.125,
                y=y,
                r=0.2, screen=framemat[5,], rads=c(0,2),
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE)
    
    
  })
  # text(x=0.4, pos=4, y=taxa.pos,
  #      labels=c("Nanno", "Foram", "Radio", "Diatom"), 
  #      cex=0.8, offset=0.75)
  
  text(x=0.13, pos=4, y=taxa.pos,
       labels=sort(c("Nannoplankton", "Foraminifera", "Radiolarians", "Diatoms")), 
       cex=0.8, offset=0.75)
  
  par(lheight=0.85)
  text(x=0.15, y=1.05, labels="Preceding\ncommunity", adj=0, cex=0.8,
       lheight=0.5)
  
  text(x=0.85, y=0.925, labels="Succeeding\ncommunity", adj=1, cex=0.8,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.55, y=1.2,
              r=0.25, screen=framemat[5,], rads=c(1.5,0.5),
              col="white", shape="circle",
              border="black", plot=TRUE)
  
  arrow.shape(x=0.45, y=1.2,
              r=0.25, screen=framemat[5,], rads=c(0.5,1.5),
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
                "Cumululative\nnovelty", "Novel\ncommunity"))
  par(lheight=1)
  
  text(x=0.575, y=0.215, adj=0.5,
       labels=bquote(bold(underline("Community"))), font=2, cex=0.8)
  text(x=0.575, y=0.15, adj=0.5,
       labels=bquote(bold(underline("classification"))), font=2, cex=0.8)
  
  par(xpd=FALSE)
  close.screen(5)
  
  dev.off()
}