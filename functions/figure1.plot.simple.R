figure1.plot.simple <- function(prob.model.list,
                                example.ssmat,
                                plot.name,
                                example.ylims){
  
  model.preds <- lapply(c(4:6), function(n){
    temp.pred <- prob.model.list$fixed.prob.models[[n]]$pred.df
    return(temp.pred)
  })
  
  rand.preds <- lapply(c(4:6), function(n){
    temp.pred <- prob.model.list$random.prob.models[[n]]$pred.df
    return(temp.pred)
  })
  
  gradient.ellipse <- function(center, a0, a1, b0, b1, n, alpha){
    
    a.grad <- seq(a0, a1, length.out=n)
    b.grad <- seq(b0, b1, length.out=n)
    
    sapply(1:n, function(n){
      
      draw.ellipse(x=center[1], y=center[2],
                   a=a.grad[n], b=b.grad[n],
                   col=rgb(1,1,1,alpha), border=NA)
      
    })
    
  }
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  if(!is.null(dev.list())){
  dev.off(which=dev.list())
  }
  close.screen(all.screens=TRUE)
  
  pdf(date.wrap(paste0("./plots/probability venn composite (",
                       plot.name, ")"), ".pdf"),
      height=5, width=10)

  split.screen(rbind(c(0.455,1,0.15,0.99), # vens
                     c(0.515,0.660,0.04,0.225), # smallplots
                     c(0.805,0.950,0.04,0.225), # smallplots
                     c(0.660,0.805,0.04,0.225), # smallplots
                     
                     c(0.455,0.99,0.01,0.99), # venn box & label
                     
                     c(0.075,0.25,0.13,0.99), # instant
                     c(0.25,0.425,0.13,0.99), # cumul
                     c(0.075,0.425,0.13,0.99))) # novel box & label
  
  screen(1)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.4,0.6,0.5)
  circle.radius <- c(0.2, 0.175)
  circle.y <- 0.55
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=3)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=3)
  
  par(xpd=NA)
  segments(x0=circle.cent + c(0, 0, 0.0225), 
           y0=rep(circle.y, 3) + c(0, 0, 0.4),
           x1=circle.cent + c(-0.2,0.2, 0.0225),
           y1=c(0, 0, 0),
           lwd=1)
  
  par(xpd=FALSE)
  rect(xleft=0.4, xright=0.6, ybottom=0.85, ytop=0.95, border=NA, col="white")
  
  # turn I circle to spatial polygons to get overlap
  I.sp <- Polygon(cbind(I$x, I$y))
  I.sp <- SpatialPolygons(list(Polygons(list(I.sp), ID = "a")))
  
  C.sp <- Polygon(cbind(C$x, C$y))
  C.sp <- SpatialPolygons(list(Polygons(list(C.sp), ID = "a")))
  
  N.sp <- raster::intersect(I.sp, C.sp)
  I.sub <- gDifference(I.sp, C.sp)
  C.sub <- gDifference(C.sp, I.sp)
  
  plot(I.sub, col="red", add=TRUE, border="black", lwd=1)
  plot(C.sub, col=cumul.col, add=TRUE, border="black", lwd=1)
  plot(N.sp, col="orange", add=TRUE, border="black", lwd=2)
  
  text(x=0.25, y=0.75,
       labels=c("Instantaneous\nnovelty (I)   "), 
       col="red", adj=1, cex=1.15)
  text(x=0.735, y=0.75,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=1.15)
  
  text(x=0.5225, y=0.915, adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=1.15)
 
  overall.means <- t(sapply(rand.preds, function(x){
    
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"u-95% CI"])
    lower <- plogis(x[1,"l-95% CI"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  text(x=c(circle.cent + c(-0.075, 0.085, 0.01)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1.1)
  
    text(x = c(circle.cent + c(-0.075, 0.085, 0.01)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), " - ",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=0.9)
  close.screen(1)
  
  sapply(2:4, function(n){
    screen(n)
    
    par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0), las=1)
    
    point.col <- c("red", cumul.col, "orange")[n-1]
    
    pred.df <- as.data.frame(model.preds[[n-1]])
    
    ylims <- c(0, ifelse(n %in% 5:6, 0.09, 0.06))
    
    text.pos <- rep(c(-0.005, 0.005), 2) * ylims[2]/0.055
    
    plot(x=NULL, y=NULL, xlab="", ylab="",
         xlim=c(0.25,5), xaxs="i", yaxs="i",
         ylim=ylims, axes=FALSE)
    
    if(n == 2){
      axis(side=2, at=c(0,0.02,0.04,0.06,0.08,0.1), cex.axis=0.8)
      axis(side=2, at=seq(0,0.1,0.01), labels=NA, tcl=-0.125)
      mtext(side=2, line=1.75, text="Probability", cex=0.8, las=0)
    } else {
      axis(side=2, at=c(0,0.02,0.04,0.06,0.08), labels=NA)
      axis(side=2, at=seq(0,0.1,0.01), labels=NA, tcl=-0.125)
    }
    
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=par("usr")[3], ytop=par("usr")[4], 
         col="white", border=NA)
    
    pred.df$lower <- plogis(pred.df[,"l-95% CI"])
    pred.df$upper <- plogis(pred.df[,"u-95% CI"])
    
    segments(x0=1:4, x1=1:4,
             y0=pred.df$lower,
             y1=pred.df$upper)
    
    points(y=plogis(pred.df$Estimate), x=1:4, bg=point.col,
           pch=21)
    
    text(x=1:4, 
         y=c(pred.df$lower[1], pred.df$upper[2], pred.df$lower[3], pred.df$upper[4]) + text.pos,
         labels=c("Diatom", "Foram", "Nanno", "Radio"), 
         adj=0.5, cex=0.8)
    
    box()
    close.screen(n)
  })
  
  screen(5)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot.new()
  text(x=relative.axis.point(0.03, "x"),
       y=relative.axis.point(0.97, "y"),
       labels="(B)", font=2)
  box()
  close.screen(5)
  
  screen(6)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  
  return.data <- identify.novel.gam(site.sp.mat = example.ssmat,
                                    alpha=0.05,
                                    metric="jaccard",
                                    site="1338",
                                    plot=FALSE,
                                    plot.data=TRUE)
  
  min.p <- return.data[[3]][-(1:10),]
  seq.p <- return.data[[2]][-(1:10),]
  save.data <- return.data
  return.data<-return.data[[1]][-(1:10),]
  
  plot(as.numeric(as.character(return.data$bins)) ~ return.data$seq.dist,
       type="n", xlim=c(0,1), yaxs="i", ylim=example.ylims,
       axes=FALSE, xlab="", ylab="", yaxt="n")
  axis(side=1, lwd=0.5, at=seq(0,0.8,0.2), mgp=c(3,0.2,0))
  lims <- par("usr")
  
  axis(side=2, las=1)
  mtext(side=2, line=2, text="Millions of years BP")
  axis(side=2, at=seq(1,65,1), tcl=-0.125, labels=NA)
  
  polygon(y=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          x=c(seq.p$lwr, rev(seq.p$upr)),
          col="grey75", border=NA)
  
  lines(x=seq.p[,4], y=as.numeric(as.character(return.data$bins)), col="grey15",
        lty="dashed")
  
  with(return.data,
       lines(x=seq.dist, y=as.numeric(as.character(bins)), lwd=0.75))
  
  with(return.data[return.data$instant & !is.na(return.data$instant),],
       points(x=seq.dist, y=as.numeric(as.character(bins)),
              pch=21, bg="red"))
  
  sub.data <- return.data[return.data$instant & !is.na(return.data$instant),]
  
  if(sum(sub.data$instant) != sum(sub.data$novel)){
  par(xpd=NA)
  with(sub.data[sub.data$cat != "novel",],
       text(x=seq.dist, y=as.numeric(as.character(bins)),
            labels="(I)", pos=4,
            col="red"))
  par(xpd=FALSE)
  }
  
  sapply(which(return.data$novel), function(x){
    segments(y0 = as.numeric(as.character(return.data$bins))[x],
             y1 = as.numeric(as.character(return.data$bins))[x],
             x1 = return.data$seq.dist[x] - (0.05 * (par("usr")[1] - par("usr")[2])),
             x0 = par("usr")[2], col="orange", lwd=2)
  })
  
  mtext(side=1, text = "Instantaneous\ndissimilarity", line=2)
  
  text(x=relative.axis.point(0.10, "x"),
       y=relative.axis.point(0.965, "y"),
       labels="(A)", font=2)
  
  close.screen(6)
  
  screen(7)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  
  plot(x=return.data$raw.min.dist,
       y=as.numeric(as.character(return.data$bins)), type="n",
       xlim=rev(c(0, 1)),
       ylim=c(lims[3], lims[4]),
       yaxs="i",
       axes=FALSE, ylab="", xlab="")
  
  polygon(y=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          x=c(min.p$lwr, rev(min.p$upr)),
          col="grey75", border=NA)
  
  lines(x=min.p[,4], y=as.numeric(as.character(return.data$bins)), col="grey15",
        lty="dashed")
  
  with(return.data,
       lines(x=raw.min.dist, y=as.numeric(as.character(bins)), lwd=1))
  
  with(return.data[return.data$cumul,],
       points(x=raw.min.dist, y=as.numeric(as.character(bins)),
              pch=21, bg=cumul.col))
  
  sub.data <- return.data[return.data$cumul,]
  
  par(xpd=NA)
  if(sum(sub.data$cumul) != sum(sub.data$novel)){
  with(sub.data[sub.data$cat != "novel",],
       text(x=seq.dist, y=as.numeric(as.character(bins)),
            labels="(C)", pos=2,
            col=cumul.col))
  }
  
  sapply(which(return.data$novel), function(x){
    segments(y0 = as.numeric(as.character(return.data$bins))[x],
             y1 = as.numeric(as.character(return.data$bins))[x],
             x0 = return.data$raw.min.dist[x] - (0.05 * (par("usr")[2] - par("usr")[1])),
             x1 =par("usr")[1], col="orange", lwd=2)
    
    points(y=as.numeric(as.character(return.data$bins))[x],
           x=par("usr")[1], pch=21, bg="orange")
    
    text(y=as.numeric(as.character(return.data$bins))[x],
         x=par("usr")[1], labels="(N)", pos=3, col="orange")
  })
  par(xpd=FALSE)
  
  axis(side=1, at=seq(0,0.8,0.2), mgp=c(3,0.2,0))
  axis(side=4, labels=NA)
  axis(side=4, at=seq(1,65,1), tcl=-0.125, labels=NA)
  
  mtext(side=1, text = "Cumulative\ndissimilarity", line=2)
  close.screen(7)
  
  screen(8)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot.new()
  box()
  close.screen(8)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}