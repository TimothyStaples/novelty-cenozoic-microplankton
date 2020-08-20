#prob.model.list <- novel.prob.list

figure1.plot <- function(prob.model.list,
                         example.ssmat,
                         plot.name){
  
  model.preds <- lapply(c(6,4,5,1,2), function(n){
    temp.pred <- prob.model.list$fixed.prob.models[[n]]$pred.df
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
  
  library(plotrix)
  library(sp)
  library(raster)
  library(rgeos)
  library(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  if(!is.null(dev.list())){
  dev.off(which=dev.list())
  }
  close.screen(all.screens=TRUE)
  
  pdf(date.wrap(paste0("./plots/probability venn composite (",
                       plot.name, 
                       ")"),
                ".pdf"),
      height=7, width=10)
  
  split.screen(rbind(c(0.5,1,0.15,0.85), # vens
                     c(0.675,0.825,0.04,0.17), # smallplots
                     c(0.525,0.675,0.04,0.17), # smallplots
                     c(0.825,0.975,0.04,0.17), # smallplots
                     c(0.60,0.75,0.84,0.97), # smallplots
                     c(0.75,0.9,0.84,0.97), # smallplots
                     
                     c(0.455,0.99,0.01,0.99), # venn box & label
                     
                     c(0.075,0.25,0.1,0.99), # instant
                     c(0.25,0.425,0.1,0.99), # cumul
                     c(0.075,0.425,0.1,0.99))) # novel box & label
  
  screen(1)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="")
  
  # top circle
  circle.cent<-c(0.42,0.58,0.5)
  circle.radius <- 0.175
  circle.y <- 0.75
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius, nv=500, 
                 border="black", col="white", lwd=3)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius, nv=500, 
                 border="black", col="white", lwd=3)
  
  par(xpd=NA)
  segments(x0=circle.cent + c(-0.15,0.15,0), 
           y0=rep(circle.y, 3)-0.45,
           x1=circle.cent + c(-0.25,0.25,0),
           y1=c(-0.2, -0.2, -0.2),
           lwd=1)
  
  segments(x0=circle.cent[1:2], 
           y0=rep(circle.y, 2),
           x1=circle.cent[1:2] + c(-0.1,0.1),
           y1=c(1.1, 1.1),
           lwd=1)
  
  par(xpd=FALSE)
  rect(xleft=0.4, xright=0.6, ybottom=0.01, ytop=0.085, border=NA, col="white")
  
  # turn I circle to spatial polygons to get overlap
  I.sp <- Polygon(cbind(I$x, I$y))
  I.sp <- SpatialPolygons(list(Polygons(list(I.sp), ID = "a")))
  
  C.sp <- Polygon(cbind(C$x, C$y))
  C.sp <- SpatialPolygons(list(Polygons(list(C.sp), ID = "a")))
  
  N.sp <- raster::intersect(I.sp, C.sp)
  I.sub <- gDifference(I.sp, C.sp)
  C.sub <- gDifference(C.sp, I.sp)
  
  Arrows(x0=circle.cent + c(-0.075,0.075,0), 
         y0=c(circle.y, circle.y, 0.585),
         x1=circle.cent + c(-0.09,0.09,0),
         y1=c(0.4825, 0.4825, 0.4675),
         lwd=1, arr.type="triangle",
         arr.length=0.15, arr.width=0.15)
  
  overlap.col <- rgb(174/255,146/255,161/255)
  plot(N.sp, col="white", add=TRUE, border="white", lwd=1.5)
  
  plot(I.sub, col="red", add=TRUE, border="red", lwd=1)
  plot(C.sub, col=cumul.col, add=TRUE, border=cumul.col, lwd=1)
  
  plot(N.sp, col=overlap.col, add=TRUE, border=overlap.col, lwd=1)
  
  text(x=0.23, y=circle.y,
       labels=c("Large\nshift in \ncommunity \ncomposition"), 
       col="red", adj=1, cex=1.15)
  text(x=0.765, y=circle.y,
       labels=c("Shift to\n a new\n community\nstate"), 
       col=cumul.col, adj=0, cex=1.15)
  
  overall.means <- sapply(model.preds[4:5], function(x){
    upper <- plogis(x$fit + 1.96*x$se.fit)
    lower <- plogis(x$fit - 1.96*x$se.fit)
    round(range(c(lower, upper))*100, 1)
  })
  
  text(x=c(circle.cent[1]-0.025,
           circle.cent[2]+0.025), 
       y= rep(circle.y, 2),
       labels=paste0(sprintf(overall.means[1,], fmt = '%#.1f'), "-",
                     sprintf(overall.means[2,], fmt = '%#.1f'),"%"),
       font=2, cex=1.15)
  
  I.sub <- elide(I.sub, shift=c(-0.05, -0.45))
  C.sub <- elide(C.sub, shift=c(0.05, -0.45))
  N.sub <- elide(N.sp, shift=c(0, -0.45))
  
  plot(N.sub, col="orange", add=TRUE, border="black", lwd=1)
  plot(I.sub, col="red", add=TRUE, border="black", lwd=1)
  plot(C.sub, col=cumul.col, add=TRUE, border="black", lwd=1)
  
  text(x=0.22, y=0.14,
       labels=c("Instantaneous\nnovelty (I)"), 
       col="red", adj=1, cex=1.15)
  text(x=0.78, y=0.14,
       labels=c("Cumulative\nnovelty (C)"), 
       col=cumul.col, adj=0, cex=1.15)
  
  text(x=0.5, y=0.05, adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=1.15)
  
  overall.means <- sapply(model.preds[c(2,1,3)], function(x){
    upper <- plogis(x$fit + 1.96*x$se.fit)
    lower <- plogis(x$fit - 1.96*x$se.fit)
    
    round(range(c(lower, upper))*100, 1)
  })
  
  text(x=c(circle.cent[1]-0.145,
           circle.cent[3],
           circle.cent[2]+0.145), y= rep(circle.y, 3)-0.45,
       labels=paste0(sprintf(overall.means[1,], fmt = '%#.1f'), "-",
                     sprintf(overall.means[2,], fmt = '%#.1f'),"%"), 
       font=2, cex=1.1)
  close.screen(1)
  
  sapply(2:6, function(n){
    screen(n)
    
    par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0), las=1)
    
    point.col <- c("orange","red",cumul.col, "red", cumul.col)[n-1]
    
    pred.df <- model.preds[[n-1]]
    
    ylims <- c(0, ifelse(n %in% 5:6, 0.09, 0.06))
    
    text.pos <- rep(c(-0.005, 0.005), 2) * ylims[2]/0.055
    
    plot(x=NULL, y=NULL, xlab="", ylab="",
         xlim=c(0.25,5), xaxs="i", yaxs="i",
         ylim=ylims, axes=FALSE)
    
    if(n %in% c(3,5)){
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
    
    pred.df$lower <- plogis(pred.df$fit - 1.96*pred.df$se.fit)
    pred.df$upper <- plogis(pred.df$fit + 1.96*pred.df$se.fit)
    
    segments(x0=1:4, x1=1:4,
             y0=pred.df$lower,
             y1=pred.df$upper)
    
    points(y=plogis(pred.df$fit), x=1:4, bg=point.col,
           pch=21)
    
    text(x=1:4, 
         y=c(pred.df$lower[1], pred.df$upper[2], pred.df$lower[3], pred.df$upper[4]) + text.pos,
         labels=c("Nanno", "Foram", "Radio", "Diatom"), 
         adj=0.5, cex=0.8)
    
    box()
    close.screen(n)
  })
  
  screen(7)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot.new()
  text(x=relative.axis.point(0.03, "x"),
       y=relative.axis.point(0.98, "y"),
       labels="(B)", font=2)
  box()
  close.screen(7)
  
  screen(8)
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
  
  ylims <- c(max(c(max(seq.p$upr, na.rm=TRUE), max(return.data$seq.dist, na.rm=TRUE))) * 1.05,
             min(c(min(seq.p$lwr, na.rm=TRUE), min(return.data$seq.dist, na.rm=TRUE))) *0.95)
  
  plot(as.numeric(as.character(return.data$bins)) ~ return.data$seq.dist,
       type="n", xlim=c(0.1,0.9), yaxs="i", ylim=rev(c(-0.25,15.5)),
       axes=FALSE, xlab="", ylab="", yaxt="n")
  axis(side=1, lwd=0.5)
  lims <- par("usr")
  
  axis(side=2, las=1)
  mtext(side=2, line=2, text="Millions of years BP")
  axis(side=2, at=seq(1,20,1), tcl=-0.125, labels=NA)
  
  polygon(y=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          x=c(seq.p$lwr, rev(seq.p$upr)),
          col="grey75", border=NA)
  
  lines(x=seq.p[,4], y=as.numeric(as.character(return.data$bins)), col="grey15",
        lty="dashed")
  
  with(return.data,
       lines(x=seq.dist, y=as.numeric(as.character(bins)), lwd=1))
  
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
  
  mtext(side=1, text = "Instantaneous\ndissimilarity", line=2.25)
  
  text(x=relative.axis.point(0.10, "x"),
       y=relative.axis.point(0.98, "y"),
       labels="(A)", font=2)
  
  close.screen(8)
  
  screen(9)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  
  ylims <- c(0.1,0.9)
  
  plot(x=return.data$raw.min.dist,
       y=as.numeric(as.character(return.data$bins)), type="n",
       xlim=rev(ylims),
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
  
  axis(side=1)
  axis(side=4, labels=NA)
  axis(side=4, at=seq(1,20,1), tcl=-0.125, labels=NA)
  
  mtext(side=1, text = "Cumulative\ndissimilarity", line=2.25)
  close.screen(9)
  
  screen(10)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot.new()
  box()
  close.screen(10)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}