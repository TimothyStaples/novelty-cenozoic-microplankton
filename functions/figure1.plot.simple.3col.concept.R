figure1.plot.concept <- function(example.ssmat,
                                plot.name,
                                example.xlims,
                                novel.lab.pos = list(3, 1, 4),
                                ts.cut=0){
  
  return.data <- identify.novel.gam(site.sp.mat = example.ssmat,
                                    alpha=0.05,
                                    metric="jaccard",
                                    site="1338",
                                    plot=FALSE,
                                    plot.data=TRUE)
  
  min.p <- return.data[[3]][-(0:ts.cut),]
  seq.p <- return.data[[2]][-(0:ts.cut),]
  save.data <- return.data
  return.data<-return.data[[1]][-(0:ts.cut),]
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  cumul.col.rgb <- col2rgb(cumul.col)/255
  
  if(!is.null(dev.list())){
  dev.off(which=dev.list())
    close.screen(all.screens=TRUE)
  }
  
  instant.point.concept <- function(example.points,
                                    nov.points, 
                                    nov.rad, axis1, axis2){
    
    plot(x=NULL, y=NULL, xlim=c(0.2,0.825), ylim=c(0.15,1), axes=FALSE,
         ylab="", xlab="")
    box()
    
    if(axis1){mtext(side=1, line=-0.3, text="nMDS 1")}
    if(axis2){mtext(side=2, line=0, text="nMDS 2", las=0)}
    
    sapply(1:length(nov.points), function(n){
    draw.circle(x=example.points[nov.points[n],1],
                y=example.points[nov.points[n],2],
                radius=nov.rad[n],
                border=NA, col=rgb(0.5,0.5,0.5,0.5))
    })
    
    lines(example.points, col="black", lwd=0.5)
    
    sapply(2:nrow(example.points), function(point){
    diff <- example.points[point,] - example.points[point-1,]
    arrows(x0=example.points[point,1],
           y0=example.points[point,2],
           x1=example.points[point,1] - (0.4*diff[1]),
           y1=example.points[point,2] - (0.4*diff[2]),
           lwd=0.5, length=0.03, angle=150)
    })
    
    sapply(1:length(nov.points), function(n){
      
      nov.point <- nov.points[n]
      
      lines(example.points[nov.point + c(-1,0),], col="red", lwd=1)
      
      points(x=example.points[nov.point-1,1], 
             y=example.points[nov.point-1,2], pch=21, cex=0.75, bg="red", lwd=0.5)
      
      # draw on axes and text to act as labels
      par(xpd=NA)
      axis(side=1, at=c(relative.axis.point(c(0.2,0.6)[n], "x"), 
                        relative.axis.point(c(0.2,0.6)[n], "x") + nov.rad[n]), tcl=0.1, labels=NA, col.ticks="red",
           line=-3, col="red", lwd=0.5)
      
      axis(side=1, at=c(relative.axis.point(0.2, "x"), 
                        relative.axis.point(0.2, "x") + dist(example.points[nov.point + c(-1,0),])), 
           tcl=0.1, labels=NA, col.ticks="red",
           line=-4, col="red", lwd=0.5)
      par(xpd=FALSE)
      
      text(x=example.points[nov.point + -1,1], y = example.points[nov.point + -1,2],
           labels=expression('I'["obs"]), col="red", cex=0.8) 
      text(x=example.points[nov.point,1], y = example.points[nov.point,2],
           labels=expression('I'["crit"]), col="red", cex=0.8) 
      
    })
    
    points(example.points[-c(nov.points, nov.points - 1),], pch=16, cex=0.75, col="black")
    points(example.points[nov.points,], pch=21, cex=0.75, col="black", bg="white",lwd=0.5)
    
    text(x=relative.axis.point(0.5, "x"),
         y=relative.axis.point(0.075, "y"), adj=0.5,
         labels=expression("I"["obs"]*" > I"["crit"]*" = Instantaneous novelty"),
         col="red", cex=0.8)
    
    box()
  }
  
  cumul.point.concept <- function(example.points,
                               nov.points, nov.rad, axis1, axis2){
  
    plot(x=NULL, y=NULL, xlim=c(0.2,0.825), ylim=c(0.15,1), axes=FALSE,
         ylab="", xlab="")
  box()
  
  if(axis1){mtext(side=1, line=0, text="nMDS 1")}
  if(axis2){mtext(side=2, line=0, text="nMDS 2", las=0)}
  
  sapply(1:length(nov.points), function(n){
    draw.circle(x=example.points[nov.points[n],1],
                y=example.points[nov.points[n],2],
                radius=nov.rad[n],
                border=NA, col=rgb(0.5,0.5,0.5,0.5))
  })
  
  lines(example.points, col="black", lwd=0.5)
  
  sapply(2:nrow(example.points), function(point){
    diff <- example.points[point,] - example.points[point-1,]
    arrows(x0=example.points[point,1],
           y0=example.points[point,2],
           x1=example.points[point,1] - (0.4*diff[1]),
           y1=example.points[point,2] - (0.4*diff[2]),
           lwd=0.5, length=0.03, angle=150)
  })
  
  min.novel <- sapply(1:length(nov.points), function(n){
    
  nov.point <- nov.points[n]
    
  novel.dist <- as.matrix(dist(example.points))
  min.novel <- which.min(novel.dist[-(nov.point:nrow(example.points)),nov.point])
  
  segments(x0=example.points[nov.point,1],
           y0=example.points[nov.point,2],
           x1=example.points[min.novel,1],
           y1=example.points[min.novel,2],
           col=cumul.col, lwd=1, lty="31")

  par(xpd=NA)
  axis(side=1, at=c(relative.axis.point(c(0.2,0.6)[n], "x"), 
                    relative.axis.point(c(0.2,0.6)[n], "x") + nov.rad[n]), tcl=0.1, labels=NA, col.ticks=cumul.col,
       line=-3, col=cumul.col, lwd=0.5)
  
  axis(side=1, at=c(relative.axis.point(0.2, "x"), 
                    relative.axis.point(0.2, "x") + dist(example.points[nov.point + c(-1,0),])), 
       tcl=0.1, labels=NA, col.ticks=cumul.col,
       line=-4, col=cumul.col, lwd=0.5)
  par(xpd=FALSE)
  
  text(x=example.points[nov.point + -1,1], y = example.points[nov.point + -1,2],
       labels=expression('C'["obs"]), col=cumul.col, cex=0.8) 
  text(x=example.points[nov.point,1], y = example.points[nov.point,2],
       labels=expression('C'["crit"]), col=cumul.col, cex=0.8) 
  
  return(min.novel)
  })
  
  points(x=example.points[min.novel,1], y=example.points[min.novel,2], 
         pch=21, cex=0.75, bg=cumul.col, lwd=0.5)

  points(example.points[-c(nov.points, min.novel),], pch=16, cex=0.75, col="black")
  points(example.points[nov.points,], pch=21, cex=0.75, col="black", bg="white",lwd=0.5)
  
  
  text(x=relative.axis.point(0.5, "x"),
       y=relative.axis.point(0.075, "y"), adj=0.5,
       labels=expression("C"["obs"]*" > C"["crit"]*" = Cumulative novelty"),
       col=cumul.col, cex=0.8)
  
  box()
  }

  pdf(date.wrap(paste0("./plots/probability venn composite (",
                       plot.name, ")"), ".pdf"),
      height=4.75, width=2.25, useDingbats = FALSE,
      colormodel = "cmyk")

  split.screen(rbind(c(0.2,0.99,0.03,0.2225),
                     c(0.2,0.99,0.2225,0.415),
                     
                     c(0.2,0.99,0.475,0.6), # instant
                     c(0.2,0.99,0.6,0.725), # cumul
                     c(0.2,0.99,0.475,0.725),  # novel box & label
                     c(0.2,0.99,0.725,0.99))) # time-series
  
  screen(1)
  par(mar=c(0,0,0,0), ps=6, tcl = -0.2, mgp=c(3,0.5,0), las=1)
  
  # example points
  example.points = rbind(c(0.5,1.15),c(0.425,0.8),c(0.47,0.7),c(0.49,0.58),
                         c(0.48,0.52), c(0.55,0.5), c(0.6,0.6),
                         c(0.65,0.65), c(0.65,0.5), c(0.7,0.55),
                         c(0.675,0.7), c(0.7,0.825), c(0.75,0.75),
                         c(0.78, 0.675), c(0.75,0.55), c(0.75, 0.375),
                         c(0.7, 0.28), c(0.675, 0.36),c(0.6,0.325),
                         c(0.5,0.35), c(0.35,0.32), c(0.3,0.78), 
                         c(0.25, 0.87), c(0.25, 0.8), c(0.15,0.75))

  instant.point.concept(example.points = example.points,
                        nov.points = nrow(example.points)- c(3, 10),
                        nov.rad = c(0.08, 0.08),
                        axis1=TRUE, axis2=TRUE)
  text(x=relative.axis.point(0.015, "x"),
       y=relative.axis.point(0.925, "y"),
       labels="(D)", font=2, adj=0)
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=6, tcl = -0.2, mgp=c(3,0.5,0), las=1)
  
  cumul.point.concept(example.points = example.points,
                        nov.point = nrow(example.points)- c(3, 10),
                        nov.rad = c(0.08, 0.08),
                        axis1=FALSE, axis2=TRUE)
  text(x=relative.axis.point(0.015, "x"),
       y=relative.axis.point(0.925, "y"),
       labels="(C)", font=2, adj=0)
  close.screen(2)
  
  screen(3)
  par(mar=c(0,0,0,0), ps=6, tcl = -0.2, mgp=c(3,0.5,0), las=1)
  
  plot(return.data$seq.dist ~ as.numeric(as.character(return.data$bins)),
       type="n", ylim=c(0,1), xlim=example.xlims,
       axes=FALSE, xlab="", ylab="", yaxt="n")
  axis(side=2, lwd=0.5, at=seq(0,0.8,0.2), mgp=c(3,0.3,0))
  lims <- par("usr")
  
  axis(side=1, las=1, mgp=c(3,-0.2,0))
  mtext(side=1, line=0.3, text="Millions of years BP")
  axis(side=1, at=seq(1,65,1), tcl=-0.1, labels=NA)
  
  with(return.data,
       lines(y=seq.dist, x=as.numeric(as.character(bins)),
             col=rgb(1,0,0,0.65)))
  
  polygon(x=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          y=c(seq.p$lwr, rev(seq.p$upr)),
          col=rgb(0.5,0.5,0.5,0.5), border=NA)
  
  sub.data <- return.data[return.data$instant & !is.na(return.data$instant),]
  
  if(sum(sub.data$instant) != sum(sub.data$novel)){
    par(xpd=NA)
    with(sub.data[sub.data$cat != "novel",],
         text(y=seq.dist, x=as.numeric(as.character(bins)),
              labels="(I)", pos=novel.lab.pos,
              col="red"))
    par(xpd=FALSE)
  }
  
  par(xpd=NA)
  legend(x=40, y=relative.axis.point(1.25, "y"),
         pch=16, pt.cex=0.35, legend="= B", bty="n",
         x.intersp=0.2, col="grey40", adj=c(0,0.3))
  par(xpd=FALSE)
  
  sapply(which(return.data$novel), function(x){
    segments(x0 = as.numeric(as.character(return.data$bins))[x],
             x1 = as.numeric(as.character(return.data$bins))[x],
             y1 = return.data$seq.dist[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
             y0 = par("usr")[4], col="orange", lwd=1)
  })
  
  with(return.data[return.data$instant & !is.na(return.data$instant),],
       points(y=seq.dist, x=as.numeric(as.character(bins)),
              pch=21, bg="red", cex=0.65, lwd=0.75))
  
  with(return.data[return.data$cat == "back" & !is.na(return.data$cat),],
       points(y=seq.dist, x=as.numeric(as.character(bins)),
              pch=16, col="grey40", cex=0.35, lwd=0.5))
  
  
  lines(y=seq.p[,4], x=as.numeric(as.character(return.data$bins)), col="black",
        lty="31", lwd=1)
  
  mtext(side=2, text = "Instantaneous\ndissimilarity", line=1, las=0)
  close.screen(3)
  
  screen(4)
  par(mar=c(0,0,0,0), ps=6, tcl = -0.2, mgp=c(3,0.5,0), las=1)
  
  plot(y=return.data$raw.min.dist,
       x=as.numeric(as.character(return.data$bins)), type="n",
       ylim=rev(c(0, 0.9)),
       xlim=c(lims[1], lims[2]),
       xaxs="i",
       axes=FALSE, ylab="", xlab="")
  
  cumul.col.rgb <- col2rgb(cumul.col)/255
  with(return.data,
       lines(y=raw.min.dist, x=as.numeric(as.character(bins)), lwd=1,
       col=rgb(cumul.col.rgb[1], cumul.col.rgb[2], cumul.col.rgb[3], 0.65)))
  
  polygon(x=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          y=c(min.p$lwr, rev(min.p$upr)),
          col=rgb(0.5,0.5,0.5,0.5), border=NA)
  
  par(xpd=NA)
  sapply(1:sum(return.data$novel), function(n){
    print(n)
    x <- which(return.data$novel)[n]
    
    segments(x0 = as.numeric(as.character(return.data$bins))[x],
             x1 = as.numeric(as.character(return.data$bins))[x],
             y0 = return.data$raw.min.dist[x] - (0.05 * (par("usr")[4] - par("usr")[3])),
             y1 = par("usr")[3], col="orange", lwd=1)
    
    points(x=as.numeric(as.character(return.data$bins))[x],
           y=par("usr")[3], pch=21, bg="orange", cex=0.65, lwd=0.75)
    
    text(x=as.numeric(as.character(return.data$bins))[x],
         y=par("usr")[3], labels="(N)", pos=novel.lab.pos, 
         col="orange", offset=0.25)
  })
  
  par(xpd=FALSE)
  
  with(return.data[return.data$cat=="back",],
       points(y=raw.min.dist, x=as.numeric(as.character(bins)),
              pch=16, col="grey40", cex=0.35, lwd=0.5))
  
  lines(y=min.p[,4], x=as.numeric(as.character(return.data$bins)), col="black",
        lty="31", lwd=1)
  
  with(return.data[return.data$cumul,],
       points(y=raw.min.dist, x=as.numeric(as.character(bins)),
              pch=21, bg=cumul.col, cex=0.65, lwd=0.75))
  

  
  
  sub.data <- return.data[return.data$cumul,]
  
  if(sum(sub.data$cumul) != sum(sub.data$novel)){
  with(sub.data[sub.data$cat != "novel",],
       text(y=raw.min.dist, x=as.numeric(as.character(bins)),
            labels="(C)", pos=novel.lab.pos,
            col=cumul.col))
  }

  axis(side=2, at=seq(0,0.8,0.2), mgp=c(3,0.3,0))

  mtext(side=2, text = "Cumululative\ndissimilarity", line=1, las=0)
  
  text(x=relative.axis.point(0.0125, "x"),
       y=relative.axis.point(0.89, "y"),
       labels="(B)", font=2, adj=0)
  
  close.screen(4)
  
  screen(5)
  par(mar=c(0,0,0,0), ps=6, tcl = -0.2, mgp=c(3,0.5,0))
  plot.new()
  box()
  close.screen(5)
  
  screen(6)
  par(mar=c(0,0,0,0), ps=6, tcl = -0.2, mgp=c(3,0.5,0))
  plot.new()
  
  circle.y <- 0.275
  circle.n <- 6
  circle.cents <- relative.axis.point(seq(0.055,0.755,len=circle.n), "x")
  arrow.offset <- seq(0.075, 0.6, len=circle.n)

  segments(x0=min(circle.cents), x1 = max(circle.cents), y0 = circle.y, y1 = circle.y, lwd=1)
  segments(y0=circle.y, y1 = circle.y, x0 = max(circle.cents), x1 = max(circle.cents)+0.1, lwd=1)
  Arrows(y0=circle.y, y1 = circle.y, x0 = max(circle.cents)+0.07, x1 = max(circle.cents)+0.08,
         arr.type="triangle", arr.width=0.125, arr.length=0.125)
  
  par(lheight=0.75)
  text(y=circle.y-0.015, x=max(circle.cents)+0.055,
       pos=4, labels="Time\nseries", cex=0.8)
  par(lheight=1)
  
  sub.dists <- return.data[as.numeric(as.character(return.data$bins)) < 30,]
    
  sub.ssmat <- example.ssmat[as.numeric(rownames(example.ssmat)) < 30, ][circle.n:1,]
  sub.dist <- as.matrix(vegdist(sub.ssmat, method="jaccard"))
  
  circle.diffs <- sapply(1:circle.n, function(n){
  
      x = relative.axis.point(seq(0.1,0.9,len=circle.n), "x")[n]
    
      return(rev(circle.cents)[1]-x)
           })

  segments(y0 = circle.y,
           y1 = circle.y - arrow.offset[2],
           x0 = circle.cents[-1],
           x1 = circle.cents[-1], col="red", lwd=1)
  
  segments(y0 = circle.y  - arrow.offset[2],
           y1 = circle.y - arrow.offset[2],
           x0 = circle.cents[-1],
           x1 = circle.cents[-1] - diff(circle.cents[1:2]) + 0.025, col="red", lwd=1)
  
  Arrows(y0 = circle.y  - arrow.offset[2],
         y1 = circle.y  - arrow.offset[2] + 0.05,
         x0 = circle.cents[-1] - diff(circle.cents[1:2]) + 0.025,
         x1 = circle.cents[-1] - diff(circle.cents[1:2]) + 0.025, 
         col="red", lwd=1,
         arr.type = "triangle", arr.length=0.1, arr.width=0.1)
  
  rect(xleft=rowMeans(cbind(circle.cents[-1], circle.cents[-length(circle.cents)]))-0.035,
       xright=rowMeans(cbind(circle.cents[-1], circle.cents[-length(circle.cents)]))+0.06,
       ybottom=circle.y  - arrow.offset[2] - 0.05, 
       ytop=circle.y  - arrow.offset[2] + 0.05,
       border=NA, col=rgb(1,1,1,1))
  
  labels=sprintf("%.2f",diag(sub.dist[-1,-ncol(sub.dist)]))
  sapply(1:length(labels), function(n){
    temp.lab <- labels[n]
    
    text(y = (circle.y  - arrow.offset[2]),
       x = (rowMeans(cbind(circle.cents[-1], circle.cents[-length(circle.cents)])))[n] +0.01,
       labels = bquote(underline(.(temp.lab))), col="red", 
       adj=c(0.5,0.5), cex=0.8)
  })
  
  par(xpd=NA)
  text(y=circle.y  - arrow.offset[2] + 0.055,
       x=mean(circle.cents), col="red", pos=1, offset=0.75,
       labels='Instantaneous dissimilarity', cex=0.8)
  
  text(x=relative.axis.point(0.015, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(A)", font=2, adj=0)
  par(xpd=FALSE)

sapply(2:length(circle.cents), function(n){
    
  segments(y0 = circle.y,
           y1 = circle.y + arrow.offset[n],
           x0 = circle.cents[n],
           x1 = circle.cents[n], col=cumul.col, lwd=1) 
    
  segments(y0 = circle.y + arrow.offset[n],
           y1 = circle.y + arrow.offset[n],
           x0 = circle.cents[1],
           x1 = circle.cents[n], col=cumul.col, lwd=1)
  
  Arrows(y0 = circle.y + arrow.offset[n],
         y1 = circle.y + arrow.offset[n] - 0.05,
         x0 = circle.cents[-1][1:(n-1)] - diff(circle.cents[1:2]),
         x1 = circle.cents[-1][1:(n-1)] - diff(circle.cents[1:2]), 
         col=cumul.col, lwd=1,
         arr.type = "triangle", arr.length=0.1, arr.width=0.1)
  
  rect(xleft = circle.cents[-1][1:(n-1)] - 0.5*diff(circle.cents[1:2]) - 0.05,
       xright = circle.cents[-1][1:(n-1)] - 0.5*diff(circle.cents[1:2]) + 0.05,
       ybottom = circle.y + arrow.offset[n] - 0.05, 
       ytop = circle.y + arrow.offset[n] + 0.05,
       border=NA, col=rgb(1,1,1,1))
  
  text.size = sub.dist[1:(n-1),n] == min(sub.dist[1:(n-1),n])
  cumul.rgb <- col2rgb(cumul.col)/255
  
  cumul.dist <- as.numeric(sub.dist[1:(n-1),n])
  
  sapply(1:length(cumul.dist), function(n1){
  
  temp.lab <- cumul.dist[n1]
  
  if(length(cumul.dist)==1 | temp.lab == min(cumul.dist)){
    
  temp.lab <- sprintf("%.2f", temp.lab)
  text(y = circle.y + arrow.offset[n],
       x = (circle.cents[-1][1:(n-1)] - 0.5*diff(circle.cents[1:2]))[n1],
       labels = bquote(underline(.(temp.lab))),
       col=rgb(cumul.rgb[1], cumul.rgb[2], cumul.rgb[3], 1),
       cex=0.8)
  } else {
    temp.lab <- sprintf("%.2f", temp.lab)
    text(y = circle.y + arrow.offset[n],
         x = (circle.cents[-1][1:(n-1)] - 0.5*diff(circle.cents[1:2]))[n1],
         labels = temp.lab,
         col=rgb(cumul.rgb[1], cumul.rgb[2], cumul.rgb[3], 0.5),
         cex=0.8)
  }
  
  })
  
})
text(y=circle.y + max(arrow.offset) - 0.055,
     x=mean(circle.cents), col=cumul.col, pos=3, offset=0.75,
     labels='Cumululative dissimilarity', cex=0.8)  

sapply(1:circle.n, function(n){
    x = circle.cents[n]
    draw.circle(y=circle.y, x=x,radius=0.04, col="white", lwd=1)
    text(y=circle.y, x=x, labels=paste0('t',n), font=2, cex=0.8)
  })
box()

close.screen(6)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}
