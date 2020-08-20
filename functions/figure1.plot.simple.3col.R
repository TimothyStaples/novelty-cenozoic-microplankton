figure1.plot.simple.3col <- function(prob.model.list,
                                example.ssmat,
                                plot.name,
                                example.xlims,
                                novel.lab.pos = list(3, 1, 4)){
  
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
    close.screen(all.screens=TRUE)
    }

  
  pdf(date.wrap(paste0("./plots/probability venn composite (",
                       plot.name, ")"), ".pdf"),
      height=10, width=5)

  # split.screen(rbind(c(0.65,1,0.15,0.99), # vens
  #                    
  #                    c(0.682,0.782,0.04,0.225), # smallplots
  #                    c(0.882,0.982,0.04,0.225), # smallplots
  #                    c(0.782,0.882,0.04,0.225), # smallplots
  #                    
  #                    c(0.65,0.99,0.01,0.99), # venn box & label
  #                    
  #                    c(0.375,0.5,0.13,0.99), # instant
  #                    c(0.5,0.625,0.13,0.99), # cumul
  #                    c(0.375,0.625,0.13,0.99),  # novel box & label
  #                    
  #                    c(0.05,0.325,0.01,0.99))) # time-series

  split.screen(rbind(c(0.15,0.95,0.05,0.375), # vens
                     
                     c(0.275,0.475,0.035,0.165), # smallplots
                     c(0.675,0.875,0.035,0.165), # smallplots
                     c(0.475,0.675,0.035,0.165), # smallplots
                     
                     c(0.15,0.95,0.01,0.4), # venn box & label
                     
                     c(0.15,0.95,0.45,0.6), # instant
                     c(0.15,0.95,0.6,0.75), # cumul
                     c(0.15,0.95,0.45,0.75),  # novel box & label
                     
                     c(0.15,0.95,0.75,0.99))) # time-series
  
  screen(1)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.4,0.6,0.5)
  circle.radius <- c(0.2, 0.175)
  circle.y <- 0.65
  
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
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  text(x=0.5225, y=circle.y+0.35, 
       adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=1)
  
  
  par(xpd=FALSE)
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
  
  text(x=0.25, y=circle.y+0.2,
       labels=c("Instantaneous\nnovelty (I)   "), 
       col="red", adj=1, cex=1)
  text(x=0.735, y=circle.y+0.2,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=1)
  
  overall.means <- t(sapply(rand.preds, function(x){
    
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  text(x=c(circle.cent + c(-0.075, 0.085, 0.01)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1)
  
    text(x = c(circle.cent + c(-0.075, 0.085, 0.01)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), " - ",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=0.8)
  close.screen(1)
  
  sapply(2:4, function(n){
    screen(n)
    
    par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0), las=1)
    
    point.col <- c("red", cumul.col, "orange")[n-1]
    
    pred.df <- as.data.frame(model.preds[[n-1]])
    
    ylims <- c(0, 0.2)
       
    plot(x=NULL, y=NULL, xlab="", ylab="",
         xlim=c(0.5,4.75), xaxs="i",
         ylim=ylims, axes=FALSE)
    
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=par("usr")[3], ytop=par("usr")[4], 
         col="white", border=NA)
    
    if(n == 2){
      axis(side=2, at=seq(0,0.2,0.05), cex.axis=0.8)
      mtext(side=2, line=1.75, text="Probability", cex=0.8, las=0)
      axis(side=4, at=seq(0,0.2,0.05), labels=NA)
    } else {
      axis(side=2, at=seq(0,0.2,0.05), labels=NA)
    }
    
    if(n == 4){
      axis(side=4, at=seq(0,0.2,0.05), labels=NA, tcl=0.25)
    }
    
    axis(side=1, at = 1:4, labels=c("D", "F", "N", "R"),
         mgp=c(3,0,0), cex.axis=0.8)
    
    pred.df$upper <- plogis(pred.df[,"fit"] + 1.96 * pred.df[, "se.fit"])
    pred.df$lower <- plogis(pred.df[,"fit"] - 1.96 * pred.df[, "se.fit"])
    
    raw.data <- prob.model.list$fixed.prob.models[[n]]$model$data
    raw.data$success <- raw.data[, c("instant", "cumul", "novel")[(n-1)]]
    raw.data$failure <- raw.data[, c("non.instant", "non.cumul", "non.novel")[(n-1)]]
    raw.data$prop <- raw.data$success / (raw.data$success + raw.data$failure)
    
    points(x = jitter(as.numeric(raw.data$taxa), amount=0.2),
           y = raw.data$prop, pch = 16, col=point.col, cex=0.5)
    
    segments(x0=1:4, x1=1:4,
             y0=pred.df$lower,
             y1=pred.df$upper)
    
    points(y=plogis(pred.df$fit), x=1:4, bg="black",
           pch=21)
    
    box()
    close.screen(n)
  })
  
  screen(5)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  text(x=relative.axis.point(0.03, "x"),
       y=relative.axis.point(0.97, "y"),
       labels="(C)", font=2)
  box()
  close.screen(5)
  
  screen(6)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0), las=1)
  
  
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
  
  plot(return.data$seq.dist ~ as.numeric(as.character(return.data$bins)),
       type="n", ylim=c(0,1), xlim=example.xlims,
       axes=FALSE, xlab="", ylab="", yaxt="n")
  axis(side=2, lwd=0.5, at=seq(0,0.8,0.2))
  lims <- par("usr")
  
  axis(side=1, las=1, mgp=c(3,0.2,0))
  mtext(side=1, line=1.25, text="Millions of years BP")
  axis(side=1, at=seq(1,65,1), tcl=-0.125, labels=NA)
  
  polygon(x=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          y=c(seq.p$lwr, rev(seq.p$upr)),
          col="grey75", border=NA)
  
  lines(y=seq.p[,4], x=as.numeric(as.character(return.data$bins)), col="grey15",
        lty="dashed")
  
  with(return.data,
       lines(y=seq.dist, x=as.numeric(as.character(bins)), lwd=0.75))
  
  with(return.data[return.data$instant & !is.na(return.data$instant),],
       points(y=seq.dist, x=as.numeric(as.character(bins)),
              pch=21, bg="red"))
  
  sub.data <- return.data[return.data$instant & !is.na(return.data$instant),]
  
  if(sum(sub.data$instant) != sum(sub.data$novel)){
  par(xpd=NA)
  with(sub.data[sub.data$cat != "novel",],
       text(y=seq.dist, x=as.numeric(as.character(bins)),
            labels="(I)", pos=novel.lab.pos[[1]],
            col="red"))
  par(xpd=FALSE)
  }
  
  sapply(which(return.data$novel), function(x){
    segments(x0 = as.numeric(as.character(return.data$bins))[x],
             x1 = as.numeric(as.character(return.data$bins))[x],
             y1 = return.data$seq.dist[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
             y0 = par("usr")[4], col="orange", lwd=2)
  })
  
  mtext(side=2, text = "Instantaneous\ndissimilarity", line=2, las=0)
  
  close.screen(6)
  
  screen(7)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0), las=1)
  
  plot(y=return.data$raw.min.dist,
       x=as.numeric(as.character(return.data$bins)), type="n",
       ylim=rev(c(0, 0.9)),
       xlim=c(lims[1], lims[2]),
       xaxs="i",
       axes=FALSE, ylab="", xlab="")
  
  polygon(x=c(as.numeric(as.character(return.data$bins)),
              rev(as.numeric(as.character(return.data$bins)))),
          y=c(min.p$lwr, rev(min.p$upr)),
          col="grey75", border=NA)
  
  lines(y=min.p[,4], x=as.numeric(as.character(return.data$bins)), col="grey15",
        lty="dashed")
  
  with(return.data,
       lines(y=raw.min.dist, x=as.numeric(as.character(bins)), lwd=1))
  
  with(return.data[return.data$cumul,],
       points(y=raw.min.dist, x=as.numeric(as.character(bins)),
              pch=21, bg=cumul.col))
  
  sub.data <- return.data[return.data$cumul,]
  
  par(xpd=NA)
  if(sum(sub.data$cumul) != sum(sub.data$novel)){
  with(sub.data[sub.data$cat != "novel",],
       text(y=raw.min.dist, x=as.numeric(as.character(bins)),
            labels="(C)", pos=novel.lab.pos[[2]],
            col=cumul.col))
  }
  
  sapply(1:sum(return.data$novel), function(n){
    x <- which(return.data$novel)[n]
    
    segments(x0 = as.numeric(as.character(return.data$bins))[x],
             x1 = as.numeric(as.character(return.data$bins))[x],
             y0 = return.data$raw.min.dist[x] - (0.05 * (par("usr")[4] - par("usr")[3])),
             y1 = par("usr")[3], col="orange", lwd=2)
    
    points(x=as.numeric(as.character(return.data$bins))[x],
           y=par("usr")[3], pch=21, bg="orange")
    
    text(x=as.numeric(as.character(return.data$bins))[x],
         y=par("usr")[3], labels="(N)", pos=novel.lab.pos[[3]][n], 
         col="orange", offset=0.25)
  })
  par(xpd=FALSE)
  
  axis(side=2, at=seq(0,0.8,0.2))

  mtext(side=2, text = "Cumulative\ndissimilarity", line=2, las=0)
  
  text(x=relative.axis.point(0.0125, "x"),
       y=relative.axis.point(0.925, "y"),
       labels="(B)", font=2, adj=0)
  
  close.screen(7)
  
  screen(8)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot.new()
  box()
  close.screen(8)
  
  screen(9)
  par(mar=c(0,0,0,0), ps=10, tcl = -0.25, mgp=c(3,0.5,0))
  plot.new()
  
  circle.y <- 0.275
  circle.n <- 6
  circle.cents <- relative.axis.point(seq(0.1,0.75,len=circle.n), "x")
  arrow.offset <- seq(0.075, 0.6, len=circle.n)

  segments(x0=min(circle.cents), x1 = max(circle.cents), y0 = circle.y, y1 = circle.y, lwd=2)
  segments(y0=circle.y, y1 = circle.y, x0 = max(circle.cents), x1 = max(circle.cents)+0.1, 
           lty="32", lwd=2)
  Arrows(y0=circle.y, y1 = circle.y, x0 = max(circle.cents)+0.09, x1 = max(circle.cents)+0.1,
         arr.type="triangle", arr.width=0.2, arr.length=0.2)
  
  par(lheight=0.75)
  text(y=circle.y, x=max(circle.cents)+0.1,
       pos=4, labels="Time\nseries")
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
           x0 = circle.cents[-1] - 0.015,
           x1 = circle.cents[-1] - 0.015, col="red", lwd=1.5)
  
  segments(y0 = circle.y  - arrow.offset[2],
           y1 = circle.y - arrow.offset[2],
           x0 = circle.cents[-1] - 0.015,
           x1 = circle.cents[-1] - diff(circle.cents[1:2]) + 0.015, col="red", lwd=1.5)
  
  Arrows(y0 = circle.y  - arrow.offset[2],
         y1 = circle.y  - arrow.offset[2] + 0.05,
         x0 = circle.cents[-1] - diff(circle.cents[1:2]) + 0.015,
         x1 = circle.cents[-1] - diff(circle.cents[1:2]) + 0.015, 
         col="red", lwd=1.5,
         arr.type = "triangle", arr.length=0.1, arr.width=0.1)
  
  rect(xleft=rowMeans(cbind(circle.cents[-1], circle.cents[-length(circle.cents)]))-0.035,
       xright=rowMeans(cbind(circle.cents[-1], circle.cents[-length(circle.cents)]))+0.035,
       ybottom=circle.y  - arrow.offset[2] - 0.05, 
       ytop=circle.y  - arrow.offset[2] + 0.05,
       border=NA, col=rgb(1,1,1,1))
  
  labels=sprintf("%.2f",diag(sub.dist[-1,-ncol(sub.dist)]))
  sapply(1:length(labels), function(n){
    temp.lab <- labels[n]
    
    text(y = (circle.y  - arrow.offset[2]),
       x = (rowMeans(cbind(circle.cents[-1], circle.cents[-length(circle.cents)])))[n],
       labels = bquote(underline(.(temp.lab))), col="red", 
       adj=c(0.5,0.5), cex=0.8)
  })
  
  par(xpd=NA)
  text(y=circle.y  - arrow.offset[2],
       x=mean(circle.cents), col="red", pos=1, offset=0.75,
       labels='Instantaneous dissimilarity ("rate-of-change")')
  
  text(x=relative.axis.point(0.015, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(A)", font=2, adj=0)
  par(xpd=FALSE)

sapply(2:length(circle.cents), function(n){
    
  segments(y0 = circle.y,
           y1 = circle.y + arrow.offset[n],
           x0 = circle.cents[n],
           x1 = circle.cents[n], col=cumul.col, lwd=1.5) 
    
  segments(y0 = circle.y + arrow.offset[n],
           y1 = circle.y + arrow.offset[n],
           x0 = circle.cents[1],
           x1 = circle.cents[n], col=cumul.col, lwd=1.5)
  
  Arrows(y0 = circle.y + arrow.offset[n],
         y1 = circle.y + arrow.offset[n] - 0.05,
         x0 = circle.cents[-1][1:(n-1)] - diff(circle.cents[1:2]),
         x1 = circle.cents[-1][1:(n-1)] - diff(circle.cents[1:2]), 
         col=cumul.col, lwd=1.5,
         arr.type = "triangle", arr.length=0.1, arr.width=0.1)
  
  rect(xleft = circle.cents[-1][1:(n-1)] - 0.5*diff(circle.cents[1:2]) - 0.035,
       xright = circle.cents[-1][1:(n-1)] - 0.5*diff(circle.cents[1:2]) + 0.035,
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
text(y=circle.y + max(arrow.offset),
     x=mean(circle.cents), col=cumul.col, pos=3, offset=0.75,
     labels='Cumululative dissimilarity ("novelness")')  

sapply(1:circle.n, function(n){
    x = circle.cents[n]
    draw.circle(y=circle.y, x=x,radius=0.05, col="white", lwd=2)
    text(y=circle.y, x=x, labels=paste0('t',n), font=2)
  })
box()

  close.screen(9)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}