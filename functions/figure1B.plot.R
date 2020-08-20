figure1B.plot <- function(prob.model.list,
                          plot.name){
  
  model.preds <- lapply(c(4:6), function(n){
    temp.pred <- prob.model.list$fixed.prob.models[[n]]$pred.df
    return(temp.pred)
  })
  
  rand.preds <- lapply(c(4:6), function(n){
    temp.pred <- prob.model.list$random.prob.models[[n]]$pred.df
    return(temp.pred)
  })
  
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

  
  pdf(date.wrap(paste0("./plots/probability venn 1B (",
                       plot.name, ")"), ".pdf"),
      height=2.5, width=2.25, colormodel="cmyk", useDingbats = FALSE)

  split.screen(rbind(c(0.01,0.99,0.05,0.99), # vens
                     
                     c(0.18,0.44,0.08,0.4), # smallplots
                     c(0.7,0.96,0.08,0.4), # smallplots
                     c(0.44,0.7,0.08,0.4))) # smallplots
                     
  screen(1)
  par(mar=c(0,0,0,0), ps=8, tcl = -0.2, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.43,0.65,0.54)
  circle.radius <- c(0.22, 0.2)
  circle.y <- 0.65
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  par(xpd=NA)
  segments(x0=circle.cent + c(0, 0, 0.0175), 
           y0=rep(circle.y, 3) + c(0, 0, 0.24),
           x1=circle.cent + c(-0.15,0.2, 0.0175),
           y1=c(0.3, 0.3, 0.3),
           lwd=1)
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  par(lheight=0.85)
  text(x=circle.cent[3]+0.02, y=circle.y+0.3, 
       adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=0.9)
  
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
  plot(N.sp, col="orange", add=TRUE, border="black", lwd=1.5)
  
  text(x=0.3, y=circle.y+0.2,
       labels=c("Instantaneous\nnovelty (I)    "), 
       col="red", adj=1, cex=0.9)
  text(x=0.75, y=circle.y+0.2,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=0.9)
  
  overall.means <- t(sapply(rand.preds, function(x){
    
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  text(x=c(circle.cent + c(-0.075, 0.085, 0.01)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=0.9)
  
    text(x = c(circle.cent + c(-0.085, 0.1, 0.01)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=0.8)
    
    text(x=relative.axis.point(0.1, "x"),
         y=relative.axis.point(0.95, "y"),
         labels="(A)", font=2, cex=0.9)
    
  close.screen(1)
  
  sapply(2:4, function(n){
    screen(n)
    
    par(mar=c(0,0,0,0), ps=8, tcl = -0.2, mgp=c(3,0.4,0), las=1)
    
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
      mtext(side=2, line=1.3, text="Probability", cex=0.8, las=0)
      axis(side=4, at=seq(0,0.2,0.05), labels=NA)
    } else {
      axis(side=2, at=seq(0,0.2,0.05), labels=NA)
    }
    
    if(n == 4){
      axis(side=4, at=seq(0,0.2,0.05), labels=NA, tcl= 0.2)
    }
    
    axis(side=1, at = 1:4, labels=c("D", "F", "N", "R"),
         mgp=c(3,-0.1,0), cex.axis=0.8)
    
    pred.df$upper <- plogis(pred.df[,"fit"] + 1.96 * pred.df[, "se.fit"])
    pred.df$lower <- plogis(pred.df[,"fit"] - 1.96 * pred.df[, "se.fit"])
    
    raw.data <- prob.model.list$fixed.prob.models[[n]]$model$data
    raw.data$success <- raw.data[, c("instant", "cumul", "novel")[(n-1)]]
    raw.data$failure <- raw.data[, c("non.instant", "non.cumul", "non.novel")[(n-1)]]
    raw.data$prop <- raw.data$success / (raw.data$success + raw.data$failure)
    
    points(x = jitter(as.numeric(raw.data$taxa), amount=0.2),
           y = raw.data$prop, pch = 16, col=point.col, cex=0.4)
    
    segments(x0=1:4, x1=1:4,
             y0=pred.df$lower,
             y1=pred.df$upper)
    
    points(y=plogis(pred.df$fit), x=1:4, bg="white",
           pch=21, cex=0.75, lwd=1.25)
    
     text(x=relative.axis.point(0.15, "x"),
          y=relative.axis.point(0.925, "y"),
          labels=paste0("(",LETTERS[c(2,4,3)][n-1],")"), font=2, cex=0.9)
    
    box()
    close.screen(n)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}