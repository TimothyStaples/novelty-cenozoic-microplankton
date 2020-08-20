figure3.plot.fixed <- function(turnover.models, plot.name,
                               left.xlim, left.ylim, right.xlim, right.ylim){
  
  
  print("Predicting from models...")
  cc.model.preds <- lapply(1:4, function(n){
    
    print(n)
    temp.model <- turnover.models[[n]]$model
    
    pred.df <- expand.grid(cat.bef=levels(temp.model@frame$cat.bef),
                           cat.aft=levels(temp.model@frame$cat.aft),
                           taxa = levels(temp.model@frame$taxa))
    pred.df$lag.scaled = 0
    
    if(n %in% 1:2){
      pred.df$edge.scaled = 0 
      pred.mat <- model.matrix(~ (cat.bef + cat.aft) * taxa + lag.scaled + edge.scaled, data=pred.df)
    } else {
      pred.mat <- model.matrix(~ (cat.bef + cat.aft) * taxa + lag.scaled, data=pred.df)  
    }
    
    glht.preds <- glht(temp.model, linfct=pred.mat)
    glht.preds <- as.data.frame(summary(glht.preds)$test[c("coefficients", 
                                                           "sigma")])
    colnames(glht.preds) <- c("fit", "se.fit")
    
    model.preds <-cbind(pred.df, glht.preds)
    model.preds$comb <- as.factor(paste0(model.preds$cat.bef, model.preds$cat.aft))
    model.preds$model <- c("ext", "orig", "emig", "immig")[n]
    
    model.preds$estimate <- plogis(model.preds$fit)
    model.preds$upper <- plogis(model.preds$fit + 1.96 * model.preds$se.fit)
    model.preds$lower <- plogis(model.preds$fit - 1.96 * model.preds$se.fit)
    
    return(model.preds)
    
  })
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  library(shape)
  library(plotrix)
  
  pdf(date.wrap(paste0("./plots/cc combined - temp vs perm (",
                       plot.name,
                       ")"),
                ".pdf"), 
      height=9, width=5, useDingbats = FALSE)
  

  framemat <- rbind(c(0.125,0.55,0.7,0.9),
                    c(0.55,0.975,0.7,0.9),
                    c(0.125,0.55,0.5,0.7),
                    c(0.55,0.975,0.5,0.7),
                    
                    c(0.125,0.55,0.25,0.45),
                    c(0.55,0.975,0.25,0.45),
                    c(0.125,0.55,0.05,0.25),
                    c(0.55,0.975,0.05,0.25),
                    
                    c(0.1,0.95,0.9,1))
  

  split.screen(framemat)
  
  # EXT/EMIG PLOTS ####
  sapply(1:4, function(n){
  
  screen(n)
  
  sub.preds.x <- cc.model.preds[[1]]
  sub.preds.x <- sub.preds.x[sub.preds.x$taxa == levels(sub.preds.x$taxa)[n],]
  sub.preds.y <- cc.model.preds[[3]]
  sub.preds.y <- sub.preds.y[sub.preds.y$taxa == levels(sub.preds.y$taxa)[n],]
  
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  plot(x=NULL, y=NULL, xlim=left.xlim, ylim=left.ylim,
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  
  abline(a=0, b=1, lwd=2, col="grey60", lty="31")

  if(n %in% 3:4){
  axis(side=1, at=seq(0,0.3,0.1), mgp=c(3,0.2,0))
    } else {
    axis(side=1, at=seq(0,0.3,0.1), mgp=c(3,0.2,0), labels=NA)  
  }
  
  if(n %in% c(1,3)){
  axis(side=2)
  } else {
    axis(side=2, labels=NA)
  }
  
  if(n==1){mtext(side=2, line=2, text="Probability of emigration", las=0, at=par("usr")[3])}
  if(n==3){mtext(side=1, line=1, text="Probability of local extinction", las=0, at=par("usr")[2])}
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(ifelse(n %in% 1:2, 0.935, 0.915), "y"),
       adj=0, labels=paste0("(", LETTERS[n], ") ",
       c("Diatoms", "Foraminifera", "Nannoplankton", "Radiolarians")[n]), font=2)
  
  segments(y0=sub.preds.y$estimate,
           y1=sub.preds.y$estimate,
           x0=sub.preds.x$upper,
           x1=sub.preds.x$lower,
           col="black", lwd=1)
  
  segments(x0=sub.preds.x$estimate,
           x1=sub.preds.x$estimate,
           y0=sub.preds.y$upper,
           y1=sub.preds.y$lower,
           col="black", lwd=1)
  
  sapply(1:dim(sub.preds.x)[1], function(x){
    print(x)
    
    aft.col = c("grey35", cumul.col, "red", "orange")[sub.preds.x$cat.aft[x]]
    bef.col = c("grey35", cumul.col, "red", "orange")[sub.preds.x$cat.bef[x]]
    shape = c("triangle", "square", "circle", "diamond")[sub.preds.x$taxa[x]]
    
    arrow.shape(x=sub.preds.x$estimate[x],
                y=sub.preds.y$estimate[x],
                r=0.175, screen=framemat[n,], rads=c(1.5, 0.5),
                col=aft.col,
                shape=shape,
                border=NA, plot=TRUE)
    
    arrow.shape(x=sub.preds.x$estimate[x],
                y=sub.preds.y$estimate[x],
                r=0.175, screen=framemat[n,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape=shape,
                col=bef.col,
                plot=TRUE, border="black", add.arrow=TRUE)
    
    arrow.shape(x=sub.preds.x$estimate[x],
                y=sub.preds.y$estimate[x],
                r=0.175, screen=framemat[n,], rads=c(0, 2),
                shape=shape,
                plot=TRUE, border="black", add.arrow=FALSE)
    
  })
    box()
  close.screen(n)
  
  })

  # ORIG/IMMIG PLOTS ####
  sapply(1:4, function(n){
  
  screen(n+4)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  sub.preds.x <- cc.model.preds[[2]]
  sub.preds.x <- sub.preds.x[sub.preds.x$taxa == levels(sub.preds.x$taxa)[n],]
  sub.preds.y <- cc.model.preds[[4]]
  sub.preds.y <- sub.preds.y[sub.preds.y$taxa == levels(sub.preds.y$taxa)[n],]
  
  plot(x=NULL, y=NULL, xlim=right.xlim, ylim=right.ylim,
       xaxs="i", yaxs="i", axes=FALSE, xlab="", ylab="")
  abline(a=0, b=1, lwd=2, col="grey60", lty="31")
  
  if(n %in% 3:4){
    axis(side=1, mgp=c(3,0.2,0))
  } else {
    axis(side=1, mgp=c(3,0.2,0), labels=NA)  
  }
  
  if(n %in% c(1,3)){
    axis(side=2)
  } else {
    axis(side=2, labels=NA)
  }
  
  if(n==1){mtext(side=2, line=2, text="Probability of immigration", las=0, at=par("usr")[3])}
  if(n==3){mtext(side=1, line=1, text="Probability of local origination", las=0, at=par("usr")[2])}
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(ifelse(n %in% 1:2, 0.935, 0.915), "y"),
       adj=0, labels=paste0("(", LETTERS[n+4], ") ",
                            c("Diatoms", "Foraminifera", "Nannoplankton", "Radiolarians")[n]), font=2)
  
  # UNCOMMENT THESE LINES TO GENERATE 95% CIS ON PLOT
  segments(y0=sub.preds.y$estimate,
           y1=sub.preds.y$estimate,
           x0=sub.preds.x$upper,
           x1=sub.preds.x$lower,
           col="black", lwd=1)
  
  segments(x0=sub.preds.x$estimate,
           x1=sub.preds.x$estimate,
           y0=sub.preds.y$upper,
           y1=sub.preds.y$lower,
           col="black", lwd=1)
  
  sapply(1:dim(sub.preds.y)[1], function(x){
    
    aft.col = c("grey35", cumul.col, "red", "orange")[sub.preds.x$cat.aft[x]]
    bef.col = c("grey35", cumul.col, "red", "orange")[sub.preds.x$cat.bef[x]]
    shape = c("triangle", "square", "circle", "diamond")[sub.preds.x$taxa[x]]
    
    
    arrow.shape(x=sub.preds.x$estimate[x],
                y=sub.preds.y$estimate[x],
                r=0.175, screen=framemat[n,], rads=c(1.5, 0.5),
                col=aft.col,
                shape=shape,
                border=NA, plot=TRUE)
    
    arrow.shape(x=sub.preds.x$estimate[x],
                y=sub.preds.y$estimate[x],
                r=0.175, screen=framemat[n,],
                rads=c(0.5, 1.5), lwd=0.5,
                shape=shape,
                col=bef.col,
                plot=TRUE, border="black", add.arrow=TRUE)
    
    arrow.shape(x=sub.preds.x$estimate[x],
                y=sub.preds.y$estimate[x],
                r=0.175, screen=framemat[n,], rads=c(0, 2),
                shape=shape,
                plot=TRUE, border="black", add.arrow=FALSE)
    
  })
  box()
  close.screen(n+4)
  })
  
  screen(9)
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot.new()
  
  par(xpd=NA)
  
  taxa.pos <- rev(seq(0.1,0.7, len=4))
  
  text(x=0.5, y=0.9, adj=0.5,
       labels=bquote(bold("Taxa groups")), font=2, cex=0.8)
  sapply(taxa.pos, function(y){
    
    arrow.shape(x=0.4,
                y=y,
                r=0.175, screen=framemat[9,], rads=c(1.5,0.5),
                col="white",
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE)

    arrow.shape(x=0.4,
                y=y,
                r=0.175, screen=framemat[9,], rads=c(0.5, 1.5),
                col="white",
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE, add.arrow=TRUE, lwd=0.5)

    arrow.shape(x=0.4,
                y=y,
                r=0.175, screen=framemat[9,], rads=c(0,2),
                shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                border="black", plot=TRUE)


  })

  text(x=0.4, pos=4, y=taxa.pos,
       labels=c("Diatoms", "Foraminifera", "Nannoplankton", "Radiolarians"),
       cex=0.8, offset=0.75)

  par(lheight=0.85)
  text(x=0.175, y=0.45, labels="Preceding\ncommunity", adj=0.8, cex=0.8,
       lheight=0.5)
  
  text(x=0.175, y=0.15, labels="Succeeding\ncommunity", adj=0.2, cex=0.8,
       lheight=0.5)
  par(lheight=1)
  
  arrow.shape(x=0.195, y=0.8,
              r=0.3, screen=framemat[9,], rads=c(1.5,0.5),
              col="white", shape="circle",
              border="black", plot=TRUE)
  
  arrow.shape(x=0.165, y=0.8,
              r=0.3, screen=framemat[9,], rads=c(0.5,1.5),
              col="white", shape="circle",
              border="black", plot=TRUE, add.arrow = TRUE)
  
  segments(x0=c(0.065, 0.05, 0.295,  0.31),
           x1=c(0.05,  0.05, 0.31, 0.31),
           y0=c(0.45, 0.45, 0.15,  0.15),
           y1=c(0.45, 0.65, 0.15,  0.65))
  
  Arrows(x0=c(0.05, 0.31),
         x1=c(0.115, 0.245),
         y0=c(0.65, 0.65),
         y1=c(0.775, 0.775),
         arr.length = 0.1, arr.width = 0.1,
         arr.type = "triangle")
  
  rect.pos <- rev(seq(0.1,0.6, len=4))
  rect.height <- 0.06
  
  rect(xleft=0.7, xright=0.73,
       ytop=rect.pos[1] + rect.height, ybottom=rect.pos[1] - rect.height, col="grey35")

  rect(xleft=0.7, xright=0.73,
       ytop=rect.pos[2] + rect.height, ybottom=rect.pos[2] - rect.height, col="red")
  
  rect(xleft=0.7, xright=0.73,
       ytop=rect.pos[3] + rect.height, ybottom=rect.pos[3] - rect.height, col=cumul.col)
  
  rect(xleft=0.7, xright=0.73,
       ytop=rect.pos[4] + rect.height, ybottom=rect.pos[4] - rect.height, col="orange")
  
  par(lheight=0.7)
  text(x=0.73, y=rect.pos, pos=4, offset=0.25, cex=0.8,
       labels=c("Background", "Instantaneous novelty",
                "Cumulative novelty", "Novel community"))
  par(lheight=1)
  
  text(x=0.875, y=0.9, adj=0.5,
       labels=bquote(bold("Community")), font=2, cex=0.8)
  text(x=0.875, y=0.775, adj=0.5,
       labels=bquote(bold("classification")), font=2, cex=0.8)
  
  par(xpd=FALSE)

  close.screen(9)
  close.screen(all.screens=TRUE)
  
  dev.off()
  
}