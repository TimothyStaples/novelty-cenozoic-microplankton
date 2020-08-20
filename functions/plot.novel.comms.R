plot.novel.comm <- function(site.sp.mat, alpha, metric, site, axis.label){
  
return.data <- identify.novel.gam(site.sp.mat, alpha, metric, site=site, plot=FALSE,
                                  plot.data=TRUE)

min.p <- return.data[[3]]
seq.p <- return.data[[2]]
save.data <- return.data
return.data<-return.data[[1]]

ylims <- c(max(c(max(seq.p$upr, na.rm=TRUE), max(return.data$seq.dist, na.rm=TRUE))) * 1.1,
           min(c(min(seq.p$lwr, na.rm=TRUE), min(return.data$seq.dist, na.rm=TRUE))) *0.9)
           
plot(return.data$seq.dist ~ 
       as.numeric(as.character(return.data$bins)), type="n",
     ylim=ylims,
     axes=FALSE, xlab="", ylab="", yaxt="n")
axis(side=2, lwd=0.5)
lims <- par("usr")

polygon(x=c(as.numeric(as.character(return.data$bins)),
            rev(as.numeric(as.character(return.data$bins)))),
        y=c(seq.p$lwr, rev(seq.p$upr)),
        col="grey75", border=NA)

lines(seq.p[,4] ~ as.numeric(as.character(return.data$bins)), col="grey15",
      lty="dashed")

with(return.data,
     lines(seq.dist ~ as.numeric(as.character(bins)), lwd=1.5))

with(return.data[return.data$instant & !is.na(return.data$instant),],
     points(seq.dist ~ as.numeric(as.character(bins)),
            pch=21, bg="red"))

sapply(which(return.data$novel), function(x){
  segments(x0 = as.numeric(as.character(return.data$bins))[x],
           x1 = as.numeric(as.character(return.data$bins))[x],
           y0 = return.data$seq.dist[x] + (0.05 * (par("usr")[3] - par("usr")[4])),
           y1 = par("usr")[3], col="orange", lwd=2)
})

segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[4], y1=par("usr")[4])
segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])
mtext(side=2, text = "Instantaneous\ndissimilarity", line=2)

ylims <- c(min(c(min(min.p$lwr, na.rm=TRUE), min(return.data$raw.min.dist, na.rm=TRUE))) *0.9,
           max(c(max(min.p$upr, na.rm=TRUE), max(return.data$raw.min.dist, na.rm=TRUE))) * 1.1)

par(xpd=NA)
legend(x=relative.axis.point(0.5, "x"),
       y=relative.axis.point(1.175, "y"),
       legend=c("Instantaneous novelty", "Cumulative novelty", "Novel community"),
       pch=21, pt.bg=c("red","skyblue", "orange"),
       xjust=0.5, y.intersp=0, bty="n", x.intersp=0.75,
       horiz=TRUE)
par(xpd=FALSE)

plot(y=return.data$raw.min.dist,
     x=as.numeric(as.character(return.data$bins)), type="n",
     ylim=ylims,
     xlim=c(lims[1], lims[2]), xaxs="i",
     axes=FALSE, ylab="", xlab="")

polygon(x=c(as.numeric(as.character(return.data$bins)),
            rev(as.numeric(as.character(return.data$bins)))),
        y=c(min.p$lwr, rev(min.p$upr)),
        col="grey75", border=NA)

lines(min.p[,4] ~ as.numeric(as.character(return.data$bins)), col="grey15",
      lty="dashed")

with(return.data,
     lines(raw.min.dist ~ as.numeric(as.character(bins)), lwd=1.5))

with(return.data[return.data$cumul,],
     points(raw.min.dist ~ as.numeric(as.character(bins)),
            pch=21, bg="skyblue"))

par(xpd=NA)
sapply(which(return.data$novel), function(x){
  segments(x0 = as.numeric(as.character(return.data$bins))[x],
           x1 = as.numeric(as.character(return.data$bins))[x],
           y0 = return.data$raw.min.dist[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
           y1 =par("usr")[4], col="orange", lwd=2)
  
  points(x=as.numeric(as.character(return.data$bins))[x],
         y=par("usr")[4], pch=21, bg="orange")
})
par(xpd=FALSE)

segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[3])
segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])

axis(side=1, mgp=c(3,0.2,0), lwd=0.5)
axis(side=1, 
     tcl=-0.125, labels=NA)
mtext(side=1, text = axis.label, line=1)

axis(side=2, lwd=0.5)
mtext(side=2, text = "Cumulative\ndissimilarity", line=2)


return(return.data)
}