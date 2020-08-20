turnover.edge.effects <- function(turnover.df){

require(merTools)
  
print("--- Testing extinction edge effects ---")
turnover.df$n.from.start <- turnover.df$n.from.start

turnover.df$log.start <- log(turnover.df$n.from.start)
turnover.df$scale.start <- scale(turnover.df$log.start)
turnover.df$scale.start2 <- scale(turnover.df$log.start^2)

turnover.df$log.end <- log(turnover.df$n.from.end)
turnover.df$scale.end <- scale(turnover.df$log.end)
turnover.df$scale.end2 <- scale(turnover.df$log.end^2)

turnover.df$taxa <- as.factor(turnover.df$taxa)

print("Running edge effect models...")

start.model <- glmer(cbind(orig, orig.rich) ~ scale.start + scale.start2 + (1|taxa), 
                     data = turnover.df[turnover.df$n.from.start > 0,], 
                     family = binomial)

end.model <- glmer(cbind(ext, ext.rich) ~ scale.end + scale.end2 + (1|taxa), 
                   data=turnover.df[!is.na(turnover.df$cat.aft),], 
                   family=binomial)

start.pred <- turnover.df[!duplicated(turnover.df$scale.start),
                          c("scale.start", "scale.start2")]
start.pred <- start.pred[order(start.pred$scale.start),]
start.pred$taxa = "a"

start.pred <- cbind(start.pred,
                    predictInterval(start.model, newdata = start.pred,
                                    n.sims = 999, which="fixed",
                                    level=0.95, include.resid.var=FALSE,
                                    type="linear.prediction"))

end.pred <- turnover.df[!duplicated(turnover.df$scale.end),
                          c("scale.end", "scale.end2")]
end.pred <- end.pred[order(end.pred$scale.end),]
end.pred$taxa = "a"

end.pred <- cbind(end.pred,
                    predictInterval(end.model, newdata = end.pred,
                                    n.sims = 999, which="fixed",
                                    level=0.95, include.resid.var=FALSE,
                                    type="linear.prediction"))

pdf(date.wrap("./plots/turnover edge effects", ".pdf"),
    height=3, width=5.5, useDingbats=FALSE)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(3,3,1,1), ps=10, tcl=-0.25,
    mgp=c(3,0.2,0), las=1)

start.pred$bin.start <- exp(start.pred$scale.start * 
                            attr(turnover.df$scale.start, "scaled:scale") +
                            attr(turnover.df$scale.start, "scaled:center")) - 1

end.pred$bin.end <- exp(end.pred$scale.end * 
                              attr(turnover.df$scale.end, "scaled:scale") +
                              attr(turnover.df$scale.end, "scaled:center")) 

plot(plogis(end.pred$fit) ~ end.pred$bin.end, 
     type = 'n', xlim = log(c(0.7,600)), ylim=c(0,0.35),
     xlab = "", ylab = "", xaxt="n", yaxt="n", xaxs="i")

per.ts.prop <- do.call("rbind", lapply(split(turnover.df, f=turnover.df$n.from.end),
                      function(x){
                        
                        data.frame(prop = mean(x$ext / (x$ext.rich + x$ext)),
                                   sample = nrow(x),
                                   x.pos = x$n.from.end[1])

                      }))

points(per.ts.prop$prop ~ log(per.ts.prop$x.pos),
       cex = log(per.ts.prop$sample+1) / max(log(per.ts.prop$sample))*1.5,
       pch=16, 
       col=rgb(0.5,0.5,0.5,0.1),
       lwd=0.1)

lines(plogis(end.pred$upr) ~ log(end.pred$bin.end), lty="dashed")
lines(plogis(end.pred$lwr) ~ log(end.pred$bin.end), lty="dashed")
lines(plogis(end.pred$fit) ~ log(end.pred$bin.end), lwd=2)

axis(side=1, at=log(c(1,10,100,500)), labels=c(1,10,100,500))
axis(side=1, at=log(c(seq(1,10,1),
                      seq(10,100,10),
                      seq(100,1000,100))), labels=NA, tcl=-0.125)
mtext(side=1, line=1.25, text="Bins from time-series end")

axis(side=2, mgp=c(3,0.5,0))
axis(side=2, at=seq(0,0.4,0.05), labels=NA, tcl=-0.125)
mtext(side=2, line=2, text="Probability", las=0)


text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(A) Local extinction", font=2, adj=0)

plot(plogis(start.pred$fit) ~ start.pred$bin.start, 
     type = 'n', xlim = log(c(0.7,600)), ylim=c(0,0.35),
     xlab = "", ylab = "", xaxt="n", yaxt="n", xaxs="i")

per.ts.prop <- do.call("rbind", lapply(split(turnover.df, f=turnover.df$n.from.start),
                                       function(x){
                                         
                                         data.frame(prop = mean(x$orig / (x$orig.rich + x$orig), na.rm=TRUE),
                                                    sample = nrow(x),
                                                    x.pos = x$n.from.start[1])
                                         
                                       }))

points(per.ts.prop$prop ~ log(per.ts.prop$x.pos - 1),
       cex = log(per.ts.prop$sample+1) / max(log(per.ts.prop$sample))*1.5,
       pch=16, 
       col=rgb(0.5,0.5,0.5,0.1),
       lwd=0.1)

with(start.pred[start.pred$bin.start > 0, ], {
lines(plogis(upr) ~ log(bin.start), lty="dashed")
lines(plogis(lwr) ~ log(bin.start), lty="dashed")
lines(plogis(fit) ~ log(bin.start), lwd=2)
      })

axis(side=1, at=log(c(1,10,100,500)), labels=c(1,10,100,500))
axis(side=1, at=log(c(seq(1,10,1),
                      seq(10,100,10),
                      seq(100,1000,100))), labels=NA, tcl=-0.125)
mtext(side=1, line=1.25, text="Bins from time-series start")

axis(side=2, mgp=c(3,0.5,0), labels=NA)
axis(side=2, at=seq(0,0.4,0.05), labels=NA, tcl=-0.125)
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(B) Local origination", font=2, adj=0)


dev.off()

}