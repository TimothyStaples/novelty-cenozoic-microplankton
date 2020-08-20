test.time.series.sensitivity <- function(all.novel.list,
                                         cut.novel.list){
  
# run models 
all.ts.models <- lapply(c("instant", "cumul", "novel"), function(x){
  
print(x)
temp.novel <- do.call("rbind", 
                        lapply(1:4, function(n){
                          
                          sub.novel <- all.novel.list[[n]]
                          
                          do.call("rbind", lapply(sub.novel$site, function(s){
                          
                          temp.df <- sub.novel$novel[[s]]
                          temp.df$bin.num <- as.numeric(as.character(temp.df$bins))
                          temp.df <- temp.df[order(temp.df$bin.num, decreasing=TRUE),]
                          temp.df$bin.n <- 1:nrow(temp.df)
                          temp.df$taxa <- c("nano", "foram", "radio", "diatom")[n]
                          return(temp.df)
                          
                        }))
  }))

temp.novel <- temp.novel[temp.novel$bin.n > 1,]
temp.novel$success <- temp.novel$cat == x
temp.novel$log.bin <- log(temp.novel$bin.n - 1) # correct for NA bin in all data

temp.model <- glmer(success ~ log.bin + I(log.bin^2) + (1|taxa/site), data=temp.novel,
                    family=binomial)
summary(temp.model)
return(temp.model)
})

cut.ts.models <- lapply(c("instant", "cumul", "novel"), function(x){  
  
  print(x)
  temp.novel <- do.call("rbind", 
                        lapply(1:4, function(n){
                          
                          sub.novel <- cut.novel.list[[n]]
                          
                          do.call("rbind", lapply(sub.novel$site, function(s){
                            
                            temp.df <- sub.novel$novel[[s]]
                            temp.df$bin.num <- as.numeric(as.character(temp.df$bins))
                            temp.df <- temp.df[order(temp.df$bin.num, decreasing=TRUE),]
                            temp.df$bin.n <- 1:nrow(temp.df)
                            temp.df$taxa <- c("nano", "foram", "radio", "diatom")[n]
                            return(temp.df)
                            
                          }))
                        }))
  
  temp.novel <- temp.novel[temp.novel$bin.n > 1,]
  temp.novel$success <- temp.novel$cat == x
  temp.novel$log.bin <- log(temp.novel$bin.n - 1) # correct for NA bin in all data
  
  temp.model <- glmer(success ~ log.bin + I(log.bin^2) + (1|taxa/site), data=temp.novel,
                      family=binomial)
  summary(temp.model)
  return(temp.model)
})

require(lme4)
require(merTools)

ts.freq.preds <- lapply(c(all.ts.models, cut.ts.models),
                        function(temp){
                          
                          pred.df <- data.frame(log.bin = seq(min(temp@frame$log.bin),
                                                              max(temp@frame$log.bin),
                                                              len=200))
                          
                          pred.df$raw.bin = exp(pred.df$log.bin)
                          pred.df$taxa <- temp@frame$taxa[1]
                          pred.df$site <- temp@frame$site[1]
                          
                          pred.df <- cbind(pred.df,
                                           sum.fit = predict(temp,
                                                             newdata=pred.df,
                                                             re.form=NA),
                                           predictInterval(temp, 
                                                           newdata = pred.df,
                                                           n.sims = 999,
                                                           which="fixed",
                                                           level=0.95,
                                                           include.resid.var=FALSE,
                                                           type="linear.prediction"))
                          
                          return(pred.df)
                        })

pdf(date.wrap("./plots/prob over ts length", ".pdf"), height=5.5, width=4.5, useDingbats=FALSE)
par(mfcol=c(3,2), mar=c(0,0,0,0), oma=c(2.5,3,2,2), tcl=-0.25, mgp=c(3,0.5,0), las=1, ps=8)

sapply(1:6, function(n){
  
  plot(x=NULL, y=NULL, xlim=c(1, 500), ylim=c(-0.005,0.225), yaxs="i",
       xlab="", ylab="", xaxt="n", yaxt="n", log="x")
  
  if(n %in% c(1,4)){temp.col = rgb(1,0,0,0.25)}
  if(n %in% c(2,5)){temp.col = rgb(0.373,0.651,0.765,0.25)}
  if(n %in% c(3,6)){
    nov.col <- col2rgb("orange")/255
    temp.col <- rgb(nov.col[1], nov.col[2], nov.col[3], 0.25)
  }
  
  if(n %in% c(1:3)){
    mtext(side=2, line=2, text="Probability", las=0)
    axis(side=2)
  } else {axis(side=2, labels=NA)}
  
  if(n %in% c(3,6)){
    axis(side=1, mgp=c(3,0.2,0))
    mtext(side=1, line=1, text="Time-series position")
  } else {axis(side=1, labels=NA)}
  
  axis(side=1, at=c(seq(1,10,1), seq(10,100,10), seq(100,500,100)),
       labels=NA, tcl=-0.125)
  
  if(n == 1){mtext(side=3, line=0.2, text="Complete time-series", font=2)}
  if(n == 4){mtext(side=3, line=0.2, text="First 5 bins removed", font=2)
             mtext(side=4, line=0.1, las=0, text="Instantaneous Novelty", font=2)}
  if(n == 5){mtext(side=4, line=0.1, las=0, text="Cumulative Novelty", font=2)}
  if(n == 6){mtext(side=4, line=0.1, las=0, text="Novel Community", font=2)}
  
  pred.df <- ts.freq.preds[[n]]  
  temp.m <- c(all.ts.models, cut.ts.models)[[n]]
  
  per.ts.prop <- tapply(temp.m@frame$success,
                        exp(temp.m@frame$log.bin),
                        function(x){sum(x) / length(x)})
  
  sample.size <- tapply(temp.m@frame$success,
                        exp(temp.m@frame$log.bin), length)
  
  points(per.ts.prop ~ as.numeric(names(per.ts.prop)),
         cex = log(sample.size+1) / max(log(sample.size))*1.5,
         pch=16, 
         col=temp.col,
         lwd=0.1)
  
  lines(x=pred.df$raw.bin,
        y=plogis(pred.df$upr),
        col="black", lwd=1, lty="31")
  
  lines(x=pred.df$raw.bin,
        y=plogis(pred.df$lwr),
        col="black", lwd=1, lty="31")
  
  lines(x=pred.df$raw.bin,
        y=plogis(pred.df$fit),
        col="black", lwd=2)
  
  text(x = 0.9,
       y = relative.axis.point(0.935, "y"),
       labels = paste0("(", LETTERS[n], ")"), adj = 0, font = 2, cex=1.15)
  
})
dev.off()

}