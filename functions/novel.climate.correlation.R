novel.climate.correlation <- function(novel.list){

zachos <- read.csv("./raw.datafiles/zachos2001.csv")
zachos$genus <- gsub(" ", "", as.character(zachos$genus))

bins <- seq(0,66, 0.1)

zachos.100k <- cut(zachos$age.mya, breaks = bins)

zachos.agg <- data.frame(bin = levels(zachos.100k),
                         d18O = as.vector(tapply(zachos$d18O,
                                       zachos.100k,
                                       mean, na.rm=TRUE)))

age.bins <- as.character(zachos.agg$bin)
zachos.agg$age <- round(as.numeric(substr(age.bins, regexpr(",", age.bins)+1,
                                                        nchar(age.bins)-1)) - 0.5*(bins[2]-bins[1]),2)

zachos.agg <- zachos.agg[order(zachos.agg$age, decreasing=TRUE),]
zachos.agg <- zachos.agg[!is.na(zachos.agg$d18O),]
zachos.agg$diff.d18O <- c(NA, diff(zachos.agg$d18O))
zachos.agg$abs.diff.d18O <- abs(zachos.agg$diff.d18O)
zachos.agg$bin.lag <- round(abs(c(NA, diff(zachos.agg$age))),2)

# model novel probability over time ####

novel.df <- do.call("rbind", lapply(1:length(novel.list), function(n){

  x <- novel.list[[n]]
  
  x.df <- do.call("rbind", x$novel)
  x.df$taxa <- c("nano", "foram", "radio", "diatom")[n]

  return(x.df)    
}))
novel.df$bin.num <- as.numeric(as.character(novel.df$bins))

novel.prop <- data.frame(novel.prop = tapply(novel.df$novel,
                                                cut(novel.df$bin.num, breaks=bins),
                                                mean))
novel.prop$age <- rowMeans(cbind(bins[-1],
                                 bins[-length(bins)]))

zachos.novel <- merge(novel.prop, zachos.agg,
                      all.x=TRUE, all.y=FALSE, sort=FALSE)

zachos.novel.long <- merge(novel.df, zachos.agg,
                           by.x = "bins", by.y="age",
                           all.x=TRUE, all.y=FALSE, sort=FALSE)
zachos.novel.long$scaled.n <- as.vector(scale(zachos.novel.long$bin.n))
zachos.novel.long$scaled.lag <- as.vector(scale(zachos.novel.long$bin.lag.x))
zachos.novel.long$taxa <- as.factor(zachos.novel.long$taxa)
zachos.novel.long$site = as.factor(zachos.novel.long$site)

# novel.trends <- gamm4(novel ~ s(bin.num, bs="cr", k=66, by=taxa) + scaled.n + scaled.lag + taxa,
#                       random= ~(1|site),
#                       family=binomial, data = zachos.novel.long)

novel.trends <- gam(novel ~ s(bin.num, bs="cr", k=66, by=taxa) + scaled.n + scaled.lag + taxa,
                   family=binomial, data = zachos.novel.long)

zachos.agg$scaled.lag <- as.vector(scale(zachos.agg$age))

temp.trends <- gam(abs.diff.d18O ~ s(age, bs="cr", k=66), data=zachos.agg)

novel.pred.df <- lapply(1:4, function(n){
  
temp.pred <-  data.frame(bin.num = seq(0,66,0.1),
                      scaled.n=0,
                      scaled.lag=0,
             taxa = levels(zachos.novel.long$taxa)[n])
  
cbind(temp.pred,
      as.data.frame(predict(novel.trends, newdata=temp.pred, se.fit=TRUE)))
})

zachos.pred.df <- data.frame(age = seq(0,66,0.1),
                             scaled.lag=0)
zachos.pred.df <- cbind(zachos.pred.df,
                       as.data.frame(predict(temp.trends, newdata=zachos.pred.df, se.fit=TRUE)))

# Direct novelty correlation ####

head(zachos.novel.long)

temp.corr <- glmer(novel ~ abs.diff.d18O * taxa + scaled.lag + scaled.n + (1|site),
                   data=zachos.novel.long, family=binomial)



}
