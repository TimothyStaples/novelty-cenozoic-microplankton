syn.test <- function(){

# Generate communities with different shapes, visualise using dissimilarity
# from newest time point.

# set sites and species numbers
sites <- rev(seq(1,20, length.out = 20))
species <- 1:30#LETTERS[1:26]
set.seed(030414)
  
base.mat <- replicate(length(species), rep(5, length(sites)))
base.mat <- base.mat + rpois(length(base.mat), 1)
dimnames(base.mat) <- list(sites, species)

# static, unchanging community
static.m.count <- base.mat
#diag(static.m.count) = 0
#static.m.count[,20:ncol(base.mat)] = base.mat[,20:ncol(base.mat)]
static.m <- prop.table(static.m.count, margin=1)

# communities that becomes more dissimilar as time goes on
linear.m.count <- base.mat
linear.m.count[lower.tri(linear.m.count)] = 0
linear.m.count[,20:ncol(base.mat)] = base.mat[,20:ncol(base.mat)]
linear.m <- prop.table(linear.m.count, margin = 1)

# community that flips  between two states
alternate.m.count <- base.mat
alternate.m.count[1:dim(alternate.m.count)[1] %% 2 == 1, 
            1:floor(dim(alternate.m.count)[2]/3)]=0
alternate.m.count[1:dim(alternate.m.count)[1] %% 2 == 0, 
            floor(dim(alternate.m.count)[2]/3)+1 :
              (floor(dim(alternate.m.count)[2]/3) *2)]=0
alternate.m.count[,20:ncol(base.mat)] = base.mat[,20:ncol(base.mat)]
alternate.m <- prop.table(alternate.m.count, margin = 1)

# stable community that diverges to another stable state, then returns
# to the original state
hump.m.count <- base.mat
hump.m.count[c(1:floor(dim(hump.m.count)[1]/3),
         (floor(dim(hump.m.count)[1]/3)*2):dim(hump.m.count)[1]),
       1:floor(dim(hump.m.count)[2]/3)]=0
hump.m.count[(floor(dim(hump.m.count)[1]/3)+1):(floor(dim(hump.m.count)[1]/3)*2-1),
       (floor(dim(hump.m.count)[2]/3)+1):(floor(dim(hump.m.count)[2]/3) *2)]=0
hump.m.count[,20:ncol(base.mat)] = base.mat[,20:ncol(base.mat)]
hump.m <- prop.table(hump.m.count, margin = 1)

# set of three stable, increasingly dissimilar communities
step.m.count <- base.mat
transitions <- floor(dim(step.m.count)[1]/3) * 1:2
step.m.count[(transitions[1]+1):transitions[2],
       (floor(dim(step.m.count)[2]/3)+1):(floor(dim(step.m.count)[2]/3)*2)]=0
step.m.count[(transitions[2]+1):dim(step.m.count)[1],
       (floor(dim(step.m.count)[2]/3)+1):(floor(dim(step.m.count)[2]/3)*3)]=0
#step.m.count[,20:ncol(base.mat)] = base.mat[,20:ncol(base.mat)]
step.m <- prop.table(step.m.count, margin = 1)

# single stable community state with two rapid, very different communities
swing.m.count <- base.mat
transitions <- floor(dim(swing.m.count)[1]/2)
swing.m.count[transitions,
          1:floor(dim(swing.m.count)[2]/2.5)] = 0
swing.m.count[transitions+1,
          (floor(dim(swing.m.count)[2]/2.5)+1):(floor(dim(swing.m.count)[2]/2.5)*2)] = 0
swing.m.count[,25:ncol(base.mat)] = base.mat[,25:ncol(base.mat)]
swing.m <- prop.table(swing.m.count, margin = 1)

mat.list <- list(static.m, hump.m, linear.m, alternate.m, step.m, swing.m)

comb.mat <- do.call("rbind", mat.list)
#comb.mat <- comb.mat[nrow(comb.mat):1,]
  
rownames(comb.mat) = nrow(comb.mat):1

# mat.list <- lapply(mat.list, function(x){
#   x <- x[nrow(x):1,]
#   return(x)}
# )

mat.list[[length(mat.list)+1]] <- comb.mat

gam.list <- lapply(mat.list, function(mat){

identify.novel.gam(mat, 0.05, "jaccard", 1, plot=FALSE, plot.data=TRUE)
  
})

# PLOT

cumul.col <- rgb(95/255,166/255,195/255)

pdf("./plots/synthetic test.pdf", height = 5, width = 10, useDingbats = FALSE)
par(mfrow=c(3,6), mar=c(0,0,0,0), oma=c(3,6,3,1), ps=10, las=1, tcl=-0.25)

# from time 1
sapply(1:(length(mat.list)-1), function(n){
  
  temp.mat <- mat.list[[n]]
#  temp.mat <- temp.mat[nrow(temp.mat):1,]
  
  plot(x=NULL, y=NULL, xlim=rev(c(1,nrow(temp.mat))), ylim=c(-0.005,1.005),
       yaxs="i", axes=FALSE)
  box()
  
  if(n==1){axis(side=2)
    mtext(side=2, las=0, line=2.5,
          text="Dissimilarity from\ntime series start", cex=0.9)
    } else {axis(side=2, labels=NA)}
  
  # from 1 dissims
  temp.dist <- as.matrix(vegdist(temp.mat, method="jaccard"))
  
  y.vals <- temp.dist[,1]
  
  lines(y=y.vals[-1], x=as.numeric(rownames(temp.mat))[-1])
  axis(side=1, mgp=c(3,0.2,0), at=c(20,15,10,5,1), labels=NA)
  
  mtext(side=3, line=0.1,
        text=c("Stable state", "Shift & return", "Slow turnover", 
               "Two state oscillation", "Progressive shift", "Two shifts & return")[n],
        font=2)
  
  text(x = relative.axis.point(0.02, "x"),
       y = relative.axis.point(0.935, "y"),
       labels = paste0("(", LETTERS[n], ")"), font=2, adj=0)
  
})

#identify
sapply(1:(length(mat.list)-1), function(n){
  
  temp.mat <- mat.list[[n]]
  temp.diss <- gam.list[[n]]
  
  plot(x=NULL, y=NULL, xlim=rev(c(1,nrow(temp.mat))), ylim=c(-0.005,1.05),
       yaxs="i", axes=FALSE, xlab="", ylab="")
  box()
  
  if(n==1){axis(side=2, at=seq(0,1,0.2), labels=c(0,0.2,0.4,0.6,0.8,""))
    mtext(side=2, las=0, line=2.5,
          text="Instantanous\ndissimilarity test", cex=0.9)
  } else {axis(side=2, labels=NA)}
  
  bins = as.numeric(rownames(temp.mat))
  
  polygon(x=c(bins, rev(bins)),
          y=c(temp.diss[[2]][,1], rev(temp.diss[[2]][,2])), 
          col="grey80", border=NA)
  lines(temp.diss[[2]][,4] ~ bins, col="grey50", lwd=2)
  
  lines(y=temp.diss[[1]]$seq.dist, x=bins)
  
  points(y=temp.diss[[1]]$seq.dist[temp.diss[[1]]$seq.p< 0.05],
         x=bins[temp.diss[[1]]$seq.p < 0.05], pch=21, bg="red",
         cex = ifelse(temp.diss[[1]]$novel[temp.diss[[1]]$seq.p< 0.05], 2, 1))
  
  axis(side=1, mgp=c(3,0.2,0), at=c(20,15,10,5,1), labels=NA)
  
  text(x = relative.axis.point(0.02, "x"),
       y = relative.axis.point(0.935, "y"),
       labels = paste0("(", LETTERS[(length(mat.list)-1) + n], ")"), font=2, adj=0)
  
})

sapply(1:(length(mat.list)-1), function(n){
  
  temp.mat <- mat.list[[n]]
  temp.diss <- gam.list[[n]]
  
  head(temp.diss)
  
  plot(x=NULL, y=NULL, xlim=rev(c(1,nrow(temp.mat))), ylim=c(-0.005,1.05),
       yaxs="i", axes=FALSE)
  box()
  
  if(n==1){axis(side=2, at=seq(0,1,0.2), labels=c(0,0.2,0.4,0.6,0.8,""))
    mtext(side=2, las=0, line=2.5,
          text="Cumulative\ndissimilarity test", cex=0.9)
  } else {axis(side=2, labels=NA)}
  
  bins = as.numeric(rownames(temp.mat))
  
  polygon(x=c(bins, rev(bins)),
          y=c(temp.diss[[3]][,1], rev(temp.diss[[3]][,2])), 
          col="grey80", border=NA)
  lines(temp.diss[[3]][,4] ~ bins, col="grey50", lwd=2)
  
  lines(y=temp.diss[[1]]$raw.min.dist, x=bins)
  
  points(y=temp.diss[[1]]$raw.min.dist[temp.diss[[1]]$min.p< 0.05],
         x=bins[temp.diss[[1]]$min.p < 0.05], pch=21, bg=cumul.col,
         cex = ifelse(temp.diss[[1]]$novel[temp.diss[[1]]$min.p< 0.05], 2, 1))
  
  axis(side=1, mgp=c(3,0.2,0), at=c(20,15,10,5,1), labels=rev(c(20,15,10,5,1)))
  
  if(n==3){mtext(side=1, at=par("usr")[2], line=1.5,
                 text="Time point", cex=0.9)}
  
  text(x = relative.axis.point(0.02, "x"),
       y = relative.axis.point(0.935, "y"),
       labels = paste0("(", LETTERS[2*(length(mat.list)-1) + n], ")"), font=2, adj=0)
  
})

dev.off()

# now the combination plot 

pdf("./plots/synthetic test (combination).pdf", height = 4, width = 6.5, useDingbats = FALSE)
par(mfcol=c(3,1), mar=c(0,0,0,0), oma=c(3,6,1,1), ps=10, las=1, tcl=-0.25)

n <- length(mat.list)
# from time 1
 temp.mat <- mat.list[[n]]
  #  temp.mat <- temp.mat[nrow(temp.mat):1,]
  
  plot(x=NULL, y=NULL, xlim=rev(c(1,nrow(temp.mat))), ylim=c(-0.01,1.005),
       yaxs="i", axes=FALSE)
  box()
  
axis(side=2)
    mtext(side=2, las=0, line=2.5,
          text="Dissimilarity from\ntime series start", cex=0.9)

  # from 1 dissims
  temp.dist <- as.matrix(vegdist(temp.mat, method="jaccard"))
  
  y.vals <- temp.dist[,1]
  
  lines(y=y.vals[-1], x=as.numeric(rownames(temp.mat))[-1])
  axis(side=1, mgp=c(3,0.2,0), at=c(1,seq(20,120,20)), labels=NA)
  

#identify

  temp.mat <- mat.list[[n]]
  temp.diss <- gam.list[[n]]
  
  head(temp.diss)
  
  plot(x=NULL, y=NULL, xlim=rev(c(1,nrow(temp.mat))), ylim=c(-0.01,1.005),
       yaxs="i", axes=FALSE)
  box()
  
  axis(side=2, at=seq(0,1,0.2), labels=c(0,0.2,0.4,0.6,0.8,""))
    mtext(side=2, las=0, line=2.5,
          text="Instantaneous\ndissimilarity test", cex=0.9)
  
  bins = as.numeric(rownames(temp.mat))
  
  polygon(x=c(bins, rev(bins)),
          y=c(temp.diss[[2]][,1], rev(temp.diss[[2]][,2])), 
          col="grey80", border=NA)
  lines(temp.diss[[2]][,4] ~ bins, col="grey50", lwd=2)
  
  lines(y=temp.diss[[1]]$seq.dist, x=bins)
  
  points(y=temp.diss[[1]]$seq.dist[temp.diss[[1]]$seq.p< 0.05],
         x=bins[temp.diss[[1]]$seq.p < 0.05], pch=21, bg="red",
         cex = ifelse(temp.diss[[1]]$novel[temp.diss[[1]]$seq.p< 0.05], 2, 1))
  
  axis(side=1, mgp=c(3,0.2,0), at=c(1,seq(20,120,20)), labels=NA)

  temp.mat <- mat.list[[n]]
  temp.diss <- gam.list[[n]]
  
  head(temp.diss)
  
  plot(x=NULL, y=NULL, xlim=rev(c(1,nrow(temp.mat))), ylim=c(-0.01,1.005),
       yaxs="i", axes=FALSE)
  box()
  
 axis(side=2, at=seq(0,1,0.2), labels=c(0,0.2,0.4,0.6,0.8,""))
    mtext(side=2, las=0, line=2.5,
          text="Cumulative\ndissimilarity test", cex=0.9)

  bins = as.numeric(rownames(temp.mat))
  
  polygon(x=c(bins, rev(bins)),
          y=c(temp.diss[[3]][,1], rev(temp.diss[[3]][,2])), 
          col="grey80", border=NA)
  lines(temp.diss[[3]][,4] ~ bins, col="grey50", lwd=2)
  
  lines(y=temp.diss[[1]]$raw.min.dist, x=bins)
  
  points(y=temp.diss[[1]]$raw.min.dist[temp.diss[[1]]$min.p< 0.05],
         x=bins[temp.diss[[1]]$min.p < 0.05], pch=21, bg=cumul.col,
         cex = ifelse(temp.diss[[1]]$novel[temp.diss[[1]]$min.p< 0.05], 2, 1))
  
  axis(side=1, at=c(1,seq(20,120,20)), labels=rev(c(1,seq(20,120,20))), mgp=c(3,0.2,0))
  
  mtext(side=1, line=1.5, text="Time point", cex=0.9)

dev.off()


}