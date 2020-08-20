cut.complete <- function(novel.list,
                         cutoff = 0.05,
                         plot.name){
  
  # aggregate sample completeness
  samp_comp_df <- do.call("rbind", lapply(1:4, function(n){
    
    test <- novel.list[[n]]$samp_comp
    test$taxa <- names(novel.list)[n]
    return(test)
    
  }))
  
  # set thresholds to 5% quantile
  thresholds <- tapply(samp_comp_df$samp.prob,
                       samp_comp_df$taxa,
                       function(x){quantile(x, 
                                            probs=cutoff, 
                                            na.rm=TRUE)})
  
  samp_comp_df$site.taxa <- paste0(samp_comp_df$site, ".",
                                   samp_comp_df$taxa)
  
  samp_comp_df$threshold <- thresholds[match(samp_comp_df$taxa,
                                             names(thresholds))]
  
  sites.to.remove <- samp_comp_df$site.taxa[samp_comp_df$samp.prob <
                                            samp_comp_df$threshold |
                                       is.na(samp_comp_df$samp.prob)]
  
  site.to.remove <- data.frame(site = substr(sites.to.remove, 1,
                                 regexpr("\\.", sites.to.remove)-1),
                           taxa = substr(sites.to.remove,
                                 regexpr("\\.", sites.to.remove)+1,
                                 nchar(sites.to.remove)))
  
  comp.novel.list <- lapply(1:length(novel.list), function(n){
    
    x <- novel.list[[n]]
    
    sub.remove <- site.to.remove[site.to.remove$taxa == names(novel.list)[n],]
    
    print(paste0(sum(!sub.remove$site %in% names(x$novel)),
                 " out of ",
                 nrow(sub.remove),
                 " already removed"))
    
    x$novel <- x$novel[which(!names(x$novel) %in%
                               as.character(sub.remove$site))]
    x$sites <- names(x$novel)
    
    return(x)
    
  })
  names(comp.novel.list) <- names(novel.list)
  
  samp_comp_df$taxa <- factor(samp_comp_df$taxa,
                              levels=c("nano", "foram", "radio", "diatom"))
  
  pdf(paste0("./plots/", plot.name), height=4, 
      width=4, useDingbats = FALSE)
  par(mar=c(2,3,1,1), ps=8, mgp=c(3,0.5,0), tcl=-0.25)
  plot(x=NULL, y=NULL, xlim=c(0.5,4.5), ylim=c(0,1), axes=FALSE)
  
  axis(side=2, at=seq(0,1,0.2), las=1)
  axis(side=2, at=seq(0,1,0.1), labels=NA, tcl=-0.125)
  mtext(side=2, line=1.75, text="Sample completeness")
  
  axis(side=1, at=c(1:4), 
       labels=c("Nanno", "Foram", "Radio", "Diatom"),
       mgp=c(3,0.2,0))
  
  segments(x0 = (1:4) - 0.5,
           x1 = (1:4) + 0.5,
           y0 = thresholds,
           y1 = thresholds,
           lwd=2.5, col="red")
  
  with(samp_comp_df[samp_comp_df$samp.prob <
                    samp_comp_df$threshold,],
  points(y = samp.prob,
         x = as.numeric(as.factor(taxa)),
         pch=16, col=rgb(0.5,0.5,0.5,0.25)))
  
  with(samp_comp_df[samp_comp_df$samp.prob >=
                      samp_comp_df$threshold,],
       points(y = samp.prob,
              x = jitter(as.numeric(as.factor(taxa)), amount=0.1),
              pch=16, col=rgb(0.5,0.5,0.5,0.25)))
  
  boxplot(samp_comp_df$samp.prob ~ samp_comp_df$taxa, add=TRUE,
          xlab="", ylab="", axes=FALSE)
  
  box()
  dev.off()
  
  return(comp.novel.list)
}
