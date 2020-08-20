test_completeness <- function(dataset, 
                              tax = "Taxon.ID",
                              age = "Age..Ma..Gradstein.et.al..2012",
                              bin.width = 0.1){
  
  dataset$age <- dataset[, age]
  dataset <- dataset[dataset$age <= 66, ]
  dataset <- droplevels(dataset[!is.na(dataset$site), ])
  
  bins<-seq(0,66,bin.width)
  
  # for each bin width
  dataset$age.bin <- cut(dataset$age, breaks = bins)
  
  age.bins <- levels(dataset$age.bin)
  mid.points <- cbind(age.bins,
                      as.numeric(substr(age.bins, regexpr(",", age.bins)+1,
                                        nchar(age.bins)-1)) - 0.5*(bins[2]-bins[1]))
  
  dataset$mid.point <- mid.points[match(dataset$age.bin, age.bins), 2]
  dataset$mid.point <- as.numeric(dataset$mid.point)

  # by site
  dd.list <- lapply(split(dataset, f=dataset$site), function(x){
  
    print(as.character(x$site[1]))  
    
    if(length(unique(x$mid.point)) < 3){
      return(data.frame(site = x$site[1],
                        samp.prob = NA))
    }
    
    temp.dd <- divDyn(x, tax = tax, bin = "mid.point", revtime = TRUE)
    
    agg.completeness <- data.frame(site = x$site[1],
               samp.prob = sum(temp.dd$t3, na.rm=T) / (sum(temp.dd$t3, na.rm=T) + 
                                                      sum(temp.dd$tPart, na.rm=T)))
    
    return(list(temp.dd, agg.completeness))
    
  })
  
  site.dd <- do.call("rbind", lapply(dd.list, function(x){x[[2]]}))
  raw.dd <- lapply(dd.list, function(x){x[[1]]})
  
  return(list(site.dd = site.dd,
              raw.dd = raw.dd))
  
  }