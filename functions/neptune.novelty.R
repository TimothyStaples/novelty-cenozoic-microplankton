neptune.novelty <- function(dataset,
                            bin.width,
                            ssmat.type = "pa",
                            taxon.res,
                            bin.cutoff,
                            taxa.cutoff,
                            novel.alpha,
                            novel.metric){
  
  # Exclude polar latitudes (is this necessary now we're treating
  # each time series separately?)
  #data.c <- subset(dataset, abs(Estimated.Paleo.Latitude)<60) 
  
  data.c <- dataset[dataset$Age..Ma..Gradstein.et.al..2012 <= 65, ]
  data.c <- data.c[!is.na(data.c$Site), ]
  
  # subset to just data with abundance information
  if(ssmat.type == "abund"){
  data.c$abund <- as.numeric(as.character(data.c$Taxon.Abundance))
  data.c <- data.c[!is.na(data.c$abund) & data.c$abund > 0, ]
  }
  
  # # get a per-site count of the number of unique age samples (slices?)
  # site.reps <- sapply(unique(data.c$Site), function(x){
  #   
  #   temp.data <- data.c[data.c$Site == x, ]
  #   sum(!duplicated(temp.data$Age..Ma..Gradstein.et.al..2012))
  #   
  # })
  # 
  # names(site.reps) <- unique(data.c$Site)
  
  # simplify some column names
  colnames(data.c)[colnames(data.c) == "Age..Ma..Gradstein.et.al..2012"] = "age"
  colnames(data.c)[colnames(data.c) == "Site"] = "site"
  colnames(data.c)[colnames(data.c) == "Latitude"] = "lat"
  colnames(data.c)[colnames(data.c) == "Longitude"] = "long"
  
  print("Creating site-species matrices...")
  # create site-species matrices using bin width
  bins<-seq(0,65,bin.width)
  
  cand.ssmats <- lapply(unique(data.c$site), function(site){
    
    create.ssmat(records = data.c,
                 site.name = site,
                 bins = bins,
                 tax.res = taxon.res,
                 type=ssmat.type)
  })
  names(cand.ssmats) <- unique(data.c$site)
  
  # remove sites that had no appropriate data
  cand.ssmats <- cand.ssmats[sapply(cand.ssmats, class)=="matrix"]
  
  # remove sites with fewer age bins or taxa than cutoff thresholds
  ssmat.dims <- t(sapply(cand.ssmats, dim))
  cand.ssmats <- cand.ssmats[ssmat.dims[,1] >= bin.cutoff & 
                             ssmat.dims[,2] >= taxa.cutoff]
  
  # use the identify.novel.gam function to find novel communities in each time-series.
  # function creates plots as it runs.
  print("Identifying novel communities...")
  cand.novel <- lapply(names(cand.ssmats), function(site){
    print(site)
      identify.novel.gam(cand.ssmats[[site]], 
                         novel.alpha, metric=novel.metric,
                         site=site)
  })
  names(cand.novel) <- names(cand.ssmats)
  
  # remove communities that couldn't be estimated due to no variation
  # in taxa
  cand.subset <- sapply(cand.novel, class) == "data.frame"
  cand.ssmats <- cand.ssmats[cand.subset]
  cand.novel <- cand.novel[cand.subset]
  
  return(list(data = data.c[data.c$site %in% names(cand.ssmats), ],
              sites = names(cand.ssmats),
              ssmats = cand.ssmats,
              novel = cand.novel))
  
}
