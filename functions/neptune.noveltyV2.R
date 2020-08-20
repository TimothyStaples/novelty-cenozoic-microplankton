neptune.novelty <- function(dataset,
                            bin.width,
                            ssmat.type = "pa",
                            taxon.res,
                            bin.cutoff,
                            taxa.cutoff,
                            novel.alpha,
                            novel.metric,
                            longhurst=FALSE,
                            rich.cutoff=0,
                            gam.max.k = -1,
                            plot = TRUE){
  
  data.c <- dataset[dataset$Age..Ma..Gradstein.et.al..2012 <= 66, ]
  data.c <- data.c[!is.na(data.c$Site), ]
  
  # subset to just data with abundance information
  if(ssmat.type == "abund"){
  data.c$abund <- as.numeric(as.character(data.c$Taxon.Abundance))
  data.c <- data.c[!is.na(data.c$abund) & data.c$abund > 0, ]
  }

  # simplify some column names
  colnames(data.c)[colnames(data.c) == "Age..Ma..Gradstein.et.al..2012"] = "age"
  colnames(data.c)[colnames(data.c) == "Site"] = "site"
  colnames(data.c)[colnames(data.c) == "Latitude"] = "lat"
  colnames(data.c)[colnames(data.c) == "Longitude"] = "long"
  
  if(longhurst){
    print("Aggregating samples into Longhurst provinces")
    require(maptools)
    require(sp)
    
    longhurst <-readShapeSpatial("./raw.datafiles/longhurst_v4_2010/Longhurst_world_v4_2010.shp")
  
    colnames(data.c)[colnames(data.c) == "site"] = "Neptune.site"
    
    data.coords <- data.c[!duplicated(data.c$Neptune.site),
                          c("long","lat", "Neptune.site")]
    coordinates(data.coords)<-c("long","lat")
    
    data.coords@data <- cbind(data.coords@data,
                              over(data.coords, longhurst))
    
    colnames(data.coords@data) <- c("Neptune.site", "site", "site.name")
    
    data.c <- merge(data.c, data.coords@data,
                    by.x="Neptune.site", by.y="Neptune.site",
                    all.x=TRUE, all.y=FALSE)
    
  }
  
  print("Testing sample completeness")
  cand.samp <- test_completeness(dataset = data.c,
                                 tax = taxon.res,
                                 age = "age",
                                 bin.width = bin.width)
  
  raw.samp <- cand.samp[[2]]
  cand.samp <- cand.samp[[1]]
  
  print("Creating site-species matrices...")
  # create site-species matrices using bin width
  bins<-seq(0,66,bin.width)
  
  cand.ssmats <- lapply(unique(data.c$site), function(site){
    print(site)
    temp <- create.ssmat(records = data.c,
                 site.name = site,
                 bins = bins,
                 tax.res = taxon.res,
                 type=ssmat.type)
    
    if(class(temp)=="matrix"){
  # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
    }
    
    return(temp)
    
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
      identify.novel.gam(site.sp.mat = cand.ssmats[[site]], 
                         alpha = novel.alpha, 
                         metric=novel.metric,
                         site=site, 
                         plot = plot,
                         gam.max.k = gam.max.k)
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
              novel = cand.novel,
              samp_comp = cand.samp,
              raw_samp = raw.samp))
  
}
