create.ssmat <- function(records,
site.name, bins, age.limits = NA, tax.res, 
type = "pa" # options are "pa" (presence absence),
            #             "relpa (presence absence relativized by sampling effort)
            #             "abund" (relative abundance on subset of sites with abundance info)
){
  
  site.df <- droplevels(records[records$site == site.name &
                                !is.na(records$age), ])
  
  if(dim(site.df)[1]==0){return(NULL)}
  
  if(!is.na(age.limits[1])){
    
    site.df <- site.df[site.df$age <= age.limits[2] &
                         site.df$age >= age.limits[1], ]
    
  }
  
  # for each bin width
  site.df$age.bin <- cut(site.df$age, breaks = bins)
  
  age.bins <- levels(site.df$age.bin)
  mid.points <- cbind(age.bins,
                      as.numeric(substr(age.bins, regexpr(",", age.bins)+1,
                                        nchar(age.bins)-1)) - 0.5*(bins[2]-bins[1]))
  
  site.df$mid.point <- mid.points[match(site.df$age.bin, age.bins), 2]
  
  # filter out single observation sites
  if(length(unique(site.df$age.bin))<=2){return(NULL)}
  
  # cut out samples without species IDs & make site-species matrix
  if(tax.res == "species"){
    site.sp <- droplevels(site.df[!is.na(site.df$species), ])
    site.sp$binom <- paste(site.sp$genus, site.sp$species)
    site.sp.mat <- as.matrix(table(site.sp$age.bin, site.sp$binom))
    rownames(site.sp.mat) <- mid.points[match(rownames(site.sp.mat),
                                              age.bins), 2]
    } else {
    
    site.sp <- droplevels(site.df[!is.na(site.df[,tax.res]), ])
    site.sp.mat <- as.matrix(table(site.sp$age.bin, site.sp[,tax.res]))
    site.sp.mat <- site.sp.mat[,colSums(site.sp.mat)>0]
    rownames(site.sp.mat) <- mid.points[match(rownames(site.sp.mat),
                                              age.bins), 2]
    
    }
  
  if(type=="pa"){
    
    site.sp.mat <- ifelse(site.sp.mat > 0, 1, 0)
    
  }
  
  if(type=="relpa"){
    
    temp <- site.df$age[site.df$mid.point== site.df$mid.point[1]]

    age.bin.rep <- tapply(site.df$age, site.df$mid.point, 
                          function(x){length(unique(x))})
    
    if(sum(names(age.bin.rep) %in% rownames(site.sp.mat)) != dim(site.sp.mat)[1]){
      print(paste0("ERROR in relativizing PA data in site",
                   site.df$site[1]))
      return(NULL)}
    
    site.sp.mat <- site.sp.mat / as.vector(age.bin.rep[match(rownames(site.sp.mat),
                                                             names(age.bin.rep))])
    
  }
  
  if(type == "abund"){
    
    site.df$abund <- as.numeric(as.character(site.df$Taxon.Abundance))
    
    abund.site <- droplevels(site.df[!is.na(site.df$abund), ])
  
    if(dim(abund.site)[1]==0){return(NULL)}  
    
    abund.site$age.mid <- mid.points[match(abund.site$age.bin,
                                           age.bins),2]
    
    site.sp.mat <- matrix(0, 
                          nrow=length(unique(abund.site$age.bin)),
                          ncol=length(unique(abund.site$Taxon.ID)))
    rownames(site.sp.mat) <- sort(unique(as.numeric(abund.site$age.mid)))
    colnames(site.sp.mat) <- unique(abund.site$Taxon.ID)
    
    for(tax in unique(abund.site$Taxon.ID)){

    temp <- abund.site[abund.site$Taxon.ID==tax,] 

    tax.mean <- tapply(temp$abund, temp$age.mid, mean)
    
    site.sp.mat[names(tax.mean),as.character(tax)] = tax.mean
      
    }
    
    site.sp.mat <- prop.table(site.sp.mat, 1)
    
  }
    
  site.sp.mat <- site.sp.mat[, colSums(site.sp.mat) > 0]
  
  return(as.matrix(unclass(site.sp.mat)))
  
}