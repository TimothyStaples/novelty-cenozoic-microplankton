local.extorig.neotoma <- function(ssmat,
                          novel,
                          site.name,
                          raw.data){
  
  
  # Check to see if age of site-species matrices increase
  # (first row is oldest)
  if(as.numeric(as.character(rownames(ssmat)[1])) < 
     as.numeric(as.character(rownames(ssmat)[dim(ssmat)[1]]))){
    
    ssmat <- ssmat[dim(ssmat)[1]:1, ]
    
  }

# create presence-absence matrix to make it easier to identify
# local extinctions and originations
bin.mat <- ifelse(ssmat > 0, 1, 0)

# create data-frame of bins to hold results
ext.df <- data.frame(bins = rownames(ssmat),
                     site = site.name,
                     n.from.start = 1:nrow(ssmat),
                     n.from.end = nrow(ssmat):1)

ext.df$bin.lag <- c(NA, abs(diff(as.numeric(as.character(ext.df$bins)))))

# cat = community category
ext.df$cat.aft <- novel$cat[match(as.character(ext.df$bins),
                                  as.character(novel$bins))]

ext.df$cat.bef <- novel$cat.bef[match(as.character(ext.df$bins),
                                  as.character(novel$bins))]


# now we run through each pair of communities and look for local extinctions
# and originations.
count.df <- do.call("rbind", lapply(1:dim(bin.mat)[1], function(bin){
  
  # if we have the edges (1st bin for causes, last bin for consequences)
  # return NAs to fit our "ext.df" object
  if(bin == 1){return(NA)}
  
  # if it's not an edge, we need to look forward in time for originations
  # and backwards in time for extinctions. We can identify our target community, 
  # either looking for the last presence (extinction) or first presence (origination),
  # based on whether we are interested in what happened before the identity community,
  # or after.
    target = c(ext = bin - 1,
               orig = bin)
  
  # for extinctions, get everything from target to end of time-series
  ext.mat <- bin.mat[target["ext"]:dim(bin.mat)[1],]

  # for originations, get everything from 1st entry to target
  orig.mat <- bin.mat[1:target["orig"],]
  
  # remove species that don't appear in the target
  ext.mat <- as.matrix(ext.mat[, bin.mat[target["ext"],]>0])
  orig.mat <- as.matrix(orig.mat[, bin.mat[target["orig"],]>0])
  
  total.n <- c(ext = dim(ext.mat)[2], orig = dim(orig.mat)[2])
  
  # now we need to look at each species present in the target community
  
  # originations are easy, we just need to find where the first 1 in each column is,
  # and see if it lines up with our target community
  orig.bin <- apply(orig.mat, 2, function(x){which(x==1)[1]}) %in% target["orig"]
  
  # immigrations are all 0 -> 1 across our two communities, minus the ones that are
  # actual local originations
  immig.bin <- apply(orig.mat, 2, function(x){x[target["orig"]] == 1 &
                                              x[target["orig"]-1] == 0})
  
  # extinctions are only a little harder. Here we look for the first 1 in the 
  # column, but reverse the column first, so we are looking for the first 1 going
  # back in time. This is the last presence of the species before it went locally
  # extinct, and then we check to see if this matches our target community by
  # re-reversing our positions, then adjusting them so they compare row numbers in
  # the overall matrix by adding the extintction target community.
  ext.bin <- (target["ext"] + dim(ext.mat)[1] - 
             apply(ext.mat, 2, function(x){which(rev(x) == 1)[1]})) %in% target["ext"]
  
  # emigrations are the same as immigations, except we're looking for 1 -> 0 patterns, and
  # our matrix starts at our target, so we can use numeric constants to pull out the preceding
  # and succeeding communities.
  emig.bin <- apply(ext.mat, 2, function(x){x[2] == 0 & x[1] == 1})
  
  return(data.frame(ext = sum(ext.bin),
                    orig = sum(orig.bin),
                    immig = sum(immig.bin) - sum(orig.bin),
                    emig = sum(emig.bin) - sum(ext.bin),
                    ext.rich = total.n[1],
                    orig.rich = total.n[2],
                    loss = sum(emig.bin),
                    gain = sum(immig.bin)))
             
  }))

ext.df <- cbind(ext.df, count.df)


return(ext.df)

}