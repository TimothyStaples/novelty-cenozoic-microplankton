identify.novel.plotting<-function(site.sp.mat, alpha, metric){
  
  require(vegan) # community disimilarity metrics
  require(mgcv) # additive modelling
  require(arm) # invlogit transformation
  
  # check to see if site species matrix rows run oldest -> youngest
  # if so, rotate site-species matrix to get oldest sites first
  if(as.numeric(as.character(rownames(site.sp.mat)[1])) < 
     as.numeric(as.character(rownames(site.sp.mat)[dim(site.sp.mat)[1]]))){
    
    site.sp.mat <- site.sp.mat[dim(site.sp.mat)[1]:1, ]
    
  }
  
  # next get distance matrix for observed communities
  site.prop.mat <- prop.table(site.sp.mat, 1)
  
  # calculate dissimilarity matrix, splitting those that can be
  # incorporated into vegdist from those that can't
  if(metric %in% c("chord", "SQchord")){
    
    require(analogue)
    
    site.dist <- as.matrix(distance(site.prop.mat,
                                    method=metric))
    
    if(metric=="chord"){site.dist <- site.dist / sqrt(2)}
    if(metric=="SQchord"){site.dist <- site.dist / 2}
    
  }
  
  if(metric == "hellinger"){
    
    hell.mat <- decostand(site.prop.mat, method="hellinger")  
    site.dist <- as.matrix(dist(hell.mat))
    site.dist <- site.dist / sqrt(2)
    
  }
  
  if(!metric %in% c("chord", "SQchord", "hellinger")){
    
    site.dist <- as.matrix(vegdist(site.prop.mat,
                                   method=metric))
  }
  
  # convert bin names to numbers to use as continuous variable
  bins <- as.numeric(rownames(site.sp.mat))
  
  # This sets the maximum knot number for the additive model. One of the biggest
  # problems with additive models is overfitting. Ideally we want penalized 
  # splines, where knot number is weighted against variance explained. This will
  # provide a good measure of local mean trends in dissimilarity. But if we want
  # this method to run on smaller time-series, setting penalized splines doesn't work
  # very well. So I've set this to set maximum knots to half the number of bins, rounded
  # down, for time-series with between 4 and 10 bins. Less than 4 will not run using
  # splines, and is set to run as a linear regression below.
  set.k <- ifelse(dim(site.sp.mat)[1] > 10,
                  -1,
                  ifelse(dim(site.sp.mat)[1] > 4,
                         floor(dim(site.sp.mat)[1])/2,
                         NA))
  
  # SEQUENTIAL DISIMILARITY ####
  
  # calculate differenced sequential dissimilarities between 
  # time t and t-1
  seq.dist <- c(NA, diag(site.dist[-1,-dim(site.dist)[2]]))
  
  # transform to remove 0s and 1s for beta regression
  seq.dist.tr <- (seq.dist * (length(seq.dist)-1) + 0.5) / length(seq.dist)
  
  # model localised trend over time and extract residuals to get dissimilarity 
  # compared to local mean.
  if(!is.na(set.k)){
    seq.gam <- gam(seq.dist.tr ~ s(bins, bs="cr", k= set.k), family=betar())
  } else{
    seq.gam <- gam(seq.dist.tr ~ bins, family=betar())
  }
  
  # MINIMUM DISIMILARITY ####
  
  # calculate minimum dissimilarity from time 1 to time t (time 1 being
  # earliest time point)
  min.dist <- sapply(1:dim(site.dist)[1], function(n){
    if(n==1){return(NA)}
    min(site.dist[n,1:n][-n])
  })
  
  # transform to remove 0s and 1s for beta regression
  min.dist.tr <- (min.dist * (length(min.dist)-1) + 0.5) / length(min.dist)
  
  # model localised trend over time 
  if(!is.na(set.k)){
    min.gam <- gam(min.dist.tr ~ s(bins, bs="cr", k= set.k), family=betar())
  } else{
    min.gam <- gam(min.dist.tr ~ bins, family=betar())
  }
  
  # COMPARING OBSERVED TO EXPECTED SEQUENTIAL DISIMILARITY ####
  
  # This process calculates the p-value of the observed disimilarity score
  # being part of the expected distribution at the point in the time-series.
  
  # convert mu & phi beta parameters to A & B to use in qbeta.
  # mu = the mean predicted dissimilarity at each time point based on the 
  #      additive model.
  # phi = a dispersion parameter, which is set globally across the model.
  #       phi is unfortunately hidden deep in the gam model as a text string.
  
  seq.mu <- c(NA, seq.gam$fitted.values)
  phi <- as.numeric(substr(seq.gam$family$family,
                           regexpr("\\(", seq.gam$family$family)+1,
                           nchar(seq.gam$family$family)-1))
  
  # shape parameters used in qbeta.
  A = seq.mu * phi
  B = (1 - seq.mu) * phi
  
  # predict 5% and 95% prediction intervals from beta distribution parameters. 
  # We use 95% rather than 97.5% because this is a one-tailed probability test.
  # We are not concerned about communities that are MORE similar than predictions.
  # This is done for each bin along the time-series.
  seq.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               seq.p = pbeta(seq.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # EXPECTED OBSERVED TO MINIMUM DISIMILARITY ####
  
  # convert mu & phi beta parameters to A & B to use in qbeta
  min.mu <- c(NA, min.gam$fitted.values)
  phi <- as.numeric(substr(min.gam$family$family,
                           regexpr("\\(", min.gam$family$family)+1,
                           nchar(min.gam$family$family)-1))
  A = min.mu * phi
  B = (1 - min.mu) * phi
  
  # predict 5% and 95% prediction intervals from beta distribution parameters
  min.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               min.p = pbeta(min.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # CREATE RETURN DATA-FRAME ####
  
  return.data <- data.frame(bins = rownames(site.sp.mat),
                            seq.dist = seq.dist,
                            seq.resid = c(NA, resid(seq.gam)),
                            raw.min.dist = min.dist,
                            min.resid = c(NA, resid(min.gam)),
                            seq.p = seq.p$seq.p,
                            min.p = min.p$min.p,
                            instant = seq.p$seq.p <= alpha,
                            cumul = min.p$min.p <= alpha,
                            novel = seq.p$seq.p <= alpha & min.p$min.p <= alpha,
                            seq.exp=seq.mu,
                            min.exp=min.mu)
  
  return.data$cat = "back"
  return.data[return.data$instant & 
                !is.na(return.data$instant) &
                !is.na(return.data$cumul) &
                !return.data$cumul, "cat"] = "instant"
  return.data[!return.data$instant & 
                !is.na(return.data$instant) &
                !is.na(return.data$cumul) &
                return.data$cumul, "cat"] = "cumul"
  return.data[return.data$instant & 
                !is.na(return.data$cumul) &
                !is.na(return.data$instant) &
                return.data$cumul, "cat"] = "novel"
  
    return(list(return.data,
                cbind(seq.p, c(NA, plogis(predict(seq.gam)))),
                cbind(min.p, c(NA, plogis(predict(min.gam))))))
}
