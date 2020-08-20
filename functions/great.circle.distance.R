# convert 2 long-lat pairs into a distance in km
gc.dist<-function(lat1, lon1, lat2, lon2){
  
  deg2rad <- function(deg) {(deg * pi) / (180)} # custom function to turn degrees to radians
  
  R<-6371 # Radius of the earth in km
  dLat<-deg2rad(lat2-lat1) # deg2rad below
  dLon<-deg2rad(lon2-lon1); 
  a<-sin(dLat/2) * sin(dLat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dLon/2) * sin(dLon/2)
  
  c = 2 * atan2(as.numeric(sqrt(a)), as.numeric(sqrt(1-a)))
  d = R * c # Distance in km
  return(d)
}