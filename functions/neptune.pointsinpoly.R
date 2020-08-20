pointsinpoly <- function(data,
                         shape,
                         shape.data){
  
  coords <- data[!duplicated(data$Site),
                 c("Site", "Longitude", "Latitude")]
  coordinates(coords) <- c("Longitude", "Latitude")

  poly.pos <- coords %over% shape
  
  data.frame(Site=coords$Site, 
             ocean=as.character(poly.pos[,shape.data]))
}