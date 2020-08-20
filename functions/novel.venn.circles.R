novel.venn.circles <- function(y.pos, x.pos, radius){
  
  I<-draw.circle(x=y.pos, y=x.pos[1], radius=radius[1], nv=500, 
                 border="black", col="white", lwd=1)
  C<-draw.circle(x=y.pos, y=x.pos[2], radius=radius[2], nv=500, 
                 border="black", col="white", lwd=1)
  
  # turn I circle to spatial polygons to get overlap
  I.sp <- Polygon(cbind(I$x, I$y))
  I.sp <- SpatialPolygons(list(Polygons(list(I.sp), ID = "a")))
  
  C.sp <- Polygon(cbind(C$x, C$y))
  C.sp <- SpatialPolygons(list(Polygons(list(C.sp), ID = "a")))
  
  N.sp <- raster::intersect(I.sp, C.sp)
  I.sub <- gDifference(I.sp, C.sp)
  C.sub <- gDifference(C.sp, I.sp)
  
  plot(I.sub, col="red", add=TRUE, border="black", lwd=1)
  plot(C.sub, col=rgb(cumul.col[1], cumul.col[2], cumul.col[3]),
       add=TRUE, border="black", lwd=1)
  plot(N.sp, col="orange", add=TRUE, border="black", lwd=1)
  
}