arrow.rectangle <- function(x,y,w,l,nsteps=100, shape="rect",
                        screen=NA, arrow.consts, add.arrow=FALSE, ...){
  
  # set aspect ratio to ensure symmetrical shapes
  plot.window<-dev.size(units="cm")
  
  if(!is.na(screen)[1]){
    plot.window <- dev.size(units="cm") * c(screen[2]-screen[1],
                                            screen[4]-screen[3])
    
  }
  
  x.axis<-par("usr")[1:2]
  y.axis<-par("usr")[3:4]
  
  # transform each point into cm scores, using the origin of the plot as our
  # 0cm point.
  x<-(x-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1]
  y<-(y-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2]
  w <-(w-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2]
  l <-(l-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1]
  
  xc <- c(x,   x-l, x-l, x,   x+l, x+l)
  yc <- c(y+w, y+w, y-w, y-w, y-w, y+w)

  shape.top <- rbind(cbind(xc,yc)[yc == max(yc) & xc==x,],
                     cbind(xc,yc)[yc == min(yc) & xc==x,])
  
#    if(add.arrow){
    
    arrow.points <- rbind(shape.top[1,] + c(0, - 0.2*arrow.consts[2]),
                          shape.top[1,] + c(0.1*arrow.consts[1],  - 0.2*arrow.consts[2]),
                          shape.top[1,] + c(0.1*arrow.consts[1],  - 0.1*arrow.consts[2]),
                          shape.top[1,] + c(0.3*arrow.consts[1],  - 0.5*arrow.consts[2]),
                          shape.top[2,] + c(0.1*arrow.consts[1],  + 0.1*arrow.consts[2]),
                          shape.top[2,] + c(0.1*arrow.consts[1],  + 0.2*arrow.consts[2]),
                          shape.top[2,] + c(0, + 0.2*arrow.consts[2]))
    
    xarrow <- c(xc[which(xc == shape.top[1,1] & yc==shape.top[1,2])],
                arrow.points[,1],
                xc[which(xc == shape.top[2,1] & yc==shape.top[2,2]):
                     which(xc == shape.top[1,1] & yc==shape.top[1,2])])
    
    yarrow <- c(yc[which(xc == shape.top[1,1] & yc==shape.top[1,2])],
                arrow.points[,2],
                yc[which(xc == shape.top[2,1] & yc==shape.top[2,2]):
                     which(xc == shape.top[1,1] & yc==shape.top[1,2])])
    
    xc <- xarrow
    yc <- yarrow
  
# }  
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...)
  # return(cbind(xc,yc))
}