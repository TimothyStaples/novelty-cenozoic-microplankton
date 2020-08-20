arrow.shape <- function(x,y,r,nsteps=100, shape="circle",
                        screen=NA, rads=c(0,2), add.arrow=FALSE, ...){
  
  # set aspet ratio to ensure symmetrical shapes
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
  
  if(shape=="circle"){
    
    if(rads[2] < rads[1]){
      
      rs <- c(seq(rads[1]*pi, 2*pi, len=nsteps*((2-rads[1]) / (2-rads[1]+rads[2]))),
              seq(0*pi, rads[2]*pi, len=nsteps*(rads[2] / (2-rads[1]+rads[2]))))
      
      width <- 2-max(rads)+min(rads)
    } else {
      rs <- seq(rads[1]*pi,rads[2]*pi,len=nsteps)  
      width <- rads[2]-rads[1]
    }
    
    xc <- x+r*cos(rs) 
    yc <- y+r*sin(rs) 
    
    if(width < 1){
      xc <- c(xc, x)
      yc <- c(yc, y)
    }
    
  }
  
  if(shape=="triangle"){
    
    xc <- c(x, x-r, x, x+r)
    yc <- c(y+r, y-r, y-r, y-r)
    
  }
  
  if(shape=="square"){
    
    xc <- c(x,   x-r, x-r, x,   x+r, x+r)
    yc <- c(y+r, y+r, y-r, y-r, y-r, y+r)
    
    
  }
  
  if(shape=="diamond"){
    
    xc <- c(x, x-1.2*r, x, x+1.2*r)
    yc <- c(y+1.2*r, y, y-1.2*r, y)
    
  }
  
  if(add.arrow){
    
    if(shape=="circle"){
      shape.top <- rbind(cbind(xc,yc)[yc == max(yc),],
                         cbind(xc,yc)[yc == min(yc),])
      
    } else {
      shape.top <- rbind(cbind(xc,yc)[yc == max(yc) & xc==x,],
                         cbind(xc,yc)[yc == min(yc) & xc==x,])
    }
    
    arrow.dist <- c(max(xc)-min(xc), max(yc)-min(yc))
    if(shape !="circle"){
      arrow.dist <- arrow.dist*c(0.45,1)
    }
    
    arrow.points <- rbind(shape.top[1,] + c(0, - 0.4*arrow.dist[2]),
                          shape.top[1,] + c(0.2*arrow.dist[1],  - 0.4*arrow.dist[2]),
                          shape.top[1,] + c(0.2*arrow.dist[1],  - 0.3*arrow.dist[2]),
                          shape.top[1,] + c(0.6*arrow.dist[1],  - 0.5*arrow.dist[2]),
                          shape.top[2,] + c(0.2*arrow.dist[1],  + 0.3*arrow.dist[2]),
                          shape.top[2,] + c(0.2*arrow.dist[1],  + 0.4*arrow.dist[2]),
                          shape.top[2,] + c(0, + 0.4*arrow.dist[2]))
    
    if(shape=="triangle"){
      arrow.points[,2] <- arrow.points[,2] - 0.2*arrow.dist[2]
    }
    
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
    
  }  
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...)
  # return(cbind(xc,yc))
}