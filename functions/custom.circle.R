custom.circle <- function(x,y,r,rads=c(0,2),nsteps=100,
                          screen=NA, ...){  
  
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
  
  # now we need to convert our points along the arc back into plot units
  xc<-(xc / plot.window[1]) * (x.axis[2]-x.axis[1]) + x.axis[1]
  yc<-(yc / plot.window[2]) * (y.axis[2]-y.axis[1]) + y.axis[1]
  
  polygon(xc,yc,...) 
  
  return(cbind(xc,yc))
}