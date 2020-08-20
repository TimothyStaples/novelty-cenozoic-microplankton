plot.bg <- function(...){
  
  rect(xleft=par("usr")[1], 
       xright=par("usr")[2], 
       ybottom=par("usr")[3], 
       ytop=par("usr")[4],
  ...)
  
  
}