cut.novel <- function(novel.list,
                      cutoff){
  lapply(novel.list,
         function(x){
           
           bin.raw <- as.numeric(as.character(x$bins))
           bin.sort <- sort(as.numeric(as.character(x$bins)))
           
           x$bin.n <- nrow(x)+1 - match(bin.raw, bin.sort)
           
           x[x$bin.n > cutoff, ]
         })
  
}
