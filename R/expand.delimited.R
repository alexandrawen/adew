#' @export
#https://www.r-bloggers.com/2012/11/expand-delimited-columns-in-r/
expand.delimited <- function(x, col1=1, col2=2, sep=",") {
  rnum <- 1
  expand_row <- function(y) {
    factr <- y[col1]
    strng <- toString(y[col2])
    expand <- strsplit(strng, sep)[[1]]
    num <- length(expand)
    factor <- rep(factr,num)
    return(as.data.frame(cbind(factor,expand),
                         row.names=seq(rnum:(rnum+num)-1)))
    rnum <- (rnum+num)-1
  }
  expanded <- apply(x,1,expand_row)
  df <- do.call("rbind", expanded)
  names(df) <- c(names(x)[col1],names(x)[col2])
  return(df)
}