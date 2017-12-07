##' Function to add transparency to a color
##'
##' @title Add transparency to a color
##' @param name Name of color
##' @param alpha Alpha level
##' @return A hex color string
##' @author Christoffer Moesgaard Albertsen
##' @importFrom grDevices col2rgb rgb
addTrans <- Vectorize(function(name,alpha=1){
    arg <- as.list(grDevices::col2rgb(name)/255)
    names(arg) <- c("red","green","blue")
    arg$alpha <- alpha
    do.call(grDevices::rgb,arg)
})



getLowerBounds <- function(parameters){
    f <- get("getLowerBounds", envir = asNamespace("stockassessment"), inherits = FALSE)
    r <- f(parameters)
    r$RE <- rep(-10,length(parameters$RE))
    return(r)
}

getUpperBounds <- function(parameters){
    f <- get("getUpperBounds", envir = asNamespace("stockassessment"), inherits = FALSE)
    r <- f(parameters)
    r$RE <- rep(10,length(parameters$RE))
    return(r)
}
