##' Function to add transparency to a color
##'
##' @title Add transparency to a color
##' @param name Name of color
##' @param alpha Alpha level
##' @return A hex color string
##' @author Christoffer Moesgaard Albertsen
##' @importFrom grDevices col2rgb rgb
addTrans <- Vectorize(function(name,alpha=1){
    if(is.na(name) | is.null(name))
        return(name)
    arg <- as.list(grDevices::col2rgb(name)/255)
    names(arg) <- c("red","green","blue")
    arg$alpha <- alpha
    do.call(grDevices::rgb,arg)
})

is.msam <- function(x){
    return(inherits(x,"msam"))
}
is.msamset <- function(x){
    return(inherits(x,"msamset"))
}

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


par2list <- function(x){
    if(!is.null(attr(x,"vdim"))){
        return(.Call("vecpar2list",x))
    }else if(!is.null(attr(x,"cdim")) & !is.null(attr(x,"rdim"))){
        return(.Call("matpar2list",x))
    }else{
        stop("Not a valid parameter array")
    }
}
    
getStockNames <- function(x){
    if(!is.msam(x))
        stop("x must be an msam object")
    stocknames <- names(x)
    if(is.null(stocknames))
        stocknames <- paste("Stock",1:length(x))
    return(stocknames)
}


cbindYearTables <- function(x){
    ylist <- lapply(x,function(y)as.integer(rownames(y)))
    ycom <- sort(unique(unlist(ylist)))
    alist <- lapply(x,colnames)
    if(!is.null(names(x))){
        nx <- names(x)
        for(i in 1:length(alist))
            alist[[i]] <- paste0(nx[i],": ",alist[[i]])
    }
    aall <- unlist(alist)
    ret <- matrix(NA,length(ycom),length(aall))
    rownames(ret)<-ycom
    colnames(ret)<-aall
    colIndx <- cumsum(c(1,unlist(lapply(alist,length))))
    for(i in 1:length(alist)){
        m <- x[[i]]
        xv <- which(ycom %in% ylist[[i]])
        ret[xv,colIndx[i]:(colIndx[i+1]-1)] <- m
    }
    return(ret)

}
