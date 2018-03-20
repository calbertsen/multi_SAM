##' Collect a number of msam objects into an msamset
##'
##' @title Collect msam objects
##' @param ... msam objects
##' @importFrom methods is
##' @return An msamset
##' @export
c.msam<-function(...){
    ret<-list(...)
    if(!all(unlist(lapply(ret,methods::is,class2="msam"))))
        stop("All arguments must be msam objects")
    class(ret)<-"msamset"
    ret
}
