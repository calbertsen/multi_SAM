##' @title Collect msam objects
##' @param ... msam objects
##' @return An msamset
##' @export
c.msam<-function(...){
  ret<-list(...)
  class(ret)<-"msamset"
  ret
}
