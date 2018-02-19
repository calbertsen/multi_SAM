##' Extract parameter estimates of msam object
##'
##' @title Extract parameter estimates of msam object
##' @param object msam object
##' @param ... extra arguments not used
##' @return a vector of parameter estimates
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats coef
##' @export
coef.msam <- function(object, ...){
    sdrep <- attr(object,"m_sdrep")
    ret <- sdrep$par.fixed
  attr(ret,"cov") <- sdrep$cov.fixed
  attr(ret,"sd") <- sqrt(diag(sdrep$cov.fixed))
  class(ret)<-"msamcoef"
  ret
}
##' Print msam coef object
##'
##' @param x msamcoef object
##' @param ... Not used
##' @author Christoffer Moesgaard Albertsen
##' @export
print.msamcoef <- function(x, ...){
    y <- as.vector(x)
    names(y) <- names(x)
    print(y)
}
