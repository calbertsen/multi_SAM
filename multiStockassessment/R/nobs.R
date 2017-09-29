##' @method nobs msamfit
##' @export
##' @inheritParams stats::nobs
##' @importFrom stats nobs
##' 
nobs.msamfit <- function(object, ...) {
    return(do.call("sum",lapply(object$obj$env$data$sam,function(s)s$nobs)))
}
