##' @method nobs msam
##' @export
##' @inheritParams stats::nobs
##' @importFrom stats nobs
##' 
nobs.msam <- function(object, ...) {
    obj <- attr(object,"m_obj")
    N0 <- lapply(obj$env$data$sam,function(s)s$nobs)
    N <- do.call("sum",N0)
    attr(N,"nall") <- N0
    return(N)
}
