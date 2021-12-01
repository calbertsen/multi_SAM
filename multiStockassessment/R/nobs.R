##' @method nobs msam
##' @export
##' @inheritParams stats::nobs
##' @importFrom stats nobs
##' 
nobs.msam <- function(object, ...) {
    obj <- attr(object,"m_obj")
    dat <- attr(object,"m_data")
    N0 <- lapply(obj$env$data$sam,function(s)s$nobs)
    N <- do.call("sum",N0) + length(dat$sharedObs$logobs)
    attr(N,"nall") <- N0
    return(N)
}
