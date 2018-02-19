##' @method logLik msam
##' @inheritParams stats::logLik
##' @importFrom stats logLik
##' @export
logLik.msam <- function(object, ...) {
    opt <- attr(object,"m_opt")
    val <-  -opt$objective
    p <- length(opt$par)
    N <- nobs(object)
    class(val) <- "logLik"
    attr(val, "nall") <- attr(N,"nall")
    attr(val, "nobs") <- as.numeric(N)
    attr(val, "df") <- p
    return(val)
}
