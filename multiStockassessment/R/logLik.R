##' @method logLik msamfit
##' @export
##' @inheritParams stats::logLik
##' 
##' 
logLik.msamfit <- function(object, ...) {
    val <-  -object$opt$objective
    p <- length(object$opt$par)
    N0 <- lapply(object$obj$env$data$sam,function(s)s$nobs)
    N <- do.call("sum",N0)
    class(val) <- "logLik"
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "df") <- p
    return(val)
}
