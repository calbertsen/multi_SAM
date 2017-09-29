##' @method AIC msamfit
##' @export
##' @inheritParams stats::AIC
AIC.msamfit <- function(object,...,k=2)
    return(object$opt$objective + k * length(object$opt$par))
