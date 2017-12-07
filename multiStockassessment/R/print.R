##' @title Print result from multisam.fit
##' @param x A msam object
##' @param ... Other parameters passes to logLik
##' @export
print.msam <- function(x, ...){
    cat("SAM model with correlated survival: log likelihood is",
        logLik(x,...),
        "Convergence",ifelse(0==attr(x,"m_opt")$convergence, "OK\n", "failed\n"))
}

