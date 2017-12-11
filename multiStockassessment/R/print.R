##' Print output from multisam.fit
##'
##' @title Print result from multisam.fit
##' @param x A msam object
##' @param ... Other parameters passes to logLik
##' @export
print.msam <- function(x, ...){
    cat("Multi-SAM model with",length(x),"stocks: log likelihood is",
        logLik.msam(x,...),
        "Convergence",ifelse(0==attr(x,"m_opt")$convergence, "OK\n", "failed\n"))
}

