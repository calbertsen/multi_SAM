##' Print output from multisam.fit
##'
##' @title Print result from multisam.fit
##' @param x A msam object
##' @param ... Other parameters passes to logLik
##' @export
print.msam <- function(x, ...){
    txt <- "Multi-SAM model with %i stock%s: log likelihood is %.4f. Convergence %s.\n"
    cat(sprintf(txt,
            length(x),
            ifelse(length(x)==1,"","s"),
            logLik(x,...),
            ifelse(0==attr(x,"m_opt")$convergence, "OK", "failed")
            ))
    invisible(x)
}

