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


##' @method print msamforecast
##' @export
print.msamforecast <- function(x, ...){
    sn <- getStockNames(attr(x,"fit"))
    cat("\n")
    for(i in 1:length(sn)){
        cat(sprintf("Multi-SAM forecast: %s\n\n",sn[i]))
        print(x[[i]])
        cat("\n")
    }
    invisible(x)
}

##' @method print msam_hcr
##' @export
print.msam_hcr <- function(x, ...){
    print(x$forecast)
}

##' @method print msam_referencepoints
##' @export
print.msam_referencepoints <- function(x, ...){
    sn <- names(x)
    cat("\n")
    for(i in 1:length(sn)){
        cat(sprintf("Multi-SAM reference points: %s\n\n",sn[i]))
        print(x[[i]])
        cat("\n")
    }
    invisible(x)
}

