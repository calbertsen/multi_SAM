##' Summary of msam object
##'
##' @param object msam object
##' @param returnList If true, a list of table is returned
##' @param digits digits for rounding output. Vector of length 3 corresponding to recruitment, ssb, and fbar.
##' @param ... not used
##' @return a summary table
##' @author Christoffer Moesgaard Albertsen
##' @method summary msam
##' @export
summary.msam <- function(object, returnList = FALSE,digits=c(0,0,3), ...){
    digits <- rep(digits,length=3)
    reclist <- rectable(object,returnList = TRUE)
    ssblist <- ssbtable(object,returnList = TRUE)
    fbarlist <- fbartable(object,returnList = TRUE)
    d <- attr(object,"m_data")$sam
    rl <- list()
    for(i in 1:length(object)){
        ret <- cbind(round(reclist[[i]],digits[1]),
                     round(ssblist[[i]],digits[2]),
                     round(fbarlist[[i]],digits[3]))
        colnames(ret)[c(1,4,7)] <- c(paste0("R(age ",d[[i]]$minAge,")"),
                                     "SSB",
                                     paste0("Fbar(",
                                            d[[i]]$fbarRange[1],
                                            "-", d[[i]]$fbarRange[2],
                                            ")")
                                     )
        rl[[i]] <- ret
    }
    names(rl) <- getStockNames(object)
    if(returnList)
        return(rl)
    return(cbindYearTables(rl))
}
