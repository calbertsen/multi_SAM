addTrans <- Vectorize(function(name,alpha = 1){
    if(is.na(name) | is.null(name))
        return(name)
    arg <- as.list(grDevices::col2rgb(name)/255)
    names(arg) <- c("red","green","blue")
    arg$alpha <- alpha
    do.call(grDevices::rgb,arg)
})

##' @method addRecruitmentCurve msam
##' @importFrom utils head
##' @importFrom stockassessment addRecruitmentCurve
##' @export
addRecruitmentCurve.msam <- function(fit,
                                    CI = TRUE,
                                    col = .plotcols.crp(length(fit)),
                                    cicol = addTrans(.plotcols.crp(length(fit)),0.3),
                                    plot = TRUE,
                                    PI = FALSE,
                                    picol = addTrans(.plotcols.crp(length(fit)),0.5),
                                    pilty = 2,
                                    ...){
    X <- summary(fit, returnList = TRUE)
    R <- unlist(lapply(X, function(x) x[, 1]))
    S <- unlist(lapply(X, function(x) x[, 4]))
    ## cf <- fit$sdrep$cov.fixed
    ## covEst <- cf[rownames(cf) %in% c("rec_pars"), colnames(cf) %in% c("rec_pars"), drop = FALSE]
    covEst <- attr(fit,"m_sdrep")$covRecPars
    m <- attr(fit,"m_obj")$env$map$rec_pars
    if(is.null(m)){
        covar <- covEst
    }else{
        covar <- covEst[m,m, drop = FALSE]
        covar[is.na(covar)] <- 0
    }

    stockIndex <- factor(gsub("(^SAM_)([[:digit:]]+)(_.+$)","\\2",colnames(covEst),perl=TRUE))
    pl <- attr(fit,"m_pl")
    rpList <- splitParameter(pl$rec_pars)
    covList <- lapply(levels(stockIndex), function(ii){
        indx <- which(stockIndex == ii)
        covar[indx,indx]
    })
    lsdlnList <- splitParameter(attr(fit,"m_pl")$logSdLogN)
    recPredSdList <- lapply(seq_along(lsdlnList), function(si){
        exp(lsdlnList[[si]][fit[[si]]$conf$keyVarLogN[1]+1])
    })
    
    doOne <- function(sii){
        srmc <- fit[[sii]]$conf$stockRecruitmentModelCode
        if(srmc[sii] %in% c(-2,-1,0, 3,62)){
            warning("addRecruitmentCurve is not implemented for time series models.")
        }
        rec_pars <- rpList[[sii]]
        covar <- covList[[sii]]
        recPredSd <- recPredSdList[[sii]]        
        srfit <- function(logssb, year = NA_real_, lastR = NA_real_) {
            v <- stockassessment:::stockRecruitmentModelR(logssb, rec_pars, 
                                                          srmc, fit[[sii]]$conf$constRecBreaks, 
                                                          year, lastR)
            val <- v$logRecruits
            g <- matrix(v$Gradient_recpars, ncol = 1)
            valsd <- as.vector(sqrt(t(g) %*% covar %*% g))
            pisig <- exp(splitParameter(pl$logSdLogN)[[sii]][fit[[sii]]$conf$keyVarLogN[1] + 
                                                             1])
            res <- exp(val)
            attr(res, "ci_low") <- exp(val - 2 * valsd)
            attr(res, "ci_high") <- exp(val + 2 * valsd)
            attr(res, "pi_low") <- exp(val - 2 * pisig)
            attr(res, "pi_high") <- exp(val + 2 * pisig)
            return(res)
        }
        getSSBtab <- function(x, year = NA) {
            tmp <- srfit(x, year)
            c(Estimate = as.vector(tmp),
              CIlow = as.vector(attr(tmp,"ci_low")),
              CIhigh = as.vector(attr(tmp,"ci_high")),
              PIlow = as.vector(attr(tmp,"pi_low")),
              PIhigh = as.vector(attr(tmp,"pi_high"))
              )
        }
        if(plot){
            ssbMax <- max(c(max(S), par("usr")[2]))
        }else{
            ssbMax <- max(S)
        }
        ssb <- seq(1e-5,
                   ssbMax, ## max(S),
                   len = 2000)
        if(srmc == 3){
            brks <- c(min(fit[[sii]]$data$years),ceiling(fit[[sii]]$conf$constRecBreaks),max(fit[[sii]]$data$years))
            labels <- levels(cut(0,brks,dig.lab=4))
            mid <- head(brks,-1) + diff(brks)/2
            tabList <- lapply(as.list(mid), function(y) sapply(log(ssb),getSSBtab, year = y))
            names(tabList) <- labels
        }else{
            tabList <- list(sapply(log(ssb), getSSBtab))
        }
        transp <- Vectorize(function(col){
            arg <- as.list(grDevices::col2rgb(col,TRUE)[,1])
            arg$alpha <- 0.1 * 255
            arg$maxColorValue <- 255
            do.call(grDevices::rgb, arg)
        })    
        if(plot){           
            if(CI)               
                for(i in 1:length(tabList)){
                    polygon(c(ssb, rev(ssb)),
                            c(tabList[[i]]["CIlow",],
                              rev(tabList[[i]]["CIhigh",])),
                            col = ifelse(length(tabList) == 1,cicol[sii],transp(i)),
                            border = NA)
                }
            if(srmc %in% c(90,91,92))
                abline(v = exp(fit[[sii]]$conf$constRecBreaks), col = "grey")
            for(i in 1:length(tabList))
                lines(ssb,tabList[[i]]["Estimate",], col = ifelse(length(tabList) == 1,col[sii],i), lwd = 3)
            if(length(tabList) > 1){
                for(i in 1:length(tabList)){
                    text(par("usr")[1],tabList[[i]]["Estimate",1],names(tabList)[[i]],pos=4, col = ifelse(length(tabList) == 1,cicol[sii],i))
                }
            }
            if(PI){
                for(i in 1:length(tabList)){               
                    lines(ssb,tabList[[i]]["PIlow",], col = ifelse(length(tabList) == 1,picol[sii],i), lwd = 3, lty = pilty)
                    lines(ssb,tabList[[i]]["PIhigh",], col = ifelse(length(tabList) == 1,picol[sii],i), lwd = 3, lty = pilty)
                }
            }
        }
        srfit
    }
    invisible(lapply(seq_along(fit), doOne))
}
