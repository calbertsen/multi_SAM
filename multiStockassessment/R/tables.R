##################################################################################################
##################################################################################################
#### Table methods modified from the stockassessment package to fit with multiStockassessment ####
##################################################################################################
##################################################################################################


##' Table helper
##'
##' @param fit msam object
##' @param what quoted name of what to extract
##' @param x rownames of table
##' @param trans function to be applied
##' @param returnList  If true, a list of matrices is returned
##' @param ... extra arguments not used
##' @return matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
tableit <-function (fit, what,
                    x=lapply(attr(fit,"m_data")$sam,function(x)x$years),
                    trans=function(x)x,
                    returnList=FALSE,
                    ...){
    UseMethod("tableit")
}

##' @rdname tableit
##' @method tableit msam
##' @export
tableit.msam <- function (fit, what,
                         x=lapply(attr(fit,"m_data")$sam,function(x)x$years),
                         trans=function(x)x,
                         returnList=FALSE,
                         ...){
    d <- attr(fit,"m_data")
    if(is.list(x)){
        xstock <- x
        x <- sort(unique(unlist(x)))
    }else{
        xstock <- replicate(length(fit),x,simplify=FALSE)
    }
    sdrep <- attr(fit,"m_sdrep")
    idx <- grep(what,names(sdrep$value))
    y<-sdrep$value[idx]
    ci<-y+sdrep$sd[idx]%o%c(-2,2)
    stockIndx <- as.numeric(factor(sub(paste0("_",what),"",names(sdrep$value)[idx])))
    ret <- list()

    for(i in unique(stockIndx)){
        m <- trans(unname(cbind(y[stockIndx==i],ci[stockIndx==i,])))
        rownames(m) <- xstock[[i]]
        colnames(m) <- c("Estimate","Low","High")
        ret[[i]] <- m
    }
    names(ret) <- getStockNames(fit)
    if(returnList)
        return(ret)
    return(cbindYearTables(ret))
}



##' SSB Table
##'
##' Table of estimated spawning stock biomasses.
##' 
##' @title SSB Table
##' @param fit msam object
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment ssbtable
##' @method ssbtable msam
##' @export
ssbtable.msam <- function(fit,...){
    ret <- tableit(fit,"logssb", trans=exp,...)
    return(ret)
}

##' TSB Table
##'
##' Table of estimated total stock biomasses.
##' 
##' @title TSB Table
##' @param fit msam object
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment tsbtable
##' @method tsbtable msam
##' @export
tsbtable.msam <- function(fit,...){
    ret <- tableit(fit,"logtsb", trans=exp,...)
    return(ret)
}



##' Fbar Table
##'
##' Table of estimated average fishing mortalities
##' 
##' @title Fbar Table
##' @param fit msam object
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment fbartable
##' @method fbartable msam
##' @export
fbartable.msam <- function(fit,...){
    ret <- tableit(fit,"logfbar", trans=exp,...)
    return(ret)
}

##' Recruitment Table
##'
##' Table of estimated recruitment to the fisheries.
##' 
##' @title Recruitment Table
##' @param fit msam object
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment rectable
##' @method rectable msam
##' @export
rectable.msam <- function(fit,...){
    ret <- tableit(fit,"logR", trans=exp,...)
    return(ret)
}


##' Catch Table
##'
##' Table of predicted catch in weight
##' 
##' @title Catch Table
##' @param fit msam object
##' @param obs.show should observed catches be included?
##' @param returnList  If true, a list of matrices is returned
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment catchtable
##' @method catchtable msam
##' @export
catchtable.msam <- function(fit, obs.show = FALSE,returnList = FALSE, ...){
    CW <- lapply(attr(fit,"m_data")$sam,function(x)x$catchMeanWeight)
    xx <- lapply(CW,function(x) as.integer(rownames(x)))
    ret <- tableit(fit, x=xx, "logCatch", trans=exp, returnList = TRUE, ...)
    if(obs.show){
        aux <- lapply(attr(fit,"m_data")$sam,function(x)x$aux)
        logobs <- lapply(attr(fit,"m_data")$sam,function(x)x$logobs)
        .goget <- Vectorize(function(y,a,s){
            ret <- exp(logobs[[s]][aux[[s]][,"fleet"]==1 & aux[[s]][,"year"]==y & aux[[s]][,"age"]==a])
            ifelse(length(ret)==0,0,ret)
        })
        sop <- sapply(1:length(fit),
                      function(i)rowSums(outer(rownames(CW[[i]]), colnames(CW[[i]]), .goget,s=i)*CW[[i]], na.rm=TRUE),
                      simplify=FALSE)
        for(i in 1:length(ret))
            ret[[i]] <- cbind(ret[[i]],"sop.catch"=sop[[i]])
    }
    if(returnList)
        return(ret)
    return(cbindYearTables(ret))
}

##' N Table
##'
##' Table of estimated numbers-at-age
##' 
##' @title N Table
##' @param fit msam object
##' @param returnList If true, a list of matrices is returned
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment ntable
##' @method ntable msam
##' @export
ntable.msam <- function(fit, returnList = FALSE, ...){
    ret <- lapply(par2list(attr(fit,"m_pl")$logN),function(x)exp(t(x)))
    d <- attr(fit,"m_data")$sam
    for(i in 1:length(d)){
        colnames(ret[[i]]) <- paste("Age",d[[i]]$minAge:d[[i]]$maxAge)
        rownames(ret[[i]]) <- d[[i]]$years
    }
    names(ret) <- getStockNames(fit)
    if(returnList)
        return(ret)
    return(cbindYearTables(ret))
}

##' F-at-age table
##'
##' Table of estimated age-wise fishing mortality
##' 
##' @title F-at-age Table
##' @param fit msam object
##' @param returnList  If true, a list of matrices is returned
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment faytable
##' @method faytable msam
##' @export
faytable.msam <- function(fit, returnList = FALSE, ...){
    d <- attr(fit,"m_data")$sam
    idx <- lapply(d,function(xx)xx$keyLogFsta[1,]+2)
    Flist <- par2list(attr(fit,"m_pl")$logF)
    ret <- sapply(1:length(Flist),function(i){
        r <- cbind(NA,exp(t(Flist[[i]])))[,idx[[i]]]
        r[,idx[[i]]==0] <- 0
        r
        },simplify=FALSE)
    for(i in 1:length(d)){
        colnames(ret[[i]]) <- paste("Age",d[[i]]$minAge:d[[i]]$maxAge)
        rownames(ret[[i]]) <- d[[i]]$years
    }
    names(ret) <- getStockNames(fit)
    if(returnList)
        return(ret)
    return(cbindYearTables(ret))
}




##' Parameter table
##'
##' Table of estimated parameters
##' 
##' @title Parameter table
##' @param fit msam object
##' @param ... Other parameters 
##' @return A matrix of estimates and confidence intervals
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment partable
##' @method partable msam
##' @export
partable.msam <- function(fit, ...){
    param <- coef(fit)
    pl <- attr(fit,"m_pl")
    stockNames <- getStockNames(fit)
    nam <- names(param)
    namstock <- rep("",length(nam))
    for(na in unique(nam)){
        pattr <- attributes(pl[[na]])
        if(!is.null(pattr$vdim)){
            sv <- rep(stockNames,times=pattr$vdim)
        }else if(!is.null(pattr$rdim) & !is.null(pattr$cdim)){
            sv <- rep(stockNames,times=pattr$rdim * pattr$cdim)
        }else{
            sv <- rep("",length(pl[[na]]))
        }
        namstock[nam==na] <- sv
    }
    nam <- paste(nam,namstock,sep="_")
    dup <- duplicated(nam)
    namadd <- rep(0, length(nam))
      for (i in 2:length(dup)) {      
        if(dup[i])namadd[i] <- namadd[i-1] + 1
    }
    nam <- paste(nam, namadd, sep = "_")
    ret <- cbind(param, attr(param,"sd"))
    ex <- exp(ret[,1])
    lo <- exp(ret[,1]-2*ret[,2])
    hi <- exp(ret[,1]+2*ret[,2])
    ret <- cbind(ret,ex,lo,hi)
    colnames(ret) <- c("par", "sd(par)", "exp(par)", "Low", "High")
    rownames(ret) <- nam
    return(ret)
}

##' Model table
##'
##' @title Model table
##' @param fits msam object
##' @param ... extra arguments
##' @return A matrix of model information
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment modeltable
##' @method modeltable msam
##' @export
modeltable.msam <- function(fits,...){
    modeltable(c(fits),...)
}

##' Model table
##'
##' @title Model table
##' @param fits msamset object
##' @param ... extra arguments
##' @return A matrix of model information
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment modeltable
##' @importFrom stats AIC pchisq
##' @method modeltable msamset
##' @export
modeltable.msamset <- function(fits,...){
    nam <- paste("M", 1:length(fits), sep="")
    logL <- sapply(fits,logLik)
     npar <- sapply(fits, function(f)attr(logLik(f),"df"))
    aic <- sapply(fits, stats::AIC)
    res <- cbind("log(L)"=logL, "#par"=npar, "AIC"=aic)
    rownames(res) <- nam
    o <- 1:length(fits)
    if(length(fits)==2){
        o <- order(npar, decreasing=TRUE)
        if(npar[o[1]]>npar[o[2]]){
            df <- npar[o[1]]>npar[o[2]]
            D <- 2*(logL[o[1]]-logL[o[2]])
            P <- 1 - stats::pchisq(D,df)
            cnam <- paste0("Pval( ",nam[o[1]]," -> ",nam[o[2]], " )")
            res <- cbind(res, c(NA, P)[o])
            colnames(res)[ncol(res)] <- cnam
        }
    }
    return(res[o,,drop=FALSE])
}



## ypr.msam <- function(fit, Flimit=2, Fdelta=0.01, aveYears=sapply(1:length(fit),function(i)min(15,length(fit[[i]]$data$years))), ageLimit=100,...){
##     tabs <- sapply(1:length(fit),function(i)stockassessment::ypr(fit[[i]],
##                                                                  Flimit=Flimit,
##                                                                  Fdelta=Fdelta,
##                                                                  aveYears=aveYears[i],
##                                                                  ageLimit=ageLimit,
##                                                                  ...),
##                    simplify=FALSE
##                    )
##     names(tabs) <- names(fit)
##     if(is.null(names(tabs)))
##         names(tabs) <- paste("Stock",1:length(tabs))
##     return(tabs)
## }
