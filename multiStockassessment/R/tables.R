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
    idx <- grep(paste0("_",what,"$"),names(sdrep$value))
    y<-sdrep$value[idx]
    ci<-y+sdrep$sd[idx]%o%c(-2,2)
    stockIndx <- as.numeric(factor(sub(paste0("_",what),"",names(sdrep$value)[idx])))
    ret <- list()

    for(i in unique(stockIndx)){
        syears <- d$sam[[i]]$years
        allyears <- sort(unique(c(syears,xstock[[i]])))
        xindx <- allyears %in% xstock[[i]]
        sindx <- allyears %in% syears[1:sum(stockIndx==i)]
        m <- matrix(NA,length(allyears),3)
        m[sindx,] <- trans(unname(cbind(y[stockIndx==i],ci[stockIndx==i,])))
        m <- m[xindx,]
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


##' Yield per recruit calculation
##'
##' @param fit msam object
##' @param Flimit Upper limit for Fbar
##' @param Fdelta increments on the Fbar axis
##' @param aveYears A list/vector of same length as the number of stocks with number of years to average over 
##' @param ageLimit Oldest age used
##' @param ... not used
##' @return A list of samypr objects
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment ypr
##' @method ypr msam
##' @export
ypr.msam <- function(fit,
                     Flimit=2,
                     Fdelta=0.01,
                     aveYears=lapply(attr(fit,"m_data")$sam,function(x)min(15,length(x$years))),
                     ageLimit=100,
                     ...){
    if(!is.list(aveYears))
        aveYears <- as.list(aveYears)    
    d <- attr(fit,"m_data")$sam
    barAges <- lapply(d,function(x)do.call(":",as.list(x$fbarRange))+(1-x$minAge))
    last.year.used <- lapply(d,function(x)max(x$years))
    idxno<-sapply(1:length(d),function(i)which(d[[i]]$years==last.year.used[[i]]),simplify = FALSE)
    F <- lapply(faytable(fit,returnList = TRUE),function(x){
        y <- t(x); y[is.na(y)]<-0; y})
    Fbar <- fbartable(fit, returnList = TRUE)
    
    sel<-function(stock){
        Sa<-rep(0,nrow(F[[stock]]))
        K<-0
        indxy <- 0:(aveYears[[stock]]-1)
        Say <- F[[stock]][,idxno[[stock]]-indxy]
        Ky<- Fbar[[stock]][idxno[[stock]]-indxy,1]
        return(rowSums(Say)/sum(Ky))
    }

    extend<-function(x,len=100){
        ret<-numeric(len)
        ret[1:length(x)]<-x
        ret[-c(1:length(x))]<-x[length(x)]
        ret
    }

    ave.sl<-sapply(1:length(fit),sel,simplify=FALSE)
    ave.sw <- sapply(1:length(d),
                     function(i)colMeans(d[[i]]$stockMeanWeight[(idxno[[i]]-aveYears[[i]]+1):idxno[[i]],,drop=FALSE]),simplify=FALSE)
    ave.cw <- sapply(1:length(d),
                     function(i)colMeans(d[[i]]$catchMeanWeight[(idxno[[i]]-aveYears[[i]]+1):(idxno[[i]]-1),,drop=FALSE]),simplify=FALSE)
    ave.pm <- sapply(1:length(d),
                     function(i)colMeans(d[[i]]$propMat[(idxno[[i]]-aveYears[[i]]+1):idxno[[i]],,drop=FALSE]),simplify=FALSE)
    ave.nm <- sapply(1:length(d),
                     function(i)colMeans(d[[i]]$natMor[(idxno[[i]]-aveYears[[i]]+1):idxno[[i]],,drop=FALSE]),simplify=FALSE)
    ave.lf <- sapply(1:length(d),
                     function(i)colMeans(d[[i]]$landFrac[(idxno[[i]]-aveYears[[i]]+1):(idxno[[i]]-1),,drop=FALSE]),simplify=FALSE)
    ave.cw.land <- sapply(1:length(d),
                     function(i)colMeans(d[[i]]$landMeanWeight[(idxno[[i]]-aveYears[[i]]+1):(idxno[[i]]-1),,drop=FALSE]),simplify=FALSE)

    N <- lapply(as.list(rep(ageLimit,length(fit))),numeric)
    N <- lapply(N,function(x){x[1] <- 1.0;x})
    M <- lapply(ave.nm,extend,len = ageLimit)
    sw <- lapply(ave.sw,extend,len = ageLimit)
    cw <- lapply(ave.cw.land,extend,len = ageLimit)
    pm <- lapply(ave.pm,extend,len = ageLimit)
    lf <- lapply(ave.lf,extend,len = ageLimit)


    doOne <- function(stock){
        deltafirst <- 0.00001
        delta <- Fdelta
        scales<-c(0, deltafirst, seq(0.01, Flimit, by=delta))
        yields<-numeric(length(scales))
        ssbs<-numeric(length(scales))
        for(i in 1:length(scales)){
            scale<-scales[i]
            F<-extend(ave.sl[[stock]]*scale,ageLimit)
            Z<-M[[stock]]+F
            for(a in 2:length(N[[stock]])){
                N[[stock]][a]<-N[[stock]][a-1]*exp(-Z[a-1])  
            }
            C<-F/Z*(1-exp(-Z))*N[[stock]]*lf[[stock]]  
            Y<-sum(C*cw[[stock]])
            yields[i]<-Y
            ssbs[i]<-sum(N[[stock]]*pm[[stock]]*sw[[stock]])
        }
        fmaxidx<-which.max(yields)
        fmax<-scales[fmaxidx]

        deltaY<-diff(yields)
        f01idx<-which.min((deltaY/delta-0.1*deltaY[1]/deltafirst)^2)+1
        f01<-scales[f01idx]
        
        f35spridx<-which.min((ssbs-0.35*ssbs[1])^2)+1
        f35<-scales[f35spridx]
    
        fbarlab <- substitute(bar(F)[X - Y], list(X = d[[stock]]$fbarRange[1], Y = d[[stock]]$fbarRange[2]))
        ret<-list(fbar=scales, ssb=ssbs, yield=yields, fbarlab=fbarlab, f35=f35, f01=f01, fmax=fmax, 
                  f35Idx=f35spridx, f01Idx=f01idx, fmaxIdx=fmaxidx)
        class(ret)<-"samypr"
        return(ret)
    }

    res <- lapply(1:length(fit),doOne)
    names(res) <- getStockNames(fit)
    class(res) <- "msamypr"
    return(res)
}



