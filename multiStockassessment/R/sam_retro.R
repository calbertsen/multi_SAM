
block <- function(...){
    requireNamespace("Matrix")
    x <- list(...)
    ii <- sapply(x, is.matrix)
    x[ii] <- lapply(x[ii],as,Class="generalMatrix")
    jj <- sapply(x, is.numeric)
    x[jj & !ii] <- lapply(x[jj & !ii],function(y)Matrix::Diagonal(length(y),y))
    n <- sapply(x,nrow)
    rn <- do.call(c,lapply(x, rownames))
    cn <- do.call(c,lapply(x, colnames))
    M <- Matrix::Matrix(0,sum(n),sum(n), dimnames = list(rn,cn))
    nc <- cumsum(c(0,n))
    for(i in 1:length(x)){
        jj <- nc[i]+1:n[i]
        M[jj,jj] <- x[[i]]
    }
    M    
}

get_triplet <- function(x){
    requireNamespace("Matrix")
    y <- as(x,"TsparseMatrix")
    if(is.null(y@uplo)){
        isU <- yv$i <= yv$j
        r <- data.frame(i = y@i[isU]+1,
                        j = y@j[isU]+1,
                        x = y@x[isU])
    }else if(y@uplo == "L"){
        r <- data.frame(i = y@j+1,
                        j = y@i+1,
                        x = y@x)
    }else{
        r <- data.frame(i = y@i+1,
                        j = y@j+1,
                        x = y@x)
    }
    attr(r,"Dims") <- y@Dim
    r
}

sparse_sym_block <- function(...){
     requireNamespace("Matrix")
   x <- list(...)
    ii <- sapply(x, is, class2 = "sparseMatrix")
    x[!ii] <- lapply(x[!ii],function(y) as(y,"sparseMatrix")) #Matrix::sparseMatrix(i=seq_along(y),j=seq_along(y),x=as.numeric(y)))
    y <- lapply(x, get_triplet)
    n <- sapply(x,nrow)
    nc <- cumsum(n)
    for(i in 2:length(y)){
        y[[i]]$i <- y[[i]]$i + nc[i-1]
        y[[i]]$j <- y[[i]]$j + nc[i-1]
    }
    yv <- do.call("rbind",y)
    #isU <- yv$i =< yv$j
    Matrix::sparseMatrix(i=yv$i, j=yv$j, x=yv$x, dims = c(sum(n),sum(n)), symmetric = TRUE)
}

makeSymPosDef <- function(x, tol = sqrt(.Machine$double.eps), warn = FALSE){
    requireNamespace("Matrix")
    ## Sym
    ##v <- 0.5 * (x + t(x))
    if(is.matrix(x))
        x <- as(x,"generalMatrix")
    v <- Matrix::symmpart(x)
    #ee <- eigen(v, symmetric = TRUE)
    #ee$values <- pmax(ee$values,1e-8 / max(ee$values))
    ##Sig1 <- ee$vectors %*% diag(x=ee$values) %*% solve(ee$vectors)
    Sc <- Matrix::expand2(Matrix::Schur(x))
    ii <- seq_len(nrow(Sc$T))
    if(all(Sc$T[cbind(ii,ii)] >= tol))
        return(v)
    if(warn){
        warning("Matrix was modified to be positive definite")
    }
    T2 <- as(Sc$T,"sparseMatrix")
    ii <- seq_len(nrow(T2))
    T2[cbind(ii,ii)] <- pmax(T2[cbind(ii,ii)],tol)
    Q <- as(Sc$Q,"sparseMatrix")
    Qi <- as(Sc$`Q.`,"sparseMatrix")
    Sig1 <- Q %*% T2 %*% Qi
    Sig1 <- Matrix::symmpart(Sig1) #0.5 * (Sig1 + Matrix::t(Sig1))
    Sig1
}

svd_solve_posdef <- function(x, eps = sqrt(.Machine$double.eps)){
    requireNamespace("Matrix")
    ss <- svd(x)
    ss$v %*% diag(1/pmax(ss$d,eps), length(ss$d), length(ss$d)) %*% t(ss$u)
}

make_quick_multisam <- function(RETRO){
    requireNamespace("Matrix")
    fitFinal <- attr(RETRO,"fit")
    framework <- .Call("getFramework", PACKAGE=fitFinal$obj$env$DLL)
    if (framework != "TMBad")
        stop("Please re-install the stockassessment and multiStockassessment packages with TMBad")
    fitList <- do.call("c",c(list(fitFinal), RETRO))
    m_obj <- multisam.fit(fitList, mohn=1, doSdreport = FALSE, run=FALSE)
    m_opt <- list(par = m_obj$par)
    ## match parameters to years
    dataYears <- fitFinal$data$years
    years <- rev(tail(dataYears,length(fitList)))
    apOrig <- do.call("rbind",lapply(fitList, function(x) sapply(split(x$opt$par,names(x$opt$par)),length)))
    ap <- split(m_opt$par,factor(names(m_opt$par),unique(names(m_opt$par))))
    parYear <- lapply(names(ap), function(nm) rep(years,times = apOrig[,nm]))
    names(parYear) <- names(ap)
    parYear <- unlist(parYear)
    info_Hx_peel <- unlist(lapply(years, function(y){ rep(max(years)-y, sum(parYear==y)) }))
    info_Hx_par <- unlist(lapply(years, function(y){ factor(names(m_opt$par[parYear==y]),unique(names(m_opt$par))) }))
    info_Hx_num <- unlist(lapply(years, function(y){ seq_along(m_opt$par[parYear==y]) }))
    j11or <- order(info_Hx_par,info_Hx_peel,info_Hx_num)
 
    all_he <- do.call(block,lapply(fitList,function(x) makeSymPosDef(x$opt$he,warn=TRUE)))[j11or,j11or]
    m_opt$he <- all_he
    attr(fitList,"m_obj") <- m_obj
    attr(fitList,"m_opt") <- m_opt
    attr(fitList,"m_pl") <- m_obj$env$parList(par=m_obj$env$last.par.best)
    class(fitList) <- c("msam", "samset","fake_msam")
    fitList
}

retro_hessian <- function(mFit, keep.diagonal = TRUE, returnSigma = TRUE){
##### Calculate hessian with correlation #####
    requireNamespace("Matrix")
    oFit <- mFit[[length(mFit)]]
    dataYears <- mFit[[1]]$data$years
    years <- rev(tail(dataYears,length(mFit)))
    m_opt <- attr(mFit,"m_opt")
    ## match parameters to years
    apOrig <- do.call("rbind",lapply(mFit, function(x) sapply(split(x$opt$par,names(x$opt$par)),length)))
    ap <- split(m_opt$par,factor(names(m_opt$par),unique(names(m_opt$par))))
    parYear <- lapply(names(ap), function(nm) rep(years,times = apOrig[,nm]))
    names(parYear) <- names(ap)
    parYear <- unlist(parYear)
    ##parYear <- unlist(lapply(ap,function(x) rep(years, each = length(x) / length(years))))
    isFirstYear <- parYear == min(parYear)
    ## Derivative of score function
    H <- m_opt$he
    ## Derivative wrt each of the old parameter sets
    Hx <- lapply(head(years,-1), function(y){
        hessian_gr(oFit$obj$gr, m_opt$par[parYear==y], method = "central", symmetrize=TRUE)
    })
    Hx <- lapply(Hx, makeSymPosDef, warn = TRUE)
    info_Hx_peel <- unlist(lapply(head(years,-1), function(y){ rep(max(years)-y, sum(parYear==y)) }))
    info_Hx_par <- unlist(lapply(head(years,-1), function(y){ factor(names(m_opt$par[parYear==y]),unique(names(m_opt$par))) }))
    info_Hx_num <- unlist(lapply(head(years,-1), function(y){ seq_along(m_opt$par[parYear==y]) }))
    ## Derivative of score function wrt new parameters (ordered by parameter then year)
    J2 <- H[!isFirstYear, !isFirstYear]
    ## Get derivative and (Gaussian) variance for data
    oldObj <- attr(mFit,"m_obj")
    dat <- oldObj$env$data
    logobsList <- lapply(dat$sam,function(x)x$logobs)
    auxList <- lapply(dat$sam,function(x){
        cbind(x$aux,fleetType = x$fleetTypes[x$aux[,"fleet"]])
    })
    auxDList <- lapply(dat$sam,function(x){
        x$auxData
    })
    logobs <- do.call("c",logobsList)
    stock <- rep(seq_along(logobsList)-1, times = sapply(logobsList,length))    
    indx <- do.call("c",lapply(logobsList,function(x)seq_along(x)-1))
    year <- do.call("c",lapply(auxList,function(x)x[,"year"]))
    dpart <- do.call("c",lapply(auxDList,function(x){
        if(ncol(x) < 5)
            return(rep(0,nrow(x)))
        x[,5]
    }))
    fleetType <- do.call("c",lapply(auxList,function(x)x[,"fleetType"]))
    fleet <- do.call("c",lapply(as.list(seq_along(auxList)),
                                function(i){
                                    prefix <- gsub(" ","0",
                                                   formatC(i,width=nchar(as.character(length(auxList)))))
                                    paste0(prefix,"_",auxList[[i]][,"fleet"])
                                }))
    age <- do.call("c",lapply(auxList,function(x)x[,"age"]))
    xtr <- ifelse(fleetType >= 80, dpart,0)
    kp <- which(!is.na(logobs) & (fleetType < 80))
    sortVals <- order(year[kp],stock[kp],fleet[kp],age[kp], xtr[kp])
    dat$fake_obs <- logobs[kp][sortVals]
    dat$fake_stock <- stock[kp][sortVals]
    dat$fake_indx <- indx[kp][sortVals]
    fake_year <- year[kp][sortVals]
    dat$useFakeObs <- 1
    ran <- c("logN", "logF", "missing","logSW", "logCW", "logitMO", "logNM", "logP","logitFseason", "shared_logFscale","shared_missingObs", "logGst", "logGtrip")
    pl <- c(list(fake_obs = dat$fake_obs), attr(mFit,"m_pl"))
    dat$fake_obs <- NULL
    map <- oldObj$env$map
    map$fake_obs <- factor(ifelse(fake_year <= min(years), NA, seq_along(pl$fake_obs)))
    Obj1 <- TMB::MakeADFun(dat,pl,map, random=ran, DLL="multiStockassessment",inner.control=list(maxit=100))
    isObs <- names(Obj1$par) == "fake_obs"
    Obj1$fn()
    Hy <- t(hessian_gr(Obj1$gr, Obj1$par, subset = which(isObs), method = "central", symmetrize = TRUE)) ## Subsets on rows, so transform to organize by columns
### The data appears multiple times for some years, use average
    diagA <- max(years) - fake_year[!is.na(map$fake_obs)] + 1
    A <- diag(sqrt(1/diagA),length(diagA))
    ## Approximate variance of data
    Vy <- svd_solve_posdef(A %*%Hy[which(isObs),which(isObs)]%*%A) #makeSymPosDef(svd_solve_posdef(A %*%Hy[which(isObs),which(isObs)]%*%A))
    ## Symmetrize for safety
    ## Vy <- 0.5 * (Vy + t(Vy))
    ## J1_2 is ordered by parameter then year
    J1_2 <- Hy[max(which(isObs)) + which(!isFirstYear),which(isObs)] # Derivative wrt new parameter then y
    J1_1 <- do.call("rbind",lapply(Hx, function(x) x-2*H[isFirstYear,isFirstYear]))
    ## Reorder to match J1_2
    j11or <- order(info_Hx_par,info_Hx_peel,info_Hx_num)
    J1_1 <- J1_1[j11or,]
    ## Prepare to reorder G
    info_H_peel <- c(info_Hx_peel[j11or],unlist(lapply(tail(years,1), function(y){ rep(max(years)-y, sum(parYear==y)) })))
    info_H_par <- c(info_Hx_par[j11or],unlist(lapply(tail(years,1), function(y){ factor(names(m_opt$par[parYear==y]),unique(names(m_opt$par))) })))
    info_H_num <- c(info_Hx_num[j11or],unlist(lapply(tail(years,1), function(y){ seq_along(m_opt$par[parYear==y]) })))
    ## Full gradient for Delta Method (ordered by parameter then year)
    G <- -svd_solve_posdef(J2) %*% (cbind(J1_1,J1_2) )
    G2 <- matrix(0,sum(isFirstYear),ncol(G))
    #G2[cbind(1:sum(isFirstYear),1:sum(isFirstYear))] <- 1
    Gx <- rbind(G,G2)
    ## Need to reorder to match fit
    Gx <- Gx[order(info_H_par,info_H_peel,info_H_num),]
    ii <- which(info_H_peel[order(info_H_par,info_H_peel,info_H_num)] == max(info_H_peel))
    Gx[cbind(ii,seq_along(ii))] <- 1
    ## Delta method to get correlation
    Vold <- svd_solve_posdef(oFit$opt$he)
    Sig1 <- Gx %*% block(Vold,Vy) %*% Matrix::t(Gx)
    ## Symmetrize for safety
    ## Sig1 <- makeSymPosDef(Sig1)
    ##Sig1 <- 0.5 * (Sig1 + t(Sig1))    
    ## if(forcePosDef){
    ##     Sig1 <- Matrix::nearPD(Sig1)$mat
    ## }
    if(keep.diagonal){
        noCorSig <- svd_solve_posdef(m_opt$he)
        D <- diag(sqrt(diag(noCorSig)), nrow(noCorSig))
        CC <- Matrix::cov2cor(Sig1)
        ##CCOld <- cov2cor(noCorSig)
        ## Insert
        ##CC[abs(CCOld) > 1e-5] <- CCOld[abs(CCOld) > 1e-5]
        Sig1 <- D %*% CC %*% D
        ## Sig1 <- makeSymPosDef(Sig1) #0.5 * (Sig1 + t(Sig1))
    }
   
    dimnames(Sig1) <- list(names(attr(mFit,"m_opt")$par),names(attr(mFit,"m_opt")$par))
 
    if(!returnSigma){
        Hes1 <- solve(Sig1)
        return(Hes1)
    }
    return(Sig1)
}

retro_hessian_RE <- function(mFit, keep.diagonal = TRUE, returnSigma = TRUE){
    requireNamespace("Matrix")
    oFit <- mFit[[length(mFit)]]
    dataYears <- mFit[[1]]$data$years
    years <- rev(tail(dataYears,length(mFit)))
    m_opt <- attr(mFit,"m_opt")
    m_obj <- attr(mFit,"m_obj")
    ## match parameters to years
    apOrig <- do.call("rbind",lapply(mFit, function(x) sapply(split(x$obj$env$last.par.best[x$obj$env$random],names(x$obj$env$last.par.best[x$obj$env$random])),length)))
    ap <- split(m_obj$env$last.par.best[m_obj$env$random],factor(names(m_obj$env$last.par.best[m_obj$env$random]),unique(names(m_obj$env$last.par.best[m_obj$env$random]))))
    parYear <- lapply(names(ap), function(nm) rep(years,times = apOrig[,nm]))
    names(parYear) <- names(ap)
    parYear <- unlist(parYear)
    ##parYear <- unlist(lapply(ap,function(x) rep(years, each = length(x) / length(years))))
    isFirstYear <- parYear == min(parYear)
    ## Derivative of score function
    ## H <- m_opt$he
    H <- m_obj$env$spHess(m_obj$env$last.par.best, random=TRUE)
    ## Derivative wrt each of the old parameter sets
    Hx <- lapply(head(years,-1), function(y){
        ##hessian_gr(oFit$obj$gr, m_opt$par[parYear==y], method = "central", symmetrize=TRUE)
        ## This is a bit difficult with different # random effects
        p0 <- oFit$obj$env$last.par.best
        np1 <- sum(parYear == y)
        H0 <- Matrix::sparseMatrix(i=1,j=1,x=0,dims = c(np1,np1))
        p1f <- mFit[[match(y,head(years,-1))]]$opt$par
        numPeel <- match(y,head(years,-1))
        ii <- unlist(apply(apOrig,2,function(x){
            rep(c(TRUE,FALSE), times = c(x[length(x)],x[[numPeel]]-x[length(x)]))            
        }, simplify=FALSE))
        p0[oFit$obj$env$lfixed()] <- p1f
        p0[!oFit$obj$env$lfixed()] <- m_obj$env$last.par.best[m_obj$env$random][parYear==y][ii]
        H0[ii,ii] <- oFit$obj$env$spHess(p0, random = TRUE)
        makeSymPosDef(H0)
    })
    info_Hx_peel <- unlist(lapply(head(years,-1), function(y){ rep(max(years)-y, sum(parYear==y)) }))
    info_Hx_par <- unlist(lapply(head(years,-1), function(y){ factor(names(m_obj$env$last.par.best[m_obj$env$random][parYear==y]),unique(names(m_obj$env$last.par.best[m_obj$env$random]))) }))
    info_Hx_num <- unlist(lapply(head(years,-1), function(y){ seq_along(m_obj$env$last.par.best[m_obj$env$random][parYear==y]) }))
    ## Derivative of score function wrt new parameters (ordered by parameter then year)
    J2 <- Matrix::symmpart(as(H[!isFirstYear, !isFirstYear],"sparseMatrix"))
    ## Get derivative and (Gaussian) variance for data
    oldObj <- attr(mFit,"m_obj")
    dat <- oldObj$env$data
    logobsList <- lapply(dat$sam,function(x)x$logobs)
    auxList <- lapply(dat$sam,function(x){
        cbind(x$aux,fleetType = x$fleetTypes[x$aux[,"fleet"]])
    })
    auxDList <- lapply(dat$sam,function(x){
        x$auxData
    })
    logobs <- do.call("c",logobsList)
    stock <- rep(seq_along(logobsList)-1, times = sapply(logobsList,length))    
    indx <- do.call("c",lapply(logobsList,function(x)seq_along(x)-1))
    year <- do.call("c",lapply(auxList,function(x)x[,"year"]))
    dpart <- do.call("c",lapply(auxDList,function(x){
        if(ncol(x) < 5)
            return(rep(0,nrow(x)))
        x[,5]
    }))
    fleetType <- do.call("c",lapply(auxList,function(x)x[,"fleetType"]))
    fleet <- do.call("c",lapply(as.list(seq_along(auxList)),
                                function(i){
                                    prefix <- gsub(" ","0",
                                                   formatC(i,width=nchar(as.character(length(auxList)))))
                                    paste0(prefix,"_",auxList[[i]][,"fleet"])
                                }))
    age <- do.call("c",lapply(auxList,function(x)x[,"age"]))
    xtr <- ifelse(fleetType >= 80, dpart,0)
    kp <- which(!is.na(logobs) & (fleetType < 80))
    sortVals <- order(year[kp],stock[kp],fleet[kp],age[kp], xtr[kp])
    dat$fake_obs <- logobs[kp][sortVals]
    dat$fake_stock <- stock[kp][sortVals]
    dat$fake_indx <- indx[kp][sortVals]
    fake_year <- year[kp][sortVals]
    dat$useFakeObs <- 1
    ran <- c("fake_obs","logN", "logF", "missing","logSW", "logCW", "logitMO", "logNM", "logP","logitFseason", "shared_logFscale","shared_missingObs", "logGst", "logGtrip")
    pl <- c(list(fake_obs = dat$fake_obs), attr(mFit,"m_pl"))
    dat$fake_obs <- NULL
    map <- oldObj$env$map
    map$fake_obs <- factor(ifelse(fake_year <= min(years), NA, seq_along(pl$fake_obs)))
    Obj1 <- TMB::MakeADFun(dat,pl,map, random=ran, DLL="multiStockassessment",inner.control=list(maxit=100))
    isObs <- names(Obj1$env$last.par[Obj1$env$random]) == "fake_obs"
    Obj1$fn()
    ##Hy <- t(hessian_gr(Obj1$gr, Obj1$par, subset = which(isObs), method = HyMethod, symmetrize = TRUE)) ## Subsets on rows, so transform to organize by columns
    Hy <- Obj1$env$spHess(Obj1$env$last.par, random=TRUE)
### The data appears multiple times for some years, use average
    diagA <- max(years) - fake_year[!is.na(map$fake_obs)] + 1
    A <- Matrix::Diagonal(x = sqrt(1/diagA),length(diagA))
    ## Approximate variance of data
    Hy2 <- (as(A %*%Hy[which(isObs),which(isObs)]%*%A,"sparseMatrix"))
    Vy <- (Matrix::solve(Matrix::Cholesky(Hy2)))
    ## J1_2 is ordered by parameter then year
    J1_2 <- Hy[max(which(isObs)) + which(!isFirstYear),which(isObs)] # Derivative wrt new parameter then y
    J1_1 <- do.call("rbind",lapply(seq_along(Hx), function(i){
        x <- Hx[[i]]
        ii <- unlist(apply(apOrig,2,function(v){
            rep(c(TRUE,FALSE), times = c(v[length(v)],v[[i]]-v[length(v)]))            
        }, simplify=FALSE))
        x[ii,ii] <- x[ii,ii]-2*H[isFirstYear,isFirstYear]
        x[,ii]
    }))
    ## Reorder to match J1_2
    j11or <- order(info_Hx_par,info_Hx_peel,info_Hx_num)
    J1_1 <- J1_1[j11or,]
    ## Prepare to reorder G
    info_H_peel <- c(info_Hx_peel[j11or],unlist(lapply(tail(years,1), function(y){ rep(max(years)-y, sum(parYear==y)) })))
    info_H_par <- c(info_Hx_par[j11or],unlist(lapply(tail(years,1), function(y){ factor(names(m_opt$par[parYear==y]),unique(names(m_opt$par))) })))
    info_H_num <- c(info_Hx_num[j11or],unlist(lapply(tail(years,1), function(y){ seq_along(m_opt$par[parYear==y]) })))
    ## Full gradient for Delta Method
    G <- -Matrix::solve(Matrix::Cholesky(J2)) %*% as(cbind(J1_1,J1_2), "sparseMatrix")
    ## G2 <- matrix(0,sum(isFirstYear),ncol(G))
    ## diag(G2[1:sum(isFirstYear),1:sum(isFirstYear)]) <- 1
    G2 <- Matrix::sparseMatrix(i=1:sum(isFirstYear),
                               j=1:sum(isFirstYear),
                               x=0,
                               dims = c(sum(isFirstYear),ncol(G)))
    Gx <- rbind(G,G2)
    ## Need to reorder to match fit
    Gx <- as(Gx[order(info_H_par,info_H_peel,info_H_num),],"sparseMatrix")
    ii <- which(info_H_peel[order(info_H_par,info_H_peel,info_H_num)] == max(info_H_peel))
    Gx[cbind(ii,seq_along(ii))] <- 1   
    ## Delta method to get correlation
    Hold <- (oFit$obj$env$spHess(oFit$obj$env$last.par.best,random=TRUE))
    Vold <- (Matrix::solve(Matrix::Cholesky(Hold)))
    combiVar <- sparse_sym_block(Vold,as(Vy,"sparseMatrix"))
    ## if(missing(subset)){
    Sig1 <- (Gx %*% combiVar) %*% Matrix::t(Gx)
    ## }else{
    ##     Sig1 <- (Gx[subset,,drop=FALSE] %*% combiVar) %*% Matrix::t(Gx)[,subset,drop=FALSE]
    ## }
    ## Symmetrize for safety
    #Sig1 <- makeSymPosDef(Sig1) 
    if(keep.diagonal){
        D <- Matrix::Diagonal(x = sqrt(Matrix::diag(Matrix::solve(Matrix::Cholesky(H)))))
        Do <- Matrix::Diagonal(x = 1.0 / sqrt(Matrix::diag(Sig1)))
        Dx <- D %*% Do
        ## CC <- Do %*% Sig1 %*% Do ##Matrix::cov2cor(Sig1)
        ## Sig1 <- D %*% CC %*% D
        Sig1 <- Dx %*% Sig1 %*% Dx ## Dx is diagonal, so Dx=t(Dx)
        ##Sig1 <- (Sig1) 
    }
   if(!returnSigma){
        Hes1 <- Matrix::solve(Matrix::Cholesky(as(Sig1,"sparseMatrix")))
        return(Hes1)
    }
    Sig1      
}


##' @export
mohn_CI <- function(fit, ...){
    UseMethod("mohn_CI")
}

##' @export
##' NOTE: Add option to work on subset with some loss of precision
mohn_CI.samset <- function(fit, addCorFix = TRUE, addCorRE = TRUE, nosim = 0, ignore.parameter.uncertainty = FALSE, ignore.re.uncertainty = FALSE, ...){
    requireNamespace("Matrix")

    call <- match.call()
    if(is.null(attr(fit,"fit")))
        stop("The samset should have a fit as attribute")
    
    fitFinal <- attr(fit,"fit")
    framework <- .Call("getFramework", PACKAGE=fitFinal$obj$env$DLL)
    if (framework != "TMBad")
        stop("Please re-install the stockassessment and multiStockassessment packages with TMBad")
    ## fitList <- do.call("c",c(list(fitFinal), fit))
    ## retroMS <- multisam.fit(fitList, mohn=1, doSdreport = FALSE)
    retroMS <- make_quick_multisam(fit)

    ## Does not allow for lag in R
    add <- 0
    nms <- c(paste("R(age ", fitFinal$conf$minAge + add, ")", sep = ""),
             "SSB",
             paste("Fbar(", fitFinal$conf$fbarRange[1], "-", fitFinal$conf$fbarRange[2], ")", sep = ""))

    ## Corrected Hessian of fixed effects
    if(!ignore.parameter.uncertainty){
        if(addCorFix){
            Sig0 <- retro_hessian(retroMS)
        }else{
            Hes <- makeSymPosDef(attr(retroMS,"m_opt")$he)
            Sig0 <- (svd_solve(Hes))
        }
    }else{
        Sig0 <- 0 * attr(retroMS,"m_opt")$he
    }
    

    obj <- attr(retroMS,"m_obj")
    par <- obj$env$last.par.best
    r <- obj$env$random
    nonr <- setdiff(seq_along(par), r)
    obj2 <- TMB::MakeADFun(obj$env$data,
                           obj$env$parameters,
                           type = "ADFun",
                           ADreport = TRUE,
                           DLL = obj$env$DLL,
                           silent = obj$env$silent)
    ## Vtheta <- Sig0
    ## hessian.random <- Hes_uu
    phi <- try(obj2$fn(par), silent=TRUE)    ## NOTE_1: obj2 forward sweep now initialized !
    if(is.character(phi) | length(phi)==0){
        stop("Error in phi")
        phi <- rep(NA,6)
    }
    ADGradForward0Initialized <- FALSE
    ADGradForward0Initialize <- function() { ## NOTE_2: ADGrad forward sweep now initialized !
        obj$env$f(par, order = 0, type = "ADGrad")
        ADGradForward0Initialized <<- TRUE
    }
    chunk <- unlist(obj$env$ADreportIndex()[c("mohnRho_rec","mohnRho_ssb","mohnRho_fbar","mohnRhoMod_rec","mohnRhoMod_ssb","mohnRhoMod_fbar")])
    w <- rep(0, length(phi))
    phiDeriv <- function(i){
        w[i] <- 1
        obj2$env$f(par, order=1, rangeweight=w, doforward=0) ## See NOTE_1
    }
    Dphi <- t( sapply(chunk, phiDeriv) )
    phi <- phi[chunk]      
    inUse <- which(apply(Dphi!=0,2,any))
    inUseR <- match(intersect(inUse,r),r)

    ## Corrected Hessian of random effects
    if(!ignore.re.uncertainty){
        if(addCorRE){
            Sig_uu_tmp <- retro_hessian_RE(retroMS)#, subset = inUseR)
            Sig_uu <- as(Sig_uu_tmp,"sparseMatrix") ## Matrix::symmpart(Sig_uu_tmp)
            ## Sig_Chol_uu <- Matrix::Cholesky(Sig_uu)
            ## Hes_uu <- Matrix::symmpart(Matrix::solve(Sig_Chol_uu))
        }else{
            Hes_uu <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
            Sig_uu <- as(Matrix::solve(Matrix::Cholesky(Hes_uu)),"sparseMatrix") #[inUseR,inUseR, drop = FALSE]
        }
    }else{
        Sig_uu <- 0 * (obj$env$spHess(obj$env$last.par.best, random = TRUE))#[inUseR,inUseR, drop = FALSE] * 0
    }
    dimnames(Sig_uu) <- list(names(attr(retroMS,"m_obj")$env$last.par.best[r]),#[inUseR]),
                             names(attr(retroMS,"m_obj")$env$last.par.best[r]))#[inUseR]))
    
    if(nosim == 0){ ## Delta method (calculated)
        Dphi.random <- as(Dphi[,r,drop=FALSE],"sparseMatrix")#[,intersect(inUse,r),drop=FALSE],"sparseMatrix")
        Dphi.fixed <- as(Dphi[,-r,drop=FALSE],"sparseMatrix")        
        ##tmp <- Matrix::solve(Hes_uu,Matrix::t(Dphi.random))
        tmp <- Sig_uu %*% Matrix::t(Dphi.random)
        ##tmp <- as.matrix(tmp)
        if(ignore.re.uncertainty){
            term1 <- matrix(0,nrow(Dphi.random),nrow(Dphi.random))
        }else{
            term1 <- Dphi.random%*%tmp ## first term.
        }
        if(ignore.parameter.uncertainty){
            term2 <- 0
        }else{
            w <- rep(0, length(par))
            if(!ADGradForward0Initialized) ADGradForward0Initialize()
            reverse.sweep <- function(i){
                ##w[intersect(inUse,r)] <- tmp[,i]
                w[r] <- tmp[,i]
                -obj$env$f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[-r]
            }
            A <- t(do.call("cbind",lapply(seq_along(phi), reverse.sweep))) + Dphi.fixed
            term2 <- (A %*% (Sig0 %*% Matrix::t(A))) ## second term
        }
        cov <- term1 + term2
        ## End of modified from TMB
        sdv <- sqrt(Matrix::diag(cov))
        tab <- list(Original = cbind("Estimate" = phi[1:3],
                                     "CI_low" = phi[1:3] - 2 * sdv[1:3],
                                     "CI_high" = phi[1:3] + 2 * sdv[1:3]),
                    Modified = cbind("Estimate" = phi[4:6],
                                     "CI_low" = phi[4:6] - 2 * sdv[4:6],
                                     "CI_high" = phi[4:6] + 2 * sdv[4:6]))
        bginfo <- list(phi = phi, cov = cov)
    }else{ ## Delta method (simulated)
        ## Calculate full precision (modified from TMB::sdreport)
         ## NOTE: TMB usually works with the Hessians, but the correlation is added on the Covariances, so it is easier to work on those
        ## tmp <- obj$env$f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=r) 
        ## A <- solve(Hes_uu, tmp)
        ## G <- Hes_uu %*% A
        ## G <- as.matrix(G) ## Avoid Matrix::cbind2('dsCMatrix','dgeMatrix')
        ## M1 <- cbind2(Hes_uu,G)
        ## M2 <- cbind2(t(G), as.matrix(t(A)%*%G)+Hes)
        ## M <- rbind2(M1,M2)
        ## M <- forceSymmetric(M,uplo="L")
        ## dn <- c(names(par)[r],names(par[-r]))
        ## dimnames(M) <- list(dn,dn)
        ## p <- invPerm(c(r,(1:length(par))[-r]))
        ## jointPrecision <- M[p,p]
        ## Find parameters/RE with influence
        if(ignore.parameter.uncertainty){
            ## Ch <- Matrix::Cholesky(Sig_uu) ## Already reduced
            ## ee <- Matrix::expand2(Ch, LDL = FALSE)
            ## if(!is.null(ee$P1)){
            ##     C0 <- ee$`P1.` %*% ee$L
            ## }else{
            ##     C0 <- ee$L
            ## }
            C0 <- try({Matrix::t(Matrix::chol(Sig_uu[inUseR,inUseR],pivot=FALSE))})
            if(is(C0,"try-error")){
                a <- svd(Sig_uu[inUseR,inUseR])
                C0 <- a$u %*% diag(sqrt(pmax(a$d,0)),length(a$d))
            }
        }else{
            ## f <- obj$env$f
            ## w <- rep(0, length(par))
            ## if(!ADGradForward0Initialized) ADGradForward0Initialize()
            ## tmp <- f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=intersect(inUse,r)) ## TMBad only !!!
            ## if(!is.matrix(tmp)) ## Happens if length(r)==1
            ##     tmp <- matrix(tmp, ncol=length(nonr) )
            ## J <- -Sig_uu %*% tmp
            ## A <- Sig_uu %*% tmp #solve(hessian.random, tmp)
            ## diag.term2 <- rowSums((A %*% Sig0)*A)
            ## hessian.random <- solve(Sig_uu)
            ## hessian.fixed <- solve(Sig0)
            ## G <- hessian.random %*% A
            ## G <- as.matrix(G) ## Avoid Matrix::cbind2('dsCMatrix','dgeMatrix')
            ## M1 <- cbind2(hessian.random,G)
            ## M2 <- cbind2(t(G), as.matrix(t(A)%*%G)+hessian.fixed )
            ## M <- rbind2(M1,M2)
            ## M <- forceSymmetric(M,uplo="L")
            ## dn <- c(names(par)[r],names(par[-r]))
            ## dimnames(M) <- list(dn,dn)
            ## p <- invPerm(c(r,(1:length(par))[-r]))
            ## jointPrecision <- M[p,p]
            ## Vfull <- svd_solve(jointPrecision)
            g_u_th <- obj$env$f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=r)
            J <- -Sig_uu %*% g_u_th
            V2 <- (J %*% Sig0 %*% Matrix::t(J))[inUseR,inUseR]
            ## V2 <- makeSymPosDef(V2)
            ## Vfull <- as(V2,"sparseMatrix")
            if(!ignore.re.uncertainty){
                #u2r <- match(intersect(inUse,r),r)
                ## ru <- r #which(inUse %in% intersect(inUse,r))
                ## Sig_uuUse <- as(Sig_uu,"generalMatrix") ## Already reduced
                ## iuu <- ru[Sig_uuUse@i+1]
                ## puu <- Sig_uuUse@p
                ## dp <- diff(puu)
                ## juu <- ru[rep(seq_along(dp),dp)]
                ## xuu <- Sig_uuUse@x
                ## isU <- iuu <= juu
                ## Vx <- Matrix::sparseMatrix(i=iuu[isU],j=juu[isU],x=xuu[isU],dims=c(length(inUse),length(inUse)), symmetric=TRUE)
                Vfull <- as(V2,"sparseMatrix") + Sig_uu[inUseR,inUseR]
            }else{
                Vfull <- V2
            }
            ## Cholesky
            ## ChVF <- Matrix::Cholesky(Vfull)
            ## Lower triangular
            ## ee <- Matrix::expand2(ChVF, LDL = FALSE)
            ##  if(!is.null(ee$P1)){
            ##     C0 <- ee$`P1.` %*% ee$L
            ## }else{
            ##     C0 <- ee$L
            ## }
            C0 <- try({Matrix::t(Matrix::chol(Vfull,pivot=FALSE))})
            if(is(C0,"try-error")){
                a <- svd(Vfull)
                C0 <- a$u %*% diag(sqrt(pmax(a$d,0)),length(a$d))
            }
        }
        ## Simulation procedure
        doOne0 <- function(sim=TRUE){
            p1 <- obj$env$last.par.best
            if(sim){
                ## if(ignore.parameter.uncertainty){
                ##     p1[r] <- p1[r] + drop(C0 %*% rnorm(length(r)))
                ## }else{
                ##     p1 <- p1 + drop(C0 %*% rnorm(length(p1)))
                ## }
                p1[intersect(inUse,r)] <- p1[intersect(inUse,r)] + as.vector(C0 %*% rnorm(length(intersect(inUse,r))))
            }
            rp0 <- obj$report(p1)
            vMod <- c(R = mean(apply((rp0$mohnRhoVec_rec),1,function(x)(x[1]-x[2])/log(10))),
                      SSB = mean(apply((rp0$mohnRhoVec_ssb),1,function(x)(x[1]-x[2])/log(10))),
                      Fbar = mean(apply((rp0$mohnRhoVec_fbar),1,function(x)(x[1]-x[2])/log(10))))
            vOld <- c(R = mean(apply(rp0$mohnRhoVec_rec,1,function(x)exp(x[1]-x[2])-1)),
                      SSB = mean(apply(rp0$mohnRhoVec_ssb,1,function(x)exp(x[1]-x[2])-1)),
                      Fbar = mean(apply(rp0$mohnRhoVec_fbar,1,function(x)exp(x[1]-x[2])-1)))
            list(Modified = vMod, Original = vOld)
        }
        estRho <- doOne0(FALSE)
        estRho_Orig <- estRho$Original
        estRho_Mod <- estRho$Modified
        simRho <- replicate(nosim,tryCatch(doOne0(TRUE), error = function(e) list(Modified = c(Fbar=NA,SSB=NA,R=NA), Original = c(Fbar=NA,SSB=NA,R=NA))), simplify = FALSE)
        simRho_Orig <- do.call("cbind",lapply(simRho, function(x) x$Original))
        simRho_Mod <- do.call("cbind",lapply(simRho, function(x) x$Modified))
        bginfo_Orig <- list(est = estRho_Orig, sim = simRho_Orig)
        bginfo_Mod <- list(est = estRho_Mod, sim = simRho_Mod)
        cil_Orig <- apply(simRho_Orig,1,quantile,prob=0.025, na.rm=TRUE)
        cih_Orig <- apply(simRho_Orig,1,quantile,prob=0.975, na.rm=TRUE)
        cil_Mod <- apply(simRho_Mod,1,quantile,prob=0.025, na.rm=TRUE)
        cih_Mod <- apply(simRho_Mod,1,quantile,prob=0.975, na.rm=TRUE)
        tab <- list(Original = cbind(Estimate = estRho_Orig,
                                     CI_low = cil_Orig,
                                     CI_high = cih_Orig),
                    Modified = cbind(Estimate = estRho_Mod,
                                     CI_low = cil_Mod,
                                     CI_high = cih_Mod))
        bginfo <- list(Original = bginfo_Orig,
                       Modified = bginfo_Mod)              
    }

    rownames(tab$Original) <- nms
    rownames(tab$Modified) <- nms

    res <- list(table = tab,
                call = call,
                bginfo = bginfo,
                mfit = retroMS)
    class(res) <- "sam_mohn"
    res    
}


##' @export
mohn_CI.sam <- function(fit, addCorrelation = TRUE, years = 5, ...){
    ## Fake retro
    getRetroInputStocks <- function(fit, n){
        f <- runwithout(fit,year=tail(fit$data$years,n),run=FALSE)
        f$opt <- list(objective=NA)
        class(f) <- "sam"
        f
    }
    fitList <- do.call("c",lapply(seq_len(years), getRetroInputStocks,fit=fit))
    attr(fitList,"fit") <- fit
    mohn_CI.samset(fitList, addCorrelation = addCorrelation, ...)
}

##' @export
print.sam_mohn <- function(x, ...){
    print(x$table)
    invisible(x)
}



##' @export
mohn_sim_CI <- function(fit, ...){
    UseMethod("mohn_sim_CI")
}


##' @export
mohn_sim_CI.sam <- function(fit, nsim, type = c("Full","Gauss","GaussF","Tail","FixF","Obs"), resampleParameters=FALSE, ncores = 1, year = 5){
    RETRO <- retro(fit, year = year, ncores = ncores)
    mohn_sim_CI(RETRO, nsim = nsim, type = type, resampleParameters = resampleParameters, ncores = ncores)
    
}

##' @export
mohn_sim_CI.samset <- function(fit, nsim, type = c("Full","Gauss","GaussF","Tail","FixF","Obs"), resampleParameters=FALSE, ncores = 1){

    mohn_2 <- function(fits, what=NULL, lag=0, modified = FALSE, ...){
        if(is.null(what)){
            what <- function(fit){
                ret <- cbind(rectable(fit,...)[,1], ssbtable(fit)[,1], fbartable(fit)[,1])
                add <- 0
                dots <- list(...)
                if(!is.null(dots$lagR)){
                    if(dots$lagR == TRUE){
                        add <- 1
                    }
                }
                colnames(ret) <- c(paste("R(age ", fit$conf$minAge + add, ")", sep = ""), "SSB",
                                   paste("Fbar(", fit$conf$fbarRange[1], "-", fit$conf$fbarRange[2], ")", sep = ""))
                ret
            }
        }
        ref <- what(attr(fits,"fit"))
        ret <- lapply(fits, what)
        if(modified){
            bias <- lapply(ret, function(x){y<-rownames(x)[nrow(x)-lag]; (log10(x)[rownames(x)==y,]-log10(ref)[rownames(ref)==y,])})
        }else{
            bias <- lapply(ret, function(x){y<-rownames(x)[nrow(x)-lag]; (x[rownames(x)==y,]-ref[rownames(ref)==y,])/ref[rownames(ref)==y,]})
        }
        colMeans(do.call(rbind,bias))
    }
    
    RETRO <- fit
    def <- mohn_2(RETRO)
    def[] <- NA
    def <- list(Original = def, Modified = def)
    if(type == "Tail"){
        doOne <- function(){
            simfit <- try({stockassessment:::addSimulatedYears(RETRO[[length(RETRO)]],rep(NA,length(RETRO)), resampleFirst = TRUE, refit = TRUE, resampleParameters = resampleParameters)})
            if(is(simfit,"try-error"))
                return(def)
            re_retro <- try({stockassessment::retro(simfit, length(RETRO), ncores=ncores)},silent=TRUE)
            if(is(re_retro,"try-error"))
                return(def)            
            list(Original = mohn_2(re_retro), Modified = mohn_2(re_retro, modified = TRUE))
        }
    }else{
        fit <- attr(RETRO,"fit")
        obj <- unserialize(serialize(fit$obj, NULL))
        ## Fix in simulation:
#### F (all but Full and Tail),
#### N (Fix for Gauss and Obs),
#### Obs (None)
        obj$env$data$simFlag <- c(ifelse(type %in% c("Full","Tail"),0,1),
                                  ifelse(type %in% c("Gauss","Obs"),1,0),0)
        obj$retape()
        
        obj$fn(fit$opt$par)
        pl <- fit$pl
        map <- fit$obj$env$map
        if(type %in% c("GaussF","Gauss")){
            C0 <- solve(obj$env$spHess(obj$env$last.par.best,random=TRUE))  ## Covariance matrix of random effects
        }
        doOne <- function(){
            ## Handle F (/N) for different cases
            ## Full: Let $simulate handle everything
            ## FixF: Let $simulate handle everything
            ## GaussF: Posterior simulation of F (simulate N and overwrite in $simulate)
            ## Gauss: Posterior simulation of F and N
            if(resampleParameters){
                
            }
            if(type %in% c("GaussF","Gauss")){
                if(resampleParameters){
                    px <- stockassessment:::rmvnorm(1, fit$opt$par, solve(fit$opt$he), TRUE)
                    obj$fn(px)
                    p0 <- obj$env$last.par
                }else{
                    p0 <- obj$env$last.par.best
                }
                p0[obj$env$random] <- stockassessment:::rmvnorm(1,p0[obj$env$random],C0, pivot = TRUE)
                pl <- obj$env$parList(par = p0)
            }else if(resampleParameters){
                px <- stockassessment:::rmvnorm(1, fit$opt$par, solve(fit$opt$he), TRUE)
                obj$fn(px)
                pl <- obj$env$parList(par = obj$env$last.par)
            }
            ##
            with.map <- intersect(names(pl), names(map))
            applyMap <- function(par.name) {
                tapply(pl[[par.name]], map[[par.name]], mean)
            }
            plOrig <- pl
            pl[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
            est <- unlist(pl)
            sval <- obj$simulate(est)
            ret <- fit$data
            nms <- intersect(names(ret), names(sval))
            ret[nms] <- sval[nms]
            ## New fit
            plOrig$missing <- NULL
            attr(plOrig,"what") <- NULL
            simfit <- try({suppressWarnings(sam.fit(ret, fit$conf, plOrig, map = map, silent = TRUE))},silent=TRUE)
            if(is(simfit,"try-error"))
                return(def)
            re_retro <- try({stockassessment::retro(simfit, length(RETRO), ncores=ncores)},silent=TRUE)
            if(is(re_retro,"try-error"))
                return(def)
            list(Original = mohn_2(re_retro), Modified = mohn_2(re_retro, modified = TRUE))
        }        
    }
    replicate(nsim, doOne())
}
