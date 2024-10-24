
block <- function(...){
    x <- list(...)
    ii <- sapply(x, is.matrix)
    x[!ii] <- lapply(x[!ii],function(y)diag(y,length(y),length(y)))
    n <- sapply(x,nrow)
    M <- matrix(0,sum(n),sum(n))
    nc <- cumsum(c(0,n))
    for(i in 1:length(x)){
        jj <- nc[i]+1:n[i]
        M[jj,jj] <- x[[i]]
    }
    M    
}

get_triplet <- function(x){
    y <- as(x,"TsparseMatrix")
    r <- data.frame(i = y@i+1,
                    j = y@j+1,
                    x = y@x)
    attr(r,"Dims") <- y@Dim
    r
}

sparse_block <- function(...){
    x <- list(...)
    ii <- sapply(x, is, class2 = "sparseMatrix")
    x[!ii] <- lapply(x[!ii],function(y)Matrix::sparseMatrix(i=seq_along(y),j=seq_along(y),x=as.numeric(y)))
    y <- lapply(x, get_triplet)
    n <- sapply(x,nrow)
    nc <- cumsum(n)
    for(i in 2:length(y)){
        y[[i]]$i <- y[[i]]$i + nc[i-1]
        y[[i]]$j <- y[[i]]$j + nc[i-1]
    }
    yv <- do.call("rbind",y)
    Matrix::sparseMatrix(i=yv$i, j=yv$j, x=yv$x, dims = c(sum(n),sum(n)))
}

makeSymPosDef <- function(x){
    ## Sym
    v <- x + t(x)
    ee <- eigen(Sig1, symmetric = TRUE)
        ee$values <- pmax(ee$values,1e-6 / max(ee$values))
        Sig1 <- ee$vectors %*% diag(x=ee$values) %*% solve(ee$vectors)
        Sig1 <- 0.5 * (Sig1 + t(Sig1))

}

retro_hessian <- function(mFit, keep.diagonal = TRUE, HyMethod = "forward", forcePosDef = TRUE, returnSigma = FALSE){
##### Calculate hessian with correlation #####
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
    Vy <- stockassessment:::svd_solve(A %*%Hy[which(isObs),which(isObs)]%*%A)
    ## Symmetrize for safety
    Vy <- 0.5 * (Vy + t(Vy))
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
    ## Full gradient for Delta Method
    G <- -stockassessment:::svd_solve(J2) %*% (cbind(J1_1,J1_2) )
    G2 <- matrix(0,sum(isFirstYear),ncol(G))
    diag(G2[1:sum(isFirstYear),1:sum(isFirstYear)]) <- 1
    Gx <- rbind(G,G2)
    ## Need to reorder to match fit
    Gx <- Gx[order(info_H_par,info_H_peel,info_H_num),]
    ## Delta method to get correlation
    Sig1 <- Gx %*% block(stockassessment:::svd_solve(oFit$opt$he),Vy) %*% t(Gx)
    ## Symmetrize for safety
    Sig1 <- 0.5 * (Sig1 + t(Sig1))    
    if(forcePosDef){
        ## ss <- Matrix::Schur(Sig1)
        ## ss@EValues <- pmax(ss@EValues,1e-6 / max(ss@EValues))
        ## Sig1 <- as(Reduce(`%*%`, Matrix::expand2(ss)), "sparseMatrix")
        ## Sig1 <- Matrix::symmpart(Sig1)
        ee <- eigen(Sig1, symmetric = TRUE)
        ee$values <- pmax(ee$values,1e-6 / max(ee$values))
        Sig1 <- ee$vectors %*% diag(x=ee$values) %*% solve(ee$vectors)
        Sig1 <- 0.5 * (Sig1 + t(Sig1))
    }
    if(keep.diagonal){
        D <- diag(sqrt(diag(svd_solve(m_opt$he))))
        CC <- cov2cor(Sig1)
        Sig1 <- D %*% CC %*% D
        Sig1 <- 0.5 * (Sig1 + t(Sig1))
    }
    
    if(!returnSigma){
        Hes1 <- solve(Sig1)
        return(Hes1)
    }
    return(Sig1)
}

retro_hessian_RE <- function(mFit, keep.diagonal = FALSE, forcePosDef = TRUE, returnSigma = TRUE){
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
        H0
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
    Hy2 <- Matrix::symmpart(as(A %*%Hy[which(isObs),which(isObs)]%*%A,"sparseMatrix"))
    Vy <- Matrix::solve(Matrix::Cholesky(Hy2))
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
                               x=rep(1,sum(isFirstYear)),
                               dims = c(sum(isFirstYear),ncol(G)))
    Gx <- rbind(G,G2)
    ## Need to reorder to match fit
    Gx <- as(Gx[order(info_H_par,info_H_peel,info_H_num),],"sparseMatrix")
    ## Delta method to get correlation
    Hold <- Matrix::symmpart(oFit$obj$env$spHess(oFit$obj$env$last.par.best,random=TRUE))
    Vold <- Matrix::solve(Matrix::Cholesky(Hold))
    combiVar <- sparse_block(Vold,as(Vy,"sparseMatrix"))
    Sig1 <- (Gx %*% combiVar) %*% Matrix::t(Gx)
    ## Symmetrize for safety
    Sig1 <- Matrix::symmpart(Sig1) 
     if(forcePosDef){
        ##ee <- eigen(Sig1, symmetric = TRUE)
        ss <- Matrix::Schur(Sig1)
        ss@T@x[] <- ss@EValues[] <- pmax(ss@T@x,1e-6 / max(ss@T@x))
        Sig1 <- as(Reduce(`%*%`, Matrix::expand2(ss)), "sparseMatrix")
        ##Sig1 <- Matrix::symmpart(Sig1)
    }
    if(keep.diagonal){
        D <- Matrix::Diagonal(x = sqrt(Matrix::diag(Matrix::solve(Matrix::Cholesky(H)))))
        Do <- Matrix::Diagonal(x = 1.0 / sqrt(Matrix::diag(Sig1)))
        Dx <- D %*% Do
        ## CC <- Do %*% Sig1 %*% Do ##Matrix::cov2cor(Sig1)
        ## Sig1 <- D %*% CC %*% D
        Sig1 <- Dx %*% Sig1 %*% Dx
        Sig1 <- Matrix::symmpart(Sig1) 
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
mohn_CI.samset <- function(fit, addCorFix = TRUE, addCorRE = TRUE, nosim = 0, ignore.parameter.uncertainty = FALSE, ...){

    call <- match.call()
    if(is.null(attr(fit,"fit")))
        stop("The samset should have a fit as attribute")
    
    fitFinal <- attr(fit,"fit")
    framework <- .Call("getFramework", PACKAGE=fitFinal$obj$env$DLL)
    if (framework != "TMBad")
        stop("Please re-install the stockassessment and multiStockassessment packages with TMBad")
    fitList <- do.call("c",c(list(fitFinal), fit))
    retroMS <- multisam.fit(fitList, mohn=1, doSdreport = FALSE)

    ## Does not allow for lag in R
    add <- 0
    nms <- c(paste("R(age ", fitFinal$conf$minAge + add, ")", sep = ""),
             "SSB",
             paste("Fbar(", fitFinal$conf$fbarRange[1], "-", fitFinal$conf$fbarRange[2], ")", sep = ""))

    ## Corrected Hessian of fixed effects
    if(!ignore.parameter.uncertainty){
        if(addCorFix){
            Sig0_tmp <- retro_hessian(retroMS, returnSigma = TRUE,  keep.diagonal = TRUE, forcePosDef = TRUE)
            Sig0 <- Matrix::symmpart(Sig0_tmp)
            ## Sig0_Chol <- Matrix::Cholesky(as(Sig0,"sparseMatrix"))
            ## Hes <- Matrix::symmpart(Matrix::solve(Sig0_Chol))
        }else{
            Hes <- attr(retroMS,"m_opt")$he
            Sig0 <- Matrix::symmpart(solve(Hes))
        }
    }else{
        Sig0 <- 0 * attr(retroMS,"m_opt")$he
    }
    
    ## Corrected Hessian of random effects
    if(addCorRE){
        Sig_uu_tmp <- retro_hessian_RE(retroMS, returnSigma = TRUE, keep.diagonal = TRUE, forcePosDef = TRUE)
        Sig_uu <- as(Sig_uu_tmp,"sparseMatrix") ## Matrix::symmpart(Sig_uu_tmp)
        ## Sig_Chol_uu <- Matrix::Cholesky(Sig_uu)
        ## Hes_uu <- Matrix::symmpart(Matrix::solve(Sig_Chol_uu))
    }else{
        obj0 <- attr(retroMS,"m_obj")
        Hes_uu <- Matrix::symmpart(obj0$env$spHess(obj0$env$last.par.best, random = TRUE))
        Sig_uu <- as(Matrix::symmpart(Matrix::solve(Matrix::Cholesky(Hes_uu))),"sparseMatrix")
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
    phi <- obj2$fn(par)
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
    if(nosim == 0){ ## Delta method (calculated)
        phi <- phi[chunk]
        Dphi.random <- as(Dphi[,r,drop=FALSE],"sparseMatrix")
        Dphi.fixed <- as(Dphi[,-r,drop=FALSE],"sparseMatrix")        
        ##tmp <- Matrix::solve(Hes_uu,Matrix::t(Dphi.random))
        tmp <- Sig_uu %*% Matrix::t(Dphi.random)
        ##tmp <- as.matrix(tmp)
        term1 <- Dphi.random%*%tmp ## first term.
        if(ignore.parameter.uncertainty){
            term2 <- 0
        }else{
            w <- rep(0, length(par))
            if(!ADGradForward0Initialized) ADGradForward0Initialize()
            reverse.sweep <- function(i){
                w[r] <- tmp[,i]
                -obj$env$f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[-r]
            }
            A <- t(do.call("cbind",lapply(seq_along(phi), reverse.sweep))) + Dphi.fixed
            term2 <- A %*% (Sig0 %*% Matrix::t(A)) ## second term
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
        inUse <- which(apply(Dphi!=0,2,any))
        if(ignore.parameter.uncertainty){
            inUse <- intersect(inUse,r)
            Ch <- Matrix::Cholesky(Sig_uu[match(inUse,r),match(inUse,r),drop=FALSE])
            C0 <- Matrix::expand2(Ch, LDL = FALSE)$`L`
        }else{
            J <- obj$env$f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=inUse)
            if(length(intersect(inUse,nonr)) > 0){
                isFc <- match(intersect(inUse,nonr),inUse)
                isFr <- match(intersect(inUse,nonr),nonr)
                J[cbind(isFr,isFc)] <- 1
            }
            V2 <- J %*% Sig0 %*% t(J)
            Vfull <- as(V2,"sparseMatrix")
            u2r <- match(intersect(inUse,r),r)
            ru <- which(inUse %in% intersect(inUse,r))
            Sig_uuUse <- as(Sig_uu[u2r,u2r,drop=FALSE],"generalMatrix")
            iuu <- ru[Sig_uuUse@i+1]
            puu <- Sig_uuUse@p
            dp <- diff(puu)
            juu <- ru[rep(seq_along(dp),dp)]
            xuu <- Sig_uuUse@x
            Vx <- Matrix::sparseMatrix(i=iuu,j=juu,x=xuu,dims=c(length(inUse),length(inUse)))
            Vfull <- Vfull + Vx
            ## Cholesky
            ChVF <- Matrix::Cholesky(Vfull)
            ## Lower triangular
            C0 <- Matrix::expand2(ChVF, LDL = FALSE)$`L`
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
                p1[inUse] <- p1[inUse] + as.vector(C0 %*% rnorm(length(inUse)))
            }
            rp0 <- obj$report(p1)
            vMod <- c(R = mean(apply((rp0$mohnRhoVec_rec),1,function(x)(x[1]-x[2])/log(10))),
                      SSB = mean(apply((rp0$mohnRhoVec_ssb),1,function(x)(x[1]-x[2])/log(10))),
                      Fbar = mean(apply((rp0$mohnRhoVec_fbar),1,function(x)(x[1]-x[2])/log(10))))
            vOld <- c(R = mean(apply(exp(rp0$mohnRhoVec_rec),1,function(x)x[1]/x[2]-1)),
                      SSB = mean(apply(exp(rp0$mohnRhoVec_ssb),1,function(x)x[1]/x[2]-1)),
                      Fbar = mean(apply(exp(rp0$mohnRhoVec_fbar),1,function(x)x[1]/x[2]-1)))
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


mohn_CI_Old <- function(fit, addCorrelation = TRUE, simDelta = 0, quantile_CI=FALSE, resampleRE = FALSE, onlyRE = FALSE, ...){
    ## Already fitted with retro
    call <- match.call()
    if(is.null(attr(fit,"fit")))
        stop("The samset should have a fit as attribute")
    if(simDelta == 0 & onlyRE)
        stop("onlyRE is only for simulation")
    fitFinal <- attr(fit,"fit")
    fitList <- do.call("c",c(list(fitFinal), fit))
    retroMS <- multisam.fit(fitList, mohn=1, doSdreport = FALSE)

    ## Does not allow for lag in R
    add <- 0
    nms <- c(paste("R(age ", fitFinal$conf$minAge + add, ")", sep = ""),
             "SSB",
             paste("Fbar(", fitFinal$conf$fbarRange[1], "-", fitFinal$conf$fbarRange[2], ")", sep = ""))

    
    if(addCorrelation){
        if(!onlyRE){
            Sig0_tmp <- retro_hessian(retroMS, returnSigma = TRUE,  keep.diagonal = TRUE, forcePosDef = TRUE)
            Sig0 <- Matrix::symmpart(Sig0_tmp)
            Sig0_Chol <- Matrix::Cholesky(as(Sig0,"sparseMatrix"))
            Hes <- Matrix::symmpart(Matrix::solve(Sig0_Chol))
        }
        if(resampleRE | onlyRE){
            Sig_uu_tmp <- retro_hessian_RE(retroMS, returnSigma = TRUE, keep.diagonal = TRUE, forcePosDef = TRUE)
            Sig_uu <- as(Sig_uu_tmp,"sparseMatrix") ## Matrix::symmpart(Sig_uu_tmp)
            Sig_Chol_uu <- Matrix::Cholesky(Sig_uu)
            Hes_uu <- Matrix::symmpart(Matrix::solve(Sig_Chol_uu))
        }
    }else{
        if(!onlyRE){
            Hes <- attr(retroMS,"m_opt")$he
            Sig0 <- Matrix::symmpart(solve(Hes))
        }
        if(resampleRE | onlyRE){
            obj0 <- attr(retroMS,"m_obj")
            Hes_uu <- Matrix::symmpart(obj0$env$spHess(obj0$env$last.par.best, random = TRUE))
            Sig_uu <- as(Matrix::symmpart(Matrix::solve(Matrix::Cholesky(Hes_uu))),"sparseMatrix")
            Sig_Chol_uu <- Matrix::Cholesky(Sig_uu)
        }
    }
    
    if(simDelta > 0){
        opt0 <- attr(retroMS,"m_opt")
        obj0 <- attr(retroMS,"m_obj")
        
        ## check rank
        ## rnk <- qr(Sig0,LAPACK=TRUE)$rank
        ## if(rnk < nrow(Sig0)){
        ##     ee <- eigen(Sig0, symmetric = TRUE)
        ##     Sig0 <- ee$vectors %*% diag(pmax(ee$values,1e-6 / max(ee$values))) %*% solve(ee$vectors)
        ##     Sig0 <- 0.5 * (Sig1 + t(Sig1))
        ## }
        ## Get error now if it won't make chol
        ## C0 <- chol(Sig0)
        isFixed <- obj0$env$lfixed()

        ## Huu <- retro_hessian_RE(retroMS)
        
        ## if(resampleRE){

        ## tmp <- obj0$env$f(obj0$env$last.par.best, order = 1, type = "ADGrad", 
        ##                   keepx = which(isFixed), keepy = which(!isFixed))
        ## A <- Matrix::solve(Hes_uu, tmp) # Sig_uu %*% tmp #
        ## G <- Hes_uu %*% A
        ## G <- Matrix::as.matrix(G)
        ## M1 <- Matrix::cbind2(Hes_uu, G)
        ## M2 <- Matrix::cbind2(Matrix::t(G), Matrix::as.matrix(Matrix::t(A) %*% G) + Hes)
        ## M <- Matrix::rbind2(M1, M2)
        ## M <- Matrix::forceSymmetric(M, uplo = "L")
        ## dn <- c(names(obj0$env$last.par.best)[!isFixed], names(obj0$env$last.par.best)[isFixed])
        ## dimnames(M) <- list(dn, dn)
        ## p <- Matrix::invPerm(c(which(!isFixed), (1:length(obj0$env$last.par.best))[which(isFixed)]))
        ## jointPrecision <- as(M[p, p],"sparseMatrix")
        ## jp_Chol <- Matrix::Cholesky(jointPrecision)
        ## jointCov <- Matrix::solve(jp_Chol)
        ## C0 <- Matrix::chol(jointCov)
        if(!onlyRE)
            C0f <- Matrix::chol(Sig0)
        if(resampleRE | onlyRE){
            C0u <- Matrix::expand2(Sig_Chol_uu, LDL = FALSE)$`L.`
            ##C0u <- Matrix::chol(Sig_uu)
        }
        #C0 <- Matrix::solve(pC0)
        ##jointCov <- Matrix::solve(jointPrecision)
        ## Siguu <- jointCov[!isFixed,!isFixed] - jointCov[!isFixed,isFixed] %*% Matrix::solve(jointCov[isFixed,isFixed]) %*% jointCov[isFixed,!isFixed]
        ## Cuu <- Matrix::chol(Siguu)
        ## }
        #C0 <- as(C0,"sparseMatrix")
        doOne0 <- function(sim=TRUE){
            ## if(sim & resampleRE){
            ##     cat("HELLO\n")
            ##     p1 <- obj0$env$last.par.best
            ##     p1[!isFixed] <- p1[!isFixed] + drop(rnorm(sum(!isFixed)) %*% Cuu)
            ## }else{
            if(sim){
                ## p0 <- stockassessment:::rmvnorm(1, opt0$par, Sig0)
                if(!onlyRE){
                    p0 <- opt0$par + drop(rnorm(length(opt0$par)) %*% C0f)
                    obj0$fn(p0)
                    p1 <- obj0$env$last.par #.best
                }else{
                    p1 <- obj0$env$last.par.best
                }
                ## ## Respect bounds
                ## lower2 <- attr(retroMS,"m_low")
                ## upper2 <- attr(retroMS,"m_high")
                ## p0 <- pmax(pmin(upper2, p0), lower2)
                ## ## bound transfIRARdist further
                ## p0[names(opt0$par)=="transfIRARdist"] <- pmin(pmax(p0[names(opt0$par)=="transfIRARdist"],-5),5)
                
                ## }else{
                ##     p0 <- opt0$par
                ## }
                
                ##p1 <- obj0$env$last.par.best + drop(rnorm(length(obj0$env$last.par.best))) %*% C0
                
                if(resampleRE | onlyRE)
                    p1[!isFixed] <- p1[!isFixed] + drop(rnorm(sum(!isFixed)) %*% C0u)
            }else{
                p1 <- obj0$env$last.par.best
            }
            ## if(sim & resampleRE){
            ##     p1[!isFixed] <- p1[!isFixed] + drop(rnorm(sum(!isFixed)) %*% Cuu)
            ## }
            rp0 <- obj0$report(p1)
            vMod <- c(R = mean(apply((rp0$mohnRhoVec_rec),1,function(x)(x[1]-x[2])/log(10))),
                      SSB = mean(apply((rp0$mohnRhoVec_ssb),1,function(x)(x[1]-x[2])/log(10))),
                      Fbar = mean(apply((rp0$mohnRhoVec_fbar),1,function(x)(x[1]-x[2])/log(10))))
            vOld <- c(R = mean(apply(exp(rp0$mohnRhoVec_rec),1,function(x)x[1]/x[2]-1)),
                      SSB = mean(apply(exp(rp0$mohnRhoVec_ssb),1,function(x)x[1]/x[2]-1)),
                      Fbar = mean(apply(exp(rp0$mohnRhoVec_fbar),1,function(x)x[1]/x[2]-1)))
            list(Modified = vMod, Original = vOld)
        }
        estRho <- doOne0(FALSE)
        estRho_Orig <- estRho$Original
        estRho_Mod <- estRho$Modified
        simRho <- replicate(simDelta,tryCatch(doOne0(TRUE), error = function(e) list(Modified = c(Fbar=NA,SSB=NA,R=NA), Original = c(Fbar=NA,SSB=NA,R=NA))), simplify = FALSE)
        simRho_Orig <- do.call("cbind",lapply(simRho, function(x) x$Original))
        simRho_Mod <- do.call("cbind",lapply(simRho, function(x) x$Modified))
        bginfo_Orig <- list(est = estRho_Orig, sim = simRho_Orig)
        bginfo_Mod <- list(est = estRho_Mod, sim = simRho_Mod)
        if(quantile_CI){
            cil_Orig <- apply(simRho_Orig,1,quantile,prob=0.025, na.rm=TRUE)
            cih_Orig <- apply(simRho_Orig,1,quantile,prob=0.975, na.rm=TRUE)
            cil_Mod <- apply(simRho_Mod,1,quantile,prob=0.025, na.rm=TRUE)
            cih_Mod <- apply(simRho_Mod,1,quantile,prob=0.975, na.rm=TRUE)
        }else{
            sdv_Orig <- apply(simRho_Orig,1,sd,na.rm=TRUE)
            cil_Orig <- estRho_Orig - 2 * sdv_Orig
            cih_Orig <- estRho_Orig + 2 * sdv_Orig
            sdv_Mod <- apply(simRho_Mod,1,sd,na.rm=TRUE)
            cil_Mod <- estRho_Mod - 2 * sdv_Mod
            cih_Mod <- estRho_Mod + 2 * sdv_Mod
        }
        tab <- list(Original = cbind(Estimate = estRho_Orig,
                                     CI_low = cil_Orig,
                                     CI_high = cih_Orig),
                    Modified = cbind(Estimate = estRho_Mod,
                                     CI_low = cil_Mod,
                                     CI_high = cih_Mod))
        bginfo <- list(Original = bginfo_Orig,
                       Modified = bginfo_Mod)
    }else{
        if(!resampleRE){
            if(addCorrelation){
                Sig_uu_tmp <- retro_hessian_RE(retroMS, returnSigma = TRUE, keep.diagonal = TRUE, forcePosDef = TRUE)
                Sig_uu <- as(Sig_uu_tmp,"sparseMatrix") ## Matrix::symmpart(Sig_uu_tmp)
                Sig_Chol_uu <- Matrix::Cholesky(Sig_uu)
                Hes_uu <- Matrix::symmpart(Matrix::solve(Sig_Chol_uu))
            }else{
                obj0 <- attr(retroMS,"m_obj")
                Hes_uu <- Matrix::symmpart(obj0$env$spHess(obj0$env$last.par.best, random = TRUE))
                Sig_uu <- as(Matrix::symmpart(Matrix::solve(Matrix::Cholesky(Hes_uu))),"sparseMatrix")
                ##Sig_Chol_uu <- Matrix::Cholesky(Sig_uu)
            }
        }## Otherwise calculated above
        
        ## Manual version of sdreport
        ## Modified from TMB
        obj <- attr(retroMS,"m_obj")
        obj2 <- TMB::MakeADFun(obj$env$data,
                          obj$env$parameters,
                          type = "ADFun",
                          ADreport = TRUE,
                          DLL = obj$env$DLL,
                          silent = obj$env$silent)
        r <- obj$env$random
        par <- obj$env$last.par.best
        ## Vtheta <- Sig0
        ## hessian.random <- Hes_uu
        phi <- obj2$fn(par)
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
        Dphi.random <- as(Dphi[,r,drop=FALSE],"sparseMatrix")
        Dphi.fixed <- as(Dphi[,-r,drop=FALSE],"sparseMatrix")        
        ##tmp <- Matrix::solve(Hes_uu,Matrix::t(Dphi.random))
        tmp <- Sig_uu %*% Matrix::t(Dphi.random)
        ##tmp <- as.matrix(tmp)
        term1 <- Dphi.random%*%tmp ## first term.

        w <- rep(0, length(par))
        if(!ADGradForward0Initialized) ADGradForward0Initialize()
        reverse.sweep <- function(i){
            w[r] <- tmp[,i]
            -obj$env$f(par, order = 1, type = "ADGrad", rangeweight = w, doforward=0)[-r]
        }
        A <- t(do.call("cbind",lapply(seq_along(phi), reverse.sweep))) + Dphi.fixed
        term2 <- A %*% (Sig0 %*% Matrix::t(A)) ## second term
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
        ## sdr <- TMB::sdreport(attr(retroMS,"m_obj"), attr(retroMS,"m_opt")$par, Hes, ...)
        ## ssdr <- summary(sdr)
        ## ssb <- ssdr[grepl("mohnRho_ssb",rownames(ssdr)),] %*% cbind(Estimate = c(1,0), CI_low = c(1,-2), CI_high = c(1,2))
        ## fbar <- ssdr[grepl("mohnRho_fbar",rownames(ssdr)),] %*% cbind(Estimate = c(1,0), CI_low = c(1,-2), CI_high = c(1,2))
        ## rec <- ssdr[grepl("mohnRho_rec",rownames(ssdr)),] %*% cbind(Estimate = c(1,0), CI_low = c(1,-2), CI_high = c(1,2))
        ## tab <- rbind(rec,ssb, fbar)
        ## bginfo <- ssdr
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
