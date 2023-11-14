
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


retro_hessian <- function(mFit, oFit){
##### Calculate hessian with correlation #####
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
        optimHess(m_opt$par[parYear==y],oFit$obj$fn,oFit$obj$gr)
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
    Hy <- optimHess(Obj1$par, Obj1$fn, Obj1$gr)
### The data appears multiple times for some years, use average
    diagA <- max(years) - fake_year[!is.na(map$fake_obs)] + 1
    A <- diag(sqrt(1/diagA),length(diagA))
    ## Approximate variance of data
    Vy <- stockassessment:::svd_solve(A %*%Hy[isObs,isObs]%*%A)
    ## J1_2 is ordered by parameter then year
    J1_2 <- Hy[max(which(isObs)) + which(!isFirstYear),isObs] # Derivative wrt new parameter then y
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
    ## Convert to Hessian for sdreport
    Hes1 <- try({solve(Sig1)},silent=TRUE)
    i <- 18
    while(is(Hes1,"try-error")){
        i <- i - 1
        Hes1 <- try({solve(Sig1 + diag(10^(-i),ncol(Sig1)))},silent=TRUE)
    }
    Hes1
}

##' @export
mohn_CI <- function(fit, ...){
    UseMethod("mohn_CI")
}

##' @export
mohn_CI.samset <- function(fit, addCorrelation = TRUE, simDelta = 0, ...){
    ## Already fitted with retro
    if(is.null(attr(fit,"fit")))
        stop("The samset should have a fit as attribute")
    fitList <- do.call("c",c(list(attr(fit,"fit")), fit))
    retroMS <- multisam.fit(fitList, mohn=1, doSdreport = FALSE)

    if(addCorrelation){
        Hes <- retro_hessian(retroMS, fitList[[length(fitList)]])
    }else{
        Hes <- attr(retroMS,"m_opt")$he        
    }

    if(simDelta > 0){
        opt0 <- attr(retroMS,"m_opt")
        Sig0 <- solve(he0)
        obj0 <- attr(retroMS,"m_obj")
        
        doOne0 <- function(sim=TRUE){
            if(sim){
                p0 <- stockassessment:::rmvnorm(1, opt0$par, Sig0, pivot = TRUE)
                p0[names(opt0$par)=="transfIRARdist"] <- pmin(pmax(p0[names(opt0$par)=="transfIRARdist"],-2),2)
            }else{
                p0 <- opt0$par
            }
            o <- capture.output(obj0$fn(p0))
            rp0 <- obj0$report()
            c(Fbar = mean(apply(exp(rp0$mohnRhoVec_fbar),1,function(x)x[1]/x[2]-1)),
              SSB = mean(apply(exp(rp0$mohnRhoVec_ssb),1,function(x)x[1]/x[2]-1)),
              R = mean(apply(exp(rp0$mohnRhoVec_rec),1,function(x)x[1]/x[2]-1)))
        }
        estRho <- doOne0(FALSE)
        simRho <- replicate(simDelta,doOne0())
        tab <- cbind(Est = estRho,
                     CI_low = apply(simRho,1,quantile,prob=0.025),
                     CI_high = apply(simRho,1,quantile,prob=0.975))
    }else{
        sdr <- sdreport(attr(retroMS,"m_obj"), attr(retroMS,"m_opt")$par, Hes, ...)
        ssdr <- summary(sdr)
        ssb <- ssdr[grepl("mohnRho_ssb",rownames(ssdr)),] %*% cbind(Est = c(1,0), CI_low = c(1,-2), CI_high = c(1,2))
        fbar <- ssdr[grepl("mohnRho_fbar",rownames(ssdr)),] %*% cbind(Est = c(1,0), CI_low = c(1,-2), CI_high = c(1,2))
        rec <- ssdr[grepl("mohnRho_rec",rownames(ssdr)),] %*% cbind(Est = c(1,0), CI_low = c(1,-2), CI_high = c(1,2))
        tab <- rbind(fbar,ssb,rec)
    }

    res <- list(table = tab,
                #tableMod = ssdr[grepl("mohnRhoMod_",rownames(ssdr)),],
                addCorrelation = addCorrelation,
                simDelta = simDelta,
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
mohn_sim_CI.sam <- function(fit, nsim, type = c("Full","Gauss","GaussF","Tail","FixF"), ncores = 1, year = 5){
    RETRO <- retro(fit, year = year, ncores = ncores)
    mohn_sim_CI(RETRO, nsim = nsim, type = type, ncores = ncores)
    
}

##' @export
mohn_sim_CI.samset <- function(fit, nsim, type = c("Full","Gauss","GaussF","Tail","FixF"), ncores = 1){

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
    def <- list(Normal = def, Modified = def)
    if(type == "Tail"){
        doOne <- function(){
            simfit <- try({stockassessment:::addSimulatedYears(RETRO[[length(RETRO)]],rep(NA,length(RETRO)), resampleFirst = TRUE, refit = TRUE)})
            if(is(simfit,"try-error"))
                return(def)
            re_retro <- try({stockassessment::retro(simfit, length(RETRO), ncores=ncores)},silent=TRUE)
            if(is(re_retro,"try-error"))
                return(def)            
            list(Normal = mohn_2(re_retro), Modified = mohn_2(re_retro))
        }
    }else{
        fit <- attr(RETRO,"fit")
        obj <- unserialize(serialize(fit$obj, NULL))
        ## Fix in simulation:
#### F (all but Full and Tail),
#### N (Fix for Gauss),
#### Obs (None)
        obj$env$data$simFlag <- c(ifelse(type %in% c("Full"),0,1),
                                  ifelse(type %in% c("Gauss"),1,0),0)
        obj$retape()
        obj$fn()
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
            if(type %in% c("GaussF","Gauss")){
                p0 <- obj$env$last.par.best
                p0[obj$env$random] <- stockassessment:::rmvnorm(1,p0[obj$env$random],C0, pivot = TRUE)
                pl <- obj$env$parList(par = p0)
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
            list(Normal = mohn_2(re_retro), Modified = mohn_2(re_retro))
        }        
    }
    replicate(nsim, doOne())
}
