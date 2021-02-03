getColMeanCovariate <- function(d, ages, agesAll){
    v <- numeric(length(agesAll))
    if(is.null(d))
        return(v)
    v[match(ages, agesAll)] <- as.vector(colMeans(d))
    v
}
getAgeRange <- function(conf){
    min(conf$minAge):max(conf$maxAge)
}

safe_values <- function(mf){
    mfl <- as.list(mf)
    r <- as.data.frame(lapply(mfl, function(xx){
        i <- is.na(xx) | is.nan(xx) | is.null(xx) | !is.finite(xx)
        if(!any(i))
            return(xx)
        if(is.factor(xx) | is.character(xx)){
            warning("Non-finite factor values replaced by 'Unknown'")
            xx <- as.character(xx)
            xx[i] <- "Unknown"
            return(factor(xx))
        }
        warning("Non-finite values replaced by zero")
        xx[i] <- 0
        
        xx
    }))
    attributes(r) <- attributes(mf)
    r
}

## J <- function(x){
##     cat(ls(env=parent.env(environment())),sep="\n")
##     x2 <- paste(as.character(match.call())[2],"J",sep="_")
##     eval(parse(text = x2))
## }
## K <- function(x){
##     x2 <- paste(as.character(match.call())[2],"K",sep="_")
##     eval(parse(text = x2))
## }


build_cov_design <- function(formula, x, type = NA){
    if(!inherits(x, "samset"))
        return("x must be a samset")

    sn <- getStockNames(x)
    ageAll <- min(unlist(lapply(x,function(y)y$conf$minAge))) : max(unlist(lapply(x,function(y)y$conf$maxAge)))
    potentialCovariateNames <- unique(unlist(lapply(x, function(y){
        i <- sapply(y$data, function(d) is.matrix(d) && ncol(d) == length(getAgeRange(y$conf)))
        names(y$data)[i]
    })))
   
    potentialCovariates <- lapply(x, function(y){
        datReduced <- y$data[potentialCovariateNames]
        names(datReduced) <- potentialCovariateNames
        v <- as.data.frame(lapply(datReduced, function(d) getColMeanCovariate(d, getAgeRange(y$conf), ageAll)))
        })

    XX <- cbind(expand.grid(Age = ageAll, Stock = sn), do.call("rbind",potentialCovariates))
    oFun <- function(x,y){
        if(is.numeric(x))
            return(x-y)
        return(paste(x,y,sep=":"))
    }
    oFunJ <- function(x,y){
        if(is.numeric(x))
            return(x)
        return(paste(x))
    }
    oFunK <- function(x,y){
        if(is.numeric(x))
            return(y)
        return(paste(y))
    }
    m0 <- outer(XX$Age, XX$Age, oFun)
    if(is.na(type)){                    # Keep lower tri
        indx <- sort(which(lower.tri(m0, diag = FALSE)))
    }else if(type == 1 || type == 2){   # Keep all but diagonal
        indx <- sort(which(lower.tri(m0, diag = FALSE) | upper.tri(m0, diag = FALSE)))
    }else{                              # Keep all
        indx <- 1:length(m0)
    }              
    designDatC <- as.data.frame(lapply(XX, function(x){
        m <- outer(x,x,oFun)
        m[indx]
    }))
    designDatJ <- as.data.frame(lapply(XX, function(x){
        m <- outer(x,x,oFunJ)
        m[indx]
    }))
    colnames(designDatJ) <- paste0(colnames(designDatJ),"_J")
    designDatK <- as.data.frame(lapply(XX, function(x){
        m <- outer(x,x,oFunK)
        m[indx]
    }))
    colnames(designDatK) <- paste0(colnames(designDatK),"_K")
    
    designDat <- cbind(designDatC, designDatJ, designDatK)
    designDat$Index <- seq_len(nrow(designDat))
    designDat$WithinStock <- unlist(lapply(strsplit(as.character(designDat$Stock),":"),function(x)as.numeric(x[1]==x[2])))

    
    ## formula <- as.formula(gsub("^(J|K)(\\()(.+)(\\))",
    ##                            "\\3_\\1",
    ##                            as.character(formula)))

    mf <- safe_values(model.frame(formula,designDat, na.action = NULL))
    tt <- terms(mf)
    mm <- model.matrix(tt,mf)
    attr(mm, "terms") <- tt
    attr(mm, "data") <- get_all_vars(formula,designDat)
    return(mm)
}
