getColMeanCovariate <- function(d, ages, agesAll){
    v <- numeric(length(agesAll))
    if(is.null(d))
        return(v)
    v[match(ages, agesAll)] <- as.vector(colMeans(d))
    v
}
getAgeRange <- function(data){
    min(data$minAgePerFleet):max(data$maxAgePerFleet)
}

safe_values <- function(mf){
    mfl <- as.list(mf)
    r <- as.data.frame(lapply(mfl, function(xx){
        i <- is.na(xx) | is.nan(xx) | is.null(xx) | !is.finite(xx)
        if(any(i))
            warning("Non-finite values replaced by zero")
        xx[i] <- 0
        xx
    }))
    attributes(r) <- attributes(mf)
    r
}

build_cov_design <- function(formula, x){
    
    if(!inherits(x, "samset"))
        return("x must be a samset")

    sn <- getStockNames(x)
    ageAll <- min(unlist(lapply(x,function(y)y$data$minAgePerFleet))) : max(unlist(lapply(x,function(y)y$data$maxAgePerFleet)))
    potentialCovariateNames <- unique(unlist(lapply(x, function(y){
        i <- sapply(y$data, function(d) is.matrix(d) && ncol(d) == length(getAgeRange(y$data)))
        names(y$data)[i]
    })))
   
    potentialCovariates <- lapply(x, function(y){
        datReduced <- y$data[potentialCovariateNames]
        names(datReduced) <- potentialCovariateNames
        v <- as.data.frame(lapply(datReduced, function(d) getColMeanCovariate(d, getAgeRange(y$data), ageAll)))
        })

    XX <- cbind(expand.grid(Age = ageAll, Stock = sn), do.call("rbind",potentialCovariates))
    oFun <- function(x,y){
        if(is.numeric(x))
            return(x-y)
        return(paste(x,y,sep=":"))
    }
    designDat <- as.data.frame(lapply(XX, function(x){
        m <- outer(x,x,oFun)
        m[lower.tri(m, diag = FALSE)]
    }))
    designDat$Index <- seq_len(nrow(designDat))
    designDat$WithinStock <- unlist(lapply(strsplit(as.character(designDat$Stock),":"),function(x)as.numeric(x[1]==x[2])))
    mf <- safe_values(model.frame(formula, data = designDat, na.action = NULL))
    tt <- terms(mf)
    mm <- model.matrix(tt,mf)
    attr(mm, "terms") <- tt
    attr(mm, "data") <- get_all_vars(formula,designDat)
    return(mm)
}
