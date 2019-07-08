#################################################################################################
#################################################################################################
#### Plot methods modified from the stockassessment package to fit with multiStockassessment ####
#################################################################################################
#################################################################################################


.plotcols <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
.plotcols.crp <- grDevices::colorRampPalette(.plotcols)

plotit <- function(fit, what, ...){
    UseMethod("plotit")
}

##' Function to actually do the plotting
##' @rdname plotit
##' @method plotit msam
##' @importFrom grDevices colorRampPalette gray
##' @importFrom graphics  grid legend lines matplot plot
##' @param fit msam fit
##' @param what variable to plot
##' @param x x-axis values 
##' @param ylab y-axis label
##' @param xlab x-axis label
##' @param ex extra years to add
##' @param trans function to transform variable to plot
##' @param add If false a new plot is created. If true everything is added to the previous plot-
##' @param ci Add confidence intervals?
##' @param cicol Color of confidence intervals
##' @param addCI Not used
##' @param drop Number of years to drop from the end
##' @param unnamed.basename not used
##' @param xlim x-axis limits. If null, the range of x is used.
##' @param ciAlpha Alpha channel value of confidence interval color
##' @param col line colors
##' @param extraLabel Not used
##' @param ... Other arguments
##' @export
plotit.msam <- function(fit, what,
                        x=lapply(attr(fit,"m_data")$sam,function(x)x$years),
                        ylab=what,
                        xlab="Years",
                        ex=numeric(0),
                        trans=function(x)x,
                        add=FALSE,
                        ci=TRUE,
                        cicol=gray(.5,alpha=.5),
                        addCI=NA,
                        drop=0,
                        unnamed.basename="current",
                        xlim=NULL,
                        ciAlpha = 0.3,
                        col = .plotcols.crp(length(fit)),
                        extraLabel=NULL,
                        ...){
    tabList <- tableit(fit=fit,what=what,x=x,trans=trans, returnList = TRUE)
    tab <- cbindYearTables(tabList)

    xTab <- as.numeric(rownames(tab))
    xAll <- min(unlist(x)):max(unlist(x))
    xindx <- which(xAll%in%xTab)
    
    y <- low <- high <- matrix(NA,length(xAll),length(fit))
    y[xindx,] <- tab[,(0:(length(tabList)-1))*3 + 1]
    low[xindx,] <- tab[,(0:(length(tabList)-1))*3 + 2]
    high[xindx,] <- tab[,(0:(length(tabList)-1))*3 + 3]
    
    didx <- 1:(nrow(tab)-drop)
 
    
    if(missing(xlim)){
        xr <- range(unlist(x))
    }else{
      xr <- xlim
    }
    x<-x[didx]
    y<-y[didx,]
    low<-low[didx,]
    high<-high[didx,]
    if(add){
        for(i in 1:length(fit))
            lines(xAll, (y[,i]), lwd=3, col=col[i],...)
    }else{
        plot(xAll, NA*xAll, xlab=xlab, ylab=ylab, type="n", lwd=3, xlim=xr,
             ylim=range(c((low),(high),0,ex),na.rm=TRUE), las=1,...)
        grid(col="black")
        for(i in 1:length(fit))
            lines(xAll, (y[,i]), lwd=3, col=col[i], ...)
    }
    if(ci){
        for(i in 1:length(fit)){
            xuse <- xAll%in%x[[i]]
            indx <- which(!is.na(low[,i]) & !is.na(high[,i]) & xuse)
            polygon(c(xAll[indx],rev(xAll[indx])), y = c((low[indx,i]),rev((high[indx,i]))), col = addTrans(col[i],ciAlpha), border=NA)
            lines(xAll, (y[,i]), lwd=3, col=col[i])
        }
    }
    stocknames <- as.expression(parse(text=gsub("[[:space:]]","~",getStockNames(fit))))
    legend("bottom",legend=stocknames, lwd=3, col=col, ncol=min(length(fit),5), bty="n", merge=TRUE)
}


##' @rdname plotit
##' @method plotit msamforecast
##' @export
plotit.msamforecast <- function(fit,
                                what,
                                x=fit$data$years,
                                ylab=what,
                                xlab="Years",
                                ex=numeric(0),
                                trans=exp,
                                add=FALSE,
                                ci=TRUE,
                                cicol=gray(.5,alpha=.5),
                                addCI=NA,
                                drop=0,
                                unnamed.basename="current",
                                xlim=NULL,
                                ...){
    xy <- unlist(lapply(fit, function(xx) unlist(lapply(xx,function(yy)yy$year))))
    thisfit <- attr(fit,"fit")
    xr <- range(attr(thisfit,"m_data")$maxYearAll,attr(thisfit,"m_data")$minYearAll, xy)
    what2 <- switch(what,
                    logfbar = "fbar",
                    logR = "rec",
                    logLagR = "rec",
                    logssb = "ssb",
                    logCatch = "catch",
                    logtsb = "tsb",
                    what)
    
    if(length(ex) == 0){
        if(ci){
            ex <- sapply(1:length(fit),function(i){
                x <- attr(fit[[i]],"tab")
                range(x[,paste(what2,"low", sep=":")],x[,paste(what2,"high", sep=":")])
            })
        }else{
            ex <- sapply(1:length(fit),function(i){
                x <- attr(fit[[i]],"tab")
                range(x[,paste(what2,"median", sep=":")])
            })         
        }
    }
    plotit(thisfit, what=what, ylab=ylab, xlab=xlab, ex=ex, trans=trans, add=add, ci=ci, cicol=cicol, drop=drop, xlim=xr,...)
}


##' @rdname addforecast
addforecast<-function(fit, what, dotcol="black", dotpch=19, dotcex=1.5, intervalcol=gray(.5,alpha=.5),...){
    UseMethod("addforecast")
}

##' @rdname addforecast
##' @method addforecast msamforecast
##' @export
addforecast.msamforecast <- function(fit, what, dotcol=.plotcols.crp(length(fit)), dotpch=19, dotcex=1.5, intervalcol=addTrans(dotcol,0.5),...){
    xlist <- lapply(fit,attr, which = "tab")
    ylist <- lapply(xlist,function(x) as.numeric(rownames(x)))
    dummy <- sapply(1:length(xlist),function(s){
        y <- ylist[[s]]
        sapply(1:length(ylist[[s]]),
               function(i){
                   x <- xlist[[s]]
                   xx<-c(x[i,paste(what,"low", sep=":")],x[i,paste(what,"high", sep=":")]);
                   units = par(c('usr', 'pin'))
                   xx_to_inches = with(units, pin[2L]/diff(usr[3:4]))
                   if(abs(xx_to_inches * diff(xx))>0.01){
                       arrows(y[i],xx[1],y[i],xx[2],lwd=3, col=intervalcol[s], angle=90, code=3, length=.1)
                   }
               }
               )
    })
    for(s in 1:length(ylist)){
        y <- ylist[[s]]
        x <- xlist[[s]]
        points(y,x[,paste(what,"median", sep=":")], pch=dotpch, cex=dotcex, col=dotcol[s])
    }
}

##' Fbar plot for msam object
##'
##' @param fit fitted msam object
##' @param partial Should F for each age in Fbar range be added?
##' @param drop Number of years to drop from the end
##' @param page List of ages to be used per stock for partial = true. Defaults to all ages used to calculate Fbar.
##' @param ... plotting arguments
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment fbarplot
##' @export
fbarplot.msam <- function(fit,partial = FALSE, drop=0, page=NULL, plot = TRUE, ...){
    args <- list(...)
    d <- attr(fit,"m_data")$sam
    fbarRange <- lapply(d,function(x){
        r <- x$fbarRange
        attr(r,"label") <- substitute(bar(F)[X-Y],list(X=r[1],Y=r[2]))
        attr(r,"page") <- r[1]:r[2]
        r
    })
    flist<-faytable(fit,returnList = TRUE)
    if(is.null(page)){
        page<-lapply(fbarRange,function(x)attr(x,"page"))
    }
    fmatsInPage <- sapply(1:length(fit),function(i)flist[[i]][,d[[i]]$minAge:d[[i]]$maxAge %in% page[[i]],drop=FALSE],simplify = FALSE)
    exx <- if(partial){unlist(fmatsInPage)}else{numeric(0)}
    fit2 <- fit
    names(fit2) <- paste0(gsub("[[:space:]]","~",getStockNames(fit)),":~",
                                      unlist(lapply(fbarRange,function(x)attr(x,"label"))))
    if(plot){
        plotit(fit2, "logfbar", ylab=expression(bar(F)), trans=exp, ex=c(ex,exx), drop=drop, add=add, ...)
        if(partial){
            if(is.null(args$col)){
                col <- .plotcols.crp(length(fit))
            }else{
                col <- args$col
            }
            for(i in 1:length(fit)){
                idxx <- 1:(length(d[[i]]$years)-drop)
                matplot(d[[i]]$years[idxx],
                        fmatsInPage[[i]][idxx,], add=TRUE, type="b", col=col[i], pch=as.character(page[[i]]))
            }
        }
    }
    invisible(list(fbarRange = fbarRange, fbarlab = expression(bar(F)),
              exx = exx, drop = drop, page = page,
              fmatsInPage = fmatsInPage, d = d, col = col,
              stocknames = names(fit2)))
}

##' @export
fbarplot.msamforecast <- function(fit,partial = FALSE, drop=0, page=NULL,...){
    fitlocal <- attr(fit,"fit")
    tmp <- fbarplot(fitlocal,partial=partial,drop=drop,page=page,plot=FALSE,...)
    fit2 <- fit
    names(attr(fit2,"fit")) <- tmp$stocknames
    cat(names(fit2))
    plotit(fit2, "logfbar", ylab=tmp$fbarlab, trans=exp, ex=tmp$exx, drop=tmp$drop, ...)
    if(partial){
        for(i in 1:length(fit)){
            idxx <- 1:(length(tmp$d[[i]]$years)-tmp$drop)
            matplot(tmp$d[[i]]$years[idxx],
                    tmp$fmatsInPage[[i]][idxx,], add=TRUE, type="b", col=tmp$col[i], pch=as.character(tmp$page[[i]]))
        }
    }
    addforecast(fit2,"fbar")
}

##' SSB plot
##'
##' @param fit msam object
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment ssbplot
##' @export
ssbplot.msam <- function(fit, ...){
    plotit(fit, "logssb", ylab="SSB", trans = exp, ...)
}

##' @export
ssbplot.msamforecast <- function(fit, ...){
    plotit(fit, "logssb", ylab="SSB", trans = exp, ...)
    addforecast(fit,"ssb")
}

##' TSB plot
##'
##' @param fit msam object
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment tsbplot
##' @export
tsbplot.msam <- function(fit,...){
    plotit(fit, "logtsb", ylab="TSB", trans=exp,...)
}

##' @export
tsbplot.msamforecast <- function(fit, ...){
    if(!all(sapply(fit,function(x) any(colnames(attr(x,"tab")) == "tsb:median"))))
        stop("The forecast was not made with addTSB = TRUE")
    plotit(fit, "logtsb", ylab="TSB", trans = exp, ...)
    addforecast(fit,"tsb")
}


##' Recruitment plot
##'
##' @param fit msam object
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment recplot
##' @export
recplot.msam <- function(fit,...){
    fit2 <- fit
    names(fit2) <- paste0(gsub("[[:space:]]","~",getStockNames(fit)),":~Age~",
                          unlist(lapply(attr(fit,"m_data")$sam,function(x)x$minAge)))
    plotit(fit2, "logR", ylab="Recruits", trans=exp,...)
}

##' @export
recplot.msamforecast <- function(fit, lagR = FALSE, ...){
    fit2 <- fit    
    fitlocal <- attr(fit,"fit")
    if(lagR){
        names(attr(fit2,"fit")) <- paste0(gsub("[[:space:]]","~",getStockNames(fitlocal)),":~Age~",
                                  unlist(lapply(attr(fitlocal,"m_data")$sam,function(x)x$minAge+1)))
        what <- "logLagR"
    }else{
        names(attr(fit2,"fit")) <- paste0(gsub("[[:space:]]","~",getStockNames(fitlocal)),":~Age~",
                                  unlist(lapply(attr(fitlocal,"m_data")$sam,function(x)x$minAge)))
        what <- "logR" 
    }
        plotit(fit2, what, trans=exp,...)
        addforecast(fit2, "rec")
}




##' Catch plot
##'
##' @param fit msam object
##' @param obs.show Show observations in plot?
##' @param drop Number of years to drop from the end
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment catchplot
##' @importFrom graphics points polygon
##' @export
catchplot.msam <- function(fit, obs.show=TRUE, drop=0, ...){
    args <- list(...)
    ctab <- catchtable(fit, returnList = TRUE,obs.show=TRUE)
    plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=drop,...)
    if(obs.show){
        if(is.null(args$col)){
            col <- .plotcols.crp(length(fit))
        }else{
            col <- args$col
        }
        for(i in 1:length(fit))
            points(as.numeric(rownames(ctab[[i]])),ctab[[i]][,4],
                   pch=4, lwd=2, cex=1.2, col=col[i])
    }
    obs <- list(x=x,y=rowSums(outer(rownames(CW), colnames(CW), Vectorize(.goget))*CW, na.rm=TRUE))
    invisible(list(drop=drop,obs=obs))
}


##' @export
catchplot.msamforecast <- function(fit, obs.show=TRUE, drop=0, ...){
    args <- list(...)
    ctab <- catchtable(attr(fit,"fit"), returnList = TRUE,obs.show=TRUE)
    plotit(fit, "logCatch", ylab="Catch", trans=exp, drop=drop,...)
    if(obs.show){
        if(is.null(args$col)){
            col <- .plotcols.crp(length(fit))
        }else{
            col <- args$col
        }
        for(i in 1:length(fit))
            points(as.numeric(rownames(ctab[[i]])),ctab[[i]][,4],
                   pch=4, lwd=2, cex=1.2, col=col[i])
    }
   addforecast(fit,"catch")
}


##' Parameter plot
##'
##' @param fit msam object
##' @param cor.report.limit Not used
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment parplot
##' @export
parplot.msam <- function(fit, cor.report.limit=0.95, ...){
    stop("Not implemented")
}


##' Stock-recruitment plot
##'
##' @param fit msam object
##' @param textcol Text color
##' @param add Should the figure be added to a current plot?
##' @param col Line colors
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment srplot
##' @export
srplot.msam <- function(fit,textcol="red",add=FALSE,
                        col = .plotcols.crp(length(fit)), ...){
    ns <- length(fit)
    X <- summary(fit, returnList = TRUE)
    n <- lapply(X,nrow)
    lag <- lapply(attr(fit,"m_data")$sam,function(x)x$minAge)
    R <- sapply(1:ns,function(i)X[[i]][(lag[[i]]+1):n[[i]],1],simplify=FALSE)
    SSB <- sapply(1:ns,function(i)X[[i]][(lag[[i]]+1):n[[i]],4],simplify=FALSE)
    y <- lapply(X,rownames)
    plot(NA,NA,xlab="SSB",ylab="R",type="n",
         xlim=range(0,max(unlist(SSB))),
         ylim=range(0,max(unlist(R))),...)    
    for(i in 1:length(fit)){
        lines(SSB[[i]],R[[i]], lwd=1, col=col[i], ...)
        text(SSB[[i]],R[[i]],labels=y[[i]][(lag[[i]]+1):n[[i]]],cex=.7,col=col[i])
    }
    stocknames <- as.expression(parse(text=gsub("[[:space:]]","~",getStockNames(fit))))
    legend("bottom",legend=stocknames, lwd=3, col=col, ncol=min(length(fit),5), bty="n", merge=TRUE)
}


##' Fit plot
##'
##' @param fit msam object
##' @param stock Stock to plot for
##' @param log Plot on log-scale?
##' @param fleets Fleets to plot
##' @param ... extra arguments for plotting
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment fitplot plotby
##' @export
fitplot.msam <- function(fit,
                         stock,
                         log=TRUE,
                         fleets=lapply(attr(fit,"m_data")$sam,function(x)unique(x$aux[,"fleet"])),
                         ...){
    d <- attr(fit,"m_data")$sam
    rep <- attr(fit,"m_rep")
    idx <- sapply(1:length(d),function(i)d[[i]]$aux[,"fleet"] %in% fleets[[i]],
                  simplify = FALSE)
    trans <- function(x)if(log){x}else{exp(x)}
    p <- sapply(1:length(d),function(i)trans(rep$predObs[[i]][idx[[i]]]),
                simplify=FALSE)
    o <- sapply(1:length(d),function(i)trans(d[[i]]$logobs[idx[[i]]]),
                simplify=FALSE)
    fl <-  sapply(1:length(d),function(i)d[[i]]$aux[idx[[i]],"fleet"],
                  simplify = FALSE)
    yy <-  sapply(1:length(d),function(i)d[[i]]$aux[idx[[i]],"year"],
                  simplify = FALSE)
    aa <- sapply(1:length(d),function(i){
        atmp <- d[[i]]$aux[idx[[i]],"age"]
        atmp[atmp < -1.0e-6] <- NA
        atmp
    }, simplify=FALSE)
    a <- paste0("a=",aa[[stock]]," ")
    f <- paste0(" f=",strtrim(attr(fit[[stock]]$data,"fleetNames")[fl[[stock]]],50))
    Year <- yy[[stock]]
    if(length(fleets)==1){
        myby <- paste(a, ":", f)
    }else{
        myby <- cbind(a,f)
    }
    stockassessment::plotby(Year, o[[stock]], y.line=p[[stock]], by=myby, y.common=FALSE, ylab="", ...)
}
