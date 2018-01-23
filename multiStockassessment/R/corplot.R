##' Correlation and partial correlation plot between survival processes
##'
##' 
##' @title Correlation and partial correlation plot between survival processes
##' @param object msam object
##' @param ... Other parameters not currently used
##' @return A matrix of correlations (lower triangular) and partial correlations (upper triangular)
##' @author Christoffer Moesgaard Albertsen
##' @importFrom graphics image axis abline box mtext text rect
##' @importFrom stats cov2cor
##' @importFrom grDevices colorRampPalette
##' @importFrom stockassessment corplot
##' @method corplot msam
##' @export
corplot.msam <- function(object,...){
    rp <- attr(object,"m_rep")
    dat <- attr(object,"m_dat")
    nc <- rp$Sigma
    m <- -solve(rp$Sigma)
    diag(m) <- -diag(m)
    npc <- stats::cov2cor(m)

    stockNames <- names(dat$sam)
    
    nAreas <- length(object)
    
    nc[lower.tri(nc)] <- npc[lower.tri(npc)]
    ala <- dat$minAgeAll:dat$maxAgeAll
    colnames(nc) <- rownames(nc) <- rep(ala,nAreas) #paste0(rep(LETTERS[1:nAreas],each=length(ala)),ala)

    ##ncCI <- getCorCI(indx)
    ageList <- lapply(dat$sam,function(x)x$minAge:x$maxAge)
    brks <- seq(-1,1,0.1)
    brks[1] <- -1.1; brks[21] <- 1.1
    graphics::image(x=1:dim(nc)[1],y=1:dim(nc)[2],nc,
          axes=FALSE,
          ylim=c(dim(nc)[1]+0.5,1-0.5),
          xlab="",ylab="",
          breaks=brks,
          col=addTrans(grDevices::colorRampPalette(c("blue","white","red"))(20),0.5))
    graphics::axis(1,labels=colnames(nc),at=1:dim(nc)[1])
    graphics::axis(2,labels=colnames(nc),at=1:dim(nc)[1])
    for(i in 1:dim(nc)[1])
        for(j in 1:dim(nc)[1])
            if(i != j & abs(nc[i,j]) > 1e-10){
                txt <- paste0(formatC(nc[i,j],format="f",digits=2),
                              ##"\n",
                              ##"(",formatC(ncCI$CIL[i,j],format="f",digits=2),"; ",
                              ##formatC(ncCI$CIU[i,j],format="f",digits=2),")",
                              ""
                              )
                                        #text(i,j,formatC(nc[i,j],format="f",digits=2))
                ## if(!(ncCI$CIL[i,j] < 0 & ncCI$CIU[i,j] > 0))
                ##     polygon(c(i+0.5,i+0.5-0.2,i+0.5,i+0.5),
                ##             c(j-0.5,j-0.5,j-0.5+0.2,j-0.5),
                ##             col="black")
                graphics::text(i,j,txt)
            }
    for(i in 1:dim(nc)[1]){
        colIndx <- as.numeric(rownames(nc)[i] %in% unlist(sapply(1:nAreas,function(aa)ageList[[aa]]))) + 1
        graphics::rect(i-0.5,i-0.5,i+0.5,i+0.5,col=c("white","black")[colIndx],border=NA)
    }
    graphics::abline(v=1:(dim(nc)[1]+1)-0.5,h=1:(dim(nc)[1]+1)-0.5,col="darkgrey")
    graphics::abline(v = (length(ala)) * (1:(nAreas-1)) + 0.5,lwd=2)
    graphics::abline(h = (length(ala)) * (1:(nAreas-1)) + 0.5,lwd=2)
    graphics::box(lwd=2)
    graphics::mtext(stockNames,side=1,at=length(ala)/2 + length(ala) * (1:nAreas -1) + 0.5, line = 3)
    graphics::mtext(stockNames,side=2,at=length(ala)/2 + length(ala) * (1:nAreas -1) + 0.5, line = 2.5)
                                        #box("figure")
                                        #mtext(paste0("(",LETTERS[which(indxUse==indx)],")"),outer=FALSE,adj=0,line=1,cex=1.5,font=2)
    colnames(nc) <- rownames(nc) <- paste(rep(stockNames,each=length(ala)),ala,sep=" age ")
    invisible(t(nc))
}

    
