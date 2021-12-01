prepareGenetics <- function(baselineFile,
                            sampleFile,
                            dataFile){

    if(missing(baselineFile) && missing(sampleFile)){
        r <- list(samples = list(),
                  Qspace = as(matrix(0,0,0),"dgTMatrix"),
                  Qtime = as(matrix(0,0,0),"dgTMatrix"),
                  QorderTime = 1L,
                  QorderSpace = 1L,
                  stock2gen = matrix(0,0,0))
        attr(r, "alleleDim") <- c(0,0) # ploidity x nloci
        attr(r, "nTrips") <- 0
        class(r) <- "msam_genetic"
        return(r)
    }
    stop("Not implemented")
        
}
                



## Allele frequency heatmap

allelefrequencyplot <- function(x,...){
   
}





## Plot PCA grouping
allelePCAplot <- function(x, ...){

}


allelePCAplot.msam_genetic <- function(x, ...){
    gg <- simplify2array(lapply(x$samples, function(x) x$alleleMatrix))
    s2g <- x$stock2gen
    gsNames <- if(!is.null(colnames(s2g))){ colnames(s2g) }else{ paste("Genetic stock",seq_len(ncol(s2g))) }
    ggStock <- lapply(x$samples, x$keyStock+1)
    pxB <- t(gen2PCA(gg))
    plot(prcomp(pxB)$x[,1:2], pch=16,col=factor(ggNames[ggStock], ggNames))
    legend("bottomright",ggNames,
           col = seq_along(ggNames),
           pch = 16)    
}
