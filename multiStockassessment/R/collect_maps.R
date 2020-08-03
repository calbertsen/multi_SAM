collect_maps <- function(x) {
    if(class(x) != "samset")
        stop("x must be a samset.")

    
    ncase <- length(x)
    newpars <- multiStockassessment:::collect_pars(x)
    singleMaps <- lapply(x,function(xx)xx$obj$env$map)
    mapNames <- lapply(singleMaps,names)
    uniqueNames <- unique(unlist(mapNames))
    newmap <- lapply(newpars[uniqueNames], function(xx){
        a <- as.numeric(seq_along(xx))
        attributes(a) <- attributes(xx)
        a
    })

    for(i in 1:ncase)
        for(nm in mapNames[[i]]){
            if(length(newmap[[nm]]) > 0){
                tmp <- splitParameter(newmap[[nm]])
                tmp[[i]] <- tmp[[i]][singleMaps[[i]][[nm]]]
                newmap[[nm]] <- combineParameter(tmp)
            }
        }
    newmap <- lapply(newmap, factor)

    return(newmap)
}
