collect_data <- function(x) {
    if(class(x) != "samset")
        stop("x must be a samset.")

    dat <- lapply(x,function(y){
        y$obj$env$data
    })
    names(dat) <- paste("Stock",1:length(dat))
    return(dat)
}
