getLowerBounds <- function(parameters){
    f <- get("getLowerBounds", envir = asNamespace("stockassessment"), inherits = FALSE)
    r <- f(parameters)
    r$RE <- rep(-10,length(parameters$RE))
    return(r)
}

getUpperBounds <- function(parameters){
    f <- get("getUpperBounds", envir = asNamespace("stockassessment"), inherits = FALSE)
    r <- f(parameters)
    r$RE <- rep(10,length(parameters$RE))
    return(r)
}
