
as.FLStock.msam <- function(mfit, unit.w = "kg", name = "", desc = "", predicted = FALSE){
    source(system.file("flr_convenience.R",package = "stockassessment"))
    toFLQ <- .SAM_FLR_HELP_toFLQ
    resize <- .SAM_FLR_HELP_resize
    na2zero <- .SAM_FLR_HELP_na2zero

    fls <- lapply(seq_along(mfit), function(i) as.FLStock.sam(mfit[[i]], name = names(mfit)[i]))
    ## Update estimated things
    for(i in seq_along(fls)){
        fit <- mfit[[i]]
        ages <- fit$conf$minAge:fit$conf$maxAge
        years <- fit$data$years
        ## Need to update catch/landings/discards if predicted = TRUE
        if(predicted){
            ## Modified from stockassessment:::getFleet
            grabFleet <- function(fit, fleet){
                fidx <- fit$data$aux[, "fleet"] == fleet
                aux <- fit$data$aux[fidx, ]
                logout <- attr(mfit,"m_rep")$predObs[[i]]
                .goget <- function(y, a) {
                    ret <- exp(logout[aux[, "year"] == y & aux[, "age"] == a])
                    ifelse(length(ret) == 0, 0, ret)
                }
                yr <- min(aux[, "year"]):max(aux[, "year"])
                ar <- min(aux[, "age"]):max(aux[, "age"])
                tmp <- outer(yr, ar, Vectorize(.goget))
                dimnames(tmp)[[1]] <- yr
                dimnames(tmp)[[2]] <- ar
                attr(tmp, "time") <- c(fit$data$sampleTimesStart[fleet], 
                                       fit$data$sampleTimesEnd[fleet])
                return(tmp)
            }
            CN_all <- simplify2array(lapply(which(fit$data$fleetTypes %in% c(0,1,7)),function(i)grabFleet(fit,i)))
            CN <- toFLQ(resize(na2zero(CN_all),ages,years))
            DN <- CN * (1 - LF)
            LN <- CN * LF
            catch.n(fls[[i]]) <- CN
            discard.n(fls[[i]]) <- DN
            landing.n(fls[[i]]) <- LN
        }
        ## Need to update F
        fixColnames <- function(x){
            colnames(x) <- gsub("Age ","",colnames(x))
            x
        }
        harvest(fls[[i]]) <- toFLQ(resize(fixColnames(faytable(mfit,returnList=TRUE)[[i]]),ages,years,FALSE), "f")
        ## Need to update SN
        stock.n(fls[[i]]) <- toFLQ(resize(fixColnames(ntable(mfit,returnList=TRUE)[[i]]),ages,years))

    }
  fls
}
