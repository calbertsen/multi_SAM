setwd(sub('script.R','',"tests/tables_onestock/script.R"))
sink('/dev/null',type='output')
capture.output({
## Tests the combination of ALN likelihood with AR correlation structure
library(stockassessment)
data(nscodData)
data(nscodConf)

nscodConf$obsLikelihoodFlag[3]<-"ALN"
nscodConf$obsCorStruct[3]<-"AR"
nscodConf$keyCorObs[3,]<-c(0,1,1,-1,-1)

par<-defpar(nscodData,nscodConf)

fit<-sam.fit(nscodData,nscodConf,par)

library(multiStockassessment,quietly=TRUE,warn.conflicts=FALSE)
cs<-suggestCorStructure(c(fit),nAgeClose=0)
mfit<-multisam.fit(c(fit),cs)
},type='message')
sink()

if(!all(abs(ssbtable(fit) - ssbtable(mfit)) < 1e-5 )){
    stop('ssbtables not equal')
}else{
    cat('ssbtable OK\n')
}

if(!all(abs(tsbtable(fit) - tsbtable(mfit)) < 1e-5 )){
    stop('tsbtables not equal')
}else{
    cat('tsbtable OK\n')
}

if(!all(abs(fbartable(fit) - fbartable(mfit)) < 1e-5 )){
    stop('fbartables not equal')
}else{
    cat('fbartable OK\n')
}

if(!all(abs(rectable(fit) - rectable(mfit)) < 1e-5 )){
    stop('rectables not equal')
}else{
    cat('rectable OK\n')
}

if(!all(abs(catchtable(fit) - catchtable(mfit)) < 1e-5 )){
    stop('catchtables not equal')
}else{
    cat('catchtable OK\n')
}

if(!all(abs(ntable(fit) - ntable(mfit)) < 1e-5 )){
    stop('ntables not equal')
}else{
    cat('ntable OK\n')
}

if(!all(abs(faytable(fit) - faytable(mfit)) < 1e-5 )){
    stop('faytables not equal')
}else{
    cat('faytable OK\n')
}

if(!all(abs(partable(fit) - partable(mfit)) < 1e-5 )){
    stop('partables not equal')
}else{
    cat('partable OK\n')
}


if(!all(abs(modeltable(fit) - modeltable(mfit)) < 1e-5 )){
    stop('modeltables not equal')
}else{
    cat('modeltable OK\n')
}

## if(!all(abs(ypr(fit) - ypr(mfit)) < 1e-6 )){
##     stop('modeltables not equal')
## }else{
##     cat('modeltable OK\n')
## }


if(!all(abs(summary(fit) - summary(mfit)) < 1e-5 )){
    stop('summaries not equal')
}else{
    cat('summary OK\n')
}
