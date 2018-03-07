success_symbol <- '\u001B[0;32m\u263A\u001B[0;0m'
failure_symbol <- '\u001B[0;31m\u2639\u001B[0;0m'
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
    stop(sprintf('ssbtables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('ssbtable OK %s\n',success_symbol))
}

if(!all(abs(tsbtable(fit) - tsbtable(mfit)) < 1e-5 )){
    stop(sprintf('tsbtables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('tsbtable OK %s\n',success_symbol))
}

if(!all(abs(fbartable(fit) - fbartable(mfit)) < 1e-5 )){
    stop(sprintf('fbartables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('fbartable OK %s\n',success_symbol))
}

if(!all(abs(rectable(fit) - rectable(mfit)) < 1e-5 )){
    stop(sprintf('rectables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('rectable OK %s\n',success_symbol))
}

if(!all(abs(catchtable(fit) - catchtable(mfit)) < 1e-5 )){
    stop(sprintf('catchtables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('catchtable OK %s\n',success_symbol))
}

if(!all(abs(ntable(fit) - ntable(mfit)) < 1e-5 )){
    stop(sprintf('ntables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('ntable OK %s\n',success_symbol))
}

if(!all(abs(faytable(fit) - faytable(mfit)) < 1e-5 )){
    stop(sprintf('faytables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('faytable OK %s\n',success_symbol))
}

if(!all(abs(partable(fit) - partable(mfit)) < 1e-5 )){
    stop(sprintf('partables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('partable OK %s\n',success_symbol))
}


if(!all(abs(modeltable(fit) - modeltable(mfit)) < 1e-5 )){
    stop(sprintf('modeltables not equal %s\n',failure_symbol))
}else{
    cat(sprintf('modeltable OK %s\n',success_symbol))
}

ypr1 <- ypr(fit)
ypr2 <- ypr(mfit)[[1]]
yprnames <- names(ypr1)

a <- capture.output(tval <- !all(abs(print(ypr1) - print(ypr2)) < 1e-6 ))
if(tval){
    stop(sprintf('ypr print table not equal %s\n',failure_symbol))
}else{
    cat(sprintf('ypr print table OK %s\n',success_symbol))
}

for(nn in yprnames){
    if(is.numeric(ypr1[[nn]])){
        tval <- !all(abs(ypr1[[nn]] - ypr2[[nn]]) < 1e-6 )
    ## }else if(is.call(ypr1[[nn]])){
    ##     tval <- deparse(ypr1[[nn]]) != deparse(ypr2[[nn]])
    }else{
        tval <- ypr1[[nn]] != ypr2[[nn]]
    }
    if(tval){
        stop(sprintf('ypr %s not equal %s\n',nn,failure_symbol))
    }else{
        cat(sprintf('ypr %s OK %s\n',nn,success_symbol))
    }
}


if(!all(abs(summary(fit) - summary(mfit)) < 1e-5 )){
    stop(sprintf('summaries not equal %s\n',failure_symbol))
}else{
    cat(sprintf('summary OK %s\n',success_symbol))
}
