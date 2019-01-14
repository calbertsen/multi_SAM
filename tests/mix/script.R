success_symbol <- '\u001B[0;32m\u263A\u001B[0;0m'
failure_symbol <- '\u001B[0;31m\u2639\u001B[0;0m'

testFun <- function(fl,mf){
    if(!isTRUE(all.equal(Reduce("+",lapply(fl,AIC)),AIC(mf)))) stop(sprintf('AIC not equal %s vs %s',round(Reduce("+",lapply(fl,AIC)),7),round(AIC(mf),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
    if(!isTRUE(all.equal(Reduce("+",lapply(fl,logLik)),logLik(mf),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(Reduce("+",lapply(fl,logLik)),7),round(logLik(mf),7),failure_symbol))
    cat(sprintf('logLik OK %s\n',success_symbol))
    if(!isTRUE(all.equal(Reduce("+",lapply(fl,nobs)),nobs(mf),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(Reduce("+",lapply(fl,nobs)),7),round(nobs(mf),7),failure_symbol))
    cat(sprintf('nobs OK %s\n',success_symbol))
    if(length(fl) == 1){
        if(!isTRUE(all.equal(Reduce("+",lapply(fl,BIC)),BIC(mf)))) stop(sprintf('BIC not equal %s vs %s',round(Reduce("+",lapply(fl,BIC)),7),round(BIC(mf),7),failure_symbol))
        cat(sprintf('BIC OK %s\n',success_symbol))
    }
    if(!isTRUE(all.equal(Reduce("+",lapply(fl,function(f)attr(logLik(f),'df'))),attr(logLik(mf),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(Reduce("+",lapply(fl,function(f)attr(logLik(f),'df'))),7),round(attr(logLik(mf),'df'),7),failure_symbol))
    cat(sprintf('logLik df OK %s\n',success_symbol))
}

setwd(sub('script.R','',"tests/ebwb_cod/script.R"))
sink('/dev/null',type='output')
capture.output({
## Tests the combination of ALN likelihood with AR correlation structure
library(stockassessment)

sprat <- fitfromweb("sp2015tmb2")
cbh <- fitfromweb("cbh2015_tmb")
wbh <- fitfromweb("wbss_herring_2017_tmb")

## Refit with current stockassessment version
sprat <- sam.fit(sprat$data,sprat$conf,defpar(sprat$data,sprat$conf))
cbh <- sam.fit(cbh$data,cbh$conf,defpar(cbh$data,cbh$conf))
wbh <- sam.fit(wbh$data,wbh$conf,defpar(wbh$data,wbh$conf))

library(multiStockassessment,quietly=TRUE,warn.conflicts=FALSE)

},type='message')
sink()



### Sprat
cat("Baltic Sprat:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(sprat),nAgeClose=0)
    mfitS<-multisam.fit(c(sprat),cs)
},type='message')
sink()
testFun(c(sprat),mfitS)

### CBH
cat("Central Baltic Herring:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(cbh),nAgeClose=0)
    mfitC<-multisam.fit(c(cbh),cs)
},type='message')
sink()
testFun(c(cbh),mfitC)


### WBH
cat("Western Baltic Herring:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(wbh),nAgeClose=0)
    mfitW<-multisam.fit(c(wbh),cs)
},type='message')
sink()
testFun(c(wbh),mfitW)

### Combined
cat("Baltic Sprat and Herring Combined:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(sprat,cbh,wbh),nAgeClose=0)
    mfitSCW<-multisam.fit(c(sprat,cbh,wbh),cs)
},type='message')
sink()
testFun(c(sprat,cbh,wbh),mfitSCW)


