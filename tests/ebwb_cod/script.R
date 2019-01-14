success_symbol <- '\u001B[0;32m\u263A\u001B[0;0m'
failure_symbol <- '\u001B[0;31m\u2639\u001B[0;0m'
setwd(sub('script.R','',"tests/ebwb_cod/script.R"))
sink('/dev/null',type='output')
capture.output({
## Tests the combination of ALN likelihood with AR correlation structure
library(stockassessment)

eb <- fitfromweb("EBC2014_final_tmb")
wb <- fitfromweb("wbcod_2014_tmb")

## Refit with current stockassessment version
eb <- sam.fit(eb$data,eb$conf,defpar(eb$data,eb$conf))
wb <- sam.fit(wb$data,wb$conf,defpar(wb$data,wb$conf))


dat <- eb$data
conf <- eb$conf
conf$corFlag <- 2
conf$keyVarObs <- matrix(c(1,2,2,2,2,2,2,
                      3,4,4,4,4,0,0,
                      5,6,6,6,6,0,0,
                      7,8,8,8,0,0,0,
                      0,9,10,10,10,10,0),5,7,byrow=TRUE)-1
par <- defpar(dat,conf)
eb2 <- sam.fit(dat,conf,par)

dat <- wb$dat
conf <- wb$conf
conf$corFlag <- 2
conf$keyVarObs <- matrix(c(0,1,2,2,2,2,2,2,
                           3,4,5,6,7,7,0,0,
                           0,0,0,8,8,9,9,0,
                           0,10,11,12,12,13,13,0),4,8,byrow=TRUE)-1
par <- defpar(dat,conf)
wb2 <- sam.fit(dat,conf,par)

library(multiStockassessment,quietly=TRUE,warn.conflicts=FALSE)

},type='message')
sink()

### EB
cat("Eastern Baltic Cod Configuration 1:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(eb),nAgeClose=0)
    mfitEB<-multisam.fit(c(eb),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(eb),AIC(mfitEB)))) stop(sprintf('AIC not equal %s vs %s',round(AIC(eb),7),round(AIC(mfitEB),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb),logLik(mfitEB),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb),7),round(logLik(mfitEB),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb),nobs(mfitEB),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb),7),round(nobs(mfitEB),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
if(!isTRUE(all.equal(BIC(eb),BIC(mfitEB)))) stop(sprintf('BIC not equal %s vs %s',round(BIC(eb),7),round(BIC(mfitEB),7),failure_symbol))
cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb),'df'),attr(logLik(mfitEB),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb),'df'),7),round(attr(logLik(mfitEB),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))

### EB
cat("Eastern Baltic Cod Configuration 2:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(eb2),nAgeClose=0)
    mfitEB2<-multisam.fit(c(eb2),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(eb2),AIC(mfitEB2)))) stop(sprintf('AIC not equal %s vs %s',round(AIC(eb2),7),round(AIC(mfitEB2),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb2),logLik(mfitEB2),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb2),7),round(logLik(mfitEB2),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb2),nobs(mfitEB2),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb2),7),round(nobs(mfitEB2),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
if(!isTRUE(all.equal(BIC(eb2),BIC(mfitEB2)))) stop(sprintf('BIC not equal %s vs %s',round(BIC(eb2),7),round(BIC(mfitEB2),7),failure_symbol))
cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb2),'df'),attr(logLik(mfitEB2),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb2),'df'),7),round(attr(logLik(mfitEB2),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))

### WB 1
cat("Western Baltic Cod Configuration 1:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(wb),nAgeClose=0)
    mfitWB<-multisam.fit(c(wb),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(wb),AIC(mfitWB)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(wb),7),round(AIC(mfitWB),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(wb),logLik(mfitWB),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(wb),7),round(logLik(mfitWB),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(wb),nobs(mfitWB),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(wb),7),round(nobs(mfitWB),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
if(!isTRUE(all.equal(BIC(wb),BIC(mfitWB)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(wb),7),round(BIC(mfitWB),7),failure_symbol))
cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(wb),'df'),attr(logLik(mfitWB),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(wb),'df'),7),round(attr(logLik(mfitWB),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))

### WB 2
cat("Western Baltic Cod Configuration 2:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(wb2),nAgeClose=0)
    mfitWB2<-multisam.fit(c(wb2),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(wb2),AIC(mfitWB2)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(wb2),7),round(AIC(mfitWB2),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(wb2),logLik(mfitWB2),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(wb2),7),round(logLik(mfitWB2),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(wb2),nobs(mfitWB2),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(wb2),7),round(nobs(mfitWB2),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
if(!isTRUE(all.equal(BIC(wb2),BIC(mfitWB2)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(wb),7),round(BIC(mfitWB),7),failure_symbol))
cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(wb2),'df'),attr(logLik(mfitWB2),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(wb2),'df'),7),round(attr(logLik(mfitWB2),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))

## EBWB 1 
cat("Eastern and Western Baltic Cod Combined 1:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(eb,wb),nAgeClose=0)
    mfitEBWB<-multisam.fit(c(eb,wb),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(eb)+AIC(wb),AIC(mfitEBWB)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(eb)+AIC(wb),7),round(AIC(mfitEBWB),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb)+logLik(wb),logLik(mfitEBWB),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb)+logLik(wb),7),round(logLik(mfitEBWB),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb)+nobs(wb),nobs(mfitEBWB),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb)+nobs(wb),7),round(nobs(mfitEBWB),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
##if(!isTRUE(all.equal(BIC(eb)+BIC(wb),BIC(mfitEBWB)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(eb)+BIC(wb),7),round(BIC(mfitEBWB),7),failure_symbol))
##cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb),'df')+attr(logLik(wb),'df'),attr(logLik(mfitEBWB),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb),'df')+attr(logLik(wb),'df'),7),round(attr(logLik(mfitEBWB),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))


## EBWB 2
cat("Eastern and Western Baltic Cod Combined 2:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(eb2,wb2),nAgeClose=0)
    mfitEB2WB2<-multisam.fit(c(eb2,wb2),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(eb2)+AIC(wb2),AIC(mfitEB2WB2)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(eb2)+AIC(wb2),7),round(AIC(mfitEB2WB2),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb2)+logLik(wb2),logLik(mfitEB2WB2),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb2)+logLik(wb2),7),round(logLik(mfitEB2WB2),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb2)+nobs(wb2),nobs(mfitEB2WB2),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb2)+nobs(wb2),7),round(nobs(mfitEB2WB2),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
##if(!isTRUE(all.equal(BIC(eb)+BIC(wb),BIC(mfitEBWB)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(eb)+BIC(wb),7),round(BIC(mfitEBWB),7),failure_symbol))
##cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb2),'df')+attr(logLik(wb2),'df'),attr(logLik(mfitEB2WB2),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb2),'df')+attr(logLik(wb2),'df'),7),round(attr(logLik(mfitEB2WB2),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))


## WBEB 1 
cat("Western and Eastern Baltic Cod Combined 1:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(wb,eb),nAgeClose=0)
    mfitWBEB<-multisam.fit(c(wb,eb),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(eb)+AIC(wb),AIC(mfitWBEB)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(eb)+AIC(wb),7),round(AIC(mfitWBEB),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb)+logLik(wb),logLik(mfitWBEB),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb)+logLik(wb),7),round(logLik(mfitWBEB),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb)+nobs(wb),nobs(mfitWBEB),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb)+nobs(wb),7),round(nobs(mfitWBEB),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
##if(!isTRUE(all.equal(BIC(eb)+BIC(wb),BIC(mfitEBWB)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(eb)+BIC(wb),7),round(BIC(mfitEBWB),7),failure_symbol))
##cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb),'df')+attr(logLik(wb),'df'),attr(logLik(mfitWBEB),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb),'df')+attr(logLik(wb),'df'),7),round(attr(logLik(mfitWBEB),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))


## WBEB 2
cat("Western and Eastern Baltic Cod Combined 1:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(wb2,eb2),nAgeClose=0)
    mfitWB2EB2<-multisam.fit(c(wb2,eb2),cs)
},type='message')
sink()   
if(!isTRUE(all.equal(AIC(eb2)+AIC(wb2),AIC(mfitWB2EB2)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(eb2)+AIC(wb2),7),round(AIC(mfitWB2EB2),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb2)+logLik(wb2),logLik(mfitWB2EB2),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb2)+logLik(wb2),7),round(logLik(mfitWB2EB2),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb2)+nobs(wb2),nobs(mfitWB2EB2),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb2)+nobs(wb2),7),round(nobs(mfitWB2EB2),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
##if(!isTRUE(all.equal(BIC(eb)+BIC(wb),BIC(mfitEBWB)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(eb)+BIC(wb),7),round(BIC(mfitEBWB),7),failure_symbol))
##cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb2),'df')+attr(logLik(wb2),'df'),attr(logLik(mfitWB2EB2),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb2),'df')+attr(logLik(wb2),'df'),7),round(attr(logLik(mfitWB2EB2),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))




cat("Eastern and Western Baltic Cod Combined and Correlated:\n")
sink('/dev/null',type='output')
capture.output({
    cs<-suggestCorStructure(c(wb2,eb2),nAgeClose=1)
    mfitWB2EB2_1 <-multisam.fit(c(wb2,eb2),cs, newtonsteps = 0)
},type='message')
sink()   
if(!(logLik(mfitWB2EB2) < logLik(mfitWB2EB2_1))) stop(sprintf("log-likelihood did not increase %s vs %s",logLik(mfitWB2EB2),logLik(mfitWB2EB2_1)))
cat(sprintf('log-likelihood increased OK %s\n',success_symbol))
if(!(attr(logLik(mfitWB2EB2),'df') < attr(logLik(mfitWB2EB2_1),'df'))) stop(sprintf("Number of parameters did not increase %s vs %s",attr(logLik(mfitWB2EB2),'df'),attr(logLik(mfitWB2EB2_1),'df')))
cat(sprintf('Number of parameters increased OK %s\n',success_symbol))
