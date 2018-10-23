success_symbol <- '\u001B[0;32m\u263A\u001B[0;0m'
failure_symbol <- '\u001B[0;31m\u2639\u001B[0;0m'
setwd(sub('script.R','',"tests/ebwb_cod/script.R"))
sink('/dev/null',type='output')
capture.output({
## Tests the combination of ALN likelihood with AR correlation structure
library(stockassessment)

wb <- fitfromweb("wbcod_2014_tmb")
eb <- fitfromweb("EBC2014_final_tmb")

library(multiStockassessment,quietly=TRUE,warn.conflicts=FALSE)

## Single wb fit
cs<-suggestCorStructure(c(wb),nAgeClose=0)
mfitWB<-multisam.fit(c(wb),cs)
## Single eb fit
cs<-suggestCorStructure(c(eb),nAgeClose=0)
mfitEB<-multisam.fit(c(eb),cs)
## Multifit 0
cs<-suggestCorStructure(c(eb,wb),nAgeClose=0)
mfitEBWB<-multisam.fit(c(eb,wb),cs)
## Multifit 1
cs<-suggestCorStructure(c(eb,wb),nAgeClose=1)
mfitEBWB1<-multisam.fit(c(eb,wb),cs)
},type='message')
sink()

### EB
cat("Eastern Baltic Cod:\n")
if(!isTRUE(all.equal(AIC(eb),AIC(mfitEB)))) stop(sprintf('AIC not equal %s vs %s %s',round(AIC(eb),7),round(AIC(mfitEB),7),failure_symbol))
cat(sprintf('AIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(logLik(eb),logLik(mfitEB),check.attributes=FALSE))) stop(sprintf('logLik not equal %s vs %s %s',round(logLik(eb),7),round(logLik(mfitEB),7),failure_symbol))
cat(sprintf('logLik OK %s\n',success_symbol))
if(!isTRUE(all.equal(nobs(eb),nobs(mfitEB),check.attributes=FALSE))) stop(sprintf('nobs not equal %s vs %s %s',round(nobs(eb),7),round(nobs(mfitEB),7),failure_symbol))
cat(sprintf('nobs OK %s\n',success_symbol))
if(!isTRUE(all.equal(BIC(eb),BIC(mfitEB)))) stop(sprintf('BIC not equal %s vs %s %s',round(BIC(eb),7),round(BIC(mfitEB),7),failure_symbol))
cat(sprintf('BIC OK %s\n',success_symbol))
if(!isTRUE(all.equal(attr(logLik(eb),'df'),attr(logLik(mfitEB),'df'),tolerance=1e-6))) stop(sprintf('logLik df not equal %s vs %s %s',round(attr(logLik(eb),'df'),7),round(attr(logLik(mfitEB),'df'),7),failure_symbol))
cat(sprintf('logLik df OK %s\n',success_symbol))
### WB
cat("Western Baltic Cod:\n")
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
## Both
cat("Eastern and Western Baltic Cod Combined:\n")
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
cat("Eastern and Western Baltic Cod Combined and Correlated:\n")
if(!(logLik(mfitEBWB) < logLik(mfitEBWB1))) stop("log-likelihood did not increase %s vs %s",logLik(mfitEBWB),logLik(mfitEBWB1))
cat(sprintf('log-likelihood increased OK %s\n',success_symbol))
if(!(attr(logLik(mfitEBWB),'df') < attr(logLik(mfitEBWB1),'df'))) stop("Number of parameters did not increase %s vs %s",attr(logLik(mfitEBWB),'df'),attr(logLik(mfitEBWB1),'df'))
cat(sprintf('Number of parameters increased OK %s\n',success_symbol))
