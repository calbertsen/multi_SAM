success_symbol <- '\u001B[0;32m\u263A\u001B[0;0m'
failure_symbol <- '\u001B[0;31m\u2639\u001B[0;0m'
#setwd(sub('script.R','',"tests/ebwb_cod/script.R"))
sink('/dev/null',type='output')
capture.output({
    library(stockassessment)

    eb <- stockassessment:::refit("EBC2014_final_tmb")

    library(multiStockassessment, ,quietly=TRUE,warn.conflicts=FALSE)

    mfitEB <- multisam.fit(c(eb))

    f1 <- modelforecast(eb, fval = rep(0.2,5), biasCorrect=FALSE)
    mf1 <- modelforecast(mfitEB, fval = rep(0.2,5), biasCorrect=FALSE)

},type='message')
sink()

if(!isTRUE(all.equal(AIC(eb),AIC(mfitEB)))){
    cat(sprintf('AIC not equal %s vs %s\n',round(AIC(eb),7),round(AIC(mfitEB),7),failure_symbol))
}else{
    cat(sprintf('AIC OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(logLik(eb),logLik(mfitEB),check.attributes=FALSE))){
    cat(sprintf('logLik not equal %s vs %s %s\n',round(logLik(eb),7),round(logLik(mfitEB),7),failure_symbol))
}else{
    cat(sprintf('logLik OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(nobs(eb),nobs(mfitEB),check.attributes=FALSE))){
    cat(sprintf('nobs not equal %s vs %s %s\n',round(nobs(eb),7),round(nobs(mfitEB),7),failure_symbol))
}else{
    cat(sprintf('nobs OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(BIC(eb),BIC(mfitEB)))){
    cat(sprintf('BIC not equal %s vs %s\n',round(BIC(eb),7),round(BIC(mfitEB),7),failure_symbol))
}else{
    cat(sprintf('BIC OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(attr(logLik(eb),'df'),attr(logLik(mfitEB),'df'),tolerance=1e-6))){
    cat(sprintf('logLik df not equal %s vs %s %s\n',round(attr(logLik(eb),'df'),7),round(attr(logLik(mfitEB),'df'),7),failure_symbol))
}else{
    cat(sprintf('logLik df OK %s\n',success_symbol))
}

tst <- all.equal(attr(f1,"shorttab"),attr(mf1[[1]],"shorttab"),tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Forecast shorttab not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Forecast shorttab OK %s\n',success_symbol))
}

tst <- all.equal(attr(f1,"tab"),attr(mf1[[1]],"tab"),tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Modelforecast tab not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Modelforecast tab OK %s\n',success_symbol))
}

tst <- all.equal(attr(f1,"shorttab"),attr(mf1[[1]],"shorttab"),tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Modelforecast fulltab not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Modelforecast fulltab OK %s\n',success_symbol))
}
