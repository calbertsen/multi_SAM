success_symbol <- '\u001B[0;32m\u263A\u001B[0;0m'
failure_symbol <- '\u001B[0;31m\u2639\u001B[0;0m'
#setwd(sub('script.R','',"tests/ebwb_cod/script.R"))
sink('/dev/null',type='output')
capture.output({
    library(stockassessment)

    eb <- stockassessment:::refit("BW_2018")

    library(multiStockassessment, ,quietly=TRUE,warn.conflicts=FALSE)

    mfitEB <- multisam.fit(c(eb))

    rp1 <- referencepoints(eb)
    mrp1 <- referencepoints(mfitEB)

},type='message')
sink()

cat("\n\nBW_2018 RW recruitment\n\n")

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

tst <- all.equal(rp1$tables$F,mrp1[[1]]$tables$F,tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint F table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint F table OK %s\n',success_symbol))
}


tst <- all.equal(rp1$tables$Yield,mrp1[[1]]$tables$Yield,tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint Yield table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint Yield table OK %s\n',success_symbol))
}


tst <- all.equal(rp1$tables$YieldPerRecruit,mrp1[[1]]$tables$YieldPerRecruit,tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint YieldPerRecruit table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint YieldPerRecruit table OK %s\n',success_symbol))
}


tst <- all.equal(rp1$tables$SpawnersPerRecruit,mrp1[[1]]$tables$SpawnersPerRecruit,tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint SpawnersPerRecruit table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint SpawnersPerRecruit table OK %s\n',success_symbol))
}


tst <- all.equal(rp1$tables$Biomass,mrp1[[1]]$tables$Biomass,tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint Biomass table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint Biomass table OK %s\n',success_symbol))
}


tst <- all.equal(rp1$tables$Recruitment,mrp1[[1]]$tables$Recruitment,tolerance=1e-5)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint Recruitment table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint Recruitment table OK %s\n',success_symbol))
}




## With recruitment model

sink('/dev/null',type='output')
capture.output({

    eb2 <- stockassessment:::refit("BW_2018", list(stockRecruitmentModelCode = 1))
    mfitEB2 <- multisam.fit(c(eb2))

    rp2 <- referencepoints(eb2)
    mrp2 <- referencepoints(mfitEB2)

},type='message')
sink()

cat("\n\nBW_2018 Ricker recruitment\n\n")

if(!isTRUE(all.equal(AIC(eb2),AIC(mfitEB2)))){
    cat(sprintf('AIC not equal %s vs %s\n',round(AIC(eb2),7),round(AIC(mfitEB2),7),failure_symbol))
}else{
    cat(sprintf('AIC OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(logLik(eb2),logLik(mfitEB2),check.attributes=FALSE))){
    cat(sprintf('logLik not equal %s vs %s %s\n',round(logLik(eb2),7),round(logLik(mfitEB2),7),failure_symbol))
}else{
    cat(sprintf('logLik OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(nobs(eb2),nobs(mfitEB2),check.attributes=FALSE))){
    cat(sprintf('nobs not equal %s vs %s %s\n',round(nobs(eb2),7),round(nobs(mfitEB2),7),failure_symbol))
}else{
    cat(sprintf('nobs OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(BIC(eb2),BIC(mfitEB2)))){
    cat(sprintf('BIC not equal %s vs %s\n',round(BIC(eb2),7),round(BIC(mfitEB2),7),failure_symbol))
}else{
    cat(sprintf('BIC OK %s\n',success_symbol))
}
if(!isTRUE(all.equal(attr(logLik(eb2),'df'),attr(logLik(mfitEB2),'df'),tolerance=1e-6))){
    cat(sprintf('logLik df not equal %s vs %s %s\n',round(attr(logLik(eb2),'df'),7),round(attr(logLik(mfitEB2),'df'),7),failure_symbol))
}else{
    cat(sprintf('logLik df OK %s\n',success_symbol))
}

tst <- all.equal(rp2$tables$F,mrp2[[1]]$tables$F,tolerance=1e-5, check.attributes = FALSE)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint F table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint F table OK %s\n',success_symbol))
}


tst <- all.equal(rp2$tables$Yield,mrp2[[1]]$tables$Yield,tolerance=1e-5, check.attributes = FALSE)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint Yield table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint Yield table OK %s\n',success_symbol))
}


tst <- all.equal(rp2$tables$YieldPerRecruit,mrp2[[1]]$tables$YieldPerRecruit,tolerance=1e-5, check.attributes = FALSE)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint YieldPerRecruit table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint YieldPerRecruit table OK %s\n',success_symbol))
}


tst <- all.equal(rp2$tables$SpawnersPerRecruit,mrp2[[1]]$tables$SpawnersPerRecruit,tolerance=1e-5, check.attributes = FALSE)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint SpawnersPerRecruit table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint SpawnersPerRecruit table OK %s\n',success_symbol))
}


tst <- all.equal(rp2$tables$Biomass,mrp2[[1]]$tables$Biomass,tolerance=1e-5, check.attributes = FALSE)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint Biomass table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint Biomass table OK %s\n',success_symbol))
}


tst <- all.equal(rp2$tables$Recruitment,mrp2[[1]]$tables$Recruitment,tolerance=1e-5, check.attributes = FALSE)
if(!isTRUE(tst)){
    cat(sprintf('Referencepoint Recruitment table not equal: %s %s\n',tst,failure_symbol))
}else{
    cat(sprintf('Referencepoint Recruitment table OK %s\n',success_symbol))
}


