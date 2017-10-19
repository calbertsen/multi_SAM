library(stockassessment)

data(nscodData)
data(nscodConf)
data(nscodParameters)
fit <- sam.fit(nscodData, nscodConf, nscodParameters,sim.condRE = FALSE)

datSim <- simulate(fit,nsim=2)
fitSim <- do.call("c",lapply(datSim,function(x)sam.fit(x,fit$conf, fit$pl)))


## Multiple stocks

library(multiStockassessment)

cs <- suggestCorStructure(fitSim,nAgeClose=0)
obj <- multisam.fit(fitSim,cs)

do.call("sum",lapply(fitSim,AIC)) - AIC(obj) < 1e-10
do.call("sum",lapply(fitSim,BIC)) - BIC(obj) < 1e-10
do.call("sum",lapply(fitSim,nobs)) - nobs(obj) < 1e-10
do.call("sum",lapply(fitSim,function(x)attr(logLik(x),"df"))) - attr(logLik(obj),"df") < 1e-10
do.call("sum",lapply(fitSim,logLik)) - as.numeric(logLik(obj)) < 1e-10

cs2 <- suggestCorStructure(fitSim,nAgeClose=3)
obj2 <- multisam.fit(fitSim,cs2)
AIC(obj); AIC(obj2)
BIC(obj); BIC(obj2)


ll1 <- logLik(obj)
ll2 <- logLik(obj2)
v <- abs(as.numeric(ll2 - ll1))
p <- abs(attr(ll2,"df") - attr(ll1,"df"))
pchisq(2*v,p,lower.tail=FALSE)


## One stock with correlation in survival

fitPC1 <- multisam.fit(c(fit),suggestCorStructure(c(fit),nAgeClose = 1,noCorInArea=FALSE))
fitPC2 <- multisam.fit(c(fit),suggestCorStructure(c(fit),nAgeClose = 2,noCorInArea=FALSE))
fitPC3 <- multisam.fit(c(fit),suggestCorStructure(c(fit),nAgeClose = 3,noCorInArea=FALSE))
fitPC4 <- multisam.fit(c(fit),suggestCorStructure(c(fit),nAgeClose = 4,noCorInArea=FALSE))
fitPC5 <- multisam.fit(c(fit),suggestCorStructure(c(fit),nAgeClose = 5,noCorInArea=FALSE))
fitPC6 <- multisam.fit(c(fit),suggestCorStructure(c(fit),nAgeClose = 6,noCorInArea=FALSE))

AIC(fitPC1)
AIC(fitPC2)
AIC(fitPC3)
AIC(fitPC4)
AIC(fitPC5)
AIC(fitPC6)
