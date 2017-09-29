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

AIC(obj)
nobs(obj)
logLik(obj)

do.call("sum",lapply(fitSim,function(x)x$opt$objective))
obj$opt$objective

cs2 <- suggestCorStructure(fitSim,nAgeClose=3)
obj2 <- multisam.fit(fitSim,cs2)
obj2$opt$objective


v <- -obj2$opt$objective + obj$opt$objective
p <- length(obj2$opt$par) - length(obj$opt$par)
1-pchisq(2*v,p)

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
