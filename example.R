library(stockassessment)

data(nscodData)
data(nscodConf)
data(nscodParameters)
## gctorture()
fit <- sam.fit(nscodData, nscodConf, nscodParameters,sim.condRE = FALSE)
## gctorture(FALSE)

datSim <- simulate(fit,nsim=2)
fitSim <- do.call("c",lapply(datSim,function(x)sam.fit(x,nscodConf, nscodParameters)))


## Multiple stocks

library(multiStockassessment)

cs <- suggestCorStructure(fitSim,nAgeClose=0)
obj <- multisam.fit(fitSim,cs)

fbarplot(obj)
ssbplot(obj)
tsbplot(obj)
recplot(obj)
catchplot(obj)
srplot(obj)
fitplot(obj,1)
fitplot(obj,2)
parplot(obj)


unique(rownames(summary(attr(obj,"m_sdrep"))))
unique(rownames(summary(fit$sdrep)))

unique(rownames(summary(attr(obj,"m_rep"))))
unique(rownames(summary(fit$rep)))


nc <- corplot(obj)

do.call("sum",lapply(fitSim,AIC)) - AIC(obj) < 1e-10
do.call("sum",lapply(fitSim,BIC)) - BIC(obj) < 1e-10
do.call("sum",lapply(fitSim,nobs)) - nobs(obj) < 1e-10
do.call("sum",lapply(fitSim,function(x)attr(logLik(x),"df"))) - attr(logLik(obj),"df") < 1e-10
do.call("sum",lapply(fitSim,logLik)) - as.numeric(logLik(obj)) < 1e-10

cs2 <- suggestCorStructure(fitSim,nAgeClose=2)
obj2 <- multisam.fit(fitSim,cs2,usePartialCors=TRUE)
AIC(obj); AIC(obj2)
BIC(obj); BIC(obj2)

nc <- corplot(obj2)


modeltable.msam <- function(object,...){
    cat("Hello\n")
}

modeltable.msamlist <- function(...){
    cat("Hello list\n")
    x <- list(...)
    lapply(x,modeltable)
}


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

names(fitPC1[[1]]$rep)

names(attr(fitPC1,"m_rep"))

