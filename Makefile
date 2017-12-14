SAM_BRANCH?=master
R?=R

all: doc updateStockassessment build check test install

doc:
	@echo "\033[0;32mUpdating documentation\033[0;0m"
	rm -f multiStockassessment/src/*.so
	$(R) -q -e 'devtools::document("multiStockassessment")'

build: doc
	@echo "\033[0;32mBuilding package\033[0;0m"
	$(R) CMD build multiStockassessment

updateStockassessment:
	@echo "\033[0;32mUpdating stockassessment from ${SAM_BRANCH} branch\033[0;0m"
	$(R) -q -e  "devtools::install_github(\"fishfollower/SAM\",subdir=\"stockassessment\",ref=\"${SAM_BRANCH}\",force=TRUE)"

check: doc build
	@echo "\033[0;32mChecking package as cran\033[0;0m"
	$(eval VERSION = $(shell $(R) -q --slave "l<-readLines(\"multiStockassessment/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	$(R) CMD check --as-cran multiStockassessment_${VERSION}.tar.gz

install: doc build
	@echo "\033[0;32mInstalling package\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"multiStockassessment/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	@echo "Version: ${VERSION}"
	$(R) CMD INSTALL multiStockassessment_${VERSION}.tar.gz

test:
	@echo "\033[0;31mNothing yet\033[0;0m"

clean_testmore_stockassessment:
	@echo "\033[0;32mCleaning testmore_stockassessment\033[0;0m"
	@rm -r -f testmore_stockassessment

prepare_testmore_stockassessment: clean_testmore_stockassessment
	@echo "\033[0;32mDownloading tests from stockassessment\033[0;0m"
	@mkdir testmore_stockassessment
	@curl -s -L https://github.com/fishfollower/SAM/archive/${SAM_BRANCH}.zip -o sam.zip > /dev/null
	@unzip sam.zip */testmore/* -d testmore_stockassessment/ > /dev/null
	@mv testmore_stockassessment/*/testmore/* testmore_stockassessment/
	@rm -r testmore_stockassessment/SAM-${SAM_BRANCH} testmore_stockassessment/revdep
	@echo "lf <- list.files('testmore_stockassessment',recursive=TRUE,full.names=TRUE,pattern='script.R'); \
	fn <- function(f){rl <- readLines(f);cat(c(paste0('setwd(sub(\'script.R\',\'\',\"',f,'\"))'),'sink(\'/dev/null\',type=\'output\')','capture.output({',rl[!grepl('cat\\\\\\\\(',rl)],'library(multiStockassessment,quietly=TRUE,warn.conflicts=FALSE)','cs<-suggestCorStructure(c(fit),nAgeClose=0)','mfit<-multisam.fit(c(fit),cs)','},type=\'message\')','sink()', \
	'if(round(AIC(fit),7) != round(AIC(mfit),7)) stop(sprintf(\'AIC not equal %s vs %s\',round(AIC(fit),7),round(AIC(mfit),7)))','cat(\'AIC OK\\\\\\\\n\')', \
	'if(round(logLik(fit),7) != round(logLik(mfit),7)) stop(sprintf(\'logLik not equal %s vs %s\',round(logLik(fit),7),round(logLik(mfit),7)))','cat(\'logLik OK\\\\\\\\n\')', \
	'if(round(nobs(fit),7) != round(nobs(mfit),7)) stop(sprintf(\'nobs not equal %s vs %s\',round(nobs(fit),7),round(nobs(mfit),7)))','cat(\'nobs OK\\\\\\\\n\')', \
	'if(round(BIC(fit),7) != round(BIC(mfit),7)) stop(sprintf(\'BIC not equal %s vs %s\',round(BIC(fit),7),round(BIC(mfit),7)))','cat(\'BIC OK\\\\\\\\n\')', \
	'if(round(attr(logLik(fit),\'df\'),7) != round(attr(logLik(mfit),\'df\'),7)) stop(sprintf(\'logLik df not equal %s vs %s\',round(attr(logLik(fit),\'df\'),7),round(attr(logLik(mfit),\'df\'),7)))','cat(\'logLik df OK\\\\\\\\n\')'), \
	sep='\n',file=f)};invisible(sapply(lf,fn))" | $(R) -q --slave

run_testmore_stockassessment: prepare_testmore_stockassessment
	@echo "\033[0;32mRunning tests from stockassessment\033[0;0m"
	@$(eval SCRIPTS := $(shell find testmore_stockassessment -name script.R))
	@for scrpt in ${SCRIPTS} ; do \
		echo "\033[1;33mRunning $$scrpt\033[0;0m"; \
		$(R) --slave -f $$scrpt ;\
	done

testmore_stockassessment: run_testmore_stockassessment

clean:
	@echo "\033[0;32mCleaning directory\033[0;0m"
	git clean -f -d

uninstall:
	@echo "\033[0;32mUninstalling package\033[0;0m"
	$(R) CMD REMOVE multiStockassessment
