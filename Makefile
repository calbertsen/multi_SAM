SAM_BRANCH?=master
R?=R
PACKAGE=multiStockassessment
VERSION = $(shell Rscript -e "l<-readLines(\"${PACKAGE}/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))")

TARBALL := ${PACKAGE}_${VERSION}.tar.gz
RFILES := $(wildcard ${PACKAGE}/R/*.R)
NAMESPACEFILE := ${PACKAGE}/NAMESPACE
RCHECK := ${PACKAGE}.Rcheck

.PHONY: all doc build check install test


all: doc updateStockassessment build check install test testmore_stockassessment readme

$(NAMESPACEFILE): $(RFILES)
	@echo "\033[0;32mUpdating documentation\033[0;0m"
	rm -f ${PACKAGE}/src/*.so
	$(R) -q -e 'roxygen2::roxygenize("${PACKAGE}")'
	@touch $@

doc: ${PACKAGE}/NAMESPACE

readme:
	@echo "\033[0;32mCreating README\033[0;0m"
	rm -r -f README_files
	$(R) -q -e 'rmarkdown::render("README.Rmd",rmarkdown::md_document(variant = "gfm"))'

$(TARBALL): $(NAMESPACEFILE)
	@echo "\033[0;32mBuilding package\033[0;0m"
	$(R) CMD build ${PACKAGE}


build: $(TARBALL)

updateStockassessment:
	@echo "\033[0;32mUpdating stockassessment from ${SAM_BRANCH} branch\033[0;0m"
	@echo "source('https://raw.githubusercontent.com/calbertsen/caMisc/master/caMisc/R/build_from_github.R'); \
	installFromGithub('fishfollower/SAM','${SAM_BRANCH}','stockassessment')" | $(R) -q --slave

$(RCHECK): $(TARBALL)
	@echo "\033[0;32mChecking package as cran\033[0;0m"
	$(R) CMD check --as-cran $(TARBALL)

check: $(RCHECK)

install: $(TARBALL)
	@echo "\033[0;32mInstalling package\033[0;0m"
	$(R) CMD INSTALL $(TARBALL)

test:
	@echo "\033[0;32mRunning tests\033[0;0m"
	@$(eval SCRIPTS := $(shell find tests -name script.R))
	@for scrpt in ${SCRIPTS} ; do \
		echo "\033[1;33mRunning $$scrpt\033[0;0m"; \
		$(R) --slave -f $$scrpt ;\
	done

clean_testmore_stockassessment:
	@echo "\033[0;32mCleaning testmore_stockassessment\033[0;0m"
	@rm -r -f testmore_stockassessment

prepare_testmore_stockassessment: clean_testmore_stockassessment
	@echo "\033[0;32mDownloading tests from stockassessment\033[0;0m"
	@mkdir testmore_stockassessment
	@curl -s -L https://github.com/fishfollower/SAM/archive/${SAM_BRANCH}.zip -o sam.zip > /dev/null
	@unzip sam.zip */testmore/* -d testmore_stockassessment/ > /dev/null
	@rm sam.zip
	@mv testmore_stockassessment/*/testmore/* testmore_stockassessment/
	@rm -r testmore_stockassessment/SAM-${SAM_BRANCH} testmore_stockassessment/revdep testmore_stockassessment/plots testmore_stockassessment/tables
	@echo "lf <- list.files('testmore_stockassessment',recursive=TRUE,full.names=TRUE,pattern='script.R'); \
	fn <- function(f){rl <- readLines(f);cat(c('success_symbol <- \'\u001B[0;32m\u263A\u001B[0;0m\'','failure_symbol <- \'\u001B[0;31m\u2639\u001B[0;0m\'',paste0('setwd(sub(\'script.R\',\'\',\"',f,'\"))'),'sink(\'/dev/null\',type=\'output\')','capture.output({',rl[!grepl('cat\\\\\\\\(',rl)],'library(multiStockassessment,quietly=TRUE,warn.conflicts=FALSE)','cs<-suggestCorStructure(c(fit),nAgeClose=0)','mfit<-multisam.fit(c(fit),cs)','},type=\'message\')','sink()', \
	'if(round(AIC(fit),7) != round(AIC(mfit),7)) stop(sprintf(\'AIC not equal %s vs %s %s\',round(AIC(fit),7),round(AIC(mfit),7),failure_symbol))','cat(sprintf(\'AIC OK %s\\\\\\\\n\',success_symbol))', \
	'if(round(logLik(fit),7) != round(logLik(mfit),7)) stop(sprintf(\'logLik not equal %s vs %s %s\',round(logLik(fit),7),round(logLik(mfit),7),failure_symbol))','cat(sprintf(\'logLik OK %s\\\\\\\\n\',success_symbol))', \
	'if(round(nobs(fit),7) != round(nobs(mfit),7)) stop(sprintf(\'nobs not equal %s vs %s %s\',round(nobs(fit),7),round(nobs(mfit),7),failure_symbol))','cat(sprintf(\'nobs OK %s\\\\\\\\n\',success_symbol))', \
	'if(round(BIC(fit),7) != round(BIC(mfit),7)) stop(sprintf(\'BIC not equal %s vs %s %s\',round(BIC(fit),7),round(BIC(mfit),7),failure_symbol))','cat(sprintf(\'BIC OK %s\\\\\\\\n\',success_symbol))', \
	'if(round(attr(logLik(fit),\'df\'),7) != round(attr(logLik(mfit),\'df\'),7)) stop(sprintf(\'logLik df not equal %s vs %s %s\',round(attr(logLik(fit),\'df\'),7),round(attr(logLik(mfit),\'df\'),7),failure_symbol))','cat(sprintf(\'logLik df OK %s\\\\\\\\n\',success_symbol))'), \
	sep='\n',file=f)};invisible(sapply(lf,fn))" | $(R) -q --slave

run_testmore_stockassessment: prepare_testmore_stockassessment
	@echo "\033[0;32mRunning tests from stockassessment\033[0;0m"
	@$(eval SCRIPTS := $(shell find testmore_stockassessment -name script.R))
	@for scrpt in ${SCRIPTS} ; do \
		echo "\033[1;33mRunning $$scrpt\033[0;0m"; \
		$(R) --slave -f $$scrpt ;\
	done

testmore_stockassessment: run_testmore_stockassessment

clean: clean_testmore_stockassessment
	@echo "\033[0;32mRemoving tar.gz files cleaning directory\033[0;0m"
	rm -f multiStockassessment_*.tar.gz
	@echo "\033[0;32mRemoving Rhistory files cleaning directory\033[0;0m"
	find . -type f -name '.Rhistory' -delete
	@echo "\033[0;32mRemoving Rcheck directory cleaning directory\033[0;0m"
	rm -f -r multiStockassessment.Rcheck


clean_hard:
	@echo "\033[0;32mHard cleaning directory\033[0;0m"
	git clean -f -d

uninstall:
	@echo "\033[0;32mUninstalling package\033[0;0m"
	$(R) CMD REMOVE multiStockassessment
