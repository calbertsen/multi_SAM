SAM_BRANCH?=master
R?=R

all: doc updateStockassessment build check test install

doc:
	@echo "\033[0;32mUpdating documentation\033[0;0m"
	$(R) -q -e 'devtools::document("multiStockassessment")'

build: doc
	@echo "\033[0;32mBuilding package\033[0;0m"
	$(R) CMD build multiStockassessment

updateStockassessment:
	@echo "\033[0;32mUpdating stockassessment from ${SAM_BRANCH} branch\033[0;0m"
	$(R) -q -e  "devtools::install_github(\"fishfollower/SAM\",subdir=\"stockassessment\",ref=\"${SAM_BRANCH}\")"

check: doc build
	@echo "\033[0;32mChecking package as cran\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"multiStockassessment/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	$(R) CMD check --as-cran multiStockassessment_${VERSION}.tar.gz

install: doc build
	@echo "\033[0;32mInstalling package\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"multiStockassessment/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	@echo "Version: ${VERSION}"
	$(R) CMD INSTALL multiStockassessment_${VERSION}.tar.gz

test:
	@echo "\033[0;31mNothing yet\033[0;0m"
