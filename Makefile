doc:
	Rscript -e 'devtools::document("multiStockassessment")'

build: doc
	R CMD build multiStockassessment

check: doc build
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"multiStockassessment/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	R CMD check --as-cran multiStockassessment_${VERSION}.tar.gz

install: doc build
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"multiStockassessment/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	R CMD INSTALL multiStockassessment_${VERSION}.tar.gz
