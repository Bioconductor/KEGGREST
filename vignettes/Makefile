
whoami=$(shell whoami)

all: KEGGREST-vignette.pdf clean

pdfclean:
	rm -f KEGGREST-vignette.pdf

clean:
	rm -rf *.tex *.bbl *.blg *.aux *.out *.log cache *.R
	rm -rf figure
	rm -rf *.md

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript


KEGGREST-vignette.pdf:
	rm -f *.R
	$(RSCRIPT) --vanilla -e "library(knitr);knit('KEGGREST-vignette.Rmd')"
	$(RSCRIPT) --vanilla -e "library(knitr);purl('KEGGREST-vignette.Rmd')"
	pandoc KEGGREST-vignette.pdf KEGGREST-vignette.md
	cp KEGGREST-vignette.md ../inst/doc
	cp KEGGREST-vignette.R ../inst/doc/KEGGREST-vignette.R.txt

	

