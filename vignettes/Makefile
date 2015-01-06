all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitrBootstrap); knit_bootstrap('GenomationManual-knitr.Rmd', show_code=TRUE, boot_style='Readable')"
clean:
	rm -f *.md