all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(methods); library(knitrBootstrap); knitrBoostrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'; if(knitrBoostrapFlag) knit_bootstrap('GenomationManual-knitr.Rmd', show_code=TRUE, boot_style='Readable') else {library(rmarkdown); render('genomation.Rmd', bootstrap_document(theme.chooser=TRUE, highlight.chooser=TRUE))}"

clean:
	rm -f *.md