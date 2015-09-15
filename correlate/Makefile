#######################
# general information #
#######################

# file:         Makefile
# created:      2015-09-14
# last update:  2015-09-15
# author(s):    Marcel Schilling <marcel.schilling@mdc-berlin.de>
# purpose:      automate correlation sub-project of lncRNA-splice project


######################################
# change log (reverse chronological) #
######################################

# 2015-09-15: switched from example alignments to real data
# 2015-09-14: added alignments example TSV file
#             initial version (building HTML report based on coverage example TSV file)


####################
# path definitions #
####################

# (absolute) path of this Makefile
MAKEFILE:=$(realpath $(lastword $(MAKEFILE_LIST)))

# (absolute) path of the directory containing this Makefile
MAKEFILE_DIRECTORY:=$(dir $(MAKEFILE))

# (absolute) path of TSV file with intron/exon/lncRNA coverage
COVERAGE_TSV:=$(MAKEFILE_DIRECTORY)coverage.example.tsv

# (absolute) path of TSV file with lncRNA:splice-site alignments
ALIGNMENTS_TSV:=$(MAKEFILE_DIRECTORY)alignments.tsv.gz

# (absolute) path of R-Markdown file used to generate the report
REPORT_RMD:=$(MAKEFILE_DIRECTORY)report.rmd

# (absolute) path of Markdown file containing the report
REPORT_MD:=$(REPORT_RMD:.rmd=.md)

# (absolute) path of CSS-style-sheet file used for the HTML version of the report
MARKDOWN_STYLESHEET:=$(MAKEFILE_DIRECTORY)marta.css

# (absolute) path of HTML file containing the report
REPORT_HTML:=$(REPORT_MD:.md=.html)


################
# make options #
################

.DELETE_ON_ERROR :
.SUFFIXES:
.SECONDARY:
SHELL:=/bin/bash


#######################
# command definitions #
#######################

# command used to run R commands
RUN_R_COMMAND:=Rscript -e


##################
# common targets #
##################

# if no target was specified, get HTML report
.PHONY : all
all : $(REPORT_HTML)


###################
# generate report #
###################

# define helper macros (see http://blog.jgc.org/2007/06/escaping-comma-and-space-in-gnu-make.html)
COMMA:=,
SPACE:=
SPACE+=

# define multi-line report parameters (must not include single-quotes)
define REPORT_PARAMS
  coverage_tsv<-"$(COVERAGE_TSV)"
  alignments_tsv<-"$(ALIGNMENTS_TSV)"
endef

# define helper variable to use multi-line variable as multi-line string
define newline


endef

# replace newline by semicolon
REPORT_PARAMS:=$(subst $(newline),;,${REPORT_PARAMS})

# knit report to Markdown passing parameters
$(REPORT_MD) : $(REPORT_RMD) $(COVERAGE_TSV) $(ALIGNMENTS_TSV) | $(dir $(REPORT_MD))
	$(RUN_R_COMMAND) '$(REPORT_PARAMS);require(knitr);knit("$<",output="$@")'

# convert Markdown to HTML
%.html : %.md $(MARKDOWN_STYLESHEET)
	$(RUN_R_COMMAND) 'require(markdown);markdownToHTML("$<","$@",stylesheet="$(word 2,$^)",options=c(markdownHTMLOptions(defaults=T),"toc"))'
