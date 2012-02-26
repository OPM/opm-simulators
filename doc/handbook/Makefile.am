# $id: $

if BUILD_HANDBOOK

#EPS_FILES = EPS/tutorial-problemconfiguration.eps 
TEXSOURCE = *.tex ./EPS/*.eps
MAINDOCS = dumux-handbook.dvi dumux-handbook.pdf dumux-handbook.ps
EXTRA_DIST = dumux-handbook.pdf $(TEXSOURCE) dumux-handbook.bib ModelDescriptions/*.tex
#Makefile.dist.in Makefile.dist.am

DOCS = $(HTDOCS) $(MAINDOCS)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/handbook

handbookdir=$(docdir)/handbook
handbook_DATA = dumux-handbook.pdf

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += dumux-handbook*.png 
CLEANFILES += $(DOCS)
CLEANFILES += $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	$CONVERT $< $@
dumux-handbook.tex: $(EPS_FILES)

dumux-handbook.pdf: $(TEXSOURCE)
	@LATEX@ dumux-handbook.tex <&-
	@BIBTEX@ dumux-handbook <&-
	@LATEX@ dumux-handbook.tex <&-
	@LATEX@ dumux-handbook.tex <&-
	@DVIPDF@ dumux-handbook <&-

#dist-hook:
#	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
#	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

EXTRA_TEXINPUTS=$(top_srcdir)

else
doc:
	echo -n ""
endif

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/latex
