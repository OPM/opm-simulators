# $id: $

if BUILD_HANDBOOK

#EPS_FILES = EPS/tutorial-problemconfiguration.eps 
TEXSOURCE = *.tex ./EPS/*.eps
MAINDOCS = dumux-handbook.dvi dumux-handbook.pdf dumux-handbook.ps
EXTRA_DIST = dumux-handbook.pdf 
#Makefile.dist.in Makefile.dist.am

if CONVERT
DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/handbook

handbookdir=$(docdir)/handbook
handbook_DATA = dumux-handbook.pdf

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += dumux-handbook*.png 
CLEANFILES += $(DOCS)
if CONVERT
CLEANFILES += $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	convert $< $@
dumux-handbook.tex: $(EPS_FILES)
endif

dumux-handbook.pdf: $(TEXSOURCE)
	latex dumux-handbook.tex
	bibtex dumux-handbook
	latex dumux-handbook.tex
	latex dumux-handbook.tex
	dvipdf dumux-handbook

#dist-hook:
#	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
#	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

EXTRA_TEXINPUTS=$(top_srcdir)

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/latex

endif
