# $id: $

if BUILD_HANDBOOK

#EPS_FILES = EPS/tutorial-problemconfiguration.eps 
TEXSOURCE = *.tex ./EPS/*.eps
MAINDOCS = ewoms-handbook.dvi ewoms-handbook.pdf ewoms-handbook.ps

DOCS = $(HTDOCS) $(MAINDOCS)

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/handbook

handbookdir=$(docdir)/handbook
handbook_DATA = ewoms-handbook.pdf

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += ewoms-handbook*.png 
CLEANFILES += $(DOCS)
CLEANFILES += $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	$CONVERT $< $@
ewoms-handbook.tex: $(EPS_FILES)

ewoms-handbook.pdf: $(TEXSOURCE)
	@LATEX@ ewoms-handbook.tex <&-
	@BIBTEX@ ewoms-handbook <&-
	@LATEX@ ewoms-handbook.tex <&-
	@LATEX@ ewoms-handbook.tex <&-
	@DVIPDF@ ewoms-handbook <&-

doc: ewoms-handbook.pdf

#dist-hook:
#	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
#	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

EXTRA_TEXINPUTS=$(top_srcdir)

else
doc:
	echo -n ""
endif
EXTRA_DIST= ewoms-handbook.pdf
include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/latex
