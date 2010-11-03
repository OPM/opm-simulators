# $id: $

if BUILD_HANDBOOK

#EPS_FILES = EPS/tutorial-problemconfiguration.eps 

DOCSOURCE = dumux-handbook.tex intro.tex propertysystem.tex navigation.html css.html

MAINDOCS = dumux-handbook.dvi dumux-handbook.pdf dumux-handbook.ps

if WML
 WMLDOCS = index.html
endif

if TEX4HT
 HTDOCS =  navigation.html css.html dumux-handbook.html
endif

EXTRA_DIST = dumux-handbook.pdf 
#Makefile.dist.in Makefile.dist.am

if BUILD_DOCS
if CONVERT
 DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif
endif


# setting like in dune-web
BASEDIR=../..
CURDIR=doc/handbook

handbookdir=$(docdir)/handbook
handbook_DATA = dumux-handbook.pdf

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.ps
CLEANFILES += *.4ct *.4tc *.css *.idv *.idx *.lg *.tid *.tmp *.tms *.xref
CLEANFILES += dumux-handbook*.html dumux-handbook*.png 
CLEANFILES += $(DOCS)
if CONVERT
CLEANFILES += $(EPS_FILES)

%.eps: $(srcdir)/%.jpg
	convert $< $@
dumux-handbook.tex: $(EPS_FILES)
endif

dumux-handbook.html:  $(DOCSOURCE) dune.cfg tex4ht.env
	$(MAKE) dumux-handbook.dvi
	$(MK4HT) htlatex dumux-handbook.tex "dune"; \
	$(MK4HT) htlatex dumux-handbook.tex "dune"; \
	$(T4HT) dumux-handbook.tex
	$(TEX) dumux-handbook.tex

dumux-handbook.pdf: 
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
