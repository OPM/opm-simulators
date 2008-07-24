# $id: $

EPS_FILES = 

DOCSOURCE = dumux-handbook.tex navigation.html css.html

MAINDOCS = dumux-handbook.dvi dumux-handbook.pdf dumux-handbook.ps

if WML
 WMLDOCS = index.html
endif

if TEX4HT
 HTDOCS =  navigation.html css.html dumux-handbook.html
endif

EXTRA_DIST = dumux-handbook.pdf Makefile.dist.in Makefile.dist.am

if BUILD_DOCS
if CONVERT
 DOCS = $(HTDOCS) $(MAINDOCS) $(WMLDOCS)
endif
endif

# setting like in dune-web
BASEDIR=../..
CURDIR=doc/handbook

docdir=$(datadir)/doc/handbook

doc_DATA = dumux-handbook.pdf

EXTRAINSTALL = dumux-handbook.pdf grid-howto*.html *.png *.css

CLEANFILES = *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.pdf *.ps
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

dist-hook:
	sed $(srcdir)/Makefile.dist.am -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.am
	sed $(srcdir)/Makefile.dist.in -e 's/Makefile\.dist/Makefile/g' > $(distdir)/Makefile.in

EXTRA_TEXINPUTS:=$(top_srcdir)

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/webstuff
include $(top_srcdir)/am/latex
