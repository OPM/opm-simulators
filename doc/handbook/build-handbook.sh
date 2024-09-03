#! /bin/sh

# this script builds the eWoms handbook from its LaTeX sources. The
# result file is called "ewoms-handbook.pdf"

latex ewoms-handbook || exit $?
bibtex ewoms-handbook || exit $?
latex ewoms-handbook || exit $?
latex ewoms-handbook || exit $?
dvipdf ewoms-handbook || exit $?
rm ewoms-handbook.dvi
