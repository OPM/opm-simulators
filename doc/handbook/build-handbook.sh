#! /bin/sh

# this script build the eWoms handbook from its LaTeX sources. The
# result file is called "ewoms-handbook.pdf"

latex ewoms-handbook
bibtex ewoms-handbook
latex ewoms-handbook
latex ewoms-handbook
dvipdf ewoms-handbook
rm ewoms-handbook.dvi
