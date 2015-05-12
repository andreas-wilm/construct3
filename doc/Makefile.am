# Copyright (C) 2001-2004 Institute of Physical Biology,
#                    Heinrich-Heine Universität Düsseldorf, Germany
#                    <construct@biophys.uni-duesseldorf.de>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# CVS $Id: Makefile.am,v 1.8 2007-10-22 10:43:18 steger Exp $



man_MANS = construct.1

DOC_FILE = manual.pdf
docdir = ${datadir}/doc/@PACKAGE@
doc_DATA =  $(DOC_FILE)

#manual.pdf: bugs.tex copyright.tex indepth.tex \
#            introduction.tex title.tex contact.tex fixme.tex \
#            installation.tex manual.tex quickstart.tex \
#            this_bib.bib this_bibstyle.bst this_style.sty
#	pdflatex manual.tex
#	bibtex manual
#	pdflatex manual.tex
#	pdflatex manual.tex

EXTRA_DIST = $(DOC_FILE) $(man_MANS)
