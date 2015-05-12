# CVS $Id: Makefile.am,v 1.7 2007-10-22 10:43:23 steger Exp $

INCLUDES = -I$(top_srcdir) -I$(includedir)
LIBS = @LIBS@ -lm
AM_CFLAGS = -Wall -ansi -pedantic @CFLAGS@
AM_LDFLAGS = @LDFLAGS@


bin_SCRIPTS = cs_add_ti_dp
bin_PROGRAMS = tinoco 
tinoco_SOURCES = tinoco.c tinoco.h

EXTRA_DIST = $(bin_SCRIPTS)
