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

# CVS $Id: Makefile.am,v 1.19 2007-09-30 13:59:25 wilm Exp $

CC=@TCL_CC@
INCLUDES = -I$(includedir) @TCL_INCLUDE_SPEC@ @TK_INCLUDE_SPEC@ @TK_XINCLUDES@
AM_CFLAGS   = -ansi -Wall @TCL_CFLAGS_OPTIMIZE@ @TCL_SHLIB_CFLAGS@ @MACOS_CFLAG@ @CFLAGS@
# @TCL_DEFS@ @TK_DEFS@
# TCL_CFLAGS_WARNING:@TCL_CFLAGS_WARNING@
LIBS     = @TCL_LIB_SPEC@ @TCL_LIBS@ @TCL_STUB_LIB_SPEC@ @TK_LIB_SPEC@ @TK_LIBS@ @TK_STUB_LIB_SPEC@
# @LIBS@
# LDFLAGS  = @TCL_LD_FLAGS@ @TCL_LDFLAGS_OPTIMIZE@ @LDFLAGS@



DRAWSTRUCT_TCL_FILES =  drawstruct.tcl draw_test.tcl mkIndex.tcl

lib_LTLIBRARIES = libdrawstructcore.la

libdrawstructcore_la_SOURCES = CS_DrawStruct.c \
	CS_DrawStruct.h \
	TkCS_DrawStruct.c \
	TkCS_DrawStruct.h \
	drawstruct.h \
	generic_utils.c \
	generic_utils.h

# "@" is needed here, otherwise it's a self reference
libdir        =  @libdir@/@PACKAGE@/drawstruct/
drawstructdir = $(libdir)


EXTRA_DIST      = $(DRAWSTRUCT_TCL_FILES)
drawstruct_DATA = $(DRAWSTRUCT_TCL_FILES)

# autogenerate tclIndex pkgIndex
install-data-hook:
	$(TCLSH) mkIndex.tcl $(drawstructdir)
