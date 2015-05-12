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


# CVS $Id: Makefile.am,v 1.35 2007-10-22 10:43:23 steger Exp $

# is added to auto_path
# construct libdir
CSLIBDIR=\"$(libdir)/@PACKAGE@\"
# drawstruct libdir (must confirm to drawstruct Makefile.am)
DSLIBDIR=\"$(libdir)/@PACKAGE@/drawstruct\"

CC=@TCL_CC@
INCLUDES = -I$(includedir) @TCL_INCLUDE_SPEC@ @TK_INCLUDE_SPEC@ @TK_XINCLUDES@
AM_CFLAGS   = -ansi -Wall @CFLAGS@ @MACOS_CFLAG@  @TCL_CFLAGS_OPTIMIZE@ -DCSLIBDIR=$(CSLIBDIR) -DDSLIBDIR=$(DSLIBDIR)
# @TK_DEFS@ @TCL_DEFS@
# TCL_CFLAGS_WARNING:@TCL_CFLAGS_WARNING@
LIBS        = @TCL_LIB_SPEC@ @TCL_LIBS@ @TK_LIB_SPEC@ @TK_LIBS@ @LIBS@
AM_LDFLAGS  = @TCL_LD_FLAGS@ @TCL_LDFLAGS_OPTIMIZE@  @LDFLAGS@

bin_PROGRAMS    =   cs_wish
cs_wish_SOURCES =   bp_prob_mat.c   bp_prob_mat.h   \
                    bmatch.c        bmatch.h        \
					callbacks.c     callbacks.h     \
					consensus.c     consensus.h     \
					cs_bmatch.c     cs_bmatch.h	    \
                    cs_dpm.c        cs_dpm.h        \
                    cs_imatch.c     cs_imatch.h     \
					cs_seqio.c      cs_seqio.h      \
					dlist.c         dlist.h         \
                    generic_utils.c generic_utils.h \
                    gui.c           gui.h           \
					if.c            if.h            \
                    imatch.c        imatch.h        \
					main.c          main.h          \
					mic.c           mic.h           \
					move_update.c   move_update.h   \
                    mwmatch.c       mwmatch.h       \
                    mx.c            mx.h            \
					opt_struct.c    opt_struct.h    \
					proj.c          proj.h          \
                    public_datatypes.h \
					seqio.c         seqio.h         \
                    seq_stat.c      seq_stat.h      \
					shift.c         shift.h         \
					stack.c         stack.h         \
                    struct.c        struct.h        \
					subopt_struct.c subopt_struct.h \
                    utils.c         utils.h         \
					stopwatch.c     stopwatch.h

