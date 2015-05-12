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

# CVS $Id: Makefile.am,v 1.28 2007-10-26 13:51:13 steger Exp $


bin_SCRIPTS = cs_fold cs_dp cs_remgaponly_cols \
              cs_proj_conv cs_shift cs_struct_displ \
              csdpbatch csfoldbatch
CS_TCL_FILES =	bi_match.tcl \
                balloon_help.tcl \
		        callbacks.tcl \
		        circles.tcl \
		        cs_proj.tcl \
                generic_utils.tcl \
		        gui.tcl \
		        initcore.tcl \
		        logogui_ns.tcl \
                minialn.tcl \
                p.tcl \
                pattern.tcl \
                rnamlout.tcl \
                struct.tcl \
                struct_aln.tcl \
		        sub_gui.tcl \
                upe.tcl \
		        utils.tcl \
		        GetOpts.tcl \
		        mkIndex.tcl
 
GPL_FILE = COPYING
EXTRA_DIST = $(bin_SCRIPTS) $(CS_TCL_FILES) 

# must be the same as in cs_wish Makefile.am for auto_index
constructdir = $(libdir)/@PACKAGE@/
construct_DATA = $(CS_TCL_FILES) $(GPL_FILE)

# autogenerate tclIndex pkgIndex
install-data-hook:
	$(TCLSH) mkIndex.tcl  $(DESTDIR)$(constructdir)
