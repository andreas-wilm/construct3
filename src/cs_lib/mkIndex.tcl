#!/bin/sh
#\
exec tclsh "$0" "$@"
######################################################################


##############################################################################
#
# mkIndex.tcl - create tclIndex and pkgIndex
#
# Copyright (C) 2001-2004 Institute of Physical Biology,
#                    Heinrich-Heine Universität Düsseldorf, Germany
#                    <construct@biophys.uni-duesseldorf.de>
#
# This file is part of ConStruct.
#
# ConStruct is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# ConStruct is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ConStruct; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
##############################################################################

#
# CVS $Id: mkIndex.tcl,v 1.12 2004-12-08 13:33:17 wilm Exp $
#


set dir [lindex $argv 0]
if { ! [file exists $dir] || ! [file isdirectory $dir]} {
	puts stderr "Missing arg: directory for index generation"
    return
}
puts "Generating tclIndex and pkgIndex.tcl in $dir"

auto_mkindex $dir *tcl
#pkg_mkIndex -verbose $rootdir *tcl *[info sharedlibextension]
pkg_mkIndex $dir *tcl *[info sharedlibextension]

