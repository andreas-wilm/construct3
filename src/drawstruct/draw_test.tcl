#!/bin/sh
#\
exec cs_wish "$0"  "$@"


##############################################################################
#
# draw_test.tcl - a simple exmaple how to use drawstruct
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
#  CVS $Id: draw_test.tcl,v 1.5 2004-05-25 13:28:45 wilm Exp $
#

package require DrawStruct 0.9

puts "THIS IS JUST A TEST OF DRAWSTRUCT:"
puts " displays a dummy structure if invoked without arguments"
puts " displays a structure from a ct-file if invoked with a ct-file as argument"

if {$argc} {
	set fname [lindex $argv 0]
	if {[file exists $fname]} {
		set result [LoadConnect $fname]
		set projname [lindex $result 0]
		set seq [lindex $result 1]
		set l_basepairs [lindex $result 2]
		set l_prob {}
		for {set i 0} {$i<[llength $l_basepairs]} {incr i} {
			if {[lindex $l_basepairs $i]>0} {
				lappend l_prob 1.0
			} else {
				lappend l_prob 0.0
			}
		}
	}
} else {
	set projname    "xb6711_768"
	set seq         "----guG-gCAucn--AUGAcGAcCuguuuu-nAAAccccnUCc---ggGg---nGGcuCugG-UCUGAugnnn-------------CGcCac-------"
	set l_basepairs [list 0 0 0 0 93 92 91 0 90 89 0 0 0 0 0 0 0 0 0 0 67 66 65 63 62 61 60 59 58 57 56 0 0 0 0 0 51 50 49 48 0 0 0 0 0 0 0 40 39 38 37 0 0 0 0 31 30 29 28 27 26 25 24 0 23 22 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 9 7 6 5 0 0 0 0 0 0 0]
	set l_prob      [list -1.000000 -1.000000 -1.000000 -1.000000 0.563499 0.577794 0.539872 -1.000000 0.414616 0.487174 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 0.306079 0.230737 0.677644 0.173839 0.128419 0.672792 0.649692 0.438619 0.572589 0.714232 0.117766 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 0.134727 0.598272 0.512213 0.069386 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 0.069386 0.512213 0.598272 0.134727 -1.000000 -1.000000 -1.000000 -1.000000 0.117766 0.714232 0.572589 0.438619 0.649692 0.672792 0.128419 0.173839 -1.000000 0.677644 0.230737 0.306079 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 0.487174 0.414616 0.539872 0.577794 0.563499 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000]
}

set header      [list {Optimal Consensus Structure} {Project: foo/bar} {May 14, 2004  11:36}]
set optimal     1
set win_pos    "+115+115"

DrawStruct::DrawStructure $seq $l_basepairs $optimal $win_pos \
	         -prob $l_prob -projname $projname -header $header

