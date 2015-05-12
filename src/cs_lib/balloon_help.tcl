##############################################################################
#
# balloon.tcl - procedures used by balloon help
#
# Copyright (C) 1996-1997 Stewart Allen
#
# This is part of vtcl source code Adapted for
# general purpose by Daniel Roche <dan@lectra.com>
#
# Support for canvas tags added by
# Andreas Wilm <wilm@biophys.uni-duesseldorf.de>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
##############################################################################
#
# Modified for use in ConStruct 3
#
###############################################################################

#
#  CVS $Id: balloon_help.tcl,v 1.8 2004-08-11 08:14:21 wilm Exp $
#

package provide BalloonHelp 1.0

# FIXME: implement clean up after window destroy

namespace eval BalloonHelp {

    variable Bulle

    ###   set_balloon
    #
    #
    proc set_balloon {wpath message {tag ""}} {
        variable Bulle

    	# tag
    	#
    	if {$tag!=""} {
    		set Bulle($tag) $message

    		$wpath bind $tag <Enter> {
    			set srctag ""
    			foreach t [%W gettags current] {
    				if {[info exists BalloonHelp::Bulle($t)]} {
    					set srctag $t
    					break
    				}
    			}
    			if {$srctag==""} {
    				p:error "upps...couldn't find srctag (window %W taglist [%W gettags current])"
    				break
    			}
    			set BalloonHelp::Bulle(set) 0
    			set BalloonHelp::Bulle(first) 1
    			set BalloonHelp::Bulle(id) \
    				[after 500 {BalloonHelp::createballoon %W $BalloonHelp::Bulle($srctag) %X %Y}]
    		}
    		$wpath bind $tag <Motion> {
    			set srctag ""
    			foreach t [%W gettags current] {
    				if {[info exists BalloonHelp::Bulle($t)]} {
    					set srctag $t
    					break
    				}
    			}
    			if {$srctag==""} {
    				p:error "upps...couldn't find srctag (window %W)"
    				break
    			}
    			if {$BalloonHelp::Bulle(set) == 0} {
    				after cancel $BalloonHelp::Bulle(id)
    				set BalloonHelp::Bulle(id) \
    					[after 500 {BalloonHelp::createballoon %W $BalloonHelp::Bulle($srctag) %X %Y}]
    			}
    		}
    		$wpath bind $tag <Leave> {
    			set BalloonHelp::Bulle(first) 0
    			BalloonHelp::kill_balloon
    		}
    		$wpath bind $tag <Button> {
    			set BalloonHelp::Bulle(first) 0
    			BalloonHelp::kill_balloon
    		}

    	# window
    	#
    	} else {
    		set Bulle($wpath) $message

    		bind $wpath <Enter> {
    			set BalloonHelp::Bulle(set) 0
    			set BalloonHelp::Bulle(first) 1
    			set BalloonHelp::Bulle(id) \
    				[after 500 {BalloonHelp::createballoon %W $BalloonHelp::Bulle(%W) %X %Y}]
    		}
    		bind $wpath <Motion> {
    			if {$BalloonHelp::Bulle(set) == 0} {
    				after cancel $BalloonHelp::Bulle(id)
    				set BalloonHelp::Bulle(id) \
    					[after 500 {BalloonHelp::createballoon %W $BalloonHelp::Bulle(%W) %X %Y}]
    			}
    		}
    		bind $wpath <Leave> {
    			set BalloonHelp::Bulle(first) 0
    			BalloonHelp::kill_balloon
    		}
    		bind $wpath <Button> {
    			set BalloonHelp::Bulle(first) 0
    			BalloonHelp::kill_balloon
    		}
    	}
    }


    ###   kill_balloon
    #
    #
    proc kill_balloon {} {
        variable Bulle

        after cancel $Bulle(id)
        if {[winfo exists .balloon] == 1} {
            destroy .balloon
        }
        set Bulle(set) 0
    }



    ###   createballoon
    #
    #
    proc createballoon {target message {cx 0} {cy 0} } {
        variable Bulle

        if {$Bulle(first) == 1 } {
            set Bulle(first) 2
    	    if { $cx == 0 && $cy == 0 } {
    	        set x [expr [winfo rootx $target] + ([winfo width $target]/2)]
    	        set y [expr [winfo rooty $target] + [winfo height $target] + 4]
    	    } else {
    	        set x [expr $cx + 4]
    	        set y [expr $cy + 4]
    	    }
    	 	if {[winfo exists .balloon]==1} {
            	destroy .balloon
        	}
    		toplevel .balloon -bg black
            wm overrideredirect .balloon 1
            label .balloon.l \
                -text $message -relief flat \
                -bg #ffffaa -fg black -padx 2 -pady 0 -anchor w
            pack .balloon.l -side left -padx 1 -pady 1
            wm geometry .balloon +${x}+${y}
            set Bulle(set) 1
        }
    }
}
# namespace eval BalloonHelp





###   example code   
#
#
# package require BalloonHelp <OR> source "./balloon.tcl"
#
# option add *highlightThickness 0
#
# label .head -text "move the pointer\nhere" -relief groove
# label .suite -text "another balloon here" -relief groove
#
# button .bye -text "Quit" -command exit
#
# pack .head -side top -fill x -pady 10 -padx 10
# pack .suite .bye -side right -pady 10 -padx 10
# pack .bye -side bottom -pady 10 -padx 10
#
# BalloonHelp::set_balloon .head "first balloon"
# BalloonHelp::set_balloon .suite "second balloon"
#
#
###
