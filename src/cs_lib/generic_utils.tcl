##############################################################################
# 
# generic_utils.tcl - generic routines
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
###############################################################################

#
#  CVS $Id: generic_utils.tcl,v 1.7 2004-08-11 08:14:22 wilm Exp $    
#



###############################################################################
###############################################################################
#
#          TCL
#
###############################################################################
###############################################################################



###   Debugger   
#
# Inspired by
# http://mini.net/tcl/1151
# and
# http://mini.net/tcl/473
#
# Invoke for online debugging
# s is an arbitrary identifier
#
# if global list DEBUGGER_SKIP includes this identifier debugging will be skipped 
#
# if global bool DEBUGGER_SKIP_ALL is 1, debugging will be completely skipped
#
proc Debugger {{s {}}} {
#############

    if {[info exists ::DEBUGGER_SKIP_ALL]} {
        if {$::DEBUGGER_SKIP_ALL==1} {
            return
        }
    }

    if {![info exists ::DEBUGGER_SKIP]} {
        set ::DEBUGGER_SKIP {}
    } elseif {[lsearch -exact $::DEBUGGER_SKIP $s]>=0} {
        return
    }


    if {[catch {lindex [info level -1] 0} caller]} {
        set caller ::
    }
    
    while 1 {
    	puts -nonewline "DEBUG(func=\"$caller\"  tag=\"$s\")\n> "
        flush stdout

    	gets stdin line

    	if {$line=="c"} {
            puts "continuing..";
            break
        }
        
        # i ?<pattern>?: info local ?<pattern>?
        # g ?<pattern>?: info globals ?<pattern>?
        #
        set argc [llength $line]
        if {$argc>2} {
            set flag "h"
        } elseif {$argc==1} {
            set flag  [lindex $line 0]
            set arg   "*"
        } else {
            set flag  [lindex $line 0]
            set arg   [lindex $line 1]
        }

        set print 0
        set local 0
        if {$flag=="l"} {
            set local 1
        } elseif {$flag=="lp"} {
            set local 1
            set print 1          
        } elseif {$flag=="gp"} {
            set print 1
        } elseif {$flag!="g"} {
            puts "Valid debugging commands are:"
            puts " c              : continue"
            puts " l ?<pattern>?  : info locals  ?<pattern>?"
            puts " g ?<pattern>?  : info globals ?<pattern>?"
            puts " lp ?<pattern>? : print local vars ?<pattern>?"
            puts " gp ?<pattern>? : print global vars ?<pattern>?"
            puts " h              : this help"        
            continue
        }
        
        if {$local} {
            set info "info locals $arg"
            set level 1
        } else {
            set info "info globals $arg"
            set level #0
        }

    	if {[catch {uplevel $info} res]} {
    		puts "   no match found"
            continue
        }
        
        # show only matches ?
        if {!$print} {
            puts "   $res"
            continue
        }
        
        # show all vars foreach match
        foreach v [split $res] {
            # array ?
            if {[uplevel $level  array exists $v]==1} {
                foreach an [uplevel $level array names $v] {
                    catch {uplevel $level  set ${v}(${an})} res
                    puts "   ${v}(${an}) = $res"
                }
            # scalar?
            } else {
                catch {uplevel $level set $v} res
                puts "   $v = $res"
            }
        }
        
    };# while
}
# Debugger




###############################################################################
###############################################################################
#
#          TK
#
###############################################################################
###############################################################################



###   SetCursor   
#
# busy|normal
#
proc SetCursor {state} {
########################

    if {$state=="busy"} {
        foreach win [winfo children .] {
            $win config -cursor {watch red white}
        }
    } elseif {$state=="normal"} {
        foreach win [winfo children .] {
            $win config -cursor {}
        }
    } else {
        puts -nonewline "ERROR: Call SetCursor with"
        puts " \"busy\" or \"normal\". (and not \"$state\")"
    }
}
# SetCursor
