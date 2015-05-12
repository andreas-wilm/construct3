# Copyright (C) 2004 Andreas Wilm <wilm@biophys.uni-duesseldorf.de>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#
# CVS $Id: p.tcl,v 1.3 2005-10-20 11:43:57 wilm Exp $
#


# -- p
#
#
# namespace for putting messages (verbose,debug,warn,error,fixme)
#


package provide p 0.1


###   namespace p
#
#
namespace eval ::p {
##################
    variable debug 0
    variable verbose 0

    ###   be_verbose
    #
    proc be_verbose {bool} {
    	variable verbose
    	set verbose $bool
    }
    # be_verbose


    ###  enable_debug
    #
    proc enable_debug {bool} {
    	variable debug
    	set debug $bool
    }
    # enable_debug


    ###   verbose
    #
    # be verbose
    #
    proc verbose {msg} {
    ####################
        variable verbose

        # if {!$verbose && !$debug}
        if {!$verbose} {
    		return
        }
        puts "$msg"
    }
    # verbose


    ###   debug
    #
    # print debug msg
    #
    proc debug {msg} {
    ##################
        variable debug

        if {!$debug} {
            return
        }
        if {[catch {lindex [info level -1] 0} caller]} {
            set caller ::
        }
        puts stdout "DEBUG($caller): $msg"
        flush stdout
    }
    # debug


    ###   error
    #
    # print error msg
    # optionally: exit with return code rc if non null
    #
    proc error {msg {rc 0}} {
    ##################

        if {[catch {lindex [info level -1] 0} caller]} {
            set caller ::
        }
        puts stderr "ERROR($caller): $msg"
        flush stderr
        if {$rc!=0} {
	        exit $rc
		}
    }
    # error


    ###   warn
    #
    # print warn msg
    #
    proc warn {msg} {
    #################

        if {[catch {lindex [info level -1] 0} caller]} {
            set caller ::
        }
        puts stdout "WARNING($caller): $msg"
        flush stdout
    }
    # warn


    ###   fixme
    #
    # print fixme msg
    #
    proc fixme {msg} {
    #################

        if {[catch {lindex [info level -1] 0} caller]} {
            set caller ::
        }
        puts stdout "PLEASE FIX ME($caller): $msg"
        flush stdout
    }
    # fixme
}
# namespace p

