##############################################################################
#
# cs_proj.tcl - procedures for project-IO
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
#  CVS $Id: cs_proj.tcl,v 1.23 2007-10-22 10:43:23 steger Exp $
#



namespace eval cs_proj {

###
#
# 1.line: header:"/// ConStruct Project-File Version <X.Y>"
# 2.line: project-name:"Project-Name: <name>"
# 3.line: alignment-file:Alignment: <file>
#
# no empty lines or comments in first three lines allowed !
# comments:" ^\ *#"
#
# entry:    begin entry
#           id: <id>
#           weight: <weight>
#           seqlen: <unaligned sequence length>
#           bpmat: <bpmat>
#           foldcmd: <foldcmd>
#           comment: <stored comment>
#           mapinfo: <mapping info (UNaligned indices!)>
#           end entry
###


# FIXME:  implement version check

# frontends :
# reading: proc read {fname}
#   returns "ERROR" on failure
#   otherwise
#   "OK"
#
# writing: proc write {fname}

# FEATURES:
#
# stores id, weight, bpmat, foldcmd and comment for each entry
# project files can be links
# all files can be relativ
#
# consider to set EXPORT_2_CSWISH to 0
#

###########################################################################

    variable proj    ;# (file|name|version|aln_path)
    variable nseq    ;# number of sequence entries
    variable seq_id  ;# sequence ids, 1-nseq
    variable weight  ;# sequence weights, 1-nseq
    variable seqlen  ;# unaligned sequence length
    variable bpmat   ;# basepair-prob-matrix file, 1-nseq
    variable foldcmd ;# rnafold cmdline foldcmderature, 1-nseq
    variable comment ;# entry comment, 1-nseq
    variable mapinfo;# mapping info, 1-nseq

    variable DEBUG 0
    variable VERSION 3.0
    variable EXPORT_2_CSWISH 1;# set to 0 if you want to use this without cs_wish
                                 ;# exports all values after read to cs_wish
    variable IMPORT_FROM_CSWISH 1;# set to 0 if you want to use this without cs_wish
                                 ;# imports all values from cs_wish before writing
    variable FILE_EXT_LIBZ   [list ".Z" ".gz"]


    #######################################################################
    ##########
    #                          MISC
    ##########
    #######################################################################



    ###   namespace mapinfo
    #
    #
    #
    namespace eval mapinfo {

    	###   parse
    	#
    	# parses and checks an mapinfo entry from project file
    	#
    	# IN:
    	#  str: mapinfo string read from project file
    	#  seqlen: optional known sequence length, enables extra checks
    	#
    	# OUT:
    	#   returns a expanded sorted list of mapinfo tokens
    	#   or an empty list on empty input or error
    	#   expaned means, ranges are transformed to single nt info tokens
    	#
    	# valid syntax:
    	#   mapinfo: <ntidx_or_range>:<ntidx_or_state> ...
    	#     where
    	#     ntidx_or_range: <int>|<int>-<int>
    	#     ntidx_or_state: p|u|<int>
    	#        where
    	#         u means nt <ntidx> is unpaired
    	#         p means nt <ntidx> is paired
    	#         <int> is the index of the bp-partner
    	#   the use of a range of course forbids the simultaneous use of ntidx
    	#
    	# example:
    	# 3:62 4:p 5:u 6-10:u 11-20:p
    	#
    	proc parse {str {seqlen -1}} {
    		########################

    		if {$str==""} {
    			return ""
    		}

    		# clean seqlen
    		if {$seqlen!=-1} {
    			if {[catch {set seqlen [expr {abs([format "%d" $seqlen])}]}]} {
    				set seqlen -1
    			}
    		}
    		# now seqlen is in any case -1 or integer > 0

    		###   precprocessing
    		#
    		# remove multiple internal spaces
    		regsub -all -- {  *} $str " " str
    		# analyze each token
    		set tokenlist {}
    		foreach token [split [string tolower $str]] {
    			if {[regexp -- {^[0-9]+[0-9\-]*:p|u|[0-9]+$} $token]} {
    				lappend tokenlist $token
    			} else {
    				p::warn "ignoring invalid mapinfo token \"$token\""
    				continue
    			}
    		}
    		if {[llength $tokenlist]==0} {
    			return {}
    		}


    		###  in depth check
    		#
    		set parsedtokenlist {}
    		foreach token $tokenlist {
    			p::debug "analyzing token $token"

    			# explicit basepair
    			#
    			if {[scan $token "%d:%d" ntidx1 ntidx2]==2} {
    				if {$ntidx1<1 || ($seqlen!=-1 && $ntidx1>$seqlen)} {
    					p::warn "ignoring mapinfo token \"$token\" (invalid index $ntidx1)"
    					continue
    				}
    				if {$ntidx2<1 || ($seqlen!=-1 && $ntidx2>$seqlen)} {
    					p::warn "ignoring mapinfo token \"$token\" (invalid index $ntidx2)"
    					continue
    				}

    				lappend parsedtokenlist "$ntidx1:$ntidx2"
    				lappend parsedtokenlist "$ntidx2:$ntidx1"
    				continue
    			}

    			# one nt
    			#
    			if {[scan $token "%d:%s" ntidx state]==2} {
    				if {$ntidx<1 || ($seqlen!=-1 &&  $ntidx>$seqlen)} {
    					p::warn "ignoring mapinfo token \"$token\" (invalid index $ntidx)"
    					continue
    				}
    				if {$state!="p" && $state!="u"} {
    					p::warn "ignoring mapinfo token \"$token\" (invalid state $state)"
    					continue
    				}

    				lappend parsedtokenlist $token
    				continue
    			}

    			# range
    			#
    			set tmp [split $token ":"]
    			set range [lindex $tmp 0]
    			set state [lindex $tmp 1]
    			if {$state!="p" && $state!="u"} {
    				p::warn "ignoring mapinfo token \"$token\" (invalid state $state)"
    				continue
    			}
    			# test for valid range
    			if {[scan $range "%d-%d" start end]!=2} {
    				p::warn "ignoring mapinfo token \"$token\" (invalid range $state)"
    				continue
    			}
    			if {$start>$end} {
    				p::warn "swapping start and end for mapinfo token $token"
    				set tmp $start
    				set start $end
    				set end $tmp
    			}
    			if {$start<1 || ($seqlen!=-1 && $end>$seqlen)} {
    				p::warn "ignoring mapinfo token \"$token\" (invalid range $range)"
    				continue
    			}

    			for {set i $start} {$i<=$end} {incr i} {
    				lappend parsedtokenlist "$i:$state"
    			}
    		}
    		return [lsort -dictionary -unique $parsedtokenlist]
    	}
    	# parse



    	###   precheck
    	#
    	# frontend to is_ functions, internal use only
    	#
    	# BEWARE: UNALIGNED INDICES
    	#
    	proc precheck {seqno ntidx} {
    		#######################

    		# mapinfo entry not found?
    		if {! [exists $seqno]} {
    			return -1
    		}
    		# request ntidx not found?
    		set info [lsearch -inline -glob $cs_proj::mapinfo($seqno) "${ntidx}:*"]
    		if {[scan $info "%d:%s" alreadyknown mapinforules] != 2} {
    			return -1
    		} else {
    			return $mapinforules
    		}
    	}
    	# precheck



    	###   exists
    	#
    	# return 1 if mapping info for a particular seqexists,0 otherwise
    	#
    	proc exists {seqno} {
    		if {[info exists cs_proj::mapinfo($seqno)]} {
    			return 1
    		} else {
    			return 0
    		}
    	}
    	# exists


    	###   is_paired
    	#
    	# OUT:
    	#  ntidx if forced set to pair with ntidx
    	#  1 if forced set to pair (p)
    	#  0 if not paired (u)
    	#  -1 if no info found
    	#
    	# BEWARE: UNALIGNED INDICES
    	#
    	proc is_paired {seqno ntidx} {
    		########################

    		set mapinforules [precheck $seqno $ntidx]
    		if {$mapinforules==-1} {
    			return -1
    		}

    		if {$mapinforules=="p" || [string is integer $mapinforules]} {
    			return 1
    		} elseif {$mapinforules=="u"} {
    			return 0
    		} else {
    			p:error "internal error: unhandled mapinforules $mapinforules"
    			return -1
    		}
    	}
    	# is_paired



    	###   is_paired_with
    	# OUT:
    	#  ntidx if forced set to pair with ntidx
    	#  0 otherwise
    	#
    	proc is_paired_with {seqno ntidx} {
    		############################

    		set mapinforules [precheck $seqno $ntidx]
    		if {$mapinforules==-1} {
    			return 0
    		}
    		if {[string is integer $mapinforules]} {
    			return $mapinforules
    		} else {
    			return 0
    		}
    	}
    	# is_paired_with



    	###   is_not_paired
    	#
    	# ~ opposite of is_paired
    	#
    	#  1 if not paired (u)
    	#  0 if forced set to pair (p or ntidx)
    	#  -1 if no info found
    	#
    	# BEWARE: UNALIGNED INDICES
    	#
    	proc is_not_paired {seqno ntidx} {
    		############################

    		set mapinforules [precheck $seqno $ntidx]
    		if {$mapinforules==-1} {
    			return -1
    		}

    		if {$mapinforules=="p" || [string is integer $mapinforules]} {
    			return 0
    		} elseif {$mapinforules=="u"} {
    			return 1
    		} else {
    			p:error "internal error: unhandled mapinforules $mapinforules"
    			return -1
    		}
    	}
    	# is_not_unpaired

    }
    # namespace eval mapinfo



    ###   print_out
    #
    # print out all namespace variables
    # debugging only
    #
    #
    proc print_out {} {
    	#############
    	variable proj
    	variable nseq
    	variable seq_id
    	variable weight
    	variable seqlen
    	variable bpmat
    	variable foldcmd
    	variable comment
    	variable mapinfo

        puts "##########"
        puts "[namespace current] variables:\n"
        if {[info exists proj(file)]} {
            puts "[namespace current]:proj(file)    = $proj(file)"
        }
        puts "[namespace current]:proj(name)     = $proj(name)"
        puts "[namespace current]:proj(version)  = $proj(version)"
        puts "[namespace current]:proj(aln_path) = $proj(aln_path)"

        puts "[namespace current]:nseq=$nseq"

        for {set n 1} {$n<=$nseq} {incr n} {
            puts "[namespace current]:seq_id($n)   = $seq_id($n)"
            puts "[namespace current]:weight($n)   = $weight($n)"
    		puts "[namespace current]:seqlen($n)   = $seqlen($n)"
            puts "[namespace current]:bpmat($n)    = $bpmat($n)"
            puts "[namespace current]:foldcmd($n)  = $foldcmd($n)"
            puts "[namespace current]:comment($n)  = $comment($n)"
            puts "[namespace current]:mapinfo($n)  = $mapinfo($n)"
        }
        puts "##########"
    }
    # print_out



    ###   unset_all
    #
    # unsets all namespace variables
    #
    proc unset_all {} {
    	#############
    	variable proj
    	variable nseq
    	variable seq_id
    	variable weight
    	variable seqlen
    	variable bpmat
    	variable foldcmd
    	variable comment
    	variable mapinfo

    	catch {unset proj}
    	catch {unset nseq}
    	catch {unset seq_id}
    	catch {unset weight}
    	catch {unset seqlen}
    	catch {unset bpmat}
    	catch {unset foldcmd}
    	catch {unset comment}
    	catch {unset mapinfo}
    }
    # unset_all



    ###########################################################################
    ##########
    #                          READ
    ##########
    ###########################################################################



    ###   rel_fileentry_2abs
    #
    # convert a relative pathname of an entry to absolut one
    # check also for project link !
    #
    proc rel_fileentry_2abs {f_proj f_entry} {
    	####################################

        while {![catch {file readlink $f_proj} result]} {
            set f_proj $result
        }
        p::debug "absolut path of proj is $f_proj"

        set old_dir [pwd]
        cd [file dirname $f_proj]
        cd [file dirname $f_entry]


        set abspath [file join [pwd] [file tail $f_entry]]
        cd $old_dir

        p::debug "absolut path of entry is $abspath"

        return $abspath
    }
    # rel_fileentry_2abs




    ###   is_comment
    #
    #
    proc is_comment {str} {
    	#################

    	if {[regexp -- {^ *\#} $str dummy]} {
            return 1
        } else {
            return 0
        }
    }
    # is_comment




    ###   read_seq_entry
    #
    # read one sequence entry
    #
    # return a list of 4 elements (on success):
    #   id:        sequence id
    #   weight:    weight
    #   bpmat:   path to bpprobmat file
    #   foldcmd: fold command
    #
    # returns empty list on failure
    #
    proc read_seq_entry {fid} {
    	#####################

    	# setup optional paramaters
    	set comment  ""
    	set mapinfo ""
    	set seqlen   "-1"

        while {[gets $fid line]>=0} {
            set line [string trim $line]

            # skip comment lines
            if {[is_comment $line]} {
                continue
            }

            # detect end of entry
            if {[regexp -- {end entry} $line dummy]} {
                if { ( ! [info exists id])     ||
                     ( ! [info exists weight]) ||
                     ( ! [info exists bpmat])  ||
                     ( ! [info exists foldcmd])     } {
                    return {}
                }
                return [list $id $weight $seqlen $bpmat $foldcmd $comment $mapinfo]
            }

            regexp -- {id: *(.*)}      $line dummy id
            regexp -- {weight: *(.*)}  $line dummy weight
    		regexp -- {seqlen: *(.*)}  $line dummy seqlen
            regexp -- {bpmat: *(.*)}   $line dummy bpmat
            regexp -- {foldcmd: *(.*)} $line dummy foldcmd
            regexp -- {comment: *(.*)} $line dummy comment
    		regexp -- {mapinfo: *(.*)} $line dummy mapinfo
        }

        p::debug "entry end not found"

        return {}
    }
    # read_seq_entry




    ###   get_projname
    #
    #
    proc get_projname {fid} {
    	###################
         variable proj

        if {[gets $fid line]<0} {
            p::error "corrupt project file"
            return "ERROR"
        }
        set line [string trim $line]
        if { ! [regexp -- {Project-Name: (.*)} $line dummy projname]} {
            p::error "can't find project-name"
            return "ERROR"
        }
        if { ! [info exists projname]} {
            p::error "couldn't find any project-name"
            return "ERROR"
        }
        set proj(name) "$projname"

        p::debug "projname=$projname"

        return "OK"
    }
    # get_projname



    ###   get_aln_path
    #
    #
    proc get_aln_path {fid} {
    	###################
         variable proj

        if {[gets $fid line]<0} {
            p::error "corrupt project file"
            return "ERROR"
        }
        set line [string trim $line]
        if { ! [regexp -- {Alignment: (.*)} $line dummy alnpath]} {
            p::error "can't find alignment path"
            return "ERROR"
        }
        if { ! [info exists alnpath]} {
            p::error "couldn't find any alignment path"
            return "ERROR"
        }
        set proj(aln_path) "$alnpath"

        p::debug "alnpath=$alnpath"

        return "OK"
    }
    # get_aln_path





    ###   check_header
    #
    # checks header of open project file with id <fid>
    # and sets proj(version)
    #
    proc check_header {fid} {
    	###################
    	variable proj

        if {[gets $fid line]<0} {
            return "ERROR"
        }
        set line [string trim $line]

        if { ! [regexp -- {/// ConStruct Project-File Version ([0-9\.]*)} \
                                                    $line dummy version]}  {
            return "ERROR"
        }
        if { ! [info exists version]} {
            return "ERROR"
        }
        set proj(version) $version

        p::debug "version=$version"

        return "OK"
    }
    # check_header





    ###   read
    #
    # reads a project file
    # returns OK on success, error message otherwise
    #
    proc read {fname {check_for_matrices 1}} {
    	####################################
        variable proj
        variable nseq
        variable seq_id
        variable weight
        variable bpmat
        variable foldcmd
    	variable comment
    	variable mapinfo
    	variable seqlen
        variable EXPORT_2_CSWISH
        variable FILE_EXT_LIBZ


        # reentrance
        unset_all

        # check existence
        #
        if { ! [file exists $fname]} {
    		set errmsg "project file \"$fname\" doesn't exist"
            p::error $errmsg
            return "$errmsg"
        }
        set fid [open "$fname" r]

        set proj(file) "$fname"

        # check header
        #
        if {[check_header $fid]=="ERROR"} {
    		set errmsg "invalid header in project file \"$fname\""
            p::error $errmsg
            close $fid
            unset_all
            return "$errmsg"
        }

        # get project name
        #
        if {[get_projname $fid]=="ERROR"} {
    		set errmsg "couldn't get project-name from project file \"$fname\""
            p::error $errmsg
            close $fid
            unset_all
            return "$errmsg"
        }

        # get alignment path
        #
        if {[get_aln_path $fid]=="ERROR"} {
    		set errmsg "couldn't get alignment path from project file \"$fname\""
            p::error §errmsg
            close $fid
            unset_all
            return "$errmsg"
        }

        # check alignment
        #
        if {[file pathtype $proj(aln_path)]=="relative"} {
            set proj(aln_path) [rel_fileentry_2abs $proj(file) $proj(aln_path)]
        }
        if { ! [file exists $proj(aln_path)]} {
    		set errmsg "alignment file $proj(aln_path) doesn't exist"
            p::error $errmsg
            close $fid
            unset_all
            return "$errmsg"
        }


        p::debug "starting to parse entries"

        # read all entries
        #
        set nof_seqs 0
        while {[gets $fid line]>=0} {

            set line [string trim $line]
            # skip comments
            if {[is_comment $line]} {
                continue
            }

            if {[regexp -- {begin entry} $line dummy]} {

                incr nof_seqs

                p::debug "parsing entry $nof_seqs"

                set seq_entry [read_seq_entry $fid]

                if {[llength $seq_entry]!=7} {
    				set errmsg "incomplete entry (no $nof_seqs) in project file \"$fname\""
                    p::error "$errmsg"
                    close $fid
                    unset_all
                    return $errmsg
                }

                # store entry
                set seq_id($nof_seqs)  [lindex $seq_entry 0]
                set weight($nof_seqs)  [expr {double([lindex $seq_entry 1])}]
    			set seqlen($nof_seqs)  [lindex $seq_entry 2]
                set bpmat($nof_seqs)   [lindex $seq_entry 3]
                set foldcmd($nof_seqs) [lindex $seq_entry 4]
                set comment($nof_seqs) [string trim [lindex $seq_entry 5]]
    			set mapinfo($nof_seqs) [mapinfo::parse [lindex $seq_entry 6] $seqlen($nof_seqs)]

                # check bpmat
                #
                if {[file pathtype $bpmat($nof_seqs)]=="relative"} {
                    set bpmat($nof_seqs) [rel_fileentry_2abs $proj(file) $bpmat($nof_seqs)]
                }
                if {$check_for_matrices} {
                    if { ! [file exists $bpmat($nof_seqs)]} {
                        set found 0
                        foreach ext $FILE_EXT_LIBZ {
                            set tmpfile    "$bpmat($nof_seqs)"
                            append tmpfile "$ext"
                            if {[file exists "$tmpfile"]} {
                                set $bpmat($nof_seqs) "$tmpfile"
                                set found 1
                            }
                        }
                        if {! $found} {
                            p::error "basepairmatrix file $bpmat($nof_seqs) doesn't exist"
                            close $fid
                            unset_all
                            return ERROR
                        }
                    }
                }
                if {! [LibZ_supported]} {
                    foreach ext $FILE_EXT_LIBZ {
                        if {[file extension $bpmat($nof_seqs)]==$ext} {
                            set msg    "basepair matrix [file tail $bpmat($nof_seqs)] seems to be compressed"
                            append msg ", but libZ support was not compiled in"
                            p::warn $msg
                        }
                    }
                }
            }
        }
        set nseq $nof_seqs

        # Success: now export to C
        if {$EXPORT_2_CSWISH} {
            Proj_Exchange  $proj(file) $proj(name) $proj(version) \
                           $proj(aln_path) $nseq                  \
                           seq_id  bpmat weight foldcmd
        }

        close $fid

        Debugger "cs_proj"

        return "OK"
    }
    # read




    ###########################################################################
    ##########
    #                          WRITE
    ##########
    ###########################################################################


    ###   write_example_entry
    #
    #
    proc write_example_entry {fid} {
    	##########################

        puts $fid "\n"
    	puts $fid "#comment lines start with a dash and are ignored"
    	puts $fid "#comments inside sequence entries have special tags (see below) and are stored"
    	puts $fid "#example entry:"
        puts $fid "#begin entry"
        puts $fid "#   id:       <string> e.g. h_SelD"
        puts $fid "#   weight:   <int/double> e.g. 0.125"
    	puts $fid "#   seqlen:   <int> e.g. 67"
        puts $fid "#   bpmat:    <file> e.g. h_SelD.dat\[.Z\] or h_SelD_dp.ps\[.Z|.gz\]"
        puts $fid "#   foldcmd:  <string> e.g. cs_rnafold -T 37 -p -d 3"
    	puts $fid "#   comment:  <string> e.g. this is a comment"
    	puts $fid "#   mapinfo: <string> e.g. 3-5:p 8-11:u 12:24"
        puts $fid "#end entry"
        puts $fid "\n"
    }
    # write_example_entry



    ###   write
    #
    #
    proc write {fname} {
    	##############
    	variable proj
    	variable nseq
    	variable seq_id
    	variable weight
    	variable bpmat
    	variable foldcmd
    	variable comment
    	variable mapinfo
    	variable seqlen
    	variable IMPORT_FROM_CSWISH

        set fid [open "$fname" w]

        puts $fid "/// ConStruct Project-File Version $proj(version)"
        puts $fid "Project-Name: $proj(name)"
        puts $fid "Alignment: $proj(aln_path)"

        write_example_entry $fid

    	# setup/check optional parameters
    	for {set n 1} {$n<=$nseq} {incr n} {
    		if { ! [info exists comment($n)]} {
    			set comment($n) ""
    		}
    		if { ! [info exists mapinfo($n)]} {
    			set mapinfo($n) ""
    		}
    		# old version don't have this entry
    		if { ! [info exists seqlen($n)]} {
    			set seqlen($n) "-1"
    		}
    	}

        # write each entry
        for {set n 1} {$n<=$nseq} {incr n} {
            puts $fid "begin entry"
            puts $fid "    id:      $seq_id($n)"
            puts $fid "    weight:  $weight($n)"
    		   puts $fid "    seqlen:  $seqlen($n)"
            puts $fid "    bpmat:   $bpmat($n)"
            puts $fid "    foldcmd: $foldcmd($n)"
    		   puts $fid "    comment: $comment($n)"
    		   puts $fid "    mapinfo: $mapinfo($n)"
            puts $fid "end entry\n"
        }
        close $fid
    }
    # write

}
### namespace eval cs_proj
