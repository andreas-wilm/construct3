##############################################################################
#
# utils.tcl - misc. routines
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
#  CVS $Id: utils.tcl,v 1.51 2007-10-22 10:43:23 steger Exp $
#

###   DumpFrontend   
#
#
proc DumpFrontend {vars_or_tags} {
    ############################
    global FILE_TYPES
    global opts

    set init_f "cs3_${vars_or_tags}.dump.txt"


    set filename [tk_getSaveFile        \
        -initialfile  $init_f           \
        -initialdir   $opts(dir,work)   \
        -filetypes    $FILE_TYPES(txt)  \
        -title        "Select dump-text file"]

    update idletasks

    if {$filename == ""} {
        puts "Save dialog aborted..."
    	set fid stdout
    	set fid_must_be_closed 0
    } else {
    	set fid [open $filename w]
    	set fid_must_be_closed 1
    }

    if {$vars_or_tags=="vars"} {
        DumpVars $fid
    } elseif {$vars_or_tags=="tags"} {
        DumpTags $fid
    } else {
        p::error "invalid arg vars_or_tags=$vars_or_tags"
        return
    }
    if {$fid_must_be_closed} {
    	close $fid
    }
}
#  DumpFrontend



###   DumpVars   
#
# Prints out all variables given in l_global
# optionally uses already opened file fid
#
proc DumpVars {{fid stdout}} {
##############

    #set l_global [list \
    #              c COLOR FONT FILE_TYPES lb_dp_pos MENUINDEX \
    #              opts scale selection seq SEQHEIGHT w ]
    set     l_global [info globals]
    set nss [list cs_proj StructAln Upe Circles]

    # global vars
    foreach v $l_global {
        global $v
        if {![info exists $v]} {
            puts $fid "Upps...$v doesn't exist"
            continue
        }
        if {[array exists $v]} {
            # puts $fid "$v (array):"
            foreach a [array names $v] {
                upvar #0 $v dummy
                puts $fid "$v\($a\) = $dummy($a)"
           }
        } else {
            # puts $fid "$v (scalar): [set $v]"
            puts $fid "$v = [set $v]"
        }
    }

    # namespace vars
    foreach ns $nss {
    	foreach v [info vars ${ns}::*] {
    		p::fixme "access $v"
    	}
    }

    flush $fid
}
# DumpVars



###   DumpTags   
#
#
proc DumpTags {{fid stdout}} {
    ########################
    global c

    foreach id [$c(dp) find withtag all] {
        puts $fid "c(dp): id=$id tags=[$c(dp) itemcget $id -tags]"
    }
    foreach id [$c(al) find withtag all] {
        puts $fid "c(al): id=$id tags=[$c(al) itemcget $id -tags]"
    }

    flush $fid
}
# DumpTags




###   ReplaceTag   
#
# general routine for replacing tags
#
proc ReplaceTag {canvas old_tag new_tag} {
    ####################################
    $canvas addtag  $new_tag withtag $old_tag
    $canvas dtag    $new_tag  $old_tag
}
# ReplaceTag




###   DumpTextWindowContent   
#
# Save the text content of an arbitrary window to a given file
#
proc DumpTextWindowContent {win filename} {
    #####################################

    set fid [open $filename w]
    foreach {key value index} [${win}.text dump 1.0 end] {
        if {$key == "text"} {
            puts -nonewline $fid "$value"
        }
    }
    close $fid
}
# DumpTextWindowContent





###   PrintPS   
#
# w: window
# save_print: "save"|"print"|"printcolor|screen"
# args: list of options for postscript command
# filename may not contain path elements
#
proc PrintPS {w filename save_print {args ""}} {
    ##########################################
    global PS_WIDTH
    global opts
    global FILE_TYPES


    set    cmd "$w postscript"
    append cmd " -colormode [string tolower $opts(print_cmd,colormode)]"
    append cmd " -pagewidth [format "%-2.1fc" $PS_WIDTH]"
    append cmd " [join $args]"


    ###   check first
    #
    if {$save_print=="print" && $opts(print_cmd,grey) == ""} {
        p::error "Printer not configured!"
        return

    } elseif {$save_print=="printcolor" && $opts(print_cmd,color) == ""} {
        p::error "Color printer not configured!"
        return

    } elseif {$save_print=="screen" && $opts(print_cmd,screen) == ""} {
        p::error "\"print to screen\" command not configured!"
        return
    }


    p::debug "cmd=$cmd"
    if {$save_print=="print"} {
        set ps [eval $cmd]
        catch {eval exec $opts(print_cmd,grey) << {$ps}} result

    } elseif {$save_print=="printcolor"} {
        set ps [eval $cmd]
        catch {eval exec $opts(print_cmd,color) << {$ps}} result

    } elseif {$save_print=="screen"} {
        set ps [eval $cmd]
        catch {eval exec $opts(print_cmd,screen) << {$ps} &} result

    } else {

    	set filename [tk_getSaveFile      \
    	     -initialfile  $filename  \
             -initialdir   $opts(dir,work) \
             -filetypes    $FILE_TYPES(ps)  \
    	     -title        "Browse PS files"]

    	update idletasks

    	if {$filename == ""} {
            puts "Postscript save dialog aborted"
            return
        }

        append cmd " -file $filename"

    	eval $cmd
    }
}
# PrintPS






###   ColorProb   
#
# Return Hex-Color according to a prob
#
proc ColorProb {prob} {
    ###################

    if {$prob==0.0} {
    	return #000000000000
    }
    if {$prob<=0.5} {
    	return [format "#ffffffff%04x" [expr {int((0.5-$prob)*131070)}]]
    }
    if {$prob<=1.0} {
    	return [format "#ffff%04x0000" [expr {int((1.0-$prob)*131070)}]]
    } else {
    	return #ffff00000000
    }
}
# ColorProb




###   Validate_MIC   
#
#
# return "abort", if mic isn't valid
# and user denied recomputation
#
# return nothing if mic is valid, or user requested recomputation
#
proc Validate_MIC {} {
    ################
    global alignment_is_modified
    global mic_is_computed_and_valid

    if {! $mic_is_computed_and_valid} {
        set    msg "Need to (re-)compute Mutual-Information-Content."
        append msg "Proceed?"

    	set answer [tk_messageBox -parent . \
                        -title   "Warning" \
                        -type	 yesno      \
                        -icon	 warning	\
                        -message "$msg"]
        if {$answer == "no"} {
            return "abort"
        }

        # recompute
        MutualInfoFrontend
    }

}
# Validate_MIC




###   GetLogFileName   
#
#
proc GetLogFileName {proj_or_seq_filename} {
    ######################################

    set log_ext ".log"
    set flog "[file rootname $proj_or_seq_filename]$log_ext"
}
# GetLogFileName



###   GetOptStructFileName   
#
#
proc GetOptStructFileName {proj_or_seq_filename} {
    ############################################

    set log_ext ".os"
    set flog "[file rootname $proj_or_seq_filename]$log_ext"
}
# GetOptStructFileName



###   PrintDebugMat   
#
# what==cons_td|cons_dpm|seq_td|mic
#
proc PrintDebugMat {what} {
    #####################
    global FILE_TYPES
    global opts
    global selection

    set fname [tk_getSaveFile  \
                      -filetypes $FILE_TYPES(txt) \
                      -title "Choose a file for dumping" \
                      -initialdir $opts(dir,work)]
    puts "FIXME(PrintDebugMat): compute mic if requested"
    if {$fname==""}  {
        return
    }
    if {$what=="cons_td"} {
        PrintCsTdProbMat $fname
    } elseif {$what=="seq_td"} {
        PrintTdProbMat $selection(seq_no) $fname
    } elseif {$what=="mic"} {
        PrintMic $fname
    } elseif {$what=="cons_dpm"} {
        PrintCsDpm $fname $opts(td,factor) $opts(mic,factor) \
                          $opts(td,threshold) $opts(mic,Colormap,0)
    } else {
        p::error "wrong arg (what=$what)"
        return
    }
}
# PrintDebugMat




###   SaveVie   
#
#
proc SaveVie {seq_arrayname {f_id stdout}} {
    ######################################
    upvar $seq_arrayname seq

    for {set i 1} {$i<=$seq(n_seq)} {incr i} {
        puts $f_id "> $seq(id,$i)"
        puts $f_id "$seq(nt,$i)\n"
    }

}
# SaveVie



###   LoadPrecompOptStructs   
#
#
proc LoadPrecompOptStructs {fname_os os_arrayname} {
    ##############################################
    upvar $os_arrayname os

    set f_id [open $fname_os r]
    set s 0
    while {[gets $f_id line]!=-1} {

    	if { ! [string length [string trim $line]]==0} {

            # store id
            #
    		if {[regexp {> (.*)$} $line dummy id]} {
                incr s
                set os(id,$s) "$id"
            # store dotbracket
            #
    		} else {

    			set os(db,$s) "$line"
            }
    	}
    }


    close $f_id
}
# LoadPrecompOptStructs





###   OutputCalcBaseOptions   
#
#
proc OutputCalcBaseOptions {{fid stdout}} {
    #####################################

    puts $fid [CalcBaseOptionsStr]
}
# OutputCalcBaseOptions




###   CalcBaseOptionsStr   
#
# Return string for calculation options
#
proc CalcBaseOptionsStr {} {
    ######################
    global opts


    set head "---   Calculation Base   ---\n"

    set func(head) "Consensus matrix function:\n"
    set func(td)   [format "%-2.3f * TD,  if TD  > %-2.3f" \
                            $opts(td,factor) $opts(td,threshold)]

    set info(td)  "(TD  = Thermodynamic Probability)"

    if {$opts(mic,use_alifoldscore)==1} {
        set func(mic)  [format "%-2.3f * RAF, if RAF > %-2.3f" \
                            $opts(mic,factor) $opts(mic,Colormap,0)]
	    set info(mic) "(RAF = RNAalifoldscore"
        if {$opts(mic,use_stacking)} {
            append info(mic) ", incl. stacking)"
        } else {
            append info(mic) ")"
        }
        set micbase " "
    } else {
        set func(mic)  [format "%-2.3f * MIC, if MIC > %-2.3f" \
                            $opts(mic,factor) $opts(mic,Colormap,0)]
        set info(mic) "(MIC = Mutual Information Content)"

        if {$opts(mic,unbiased)} {
            set micbase "Unbiased probability estimation (Chiu & Kolodziejczak, 1991)\n"
        } else {
            set micbase "Maximum likelihood estimation (Gutell et al., 1992)\n"
        }
        if {$opts(mic,bit)} {
            append micbase "with log_2 (Schneider et al., 1986)\n"
        } else {
            append micbase "with log_e (Chiu & Kolodziejczak, 1991)\n"
        }
        if {$opts(mic,pair_entropy_norm)} {
            append micbase "Using pair entropy normalization (Martin et al., 2005)"
        }
    }
    set str "$head"
    append str "\n$func(head)"
    append str "$func(td) $info(td)\n      +\n$func(mic) $info(mic)\n"
    append str "\n$micbase"
    append str "\n\n"

    return $str
}
# CalcBaseOptionsStr



###   MapNtIndex   
#
# Maps an index from an unaligned sequence to the
# corresponding aligned one vice versa
#
# Returns:
#   the mapped index
#   -1 on error
#
proc MapNtIndex {alignedseq ntidx direction} {
    ########################################

    if {$direction!="ALIGNED_TO_DEGAPPED" && $direction!="DEGAPPED_TO_ALIGNED"} {
    	p:error "wrong arg"
    	return -1
    }
    set result -1
    set gapcount 0
    set ntcount 0

    for {set i 0} {$i<[string length $alignedseq]} {incr i} {
    	if {[IsGap [string index $alignedseq $i]]} {
    		incr gapcount
    	} else {
    		incr ntcount
    	}
    	if {$direction=="ALIGNED_TO_DEGAPPED" && $ntidx==[expr {$i+1}]} {
    		# count == number of gaps
    		set result [expr {$ntidx-$gapcount}]
    		break
    	}
    	if {$direction=="DEGAPPED_TO_ALIGNED" && $ntidx==$ntcount} {
    		# count == number of nts
    		set result [expr {$i+1}]
    		break
    	}
    }
    return $result
}
# MapNtIndex




# ---   GzipIsInstalled
#
# <FIXME:shortdescription>
#
# IN:
# OUT:
# SIDEEFFECTS:
# NOTES:
#
proc GzipIsInstalled {} {
    ###################

    if {[catch {exec gzip -h}]} {
        return 0
    } else {
        return 1
    }
}
# GzipIsInstalled
