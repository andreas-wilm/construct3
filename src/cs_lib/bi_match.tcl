#############################################################################
#
# bi_match.tcl - bmatch and imatch procedures
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
#  CVS $Id: bi_match.tcl,v 1.28 2007-10-22 10:43:22 steger Exp $
#


###   Imatch
#
#
#
#
#
proc Imatch {display_method} {
    ########################
    global opts      ;# r
    global seq       ;# r
    global selection ;# r


    if {$display_method!="circles" && $display_method!="struct_aln" && $display_method!="struct_logo"} {
        p::error "Invalid display_method \"$display_method\""
        return
    }

    if {$opts(displ,use_cons_seq)==1} {
        set seq_nt  [Get_ConsSeq]
        set seq_id "Consensus Sequence"
    } else {
        set seq_nt  $seq(nt,$selection(seq_no))
        set seq_id  $seq(id,$selection(seq_no))
    }


    if {$opts(mic,factor)>0.0} {
        if {[Validate_MIC]=="abort"} {
            return
        }
    }

    SetCursor busy


    #####
    #
    p::debug "Invoking \"CS_Imatch  $seq_nt  bpp $opts(td,factor) $opts(mic,factor) $opts(td,threshold) $opts(mic,Colormap,0) $opts(struct,allow_single_bp)\""
    set imatch_result [CS_Imatch $seq_nt  bpp  \
    					   $opts(td,factor)    $opts(mic,factor) \
    					   $opts(td,threshold) $opts(mic,Colormap,0)    \
                           $opts(struct,allow_single_bp)]
    #
    #####
    Debugger "bi_match"

    if {$display_method=="circles"} {
        set header_l1 "Imatch Circles of $cs_proj::proj(name)/$seq_id"
        set header_l2 ""

        # display structure as circles plot
        #
        CirclesFrontend bpp  $seq_id  $seq_nt "i" $header_l1 $header_l2
		  
    } elseif {$display_method=="struct_logo"} {
        StructLogoFrontend bpp

    } else {
    	set tertiary 1
        StructAlnFrontend bpp $tertiary
    }

    # generate window with ipaired nucleotides
    #
    Gen_IBMatch_BpWin "imatch" "$imatch_result"

    SetCursor normal
}
# GenImatch




####   Bmatch
#
#
# Run bmatch-algorithm, see TkCS_Bmatch.c and CS_Bmatch.c
#
proc Bmatch {display_method} {
    #######################
    global opts      ;# r
    global seq       ;# r
    global selection ;# r


    if {$display_method!="circles" && $display_method!="struct_aln" && $display_method!="struct_logo"} {
        p::error "Invalid display_method \"$display_method\""
        return
    }


    if {$opts(displ,use_cons_seq)==1} {
        set seq_nt  [Get_ConsSeq]
        set seq_id "Consensus Sequence"
    } else {
        set seq_nt  $seq(nt,$selection(seq_no))
        set seq_id  $seq(id,$selection(seq_no))
    }


    if {$opts(mic,factor)>0.0} {
        if {[Validate_MIC]=="abort"} {
            return
        }
    }

    SetCursor busy



    #####
    #
    set dmsg "Invoking \"CS_Bmatch  $seq_nt  bpp"
    append dmsg " $opts(td,factor) $opts(mic,factor)"
    append dmsg " $opts(td,threshold) $opts(mic,Colormap,0)\""
    p::debug  $dmsg
    set bmatch_result [CS_Bmatch $seq_nt  bpp \
    					   $opts(td,factor)    $opts(mic,factor) \
    					   $opts(td,threshold) $opts(mic,Colormap,0)]
    #
    #####


    if {$display_method=="circles"} {
        set header_l1 "Bmatch Circles of $cs_proj::proj(name)/$seq_id"
        set header_l2 ""

        # display structure as circles plot
        #
        CirclesFrontend bpp  $seq_id  $seq_nt "b" $header_l1 $header_l2
    } elseif {$display_method=="struct_logo"} {
        StructLogoFrontend bpp

    } else {
    	set tertiary 1
        StructAlnFrontend bpp $tertiary
    }

    # generate window with ipaired nucleotides
    #
    Gen_IBMatch_BpWin "bmatch" "$bmatch_result"

    SetCursor normal
}
# GenBmatch




###   Gen_IBMatch_BpWin
#
# Popup window with i/b/matched bp's
# match_method: imatch|bmatch
#
#
proc Gen_IBMatch_BpWin {match_method match_result_str} {
    ##################################################


    if {$match_method=="imatch"} {
        set title_str "Imatched Basepairs"
        set w .ibp
    } elseif {$match_method=="bmatch"}  {
    	set title_str "Bmatched Basepairs"
    	set w .bbp
    } else {
    	p::error "Invalid match_method \"$match_method\""
        return
    }


    if [winfo exists $w] { destroy $w }
    toplevel 	$w

   	wm title    $w "$title_str"
  	wm iconname $w "$title_str"


    frame  $w.buttons
    pack   $w.buttons         -side bottom -fill x -pady 2m

    button $w.buttons.dismiss -text Dismiss -command "destroy $w"
    button $w.buttons.save    -text Save    -command "SaveBplistW  $match_method $w"
    button $w.buttons.print   -text Print   -command "PrintBplistW $match_method $w"
    pack   $w.buttons.save $w.buttons.print $w.buttons.dismiss -side left -expand 1

    text $w.text -yscrollcommand "$w.scroll set" -setgrid true \
                 -width 78 -height 32 -wrap none; # -wrap word
    scrollbar $w.scroll -command "$w.text yview"
    pack $w.scroll -side right -fill y
    pack $w.text   -expand yes -fill both


    # Set up display styles
    #
    $w.text tag configure bold    -font {Courier 12 bold italic}
    $w.text tag configure big     -font {Courier 14 bold}
    $w.text tag configure verybig -font {Helvetica 24 bold}
    $w.text tag configure right   -justify right
    $w.text tag configure center  -justify center
    $w.text tag configure margins -lmargin1 12m -lmargin2 6m -rmargin 10m
    $w.text tag configure spacing -lmargin1 12m -lmargin2 6m -rmargin 10m \
    	-spacing1 10p -spacing2 2p 

    $w.text insert end "$title_str" "big center"
    $w.text insert end "\n"
    $w.text insert end  "$match_result_str" margins
    $w.text insert end "\n"
}
# Gen_IBMatch_BpWin





####   SaveBplistW
#
#
# Save bp-list from imatch/bmatch
#   $matchname: 'imatch' or 'bmatch'
#   $w:    window with text to save
#
#
proc SaveBplistW {matchname w} {
    #########################
    global FILE_TYPES
    global opts

    if {$matchname!="imatch" && $matchname!="bmatch"} {
        p::error "Invalid match method \"$matchname\""
        return
    }

    set init_file "${cs_proj::proj(name)}_${matchname}.txt"

    set filename [tk_getSaveFile \
                   -filetypes $FILE_TYPES(txt) \
                   -initialfile "$init_file" \
                   -initialdir  $opts(dir,work) \
                   -title       "Save $matchname File"]
    update idletasks
    if {$filename == ""} {
    	return
    }


    DumpTextWindowContent $w $filename
}
# SaveBplistW
