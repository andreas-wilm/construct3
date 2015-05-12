##############################################################################
#
# sub_gui.tcl - some extra popup gui stuff
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
#  CVS $Id: sub_gui.tcl,v 1.39 2007-10-26 13:51:13 steger Exp $    
#

###   ChangeProbThreshold   
#
#
proc ChangeProbThreshold {menuW} {
#########################
    global opts
    global w

    set f .limit
    if {[winfo exists $f]==0} {
        set oldFocus [focus]
        set row 0
        toplevel    $f
        wm title    $f "Change thermodynamic probability threshold"
        wm iconname $f ChangeProbThreshold

        set    info_txt "Consensus basepairs with a probabilities below treshold will be ignored."
        append info_txt "Treshold value must be between 0.0 and 1.0."
        append info_txt "\n\n(Tip: Clicking a consensus basepair with middle mouse button outputs probability to console)\n"

        label $f.message  -wraplength 3i  -justify left -text "$info_txt"
        grid  $f.message  -row $row; incr row
        label $f.leer$row -text " "
        grid  $f.leer$row -row $row; incr row

        frame $f.val    -bd 2
        label $f.val.label -text "Threshold = "
        entry $f.val.entry -relief sunken -width 5
        grid  $f.val.label -column 0 -row 0 -sticky w
        grid  $f.val.entry -column 1 -row 0 -sticky e
        grid  $f.val -row $row; incr row
        bind  $f.val.entry <Return> "SetProbThreshold $f.val.entry $menuW"
        $f.val.entry insert 0 [format "%-2.3f" $opts(td,threshold)]

        label  $f.leer$row    -text " "
        grid   $f.leer$row -row $row; incr row
        frame  $f.buttons
        grid   $f.buttons  -row $row; incr row
        button $f.buttons.dismiss -text Dismiss  -command "CancelEntry  $f $oldFocus"
        button $f.buttons.accept  -text Return   -command "SetProbThreshold $f.val.entry $menuW; CancelEntry $f $oldFocus"
        grid   $f.buttons.accept  -column 0 -row 0
        grid   $f.buttons.dismiss -column 1 -row 0
    }

    set pos [winfo pointerxy $w(dp)]
    wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
    raise $f
    grab  $f
    focus $f
}
# ChangeThreshold




###   CancelEntry   
#
proc CancelEntry {f oldFocus} {
################################
    focus        $oldFocus
    grab release $f
    destroy      $f
}
# CancelEntry




###   SetProbThreshold   
#
proc SetProbThreshold {entryW menuW} {
######################
   global opts
   global MENUINDEX

    # check value
    set x [expr {0.001*int([$entryW get]*1000.)}]
    set x [expr {($x>1.0) ? 1.0 : $x}]
    set x [expr {($x<0.0) ? 0.0 : $x}]

    $entryW delete 0 end
    $entryW insert 0 [format "%-2.3f" $x]

    set opts(td,threshold) $x

    $menuW entryconfigure $MENUINDEX(CalcBase,ThreshTd) \
           -label [format "  TD Threshold = %-2.3f" $opts(td,threshold)]
}
# SetProbThreshold




###   ShowHelpHint   
#
#
#
proc ShowHelpHint {} {
###################

    set title "Help Hint"
    set msg   "The manual is provided as separate postscript file."

    tk_messageBox  -title $title  -type ok \
                   -icon  info    -message  $msg

}
# ShowHelpHint



###   ShowVersion   
#
#
#
proc ShowVersion {} {
###################

    set title "Version"
    set msg   [Version]

    tk_messageBox  -title $title  -type ok \
                   -icon  info    -message  $msg

}
# ShowVersion



###   ChooseProject   
#
# Get filename of ProjectFile from user and startup.
#
proc ChooseProject {} {
##################
    global FILE_TYPES
    global opts

    set projfile [tk_getOpenFile  \
                      -filetypes $FILE_TYPES(proj) \
                      -title "Choose a project" \
                      -initialdir $opts(dir,work)]


    if {$projfile==""}  {
        return ERROR
    }
    if {![file exists $projfile]} {
        return ERROR
    }

    set opts(dir,work) [file dirname $projfile]

    if {[LoadProject $projfile]!="success"} {
            set msg "Failed to load project-file \"$projfile\"!"
            tk_messageBox -message $msg -type ok -title "Error"
    }
}
# ChooseProject





###   SaveAln   
#
# Saves the modified alignment as vienna
# Pops up a frontend if no filename is specified
#
#
proc SaveAln {{fname ""}} {
############
    global seq
    global opts
    global FILE_TYPES
    global alignment_is_modified

    set init_file "${cs_proj::proj(name)}_construct.vie"

    if {$fname==""} {
        set fname [tk_getSaveFile -title "Save Alignment" \
                     -filetypes   $FILE_TYPES(seqout) \
                     -initialdir  $opts(dir,work) \
                     -initialfile $init_file]
        if {$fname==""} {
            return
        }
    }

    set opts(dir,work) [file dirname $fname]

    # determine format
    #
    set match 0
    set out_ext [file extension $fname]

    # intercept vienna:
    # unsupported by seqio, therefore done by tcl
    if {$out_ext!=".vie"} {
        foreach type $FILE_TYPES(seqout) {
            set name [lindex $type 0]
            set ext  [lindex $type 1]
            if {$out_ext==$ext} {
                set match 1
                break
            }
        }
        if {$match==1} {
            ConvertAln $cs_proj::proj(aln_path) $fname $out_ext
            set alignment_is_modified 0

            return ;# !!!
        } else {
            p::warn "Couldn't determine format, using vienna"
        }
    }

    set fid [open $fname w+]
    # FIXME: fix seqio for use of vienna
    SaveVie seq $fid
    close $fid

    set alignment_is_modified 0
}
# SaveAln





###   MIC_ShowColorMapWin   
#
# Popup a window to select color limits for GenStatInfo
# proc StatInfoShowColorW {}
#
proc MIC_ShowColorMapWin {} {
########################
    global opts ;# r
    global w    ;# r


    set opts(mic,oldlimit,0)    $opts(mic,limit,0)
    set opts(mic,oldlimit,1)    $opts(mic,limit,1)
    set opts(mic,oldColormap,0) $opts(mic,Colormap,0)
    set opts(mic,oldColormap,1) $opts(mic,Colormap,1)

    set w(mic_colormap) .colorsel
    catch {destroy $w(mic_colormap)}
    toplevel $w(mic_colormap)
    wm title $w(mic_colormap)    "Colormapping"
    wm iconname $w(mic_colormap) "Colormapping"
    wm resizable $w(mic_colormap) 0 0

    label $w(mic_colormap).l -width 20 -text "Colormapping"

    canvas $w(mic_colormap).ct -width 130 -height 2
    pack $w(mic_colormap).l $w(mic_colormap).ct -side top

    $w(mic_colormap).ct create line 0 1 130 1 -width 2

    frame $w(mic_colormap).buttons
    canvas $w(mic_colormap).cb -width 130 -height 1
    pack $w(mic_colormap).buttons $w(mic_colormap).cb -side bottom  -pady 1m

    $w(mic_colormap).cb create line 0 1 130 1 -width 2

    button $w(mic_colormap).buttons.apply -text Apply   -command "MIC_ApplyWin"
    button $w(mic_colormap).buttons.cancel -text Cancel -command "MIC_CancelWin"
    pack $w(mic_colormap).buttons.apply $w(mic_colormap).buttons.cancel -side left -expand 1 -padx 4

    set cl $w(mic_colormap).cl; # left canvas, for adjusting left distance to border
    set c0 $w(mic_colormap).c0; # left caption
    set c1 $w(mic_colormap).c1; # left rectangle, both arrows and the bars
    set c2 $w(mic_colormap).c2; # right rectangle incl. caption

    canvas $cl -width 2  -height 140
    canvas $c0 -width 22 -height 140
    canvas $c1 -width 45 -height 140
    canvas $c2 -width 53 -height 140
    pack $cl $c0 $c1 $c2 -side left -fill x -expand true

    set plotFont {Helvetica 10}


    # left caption
    #
    foreach point {{2 7 0.00} {2 39 0.25} {2 71 0.50} {2 103 0.75} {2 135 1.00}} {
    	$c0 create text  [lindex $point 0] [lindex $point 1] -text [lindex $point 2] -anchor w -font $plotFont
    }

    # left rectangle
    #
    $c1 create rectangle 4 6 25 136
    for {set i 0} {$i<=128} {incr i} {
    	$c1 create line 5 [expr {7+$i}] 25 [expr {7+$i}] -fill [GetRgb [expr {double($i)/128.0}] 0 1] -tags recl
    }

    # both rectangles for shadowing
    #
    $c1 create rectangle 5   6 25 [expr {$opts(mic,limit,0)-1}] -fill white -outline "" -tags rec0
    $c1 create rectangle 5 136 25 [expr {$opts(mic,limit,1)+1}] -fill black -outline "" -tags rec1
    $c1 create line 4 6 25 6

    # both black arrows
    #
    $c1 create polygon 27 $opts(mic,limit,0) 33 \
    	[expr {$opts(mic,limit,0)-4}] 33 [expr {$opts(mic,limit,0)+4}] -fill black -tags arrow0
    $c1 create polygon 27 $opts(mic,limit,1) 33 \
    	[expr {$opts(mic,limit,1)-4}] 33 [expr {$opts(mic,limit,1)+4}] -fill black -tags arrow1
    $c1 create line 33 $opts(mic,limit,0) 45 7 -tags line0
    $c1 create line 33 $opts(mic,limit,1) 45 135 -tags line1

    # right caption
    #
    label $c2.l0 -font $plotFont -text [format "%.2f" [expr {($opts(mic,limit,0)-7.)/128.0}]]
    label $c2.l1 -font $plotFont -text [format "%.2f" [expr {($opts(mic,limit,1)-7.)/128.0}]]
    label $c2.lm -font $plotFont -text [format "%.2f" [expr {($opts(mic,limit,1)+$opts(mic,limit,0)-14)/256.0}]]
    $c2 create window 0   7 -window $c2.l0 -anchor w
    $c2 create window 0 135 -window $c2.l1 -anchor w
    $c2 create window 0  71 -window $c2.lm -anchor w

    # right rectangle
    #
    $c2 create rectangle 28 6 49 136
    for {set i 0} {$i<=128} {incr i} {
    	$c2 create line 29 [expr {7+$i}] 49 [expr {7+$i}] -fill [GetRgb [expr {double($i)/128.0}] 0 1]
    }

    # bindings
    #
    $c1 bind arrow0 <1>               "MIC_ColorDown  $c1 %x %y 0"
    $c1 bind arrow1 <1>               "MIC_ColorDown  $c1 %x %y 1"
    $c1 bind arrow0 <ButtonRelease-1> "MIC_ColorLeave $c1 0"
    $c1 bind arrow1 <ButtonRelease-1> "MIC_ColorLeave $c1 1"
    $c1 bind arrow0 <B1-Motion>       "MIC_ColorMove  $c1 $c2 %x %y 0"
    $c1 bind arrow1 <B1-Motion>       "MIC_ColorMove  $c1 $c2 %x %y 1"
}
# MIC_ShowColorMapWin




###   MIC_ColorDown   
#
# proc StatInfoColorDown {win x y z}
#
proc MIC_ColorDown {win x y z} {
##############################
    global opts

    $win dtag selected
    $win addtag selected withtag current
    $win raise current

    set opts(mic,limit,$z) $y
}
# MIC_ColorDown




###   MIC_ColorMove   
#
# proc StatInfoColorMove {c1 c2 x y z}
#
#
proc MIC_ColorMove {c1 c2 x y z} {
################################
    global opts

    if {$y >= 7 && $y < 136} {
    	if {($z==0 && $y<[expr {$opts(mic,limit,1)-2}])  \
                              ||                       \
            ($z==1 && $y>[expr {$opts(mic,limit,0)+2}])} {
    		$c2.l$z configure -text [format "%.2f" [expr {($y-7.)/128.}]]
    		if {$z==0} {
    			$c1 coords line0  33 $y 45 7
    			$c1 coords arrow0 27 $y 33 [expr {$y-4}] 33 [expr {$y+4}]
    			$c1 coords rec0    5  6 25 $y
    		} else {
    			$c1 coords line1  33 $y 45 135
    			$c1 coords arrow1 27 $y 33 [expr {$y-4}] 33 [expr {$y+4}]
    			$c1 coords rec1    5 [expr {$y+1}] 25 136
    		}
    		set opts(mic,limit,$z) $y
    		$c2.lm configure -text [format "%.2f" [expr {($opts(mic,limit,1)+$opts(mic,limit,0)-14)/256.0}]]
    	}
    }
}
# MIC_ColorMove




###   MIC_ColorLeave   
#
# proc StatInfoColorLeave {c1 z}
#
proc MIC_ColorLeave {c1 z} {
###########################
    global opts  ;#w

    $c1 dtag selected

    set opts(mic,Colormap,$z) [expr {($opts(mic,limit,$z)-7.0)/128.0}]
}
# MIC_ColorLeave




###   MIC_CancelWin   
#
# proc StatInfoCancelW {}
#
proc MIC_CancelWin {} {
#####################
    global opts  ;# w
    global w     ;# r

    set opts(mic,limit,0) $opts(mic,oldlimit,0)
    set opts(mic,limit,1) $opts(mic,oldlimit,1)
    set opts(mic,Colormap,0) $opts(mic,oldColormap,0)
    set opts(mic,Colormap,1) $opts(mic,oldColormap,1)

    destroy $w(mic_colormap)
}
# MIC_CancelWin




###   MIC_ApplyWin   
#
# proc StatInfoApplyW {}
#
proc MIC_ApplyWin {} {
#################
    global opts
    global mic_is_computed_and_valid
    global w
    global MENUINDEX

    set mic_thresh_label $w(dp).f_menu.calc_base.menu

    set opts(mic,oldlimit,0)	 $opts(mic,limit,0)
    set opts(mic,oldlimit,1)	 $opts(mic,limit,1)
    set opts(mic,oldColormap,0)  $opts(mic,Colormap,0)
    set opts(mic,oldColormap,1)  $opts(mic,Colormap,1)

    set mic_is_computed_and_valid 0

    # update menu label
    $mic_thresh_label entryconfigure $MENUINDEX(CalcBase,ThreshMic) \
           -label [format "  MIC Threshold = %-2.3f" $opts(mic,Colormap,0)]


    # and recompute mic
    MutualInfoFrontend

}
# MIC_ApplyWin




###   ChangeCalcBaseFactors   
#
proc ChangeCalcBaseFactors {menupath} {
##########################
    global opts
    global w


    set f .calc_base_factors
    if {[winfo exists $f]==0} {
        set oldFocus [focus]
        set row 0
        toplevel    $f
        wm title    $f "Calculation Factors"
        wm iconname $f ChangeCalcBaseFactors

        set info_txt "Change weighting factors for structure prediction"
        append info_txt "\n\t0.0 =< factor =< 1.0"

        label $f.message  -wraplength 3i  -justify left -text "$info_txt"
        grid  $f.message  -row $row; incr row
        label $f.leer$row -text " "
        grid  $f.leer$row -row $row; incr row

        frame $f.val    -bd 2
            label $f.val.label_td -text "Thermodyn.Basepair-Prob. = "
            entry $f.val.entry_td -relief sunken -width 5
            grid  $f.val.label_td -column 0 -row 0 -sticky w
            grid  $f.val.entry_td -column 1 -row 0 -sticky e
            grid  $f.val -row $row; incr row
            bind  $f.val.entry_td <Return> "SetCalcFactor  td  $f.val.entry_td  $f.val.entry_mic"
            bind  $f.val.entry_td <Leave>  "SetCalcFactor  td  $f.val.entry_td  $f.val.entry_mic"
            $f.val.entry_td insert 0 [format "%-2.3f" $opts(td,factor)]


            label $f.val.label_mic -text "Mutual-Info.-Content = "
            entry $f.val.entry_mic -relief sunken -width 5
            grid  $f.val.label_mic -column 0 -row 1 -sticky w
            grid  $f.val.entry_mic -column 1 -row 1 -sticky e
            grid  $f.val -row $row; incr row
            bind  $f.val.entry_mic <Return> "SetCalcFactor mic $f.val.entry_td  $f.val.entry_mic"
            bind  $f.val.entry_mic <Leave>  "SetCalcFactor mic $f.val.entry_td  $f.val.entry_mic"
            $f.val.entry_mic insert 0 [format "%-2.3f" $opts(mic,factor)]


        label  $f.leer$row    -text " "
        grid   $f.leer$row -row $row; incr row
        frame  $f.buttons
        grid   $f.buttons  -row $row; incr row

        button $f.buttons.accept  -text "Return" \
               -command "SetCalcFactor  td  $f.val.entry_td  $f.val.entry_mic; \
                         ApplyCalcFactor $f.val.entry_td $f.val.entry_mic $menupath; \
                         CancelEntry $f $oldFocus;"
        button $f.buttons.dismiss -text "Dismiss" \
               -command "CancelEntry  $f $oldFocus"
        #
        grid   $f.buttons.dismiss -column 1 -row 0
        grid   $f.buttons.accept  -column 0 -row 0
    }

    set pos [winfo pointerxy $w(dp)]
    wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
    raise $f
    grab  $f
    focus $f
}
# ChangeCalcBaseFactors





###   SetCalcFactor   
#
# what: mic|td
#
proc SetCalcFactor {what  td_entry  mic_entry} {
##################

    if {$what=="td"} {
        set current_entry $td_entry
        set chained_entry $mic_entry
    } elseif {$what=="mic"} {
        set current_entry $mic_entry
        set chained_entry $td_entry
    } else {
        p::error "Invalid option \"$what\""
        return
    }

    # check value
    set x [expr {0.001*int([$current_entry get]*1000.)}]
    set x [expr {($x>1.0) ? 1.0 : $x}]
    set x [expr {($x<0.0) ? 0.0 : $x}]
    #
    $current_entry delete 0 end
    $current_entry insert 0 [format "%-2.3f" $x]

    $chained_entry delete 0 end
    $chained_entry insert 0 [format "%-2.3f" [expr {1.0-$x}]]
}
# SetCalcFactor




###   ApplyCalcFactor   
#
#
proc ApplyCalcFactor {td_entry mic_entry calcfunc_menu} {
####################
    global opts
    global MENUINDEX

    set opts(mic,factor) [format "%-2.3f" [$mic_entry get]]
    set opts(td,factor)  [format "%-2.3f" [$td_entry get]]

    # update menu label
    $calcfunc_menu entryconfigure $MENUINDEX(CalcBase,Function) \
                   -label [format "  Function = %-2.3f*TD + %-2.3f*MIC" $opts(td,factor) $opts(mic,factor)]
}
# ApplyCalcFactor



###   ChangeSuboptStructNum   
#
#
proc ChangeSuboptStructNum {menuW} {
#########################
    global opts
    global w



    set f .suboptstructs
    if {[winfo exists $f]==0} {
        set oldFocus [focus]
        set row 0
        toplevel    $f
        wm title    $f "Change number of suboptimal structures"
        wm iconname $f ChangeSuboptStructNum

        set    info_txt "Change number of suboptimal structures"
        append info_txt "which will be computed"

        label $f.message  -wraplength 3i  -justify left -text "$info_txt"
        grid  $f.message  -row $row; incr row
        label $f.leer$row -text " "
        grid  $f.leer$row -row $row; incr row

        frame $f.val    -bd 2
        label $f.val.label -text "Number = "
        entry $f.val.entry -relief sunken -width 5
        grid  $f.val.label -column 0 -row 0 -sticky w
        grid  $f.val.entry -column 1 -row 0 -sticky e
        grid  $f.val -row $row; incr row
        bind  $f.val.entry <Return> "SetSuboptStructNum $f.val.entry $menuW"
        $f.val.entry insert 0 $opts(suboptstruct,number)

        label  $f.leer$row    -text " "
        grid   $f.leer$row -row $row; incr row
        frame  $f.buttons
        grid   $f.buttons  -row $row; incr row
        button $f.buttons.dismiss -text Dismiss  -command "CancelEntry  $f $oldFocus"
        button $f.buttons.accept  -text Return   -command "SetSuboptStructNum $f.val.entry $menuW; CancelEntry $f $oldFocus"
        grid   $f.buttons.accept  -column 0 -row 0
        grid   $f.buttons.dismiss -column 1 -row 0
    }

    set pos [winfo pointerxy $w(dp)]
    wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
    raise $f
    grab  $f
    focus $f
}
# ChangeSuboptStructNum



###   SetSuboptStructNum   
#
proc SetSuboptStructNum {entryW menuW} {
#######################
   global opts
   global MENUINDEX

    # check value
    set x [format "%d" [$entryW get]]

    $entryW delete 0 end
    $entryW insert 0 $x

    set opts(suboptstruct,number) $x

    set txt [format "Number of suboptimal structures = %d" $opts(suboptstruct,number)]
    $menuW entryconfigure $MENUINDEX(Options,SuboptStructNum) -label $txt
}
# SetSuboptStructNum




###   StructLogoFrontend
#
#
proc StructLogoFrontend {bpp_arrayname {tertiary 0}} {
#######################
    global opts
    global seq
    global mic_is_computed_and_valid ;# r

    upvar $bpp_arrayname bpp

    SetCursor busy

    set csseq [Get_ConsSeq]
    LogoAln seq $csseq bpp $opts(dir,work) ${cs_proj::proj(name)}
	 
    SetCursor normal
}



###   LogoAln
#
#
proc LogoAln {alnseq_array cons_seq bpp_array workdir title} {
############
    upvar $alnseq_array seq
    upvar $bpp_array    bpp

    if {[catch {package require TclCurl}]} {
        set msg "You need the TclCurl for this feature."
        tk_messageBox -message $msg -type ok -title "Error"
        return
    }

    set bplist {}
    for {set i 1} {$i<=[string length $cons_seq]} {incr i} {
    	lappend bplist $bpp($i,partner)
    }
	set struct_alnnum [StructTransform $bplist $seq(aln_len)]
	 
	logo::frontEnd $struct_alnnum  $title $workdir 
	
# Hier gehoert jetzt eine Prozedur (inkl. tk-Widget hin, 
#   die die Benutzereingaben abfragt und in den folgenden 
#   String einbaut (siehe /usr/share/doc/tk8.4/examples/widget und 
#   sub_gui.tcl:ChangeCalcBaseFactors)
#   Fuer Angabe von zwei Teilstrings, um z.B. eine Helix 
#   aus der Mitte der Struktur zu beurteilen, muss das 
#   Logo-Alignment aus dem kompletten Alignment zusammengestellt
#   werden!
#    set    logoInput "ntprobs=0.25++0.25++0.25++0.25&probAU=1.0&probCG=1.0&probGU=1.0&strucseq=Y&assignment=Y&logotype=2&startpos=1&poslabel=Y&TorUs=N&comments="
#    append logoInput "%3E+";    # "> "
#    append logoInput [StructTransform $bplist $seq(aln_len)]
#    append logoInput "%0D%0A";  # "<cr>"
#    for {set s 1} {$s<=$seq(n_seq)} {incr s} {
#        append logoInput "%3E+"
#        append logoInput [string toupper $seq(nt,$s)]
#        append logoInput "%0D%0A"
#    }
#	puts "DEBUG: logoInput $logoInput"

}





###   ShowSeqSearchWin   
#
proc ShowSeqSearchWin {} {
######################
    global w
    global seq
    global seqsearch_re

    set    help_txt "Regular Expressions allow very powerful string matching."
    append help_txt "\nThe preinserted regular expression describes"
    append help_txt " extra stable tetraloops GNRA and UNCG.\n\n"


    set f .seqsearch
    if {[winfo exists $f]==0} {
        set oldFocus [focus]
        set row 0
        toplevel    $f
        wm title    $f "Sequence Search"
        wm iconname $f SeqSearch

        set    info_txt "Search in unaligned sequences\n"
        append info_txt "(using regular expressions)"


        label  $f.message  -wraplength 3i  -justify left -text "$info_txt"
        grid   $f.message  -row $row; incr row
        button $f.help     -text "Help" \
                           -command "RegExpHelpWindow \"$help_txt\""
        grid   $f.help -row $row; incr row
        label  $f.leer$row -text " "
        grid   $f.leer$row -row $row; incr row

        frame  $f.val    -bd 2
        label  $f.val.l_for  -text "Search for "
        entry  $f.val.entry  -relief sunken
        #
        label   $f.val.l_in   -text "Search in "
        listbox $f.val.idlist -selectmode extended
        #
        $f.val.idlist insert end all
        for {set i 1} {$i<=$seq(n_seq)} {incr i} {
            $f.val.idlist insert end "$seq(id,$i)"
        }
        $f.val.idlist selection set 0;# first index : all
        $f.val.entry insert 0 $seqsearch_re(eshp_re)
        #
        button $f.val.color  -text "color" -bg $seqsearch_re(cur_col) \
                             -command "SetReMatchColor $f.val.color"

        button $f.val.find   -text "Find & Mark" \
                             -command "SeqSearchFrontend $f.val.entry $f.val.idlist"
        button $f.val.clear  -text "Clear All" \
                             -command {PatternMatchDeleteAll}
        button $f.val.exit   -text "Return" \
                             -command "CancelEntry  $f $oldFocus"
        #
        grid  $f.val.l_for  -column 0 -row $row -sticky w
        grid  $f.val.entry  -column 1 -row $row -sticky e
        incr row
        #
        grid  $f.val.l_in   -column 0 -row $row -sticky w
        grid  $f.val.idlist -column 1 -row $row -sticky e
        incr row
        #
        grid  $f.val.color -column 0 -row $row -sticky w
        grid  $f.val.find  -column 1 -row $row -sticky we
        incr row
        grid  $f.val.clear -column 1 -row $row -sticky we
        incr row
        grid  $f.val.exit  -column 1 -row $row -sticky we
        incr row
        #
        grid  $f.val
        #
        bind  $f.val.entry <Return> "SeqSearchFrontend $f.val.entry $f.val.idlist"
    }

    set pos [winfo pointerxy $w(dp)]
    wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
    raise $f
    # no grab:  alignment scrolling is allowed
    focus $f
}
# ShowSeqSearchWin



###   SetReMatchColor   
#
#
proc SetReMatchColor {{widget ""}} {
####################
    global seqsearch_re


    set color [tk_chooseColor -initialcolor $seqsearch_re(cur_col) -title "Choose color"]
    if {$color==""} {
        return
    }
    set seqsearch_re(cur_col) $color
    if {$widget!=""} {
        $widget configure -bg $color
    }
}
# SetReMatchColor




###   SeqSearchFrontend   
#
#
proc SeqSearchFrontend {re_widget id_widget} {
######################
    global seq
    global seqsearch_re


    set re [$re_widget get]
    if {$re==""} {
        return
    }

    set seq_idxs [$id_widget curselection]
    if {![llength $seq_idxs]} {
        return
    } elseif {[lsearch -exact $seq_idxs 0]!=-1} {
        set seq_idxs {}
        for {set i 1} {$i<=$seq(n_seq)} {incr i} {
            lappend seq_idxs $i
        }
    }

    SetCursor busy


    set seqsearch_re(cur_re) $re
    lappend seqsearch_re(color_list)   $seqsearch_re(cur_col)
    lappend seqsearch_re(re_list)      $re
    lappend seqsearch_re(seqidx_lists) $seq_idxs

    p::debug "seqsearch_re(color_list)   = $seqsearch_re(color_list)"
    p::debug "seqsearch_re(re_list)      = $seqsearch_re(re_list)"
    p::debug "seqsearch_re(seqidx_lists) = $seqsearch_re(seqidx_lists)"

    SeqSearchAndMark $seq_idxs $re $seqsearch_re(cur_col)

    SetCursor normal
}
# SeqSearchFrontend






###   SeqSearchAndMark   
#
#
# seq_idxs: list of indices for seqs to mark
# re: list of different RegExps to apply
# col: corrsponding colors to re_list
#
proc SeqSearchAndMark {seq_idxs re col} {
######################
    global seq
    global seqsearch_re

    if {![info exists seqsearch_re]} {
        return
    }

    foreach s $seq_idxs {
        p::debug "searching in $seq(nt,$s)"
        set idxl [SeqSearch $seq(nt,$s) $re]
        if {[llength $idxl]==0} {
            continue
        }

        foreach idx $idxl {
            set start [lindex $idx 0]
            set end   [lindex $idx 1]
            set match_str    "match in $seq(id,$s) at"
            append match_str " [expr {$start+1}]-[expr {$end+1}]:"
            append match_str " [string range $seq(nt,$s) $start $end]"
            PatternMatchMark $s $col [expr {$start+1}] [expr {$end+1}]

            p::verbose "$match_str"
        }
    }
}
# SeqSearchAndMark




###   SeqSearch   
#
# lookup regexp <re>
# in dealigned seq and return match start and end indices (zero-offset)
# in aligned sequence as list of list
# or empty list on mismatch
#
proc SeqSearch {alnseq re} {
##########################

    set nalseq [Degap $alnseq]

    set offset 0
    set raw_idxs_l {}
    while {[regexp -indices -nocase -start $offset -- $re $nalseq idxl]} {
        p::debug "idxl = $idxl"
        lappend raw_idxs_l $idxl
        set offset [expr {[lindex $idxl 0]+1}]
        p::debug "offset = $offset"
        set highest_idx [lindex $idxl 1];# for lasy evaluation of gapshift
    }

    if {![llength $raw_idxs_l]} {
        return {}
    }


    ###   precompute gap shift til highest index
    #
    set introd_gaps 0
    set j 0
    for {set i 0} {$i<[string length $alnseq]} {incr i} {
        if {[IsGap [string index $alnseq $i]]} {
            incr introd_gaps;
        } else {
            set gap_shift($j) $introd_gaps
            p::debug "gap_shift($j)=$introd_gaps"
            incr j
            # be lazy....
            if {[expr {$j-1}]>$highest_idx} {
                break
            }
        }
    }


    set ret_idx_l {}
    foreach idxl $raw_idxs_l {
        set aln_start [expr {[lindex $idxl 0]+$gap_shift([lindex $idxl 0])}]
        set aln_end   [expr {[lindex $idxl 1]+$gap_shift([lindex $idxl 1])}]
        p::debug " idxl=$idxl -real-> $aln_start $aln_end [string range $alnseq $aln_start $aln_end]"

        lappend ret_idx_l [list $aln_start $aln_end]
    }
    return $ret_idx_l
}
# SeqSearchFrontend




###   PatternMatchMark   
#
#
proc PatternMatchMark {seq_idx col nt_idx1 nt_idx2} {
######################
    global c

    set c1 [GetNtAliCoords $seq_idx $nt_idx1]
    set c2 [GetNtAliCoords $seq_idx $nt_idx2]

    set cell_width  [GetNtAliCellDim width]
    set cell_height [GetNtAliCellDim height]

    set x1 [expr {[lindex $c1 0]-$cell_width/2}]
    set y1 [expr {[lindex $c1 1]-$cell_height/2}]

    set x2 [expr {[lindex $c2 0]+$cell_width/2}]
    set y2 [expr {[lindex $c1 1]+$cell_height/2}]

    set box [$c(al) create rectangle ${x1}c ${y1}c ${x2}c ${y2}c \
                                     -tags "REMatchBox REMatchBox_seq_${seq_idx}" \
                                     -fill $col -outline "" ]

    $c(al) lower $box

}
# PatternMatchMark




###   PatternSearchActive   
#
#
proc PatternSearchActive {} {
########################
    global c

    if {[llength [$c(al) find withtag REMatchBox]]} {
        return 1
    } else {
        return 0
    }
}
# PatternSearchActive




###   PatternMatchDelete   
#
# remove box from regexps matched aln-nts
# optionally: one for one sequeunce
#
proc PatternMatchDelete {{seq_idx -1}} {
#######################
   global c

   if {$seq_idx==-1} {
       $c(al) delete REMatchBox
    } else {
       $c(al) delete REMatchBox_seq_${seq_idx}
    }
}
# PatternMatchDelete




###   PatternMatchDeleteAll   
#
# remove all boxes and the saved regexps, seqlists etc.
#
proc PatternMatchDeleteAll {} {
#######################
   global c
   global seqsearch_re

    PatternMatchDelete

    unset seqsearch_re(re_list)
    unset seqsearch_re(color_list)
    unset seqsearch_re(seqidx_lists)

}
# PatternMatchDelete



###   RegExpHelpWindow   
#
#
proc RegExpHelpWindow {extra_txt} {
######################

    set w      .regexp_help
    set header "Regular Expression Help"

    # escape brackets and curlybraces
    set re_help_txt "The following text is extracted from Tcl/TK 8.3 re_syntax documentation:


Syntax of Tcl regular expressions.


DESCRIPTION
A regular expression describes strings of characters. It's a pattern that matches certain strings and doesn't match others.


DIFFERENT FLAVORS OF REs
Regular expressions (``RE''s), as defined by POSIX, come in two flavors: extended REs (``EREs'') and basic REs (``BREs''). EREs are roughly those of the traditional egrep, while BREs are roughly those of the traditional ed. This implementation adds a third flavor, advanced REs (``AREs''), basically EREs with some significant extensions.

This manual page primarily describes AREs. BREs mostly exist for backward compatibility in some old programs; they will be discussed at the end. POSIX EREs are almost an exact subset of AREs. Features of AREs that are not present in EREs will be indicated.


REGULAR EXPRESSION SYNTAX
Tcl regular expressions are implemented using the package written by Henry Spencer, based on the 1003.2 spec and some (not quite all) of the Perl5 extensions (thanks, Henry!). Much of the description of regular expressions below is copied verbatim from his manual entry.

An ARE is one or more branches, separated by `|', matching anything that matches any of the branches.

A branch is zero or more constraints or quantified atoms, concatenated. It matches a match for the first, followed by a match for the second, etc; an empty branch matches the empty string.

A quantified atom is an atom possibly followed by a single quantifier. Without a quantifier, it matches a match for the atom. The quantifiers, and what a so-quantified atom matches, are:

*
a sequence of 0 or more matches of the atom

+
a sequence of 1 or more matches of the atom

?
a sequence of 0 or 1 matches of the atom

\{m\}
a sequence of exactly m matches of the atom

\{m,\}
a sequence of m or more matches of the atom

\{m,n\}
a sequence of m through n (inclusive) matches of the atom; m may not exceed n

*? +? ?? \{m\}? \{m,\}? \{m,n\}?
non-greedy quantifiers, which match the same possibilities, but prefer the smallest number rather than the largest number of matches (see MATCHING)

The forms using \{ and \} are known as bounds. The numbers m and n are unsigned decimal integers with permissible values from 0 to 255 inclusive.

An atom is one of:

(re)
(where re is any regular expression) matches a match for re, with the match noted for possible reporting

(?:re)
as previous, but does no reporting (a ``non-capturing'' set of parentheses)

()
matches an empty string, noted for possible reporting

(?:)
matches an empty string, without reporting

\[chars\]
a bracket expression, matching any one of the chars (see BRACKET EXPRESSIONS for more detail)

.
matches any single character

k
(where k is a non-alphanumeric character) matches that character taken as an ordinary character, e.g.  matches a backslash character

c
where c is alphanumeric (possibly followed by other characters), an escape (AREs only), see ESCAPES below

\{
when followed by a character other than a digit, matches the left-brace character `\{'; when followed by a digit, it is the beginning of a bound (see above)

x
where x is a single character with no other significance, matches that character.

A constraint matches an empty string when specific conditions are met. A constraint may not be followed by a quantifier. The simple constraints are as follows; some more constraints are described later, under ESCAPES.

^
matches at the beginning of a line

$
matches at the end of a line

(?=re)
positive lookahead (AREs only), matches at any point where a substring matching re begins

(?!re)
negative lookahead (AREs only), matches at any point where no substring matching re begins

The lookahead constraints may not contain back references (see later), and all parentheses within them are considered non-capturing.

An RE may not end with `\\'.
"
    if {![winfo exists $w]} {
    	toplevel 	$w
    	wm title    $w $header
    	wm iconname $w $header

    	frame     $w.buttons
    	pack      $w.buttons    -side bottom  -fill x
    	button    $w.buttons.ok -text "OK"  -command "destroy $w"
    	pack      $w.buttons.ok -side left  -expand 1
    	scrollbar $w.scrolly    -command "$w.text yview"
    	scrollbar $w.scrollx    -orient horiz -command "$w.text xview"
    	pack      $w.scrolly    -side right   -fill y
    	pack      $w.scrollx    -side bottom  -fill x
    	text      $w.text       -height 40 -width  90 \
    	                        -yscrollcommand "$w.scrolly set"   \
    	                         -xscrollcommand "$w.scrollx set"   \
    	                         -wrap none                         \
    	                         -setgrid 1
    	pack      $w.text -side left  -expand 1  -fill both

        $w.text config -state normal
        $w.text delete 1.0 end

        $w.text insert end "$extra_txt"
        $w.text insert end "\n$re_help_txt"
        $w.text config -state disabled

    }
}
# RegExpHelpWindow
