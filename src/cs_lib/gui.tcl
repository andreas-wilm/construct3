##############################################################################
#
# gui.tcl - gui procedures
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
#  CVS $Id: gui.tcl,v 1.73 2007-10-26 13:51:13 steger Exp $
#


###   CleanGui   
#
# for reentrance
# safefail killing of all items
#
#
proc CleanGui {} {
    ############
    global w
    global c


    # kill aln.win.
    catch {destroy $w(al)}


    # kill seq labels
    catch {$c(seq_top)   delete DpSeqLabel}
    catch {$c(seq_left)  delete DpSeqLabel}
    catch {$c(seq_right) delete DpSeqLabel}
    catch {$c(seq_below) delete DpSeqLabel}


    # kill dotplot items separately
    catch {$c(dp) delete LAYER_ConsGaps}

    catch {$c(dp) delete DpSeqLabel}
    catch {$c(dp) delete DpNum}


    catch {$c(dp) delete AlNum} ;# numbering labels
    catch {$c(dp) delete AlNt} ;# nts
    catch {$c(dp) delete  Gap}

    catch {$c(dp) delete DpBp}
    catch {$c(dp) delete DpConsGap}
    catch {$c(dp) delete DpConsBp}


    catch {DeleteInfoDp}

}
# CleanGui




###   CreateDpWindow   
#
#   create the main dot plot window
#   including the sequence windows above, below, left, and right of the
#   window dot plot and the numbering window for the upper and left sequence
#
#
proc CreateDpWindow {} {
###################

    global c ;# w
    global w ;# w
    global SEQHEIGHT COLOR  FONT
    global lb_dp_pos ;# w
    global IMG

    set w(dp)  ".dotplot"
    set c(dp)  $w(dp).c_dp

    if {[winfo exists $w(dp)]} {
        destroy $w(dp)
    }
    toplevel    $w(dp)

    wm title    $w(dp) "ConStructDP"
    wm iconname $w(dp) "ConStructDP"
    wm withdraw .
    wm geometry $w(dp) +175+50
    wm protocol $w(dp) WM_DELETE_WINDOW "B4Exit $w(dp)"


    set   f_dp  $w(dp).f_dp
    frame $f_dp -relief raised -bd 2
    set   c(dp) $f_dp.c_dp
    grid  $f_dp -row 1 -column 0 -sticky nwse

    set   lb_dp_pos  $f_dp.lb_dp_pos_marker
    label $lb_dp_pos -bd 2 -width 10 -anchor w \
                     -font [format $FONT(dp_nm,format) $FONT(dp_nm,pt)]
    grid  $lb_dp_pos -column 3 -columnspan 2 -row 3 -sticky nsew

    # zoom button
    set zoom_balloon_text   "Zoom Dotplot in/out"
    button $f_dp.b_zoom -image $IMG(zoom) -relief raised
    grid   $f_dp.b_zoom -column 6 -columnspan 2 -row 3 -sticky nsew
    bind   $f_dp.b_zoom <1> "ZoomDotplot 2.0"
    bind   $f_dp.b_zoom <3> "ZoomDotplot 0.5"
    BalloonHelp::set_balloon $f_dp.b_zoom "$zoom_balloon_text"


    set c(dp_num_above)  $f_dp.c_dp_numbering_above
    set c(dp_num_left)   $f_dp.c_dp_numbering_left
    set c(seq_top)       $f_dp.c_seq_top
    set c(seq_below)     $f_dp.c_seq_below
    set c(seq_left)      $f_dp.c_seq_left
    set c(seq_right)     $f_dp.c_seq_right

    canvas $c(dp_num_above) -width "$c(size)c" \
                            -height ${SEQHEIGHT}c -borderwidth 2 \
                            -scrollregion "0c 0c $c(size) ${SEQHEIGHT}c"
    canvas $c(dp_num_left)  -width  ${SEQHEIGHT}c -borderwidth 2 \
                            -height "$c(size)c" \
                            -scrollregion "0c 0c $c(size) ${SEQHEIGHT}c"

    grid   $c(dp_num_above) -column 5 -row 3 -sticky wens
    grid   $c(dp_num_left)  -column 3 -row 5 -sticky wens

    canvas $c(seq_top) -width "$c(size)c" \
             -height ${SEQHEIGHT}c \
      	     -scrollregion "0c 0c $c(size) ${SEQHEIGHT}c" \
             -relief sunken -borderwidth 2
    canvas $c(seq_below) -width "$c(size)c" \
             -height ${SEQHEIGHT}c \
      	     -scrollregion "0c 0c $c(size) ${SEQHEIGHT}c" \
             -relief sunken -borderwidth 2
    canvas $c(seq_left) -width ${SEQHEIGHT}c \
             -height "$c(size)c" \
      	     -scrollregion "0c 0c $c(size) ${SEQHEIGHT}c" \
             -relief sunken -borderwidth 2
    canvas $c(seq_right) -width ${SEQHEIGHT}c \
             -height "0.0c" \
      	     -scrollregion "0c 0c $c(size) ${SEQHEIGHT}c" \
             -relief sunken -borderwidth 2

    grid $c(seq_top)   -column 5 -row 4 -sticky ew
    grid $c(seq_below) -column 5 -row 6 -sticky ew
    grid $c(seq_left)  -column 4 -row 5 -sticky ns
    grid $c(seq_right) -column 6 -row 5 -sticky ns

    scrollbar $f_dp.hscroll -orient horiz -command "DpScroll x"
    scrollbar $f_dp.vscroll -orient vert  -command "DpScroll y"
    grid      $f_dp.hscroll -column 5 -row 7 -sticky ew
    grid      $f_dp.vscroll -column 7 -row 5 -sticky ns

    canvas $c(dp) -width "0.0c" -height "0.0c" \
      	      -scrollregion "0c 0c $c(size) $c(size)" \
      	      -xscrollcommand "DpDrag $f_dp.hscroll x" \
      	      -yscrollcommand "DpDrag $f_dp.vscroll y" \
      	      -relief sunken \
      	      -borderwidth 0

    grid $c(dp) -column 5 -row 5 -sticky nwse

    $c(dp) create rectangle 0.0c 0.0c $c(size)c $c(size)c \
      	        -fill $COLOR(CDpBg) -tags "DPBorderRect"

    $c(dp) create line 0.0c 0.0c "$c(size)c" "$c(size)c"

    grid rowconfigure    $w(dp) 1 -weight 1	      ;# allow adjusting to window size
    grid columnconfigure $w(dp) 0 -weight 1
    grid rowconfigure    $f_dp  5 -weight 1
    grid columnconfigure $f_dp  5 -weight 1


    $c(dp) create rectangle 0c 0c -1c -1c -tags MirroredPosRectangle
}
# CreateDpWindow




###   CreateAlWindow   
#
# Creates the alignment window
#
#
proc CreateAlWindow {nseq aln_len} {
    ##############################
    global c
    global w
    global scale
    global c3
   #c3 ist separat scrollbares alignmentfenster fuer das 3' Ende#
    set w(al)    .alignment
    set p    $w(al).paned 
    set l    $w(al).paned.links 
    set r    $w(al).paned.rechts
    set s    $w(al).paned.seqlab
    set c(al)           $l.c_al
    set c3(al)          $r.c3_al
    set c(al_sn)        $s.c_al_sname
    set c(al_nm)        $l.c_al_numbering
    set c3(al_nm)       $r.c3_al_numbering
    set c(al_nm_dummy)  $s.c_al_numb_dummy

    if {[winfo exists $w(al)]} {
        destroy $w(al)
    }
    toplevel    $w(al)
    panedwindow $p  -orient horizontal -opaqueresize 1 -borderwidth 1 -sashpad 1 -showhandle false
    grid $p -column 0 -columnspan 3 -row 3 -rowspan 3   -sticky nwes 

    grid columnconfigure $w(al) 0 -weight 1
    grid columnconfigure $w(al) 1 -weight 1
    grid columnconfigure $w(al) 2 -weight 1 
    grid rowconfigure    $w(al) 3 -weight 1

    canvas $s; #-width 3c -height [expr {5*$scale(al)}]c 
    canvas $l; #-width  "$c(size)c" -height [expr {7*$scale(al)}]c 
    canvas $r; #-width  "$c(size)c" -height [expr {7*$scale(al)}]c

    grid $s -column 0 -row 2 -rowspan 3 -sticky nwes 
    grid $l -column 1 -row 2 -rowspan 3 -sticky nwes 
    grid $r -column 2 -row 2 -rowspan 3 -sticky nwes 

    wm title    $w(al)   "ConStruct alignment: $cs_proj::proj(name)"
    wm iconname $w(al)   "ConStruct alignment: $cs_proj::proj(name)"

    # FIXME: wm geometry $w(al)   $pos_w(align)
    wm protocol $w(al)   WM_DELETE_WINDOW "B4Exit $w(al)"
#	 
    scrollbar $l.hscroll     -orient horiz -command "AlScroll x"
    scrollbar $r.h3scroll    -orient horiz -command "DreiAlScroll x"
    scrollbar $s.hsnscroll   -orient horiz -command "SnScroll x"
    scrollbar $w(al).vscroll -orient vert  -command "AlScroll y" 
    grid      $w(al).vscroll -column 3 -row 3 -sticky nes 
 
 
    set scrx_min 0
    set scrx_max [expr {$scale(al)*$aln_len}]
    set scry_min 0
    set scry_max [expr {$scale(al)*$nseq}]
    set scrxsn_max [expr {$scale(al)*10}]
    set listSCR [list ${scrx_min}c ${scry_min}c ${scrx_max}c ${scry_max}c]

    canvas $c(al_sn)        -width  3c          -height [expr {5*$scale(al)}]c                        -scrollregion [list ${scrx_min}c ${scry_min}c ${scrxsn_max}c ${scry_max}c] -xscrollcommand "$s.hsnscroll set"   
    canvas $c(al_nm_dummy)  -width  3c          -height ${scale(al)}c 
    canvas $c(al_nm)        -width  $c(size)c   -height ${scale(al)}c   -relief sunken -borderwidth 1 -scrollregion $listSCR -xscrollcommand "$l.hscroll set"   -yscrollcommand "$w(al).vscroll set"
    canvas $c3(al_nm)       -width  "$c(size)c" -height ${scale(al)}c   -relief sunken -borderwidth 1 -scrollregion $listSCR -xscrollcommand "$r.h3scroll set"  -yscrollcommand "$w(al).vscroll set"
    canvas $c(al)           -width  "$c(size)c" -height [expr {5*$scale(al)}]c \
                                                                        -relief sunken -borderwidth 1 -scrollregion $listSCR -xscrollcommand "$l.hscroll set"   -yscrollcommand "$w(al).vscroll set" 
    canvas $c3(al)          -width  "$c(size)c" -height [expr {5*$scale(al)}]c \
                                                                        -relief sunken -borderwidth 1 -scrollregion $listSCR -xscrollcommand "$r.h3scroll set"  -yscrollcommand "$w(al).vscroll set"

    $p add  $s -minsize 30				  
    $p add  $l -minsize 40
    $p add  $r -minsize 40	

    pack  $l.hscroll      -fill both -side bottom
    pack  $r.h3scroll     -fill both -side bottom
    pack  $s.hsnscroll    -fill both -side bottom
     	 
    pack  $c(al_nm_dummy) -fill both -side top
    pack  $c(al_nm)       -fill both -side top
    pack  $c3(al_nm)      -fill both -side top
    pack  $c(al_sn)       -fill both -expand true
    pack  $c(al)          -fill both -expand true
    pack  $c3(al)         -fill both -expand true
}
# CreateAlWindow





###   CreateDpMenu   
#
#
#
proc CreateDpMenu {} {
    ################
    global w            ;# r
    global MENUINDEX    ;# w
    global opts         ;# r

    set column 0

    set    f_menu $w(dp).f_menu
    frame $f_menu -relief raised -bd 2
    grid  $f_menu -row 0 -sticky ew


    #####   File
    #
    #
    set midx 1
    set MENUINDEX(File,OpenProject)	   $midx;incr midx
    set MENUINDEX(File,SaveAlignment)  $midx;incr midx
    incr midx;# separator
    set MENUINDEX(File,PrintDp)        $midx;incr midx
    incr midx;# separator
    set MENUINDEX(File,Exit)           $midx;incr midx
    #
    set m $f_menu.file.menu
    menubutton $f_menu.file  -text "File"  -menu $m
    menu $m
    $m add command  -label "Open Project"   -command "ChooseProject"
    $m add command  -label "Save Alignment" -command "SaveAln"
    #
    $m add separator
    #
    $m add cascade  -label "Print Dotplot" -menu $m.printdp
    #
    $m add cascade  -label "Print Alignment" -menu $m.printaln
    #
    $m add separator
    #
    $m add command  -label "Exit"          -command "B4Exit $w(dp)"

    grid $f_menu.file   -column $column -row 0 -sticky w


    ###  submenu (cascade) printdp
    #
    menu $m.printdp
    $m.printdp add command -label "to printer"        -command "PrintDotplot print"
    $m.printdp add command -label "to color printer"  -command "PrintDotplot printcolor"
    $m.printdp add command -label "to file"           -command "PrintDotplot save"
    $m.printdp add command -label "to screen"         -command "PrintDotplot screen"
    $m.printdp add cascade -label "Options" -menu $m.dpp_options

    ###  subsubmenu (cascade) dpp_options
    #
    menu $m.dpp_options
    $m.dpp_options add checkbutton -label "Remove sequence selection" \
    	-variable opts(print_dp,rem_selection)
    $m.dpp_options add checkbutton -label "Add project" \
    	-variable opts(print_dp,add_project)
    $m.dpp_options add checkbutton -label "Add sequence" \
    	-variable opts(print_dp,add_sequence)
    $m.dpp_options add checkbutton -label "Add numbering" \
    	-variable opts(print_dp,add_numbering)

    ###  submenu (cascade) printaln
    #
    menu $m.printaln
    $m.printaln add command -label "to printer"        -command "PrintAlignment print"
    $m.printaln add command -label "to color printer"  -command "PrintAlignment printcolor"
    $m.printaln add command -label "to file"           -command "PrintAlignment save"
    $m.printaln add command -label "to screen"         -command "PrintAlignment screen"
#        $m.printaln add cascade -label "Options" -menu $m.dpp_options
    #
#            ###  subsubmenu (cascade) dpp_options
#            
#            menu $m.dpp_options
#            $m.dpp_options add checkbutton -label "Remove sequence selection" \
#                                           -variable opts(print_dp,rem_selection)
#            $m.dpp_options add checkbutton -label "Add project" \
#                                           -variable opts(print_dp,add_project)
#            $m.dpp_options add checkbutton -label "Add sequence" \
#                                           -variable opts(print_dp,add_sequence)
#            $m.dpp_options add checkbutton -label "Add numbering" \
#                                           -variable opts(print_dp,add_numbering)
#

    incr column



    #####   Structure
    #
    #
    set m $f_menu.structure.menu
    menubutton $f_menu.structure  -text "Structure_Prediction"  -menu $m
    menu $m
    #
    $m add cascade  -label "Optimal"    -menu $m.opt
    $m add cascade  -label "Suboptimal" -menu $m.subopt
    $m add separator
    $m add cascade  -label "Tertiary"   -menu $m.tertiary

    grid $f_menu.structure   -column $column -row 0 -sticky w


    ###  submenu (cascade) optimal
    #
    menu $m.opt
    #
    $m.opt add command -label   "Draw Structure" \
    	-command "OptimalConStruct drawstruct"
    $m.opt add command -label   "Circles" \
    	-command "OptimalConStruct circles"
    $m.opt add command -label   "Structural Alignment" \
    	-command "OptimalConStruct struct_aln"
    $m.opt add command -label   "Request Structure Logo" \
    	-command "OptimalConStruct struct_logo"


    ###  submenu (cascade) suboptimal
    #
    menu $m.subopt
    #
    $m.subopt add command -label    "Draw Structure" \
    	-command  "SuboptimalConStruct drawstruct"
    $m.subopt add command -label    "Circles" \
    	-command  "SuboptimalConStruct circles"
    $m.subopt add command -label    "Structural Alignment" \
    	-command  "SuboptimalConStruct struct_aln"
    $m.subopt add command -label   "Request Structure Logo" \
    	-command "SuboptimalConStruct struct_logo"



    ###  submenu (cascade) tertiary
    #
    menu $m.tertiary
    #
    $m.tertiary add cascade -label "Pseudoknots (Imatch)" -menu $m.tertiary.imatch
    $m.tertiary add cascade -label "Basetriples (Bmatch)" -menu $m.tertiary.bmatch
    menu $m.tertiary.imatch
    $m.tertiary.imatch add command -label "Circles"                -command "Imatch circles"
    $m.tertiary.imatch add command -label "Structural Alignment"   -command "Imatch struct_aln"
    $m.tertiary.imatch add command -label "Request Structure Logo" -command "Imatch struct_logo"
    menu $m.tertiary.bmatch
    $m.tertiary.bmatch add command -label "Circles"                -command "Bmatch circles"
    $m.tertiary.bmatch add command -label "Structural Alignment"   -command "Bmatch struct_aln"
    $m.tertiary.bmatch add command -label "Request Structure Logo" -command "Bmatch struct_logo"



    incr column




    ### Calculation Base
    #
    # FIXME: some unused
    #
    set midx 1
    incr midx;# title
    set MENUINDEX(CalcBase,ComputeMic)      $midx;incr midx
    set MENUINDEX(CalcBase,HideMic)         $midx;incr midx
    set MENUINDEX(CalcBase,Colormapping)    $midx;incr midx
    incr midx;# separator
    set MENUINDEX(CalcBase,use_alifoldscore) $midx;incr midx
    set MENUINDEX(CalcBase,use_mic)          $midx;incr midx
    incr midx;# separator
    incr midx;# title    
    set MENUINDEX(CalcBase,use_stacking)     $midx;incr midx
    incr midx;# separator
    incr midx;# title    
    set MENUINDEX(CalcBase,LogE)            $midx;incr midx
    set MENUINDEX(CalcBase,Log2)            $midx;incr midx
    set MENUINDEX(CalcBase,Unbiased)        $midx;incr midx
    set MENUINDEX(CalcBase,MLM)             $midx;incr midx
    set MENUINDEX(CalcBase,PairEntropyNorm) $midx;incr midx
    incr midx;# separator
    incr midx;# dummy
    set MENUINDEX(CalcBase,ThreshTd)       $midx;incr midx
    set MENUINDEX(CalcBase,ThreshMic)      $midx;incr midx
    set MENUINDEX(CalcBase,Function)       $midx;incr midx
    incr midx;# separator
    set MENUINDEX(CalcBase,Output)         $midx;incr midx


    set m $f_menu.calc_base.menu
    menubutton $f_menu.calc_base  -text "Calculation_Base"  -menu $m
    menu $m
    #
    #
    $m add command -label "Covariation:"
    $m add command -label "  Compute/Show Covariation" -command "MutualInfoFrontend"
    $m add command -label "  Hide Covariation" -command "DeleteInfoDp"
    $m add command -label "  Treshold (Colormapping)" -command "MIC_ShowColorMapWin"
    #
    $m add separator
    $m add radiobutton -label "  Use RNAalifold score"  \
        -variable opts(mic,use_alifoldscore) -value 1 -command "MIC_SetOpts use_alifoldscore 1"
    $m add radiobutton -label "  Use Mutual Information (MI)" \
        -variable opts(mic,use_alifoldscore) -value 0 -command "MIC_SetOpts use_alifoldscore 0"
    #
    $m add separator
    $m add command -label "RNAalifold Options:"
    $m add checkbutton -label "  Use Stacking for RNAalifold" \
        -variable opts(mic,use_stacking) -onvalue 1 -offvalue 0 \
        -command "MIC_SetOpts opts(mic,use_stacking) $opts(mic,use_stacking)"
    #
    $m add separator
    #
    $m add command -label "MI Options:"
    $m add radiobutton -label "  Use log_e for prob. estimation" \
        -variable opts(mic,bit) -value 0 -command "MIC_SetOpts bit 0"
    $m add radiobutton -label "  Use log_2 for prob. estimation"  \
        -variable opts(mic,bit) -value 1 -command "MIC_SetOpts bit 1"
    $m add radiobutton -label "  Use unbiased prob. estimation" \
        -variable opts(mic,unbiased) -value 1 -command "MIC_SetOpts unbiased 1"
    $m add radiobutton -label "  Use max likelihood estimation" \
        -variable opts(mic,unbiased) -value 0 -command "MIC_SetOpts unbiased 0"
    $m add checkbutton -label "  Pair-Entropy Normalization" \
        -variable opts(mic,pair_entropy_norm) -onvalue 1 -offvalue 0 \
        -command "MIC_SetOpts opts(mic,pair_entropy_norm) $opts(mic,pair_entropy_norm)"
    
    #
    $m add separator
    #
    $m add command -label   "Combined Dotplot:"
    $m add command -label   [format "  TD Threshold = %-2.3f" $opts(td,threshold)] \
                   -command "ChangeProbThreshold $m"
    $m add command -label   [format "  MI Threshold = %-2.3f" $opts(mic,Colormap,0)] \
                   -command "MIC_ShowColorMapWin"
    $m add command -label   [format "  Function = %-2.3f*TD + %-2.3f*MIC" $opts(td,factor) $opts(mic,factor)] \
                   -command "ChangeCalcBaseFactors $m"
    #
    $m add separator
    #
    $m add command -label   "Output to console" -command "OutputCalcBaseOptions"


    grid $f_menu.calc_base   -column $column -row 0 -sticky w

    incr column



    #####   Alignment
    #
    #
    set m $f_menu.alignment.menu
    menubutton $f_menu.alignment  -text "Alignment"  -menu $m
    menu $m

    $m add command  -label "Seq. Search"     -command ShowSeqSearchWin
    $m add separator
    $m add command  -label "Map Nt to Helix" -command {AlnMapNtToHelix show}
    $m add command  -label "Clear"           -command {AlnMapNtToHelix clear}
    #$m add separator
    # AW: need structure for it $m add command  -label "Request Logo" -command {StructLogoFrontend}
### KL
    $m add separator
    $m add command  -label "Minialignment"  -command {MiniAl} -state disabled
    # AW: tmp disabled for paper release
### KL
    #
    grid $f_menu.alignment   -column $column -row 0 -sticky w

    incr column




    #####  Options
    #
    #
    set midx 1
    set MENUINDEX(Options,UseConSeq)       $midx;incr midx
    set MENUINDEX(Options,SingleBp)        $midx;incr midx
    set MENUINDEX(Options,RemoveGaps)      $midx;incr midx
    incr midx; # separator
    set MENUINDEX(Options,Highlight)       $midx;incr midx
    #
    set MENUINDEX(Options,ClearMIC)        $midx;incr midx
    incr midx; # separator
    #
    set MENUINDEX(Options,ConsBp)          $midx;incr midx
    set MENUINDEX(Options,Rectangle)       $midx;incr midx
    incr midx; # separator
    set MENUINDEX(Options,SuboptStructNum) $midx;incr midx
    incr midx; # separator
    set MENUINDEX(Options,UseMapForStructAln) $midx;incr midx
    set MENUINDEX(Options,UseMapInDotplot)    $midx;incr midx

    #
    set m $f_menu.options.menu
    menubutton $f_menu.options  -text "Options"  -menu $m
    menu $m
    #
    $m add checkbutton -label "Use Cons.Seq. for Struct.Prediction" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,use_cons_seq)
    $m add checkbutton -label "Remove Gaps for DrawStructure" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,remove_gaps)
    $m add checkbutton -label "Allow single basepairs"  -onvalue 1 -offvalue 0 \
                       -variable opts(struct,allow_single_bp)
    #
    $m add separator
    #
    $m add checkbutton -label "Show Consensus Basepairs" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,cons_bp) \
                       -command {DpShowHide DpConsBp $opts(displ,cons_bp)}
    $m add checkbutton -label "Show Gaps" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,gaps) \
                       -command {DpShowHide DpConsGap $opts(displ,gaps)}
    $m add checkbutton -label "Show Basepairs" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,bps) \
                       -command {DpShowHide DpBp $opts(displ,bps)}
    $m add checkbutton -label "Show mirrored rectangle" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,show_rectangle)
    $m add checkbutton -label "Highlight Consensus-Nts" -onvalue 1 -offvalue 0 \
                       -variable opts(displ,highlight_nt_from_csbp)

    #
    $m add separator
    #
    set txt [format "Number of suboptimal structures = %d" $opts(suboptstruct,number)]
    $m add command -label $txt -command "ChangeSuboptStructNum $m"

    #
    $m add separator
    #
    $m add checkbutton -label "Show Mapping Info in Struct.Aln." \
    	-variable opts(mapping,use_in_structaln) -onvalue 1 -offvalue 0 \
    	-command {set StructAln::displ_mapviols $opts(mapping,use_in_structaln)}
    $m add checkbutton -label "Use Mapping Info in Dotplot" \
    	-variable opts(mapping,use_in_dotplot) -onvalue 1 -offvalue 0 \
    	-command {DotplotMappingFrontend}
    #
    $m add separator
    #
    $m add checkbutton -label "Show Seq. Stats in Struct.Aln." \
    	-variable opts(structaln,show_seq_stat) -onvalue 1 -offvalue 0
    $m add checkbutton -label "Show Struct. Stats in Struct.Aln." \
    	-variable opts(structaln,show_struct_stat) -onvalue 1 -offvalue 0
    $m add checkbutton -label "Show Pattern Stats in Struct.Aln." \
    	-variable opts(structaln,show_patt_stat) -onvalue 1 -offvalue 0



    grid $f_menu.options   -column $column -row 0 -sticky e


    incr column



    #####    Debug
    #
    #
    if {$opts(debug)} {
        set m $f_menu.debug.menu
        menubutton $f_menu.debug  -text "Debug"  -menu $m
        menu $m
        $m add command  -label "Dump Vars (Tcl)"   -command "DumpFrontend vars"
        $m add command  -label "Dump Tags"         -command "DumpFrontend tags"

        $m add checkbutton -label "Online Debugging" -onvalue 0 -offvalue 1 \
                           -variable DEBUGGER_SKIP_ALL

        $m add cascade  -label "Print Matrix" -menu $m.print_mat


        grid $f_menu.debug   -column $column -row 0 -sticky e
        grid columnconfigure $f_menu $column -weight 1

            # submenu print matrix
            menu $m.print_mat
            $m.print_mat add command  -label "Print Cons.TD-Matrix" \
                                      -command "PrintDebugMat cons_td"
            $m.print_mat add command  -label "Print TD-Matrix for selected sequence" \
                                      -command "PrintDebugMat seq_td"
            $m.print_mat add command  -label "Print MIC" \
                                      -command "PrintDebugMat mic"
            $m.print_mat add command  -label "Print Combined Dp-Matrix" \
                                      -command "PrintDebugMat cons_dpm"

        incr column
    }




    #####   Help
    #
    #
    set MENUINDEX(Help,Help)    1
    set MENUINDEX(Help,Version)  2
    #
    set m $f_menu.help.menu
    menubutton $f_menu.help  -text "Help"  -menu $m
    menu $m

    $m add command  -label "Manual"    -command {ShowHelpHint}
    $m add command  -label "Version" -command {ShowVersion}

    grid $f_menu.help   -column $column -row 0 -sticky e
    grid columnconfigure $f_menu $column -weight 1


    incr column
}
# CreateDpMenu





###   CreateAlnNtsAndLabels   
#
#
#
proc CreateAlnNtsAndLabels {} {
    #########################

    global scale
    global FONT
    global selection
    global COLOR
    global selection
    global c
    global c3
    global seq ;# r


    ### Setup color
    #
    for {set s 1} {$s<=$seq(n_seq)} {incr s} {
        if {$s == $selection(seq_no)} {
            set col($s) $COLOR(SelSeq)
        } else {
            set col($s) $COLOR(UnselSeq)
        }
    }


    ### Create numbering labels
    #
    for {set i 10} {$i<=$seq(aln_len)} {incr i +10} {
        set coords [GetNtAliCoords 1 $i]
        set x1 [lindex $coords 0]
        set y1 [lindex $coords 1]
        $c(al_nm) create text "${x1}c" "${y1}c" -tags AlNum     \
                        -text $i  -font $FONT(al_sn)  -anchor center
        $c3(al_nm) create text "${x1}c" "${y1}c" -tags AlNum    \
                        -text $i  -font $FONT(al_sn)  -anchor center		
    }


    ###  place nts and seq labels
    #
    for {set s 1} {$s<=$seq(n_seq)} {incr s} {

        for {set n 1} {$n<=$seq(aln_len)} {incr n} {

            set seqchar [string index $seq(nt,$s) [expr {$n-1}]]
            set coords [GetNtAliCoords $s $n]
            set x1 [lindex $coords 0]
            set y1 [lindex $coords 1]
            if { ! [IsGap $seqchar] } {
                $c(al) create text "${x1}c" "${y1}c" \
                	-tags "AlnNt  AlnNt_seq_${s}_nt_${n}" \
                    -anchor center -text $seqchar -font $FONT(al)
				$c3(al) create text "${x1}c" "${y1}c" \
                	-tags "AlnNt  AlnNt_seq_${s}_nt_${n}" \
                    -anchor center -text $seqchar -font $FONT(al)
            } else {
                $c(al) create text "${x1}c" "${y1}c" \
                	-tags "Gap Gap_seq_${s}_pos_${n}" \
                    -anchor center  -text $seqchar -font $FONT(al)
				$c3(al) create text "${x1}c" "${y1}c" \
                	-tags "Gap Gap_seq_${s}_pos_${n}" \
                    -anchor center  -text $seqchar -font $FONT(al)
            }
        }


        # Create sequence labels
        #
    	set balloonhelp "$cs_proj::comment($s)"
        label $c(al_sn).l_sn${s}  -relief raised  -bd 2               \
                                  -text $seq(id,$s) -font $FONT(al_sn)\
                                  -fg black -bg $col($s) -anchor w 
    	if {$balloonhelp!=""} {
    		BalloonHelp::set_balloon $c(al_sn).l_sn${s} "$balloonhelp"
    	}
        if {$col($s)==$COLOR(SelSeq)} {
            $c(al_sn).l_sn${s} config -fg white
        }
        $c(al_sn) create window 0c [expr {$scale(al)*($s)-$scale(al)/2}]c \
                          -width 5c -anchor w -window $c(al_sn).l_sn${s}
    }

    $c(al) yview

}
# CreateAlnNtsAndLabels





###   CreateLayers   
#
# Create some dummy items in dotplot
# all other items gaps, cons-gaps, bps, cons-bps
# are stacked relative to these
#
proc CreateLayers {} {
    ################
    global c     ;# r
    global seq   ;# r

    # create a dummy object foreach possible gap number
    #
    for {set i 0} {$i<=$seq(n_seq)} {incr i} {
        $c(dp) create line 0c -1.5c 0c -1.0c -fill white -width 1 \
                                          -tags "LAYER_ConsGaps_$i LAYER_ConsGaps"

    }

    $c(dp) create line 0c -1.5c 0c -1.0c -fill white -width 1 \
                                      -tags "LAYER_DpBp"
    $c(dp) create line 0c -1.5c 0c -1.0c -fill white -width 1 \
                                      -tags "LAYER_ConsBP"
}
# CreateLayers




###   CreateNewDpBp   
#
# called from core (knows the coords)
#
proc CreateNewDpBp {seq_idx nt_i nt_j x1 y1 x2 y2} {
    #############################################
    global c
    global COLOR

    set tags "DpBp DpBp_seq_${seq_idx} DpBp_seq_${seq_idx}_nti_${nt_i}_ntj_${nt_j}"

    set bp [$c(dp) create rectangle $x1 $y1 $x2 $y2 \
                            -activefill $COLOR(SelBp)  \
                            -fill $COLOR(UnselBp) \
                            -outline "" -tags $tags]
    return $bp
}
# CreateNewDpBp




###   CreateNewConsDpBp   
#
# called from core (knows the coords)
#
proc CreateNewConsDpBp {fill_color nt_i nt_j x1 y1 x2 y2} {
    #####################################################
    global c
    global COLOR

    set tags "DpConsBp DpConsBp_nti_${nt_i}_ntj_${nt_j}"

    $c(dp) create rectangle $x1 $y1 $x2 $y2 \
                            -fill $fill_color \
                            -outline "" -tags $tags

}
# CreateNewConsDpBp




###   RaiseDpBpsToLayer   
#
#
#
proc RaiseDpBpsToLayer {} {
    #####################
    global c

    $c(dp) raise DpBp LAYER_DpBp
}
# RaiseDpBpsToLayer




###   ResizeConsDp   
#
#
#
proc ResizeConsDp {nt_i nt_j x1 y1 x2 y2} {
##################
    global c

    set tag "DpConsBp_nti_${nt_i}_ntj_${nt_j}"
    $c(dp) coords $tag  $x1 $y1 $x2 $y2
}
# ResizeConsDp




###   DeleteConsBp   
#
#
#
proc DeleteConsBp {nt_i nt_j} {
    #########################
    global c

    set tag "DpConsBp_nti_${nt_i}_ntj_${nt_j}"
    $c(dp) delete $tag
}
# DeleteConsBp




###   CreateDummyConsGap   
#
# A dummy consensus gap
# This is for speedup:
# A move update doesn't need to create a new object an stack it correctly
# Instead it just recoords the dummy poylygon
#
# WARNING: 3 (not 2) coords in all unpatched tk versions (<8.4.2 <8.3.5?) for
#  postscript generation needed
# (BUG in tkCanvPoly.c:PolygonToPostscript, see
#  http://sourceforge.net/tracker/index.php?func=detail&aid=734498&group_id=12997&atid=112997
#  http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=192302
# )
#
proc CreateDummyConsGap {nt_idx} {
    ############################
    global c

    set tags "DpConsGap DpConsGap_nt_${nt_idx}"
    $c(dp) create polygon 0 0 0 0 0 0 -tags $tags
}
# CreateDummyConsGap




###   CreateConsGap   
#
# A consensus gap
#
#
proc CreateConsGap {nt_idx fill_color coord_list} {
    #############################################
    global c

    set tags "DpConsGap DpConsGap_nt_${nt_idx}"
    $c(dp) create polygon $coord_list -fill $fill_color -outline "" -tags $tags
}
# CreateConsGap




###   RecoordConsGap   
#
# Update a consensus gap after move
# if coord_list is omitted, the gap is converted to
# a dummy  (see CreateDummyConsGap)
#
proc RecoordConsGap {nt_idx {coord_list "DUMMY"}} {
    #############################################
    global c

    set tag "DpConsGap_nt_${nt_idx}"
    if {$coord_list=="DUMMY"} {
        $c(dp) coords $tag 0 0 0 0 0 0
    } else {
        $c(dp) coords $tag $coord_list
    }
}
# RecoordConsGap




###   RecolorConsGap   
#
# Update color of a consensus gap after move
#
proc RecolorConsGap {nt_idx fill_color} {
    ###################################
    global c

    set tag "DpConsGap_nt_${nt_idx}"

    $c(dp) itemconfigure $tag -fill $fill_color
}
# RecolorConsGap




###   RecolorConsDp   
#
#
#
proc RecolorConsDp {nt_i nt_j fill_color} {
    #####################################
    global c

    set tag "DpConsBp_nti_${nt_i}_ntj_${nt_j}"
    $c(dp) itemconfigure $tag -fill $fill_color
}
# RecolorConsDp




###   CreateToolButtons   
#
#
#
proc CreateToolButtons {} {
    #####################
    global IMG
    global w


    set move_balloon_text "Move selected sequence slice left/right"
    set gapins_balloon_text "Move consecutive succession of nucleotides left/right"


    ##########   Alignment window first
    #
    #

    ### Setup frames in alignment window
    #
    set    a_tl $w(al).f_tl
    if {[winfo exists $a_tl]} {
        destroy $a_tl
    }

    frame $a_tl -bd 2 ;# tool buttons
    grid  $a_tl -in $w(al)  -column 1 -row 0 -sticky w


    ### Setup frames for images inside top frame
    #
    frame $a_tl.f_move   -bd 2
    frame $a_tl.f_gapi   -bd 2
    grid  $a_tl.f_move   -column 0 -row 0 -sticky ew
    grid  $a_tl.f_gapi   -column 1 -row 0 -sticky ew



    ### Display and bind the stuff
    #

    button $a_tl.f_move.b_move -image $IMG(move)
    bind   $a_tl.f_move.b_move <1> {
        SetCursor busy
        RequestToMove -1
        SetCursor normal
    }
    bind   $a_tl.f_move.b_move <3> {
        SetCursor busy
        RequestToMove 1
        SetCursor normal
    }
    BalloonHelp::set_balloon $a_tl.f_move.b_move "$move_balloon_text"

    
    button $a_tl.f_gapi.b_gapi -image $IMG(gapins)
    bind   $a_tl.f_gapi.b_gapi <1> {
        SetCursor busy
        EmulateGapInsertion -1
        SetCursor normal
    }
    bind   $a_tl.f_gapi.b_gapi <3> {
        SetCursor busy
        EmulateGapInsertion 1
        SetCursor normal
    } 
    BalloonHelp::set_balloon $a_tl.f_gapi.b_gapi "$gapins_balloon_text"

    grid $a_tl.f_move.b_move  -column 0 -row 1 -sticky ns
    grid $a_tl.f_gapi.b_gapi  -column 0 -row 1 -sticky ns

}
# CreateToolButtons




###   CreateDpSeqLabels   
#
#
# Generates the sequence labels positioned beside the dotplot
#
proc CreateDpSeqLabels {seq} {
    ########################
    global c     ;# r
    global FONT  ;# w
    global scale ;# r


    set FONT(dp_seq,pt) [expr {int([expr {$scale(dp)*100.0/2.54*10}])}]
    if {$FONT(dp_seq,pt) < 10} {
        set FONT(dp_seq,pt) 10
    }
    set FONT(dp_seq) [format $FONT(dp_seq,format) $FONT(dp_seq,pt)]

    for {set i 1} {$i<=[string length $seq]} {incr i +1} {
        set nt [string index $seq [expr {$i-1}]]

        set coords [GetDpSeqCoords $i horiz]

        $c(seq_top)   create text $coords \
                       -text "$nt"\
                       -font $FONT(dp_seq) \
                       -anchor center \
                       -tags "DpSeqLabel DpSeqLabel_ntidx_$i"

        $c(seq_below) create text $coords \
                       -text "$nt"\
                       -font $FONT(dp_seq) \
                       -anchor center \
                       -tags "DpSeqLabel DpSeqLabel_ntidx_$i"

        set coords [GetDpSeqCoords $i vert]

        $c(seq_left) create text $coords \
                       -text "$nt"\
                       -font $FONT(dp_seq) \
                       -anchor center \
                       -tags "DpSeqLabel DpSeqLabel_ntidx_$i"

        $c(seq_right) create text $coords \
                       -text "$nt"\
                       -font $FONT(dp_seq) \
                       -anchor center \
                       -tags "DpSeqLabel DpSeqLabel_ntidx_$i"

    }
}
# CreateDpSeqLabels




###   GenDpNumLabels   
#
# FIXME: allow calling me after zoom
#
#
proc GenDpNumLabels {} {
    ##################
    global c     ;# r
    global seq   ;# r
    global FONT  ;# w
    global SEQHEIGHT ;# r
    global scale ;# r


    set viewable_dp_area [expr {$c(size)/$c(virtual_size)}]
    set viewable_nts     [expr {$seq(aln_len)*$viewable_dp_area}]
    set tickwidth        [expr {int($viewable_nts/4)}]
    set digits           [string length [expr {int($viewable_nts)}]]
    set digits           [expr {$digits-1}]
    set nulls ""
    for {set i 1} {$i<=$digits} {incr i} {
        set nulls ${nulls}0
    }
    set nulls1 [string range $nulls 1 end]
    set l_rangetest [format "{5%s 2%s} {2%s 1%s} {1%s 5%s}" \
                            $nulls $nulls $nulls $nulls $nulls $nulls1]
    foreach test $l_rangetest {
        if {$viewable_nts > [lindex $test 0]} {
            set tickwidth [lindex $test 1]
            break
        }
    }
    #puts "TMP_DEBUG(GenDpNumLabels): viewable_nts=$viewable_nts"
    #puts "TMP_DEBUG(GenDpNumLabels): tickwidth=$tickwidth"


    $c(dp_num_above) delete DpNum
    $c(dp_num_left)  delete DpNum

    set FONT(dp_nm) [format $FONT(dp_nm,format) $FONT(dp_nm,pt)]

    # numbering labels
    #
    for {set i $tickwidth} {$i<=$seq(aln_len)} {incr i +$tickwidth} {
        set  NtDpX    [expr {($i-0.5)*$scale(dp)}]
        set  NtDpY    [expr {$SEQHEIGHT/2}]
        set  C_NtDp_i [list [expr {${NtDpX}*$c(virtual_size)/$c(size)}]c ${NtDpY}c]
        set  C_NtDp_j [list ${NtDpY}c [expr {${NtDpX}*$c(virtual_size)/$c(size)}]c]

        # numbering labels i
        $c(dp_num_above) create text $C_NtDp_i \
                -text $i -font $FONT(dp_nm) \
                -anchor center -tags "DpNum"

        # numbering labels j
        $c(dp_num_left) create text $C_NtDp_j \
                -text $i -font $FONT(dp_nm) \
                -anchor center -tags "DpNum"
    }
}
# GenDpNumLabels




###   DeleteInfoDp   
#
#
proc DeleteInfoDp {} {
    #################
    global c ;# r

    # FIXME: hide first and delete only when options have changed
    foreach ip [$c(dp) find withtag MutInfoPair] {
        $c(dp) delete $ip
    }
}
# DeleteInfoDp




###   MutualInfoFrontend   
#
#
#
proc MutualInfoFrontend {} {
    #######################
    global c           ;# r
    global opts        ;# r
    global mic_is_computed_and_valid ;# r


    set mi_tag "MutInfoPair"
    SetCursor busy

    DeleteInfoDp

    if { ! $mic_is_computed_and_valid } {
        if {$opts(mic,use_alifoldscore)==1} {
            p::verbose "Computing AlifoldScore"
        } else {
            p::verbose "Computing MutualInfoContent"
        }

        ComputeMutualInfo $opts(mic,unbiased)           \
                          $opts(mic,bit)                \
                          $opts(mic,pair_entropy_norm)  \
                          $opts(mic,use_alifoldscore)   \
                          $opts(mic,use_stacking)

        set mic_is_computed_and_valid 1
    }

    CreateInfoDp $mi_tag $opts(mic,Colormap,0) $opts(mic,Colormap,1)

    update

    SetCursor normal
}
#




###   GetNtAliCoords   
#
# returns dimension in c without unit !
# -> append c if don't need to do expressions
#
proc GetNtAliCoords {seq_no nt_no} {
    ##############################
    global scale

    set NTAliX  [expr {($nt_no-0.5)*$scale(al)}]
    set NTAliY  [expr {$scale(al)*$seq_no - $scale(al)/2}]

    return [list ${NTAliX} ${NTAliY}]
}
# GetNtAliCoords



###   GetNtAliCellDim   
#
# width_or_height = width|height
# returns dimension in c without unit !
# -> append c if don't need to do expressions
#
proc GetNtAliCellDim {width_or_height} {
###################
    global scale

    set nt_no  1;# dummy
    set seq_no 1;# dummy

    if {$width_or_height=="width"} {
        set x1  [expr {($nt_no-0.5)*$scale(al)}]
        incr nt_no
        set x2  [expr {($nt_no-0.5)*$scale(al)}]
        return [expr {$x2-$x1}]
    } elseif {$width_or_height=="height"} {
        set y1  [expr {$scale(al)*$seq_no - $scale(al)/2}]
        incr seq_no
        set y2  [expr {$scale(al)*$seq_no - $scale(al)/2}]
        return [expr {$y2-$y1}]
    } else {
        p::error "invalid arg \"$width_or_height\", must be \"width\" or \"height\""
        return 0.0
    }
}
# GetNtAliCellDim



###   GetDpSeqCoords   
#
# nt_no   = nt number
# oritent = horiz|vert
#
proc GetDpSeqCoords {nt_no orient} {
    ##############################
    global scale
    global SEQHEIGHT

    set dp_x [expr {($nt_no-0.5)*$scale(dp)}]
    set dp_y [expr {$SEQHEIGHT/2}]

    if {$orient=="horiz"} {
        return [list ${dp_x}c ${dp_y}c]
    } elseif {$orient=="vert"} {
        return [list ${dp_y}c ${dp_x}c]
    } else {
        p::error "orient must be one of horiz|vert"
        return [list 0 0]
    }

}
# GetDpSeqCoords




###   ChangeMenuesState   
#
# Change cs_dp menus (dependent on successfull project loading)
#
proc ChangeMenuesState {state} {
    ##########################
    global w
    global MENUINDEX


    # all: structure prediction
    set menub  $w(dp).f_menu.structure
    $menub configure -state $state

    # all: alignment
    set menub  $w(dp).f_menu.alignment
    $menub configure -state $state

    # mic computation
    set m $w(dp).f_menu.calc_base.menu
    $m entryconfigure $MENUINDEX(CalcBase,ComputeMic) -state $state

    # colormapping
    set m $w(dp).f_menu.calc_base.menu
    $m entryconfigure $MENUINDEX(CalcBase,Colormapping) -state $state
    $m entryconfigure $MENUINDEX(CalcBase,ThreshMic)    -state $state

    # save alignment
    set m $w(dp).f_menu.file.menu
    $m entryconfigure $MENUINDEX(File,SaveAlignment) -state $state

    # print dotplot
    set m $w(dp).f_menu.file.menu
    $m entryconfigure $MENUINDEX(File,PrintDp) -state $state

}
# ChangeMenuesState





###   InitImages   
#
# see
# http://aspn.activestate.com/ASPN/Cookbook/Tcl/Recipe/117247
# http://mini.net/tcl/2897
#
# proc inlineGIF {img {name ""}} {
#     package require base64
#     set f [open $img]
#     fconfigure $f -translation binary
#     set data [base64::encode [read $f]]
#     close $f
#     if {[llength [info level 0]] == 2} {
# 	    # base name on root name of the image file
# 	    set name [file root [file tail $img]]
#     }
#     return "image create photo [list $name] -data {\n$data\n}"
# }
#
#
proc InitImages {} {
    ##############
    global IMG;# r

    set IMG(construct) [image create photo -data "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"]


    set IMG(zoom)      [image create photo -data "R0lGODdhHAAcAPIAAAAAAAAAgACAgICAgMDAwAAA/wD//////ywAAAAAHAAcAAADjEi63P4QyAmgVfLofeR1wEYBglZ92CF4jCSsaNhZJEyr5wa+5xPmmogNhFMEhcPWzPgh9ZTPZnLhRFGnGCxNW7USuo1aFDnO9gYNtHk7HLjdV238Cc+uyi1e1HUP4ON3Ey93Bgd+TTWDegABBYUFf1QjPYyOhixeYY2FAJCZIJ0/nz4ABpFeoaM0qpkJADs="]
    set IMG(move)      [image create photo -data "R0lGODdhHAAcAPEAAMDAwP///wAAAICAgCwAAAAAHAAcAAACVYSPqcvtD6OctNo7A9M0CCZwkfctQggFZKmcoZio7qnW9fwiMs73rmrwCYXBoXFmGAyOQyVCiVNKpdEBA3oCnZyPJSuxnSxBVkpZccao1+y2+w2XFAAAOw=="]
    set IMG(gapins)    [image create photo -data "R0lGODdhHAAcAPEAAMDAwP///wAAAPoAACwAAAAAHAAcAAACcISPqcvtG2IwcoL6gNjC8K5x2UgyUfkMA+qorCKp0muotk3X95rrbu/jAQG/IVFoTF4oxkEE+TgdpD8pxTLdWALa46X7FWA1njLZBUqbzeoTOpKGIz7ljZcesh/UGjmavUYmGKIjyAeSESDTVKNEUgAAOw=="]
    set IMG(reset)     [image create photo -data "R0lGODdhHAAcAPEAAMDAwP///wAAAICAgCwAAAAAHAAcAAACZ4SPqcvtD6OcKFAXxGVZ79p5XxAKAommFlea7vs28EyLSo2bzTDk9MOD8YbEQSSo+yR6NqWh51QYo1GUIQXAaq0HlyUGdgFiXZLmtEJ/Sd1bUjwWw8e38ZqNOnObdHpYR8ZHNUjoUAAAOw=="]



}
# InitImages




###   PrintAlignment   
#
# where: printer|color|file|screen
#
proc PrintAlignment {where} {
###################
    global c

    set plotfile_ext  "_aln.ps"
    set projname      $cs_proj::proj(name)
    set psfile        "${projname}${plotfile_ext}"
    set extra_ps_args ""

p::fixme "implement: remove selection"
p::fixme "implement: view port"
p::fixme "implement: numbering"

    eval PrintPS $c(al) $psfile $where $extra_ps_args

}
# PrintAlignment




###   PrintDotplot   
#
# where: printer|color|file|screen
#
#
proc PrintDotplot {where} {
    #####################
    global c
    global SEQHEIGHT
    global FONT
    global selection
    global opts



    set plotfile_ext  "_dp.ps"
    set projname      $cs_proj::proj(name)
    set psfile        "${projname}${plotfile_ext}"
    set extra_ps_args ""




    # actual viewport
    set view(x1) [expr {$c(virtual_size)*[lindex [$c(dp) xview] 0]}]
    set view(y1) [expr {$c(virtual_size)*[lindex [$c(dp) yview] 0]}]
    set view(x2) [expr {$c(virtual_size)*[lindex [$c(dp) xview] 1]}]
    set view(y2) [expr {$c(virtual_size)*[lindex [$c(dp) yview] 1]}]


    # coords for searching horizontally
    set find(h,x1) $view(x1)
    set find(h,y1) 0
    set find(h,x2) $view(x2)
    set find(h,y2) $SEQHEIGHT

    # coords for searching vertically
    set find(v,x1) 0
    set find(v,y1) $view(y1)
    set find(v,x2) $SEQHEIGHT
    set find(v,y2) $view(y2)

    set actualseqheight [expr {$FONT(dp_seq,pt)*2.54/100.0/10.0}]
    set actuallabheight [expr {$FONT(dp_nm,pt)*2.54/100.0/10.0}]
    set moveseqwidth1   [expr {($SEQHEIGHT/2)+$actualseqheight}]
    set moveseqwidth2   [expr {($SEQHEIGHT/2)-$actualseqheight}]
    set movelabwidth1   [expr {$actualseqheight+$actuallabheight}]
    set ps_offset       [expr {2*$actualseqheight}]
    set lab_size(o)     [expr {2*$actuallabheight}]
    set lab_size(l)     [expr {2*$actuallabheight}]



    ### remove selection
    #
    if {$opts(print_dp,rem_selection)} {
        set current_selection $selection(seq_no)
        SeqSelectionChange "DESELECT"
    }

    ### print numbering
    #
    #
    if {$opts(print_dp,add_numbering)} {
        foreach nt [$c(dp_num_above) find overlapping \
                                     ${find(h,x1)}c  ${find(h,y1)}c ${find(h,x2)}c ${find(h,y2)}c] {
            if {[lsearch -exact [$c(dp_num_above) gettags $nt] "DpNum"] > -1} {
                $c(dp) create text [$c(dp_num_above) coords $nt] \
                    -text [lindex [$c(dp_num_above) itemconfigure $nt -text] end] \
                    -font [lindex [$c(dp_num_above) itemconfigure $nt -font] end] \
                    -anchor center \
                    -tags "DpPrint DpNtsPrint DpNtsPrint_Numabove"
            }
        }
        foreach nt [$c(dp_num_left) find overlapping \
                                    ${find(v,x1)}c  ${find(v,y1)}c ${find(v,x2)}c ${find(v,y2)}c] {
    		if {[lsearch -exact [$c(dp_num_left) gettags $nt] "DpNum"] > -1} {
                $c(dp) create text [$c(dp_num_left) coords $nt] \
                    -text [lindex [$c(dp_num_left) itemconfigure $nt -text] end] \
                    -font [lindex [$c(dp_num_left) itemconfigure $nt -font] end] \
                    -anchor center \
                    -tags "DpPrint DpNtsPrint DpNtsPrint_Numleft"
    		}
        }
        # move the stuff to the right position
        $c(dp) move DpNtsPrint_Numabove 0c [expr {$view(y1)-$moveseqwidth1-$movelabwidth1}]c
        $c(dp) move DpNtsPrint_Numleft  [expr {$view(x1)-$moveseqwidth1-$movelabwidth1}]c 0c
    }


    ### print sequence
    #
    if {$opts(print_dp,add_sequence)} {

        # update labels to consensus sequence
        if {$opts(displ,use_cons_seq)} {
            UpdateDpSeqLabels [Get_ConsSeq]
        }

        foreach seq "seq_top seq_below" {
            foreach nt [$c($seq) find overlapping \
                                 ${find(h,x1)}c  ${find(h,y1)}c ${find(h,x2)}c ${find(h,y2)}c] {
                if {[lsearch -exact [$c($seq) gettags $nt] "DpSeqLabel"] > -1} {
                    $c(dp) create text [$c($seq) coords $nt] \
                        -text [lindex [$c($seq) itemconfigure $nt -text] end] \
       				    -font [lindex [$c($seq) itemconfigure $nt -font] end] \
       				    -anchor center \
       				    -tags "DpPrint DpNtsPrint DpNtsPrint_$seq"
                }
            }
        }

        foreach seq "seq_left seq_right" {
            foreach nt [$c($seq) find overlapping \
                                 ${find(v,x1)}c  ${find(v,y1)}c ${find(v,x2)}c ${find(v,y2)}c] {
                if {[lsearch -exact [$c($seq) gettags $nt] "DpSeqLabel"] > -1} {
                    $c(dp) create text [$c($seq) coords $nt] \
                        -text [lindex [$c($seq) itemconfigure $nt -text] end] \
       				    -font [lindex [$c($seq) itemconfigure $nt -font] end] \
       				    -anchor center \
       				    -tags "DpPrint DpNtsPrint DpNtsPrint_$seq"
                }
            }
        }
        # move the stuff to the right position
       $c(dp) move DpNtsPrint_seq_top   0c [expr {$view(y1)-$moveseqwidth1}]c
       $c(dp) move DpNtsPrint_seq_below 0c [expr {$view(y2)-$moveseqwidth2}]c
       $c(dp) move DpNtsPrint_seq_left  [expr {$view(x1)-$moveseqwidth1}]c 0c
       $c(dp) move DpNtsPrint_seq_right [expr {$view(x2)-$moveseqwidth2}]c 0c

    }




    ### print proj
    #
    if  {$opts(print_dp,add_project)} {
        set txt "$projname\n[clock format [clock seconds] -format "%b %d, %Y  %R"]"
        $c(dp) create text [expr {$view(x1)+0.5}]c [expr {$view(y2)-1.0}]c \
                           -text $txt -anchor w -tags "DpPrint DpPrintProjname"
    }


    # hide the basepair position rectangle
    $c(dp) lower MirroredPosRectangle


    set extra_ps_args    "-x [expr {$view(x1)-$ps_offset-$lab_size(l)}]c"
    append extra_ps_args " -width  [expr {$view(x2)-$view(x1)+2*$ps_offset+$lab_size(l)}]c"
    append extra_ps_args " -y [expr {$view(y1)-$ps_offset-$lab_size(o)}]c"
    append extra_ps_args " -height [expr {$view(y2)-$view(y1)+2*$ps_offset+$lab_size(o)}]c"



    eval PrintPS $c(dp) $psfile $where $extra_ps_args

    # restore selection
    if {$opts(print_dp,rem_selection)} {
        SeqSelectionChange $current_selection
    }


    # raise the basepair position rectangle
    $c(dp) raise MirroredPosRectangle

    # cleanup
    $c(dp) delete DpPrint
}
# PrintDotplot






###   AlnMapNtToHelix   
#
#
#
proc AlnMapNtToHelix {clear_or_show} {
    ################################
    global seq
    global optstruct
    global c

    if { ! $optstruct(loaded)} {
        p::error "no precomputed optimal structures loaded"
        return
    }
    if {$clear_or_show!="show" && $clear_or_show!="clear"} {
        p::error "wrong arg \"$clear_or_show\", must be clear or show"
        return
    }

    SetCursor busy


    for {set s 1} {$s<=$seq(n_seq)} {incr s} {

        if {$clear_or_show=="show"} {
            set struct [DotbracketToAlNum $optstruct(db,$s)]
        }

        set ngaps 0
        for {set n 1} {$n<=$seq(aln_len)} {incr n} {
            set nt [string index $seq(nt,$s) [expr {$n-1}]]
            if {$clear_or_show=="show"} {

                if {! [IsGap $nt]} {
                    set char [string index $struct [expr {$n-1-$ngaps}]]
                } else {
                    incr ngaps
                    continue
                }
            } else {
                set char $nt
            }
            $c(al) itemconfigure AlnNt_seq_${s}_nt_${n} -text $char
        }
    }

    SetCursor normal
}
# AlnMapNtToHelix



###   DotplotMappingFrontend
#
# reads opts(mapping,use_in_dotplot) and decides if
# mapping info should be merge in dotplot or removed
#
proc DotplotMappingFrontend {} {
    #######################
    global dp_hidden_objects ;# w
    global opts              ;# r

    if {$opts(mapping,use_in_dotplot)==0} {
    	RemoveMappingFromDotplot $dp_hidden_objects(mapping)
    	set dp_hidden_objects(mapping) {}

    } else {
    	set dp_hidden_objects(mapping) [MergeMappingInDotplot]
    }
}
# DotplotMappingFrontend



###   MergeMappingInDotplot
#
# Merges Mapping Info in dotplot, i.e.
# Remove amibious basepairs from dotplot where
#   nucleotides are mapped to be unpaired or
#   set to pair with one particular nt
#
proc MergeMappingInDotplot {} {
    #########################
    global seq               ;# r

    p::verbose "Merging mapping info in dotplot"

    set hiddenobjects {}
    for {set s 1} {$s<=$seq(n_seq)} {incr s} {
    	if { ! [cs_proj::mapinfo::exists $s]} {
    		p::debug "skipping seq $s, no info found"
    		continue
    	}
    	# FIXME: which is more work: find all dpbp tags or iterate over each info entry?
    	p::debug "merging mapping info for seq $s"
    	set degapindex 0

    	for {set n 1} {$n<=$seq(aln_len)} {incr n} {
    		if {! [IsGap [string index $seq(nt,$s) [expr {$n-1}]]]} {
    			incr degapindex
    		}

    		set notpaired [cs_proj::mapinfo::is_not_paired $s $degapindex]
    		if {$notpaired==-1} {
    			# don't know anything about this one
    			continue

    		} elseif {$notpaired==1} {
    			# should not pair: hide all possible bps
    			p::debug "hiding all dpbps for seq $s where nt $n participates"
    			HideThisNtDpBps $s $n
    			lappend hiddenobjects [join [HideThisNtDpBps $s $n]]

    		} else {
    			# should pair
    			set degappairedwith  [cs_proj::mapinfo::is_paired_with $s $degapindex]
    			if {$degappairedwith==0} {
    				continue ;# don't know partner
    			}
    			set pairedwith [MapNtIndex $seq(nt,$s) $degappairedwith DEGAPPED_TO_ALIGNED]
    			if {$pairedwith==-1} {
    				p::error "Upss...MapNtIndex failed for seq(nt,$s) ntidx $degappairedwith"
    			}

    			# hide all possible bps and reactivate the set one
    			#
    			lappend hiddenobjects [join [HideThisNtDpBps $s $n]]
    			# nti>ntj treat only this one
    			if {$n>$pairedwith} {
    				set goodtag DpBp_seq_${s}_nti_${n}_ntj_${pairedwith}
    			} else {
    				set goodtag DpBp_seq_${s}_nti_${pairedwith}_ntj_${n}
    			}
    			DpShowHide $goodtag 1
    			set replidx [lsearch -exact $hiddenobjects $goodtag]
    			if {$replidx==-1} {
    				p::verbose "$goodtag should pair but was not predicted"
    			} else {
    				set hiddenobjects [lreplace $hiddenobjects $replidx $replidx]
    			}
    		}
    	}
    }
    return $hiddenobjects
}
# MergeMappingInDotplot



###   RemoveMappingFromDotplot
#
# 
#
proc RemoveMappingFromDotplot {dptags} {
    ################

    p::verbose "Removing mapping info from dotplot"

    foreach id $dptags {
    	p::debug "showing $id"
    	DpShowHide $id 1
    }
}
# RemoveMappingFromDotplot
