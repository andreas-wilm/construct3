##############################################################################
#
# callbacks.tcl - procedures for handling bindings and actions invoked
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
#  CVS $Id: callbacks.tcl,v 1.57 2007-10-26 13:51:13 steger Exp $
#


###   BindAlignment   
#
#
proc BindAlignment {num_seq} {
    ########################
    global c         ;# r
	 global c3
    global lastfocus ;# rw
    

    ##########   Alignment
    #
    #

    # allow up and down
    #
    foreach ca [list $c(al) $c3(al) $c(al_sn) ] {
        bind $ca <Enter> {
            set lastfocus [focus]
            focus %W
        }
        bind $ca <KeyPress-Up> {
            if {$selection(seq_no)!=1} {
                #SeqSelectionChange [expr {$selection(seq_no)-1}]
                set newseq [expr {$selection(seq_no)-1}]
                event generate $c(al_sn).l_sn$newseq <1>
                MoveAlnView $newseq 1 ${seq(aln_len)}
            }
        }
        bind $ca <KeyPress-Down> {
            if {$selection(seq_no)<$seq(n_seq)} {
            #SeqSelectionChange [expr {$selection(seq_no)+1}]
                set newseq [expr {$selection(seq_no)+1}]
                event generate $c(al_sn).l_sn$newseq <1>
                MoveAlnView $newseq 1 ${seq(aln_len)}
            }
        }
    }

    
    ### Sequence labels-> Selection of  another sequence
    #
    #
    for {set s 1} {$s<=$num_seq} {incr s} {
        bind $c(al_sn).l_sn${s} <1> "SeqSelectionChange $s"
    }

    ### Nt selection
    #
    $c(al) bind AlnNt <1> {NtSelectionChange}

    $c(al) bind AlnNt <Control-1> {SuperSelectionChange}



    # mouse over aligned nucleotide: fill red
    #
    $c(al) bind AlnNt <Any-Enter> {
        $c(al) itemconfigure current -fill $COLOR(ActiveNt)

        # colorize corresponding dp bp
        #
        set seq_no -1; set nt_no -1;
        foreach tag [$c(al) gettags current] {
            if {[regexp {AlnNt_seq_([0-9]*)_nt_([0-9]*)} $tag \
                                           dummy seq_no nt_no]} {
                break
            }
        }
        Highlight_DpBp    $seq_no $nt_no on
        Highlight_DpSeqNt $nt_no  on
    }


    # mouse leaves nucleotide : change to old color
    #
    $c(al) bind AlnNt <Any-Leave> {
        $c(al) itemconfigure current -fill $COLOR(UnactiveNt)
        
        # decolorize corresponding dp bp
        #
        set seq_no -1; set nt_no -1;
        foreach tag [$c(al) gettags current] {
            if {[regexp {AlnNt_seq_([0-9]*)_nt_([0-9]*)} $tag \
                                           dummy seq_no nt_no]} {
                break
            }
        }
        Highlight_DpBp    $seq_no $nt_no off
        Highlight_DpSeqNt $nt_no  off

        # reset to selection color if this one was selected
        if {$selection(seq_no)==$seq_no} {
            $c(dp) itemconfigure DpBp_seq_${seq_no} -fill $COLOR(SelSeq)
        }
    }


    # Special Mouse binding
    # we cannot use the generic MouseScroll function
    # since we always need to move 2 canvas
    bind $c(al) <Button-5>         {
        $c(al)    yview scroll 5 units
        $c(al_sn) yview scroll 5 units
        $c3(al)   yview scroll 5 units
    }
    bind $c(al) <Button-4>         {
        $c(al)    yview scroll -5 units
        $c(al_sn) yview scroll -5 units
        $c3(al)   yview scroll -5 units
    }
    bind $c(al) <Shift-Button-5>   {
        $c(al)    xview scroll 5 units
        $c(al_nm) xview scroll 5 units
    }
    bind $c(al) <Shift-Button-4>   {
        $c(al)    xview scroll -5 units
        $c(al_nm) xview scroll -5 units
    }
    bind $c3(al) <Button-5>         {
        $c(al)    yview scroll 5 units
        $c(al_sn) yview scroll 5 units
        $c3(al)   yview scroll 5 units
    }
    bind $c3(al) <Button-4>         {
        $c(al)    yview scroll -5 units
        $c(al_sn) yview scroll -5 units
        $c3(al)   yview scroll -5 units
    }
    bind $c3(al) <Shift-Button-5>   {
        $c3(al)    xview scroll 5 units
        $c3(al_nm) xview scroll 5 units
    }
    bind $c3(al) <Shift-Button-4>   {
        $c3(al)    xview scroll -5 units
        $c3(al_nm) xview scroll -5 units
    }
    
    #puts "TMP_DEBUG(BindAlignment): mouse bind"
    #$c(al) bind AlnNt <1> {puts "c(al) mouse position %x %y"}
	 
	 
	##kl 
	 
   ### Nt selection
    #
    $c3(al) bind AlnNt <1> {NtSelectionChange}

    $c3(al) bind AlnNt <Control-1> {SuperSelectionChange}



    # mouse over aligned nucleotide: fill red
    #
	
    $c3(al) bind AlnNt <Any-Enter> {
        $c3(al) itemconfigure current -fill $COLOR(ActiveNt)
       
        # colorize corresponding dp bp
        #
        set seq_no -1; set nt_no -1;
        foreach tag [$c3(al) gettags current] {
            if {[regexp {AlnNt_seq_([0-9]*)_nt_([0-9]*)} $tag \
                                           dummy seq_no nt_no]} {
                break
            }
        }
        Highlight_DpBp    $seq_no $nt_no on
        Highlight_DpSeqNt $nt_no  on
    }


    # mouse leaves nucleotide : change to old color
    #
    $c3(al) bind AlnNt <Any-Leave> {
        $c3(al) itemconfigure current -fill $COLOR(UnactiveNt)

        # decolorize corresponding dp bp
        #
        set seq_no -1; set nt_no -1;
        foreach tag [$c3(al) gettags current] {
            if {[regexp {AlnNt_seq_([0-9]*)_nt_([0-9]*)} $tag \
                                           dummy seq_no nt_no]} {
                break
            }
        }
        Highlight_DpBp    $seq_no $nt_no off
        Highlight_DpSeqNt $nt_no  off

        # reset to selection color if this one was selected
        if {$selection(seq_no)==$seq_no} {
            $c(dp) itemconfigure DpBp_seq_${seq_no} -fill $COLOR(SelSeq)
        }
    }


#    # Special Mouse binding
 #   # we cannot use the generic MouseScroll function
  #  # since we always need to move 2 canvas
   # bind $c3(al) <Button-5>         {
    #    $c3(al)    yview scroll 5 units
     #   $c(al_sn) yview scroll 5 units
#    }
 #   bind $c3(al) <Button-4>         {
  #      $c3(al)    yview scroll -5 units
   #     $c(al_sn) yview scroll -5 units
    #}
#    bind $c3(al) <Shift-Button-5>   {
 #       $c3(al)    xview scroll 5 units
  #      $c(al_nm) xview scroll 5 units
   # }
    #bind $c3(al) <Shift-Button-4>   {
#        $c3(al)    xview scroll -5 units
 #       $c(al_nm) xview scroll -5 units
  #  }


    #puts "TMP_DEBUG(BindAlignment): mouse bind"
    #$c(al) bind AlnNt <1> {puts "c(al) mouse position %x %y"}
   
	
	##kl
	
	
}
# BindAlignment




###   BindDotplot   
#
#
# FIXME: replace tag regexp search with scan (faster!)
#
proc BindDotplot {} {
    ###############
    global c         ;# r
    global c3 
    global opts      ;# r

    # foreach ca [list $c(al) $c3(al) ] 
	 
    ### display bp position
    #
    bind $c(dp) <Motion>    "MousePos_to_BpIndex %x %y"; # show mouse position in nt coordinates
    bind $c(dp) <Any-Leave> "$c(dp) coords Showrec 0 0 -1 -1"



    ###########   DpBp
    #
    #
    #

    $c(dp) bind DpBp <Any-Enter> {

        # highlight corresponding alignment nt
        #
        set seq_no -1; set bp_idx -1;

        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpBp_seq_([0-9]*)_nti_([0-9]*)_ntj_([0-9]*)}  $tag \
                                                    dummy seq_no nt_i nt_j]} {
                break
            }
        }
        Highlight_AlnNt_From_DpBp $seq_no $nt_i $nt_j on

        catch {unset tag dummy seq_no nt_i nt_j}
    }


    $c(dp) bind DpBp <Any-Leave> {

        # de-highlight corresponding alignment nt
        #
        set seq_no -1; set bp_idx -1;
        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpBp_seq_([0-9]*)_nti_([0-9]*)_ntj_([0-9]*)}  $tag \
                                                    dummy seq_no nt_i nt_j]} {
                break
            }
        }
        Highlight_AlnNt_From_DpBp $seq_no $nt_i $nt_j off

        catch {unset tag dummy seq_no nt_i nt_j}
    }


    $c(dp) bind DpBp <3> {
        set seq_no -1; set nt_i -1; set nt_j -1;
        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpBp_seq_([0-9]*)_nti_([0-9]*)_ntj_([0-9]*)}  $tag \
                                                    dummy seq_no nt_i nt_j]} {
                break
				 }     
        }
        MoveAlnView $seq_no $nt_j $nt_i
        catch {unset tag dummy seq_no nt_i nt_j}
    }
	 
    $c(dp) bind DpBp <1> {
        set seq_no -1; set nt_i -1; set nt_j -1;
        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpBp_seq_([0-9]*)_nti_([0-9]*)_ntj_([0-9]*)}  $tag \
                                                    dummy seq_no nt_i nt_j]} {
                break
				 }     
        }
        MoveAlnView $seq_no $nt_j $nt_i
        catch {unset tag dummy seq_no nt_i nt_j}
    }
	 
	 

    ###########   DpConsBp
    #
    #
    #

    ###   Highlight

    $c(dp) bind DpConsBp <Any-Enter> {
        if {$opts(displ,highlight_nt_from_csbp)} {
            foreach tag [$c(dp) gettags current] {
                if {[regexp {DpConsBp_nti_([0-9]*)_ntj_([0-9]*)} $tag \
                                                     dummy nt_i nt_j]} {
                    break
                }
            }
            Highlight_AlnNt_From_DpConsBp $nt_i $nt_j on
            catch {unset tag dummy nt_i nt_j}
        }
    }

    $c(dp) bind DpConsBp <Any-Leave> {
        if {$opts(displ,highlight_nt_from_csbp)} {
            foreach tag [$c(dp) gettags current] {

                if {[regexp {DpConsBp_nti_([0-9]*)_ntj_([0-9]*)} $tag \
                                                     dummy nt_i nt_j]} {
                    break
                }
            }
            Highlight_AlnNt_From_DpConsBp $nt_i $nt_j off
            catch {unset tag dummy nt_i nt_j}
        }
    }


    ###   Move Alignment window

    $c(dp) bind DpConsBp <1> {
        set nt_i -1; set nt_j -1;
        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpConsBp_nti_([0-9]*)_ntj_([0-9]*)} $tag \
                                                      dummy nt_i nt_j]} {
                break
            }
        }
        MoveAlnView "consensus" $nt_j $nt_i
        catch {unset tag dummy nt_i nt_j}
    }

    ###   Info

    $c(dp) bind DpConsBp <2> {
        set nt_i -1; set nt_j -1;
        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpConsBp_nti_([0-9]*)_ntj_([0-9]*)} $tag \
                                                      dummy nt_i nt_j]} {
                break
            }
        }
        puts [ConsBpInfo $nt_i $nt_j]
        catch {unset tag dummy nt_i nt_j}
    }


    MouseScroll $c(dp)

   ##kl
    ####### Get_MiniAlNuc_from_Dp
    #
    # gets ConsensusBasePairs from Dp (with right mouse-Click)
    # as start and end points
    # for extracting a Minialignment
    #
	global minial

	set minial(nuc1) "1st BP 5'"
	set minial(nuc2) "2nd BP 5'"
    set minial(nuc3) "1st BP 3'"
	set minial(nuc4) "2nd BP 3'"
	
	
	
    $c(dp) bind DpConsBp <3> {
        set nt_i -1; set nt_j -1;
        foreach tag [$c(dp) gettags current] {
            if {[regexp {DpConsBp_nti_([0-9]*)_ntj_([0-9]*)} $tag dummy nt_i nt_j]} {
                break
            }
        }

        # second nt selected
        if {$minial(nuc4) == $minial(nuc3) &&  \
            $minial(nuc3) != "1st BP 3'"        } {
            if {$minial(nuc4) < $nt_i } {
                set dummy        $minial(nuc4)
			    set minial(nuc4) $nt_i 
                set minial(nuc3) $dummy 
				
	   	        set dummy        $minial(nuc1)
			    set minial(nuc1) $nt_j
			    set minial(nuc2) $dummy 
			    MiniAl
            } else {
                set minial(nuc3) $nt_i 
			    set minial(nuc2) $nt_j
			    MiniAl
            }
        # first bp selected
        } else {
            set minial(nuc4) $nt_i
            set minial(nuc1) $nt_j  
            set minial(nuc3) $nt_i
            set minial(nuc2) $nt_j
            MiniAl
        }
 	    MoveAlnView "consensus" $nt_j $nt_i
        catch {unset tag dummy nt_i nt_j}
	}
  ##kl

}
# BindDotplot

###   Highlight_DpSeqNt   
#
#
# Highlights the Nts beside Dotplot
# state: on|off
#
proc Highlight_DpSeqNt {nt_no state} {
######################
    global COLOR  ;# r
    global c      ;# r

    if {$state=="off"} {
        set color $COLOR(UnactiveNt)
    } else {
        set color $COLOR(ActiveNt)
    }
    foreach canv [list $c(seq_top) $c(seq_below) $c(seq_left) $c(seq_right)] {
        $canv itemconfigure DpSeqLabel_ntidx_${nt_no} -fill $color
    }
}
# Highlight_DpSeqNt




###    SeqSelectionChange   
#
# new_sel_seq = index of selected sequence
# if DESELECT then deselect
#
proc SeqSelectionChange {new_sel_seq} {
#######################
    global c            ;# r
    global c3           ;# r
    global COLOR        ;# r
    global selection    ;# rw
    global seq          ;# r


    # sequence already selected ?
    if {$selection(seq_no)==$new_sel_seq} {
        return
    }

    # deselect ?
    if {$new_sel_seq=="DESELECT"} {
        # deselect nt in old seq
        $c(al) delete AlnSelBox

        # change sequence label color
        #
        $c(al_sn).l_sn${selection(seq_no)} config -fg black -bg $COLOR(UnselSeq)

        # change dotplot basepair color
        #
        $c(dp) itemconfigure DpBp_seq_${selection(seq_no)} -fill $COLOR(UnselSeq)

        # reset nt selection
        #
        set selection(seq_no)      -1
        set selection(nt_low)      -1
        set selection(nt_high)     -1
        set selection(super_sel) 0
        return
    }



    if {$selection(seq_no)==-1} {
        set  selection(seq_no) $new_sel_seq
    } else {
        # deselect nt in old seq
        $c(al) delete AlnSelBox
    }

    # change sequence label color
    #
    $c(al_sn).l_sn${selection(seq_no)} config -fg black -bg $COLOR(UnselSeq)
    $c(al_sn).l_sn${new_sel_seq} config -fg white -bg $COLOR(SelSeq)


    # change dotplot basepair color
    #
    $c(dp) itemconfigure DpBp_seq_${selection(seq_no)} -fill $COLOR(UnselSeq)
    $c(dp) itemconfigure DpBp_seq_${new_sel_seq} -fill $COLOR(SelSeq)




    set selection(seq_no) $new_sel_seq

    ### raise selected seq dotplot bps
    #
    $c(dp) raise DpBp_seq_${new_sel_seq} DpBp

    # reset nt selection
    #
    set selection(nt_low)      -1
    set selection(nt_high)     -1


    UpdateDpSeqLabels $seq(nt,$new_sel_seq)

    set selection(super_sel) 0
}
# SeqSelectionChange





###   MoveAlnView   
#
# Move Alignment Window to Dotplot-clicked Nucleotide
#
#
proc MoveAlnView {seq_idx nt_idx nt3_idx} {
################
    global c         ;# r
    global c3 
    global seq       ;# r
    global selection ;# r

    set local_debug 0

    if {$local_debug} {
        puts "DEBUG(MoveAlnView): seq_idx=$seq_idx nt_idx=$nt_idx"
    }

    if {$seq_idx=="consensus"} {
        set seq_idx $selection(seq_no)
        if {$seq_idx==-1} {set seq_idx 1}
    }


    # use Coord instead of tag first, since we don't know if max elems are gaps

    set tmp_coords  [GetNtAliCoords 1 1]
    set nt_coords(min,x) [lindex $tmp_coords 0]
    set nt_coords(min,y) [lindex $tmp_coords 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): nt_coords(min,x)=$nt_coords(min,x) nt_coords(min,y)=$nt_coords(min,y)"
    }
    #
    set tmp_coords  [GetNtAliCoords ${seq(n_seq)} ${seq(aln_len)}]
    set nt_coords(max,x) [lindex $tmp_coords 0]
    set nt_coords(max,y) [lindex $tmp_coords 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): nt_coords(max,x)=$nt_coords(max,x) nt_coords(max,y)=$nt_coords(max,y)"
    }
    #
    set tmp_coords  [GetNtAliCoords ${seq_idx} ${nt_idx}]
    set nt_coords(cur,x) [lindex $tmp_coords 0]
    set nt_coords(cur,y) [lindex $tmp_coords 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): nt_coords(cur,x)=$nt_coords(cur,x) nt_coords(cur,y)=$nt_coords(cur,y)"
    }
	 #
	 set nt_coords(rel,x)  [expr {$nt_coords(cur,x)/$nt_coords(max,x)}]
    set nt_coords(rel,y)  [expr {$nt_coords(cur,y)/$nt_coords(max,y)}]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): nt_coords(rel,x)=$nt_coords(rel,x) nt_coords(rel,y)=$nt_coords(rel,y)"
    }

	 
	 
	 
	 #kl
	 
    set tmp3_coords  [GetNtAliCoords ${seq_idx} ${nt3_idx}]
    set nt_coords(cur3,x) [lindex $tmp3_coords 0]
    set nt_coords(cur3,y) [lindex $tmp3_coords 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): nt_coords(cur3,x)=$nt_coords(cur3,x) nt_coords(cur3,y)=$nt_coords(cur3,y)"
    }
    #
    set nt_coords(rel3,x)  [expr {$nt_coords(cur3,x)/$nt_coords(max,x)}]
    set nt_coords(rel3,y)  [expr {$nt_coords(cur3,y)/$nt_coords(max,y)}]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): nt_coords(rel,x)=$nt_coords(rel,x) nt_coords(rel,y)=$nt_coords(rel,y)"
    }

    #kl





    # rel. visibility coords of alignment canvas
    #
    set al_vis_start(x) [lindex [$c(al) xview] 0]
    set al_vis_end(x)   [lindex [$c(al) xview] 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): al_vis_start(x)=$al_vis_start(x) al_vis_end(x)=$al_vis_end(x)"
    }
    #
    set al_vis_start(y) [lindex [$c(al) yview] 0]
    set al_vis_end(y)   [lindex [$c(al) yview] 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): al_vis_start(y)=$al_vis_start(y) al_vis_end(y)=$al_vis_end(y)"
    }


    # move x
    #
    if {$nt_coords(rel,x)<$al_vis_start(x) || $nt_coords(rel,x)>$al_vis_end(x)} {
        if {$local_debug} {puts "nt NOT VISIBLE in x view"}
        set al_width        [expr {$al_vis_end(x)-$al_vis_start(x)}]
        set fract_off_left  [expr {$nt_coords(rel,x)-$al_width/2}]
        $c(al)    xview moveto $fract_off_left
        $c(al_nm) xview moveto $fract_off_left
    } else {
        if {$local_debug} {puts "nt VISIBLE in x view"}
    }

    # move y
    #
    if {$nt_coords(rel,y)<$al_vis_start(y) || $nt_coords(rel,y)>$al_vis_end(y)} {
        if {$local_debug} {puts "nt NOT VISIBLE in y view"}
        set al_height      [expr {$al_vis_end(y)-$al_vis_start(y)}]
        set fract_off_top  [expr {$nt_coords(rel,y)-$al_height/2}]
        $c(al)    yview moveto $fract_off_top
        $c(al_sn) yview moveto $fract_off_top
    } else {
        if {$local_debug} {puts "nt VISIBLE in y view"}
    }



  #kl
  # rel. visibility coords of alignment canvas
    #
    set al_vis_start(x) [lindex [$c3(al) xview] 0]
    set al_vis_end(x)   [lindex [$c3(al) xview] 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): al_vis_start(x)=$al_vis_start(x) al_vis_end(x)=$al_vis_end(x)"
    }
    #
    set al_vis_start(y) [lindex [$c3(al) yview] 0]
    set al_vis_end(y)   [lindex [$c3(al) yview] 1]
    if {$local_debug} {
        puts "DEBUG(MoveAlnView): al_vis_start(y)=$al_vis_start(y) al_vis_end(y)=$al_vis_end(y)"
    }


    # move x
    #
    if {$nt_coords(rel3,x)<$al_vis_start(x) || $nt_coords(rel3,x)>$al_vis_end(x)} {
        if {$local_debug} {puts "nt NOT VISIBLE in x view"}
        set al_width        [expr {$al_vis_end(x)-$al_vis_start(x)}]
        set fract_off_left3  [expr {$nt_coords(rel3,x)-$al_width/2}]
        $c3(al)    xview moveto $fract_off_left3
        $c3(al_nm) xview moveto $fract_off_left3
    } else {
        if {$local_debug} {puts "nt VISIBLE in x view"}
    }

    # move y
    #
    if {$nt_coords(rel3,y)<$al_vis_start(y) || $nt_coords(rel3,y)>$al_vis_end(y)} {
        if {$local_debug} {puts "nt NOT VISIBLE in y view"}
        set al_height      [expr {$al_vis_end(y)-$al_vis_start(y)}]
        set fract_off_top3  [expr {$nt_coords(rel3,y)-$al_height/2}]
        $c3(al)    yview moveto $fract_off_top3
        $c(al_sn) yview moveto $fract_off_top3
    } else {
        if {$local_debug} {puts "nt VISIBLE in y view"}
    }
    #kl

}
# MoveAlnView





###    NtSelectionChange   
#
# New nucleotides have been selected in alignment
#
proc NtSelectionChange {} {
######################
    global selection;# rw
    global c        ;# r
    global c3

    # which nt and sequence?
    #
    foreach ca [list $c(al) $c3(al)] {
        foreach tag [$ca gettags current] {
            if {[regexp {.*_seq_([0-9]+)_nt_([0-9]+)} $tag \
                           dummy new_seq_idx new_nt_idx]} {

                # new sequence selected ?
                #
                if {$new_seq_idx != $selection(seq_no)} {
                    SeqSelectionChange $new_seq_idx
                }
                break
            }
        }
    }

    # second nt selected
    #
    if { ($selection(nt_low)==$selection(nt_high)) && ($selection(nt_low)!=-1)} {

        if {$selection(nt_low) < $new_nt_idx} {
            set selection(nt_high)  $new_nt_idx
        } else {
            set selection(nt_low)   $new_nt_idx
        }
        BoxNtSelection $selection(seq_no) $selection(seq_no) \
                       $selection(nt_low) $selection(nt_high) "add"

    # first nt selected
    #
    } else {

        # decolorize old selection
        BoxNtSelection $selection(seq_no) $selection(seq_no) \
                       $selection(nt_low) $selection(nt_high) "remove"

        set selection(nt_low) $new_nt_idx
        set selection(nt_high) $new_nt_idx
        BoxNtSelection $selection(seq_no) $selection(seq_no) \
                       $selection(nt_low) $selection(nt_high) "add"

    }
}
# NtSelectionChange



###   EmulateGapInsertion
#
# move consecutive succession of nucleotides to next gap
#
proc EmulateGapInsertion {direction} {
########################
    global selection
    global seq
    
    if {$selection(nt_low)==-1} {
        puts "No nucleotide selected!"
        return
    }

    if {$selection(super_sel)!=0} {
        puts "This works for single sequence selection only"
        return
    }


    if {$direction==-1} {
        set fixed_idx $selection(nt_high)
        set step_idx  $selection(nt_low)
    } elseif {$direction==1} {
        set fixed_idx $selection(nt_low)
        set step_idx  $selection(nt_high)
        
    } else {
        p:error "wrong direction \"$direction\""
    }

    set this_seq $seq(nt,$selection(seq_no))

    # grow selection to next gap
    while { $step_idx!=0 && $step_idx!=$seq(aln_len) } {
        set this_nt [string index $this_seq [expr {$step_idx-1+$direction}]]
        
        if {[IsGap $this_nt]} {
            
            if {$step_idx>$fixed_idx} {
                set selection(nt_low) $fixed_idx
                set selection(nt_high) $step_idx
            } else {
                set selection(nt_low) $step_idx
                set selection(nt_high) $fixed_idx
            }
            
            BoxNtSelection $selection(seq_no) $selection(seq_no) \
                $selection(nt_low) $selection(nt_high) "add"

            RequestToMove $direction

            break
        }
        incr step_idx $direction
    }
}
# EmulateGapInsertion



###   SuperSelectionChange   
#
#
proc SuperSelectionChange {} {
#########################
    global selection;# rw
    global c        ;# r
    global c3


    # paranoia: no sequence preselected ? return
    if {$selection(seq_no)==-1} {
        p::debug "new sequence preselected ! falling back"
        return
    }

    # which nt and sequence?
    #
    foreach tag [$c(al) gettags current] {
        if {[regexp {.*_seq_([0-9]+)_nt_([0-9]+)} $tag \
                           dummy new_seq_idx new_nt_idx]} {

            # return if same sequence has been selected
            if {$new_seq_idx == $selection(seq_no)} {
                p::debug "same sequence! returning"
                return
            }
            break
        }
    }
	 
	 # which nt and sequence? c3(al)
    #
    foreach tag [$c3(al) gettags current] {
        if {[regexp {.*_seq_([0-9]+)_nt_([0-9]+)} $tag \
                           dummy new_seq_idx new_nt_idx]} {

            # return if same sequence has been selected
            if {$new_seq_idx == $selection(seq_no)} {
                p::debug "same sequence! returning"
                return
            }
            break
        }
    }
    #puts "DEBUG(SuperSelectionChange): new_seq_idx = $new_seq_idx"
    #puts "DEBUG(SuperSelectionChange): new_nt_idx  = $new_nt_idx"
    #puts "DEBUG(SuperSelectionChange): new_seq_idx = $new_seq_idx"
    #puts "DEBUG(SuperSelectionChange): selection(seq_no)  = $selection(seq_no)"
    #puts "DEBUG(SuperSelectionChange): old selection(nt_low)  = $selection(nt_low)"
    #puts "DEBUG(SuperSelectionChange): old selection(nt_high) = $selection(nt_high)"



    if {$selection(nt_low) < $new_nt_idx} {
        set selection(nt_high)  $new_nt_idx
    } else {
        set selection(nt_low)   $new_nt_idx
    }
    BoxNtSelection $new_seq_idx $selection(seq_no) $selection(nt_low) $selection(nt_high) "add"

    #puts "DEBUG(SuperSelectionChange): new selection(nt_low)  = $selection(nt_low)"
    #puts "DEBUG(SuperSelectionChange): new selection(nt_high) = $selection(nt_high)"


    # update the list of selected sequences
    #
    set selection(super_sel) 1
    if {$new_seq_idx > $selection(seq_no)} {
        set start $selection(seq_no)
        set end   $new_seq_idx
    } else {
        set start $new_seq_idx
        set end   $selection(seq_no)
    }
    set selection(super_sel,seq_indices) {}
    for {set i $start} {$i<=$end} {incr i} {
        lappend selection(super_sel,seq_indices) $i
    }
}
# SuperSelectionChange




###   BoxNtSelection   
#
# Draws a box around selected nucleotides in alignment window
#
#
proc BoxNtSelection {seq_idx1 seq_idx2 nt_idx1 nt_idx2 action} {
###################
    global c        ;# r
    global c3       ;# r	 
    global COLOR    ;# r
    global seq;     ;# r

    if {$action=="remove"} {
        $c(al)  delete AlnSelBox
        $c3(al) delete AlnSelBox
        return
    }

    if {$seq_idx1>$seq_idx2} {
        set dummy $seq_idx1
        set seq_idx1 $seq_idx2
        set seq_idx2 $dummy
    }

    # if multiple sequences are selected,
    # none of them may have a gap at either end
    # (this is implicitly the same case for single selected sequences)
    for {set s $seq_idx1} {$s<=$seq_idx2} {incr s} {
        set seqslice [string range $seq(nt,$s) \
                     [expr {$nt_idx1-1}] [expr {$nt_idx2-1}]]
        if {[IsGap [string index $seqslice 0]]
                  ||
            [IsGap [string index $seqslice end]]} {
                puts "Selection for sequence \"$seq(id,$s)\" has a gap at one of it's ends."
                puts "Please exclude it for this multiple-selection"
                return
        }
    }


    set coords_1  [$c(al) bbox AlnNt_seq_${seq_idx1}_nt_${nt_idx1}]
    set coords_2  [$c(al) bbox AlnNt_seq_${seq_idx2}_nt_${nt_idx2}]

    set x1_coord_list [list [lindex $coords_1 0] [lindex $coords_2 0]]
    set x2_coord_list [list [lindex $coords_1 2] [lindex $coords_2 2]]
    set y1_coord_list [list [lindex $coords_1 1] [lindex $coords_2 1]]
    set y2_coord_list [list [lindex $coords_1 3] [lindex $coords_2 3]]

    set x1 [lindex [lsort -real $x1_coord_list] 0]
    set x2 [lindex [lsort -real $x2_coord_list] 1]
    set y1 [lindex [lsort -real $y1_coord_list] 0]
    set y2 [lindex [lsort -real $y2_coord_list] 1]


    foreach ca [list $c(al) $c3(al)] {
        set box [$ca create rectangle $x1 $y1 $x2 $y2  -tags AlnSelBox \
                                  -fill $COLOR(SelNtBox) -outline "" ]
        $ca lower $box
        $ca lower REMatchBox
    }
}
# BoxNtSelection





###   AdjustAlnCanvas   
#
# Must be called b4 Sequences are changed !
#
# returns the x-coords difference (for AlnSelBox movement)
#
proc AdjustAlnCanvas {seq_no nt_low nt_high offset} {
####################
    global c      ;# r
    global c3      ;# r
    global seq    ;# r

    ###   calc coord difference
    #     by getting the coord from nt and the gap beside
    #
    #
    if {$offset==-1} {
        set old_gap_tag "Gap_seq_${seq_no}_pos_[expr {$nt_low-1}]"
        set new_gap_tag "Gap_seq_${seq_no}_pos_${nt_high}"

        set leader $nt_low
        set trail  $nt_high
    } elseif {$offset==1} {
        set old_gap_tag "Gap_seq_${seq_no}_pos_[expr {$nt_high+1}]"
        set new_gap_tag "Gap_seq_${seq_no}_pos_${nt_low}"

        set leader $nt_high
        set trail $nt_low
    } else {
        p::error "Invalid offset $offset"
        return
    }



    set x_diff [GetNtAliCellDim width]
    if {$offset==-1} {
        set x_diff [expr {$x_diff*-1.0}]
    }


    ###   move and re-tag the block nts/gaps
    #
    #     order differs dependent on offset (otherwise tags will be duplicated)
    #
    if {$offset==-1} {
        for {set i $nt_low} {$i<=$nt_high} {incr i +1} {
            # tags must be set according to current residue inside moved block
            set residue [string index $seq(nt,$seq_no) [expr {$i-1}]]
            if {[IsGap $residue]} {
                set old_tag "Gap_seq_${seq_no}_pos_$i"
                set new_tag "Gap_seq_${seq_no}_pos_[expr {$i+$offset}]"
            } else {
                set old_tag "AlnNt_seq_${seq_no}_nt_$i"
                set new_tag "AlnNt_seq_${seq_no}_nt_[expr {$i+$offset}]"
            }
            $c(al)  move $old_tag "${x_diff}c" 0
            $c3(al) move $old_tag "${x_diff}c" 0
            ReplaceTag $c(al)  $old_tag $new_tag
            ReplaceTag $c3(al) $old_tag $new_tag
        }
    } else {
        for {set i $nt_high} {$i>=$nt_low} {incr i -1} {
            # tags must be set according to current residue inside moved block
            set residue [string index $seq(nt,$seq_no) [expr {$i-1}]]
            if {[IsGap $residue]} {
                set old_tag "Gap_seq_${seq_no}_pos_$i"
                set new_tag "Gap_seq_${seq_no}_pos_[expr {$i+$offset}]"
            } else {
                set old_tag "AlnNt_seq_${seq_no}_nt_$i"
                set new_tag "AlnNt_seq_${seq_no}_nt_[expr {$i+$offset}]"
            }
            $c(al)  move $old_tag "${x_diff}c" 0
            $c3(al) move $old_tag "${x_diff}c" 0
            ReplaceTag $c(al)  $old_tag $new_tag
            ReplaceTag $c3(al) $old_tag $new_tag
        }
    }


    ###   flip and re-tag overwritten gap
    #
    set gap_offset  [expr {-$x_diff*($nt_high-$nt_low+1)}]
    $c(al)  move       $old_gap_tag "${gap_offset}c" 0
    $c3(al) move       $old_gap_tag "${gap_offset}c" 0
    ReplaceTag $c(al)  $old_gap_tag $new_gap_tag
    ReplaceTag $c3(al) $old_gap_tag $new_gap_tag
    return $x_diff
}
# AdjustAlnCanvas



###   RequestToMove   
#
# Request for alignment move/modification
#
proc RequestToMove {offset} {
##################
    global selection                 ;# rw
    global seq                       ;# rw
    global alignment_is_modified     ;# w
    global mic_is_computed_and_valid ;# r
    global c                         ;# r
    global c3                        ;# r
    global opts                      ;# r
    global seqsearch_re              ;# r


    if {$offset==-1} {
        set ori "left"
    } elseif {$offset==1} {
        set ori "right"
    } else {
        p::error "Invalid offset $offset"
        return
    }


    # no selection at all ?
    #
    if {$selection(seq_no)==-1 || $selection(nt_low)==-1} {
        puts "No nucleotide range selected!"
        return
    }


    if {$selection(super_sel)==1} {
        set sel_seqs $selection(super_sel,seq_indices)
    } else {
        set sel_seqs $selection(seq_no)
    }
    p::debug "sel_seqs = $sel_seqs"

    # check if movement is valid for each selected sequence
    # before doing anything
    #
    foreach move_seq $sel_seqs {
        if { ! [MoveRequestIsOk $move_seq \
                    $selection(nt_low) $selection(nt_high) $ori]} {
            return
        }
    }


    foreach move_seq $sel_seqs {
        set xdiff [MoveAsRequested $move_seq \
                       $selection(nt_low) $selection(nt_high) $offset]
    }

    # now it's save to update sequences in tcl
    #
    Seq_Exchange seq C2TCL


    ###   move the selection box
    #
    $c(al)  move AlnSelBox "${xdiff}c" 0
    $c3(al) move AlnSelBox "${xdiff}c" 0

    # delete re-match from this seq and redraw
    if {[PatternSearchActive]} {
        foreach move_seq $sel_seqs {
            PatternMatchDelete $move_seq

            for {set i 0} {$i<[llength $seqsearch_re(re_list)]} {incr i} {
                set re [lindex $seqsearch_re(re_list)      $i]
                set col [lindex $seqsearch_re(color_list)   $i]
                set seq_list [lindex $seqsearch_re(seqidx_lists) $i]

                # moved seq inside this applied re seq list
                if {[lsearch -exact $seq_list $move_seq]!=-1} {
                    SeqSearchAndMark   $move_seq $re $col
                }
            }
        }
    }

    # update dotplot sequence labels
    #
    # FIXME: select one and reset fully
    UpdateDpSeqLabels $seq(nt,$selection(seq_no))

    # only hide when requested
    # don't explicitly show when already drawn
    if {! $opts(displ,cons_bp)} {
        DpShowHide DpConsBp 0
    }
    if {! $opts(displ,bps)} {
        DpShowHide DpBp 0
    }
    if {! $opts(displ,gaps)} {
        DpShowHide Gap 0
    }




    # FIXME: needed ?
    update idletasks


    # update selection
    #
    incr selection(nt_low)  $offset
    incr selection(nt_high) $offset


    set  mic_is_computed_and_valid  0
    set  alignment_is_modified      1
}
# RequestToMove




###   MoveRequestIsOk   
#
#
# ori=left|right
#
proc MoveRequestIsOk {seq_no nt_low nt_high ori} {
####################
    global seq

    ###   If neighbour to requested offset is not a gap
    #     or if selection is at either end then return 0
    #
    if {$ori=="right"} {
        if {$nt_high==$seq(aln_len)} {
            p::verbose "You cannot move this block $ori, since it is already at end!"
            return 0
        }
        set neighbour_idx [expr {$nt_high-0}]

    # left
    } else {
        if {$nt_low==1} {
            p::verbose "You cannot move this block $ori, since it is already at beginning!"
            return 0
        }
        set neighbour_idx [expr {$nt_low-2}]
    }

    set neighbour_nt  [string index \
                       $seq(nt,$seq_no) $neighbour_idx]
    if { ! [IsGap $neighbour_nt]} {
        p::verbose "You cannot move this block $ori, since there is no gap!"
        return 0
    }

    return 1
}
# MoveRequestIsOk




###   MoveAsRequested   
#
# Returns the x coordinate difference (for AlnSelBox movement)
#
# offset = -1|1
#
proc MoveAsRequested {seq_no nt_low nt_high offset} {
####################
    global seq

    set debug_seq_slice [string range $seq(nt,$seq_no) \
                                    [expr {$nt_low-1}] [expr {$nt_high-1}]]
    p::debug   "moving $seq(id,$seq_no)"
    p::debug   "seq_no=$seq_no offset=$offset"
    p::debug   "nt_low=$nt_low   nt_high=$nt_high"
    p::debug   "that is \"$debug_seq_slice\""
    p::verbose "Moving $debug_seq_slice of $seq(id,$seq_no)"

    p::debug "MoveNtSelection  $seq_no $nt_low $nt_high $offset"

    set modseqslice [MoveNtSelection  $seq_no \
                                      $nt_low $nt_high $offset]

    p::debug "Calling: AdjustAlnCanvas  $seq_no $nt_low $nt_high $offset"

    set xdiff [AdjustAlnCanvas  $seq_no $nt_low $nt_high $offset]

    p::debug "xdiff=$xdiff"

    Debugger "callbacks"

    return $xdiff

}
# MoveAsRequested




###   UpdateDpSeqLabels   
#
# seq_nts maybe a slide
#
proc UpdateDpSeqLabels {seq_nts} {
######################
    global c     ;# r

    for {set n 1} {$n<=[string length $seq_nts]} {incr n} {

        set nt [string index $seq_nts [expr {$n-1}]]

        $c(seq_top)   itemconfigure DpSeqLabel_ntidx_$n -text "$nt"
        $c(seq_below) itemconfigure DpSeqLabel_ntidx_$n -text "$nt"
        $c(seq_left)  itemconfigure DpSeqLabel_ntidx_$n -text "$nt"
        $c(seq_right) itemconfigure DpSeqLabel_ntidx_$n -text "$nt"
    }
}
# UpdateDpSeqLabels




###   ZoomDotplot   
#
#
proc ZoomDotplot {scale_fac} {
################
    global c           ;# w
    global scale       ;# w
    global SEQHEIGHT   ;# r
    global FONT        ;# w
    global seq         ;# r

    p::debug "old scale(dp)=$scale(dp)"
    p::debug "old c(virtual_size)=$c(virtual_size)"


    SetCursor busy

    set c(virtual_size)  [expr {$c(virtual_size)*$scale_fac}]
    set scale(dp)        [expr {$scale(dp)*$scale_fac}]

    # remember the visible scrollregions
    #
    set dp_xview [$c(dp) xview]
    set vis_scr_reg(x0) [lindex $dp_xview 0]
    set vis_scr_reg(x1) [lindex $dp_xview 1]
    set dp_yview [$c(dp) yview]
    set vis_scr_reg(y0) [lindex $dp_yview 0]
    set vis_scr_reg(y1) [lindex $dp_yview 1]



    foreach id [$c(dp) find all] {
        $c(dp) scale ${id} 0 0 $scale_fac $scale_fac
    }
    $c(dp) config -scrollregion [list 0c 0c ${c(virtual_size)}c ${c(virtual_size)}c]


    # restore the visible scrollregion
    #
    set move_to_x [expr {($vis_scr_reg(x0)+$vis_scr_reg(x1))/2.0}]
    set move_to_x [expr {$move_to_x - ($vis_scr_reg(x1)-$vis_scr_reg(x0)) /2.0/$scale_fac}]
    set move_to_y [expr {($vis_scr_reg(y0)+$vis_scr_reg(y1))/2.0}]
    set move_to_y [expr {$move_to_y - ($vis_scr_reg(y1)-$vis_scr_reg(y0)) /2.0/$scale_fac}]
    # puts "TMP_DEBUG(ZoomDotplot): move_to_x=$move_to_x"
    # puts "TMP_DEBUG(ZoomDotplot): move_to_y=$move_to_y"

    $c(dp) xview moveto $move_to_x	;# left border of visible part
    $c(dp) yview moveto $move_to_y	;# top  border of visible part

    set pt_sh [expr {int($SEQHEIGHT * 100.0/2.54*10)}]
    set pt_sc [expr {int($scale(dp) * 100.0/2.54*10)}]
    if {$pt_sc < $pt_sh} {
        set FONT(dp_seq,pt) $pt_sc
    } else {
        set FONT(dp_seq,pt) $pt_sh
    }
    set FONT(dp_seq) [format $FONT(dp_seq,format) $FONT(dp_seq,pt)]


    set scr_reg [list 0c 0c ${c(virtual_size)}c ${c(virtual_size)}c]
    # puts "TMP_DEBUG(ZoomDotplot): scr_reg=$scr_reg"

    foreach id [$c(dp_num_above) find withtag DpNum] {
        $c(dp_num_above) scale $id 0 0 $scale_fac 1
        $c(dp_num_above) itemconfig $id -font $FONT(dp_nm)
    }
    $c(dp_num_above) config -scrollregion $scr_reg

    foreach id [$c(dp_num_left) find withtag DpNum] {
        $c(dp_num_left) scale $id 0 0 1 $scale_fac
        $c(dp_num_left) itemconfig $id -font $FONT(dp_nm)
    }
    $c(dp_num_left) config -scrollregion $scr_reg


    foreach canvas "$c(seq_top) $c(seq_below)" {
        foreach id [$canvas find withtag DpSeqLabel] {
            $canvas scale $id 0 0 $scale_fac 1
            $canvas itemconfig $id -font $FONT(dp_seq)
        }
        $canvas config -scrollregion $scr_reg
    }


    foreach canvas [list $c(seq_left) $c(seq_right)] {
        foreach id [$canvas find withtag DpSeqLabel] {
            $canvas scale $id 0 0 1 $scale_fac
            $canvas itemconfig $id -font $FONT(dp_seq)
        }
        $canvas config -scrollregion $scr_reg
    }

    p::debug "new scale(dp)=$scale(dp)"
    p::debug "new c(virtual_size)=$c(virtual_size)"

    SetCursor normal
}
# ZoomDotplot




###   MousePos_to_BpIndex   
#
# Converts mouse position in dotplot to a basepair index
# and displays it
# if the mirrored rectangle option is 1, display this also
#
proc MousePos_to_BpIndex {x y} {
########################
    global scale     ;# r
    global c         ;# r
    global lb_dp_pos ;# r
    global opts      ;# r


    set scale_pix [winfo fpixels $c(dp) ${scale(dp)}c]

    # calculate nt indices
    set i [expr {int(ceil([$c(dp) canvasx $x]/$scale_pix))}]
    set j [expr {int(ceil([$c(dp) canvasy $y]/$scale_pix))}]

    $lb_dp_pos config -text "BP $i:$j"


    # show the rectangle
    #
    if {$opts(displ,show_rectangle)} {
        $c(dp) coords MirroredPosRectangle [expr {($j-1)*$scale(dp)}]c \
                                           [expr {($i-1)*$scale(dp)}]c \
                                           [expr {$j*$scale(dp)}]c \
                                           [expr {$i*$scale(dp)}]c
        $c(dp) raise MirroredPosRectangle
    } else {
        $c(dp) coords MirroredPosRectangle 0c 0c -1c -1c
    }
}
# MousePos_to_BpIndex




###   DpScroll    
#
# routine used to scroll dotplotseq including the surrounding sequence and
# numbering windows
#
#
proc DpScroll {xy command pos {pages ""}} {
#############
   global c

   if {$xy == "x"} {
      set l_canvas [list $c(dp) $c(dp_num_above) $c(seq_top) $c(seq_below)]
   } else {
      set l_canvas [list $c(dp) $c(dp_num_left) $c(seq_left) $c(seq_right)]
   }
   foreach canvas $l_canvas {
        eval $canvas ${xy}view $command $pos $pages
        # puts "DEBUG(DpScroll): $canvas ${xy}view $command $pos $pages"
   }
}
# DpScroll




###   DpDrag   
#
# routine used to scroll dot plot including the surrounding sequence and
# numbering windows
#
proc DpDrag {sb_name xy lo hi} {
###########
   global c

   if {$xy == "x"} {
      $sb_name set $lo $hi
      $c(dp_num_above)  xview moveto $lo
      $c(seq_top)       xview moveto $lo
      $c(seq_below)     xview moveto $lo
   } else {
      $sb_name set $lo $hi
      $c(dp_num_left)   yview moveto $lo
      $c(seq_left)      yview moveto $lo
      $c(seq_right)     yview moveto $lo
   }
}
# DpDrag






###  SnScroll    
#
# routine used to scroll Seqname window
#
proc SnScroll {xy command pos {pages ""}} {
#############
    global c 
    global c3 

    set l_canvas [list $c(al_sn)]

    foreach canvas $l_canvas {
        eval $canvas ${xy}view $command $pos $pages
    }
}
# SnScroll


###  AlScroll    
#
# routine used to scroll alignment window
#
proc AlScroll {xy command pos {pages ""}} {
#############
    global c ;# r
    global c3 

    if {$xy == "x"} {
        set l_canvas [list $c(al) $c(al_nm)]
    } else {
        set l_canvas [list $c(al) $c3(al) $c(al_sn)]
    }
    foreach canvas $l_canvas {
        eval $canvas ${xy}view $command $pos $pages
    }
}
# AlScroll

#kl
###  DreiAlScroll    
#
# routine used to scroll alignment window
#
proc DreiAlScroll {xy command pos {pages ""}} {
#################
    global c3 ;# r

    if {$xy == "x"} {
        set l_canvas [list $c3(al) $c3(al_nm)]
    } else {
        set l_canvas [list $c3(al) $c(al_sn)]
    }
    foreach canvas $l_canvas {
        eval $canvas ${xy}view $command $pos $pages
    }
}
# DreiAlScroll


###  MiniAlScroll    
#
# routine used to scroll Minialignment window
#
proc MiniAlScroll {xy command pos {pages ""}} {
#################
    global g ;# r

    if {$xy == "x"} {
        set l_canvas [list $g(al) $g(al_nm)]
    } else {
        set l_canvas [list $g(al) $g(al_sn)]
    }
    foreach canvas $l_canvas {
        eval $canvas ${xy}view $command $pos $pages
    }
}
# MiniAlScroll
#kl


###   MIC_SetOpts   
#
#
proc MIC_SetOpts {var_name var_val} {
################
    global opts
    global mic_is_computed_and_valid

    set  opts(mic,$var_name)         $var_val

    if {$mic_is_computed_and_valid} {
        DeleteInfoDp
    }
    set  mic_is_computed_and_valid  0

}
# MIC_SetOpts




###   DpShowHide   
#
# Generic function for displaying/hiding dotplot items
# status: 0|1
#
proc DpShowHide {tag status} {
###############
    global c

    if {$status==0} {
        set state "hidden"
    } else {
        set state "normal"
    }
    $c(dp) itemconfigure $tag -state $state
}
# DpShowHide




###   MouseScroll   
# Generic mouse scrolling function for UNIX systems
# taken from http://koala.ilog.fr/colas/mouse-wheel-scroll/#tcl (Cola Nahaboo)
#
#
proc MouseScroll {bindtag} {
################

    # horizontal
    bind $bindtag <Button-5>         [list %W yview scroll 5 units]
    bind $bindtag <Button-4>         [list %W yview scroll -5 units]
    # vertikal with shift modifier
    bind $bindtag <Shift-Button-5>   [list %W xview scroll 5 units]
    bind $bindtag <Shift-Button-4>   [list %W xview scroll -5 units]
}
# MouseScroll



###   HideThisNtDpBps
#
# Hide all Dotplot Basepairs in which a particular nt participates
# Returns list of tags of hidden objects
#
proc HideThisNtDpBps {seqidx ntidx} {
####################
    global c

    set hiddenobjects {}
    foreach id [$c(dp) find withtag DpBp_seq_${seqidx}] {
    	foreach tag [$c(dp) gettags $id] {
    		if {[scan $tag "DpBp_seq_${seqidx}_nti_%d_ntj_%d" i j]==2} {
    			if {$i==$ntidx || $j==$ntidx} {
    				p::debug "hiding $tag"
    				DpShowHide $tag 0
    				lappend hiddenobjects $tag
    			}
    		}
    	}
    }
    return $hiddenobjects
}
# HideThisNtDpBps
