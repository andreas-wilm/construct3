##############################################################################
#
# circles.tcl - procedures to display rna-structures as circles-plot
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
#  CVS $Id: circles.tcl,v 1.14 2007-10-22 10:43:22 steger Exp $
#

##########
#
# Use CirclesFrontend (struct.tcl) your your convenience
# It wraps some things up
#
##########


package provide Circles 0.1


namespace eval Circles {


    # FIXME: Doku
    #        merge some parts with drawstruct ?
    #
    #
    variable CCircles
    variable alignlength    ;# set in Circles
    variable IMG_ZOOM ""    ;# filename for zoom image
    variable bplist; # list of basepair indices
    variable seqnt
    variable seqid

    ###   Init
    #
    #
    #
    proc Init {} {
    	########
    	variable CCircles
    	variable bplist
    	variable seqnt
    	variable seqid

    	set CCircles(Displ,bg_col)	    grey85;	# background of Circles plot
    	set CCircles(Displ,sequence)	normal;	# print/show the nts
    	set CCircles(SeqRanges)	        {0,0};	# print all nts (not only certain ranges)
    	set CCircles(Displ,numbers)	    normal;	# print/show the numbers
    	set CCircles(NumberSteps)	    10;    	# print numbering labels in steps of 10
    	set CCircles(Displ,mapping)	    normal;	# print/show mapping data
    	set CCircles(Displ,header)	    normal;	# print/show headers (title)
    	set CCircles(Zoom)		        1.;    	# start with no zoom
    	set CCircles(VIEW_X1)		    [expr {([lindex [wm maxsize .] 0]>[lindex [wm maxsize .] 1]) ? [expr [lindex [wm maxsize .] 1]*.6] : [expr [lindex [wm maxsize .] 0]*.6]}]
    	# size of plots in canvas windows; 400 pixels is fine
    	set CCircles(VIEW_Y1)		    $CCircles(VIEW_X1);	# dimish in case the plots are too large for your monitor
    	set CCircles(Font_Header)	    [DetermineFontSize [expr int($CCircles(VIEW_X1)/4)]]

    	set bplist {}
    	set seqnt "internal error"
    	set seqid "internal error"

    	### check if needed external commands exist
    	#
    	set external_cmds [list PrintPS ColorProb]

    	foreach cmd $external_cmds {
    		if {[lsearch -exact [info commands] $cmd]==-1} {
    			#auto_reset;# paranoia...no: bullshit all packages and their inits are lost
    			if {[auto_load $cmd]==0} {
    				puts "ERROR([namespace current]): cannot find required procedure \"$cmd\""
    				puts "ERROR([namespace current]): You must define it first"
    				return ERROR
    			}
    		}
    	}
    	return OK
    }
    # Init



    ###   Circles
    #
    #   Create a "circles"-like plot with two title lines from $LINE1 and $LINE2
    #   The X11 header is $HEADER.
    #   SEQ:  Array of nt chars from 1 to $SEQLEN
    #   PAIR: Array of int; nt i (==1 to $SEQLEN) pairs with nt j when $PAIR($I)==$J
    #   PROB: Array of real; base pairing probabilities
    #   If "show" is element of list $PltOptions then nt chars are written.
    #   Canvas size is $VIEW_X1,$VIEW_Y1
    #   For writing the title, font
    #   $Font_CirclesHeader==-adobe-times-bold-r-normal-*-*-80-*-*-*-*-*-*
    #   with fontsize 14 is used.
    #   If run from imatch, $imatch=1.
    #
    proc Circles {header line1 line2 sid s p w win_pos {doMapping 0} {line 0} {map 0} {imatch 0}} {
    	#####################################################################################
    	variable CCircles
    	variable alignlength
    	variable bplist
    	variable seqnt
    	variable seqid

    	upvar $p pair
    	upvar $w prob


    	set seqnt "$s"
    	set seqid "$sid"
    	set seqlen [string length $s]

    	for {set i 0} {$i<$seqlen} {incr i} {
    		set seq([expr {$i+1}]) [string index $s $i]
    	}
    	set alignlength $seqlen

    	set bplist {}
    	for {set i 1} {$i<=[string length $s]} {incr i} {
    		lappend bplist $pair($i)
    	}


    	if {$doMapping==1} {
    		upvar $map mapData
    	}


    	if {$imatch==1} {
    		set cw ".icircles"
    	} elseif {$imatch==2} {
    		set cw ".bcircles"
    	} else {
    		set cw ".circles"
    	}

    	OpenCirclesWindow $cw $win_pos $header $imatch pair

    	set pi         3.1415927
    	set 90        [expr {$pi/180.* 90.        }]
    	set 180       [expr {$pi/180.*180.        }]
    	set angle     [expr {$pi/180.*355./$seqlen}]
    	set x0        [expr {$CCircles(VIEW_X1)/2.}]; # center of circle
    	set y0        [expr {$CCircles(VIEW_Y1)/2.}]

    	set headersize 14;	# only a guess
    	set diminish    0
    	if {$CCircles(Displ,header)=="normal"} {
    		$cw.canvas create text $x0 0 \
    			-text $line1 -font $CCircles(Font_Header) \
    			-anchor n -justify center -tag header
    		$cw.canvas create text $x0 $headersize \
    			-text $line2 -font $CCircles(Font_Header) \
    			-anchor n -justify center -tag header
    		incr diminish [expr {2*$headersize}]
    	}

    	set ticklength [expr {int($x0/3.*sin($angle))}]
    	set ticklength [expr {($ticklength>1) ? $ticklength : 2}]
    	#  __
    	# \  \__
    	#  \    \__
    	#   \      \
    	#   /    __/  ^
    	#  /  __/     | arrow3; distance from the outside edge of the line to the trailing points
    	# /__/        v
    	#
    	#    <------> arrow1; distance along the line from the neck of the arrowhead to its tip
    	# <---------> arrow2; distance along the line from the trailing points of the arrowhead to the tip
    	#
    	if {$line>0} {
    		# set arrow1 [expr {([expr $ticklength*.8]> 8) ? [expr $ticklength*.8] : 8}]
    		# set arrow2 [expr {(      $ticklength    >10) ?       $ticklength     :10}]
    		# set arrow3 [expr {([expr $ticklength*.3]> 3) ? [expr $ticklength*.3] : 3}]
    		# set diminish [expr $diminish+$line*1.5*$arrow2]
    		set arrow1 [expr {$ticklength*1.8}]
    		set arrow2 [expr {$ticklength*1.8}]
    		set arrow3 [expr {$ticklength* .4}]
    		set arrowshape [list $arrow1 $arrow2 $arrow3]
    	  # set diminish [expr {$diminish+$line*2.*$arrow2}]
    		set diminish [expr {$diminish+$line*$arrow2/2.}]
    	}
    	if {$CCircles(Displ,sequence)=="normal"} {
    		set show 7
    	} else {
    		set show 5
    	}
    	set diminish   [expr {$diminish+$show*$ticklength}]
    	set radius     [expr {$x0-$diminish}]
    	set font       [DetermineFontSize [expr {$ticklength*20}]]

    	Kreisbogen $cw.canvas 0. 0. $x0 $y0 $radius 0.0 [expr {($seqlen-1)*$angle}] black
    	Kreisticks $cw.canvas 0. 0. $x0 $y0 $radius 0   [expr {$seqlen-1}] $angle $ticklength seq $font pair
    	for {set i 1} {$i<=$line} {incr i} {
    		set list [array names mapData "$i,\[0-9\]*"]
    		foreach j $list {
    			regexp "\[0-9\]*,(\[0-9\]*)" $j schrott k
    			incr k -1;
    			set dummy   [expr {$radius + $show*$ticklength + ($i-1)*$arrow2}]
    			set dummy1  [expr {$radius + $show*$ticklength +  $i   *$arrow2}]
    			set dummy_x [expr {cos($k*$angle)}]
    			set dummy_y [expr {sin($k*$angle)}]
    			$cw.canvas create line [expr {$x0+$dummy *$dummy_x}] \
    				[expr {$y0+$dummy *$dummy_y}] \
    				[expr {$x0+$dummy1*$dummy_x}] \
    				[expr {$y0+$dummy1*$dummy_y}] \
    				-fill       $mapData($i,color) \
    				-width      1.0 \
    				-arrow      first  \
    				-arrowshape $arrowshape \
    				-tag        mapping
    		}
    	}
    	if {$imatch != 2} {
    		set dummyseqlen $seqlen
    	} else {
    		set dummyseqlen [expr {$seqlen*2}]
    	}
    	for {set n 1} {$n<=$dummyseqlen} {incr n} {
    		if {$imatch != 2} {
    			set i $n
    			set j $pair($i)
    		} else {
    			if {$n > $seqlen} {
    				set m [expr {$n-$seqlen}]
    			} else {
    				set m $n
    			}
    			if {$m<$pair($n) || $pair($n)==0} {
    				set i $m
    				set j $pair($n)
    			} else {
    				set j $m
    				set i $pair($n)
    			}
    		}
    		if {$pair($i)>0 && $i<$j} {
    			#	set j $pair($i)
    			set alphaI [expr {$angle*($i-1)}]
    			set alphaJ [expr {$angle*($j-1)}]
    			if {$prob($i)==0.} {
    				set color [format "#000000000000"];   # schwarz
    			} else {
    				set color [ColorProb $prob($i)]
    			}
    			if {[expr {abs($180-$alphaJ+$alphaI)}]<[expr {$pi/180.*2.}]} {
    				# puts "TMP_DEBUG(): Straight line between $i and $j"
    				$cw.canvas create line [expr {$x0+$radius*cos($alphaI)}] [expr {$CCircles(VIEW_Y1)-$y0+$radius*sin($alphaI)}] \
    					[expr {$x0+$radius*cos($alphaJ)}] [expr {$CCircles(VIEW_Y1)-$y0+$radius*sin($alphaJ)}] \
    					-fill $color -width 2. -tag line_${i}_$pair($i)
    				$cw.canvas addtag line withtag line_${i}_$pair($i)
    				$cw.canvas bind line_${i}_$pair($i) <Button-1> "[namespace code ChangeBpState] %x %y %W"
    			} else {
    				set gammaI [expr {$alphaI-$90}]
    				set gammaJ [expr {$alphaJ+$90}]
    				set x1     [expr {$radius*cos($alphaI)}]
    				set y1     [expr {$radius*sin($alphaI)}]
    				set x2     [expr {$radius*cos($alphaJ)}]
    				set y2     [expr {$radius*sin($alphaJ)}]
    				set alpha  [expr {$alphaI+$90}]
    				set beta   [expr {$alphaJ-$90}]
    				set nenner [expr {cos($alpha)-cos($beta)}]
    				if {$nenner<0.001} {
    					set nenner  [expr {sin($alpha)-sin($beta)}]
    					set radius2 [expr {($y2-$y1)/$nenner}]
    				} else {
    					set radius2 [expr {($x2-$x1)/$nenner}]
    				}
    				set x02 [expr {$x1+$radius2*cos($alpha)}]
    				set y02 [expr {$y1+$radius2*sin($alpha)}]
    				if {[expr {$alphaJ-$alphaI}]>$180} {
    					Kreisbogen $cw.canvas $x02 $y02 $x0 $y0 $radius2 $gammaI $gammaJ $color $i $pair($i)
    				} else {
    					Kreisbogen $cw.canvas $x02 $y02 $x0 $y0 $radius2 $gammaJ $gammaI $color $i $pair($i)
    				}
    			}
    		}
    	}
    }
    # Circles




    ###   Kreisbogen
    #
    #   Plot an arc with RADIUS from ANGLE0 to ANGLE9 (in rad)
    #   at center X0+X0FF,Y0+YOFF using COLOR in pathname W
    #
    proc Kreisbogen {w x0 y0 xoff yoff radius angle0 angle9 color {ipair 0} {jpair 0}} {
    	##############################################################################

    	set x0 [expr {$x0+$xoff}]
    	set y0 [expr {$y0+$yoff}]
    	set pi 3.1415927
    	set deltaAngle [WDomain [expr {($angle9 - $angle0)*180./$pi}]]

    	if {$deltaAngle==0.} {
    		set grad 360
    	} else {
    		set grad [expr {int($deltaAngle)}]
    	}
    	set coords ""
    	for {set i 0} {$i<=$grad} {incr i} {
    		set ii [expr {$i*$pi/180.}]
    		lappend coords [expr {$x0+$radius*cos($ii+$angle0)}] [expr {$y0+$radius*sin($ii+$angle0)}]
    	}
    	lappend coords [expr {$x0+$radius*cos(    $angle9)}] [expr {$y0+$radius*sin(    $angle9)}]
    	if {$ipair > 0 && $jpair > 0 } {
    		eval $w create line $coords -fill $color -width 2 -tag line_${ipair}_${jpair}
    		$w addtag line withtag line_${ipair}_${jpair}
    		$w bind line_${ipair}_${jpair} <Button-1> "[namespace code ChangeBpState] %x %y %W"
    	} else {
    		eval $w create line $coords -fill $color -width 2;	# -smooth true -splinesteps 5
    	}
    }
    # Kreisbogen




    ###   Kreisticks
    #
    #   Plot polar ticks with RADIUS from BEGIN to END with ANGLE (in rad)
    #   at center X0+X0FF,Y0+YOFF in pathname W
    #   If $CCircles(Displ,sequence)=="normal" (not "hidden") then nt chars
    #   are written in between indices #i,#j from the list $CCircles(SeqRanges)
    #   If $CCircles(Displ,numbers)=="normal" (not "hidden") then ticks are numbered
    #   in steps of $CCircles(NumberSteps); otherwise a step of 10 is used
    #
    proc Kreisticks {w x0 y0 xoff yoff radius begin end angle length s font p} {
    	######################################################################

    	variable CCircles
    	#
    	upvar $s seq
    	upvar $p pair

    	set y0 [expr {$yoff-$y0}]
    	set x0 [expr {$x0+$xoff}]
    	set coords ""
    	for {set i $begin} {$i<=$end} {incr i} {
    		set dummy_x [expr {cos($i*$angle)}]
    		set dummy_y [expr {sin($i*$angle)}]
    		$w create line  [expr {$x0+$radius*$dummy_x}]           [expr {$y0+$radius*$dummy_y}]	\
    			[expr {$x0+($radius+$length)*$dummy_x}] [expr {$y0+($radius+$length)*$dummy_y}] \
    			-fill black -width 1.
    	}
    	set distance 0
    	if {$CCircles(Displ,sequence)=="normal"} {
    		set distance 3
    		foreach val $CCircles(SeqRanges) {
    			set dummy [split $val ","]
    			set beg1  [expr {[lindex $dummy 0]-1}]; if {$beg1==-1} {set beg1 0   }
    			set end1  [expr {[lindex $dummy 1]-1}]; if {$end1==-1} {set end1 $end}
    			for {set i $beg1} {$i<=$end1} {incr i} {
    				if {$seq([expr {1+$i}])!="X" && $seq([expr {1+$i}])!="x"} {
    					set dummy_x [expr {cos($i*$angle)}]
    					set dummy_y [expr {sin($i*$angle)}]
    					set j [expr {1+$i}]
    					$w create text [expr {$x0+($radius+2*$length)*$dummy_x}] \
    						[expr {$y0+($radius+2*$length)*$dummy_y}] \
    						-text "$seq($j)" \
    						-anchor center -fill black -font $font -tag sequence
    				}
    			}
    		}
    	}
    	if {$CCircles(Displ,numbers)=="normal"} {
    		set step $CCircles(NumberSteps)
    		$w create line	[expr {$x0+($radius+$distance*$length)}]     [expr {$y0}] \
    			[expr {$x0+($radius+($distance+2)*$length)}] [expr {$y0}] \
    			-fill black  -width 1.0  -tag numbers
    		$w create text	[expr {$x0+($radius+($distance+3)*$length)}] [expr {$y0}] \
    			-text 1  -anchor center  -fill black  -tag numbers  -font $font
    		for {set i [expr {$step-1}]} {$i<=$end} {set i [expr {$i+$step}]} {
    			set dummy_x [expr {cos($i*$angle)}]
    			set dummy_y [expr {sin($i*$angle)}]
    			$w create line	[expr {$x0+($radius+ $distance*$length)*$dummy_x}] \
    				[expr {$y0+($radius+ $distance   *$length)*$dummy_y}] \
    				[expr {$x0+($radius+($distance+2)*$length)*$dummy_x}] \
    				[expr {$y0+($radius+($distance+2)*$length)*$dummy_y}] \
    				-fill black  -width  1.0  -tag numbers
    			$w create text	[expr {$x0+($radius+($distance+3)*$length)*$dummy_x}] \
    				[expr {$y0+($radius+($distance+3)*$length)*$dummy_y}] \
    				-text [expr {1+$i}] \
    				-anchor center  -fill black  -tag numbers  -font $font
    		}
    	}
    }
    # Kreisticks


    ###   ChangeBpState
    #
    #
    proc ChangeBpState {x y w} {
    	######################

    	set x [$w canvasx $x]
    	set y [$w canvasy $y]

    	set dummy [$w gettags [$w find closest $x $y]]

    	regexp {line_([0-9]+_[0-9]+) } $dummy match nt
    	regexp {(line_[0-9]+_[0-9]+) } $dummy match nt_tag

    	if {[info exists nt_tag] && [info exists nt]} {
    		if { [$w itemcget $nt_tag -dash] == "" } {
    			$w itemconfigure line_$nt -dash {6}
    			$w dtag line_$nt line
    		} else {
    			$w itemconfigure line_$nt -dash {}
    			$w addtag line withtag line_$nt
    		}
    	}
    }
    # ChangeBpState




    ###   WDomain
    #
    #   Takes any angle in degrees and returns the equivalent angle
    #   in as an angle between 0 and 360.
    #
    proc WDomain {angle} {
    	######################
    	set mult [expr {int(abs($angle/360.0))}]
    	if {$angle >= 0.0} {
    		return [expr {$angle -  $mult       *360.0}]
    	} else {
    		return [expr {$angle + ($mult + 1.0)*360.0}]
    	}
    }
    # WDomain




    ###   PrintPlots
    #
    #	print plots
    #
    proc PrintPlots {w w1 w2 datei ext1 ext2 save_print} {
    	################################################

    	# cursor disable blue

    	if {$w1!="none"} {
    		itemChangeState $w.$w1 print normal
    		PrintPS $w.$w1 [format "%s_%s.eps" $datei $ext1] $save_print
    		itemChangeState $w.$w1 print hidden
    	}
    	if {$w2!="none"} {
    		itemChangeState $w.$w2 print normal
    		PrintPS $w.$w2 [format "%s_%s.eps" $datei $ext2] $save_print
    		itemChangeState $w.$w2 print hidden
    	}

    	# cursor enable
    }
    # PrintPlots




    ###   itemChangeState
    #
    #
    proc itemChangeState {w tag state} {
    	##############################
    	foreach id [$w find withtag $tag] {
    		$w itemconfigure $id -state $state
    	}
    }
    # itemChangeState




    ###   OpenCirclesWindow
    #
    #
    #
    proc OpenCirclesWindow {w win_pos datei imatch p} {
    	############################################
    	variable CCircles
    	variable IMG_ZOOM

    	upvar $p pair

    	if {[winfo exists $w]} {
    		destroy $w
    	}
    	toplevel	$w
    	wm title	$w "View of [file tail $datei]"
    	wm iconname	$w          [file tail $datei]
    	wm geometry  $w $win_pos

    	canvas	$w.canvas -width  $CCircles(VIEW_X1)	\
    		-height         $CCircles(VIEW_Y1)   \
    		-xscrollcommand "$w.hscroll set" \
    		-yscrollcommand "$w.vscroll set" \
    		-scrollregion   "0 0 $CCircles(VIEW_X1) $CCircles(VIEW_Y1)"  \
    		-relief         raised  		 \
    		-background     $CCircles(Displ,bg_col) \
    		-borderwidth    1


    	### Buttons
    	#
    	#

    	frame	$w.buttons	-bd 2 -relief raised
    	#   button	$w.buttons.dismiss	-text Dismiss		-command "destroy $w"
    	#   button	$w.buttons.save		-text Save      	-command "PrintPlots $w canvas none [file rootname $datei] circles x save"
    	#   button	$w.buttons.print	-text Print     	-command "PrintPlots $w canvas none [file rootname $datei] circles x print"
    	#   button	$w.buttons.printcolor	-text PrintColor     	-command "PrintPlots $w canvas none [file rootname $datei] circles x printcolor"
    	if {$IMG_ZOOM!=""} {
    		button	$w.buttons.zoom -image $IMG_ZOOM
    	} else {
    		button	$w.buttons.zoom -text "Zoom"
    	}


    	### Menu
    	#
    	#
    	frame $w.f_menu	-bd 2 -relief raised


    	### File menu
    	#
    	#
    	menubutton	$w.f_menu.file	-text File \
    		-menu $w.f_menu.file.menu

    	set m	$w.f_menu.file.menu

    	menu $m

    	if {$imatch != 2} {
    		#       $m add command	-label " Write SQU file " \
    						   #		                -command "WriteSquFile $imatch [file rootname $datei]"
    		#		$m add sep
    	}
    	$m add cascade	-label "Save... " -menu $m.save
    	$m add cascade	-label "Print..." -menu $m.print
    	$m add sep
    	$m add command	-label " Close " -command "destroy $w"


    	###   Submenu Save

    	menu $m.save
    	$m.save add command	-label "Consensus" -command "[namespace code SaveCons] cs-consensus"

    	$m.save add cascade	-label "Connect..." -menu $m.ct 
    	menu  $m.ct
    	$m.ct add command -label "Zuker" -command "[namespace code SaveCons] connect zuker"
    	$m.ct add command -label "Gcg" -command "[namespace code SaveCons] connect gcg"
    	$m.ct add command -label "Rnaviz" -command "[namespace code SaveCons] connect rnaviz"


    	
    	# cannot store triples !
    	#if {$imatch!=1 && $imatch!=2} {
    	#	for {set i 0} {$i<3} {incr i} {
    	#		$m.ct entryconfigure $i -state disabled
    	#	}
    	#}
    	$m.save add command	-label "Rnaml" -command "[namespace code SaveCons] rnaml"
    	$m.save add command	-label "Stockholm" -command "[namespace code SaveCons] stockholm"

    	
    	###   Submenu print

    	menu $m.print
    	$m.print add command \
    		-label "to printer " \
    		-underline 3 \
    		-accelerator "Ctrl+p"	\
    		-command "[namespace code PrintPlots] $w canvas none [file tail [file rootname $datei]] circles x print"
    	bind $w <Control-p> "[namespace code PrintPlots] $w canvas none [file tail [file rootname $datei]] circles x print"
    	$m.print add command \
    		-label "to color printer " \
    		-command "[namespace code PrintPlots] $w canvas none [file tail [file rootname $datei]] circles x printcolor"
    	$m.print add command \
    		-label "to file " \
    		-command "[namespace code PrintPlots] $w canvas none [file tail [file rootname $datei]] circles x save"
    	$m.print add command \
    		-label "to screen " \
    		-command "[namespace code PrintPlots] $w canvas none [file tail [file rootname $datei]] circles x screen"


    	### Display menu
    	#
    	#
    	menubutton	$w.f_menu.displ	-text "Display" \
    		-menu $w.f_menu.displ.menu

    	set m	$w.f_menu.displ.menu
    	menu		$m

    	foreach displmode {sequence numbers mapping header} {
    		$m add check    -label [format " %s%s " [string toupper [string index $displmode 0]] [string range $displmode 1 end]]	\
    			-variable Circles::CCircles(Displ,$displmode) \
    			-onvalue  normal \
    			-offvalue hidden \
    			-command  "[namespace  code CirclesItemConfig] $w.canvas $displmode"
    	}
    	$m add sep
    	$m add cascade	 -label " Background " \
    		-menu $w.f_menu.displ.menu.bg

    	menu $m.bg
    	$m.bg add radio -label     " white " \
    		-variable  Circles::CCircles(Displ,bg_col) \
    		-value     white \
    		-command   "[namespace code CirclesItemConfig] $w.canvas white"

    	$m.bg add radio -label     " grey " \
    		-variable  Circles::CCircles(Displ,bg_col) \
    		-value     grey85 \
    		-command   "[namespace code CirclesItemConfig] $w.canvas grey85"


    	### Options menu
    	#
    	#
    	menubutton $w.f_menu.opt	-text "Options" \
    		-menu $w.f_menu.opt.menu

    	set m	$w.f_menu.opt.menu
    	menu		$m
    	$m add command	-label		" Get sequence ranges "	\
    		-command	"[namespace code StructSeqRanges] $m"
    	$m add command	-label		[format " Numbering step = %d" $CCircles(NumberSteps)]	\
    		-command	"[namespace code StructNumSteps] $m"

    	####################################

    	scrollbar $w.hscroll -orient horiz -command "$w.canvas xview"
    	scrollbar $w.vscroll -orient vert  -command "$w.canvas yview"

    	grid                 $w.canvas		-column 1 -row 2		-sticky wens
    	grid columnconfigure $w 1 -weight 1
    	grid rowconfigure    $w 2 -weight 1

    	grid 		$w.f_menu		-column 1 -row 1 -columnspan 2	-sticky we
    	grid			$w.f_menu.file		-column 1 -row 1
    	grid			$w.f_menu.displ		-column 2 -row 1
    	grid			$w.f_menu.opt		-column 3 -row 1

    	grid			$w.buttons		-column 1 -row 4 -columnspan 2	-sticky we
    	#   grid			$w.buttons.dismiss	-column 1 -row 1 				-padx 32
    	#   grid			$w.buttons.save		-column 2 -row 1 				-padx 16
    	#   grid			$w.buttons.print	-column 3 -row 1 				-padx 16
    	#   grid			$w.buttons.printcolor	-column 4 -row 1 				-padx 16
    	grid			$w.buttons.zoom		-column 5 -row 1 				-padx 16

    	grid			$w.hscroll		-column 1 -row 3 -sticky ew
    	grid			$w.vscroll		-column 2 -row 2 -sticky wns


    	bind $w.canvas	<3>		"$w.canvas scan mark   %x %y"
    	bind $w.canvas	<B3-Motion>	"$w.canvas scan dragto %x %y"

    	bind $w.buttons.zoom	<Button-3>	"[namespace code ZoomCircleCanvas] $w.canvas 0.5"
    	bind $w.buttons.zoom	<Button-1>	"[namespace code ZoomCircleCanvas] $w.canvas 2.0"

    	$w.canvas create rect -1600 -1280 1600 1280 -fill $CCircles(Displ,bg_col) -tag Cbackgnd

    }
    # OpenCirclesWindow




    ###   SaveCons
    #
    # save consensus
    #
    proc SaveCons {structformat {ctsubformat "zuker"}} {
    	#############
    	variable seqnt
    	variable seqid
    	variable bplist

    	SaveStructFrontend $seqid $seqnt $bplist \
    		$structformat $ctsubformat
    }
    # SaveCons



    ###   StructSeqRanges
    #
    #
    #
    proc StructSeqRanges {menuW} {
    	########################
    	variable CCircles
    	variable alignlength

    	set f .seqRanges
    	if {[winfo exists $f]==0} {
    		set oldFocus [focus]
    		set row 0
    		toplevel    $f
    		wm title    $f SeqRanges
    		wm iconname $f SeqRanges
    		set    txt "Set ranges for sequence labeling\n"
    		append txt "\nGive two comma-separated numbers for each range.\n"
    		append txt "Individual ranges are separated by blanks.\n"
    		append txt "\nF.e.: a,b c,d\n1 =< a < b =< $alignlength\n"
    		append txt "0,0 is equivalent to 1,$alignlength"
    		label $f.message	-wraplength 5i	\
    			-justify left -text "$txt"
    		grid	  $f.message	-row $row; incr row
    		label	  $f.leer$row	-text " "
    		grid	  $f.leer$row	-row $row; incr row
    		frame	  $f.a 		-bd 2
    		label	  $f.a.label1	-text "Ranges = "
    		entry	  $f.a.entry	-relief	sunken	\
    			-width	10	\
    			-xscrollcommand "$f.a.xscroll set"
    		scrollbar	  $f.a.xscroll	-relief	sunken	\
    			-orient horiz	\
    			-command "$f.a.entry xview"
    		grid	  $f.a.label1	-column 0 -row 0 -sticky w
    		grid	  $f.a.entry	-column 1 -row 0 -sticky we
    		grid	  $f.a.xscroll	-column 1 -row 1 -sticky we
    		grid	  $f.a		-row $row; incr row
    		bind	  $f.a.entry <Return> "GetSeqRanges $f.a.entry $menuW"
    		eval set x [format "\$%s" CCircles(SeqRanges)]
    		$f.a.entry insert 0 [format "%s" $x]
    		label	  $f.leer$row	-text " "
    		grid        $f.leer$row	-row $row; incr row
    		frame       $f.buttons
    		grid        $f.buttons	-row $row; incr row
    		button      $f.buttons.dismiss	-text Dismiss	-command "Circles::CancelEntry $f $oldFocus"
    		button      $f.buttons.accept	-text Return	-command "Circles::GetSeqRanges $f.a.entry $menuW; Circles::CancelEntry $f $oldFocus"
    		grid        $f.buttons.dismiss	-column 0 -row 0
    		grid        $f.buttons.accept	-column 1 -row 0
    	}
    	set pos [winfo pointerxy $menuW]
    	wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
    	raise $f
    	grab  $f
    	focus $f
    }
    # StructSeqRanges




    ###   GetSeqRanges
    #
    #
    #
    proc GetSeqRanges {entryW menuW} {
    	###########################
    	variable CCircles
    	variable alignlength

    	set x [$entryW get]
    	foreach val $x {
    		set dummy [split $val ","]
    		set beg  [expr {[lindex $dummy 0]}]; set beg [expr {($beg<1) ? 1 : $beg}]
    		set beg  [expr {($beg>$alignlength) ? $alignlength: $beg}]
    		set end  [expr {[lindex $dummy 1]}]; set end [expr {($end<1) ? $alignlength: $end}]
    		set end [expr {($end>$alignlength) ? $alignlength: $end}]
    		if {$beg>$end} {
    			set dummy $beg
    			set beg   $end
    			set end   $dummy
    		}
    		lappend y "$beg,$end"
    	}
    	$entryW delete 0 end
    	$entryW insert 0 [format "%s" $y]
    	set CCircles(SeqRanges) $y

    	NotifyReload
    }
    # GetSeqRanges




    ###   StructNumSteps
    #
    #
    proc StructNumSteps {menuW} {
    	######################
    	variable CCircles
    	variable alignlength

    	set f .numSteps
    	if {[winfo exists $f]==0} {
    		set oldFocus [focus]
    		set row 0
    		toplevel    $f
    		wm title    $f NumSteps
    		wm iconname $f NumSteps
    		label   $f.message    -wraplength 5i  \
    			-justify left	\
    			-text "Set numbering steps\n\n1 nts =< step size =< $alignlength nts"
    		grid	  $f.message	-row $row; incr row
    		label	  $f.leer$row	-text " "
    		grid	  $f.leer$row	-row $row; incr row
    		frame	  $f.a 		-bd 2
    		label	  $f.a.label1	-text "Step size = "
    		entry	  $f.a.entry	-relief sunken -width 4
    		label	  $f.a.label2	-text " nts"
    		grid	  $f.a.label1	-column 0 -row 0 -sticky w
    		grid	  $f.a.entry	-column 1 -row 0 -sticky e
    		grid	  $f.a.label2	-column 2 -row 0 -sticky w
    		grid	  $f.a		-row $row; incr row
    		bind	  $f.a.entry <Return> "Circles::GetNumSteps $f.a.entry $menuW"
    		eval set x [format "\$%s" CCircles(NumberSteps)]
    		$f.a.entry insert 0 [format "%d" $x]
    		label	  $f.leer$row	-text " "
    		grid        $f.leer$row	-row $row; incr row
    		frame       $f.buttons
    		grid        $f.buttons	-row $row; incr row
    		button      $f.buttons.dismiss	-text Dismiss	-command "Circles::CancelEntry $f $oldFocus"
    		button      $f.buttons.accept	-text Return	-command "Circles::GetNumSteps $f.a.entry $menuW; Circles::CancelEntry  $f $oldFocus"
    		grid        $f.buttons.dismiss	-column 0 -row 0
    		grid        $f.buttons.accept	-column 1 -row 0
    	}
    	set pos [winfo pointerxy $menuW]
    	wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
    	raise $f
    	grab  $f
    	focus $f
    }
    # StructNumSteps




    ###   GetNumSteps
    #
    #
    proc GetNumSteps {entryW menuW} {
    	####################### ###
    	variable CCircles
    	variable alignlength

    	set x [expr {int([$entryW get])}]
    	set x [expr {($x>$alignlength) ? $alignlength : $x}]
    	set x [expr {($x<1           ) ? 1            : $x}]
    	$entryW delete 0 end
    	$entryW insert 0 [format "%d" $x]
    	set CCircles(NumberSteps) $x
    	$menuW entryconfigure 2 -label [format " Numbering step = %d" $CCircles(NumberSteps)]

    	NotifyReload
    }
    # GetNumSteps



    ###   ZoomCircleCanvas
    #
    #
    proc ZoomCircleCanvas {w factor} {
    	############################
    	variable CCircles

    	if {[expr {$CCircles(Zoom)*$factor}] > 0.45} {
    		set CCircles(Zoom) [expr {$CCircles(Zoom)*$factor}]
    		set dummy [$w xview]; set x0 [lindex $dummy 0]; set x1 [lindex $dummy 1];	# remember the visible scrollregion
    		set dummy [$w yview]; set y0 [lindex $dummy 0]; set y1 [lindex $dummy 1]

    		foreach id [$w find all] {;							# scale the canvas
    			$w scale $id 0. 0. $factor $factor
    		}

    		$w config -scrollregion [list 0. 0. [expr {$CCircles(VIEW_X1)*$CCircles(Zoom)}]	\
    									 [expr {$CCircles(VIEW_Y1)*$CCircles(Zoom)}]];	# new size of canvas

    		# ($a0+$a1)/2.      = center of old scrollregion
    		# ($a1-$a0)/$factor = span of new visible scrollregion
    		$w xview moveto [expr {($x0+$x1)/2.-($x1-$x0)/2./$factor}]; # left border of visible part
    		$w yview moveto [expr {($y0+$y1)/2.-($y1-$y0)/2./$factor}]; # top  border of visible part
    	}
    }
    # ZoomCircleCanvas




    ###   CirclesItemConfig
    #
    proc CirclesItemConfig {w tag} {
    	##########################
    	variable CCircles

    	if {$tag == "white" || $tag == "grey85"} {
    		$w config -background $CCircles(Displ,bg_col)
    		$w itemconfig Cbackgnd -fill $CCircles(Displ,bg_col)
    	} else {
    		$w itemconfig $tag -state $CCircles(Displ,$tag)
    	}
    }
    # CirclesItemConfig



    ### WriteStructure
    #
    #
    proc WriteStructure {w} {
    	##################
    	variable alignlength

    	set filename [tk_getSaveFile        \
    					  -title        "Select structure file"]

    	update idletasks

    	if {$filename == ""} {
    		puts "Save dialog aborted: writing to stdout"
    		set fid stdout
    	} else {
    		set fid [open $filename w+]
    	}



    	for {set i 1} {$i <= $alignlength} {incr i} {
    		set pair($i) 0
    	}

    	puts -nonewline $fid "> "
    	set Header [$w find withtag header]
    	foreach i $Header {
    		puts -nonewline $fid  [$w itemcget $i -text]
    	}
    	puts $fid " "

    	set AllId [$w find withtag line]
    	foreach id $AllId {
    		set Tags [$w gettags $id]
    		regexp {line_([0-9]+)_([0-9]+) } $Tags match nti ntj
    		set pair($ntj) $nti
    		set pair($nti) $ntj
    	}
    	set structKlammer ""
    	for {set i 1} {$i <= $alignlength} {incr i} {
    		if {$pair($i)>0} {
    			if {$pair($i)>$i} {
    				set structKlammer "${structKlammer}("
    			} else {
    				set structKlammer "${structKlammer})"
    			}
    		} else {
    			set structKlammer "${structKlammer}."
    		}
    	}
    	puts -nonewline $fid $structKlammer

    	if {$filename != ""} {close $fid}
    }
    # WriteStructure





    ###   DetermineFontSize
    #
    #
    #
    proc DetermineFontSize {fontsize} {
    	#############################

    	set fontsize [expr {int($fontsize)}]
    	set server [winfo server .]


    	if {[lsearch -exact $server DECWINDOWS] > -1 } {;	# VMS interpoliert nicht
    		if {[lsearch -exact $server X11R0] > -1 } {
    			if {$fontsize < 90} {
    				set fontsize	80
    			} elseif {$fontsize < 110} {
    				set fontsize 100
    			} elseif {$fontsize < 130} {
    				set fontsize 120
    			} elseif {$fontsize < 160} {
    				set fontsize 140
    			} elseif {$fontsize < 210} {
    				set fontsize 180
    			} else {
    				set fontsize 240
    			}
    		} else {
    			if {$fontsize < 25} {
    				set fontsize 20
    			} elseif {$fontsize < 55} {
    				set fontsize 50
    			} elseif {$fontsize < 65} {
    				set fontsize 60
    			} elseif {$fontsize < 75} {
    				set fontsize 70
    			} elseif {$fontsize < 85} {
    				set fontsize 80
    			} elseif {$fontsize < 95} {
    				set fontsize 90
    			} elseif {$fontsize < 105} {
    				set fontsize 100
    			} elseif {$fontsize < 115} {
    				set fontsize 110
    			} elseif {$fontsize < 125} {
    				set fontsize 120
    			} elseif {$fontsize < 135} {
    				set fontsize 130
    			} elseif {$fontsize < 145} {
    				set fontsize 140
    			} elseif {$fontsize < 155} {
    				set fontsize 150
    			} elseif {$fontsize < 165} {
    				set fontsize 160
    			} elseif {$fontsize < 175} {
    				set fontsize 170
    			} elseif {$fontsize < 185} {
    				set fontsize 180
    			} elseif {$fontsize < 195} {
    				set fontsize 190
    			} elseif {$fontsize < 215} {
    				set fontsize 200
    			} elseif {$fontsize < 235} {
    				set fontsize 230
    			} elseif {$fontsize < 245} {
    				set fontsize 240
    			} elseif {$fontsize < 265} {
    				set fontsize 250
    			} elseif {$fontsize < 290} {
    				set fontsize 280
    			} elseif {$fontsize < 325} {
    				set fontsize 300
    			} elseif {$fontsize < 355} {
    				set fontsize 350
    			} elseif {$fontsize < 380} {
    				set fontsize 360
    			} elseif {$fontsize < 440} {
    				set fontsize 400
    			} elseif {$fontsize < 600} {
    				set fontsize 480
    			} else {
    				set fontsize 720
    			}
    		}
    	} else {
    		if {$fontsize < 10} {set fontsize 10};		# that's it
    	}

    	return [format "-adobe-courier-bold-r-*-*-*-%d-*-*-*-*-*-*" $fontsize]

    }
    # DetermineFontSize


    ###   CancelEntry
    #
    #
    proc CancelEntry {f oldFocus} {
    	#########################
    	
    	focus        $oldFocus
    	grab release $f
    	destroy      $f

    }
    # CancelEntry



    ###   NotifyReload
    #
    #
    proc NotifyReload {} {
    	################
    	
    	set title "Note"
    	set msg   "Close the Circles window and open it again to apply made changes !"
    	tk_messageBox  -title  "$title" -message "$msg" \
    		-icon error  -type  ok
    }
    # NotifyReload


}
# namespace eval circles
