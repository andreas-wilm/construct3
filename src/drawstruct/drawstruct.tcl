##############################################################################
#
# drawstruct.tcl - routines for squggles like display of rna structures
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
##############################################################################

#
#  CVS $Id: drawstruct.tcl,v 1.23 2007-10-22 10:43:23 steger Exp $
#

package provide DrawStruct 0.9

namespace eval DrawStruct {

    ### Namespace wide variables
    #
    # initialized once on namespace evaluation
    #
    #
    variable debug_level 0
    #
    variable DStr
    variable DS_DisplOpt
    #
    #
    variable print_pswidth
    #
    variable tmp_dir
    #
    variable zoom_img      ""
    variable reset_img     ""


    # load c-core
    package require drawstructcore



    ###   Init
    #
    # Called only once after "package require"
    #
    proc Init {} {
		########
        variable print_pswidth
        variable tmp_dir


        set tmp_dir /tmp
        set print_pswidth 16

        ### check if needed external commands exist
        #
        set external_cmds [list PrintPS ColorProb MouseScroll]

        foreach cmd $external_cmds {
            if {[lsearch -exact [info commands] $cmd]==-1} {
                auto_reset;# paranoia
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


    ###   DrawStructure
    #
    # Core Routine
    #
    # Arguments:
    #   seq : sequence nt's
    #   l_basepairs : list of basepairs
    #   optimalstruct : 1 if provided structure is optimal, 0 if suboptimal
    #   win_pos : window geometry
    #   args:
    #       following optionals:
    #            -header : list of sequence header displayed in canvas
    #            -prob : list of basepairing probability
    #            -projname : project name
    #
    #
    #
    proc DrawStructure {seq  l_basepairs  optimal_struct  win_pos  args} {
		###############################################################
         variable DStr
         variable debug_level

       if {$debug_level > 0} {
            puts "DrawStruct::DrawStructure:  seq:            $seq"
            puts "DrawStruct::DrawStructure:  l_basepairs:    $l_basepairs"
            puts "DrawStruct::DrawStructure:  optimal_struct: $optimal_struct"
            puts "DrawStruct::DrawStructure:  args:           $args"
        }


        # find a new unique window number
        #
        set m 0
        while {[winfo exists .w_drawstruct_$m]} {
			incr m
		}


        WinSetup $m

        set DStr($m,l_basepairs) $l_basepairs
        set DStr($m,seqlen) [string length $seq]
        set DStr($m,seq) $seq
        wm geometry $DStr($m,win) $win_pos

        ParseArgs $m $DStr($m,seqlen) $args



        canvas $DStr($m,can)  -xscrollcommand "$DStr($m,win).hscroll set" \
                              -yscrollcommand "$DStr($m,win).vscroll set" \
                              -closeenough 8.0 \
                              -relief raised \
                              -borderwidth 1

        grid $DStr($m,can)    -column 1 -row 2 -sticky wens
        grid columnconfigure  $DStr($m,win) 1 -weight 1
        grid rowconfigure     $DStr($m,win) 2 -weight 1

        frame $DStr($m,win).f_menu -bd 2 -relief raised
        frame $DStr($m,win).f_tools -bd 2 -relief raised


        ### File menu
        #
        CreateFileMenu $m $optimal_struct


        ### Display menu
        #
        CreateDisplayMenu $m


        ### Buttons
        #
        CreateButtons $m


       scrollbar $DStr($m,win).hscroll -orient horiz -command "$DStr($m,can) xview"
       scrollbar $DStr($m,win).vscroll -orient vert  -command "$DStr($m,can) yview"

       grid $DStr($m,win).f_menu           -column 1 -row 1 -columnspan 2 -sticky we
       grid $DStr($m,win).f_menu.file      -column 1 -row 1 -sticky w
       grid $DStr($m,win).f_menu.displ     -column 2 -row 1 -sticky e

       grid $DStr($m,win).f_tools          -column 3 -row 1 -rowspan 3 -sticky ns
       grid $DStr($m,win).f_tools.b_zoom   -column 3 -row 1
       grid $DStr($m,win).f_tools.b_reset  -column 3 -row 2

       grid $DStr($m,win).hscroll          -column 1 -row 3 -sticky ew
       grid $DStr($m,win).vscroll          -column 2 -row 2 -sticky wns

       CreateStructure $m

       MouseScroll $DStr($m,can)

    }
	# DrawStructure



    ###   ParseArgs
    #
    # Process provided arguments
    #
    #
    proc ParseArgs {win_no seqlen args} {
		##############################
        variable DStr

        # örks... args is now a list of lists of lists of
        set args [lindex $args 0]

        if {[llength $args] > 0} {

            # Extract header
            #
			if {[set headerstart [lsearch $args "-header"]] > -1} {
				set l_headers [lindex $args [expr {$headerstart+1}]]
				set DStr($win_no,headline,1) [lindex $l_headers 0]
				set DStr($win_no,headline,2) [lindex $l_headers 1]
				set DStr($win_no,headline,3) [lindex $l_headers 2]
			}

            # Extract BP-propability
            #
           if {[set probstart [lsearch $args "-prob"]] > -1} {
               set DStr($win_no,l_prob) [lindex $args [expr {$probstart+1}]]
           }

            # Extract Proj-Name
            #
            if {[set namestart [lsearch $args "-projname"]] > -1} {
               set DStr($win_no,name) [lindex $args [expr {$namestart+1}]]
           }
       }

        ### setup missing essentials
        #
		if {![info exists DStr($win_no,name)]} {
		   set DStr($win_no,name) "$win_no"
		   puts "No structure name provided to DrawStructure!"
		}
        if {![info exists DStr($win_no,l_prob)]} {
            for {set i 0} {$i < $seqlen} {incr i} {
				lappend DStr($win_no,l_prob) 1.0
			}
            puts "No base-pair probabilities provided to DrawStructure!"
        }

    }
	# ParseArgs




    ###   CreateFileMenu
    #
    # Guess what... creates the file menu
    #
    #
    proc CreateFileMenu {m optimal_struct} {
		##################################
        variable DStr

        menubutton $DStr($m,win).f_menu.file    \
                       -text File                    \
                       -underline 0                \
                       -menu $DStr($m,win).f_menu.file.menu

        set thismenu $DStr($m,win).f_menu.file.menu
		menu $thismenu

        #        $DStr($m,win).f_menu.file.menu add command    \
        #                    -label "Open SquFile "                \
        #                    -underline 0
#        if {$optimal_struct} {
#               $DStr($m,win).f_menu.file.menu add command    \
#                           -label "Write SQU file"                \
#                           -command "WriteSquFile 0 $DStr($m,name)"
#        } else {
#            $DStr($m,win).f_menu.file.menu add command    \
#                           -label "Write SQU file"                \
#                           -command "WriteSubSquFile $DStr($m,name) $m"
#        }
#        $DStr($m,win).f_menu.file.menu add sep
		$thismenu add cascade -menu $thismenu.save -label "Save..."
        $thismenu add cascade -menu $thismenu.print -label "Print..."

        $thismenu add sep
        $thismenu add command -command "destroy $DStr($m,win)" -label "Close "


		set thismenu $DStr($m,win).f_menu.file.menu.save
		menu $thismenu
		$thismenu add command -label "Consensus" \
			-command "[namespace code SaveCons] $m cs-consensus"
		$thismenu add cascade -label "Connect..." -menu $thismenu.ct

		menu  $thismenu.ct
		$thismenu.ct add command -label "Zuker" \
			-command "[namespace code SaveCons] $m connect zuker"
		$thismenu.ct add command -label "Gcg" \
			-command "[namespace code SaveCons] $m connect gcg"
		$thismenu.ct add command -label "Rnaviz" \
			-command "[namespace code SaveCons] $m connect rnaviz"

		$thismenu add command -label "Rnaml" \
			-command "[namespace code SaveCons] $m rnaml"
		$thismenu add command -label "Stockholm" \
			-command "[namespace code SaveCons] $m stockholm"
		
		

		set thismenu $DStr($m,win).f_menu.file.menu.print
        menu $thismenu

        $thismenu add command -label "to printer " \
                   -underline 3 -accelerator "Ctrl+p" \
                   -command "[namespace code PrintStructWin] print $m"
        bind $DStr($m,win) <Control-p> "[namespace code PrintStructWin] printer $m"
        $thismenu add command -label "to color printer "        \
                   -command "[namespace code PrintStructWin] printcolor $m"
        $thismenu add command -label "to file "            \
                   -command "[namespace code PrintStructWin] save $m"
        $thismenu add command -label "to screen "            \
                   -command "[namespace code PrintStructWin] screen $m"

    }
	# CreateFileMenu





    ###   CreateDisplayMenu
    #
    # Creates the display option
    # Values are stored in global array DS_DisplOpt
    #
    proc CreateDisplayMenu {m} {
		######################
          variable DStr
          variable DS_DisplOpt

        menubutton $DStr($m,win).f_menu.displ \
                   -text "Display" \
                   -menu $DStr($m,win).f_menu.displ.menu \
                   -underline 0

        menu $DStr($m,win).f_menu.displ.menu

        foreach displmode {"sequence" "numbers" "header"} {
			set label [format "%s%s " \
                               [string toupper [string index $displmode 0]] \
                               [string range $displmode 1 end]]
			$DStr($m,win).f_menu.displ.menu add check \
				-label $label \
				-variable [namespace current]::DS_DisplOpt($m,${displmode}_state) \
				-onvalue normal -offvalue hidden \
				-command "[namespace code ChangeItemConfig] $displmode $m"
           }

           $DStr($m,win).f_menu.displ.menu add sep

           $DStr($m,win).f_menu.displ.menu add cascade        \
                   -label "Background "            \
                   -menu $DStr($m,win).f_menu.displ.menu.bg

           $DStr($m,win).f_menu.displ.menu add cascade \
                   -label "Probability " \
                   -menu $DStr($m,win).f_menu.displ.menu.prob

           $DStr($m,win).f_menu.displ.menu add cascade \
                   -label "Single bp " \
                   -menu $DStr($m,win).f_menu.displ.menu.singles

           $DStr($m,win).f_menu.displ.menu add cascade \
                   -label "Loop mode " \
                   -menu $DStr($m,win).f_menu.displ.menu.loops

           $DStr($m,win).f_menu.displ.menu add cascade \
                   -label "Char mode " \
                   -menu $DStr($m,win).f_menu.displ.menu.chars

           $DStr($m,win).f_menu.displ.menu add command \
                   -label "Line widths " \
                   -command "[namespace code GetLineWidths] $m"

       menu    $DStr($m,win).f_menu.displ.menu.bg

           $DStr($m,win).f_menu.displ.menu.bg add radio \
                   -label "black" \
                   -variable [namespace current]::DS_DisplOpt($m,bg_col) \
                   -value black \
                   -command "[namespace code ChangeItemConfig] black $m"

           $DStr($m,win).f_menu.displ.menu.bg add radio \
                   -label "white" \
                   -variable [namespace current]::DS_DisplOpt($m,bg_col) \
                   -value white \
                   -command "[namespace code ChangeItemConfig] white $m"

       menu    $DStr($m,win).f_menu.displ.menu.prob

           $DStr($m,win).f_menu.displ.menu.prob add radio \
                   -label "colored" \
                   -variable [namespace current]::DS_DisplOpt($m,prob_col) \
                   -value colored \
                   -command "[namespace code CreateStructure] $m"

           $DStr($m,win).f_menu.displ.menu.prob add radio \
                   -label "greyscaled" \
                   -variable [namespace current]::DS_DisplOpt($m,prob_col) \
                   -value greyscaled \
                   -command "[namespace code CreateStructure] $m"

       menu    $DStr($m,win).f_menu.displ.menu.singles

           $DStr($m,win).f_menu.displ.menu.singles add radio \
                   -label "off" \
                   -variable [namespace current]::DS_DisplOpt($m,singles) \
                   -value 0 \
                   -command "[namespace code CreateStructure] $m"

           $DStr($m,win).f_menu.displ.menu.singles add radio \
                   -label "on" \
                   -variable [namespace current]::DS_DisplOpt($m,singles) \
                   -value 1 \
                   -command "[namespace code CreateStructure] $m"

       menu    $DStr($m,win).f_menu.displ.menu.loops

           $DStr($m,win).f_menu.displ.menu.loops add    radio \
                   -label "straight" \
                   -variable [namespace current]::DS_DisplOpt($m,loopmode) \
                   -value straight \
                   -command "[namespace code CreateStructure] $m"

           $DStr($m,win).f_menu.displ.menu.loops add    radio \
                   -label "symmetric" \
                   -variable [namespace current]::DS_DisplOpt($m,loopmode) \
                   -value symmetric \
                   -command "[namespace code CreateStructure] $m"

       menu    $DStr($m,win).f_menu.displ.menu.chars

           $DStr($m,win).f_menu.displ.menu.chars add    radio \
                   -label "uppercase" \
                   -variable [namespace current]::DS_DisplOpt($m,charmode) \
                   -value uppercase \
                   -command "[namespace code CreateStructure] $m"

           $DStr($m,win).f_menu.displ.menu.chars add    radio \
                   -label "lowercase" \
                   -variable [namespace current]::DS_DisplOpt($m,charmode) \
                   -value lowercase \
                   -command "[namespace code CreateStructure] $m"

    }
	# CreateDisplayMenu



    ###   CreateButtons
    #
    # Creates the Tool Buttons, with
    #        images, if zoom_img and reset_img are set
    #        text, otherwise
    #
    proc CreateButtons {no} {
		###################
          variable DStr
          variable zoom_img
          variable reset_img

        ### zoom
        #
       if {$zoom_img!=""} {
            button $DStr($no,win).f_tools.b_zoom \
                     -image $zoom_img
        } else {
            button $DStr($no,win).f_tools.b_zoom \
                     -text "Zoom"
        }

        ### reset
        #
       if {$reset_img!=""} {
           button $DStr($no,win).f_tools.b_reset \
                    -image $reset_img \
                    -command "[namespace code Reset] $no"
        } else {
           button $DStr($no,win).f_tools.b_reset \
                  -text "Reset" \
                     -command "[namespace code Reset] $no"
        }

    }
	# CreateButtons




    ###   CreateStructure
    #
    #
    #
    proc CreateStructure {m} {
		####################
         variable DStr
         variable debug_level
         variable   DS_DisplOpt

        set bond_width $DStr($m,bond_width);    # 4;    # 1
        set rect_width $DStr($m,rect_width);    # 2;    # 1
        set back_width $DStr($m,back_width);    # 2;    # 1

        $DStr($m,can) delete all

        if {$debug_level>4} {
            puts "DrawStruct::CreateStructure:  all destroyed in DStr($m,can)=$DStr($m,can)"
        }

        #####   canvas-scale computing
        #

        GetXYPlot $m

        scan [lindex $DStr($m,l_xyplot) 0] "%f %f %f %f %f %f %f %f %f %f %f %f" \
                               minx miny maxx maxy \
                               label1_x label1_y label2_x label2_y label3_x label3_y \
                               char_size char_width

        set xdraw_size [expr {$maxx-$minx}]    ;# dx
        set ydraw_size [expr {$maxy-$miny}]    ;# dy

        if {$debug_level>2} {
            puts "DrawStruct::CreateStructure startup"
            puts "\tmaxx=$maxx  minx=$minx"
            puts "\tmaxy=$maxy  miny=$miny"
            puts "\tdx=$xdraw_size  dy=$ydraw_size"
        }


        set max_draw_size  [expr {($xdraw_size>$ydraw_size)?$xdraw_size:$ydraw_size}]
        set fpixels        [winfo fpixels $DStr($m,win) $DStr($m,csize)c]
        set csize          [expr {$fpixels/$max_draw_size}]
        set DStr($m,scale) [expr {$fpixels/$max_draw_size}]

        #vputs [format "fontsize = %d pixel" [expr {int($char_size*$DStr($m,scale))}]] 5
        set pixel2point [winfo fpixels $DStr($m,win) 0.424c]
        if {$debug_level>3} {
            puts "DrawStruct::CreateStructure 12 point = $pixel2point pixel"
        }
        set pt_dstr [expr {int($char_size*$DStr($m,scale)*120./$pixel2point)}]
        if {$debug_level>2} {
            puts "DrawStruct::CreateStructure fontsize = $pt_dstr point"
        }
        set ft_dstr [DetermineFontSize $pt_dstr ]

        $DStr($m,can) config -width [expr {$csize*$xdraw_size}] \
                       -height [expr {$csize*$ydraw_size}] \
                       -scrollregion "[expr {$DStr($m,scale)*$minx}] [expr {$DStr($m,scale)*$miny}] \
                                           [expr {$DStr($m,scale)*$maxx}] [expr {$DStr($m,scale)*$maxy}]"


       $DStr($m,can) create rectangle $minx $miny $maxx $maxy \
                       -fill $DS_DisplOpt($m,bg_col) \
                       -tags ccolor

       ### Headlines Schreiben ###############################################
       if {$debug_level>2} {
            puts "DrawStruct::CreateStructure Writing headlines"
       }

       if {[info exists DStr($m,headline,1)]} {

           $DStr($m,can) create text $label1_x $label1_y \
                       -text $DStr($m,headline,1) \
                       -anchor n \
                       -fill $DS_DisplOpt($m,header_col) \
                       -font $ft_dstr \
                       -tags "label header" \
                       -state $DS_DisplOpt($m,header_state)

           $DStr($m,can) create text $label2_x $label2_y \
                       -text $DStr($m,headline,2) \
                       -anchor n \
                       -fill $DS_DisplOpt($m,header_col) \
                       -font $ft_dstr \
                       -tags "label header" \
                       -state $DS_DisplOpt($m,header_state)

           $DStr($m,can) create text $label3_x $label3_y \
                       -text $DStr($m,headline,3) \
                       -anchor n \
                       -fill $DS_DisplOpt($m,header_col) \
                       -font $ft_dstr \
                       -tags "label header" \
                       -state $DS_DisplOpt($m,header_state)
       }

       ################# Create Structure ###############################################

       for {set i 0} {$i<$DStr($m,seqlen)} {incr i} {

           scan [lindex $DStr($m,l_xyplot) [expr {1+$i}]]        "%d %d %s %d %f %f %f %f %f %f %f %f %f %f" \
                       nt($i) j($i) nuc($i) tag($i) nuc_x($i) nuc_y($i) \
                       back_x1($i) back_y1($i) back_x2($i) back_y2($i) \
                       bond_x1($i) bond_y1($i) bond_x2($i) bond_y2($i)

           scan [lindex $DStr($m,l_xyplot) [expr {1+$i}]]        "%*d %*d %*s %*d %*f %*f %*f %*f %*f %*f %*f \
                                                                           %*f %*f %*f %f %f %f %f %f %f %d" \
                       l_line_x1($i) l_line_y1($i) l_line_x2($i) l_line_y2($i) \
                       l_pos_x($i) l_pos_y($i) label($i)

       }

       if         {$DS_DisplOpt($m,charmode)=="uppercase"} {
           for {set i 0} {$i<$DStr($m,seqlen)} {incr i} {
               set nuc($i) [string toupper $nuc($i)]
           }
       } elseif    {$DS_DisplOpt($m,charmode)=="lowercase"} {
           for {set i 0} {$i<$DStr($m,seqlen)} {incr i} {
               set nuc($i) [string tolower $nuc($i)]
           }
       }

   #    set DStr($m,col1) [lindex $DStr($m,bg) 1]
   #    set DStr($m,col2) [lindex $DStr($m,bg) 2]

       for {set i 0} {$i<$DStr($m,seqlen)} {incr i} {

           if {$nuc($i) != "X"} {

               set prob [lindex $DStr($m,l_prob) $i]

                   if {$DS_DisplOpt($m,prob_col) == "colored"} {

                   #### white - yellow - red / bg: black ####

                       set DStr($m,col) [ColorProb $prob]

                   } elseif {$DS_DisplOpt($m,bg_col) == "black"} {

                       #### darkgrey -> white / bg: black ####

                       if {$prob > 0} {
                           set Gr [format "%02x" [expr {int((255-(1-[lindex $DStr($m,l_prob) $i])*220))}]]
                           set DStr($m,col) "#$Gr$Gr$Gr"
                       } else {
                           set DStr($m,col) "#FFFFFF"
                       }

                   } elseif {$DS_DisplOpt($m,bg_col) == "white"} {

                       #### lightgrey -> black / bg: white ####

                       if {$prob > 0} {
                           set Gr [format "%02x" [expr {int((1-[lindex $DStr($m,l_prob) $i])*200)}]]
                           set DStr($m,col) "#$Gr$Gr$Gr"
                       } else {
                           set DStr($m,col) "#FFFFFF"
                       }
                   }

                   if {($tag($i)<0)} {            ;# STAMMNUKLEOTIDE

                       if {($nt($i) == $DStr($m,seqlen))} {

                           $DStr($m,can) create line     $bond_x1($i) $bond_y1($i)    \
                                                       $bond_x2($i) $bond_y2($i)    \
                               -width $bond_width        \
                               -fill $DStr($m,col)        \
                               -tags "l_nt[expr {$i+1}] bo_nt$i bonds"

                           $DStr($m,can) create text     $nuc_x($i) $nuc_y($i)        \
                               -text $nuc($i)            \
                               -fill $DS_DisplOpt($m,sequence_col) \
                               -font $ft_dstr            \
                               -tags "ft$i sequence"        \
                               -state $DS_DisplOpt($m,sequence_state)

                           $DStr($m,can) create rectangle     $back_x1($i) $back_y1($i)    \
                                                               $back_x1($i) $back_y1($i)    \
                               -width $rect_width        \
                               -outline khaki2            \
                               -tags "stem[expr {-$tag($i)}] nt[expr {$i+1}]"

                       } else {

                           $DStr($m,can) create line    $back_x1($i) $back_y1($i)    \
                                                       $back_x2($i) $back_y2($i)    \
                               -width $back_width        \
                               -fill $DS_DisplOpt($m,backs_col)    \
                               -tags "l_nt[expr {$i+1}] ba_nt$i back[expr {$i+1}] backs"

                           $DStr($m,can) create line    $bond_x1($i) $bond_y1($i)    \
                                                       $bond_x2($i) $bond_y2($i)    \
                               -width $bond_width        \
                               -fill $DStr($m,col)        \
                               -tags "l_nt[expr {$i+1}] bo_nt$i bonds"

                           $DStr($m,can) create text    $nuc_x($i) $nuc_y($i)        \
                               -text $nuc($i)            \
                               -fill $DS_DisplOpt($m,sequence_col) \
                               -font $ft_dstr            \
                               -tags "ft$i sequence"        \
                               -state $DS_DisplOpt($m,sequence_state)

                           $DStr($m,can) create rectangle    $back_x1($i) $back_y1($i)    \
                                                               $back_x1($i) $back_y1($i)    \
                               -width $back_width        \
                               -outline khaki2            \
                               -tags "stem[expr {-$tag($i)}] nt[expr {$i+1}]"
                       }

                   } elseif {($tag($i)>0)} { ;# HELIX ohne Stammnukleotide

                       if {($nt($i) == $DStr($m,seqlen))} {

                       $DStr($m,can) create text     $nuc_x($i) $nuc_y($i)        \
                               -text $nuc($i)            \
                               -fill $DS_DisplOpt($m,sequence_col)    \
                               -font $ft_dstr            \
                               -tags "ft$i sequence"        \
                               -state $DS_DisplOpt($m,sequence_state)

                       $DStr($m,can) create rectangle    $back_x1($i) $back_y1($i)    \
                                                               $back_x1($i) $back_y1($i)    \
                               -width $back_width        \
                               -outline khaki2            \
                               -tags "helix$tag($i) nt[expr {$i+1}]"

                       $DStr($m,can) create line    $bond_x1($i) $bond_y1($i)    \
                                                       $bond_x2($i) $bond_y2($i)    \
                               -width $bond_width        \
                               -fill $DStr($m,col)        \
                               -tags "l_nt[expr {$i+1}] bo_nt$i bonds"
                       } else {

                           $DStr($m,can) create line    $back_x1($i) $back_y1($i)    \
                                                       $back_x2($i) $back_y2($i)    \
                               -width $back_width        \
                               -fill $DS_DisplOpt($m,backs_col)    \
                               -tags "l_nt[expr {$i+1}] ba_nt$i back[expr {$i+1}] backs"

                           $DStr($m,can) create line    $bond_x1($i) $bond_y1($i)    \
                                                       $bond_x2($i) $bond_y2($i)    \
                               -width $bond_width        \
                               -fill $DStr($m,col)        \
                               -tags "l_nt[expr {$i+1}] bo_nt$i bonds"

                           $DStr($m,can) create text    $nuc_x($i) $nuc_y($i)        \
                               -text $nuc($i)            \
                               -fill $DS_DisplOpt($m,sequence_col) \
                               -font $ft_dstr            \
                               -tags "ft$i sequence"        \
                               -state $DS_DisplOpt($m,sequence_state)

                           $DStr($m,can) create rectangle    $back_x1($i) $back_y1($i)    \
                                                               $back_x1($i) $back_y1($i)    \
                               -width $back_width        \
                               -outline khaki2            \
                               -tags "helix$tag($i) nt[expr {$i+1}]"
                       }

                   } elseif {($tag($i) == 0)} {;    # LOOPS (freie Nukleotide)

                       if {($nt($i) == $DStr($m,seqlen))} {

                       $DStr($m,can) create text    $nuc_x($i) $nuc_y($i)        \
                               -text $nuc($i)            \
                               -fill white            \
                               -font $ft_dstr            \
                               -tags "ft$i sequence"        \
                               -state $DS_DisplOpt($m,sequence_state)

                       $DStr($m,can) create rectangle    $back_x1($i) $back_y1($i)    \
                                                               $back_x1($i) $back_y1($i)    \
                               -width $back_width        \
                               -outline khaki2            \
                               -tags nt[expr {$i+1}]

                       } else {

                           $DStr($m,can) create line    $back_x1($i) $back_y1($i)    \
                                                       $back_x2($i) $back_y2($i)    \
                               -width $back_width        \
                               -fill $DS_DisplOpt($m,backs_col)    \
                               -tags "l_nt[expr {$i+1}] ba_nt$i back[expr {$i+1}] backs"

                           $DStr($m,can) create text     $nuc_x($i) $nuc_y($i)        \
                               -text $nuc($i)            \
                               -fill $DS_DisplOpt($m,sequence_col) \
                               -font $ft_dstr            \
                               -tags "ft$i sequence"        \
                               -state $DS_DisplOpt($m,sequence_state)

                           $DStr($m,can) create rectangle    $back_x1($i) $back_y1($i)    \
                                                               $back_x1($i) $back_y1($i)    \
                               -width $back_width        \
                               -outline khaki2            \
                               -tags nt[expr {$i+1}]
                       }
                   }

                   if {$label($i) > 0} {

                       $DStr($m,can) create line    $l_line_x1($i) $l_line_y1($i)    \
                                                   $l_line_x2($i) $l_line_y2($i)    \
                               -width 1            \
                               -fill $DS_DisplOpt($m,numbers_col)    \
                               -tags numbers \
                               -state $DS_DisplOpt($m,numbers_state)
                       $DStr($m,can) create text    $l_pos_x($i) $l_pos_y($i)    \
                               -text $label($i)        \
                               -fill $DS_DisplOpt($m,numbers_col)    \
                               -font $ft_dstr            \
                               -tags "numbers label"        \
                               -state $DS_DisplOpt($m,numbers_state)
                   }


               ############################################################################

                   $DStr($m,can) bind    stem$tag($i) <Any-Enter> "[namespace code StemEnter]"
                   $DStr($m,can) bind    stem$tag($i) <Any-Leave> "[namespace code StemLeave]"
                   $DStr($m,can) bind  helix$tag($i) <Any-Enter> "[namespace code HelixEnter]"
                   $DStr($m,can) bind  helix$tag($i) <Any-Leave> "[namespace code HelixLeave]"

               } else {

                   $DStr($m,can) itemconfig ba_nt$nt([expr {$i-2}]) -state hidden
               }

           }

           $DStr($m,can) bind selected <ButtonPress-1> "[namespace code ButtonPress] %x %y"
           bind $DStr($m,can) <B1-Motion> "[namespace code ButtonMove] %x %y"
           bind $DStr($m,can) <ButtonRelease-1> "[namespace code ButtonLeave] %x %y"

           bind $DStr($m,can) <3> "$DStr($m,can) scan mark %x %y"
           bind $DStr($m,can) <B3-Motion> "$DStr($m,can) scan dragto %x %y"

           bind $DStr($m,win).f_tools.b_zoom <Button-3> \
                 "[namespace code ZoomCanvas] $m 0.5 $minx $miny $maxx $maxy"
           bind $DStr($m,win).f_tools.b_zoom <Button-1> \
                 "[namespace code ZoomCanvas] $m 2 $minx $miny $maxx $maxy"
   #        bind $DStr($m,can) <ButtonPress-2> "DumpMousePos %x %y"

       $DStr($m,can) lower bonds backs

       foreach id [$DStr($m,can) find all] {
           $DStr($m,can) scale $id 0 0 $DStr($m,scale) $DStr($m,scale)
       }

        if {$debug_level>3} {
            puts "DrawStruct::CreateStructure done"
        }

   }
	# CreateStructure





    ############################################################################
    ###            #############################################################
    ###   EVENTS   #############################################################
    ###            #############################################################
    ############################################################################



   ###   StemEnter
   #
   #
   #
   proc StemEnter {} {        ;# STAMMNUKLEOTIDE
	   #############
         variable DStr
        variable   DS_DisplOpt

       set m [GetFocusedWindowNo]

       #vputs "StemEnter $m" 4

       set stem_tag1 [lindex [$DStr($m,can) gettags current] 0]
       #vputs $stem_tag1 6
       set stem_tag2 [$DStr($m,can) find withtag $stem_tag1]
       #vputs $stem_tag2 6
       set lack ""
       foreach element $stem_tag2 {
           lappend lack [$DStr($m,can) gettags $element]
       }
       set ntnr1 [lindex [lindex $lack 0] 1]
       set ntnr2 [lindex [lindex $lack 1] 1]
       #vputs "$ntnr1 $ntnr2" 6
       scan $ntnr1 "nt%d" ntnr1
       scan $ntnr2 "nt%d" ntnr2
       #set ntnr1 [expr $ntnr1] ;# wat soll dat?
       #set ntnr2 [expr $ntnr2]
       ##vputs "$ntnr1 $ntnr2" 6

       for {set i $ntnr1} {$i<=$ntnr2} {incr i} {

           if {$i==$ntnr2} {
               set lack [$DStr($m,can) find withtag nt$i]
               # #vputs ">$nt($i)" 6
               # #vputs ">$lack" 6
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }

           } else {

               set lack [$DStr($m,can) find withtag nt$i]
               # #vputs ">>$nt($i)" 6
               # #vputs ">>$lack" 6
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }
               set linien [$DStr($m,can) find withtag back$i]
               # #vputs "$linien"    6
               foreach element $linien {
                   $DStr($m,can) itemconfigure $element -fill $DS_DisplOpt($m,hlback_col)
                   $DStr($m,can) addtag selected withtag $element
               }
           }
       }
   }
	# StemEnter


	
   ###   StemLeave
   #
   #
   #
   proc StemLeave {} {
	   #############
         variable DStr
        variable   DS_DisplOpt

       set m [GetFocusedWindowNo]

       #vputs "StemLeave $m" 4

       set stem_tag1 [lindex [$DStr($m,can) gettags current] 0]
       #vputs $stem_tag1 6
       set stem_tag2 [$DStr($m,can) find withtag $stem_tag1]
       #vputs $stem_tag2 6
       set lack ""
       foreach element $stem_tag2 {
           lappend lack [$DStr($m,can) gettags $element]
       }
       set ntnr1 [lindex [lindex $lack 0] 1]
       set ntnr2 [lindex [lindex $lack 1] 1]
       #vputs "$ntnr1 $ntnr2" 6
       scan $ntnr1 "nt%d" ntnr1
       scan $ntnr2 "nt%d" ntnr2
       #set ntnr1 [expr $ntnr1] ;# wat soll dat?
       #set ntnr2 [expr $ntnr2]
       ##vputs "$ntnr1 $ntnr2" 6

       for {set i $ntnr1} {$i<=$ntnr2} {incr i} {

           if {$i==$ntnr2} {
               set lack [$DStr($m,can) find withtag nt$i]
               #vputs ">$lack" 6
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }

           } else {

               set lack [$DStr($m,can) find withtag nt$i]
               #vputs ">>$lack" 6
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }
               set linien [$DStr($m,can) find withtag back$i]
               #vputs "$linien"     6
               foreach element $linien {
                   $DStr($m,can) itemconfigure $element -fill $DS_DisplOpt($m,backs_col)
                   $DStr($m,can) addtag selected withtag $element
               }
           }
       }
   }
	# StemLeave


   ###   HelixEnter
   #
   #
   #
   proc HelixEnter {} {    ;# HELIX ohne Stammnukleotide
	   ##############
         variable DStr
          variable   DS_DisplOpt

       set m [GetFocusedWindowNo]

       #vputs "HelixEnter $m" 4

       set helix_tag1 [lindex [$DStr($m,can) gettags current] 0]
       #vputs $helix_tag1 6
       scan $helix_tag1 "helix%d" lack
       #vputs $lack 6

       set helix_tag2 [$DStr($m,can) find withtag stem$lack]
       set lack ""
       foreach element $helix_tag2 {
           lappend lack [$DStr($m,can) gettags $element]
       }
       set ntnr1 [lindex [lindex $lack 0] 1]
       set ntnr2 [lindex [lindex $lack 1] 1]
       scan $ntnr1 "nt%d" ntnr1
       scan $ntnr2 "nt%d" ntnr2

       for {set i $ntnr1} {$i<=$ntnr2} {incr i} {

           if {$i==$ntnr2} {
               set lack [$DStr($m,can) find withtag nt$i]
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }

           } else {

               set lack [$DStr($m,can) find withtag nt$i]
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }
               set linien [$DStr($m,can) find withtag back$i]
               foreach element $linien {
                   $DStr($m,can) itemconfigure $element -fill $DS_DisplOpt($m,hlback_col)
                   $DStr($m,can) addtag selected withtag $element
               }
           }
       }
   }
	# HelixEnter


	
   ###   HelixLeave
   #
   #
   #
   #
   proc HelixLeave {} {
	   ##############
         variable DStr
          variable   DS_DisplOpt

       set m [GetFocusedWindowNo]

       #vputs "HelixLeave $m" 4

       set helix_tag1 [lindex [$DStr($m,can) gettags current] 0]
       scan $helix_tag1 "helix%d" lack

       set helix_tag2 [$DStr($m,can) find withtag stem$lack]
       set lack ""
       foreach element $helix_tag2 {
           lappend lack [$DStr($m,can) gettags $element]
       }
       set ntnr1 [lindex [lindex $lack 0] 1]
       set ntnr2 [lindex [lindex $lack 1] 1]
       scan $ntnr1 "nt%d" ntnr1
       scan $ntnr2 "nt%d" ntnr2
       #set ntnr1 [expr $ntnr1] ;# wat soll dat?
       #set ntnr2 [expr $ntnr2]

       for {set i $ntnr1} {$i<=$ntnr2} {incr i} {

           if {$i==$ntnr2} {
               set lack [$DStr($m,can) find withtag nt$i]
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }

           } else {

               set lack [$DStr($m,can) find withtag nt$i]
               foreach element $lack {
                   $DStr($m,can) addtag selected withtag $element
               }
               set linien [$DStr($m,can) find withtag back$i]
               foreach element $linien {
                   $DStr($m,can) itemconfigure $element -fill $DS_DisplOpt($m,backs_col)
                   $DStr($m,can) addtag selected withtag $element
               }
           }
       }
   }
	# HelixLeave


   ###   ButtonPress
   #
    #
    #
   proc ButtonPress {x y} {
	   ##################
         variable DStr
        variable   DS_DisplOpt

       set m [GetFocusedWindowNo]

       #vputs "ButtonPress $m" 4

       #vputs "\n\n[lrange [$DStr($m,can) gettags current] 0 1]\n" 6
       $DStr($m,can) delete LINE
       if {[scan [lrange [$DStr($m,can) gettags current] 0 1] "helix%d nt%d" DStr($m,ntA) DStr($m,ntB)] == 2} {

           #vputs "$DStr($m,ntA) $DStr($m,ntB) [lrange [$DStr($m,can) gettags current] 0 1]" 6
           set DStr($m,PktAx) 0
           set DStr($m,PktAy) 0
           set DStr($m,PktBx) 0
           set DStr($m,PktBy) 0

           set DStr($m,stem_bp) [$DStr($m,can) find withtag stem$DStr($m,ntA)]

           set stem_a [lindex $DStr($m,stem_bp) 0]
           scan [lindex [$DStr($m,can) gettags $stem_a ] 1] "nt%d" stem_a
           set d_stem_a [expr {$DStr($m,ntB)-$stem_a}]

           set stem_b [lindex $DStr($m,stem_bp) 1]
           scan [lindex [$DStr($m,can) gettags $stem_b ] 1] "nt%d" stem_b
           set d_stem_b [expr {-($DStr($m,ntB)-$stem_b)}]

           set DStr($m,ntA) [expr {($d_stem_a<$d_stem_b)?$stem_a:$stem_b}]
           set coord_PktA [$DStr($m,can) coords nt$DStr($m,ntA)]
           #vputs "\nPunkt A: $coord_PktA" 6

           set coord_PktB [$DStr($m,can) coords current]
           #vputs "Punkt B: $coord_PktB" 6

           scan [lrange $coord_PktA 0 1] "%f %f" DStr($m,PktAx) DStr($m,PktAy)
           scan [lrange $coord_PktB 0 1] "%f %f" DStr($m,PktBx) DStr($m,PktBy)


           #vputs "line_width = $DStr($m,line_width)" 6
           $DStr($m,can) create line $DStr($m,PktAx) $DStr($m,PktAy)    \
                           [$DStr($m,can) canvasx $x] [$DStr($m,can) canvasy $y] \
                   -width $DStr($m,line_width)        \
                   -fill    $DS_DisplOpt($m,hlback_col)    \
                   -tags    LINE
       } else {
           puts "ARGH! DrawStruct::ButtonPress No Stem-Basepair selected"

           set DStr($m,PktAx) nix
           set DStr($m,PktAy) nix
           set DStr($m,PktBx) nix
           set DStr($m,PktBy) nix
           set DStr($m,ntA) nix
           set DStr($m,ntB) nix
           return break
       }
   }
	# ButtonPress



    ###   ButtonMove
    #
    #
    #
   proc ButtonMove {x y} {
	   #################
         variable DStr
          variable DS_DisplOpt

       set m [GetFocusedWindowNo]

       #vputs "ButtonMove $m" 4

       if {[info exists DStr($m,PktAy)] != 1 } {

           return break

       } else {

           if {($DStr($m,PktAx)=="nix")||($DStr($m,PktAy)=="nix")||($DStr($m,PktBx)=="nix")|| \
               ($DStr($m,PktBy)=="nix")||($DStr($m,ntA)=="nix")||($DStr($m,ntB)=="nix")} {

               return break

           } else {
                       $DStr($m,can) delete LINE
                       $DStr($m,can) create line $DStr($m,PktAx) $DStr($m,PktAy) \
                                                       [$DStr($m,can) canvasx $x] [$DStr($m,can) canvasy $y] \
                       -width $DStr($m,line_width)    \
                       -fill $DS_DisplOpt($m,hlback_col)    \
                       -tags LINE

   #            $DStr($m,can) create line $DStr($m,PktAx) $DStr($m,PktAy) $x $DStr($m,PktAy) \ ;###############
   #                    -width .3 \                                ;# Hilfslinien #
   #                -fill white \                                ;###############
   #                -tags LINE
   #            $DStr($m,can) create line ${DStr($m,PktAx)} ${DStr($m,PktAy)} \
   #                                                    ${DStr($m,PktBx)} ${DStr($m,PktAy)} \
   #                    -width .3 \
   #                -fill white \
   #                -tags LINE
   #            $DStr($m,can) create line ${DStr($m,PktBx)} ${DStr($m,PktBy)} \
   #                                                    ${DStr($m,PktBx)} ${DStr($m,PktAy)} \
   #                    -width .3 \
   #                -fill green \
   #                -tags LINE
   #            $DStr($m,can) create line $x ${DStr($m,PktAy)} $x $y \
   #                    -width .3 \
   #                -fill skyblue \
   #                -tags LINE
           }
       }
   }
	# ButtonMove



    ###   ButtonLeave
    #
    #
    #
    proc ButtonLeave {x y} {
		##################
         variable DStr

       set m [GetFocusedWindowNo]

       #vputs "ButtonLeave $m" 4

       if {[info exists DStr($m,PktAy)] != 1 } {

           return break

       } else {

           if {($DStr($m,PktAx)=="nix")||($DStr($m,PktAy)=="nix")||($DStr($m,PktBx)=="nix")||($DStr($m,PktBy)=="nix")||($DStr($m,ntA)=="nix")||($DStr($m,ntB)=="nix")} {

               return break

           } else {

               #vputs "\n    x / y\nA: $DStr($m,PktAx) / $DStr($m,PktAy) \nB: $DStr($m,PktBx) / $DStr($m,PktBy) \nC: [$DStr($m,can) canvasx $x] / [$DStr($m,can) canvasy $y]" 6

               set DStr($m,PktCx) 0
               set DStr($m,PktCy) 0
               set DStr($m,PktCx) [$DStr($m,can) canvasx $x]
               set DStr($m,PktCy) [$DStr($m,can) canvasy $y]

               set pi 3.141592654
               set Seite_a [expr {sqrt(pow($DStr($m,PktCx)-$DStr($m,PktBx),2)+pow($DStr($m,PktBy)-$DStr($m,PktCy),2))}]  ;# Dreieckseiten    #
               set Seite_b [expr {sqrt(pow($DStr($m,PktAx)-$DStr($m,PktCx),2)+pow($DStr($m,PktCy)-$DStr($m,PktAy),2))}]  ;# berechnen         #
               set Seite_c [expr {sqrt(pow($DStr($m,PktBx)-$DStr($m,PktAx),2)+pow($DStr($m,PktAy)-$DStr($m,PktBy),2))}]  ;#(nach Standard-Dreieck)#

                # FIXME prevent divide by zero
                if {$DStr($m,PktAx)==$DStr($m,PktBx)} {
                    return
                }


               set baseline_b [expr {$DStr($m,PktBx)-$DStr($m,PktAx)}]   ;###############
               set lot_b [expr {$DStr($m,PktAy)-$DStr($m,PktBy)}]     ;# Hilfslinien #
               set baseline_c [expr {$DStr($m,PktCx)-$DStr($m,PktAx)}]   ;# berechnen   #
               set lot_c [expr {$DStr($m,PktAy)-$DStr($m,PktCy)}]     ;###############

               ###### Winkel zwischen x-Achse (->) und Helix ####################
               set cos_phi [expr {(pow($Seite_c,2)+pow($baseline_b,2)-pow($lot_b,2))/(2*$Seite_c*$baseline_b)}]
               set phi [expr {acos($cos_phi)*180/$pi}]
               #vputs "Phi= $phi" 6

               ###### Winkel zwischen x-Achse (->) und Gerade ###################
               set cos_epsilon [expr {(pow($Seite_b,2)+pow($baseline_c,2)-pow($lot_c,2))/(2*$Seite_b*$baseline_c)}]
               set epsilon [expr {acos($cos_epsilon)*180/$pi}]
               #vputs "Epsilon= $epsilon" 6

               ###### Drehwinkel (zur Ueberpruefung von Delta, s.u.) ############
               set cos_alpha [expr {(pow($Seite_b,2)+pow($Seite_c,2)-pow($Seite_a,2))/(2*$Seite_b*$Seite_c)}]
               set alpha [expr {acos($cos_alpha)*180/$pi}]
               #vputs "Alpha= $alpha\n" 6

               ###### Fallunterscheidungen ######################################

               if {$baseline_b>=0    &&    $lot_b>=0} {        ;### Helix im 1.Quadrant

                   if {$baseline_c>0 && $lot_c>0} {
                       set delta [expr {$phi-$epsilon}]
                   } elseif {$baseline_c>0 && $lot_c<0} {
                       set delta [expr {$phi+$epsilon}]
                   } elseif {$baseline_c<0 && $lot_c>0} {
                       set delta [expr {$phi-$epsilon}]
                   } elseif {$baseline_c<0 && $lot_c<0} {
                       set delta [expr {$phi+$epsilon}]
                   }

               } elseif {$baseline_b<0 && $lot_b>=0} {;    ### Helix im 2.Quadrant

                   if {$baseline_c>=0 && $lot_c>=0} {
                       set delta [expr {$phi-$epsilon}]
                   } elseif {$baseline_c>0 && $lot_c<0} {
                       set delta [expr {$phi+$epsilon}]
                   } elseif {$baseline_c<0 && $lot_c>0} {
                       set delta [expr {$phi-$epsilon}]
                   } elseif {$baseline_c<0 && $lot_c<0} {
                       set delta [expr {$phi+$epsilon}]
                   }

               } elseif {$baseline_b<=0 && $lot_b<0} {    ;### Helix im 3.Quadrant

                   if {$baseline_c>0 && $lot_c>0} {
                       set delta [expr {-($phi+$epsilon)}]
                   } elseif {$baseline_c>0 && $lot_c<0} {
                       set delta [expr {$epsilon-$phi}]
                   } elseif {$baseline_c<0 && $lot_c>0} {
                       set delta [expr {-($phi+$epsilon)}]
                   } elseif {$baseline_c<0 && $lot_c<0} {
                       set delta [expr {$epsilon-$phi}]
                   }

               } elseif {$baseline_b>0 && $lot_b<0} {    ;### Helix im 4.Quadrant

                   if {$baseline_c>0 && $lot_c>0} {
                       set delta [expr {-($phi+$epsilon)}]
                   } elseif {$baseline_c>0 && $lot_c<0} {
                       set delta [expr {$epsilon-$phi}]
                   } elseif {$baseline_c<0 && $lot_c>0} {
                       set delta [expr {-($phi+$epsilon)}]
                   } elseif {$baseline_c<0 && $lot_c<0} {
                       set delta [expr {$epsilon-$phi}]
                   }
               }

               if {$delta>180} {
                   set delta [expr {$delta-360}]

               } elseif {$delta<-180} {
                   set delta [expr {$delta+360}]
               }

               #vputs "Delta= $delta" 6

               set DStr($m,PktAx) nix
               set DStr($m,PktAy) nix
               set DStr($m,PktBx) nix
               set DStr($m,PktBy) nix
               set DStr($m,PktCx) nix
               set DStr($m,PktCy) nix

               scan [lindex [$DStr($m,can) gettags [lindex $DStr($m,stem_bp) 0]] 1] "nt%d" stempair1
               scan [lindex [$DStr($m,can) gettags [lindex $DStr($m,stem_bp) 1]] 1] "nt%d" stempair2

               #vputs "stempair1=$stempair1, stempair2=$stempair2" 6

               if {![info exists DStr($m,pivot)]} {
                   set DStr($m,pivot) ""
               }

               set pivotstart [lsearch $DStr($m,pivot) $stempair1]

               #vputs "pivotstart=$pivotstart" 6

               if {$pivotstart > -1} {
                   set deltastart [expr {$pivotstart+2}]
                   set olddelta [lindex $DStr($m,pivot) $deltastart]
                   set newdelta [expr {$olddelta+$delta}]
                   set DStr($m,pivot) [lreplace $DStr($m,pivot) $deltastart $deltastart $newdelta]
               } else {
                   lappend DStr($m,pivot) $stempair1 $stempair2 $delta
               }

               #vputs "DStr($m,pivot)=$DStr($m,pivot) " 6

               #vputs "$DS_DisplOpt($m,sequence_state) $DS_DisplOpt($m,numbers_state) $DS_DisplOpt($m,header_state)" 6

               CreateStructure $m
           }
       }
   }
	# ButtonLeave

    ############################################################################
    ###                #########################################################
    ###   END EVENTS   #########################################################
    ###                #########################################################
    ############################################################################





    ############################################################################
    ###          ###############################################################
    ###   MISC   ###############################################################
    ###          ###############################################################
    ############################################################################






    ###   WinSetup 
    #
    # Setup called once for a new drawstruct window
    #
    proc WinSetup {win_no} {
		#################

        variable DStr
        variable DS_DisplOpt

        if {![info exists DStr($win_no,csize)]} {
            set DStr($win_no,csize) 20
        }


         set DS_DisplOpt($win_no,bg_col)         black
          set DS_DisplOpt($win_no,header_col)     white
          set DS_DisplOpt($win_no,backs_col)      white
          set DS_DisplOpt($win_no,sequence_col)   white
          set DS_DisplOpt($win_no,numbers_col)    white
          set DS_DisplOpt($win_no,hlback_col)     red
          set DS_DisplOpt($win_no,sequence_state) normal
          set DS_DisplOpt($win_no,numbers_state)     normal
          set DS_DisplOpt($win_no,header_state)     normal
          set DS_DisplOpt($win_no,singles)             0
          set DS_DisplOpt($win_no,loopmode)         straight
          set DS_DisplOpt($win_no,prob_col)         colored
          set DS_DisplOpt($win_no,charmode)         uppercase

          set DStr($win_no,bond_width)                2;    # 4    # "hydrogen" bonds
          set DStr($win_no,rect_width)                1;    # 2    # mark for stem nts (for selection)
          set DStr($win_no,back_width)                0;    # 2    # backbone
          set DStr($win_no,line_width)                1;    # 3    # width of border and angle



        set DStr($win_no,win) .w_drawstruct_$win_no
       set DStr($win_no,can) $DStr($win_no,win).c_dstr

       toplevel    $DStr($win_no,win)
       wm title       $DStr($win_no,win) "Draw Structure $win_no"
       wm iconname    $DStr($win_no,win) "Draw Structure $win_no"
       wm protocol $DStr($win_no,win) WM_DELETE_WINDOW "[namespace code B4DelWindow] $DStr($win_no,win) $win_no"

       set DStr($win_no,vir_csize)        $DStr($win_no,csize)

    }
	# WinSetup




   ###   GetLineWidths
   #
   #
   proc GetLineWidths {m} {
	   ##################
         variable DStr

       set old_bond_width    $DStr($m,bond_width)
       set old_back_width    $DStr($m,back_width)

       set f .lineWidths
       if {[winfo exists $f]==0} {
           set oldFocus [focus]
           set row 0
           toplevel        $f
           wm title        $f LineWidths
           wm iconname $f LineWidths
           label            $f.message    -wraplength 5i  \
                                           -justify left    \
                                           -text "Set line widths\n\n1 =< widths =< 5"
           grid            $f.message    -row $row; incr row
           label           $f.leer$row    -text " "
           grid            $f.leer$row    -row $row; incr row
           frame           $f.a            -bd 2
           label           $f.a.label    -text "bond width = "
           entry           $f.a.entry    -relief sunken -width 2
           grid            $f.a.label    -column 0 -row 0 -sticky w
           grid            $f.a.entry    -column 1 -row 0 -sticky e
           grid            $f.a            -row $row; incr row
           bind            $f.a.entry <Return> "[namespace code GetWidthsVal] $f $m"
           frame           $f.b            -bd 2
           label           $f.b.label    -text "backbone width = "
           entry           $f.b.entry    -relief sunken -width 2
           grid            $f.b.label    -column 0 -row 0 -sticky w
           grid            $f.b.entry    -column 1 -row 0 -sticky e
           grid            $f.b            -row $row; incr row
           bind            $f.b.entry <Return> "[namespace code GetWidthsVal] $f $m"
           $f.a.entry insert 0 [format "%-2d" $DStr($m,bond_width)]
           $f.b.entry insert 0 [format "%-2d" $DStr($m,back_width)]
           label           $f.leer$row    -text " "
           grid            $f.leer$row    -row $row; incr row
           frame           $f.buttons
           grid            $f.buttons    -row $row; incr row
           button          $f.buttons.dismiss    -text Dismiss \
                        -command "[namespace code CancelEntry] $f $oldFocus"
           button        $f.buttons.accept        -text Return \
                        -command "[namespace code ModWidths] $old_bond_width $old_back_width $f $m $oldFocus"
           grid            $f.buttons.dismiss    -column 0 -row 0
           grid            $f.buttons.accept        -column 1 -row 0
       }
       set pos [winfo pointerxy $DStr($m,win)]
       wm geometry $f [format "+%d+%d" [lindex $pos 0] [lindex $pos 1]]
       raise $f
       grab  $f
       focus $f

   }
	# GetLineWidths



   ###   ModWidths
   #
   #
   #
   proc ModWidths {old_bond_width old_back_width f m oldFocus} {
	   ######################################################
         variable DStr

       GetWidthsVal $f $m
       CancelEntry $f $oldFocus
       if {$old_bond_width!=$DStr($m,bond_width)    ||\
            $old_back_width!=$DStr($m,back_width)} {
            CreateStructure $m
       }
   }
	# ModWidths



   ###   GetWidthsVal
   #
   #
   #
   proc GetWidthsVal {entryW m} {
	   ########################
         variable DStr

      foreach entry [list a b] {
           set x [${entryW}.${entry}.entry get]
          set x [expr {($x>5) ? 5 : $x}]
          set x [expr {($x<1) ? 1 : $x}]
          ${entryW}.${entry}.entry delete 0 end
          ${entryW}.${entry}.entry insert 0 [format "%-2d" $x]
          if {$entry=="a"} {
               set DStr($m,bond_width) $x
           }
          if {$entry=="b"} {
               set DStr($m,back_width) $x
           }
       }
   }
	# GetWidthsVal





   ###   GetFocusedWindowNo
   #
   #
   #
   proc GetFocusedWindowNo {} {
	   #####################
       variable debug_level

       regexp {.*([0-9]+)$} [focus] dummy dswnr
        if {![info exists dswnr]} {
            puts "ARRGH! DrawStruct::GetFocusedWindowNo couldn't screen focused windows"
            puts "Focus DrawStruct window to fix me!"
            set dswnr -1
        }

       return $dswnr
   }
	# GetFocusedWindowNo






   ############################################################################################
   ### PROZEDUREN #############################################################################
   ############################################################################################
   ################################################################################


   ###   GetXYPlot
   #
   #
   #
   proc GetXYPlot {m {backbone 0}} {
	   ###########################
        variable DStr
        variable DS_DisplOpt

        # this is shit, but otherwise
        # I wasn't able to get these vars from core
        # AW
        variable ds_core_nucrange
        variable ds_core_single
        variable ds_core_loopmode
        variable ds_core_pivot
        variable ds_core_text

        # accessed from core: don't change var-names
        set ds_core_nucrange [list 0 0]
        set ds_core_single   $DS_DisplOpt($m,singles)
        set ds_core_loopmode $DS_DisplOpt($m,loopmode)
        if {[array names DStr "$m,pivot"]=="$m,pivot"} {
           set ds_core_pivot $DStr($m,pivot)
           # set ds_core_pivot  -pivot [list 3 255 -53.9937439909]
        } else {
            catch {unset ds_core_pivot}
        }
        if {[info exists DStr($m,headline,1)]} {
            set ds_core_text 1
        } else {
            set ds_core_text 0
        }

       set DStr($m,l_xyplot) [cs_drawstructure $DStr($m,seq) \
                                               $DStr($m,l_basepairs) \
                                               $DS_DisplOpt($m,loopmode)]
   }
   # GetXYPlot


    ###   PrintStructWin
    #
    #
    # mode:
    #   printer
    #   printercolor
    #   file
    #   screen
    #
    #
    proc PrintStructWin {mode m} {
		########################
        variable DStr
        variable tmp_dir
        variable print_pswidth


        if {[info exists DStr($m,name)]} {
            # FIXME remove: since PrintPS uses it's own workdir
            # set psfilename "${tmp_dir}/$DStr($m,name)_os.ps"
            set psfilename "$DStr($m,name)_os.ps"
        } else {
            # FIXME remove: since PrintPS uses it's own workdir
            #set psfilename "${tmp_dir}/my_os.ps"
            set psfilename "dummy_os.ps"
        }



        if {[$DStr($m,can) cget -height]>[$DStr($m,can) cget -width]} {; # immer die laengere Seite des Canvas
            set args  "-pageheight ${print_pswidth}c"
        } else {
            #set args "-pagewidth [expr {sqrt(2.)*${print_pswidth}}]c -rotate true"
            set args "-rotate true"
        }

        PrintPS  $DStr($m,can) $psfilename $mode $args

    }
    # end PrintStructWin


	
    ###   ZoomCanvas
    #
    #
    #
    proc ZoomCanvas {m scale_fac minx miny maxx maxy} {
		#############################################
        variable DStr
        variable debug_level

       if {$DStr($m,vir_csize) == $DStr($m,csize) && $scale_fac<1} {
           return break
       } else {
           set DStr($m,vir_csize) [expr {$DStr($m,vir_csize)*$scale_fac}]
           set dstr_scale             [expr {$DStr($m,scale)*$scale_fac}]
           set DStr($m,scale)        $dstr_scale

           foreach id [$DStr($m,can) find all] {
               $DStr($m,can) scale $id 0 0 $scale_fac $scale_fac
           }

           $DStr($m,can) config -scrollregion "[expr {$dstr_scale*$minx}] [expr {$dstr_scale*$miny}] \
                                                [expr {$dstr_scale*$maxx}] [expr {$dstr_scale*$maxy}]"

           if {$debug_level>4} {
                puts "DrawStruct::ZoomCanvas: DStr($m,vir_csize) = $DStr($m,vir_csize)"
            }

           set dstr_pt_old [lindex [split [lindex [$DStr($m,can) itemconfig ft1 -font] 4] -] 8]
           set dstr_pt     [expr {int($dstr_pt_old*$scale_fac)}]
           set dstr_ft     [DetermineFontSize $dstr_pt ]
           foreach id [$DStr($m,can) find withtag label] {
               $DStr($m,can) itemconfig $id -font $dstr_ft
           }
           for {set i 0} {$i<$DStr($m,seqlen)} {incr i} {
               $DStr($m,can) itemconfig ft$i -font $dstr_ft
           }
       }
   }
	# ZoomCanvas




    ###   ChangeItemConfig   ###################################################
   #
    #
    #
    proc ChangeItemConfig {tag  m} {
		##########################
         variable DStr
         variable debug_level
         variable DS_DisplOpt

        if {$debug_level>4} {
            puts "DrawStruct::ChangeItemConfig:   tag = $tag"
            puts "DrawStruct::ChangeItemConfig:   m = $m"
        }

       if {$tag == "white"} {
           $DStr($m,can) itemconfigure ccolor -fill $DS_DisplOpt($m,bg_col)
           if {$DS_DisplOpt($m,prob_col) != "colored"} {
               CreateStructure $m
           }
           foreach id {"header" "numbers" "sequence" "backs"} {
               $DStr($m,can) itemconfigure $id -fill black
               set DS_DisplOpt($m,[format "%s_col" $id]) black
           }

       } elseif {$tag == "black"} {
           $DStr($m,can) itemconfigure ccolor -fill $DS_DisplOpt($m,bg_col)
           if {$DS_DisplOpt($m,prob_col) != "colored"} {
               CreateStructure $m
           }
           foreach id {"header" "numbers" "sequence" "backs"} {
               $DStr($m,can) itemconfigure $id -fill white
               set DS_DisplOpt($m,[format "%s_col" $id]) white
           }

       } else {
           $DStr($m,can) itemconfigure $tag -state $DS_DisplOpt($m,${tag}_state)
       }

   }
	# ChangeItemConfig




   ###   Reset
   #
    # Discards all made changes to structure and recreates it
    #
   proc Reset {m} {
	   ##########
         variable DStr

       set DStr($m,pivot) ""
       if {[info exists DStr($m,pivot)]} {
           unset DStr($m,pivot)
           CreateStructure $m
       }
   }
	# Reset




    ### CancelEntry
    #
    #
    proc CancelEntry {f oldFocus} {
		#########################

        focus        $oldFocus
        grab release $f
        destroy      $f

    }
	# CancelEntry




    ###   DumpMousePos
    #
    # report mouse position in actual canvas to stdout
    #
   proc DumpMousePos {x y} {
	   ###################
         variable DStr

        set m [GetFocusedWindowNo]
        set cx [$DStr($m,can) canvasx $x]
        set cy [$DStr($m,can) canvasy $y]
        puts "mouse position in canvas $m: x=$cx y=$cy"

   }
	# DumpMousePos




    ###   B4DelWindow
    #
    # Unsets all variables associated with 2 be destroyed window
    #
    proc B4DelWindow {win_path win_no} {
		##############################
          variable DStr
          variable debug_level

        set oldsize [array size DStr]
        array unset DStr "$win_no,*"
        set destr_elems [expr {$oldsize-[array size DStr]}]
        if {$debug_level} {
            puts "DrawStruct::B4DelWindow:   $destr_elems associated vars unset"
        }

        destroy $win_path
        if {$debug_level} {puts "DrawStruct::B4DelWindow:   $win_path destroyed"}

    }
	# B4DelWindow




    ###   DetermineFontSize
    #
    #
    #
    proc DetermineFontSize {fontsize} {
		############################
        # global : export 2 sd_procs

        set fontsize [expr {int($fontsize)}]
        set server [winfo server .]

        #vputs "DetermineFontSize:" 5
        #vputs "fontsize = >$fontsize<" 5
        #vputs "server    = >$server<" 5

        if {[lsearch -exact $server DECWINDOWS] > -1 } {;    # VMS interpoliert nicht
            #vputs "DECWINDOWS" 5
               if {[lsearch -exact $server X11R0] > -1 } {
               #vputs "X11R0" 5
               if {$fontsize < 90} {
                      set fontsize    80
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
            if {$fontsize < 10} {set fontsize 10};        # kleiner mag X11/IRIX nicht
        }
        #vputs "fontsize: $fontsize" 5
        #    if {$fontsize < 50} {set fontsize 50};    # kleiner mag X11/IRIX nicht

        return [format "-adobe-courier-bold-r-*-*-*-%d-*-*-*-*-*-*" $fontsize]

    }
	# DetermineFontSize



	###   SaveCons
	#
	# save consensus
	#
	proc SaveCons {winnum structformat {ctsubformat "zuker"}} {
		#############
		variable DStr

		set seq $DStr($winnum,seq)
		set bplist $DStr($winnum,l_basepairs)
		set id $DStr($winnum,name)
		
		SaveStructFrontend $id $seq $bplist \
			$structformat $ctsubformat
	}
	# SaveCons




    ###   LoadSquFile
    #
    # CURRENTLY NOT USED (and not yet namespace ported)!
    #
#    proc LoadSquFile {squfile} {
#
#        # work: use variables instead of globals
#        global DStr
#        global schmockseq l_basepairs l_prob
#        global headline1 headline2 headline3
#
#        set DStr(directory) [file dirname $squfile]
#
#
#        set schmockseq ""
#        set fid_squ [open $squfile r]
#
#       set clockval [clock format [clock seconds] -format "%b %d, %Y    %R"]
#       set  headline1 [format "PlotStructure of: %s; %s" $squfile $clockval]
#       gets $fid_squ headline2                ;# reading header "FOLD of: ..."
#       gets $fid_squ line                    ;# empty line
#       gets $fid_squ line                    ;# empty line
#       gets $fid_squ headline3                ;# reading header "Length: ... Energy: ..."
#       gets $fid_squ line                    ;# two dots == separator
#       #vputs "header1 = >$headline1<" 5
#       #vputs "header2 = >$headline2<" 5
#       #vputs "header3 = >$headline3<" 5
#       #vputs $line 5
#
#       set probfound 0
#       while {[gets $fid_squ line] != -1} {
#           set numvars [scan $line "%d %s %*d %*d %d %*d %s %s" i nt pairs prob_dummy text]
#           if    { $numvars==5 } {;    # $text=="cleaned"
#                   lappend l_prob 0.0
#           } elseif    { $numvars==4 } {
#               if {[string first "*" $prob_dummy]>=0} {
#                   lappend l_prob 0.0
#               } else {
#                   lappend l_prob $prob_dummy
#                   set probfound 1
#               }
#           } elseif { $numvars==3 } {
#                   lappend l_prob 0.0
#           } else {
#               puts "ERROR in LoadSquFile: >$line<\n"
#               close $fid_squ
#               return
#           }
#           set schmockseq "$schmockseq$nt"
#           lappend l_basepairs $i $pairs
#       }
#       return $probfound
#
#   }
	# LoadSquFile



    ############################################################################
    ###              ###########################################################
    ###   END MISC   ###########################################################
    ###              ###########################################################
    ############################################################################


}
# namespace eval DrawStruct
