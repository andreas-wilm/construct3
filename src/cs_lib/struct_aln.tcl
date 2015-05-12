##############################################################################
#
# struct_aln.tcl - procedures to create a nice structural alignment
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
#  CVS $Id: struct_aln.tcl,v 1.40 2007-10-22 10:43:23 steger Exp $
#

package provide StructAln 0.1

namespace eval StructAln {

    # constants and global declarations
    variable LINE_LENGTH          100
    variable EMPTY_LINES          7
    variable CSTRUCT_HEADER       "Cons.Struct"
    variable CSEQ_HEADER          "Cons.Seq"
    variable PS_WIDTH             16
    variable NUMBERING_BLOCK_TAG  "StructAln_NextBlock"
    variable MAP_VIOL_TAG         "MapViolation"
    # Background colors for nucleotides in ConsensusAlignment
    variable COLOR

    variable Offset

    variable cs_seq       ;# consensus-sequence

    variable prob_est     ;# bit        ;# bool: 1->log_2 / 0->log_e
    # bool: 1->unbiased / 0->mlm prob estimator
    variable work_dir     "[pwd]"
    variable project_name "DUMMY"

    variable struct_repr    ;# (dotbracket|al_num)

    variable nmapviols      ;# number of mappingviolations
    variable displ_mapviols ;# bool

    variable filetypestxt

    # sum of pairs
    variable sop_score      ;# absolute sop score
    variable sop_col_score  ;# 1<=i<=aln_len: sop score pre column
    variable sop_col_max    ;# max sop score for a column
    variable sop_max        ;# max sop score for this alignment length

    variable bplist         ;# a list of basepair partner indices

    # the main arrays seq and bpp
    # have to be handed over by each routine separately,
    # since "Run" cannot upvar it and make it namespace wide visible as "variable"

    # seq ;# array with indices:
    #    aln_len (alignment length == length of all seqs)
    #    n_seq  (number of sequences)
    #    nt,<i> (where 1<=i<=n_seq   sequence nucleotides)
    #    id,<i> (where 1<=i<=n_seq   sequence ids)
    # bpp ;# array with indices:
    #    <i>,prob    (where 1<=i<=aln_len   basepair probability)
    #    <i>,partner (where 1<=i<=aln_len   basepair partner)



    ###   Init
    #
    #
    proc Init {} {
    	########
    	variable Font
    	variable COLOR
    	variable Offset
    	variable filetypestxt
    	variable sop_score
    	variable sop_col_score
    	variable sop_col_max
    	variable sop_max
    	variable displ_mapviols
    	variable nmapviols
    	variable bplist
        
    	set bplist {}
    	
    	set displ_mapviols 1
    	set nmapviols      0
        
    	set CourierMediumNormal "-*-helvetica-bold-r-normal-*-12-120-*-*-*-*-*-*;"

    	set COLOR(NoBp)	white;           # "#0000ffff0000"; # hellgruen
    	set COLOR(CsBp)	"#ffffcccccccc"; # faded pink;      # "#00008000ffff";  # hellblau
    	set COLOR(Bp)   "#ffff66660000"; # orange;          # "#d000f000f000";  # hellgrau
    	set COLOR(SS)   "#bf00ffffbf00"; # "#3333cccc3333"; # spring green;     # "#000000000000";  # white
    	set COLOR(MapViol) "#bfbf3e3effff";# DarkOrchid1

    	set Font(Name)   $CourierMediumNormal
    	set Font(Width)  [expr {-1+int([font measure $Font(Name) G        ])}]
    	set Font(Height) [expr {-1+    [font metrics $Font(Name) -linespace]} ]

    	# In case of proportional fonts the chars should be centered
    	# thus create an offset table
    	#
    	foreach char {A C G U a c g u . ( ) T X N t x n - B D E F H I J K L M O P Q R S V W Y Z b d e f h i j k l m o p q r s v w y z} {
    		set Offset($char) [expr {int(($Font(Width) - [font measure $Font(Name) $char])/2.)}]
    	}

    	set sop_score     0
    	catch {array unset sop_col_score}
    	set sop_col_max   0
    	set sop_max       0

    	### check if needed external commands exist
    	#
    	set external_cmds [list PrintPS ColorProb IsBasepair MouseScroll SumOfPairs \
    						   SumOfPairs CalcBaseOptionsStr Get_ConsBpProb PwIdent \
    						   StructTransform HelixEnd Debugger]

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

    	set filetypestxt {
    		{{Text Files}  {.txt}}
    		{{All Files}       *}
    	}

    	return OK
    }
    # Init




    ###   Run
    #
    # seq must be an array as described above
    # cons_seq = consensus_sequence
    # title = name of this project
    #
    #
    #
    proc Run {alnseq_array cons_seq bpp_array unbiased bit title {show_seq_stats 1} {show_struct_stats 1} {show_patt_stats 0}} {
    	###########################################################
    	variable cs_seq
    	variable CSEQ_HEADER
    	variable prob_est
    	variable struct_repr
    	variable bplist
        
    	upvar $alnseq_array seq
    	upvar $bpp_array    bpp

    	set w_name_aln  ".struct_aln"
    	set w_name_stat ".struct_stat"

    	set cs_seq(id) "$CSEQ_HEADER"
    	set cs_seq(nt) "$cons_seq"

    	set bplist {}
    	for {set i 1} {$i<=[string length $cons_seq]} {incr i} {
    		lappend bplist $bpp($i,partner)
    	}
    	set struct_repr(dotbracket) [StructTransform $bplist $seq(aln_len) "dotbracket"]
    	set struct_repr(al_num) [StructTransform $bplist $seq(aln_len)]

    	set prob_est(bit)      $bit
    	set prob_est(unbiased) $unbiased

    	SetRowheaderFontWidth seq

    	Debugger "struct_aln"

    	CalcSumOfPairs seq


    	
# GST
#puts "tk scaling 1:   [tk scaling]"
set tk_scale [tk scaling]
tk scaling 1.0
#puts "tk scaling 2:   [tk scaling]"
# GST

        set can [OpenStructAlignWindow $w_name_aln  $title seq bpp]
    	CreateNiceAlignment $can $title seq bpp

# GST
#puts "tk scaling 3:   [tk scaling]"
tk scaling $tk_scale
#puts "tk scaling 4:   [tk scaling]"
# GST

    	### do structure statistics
    	#

        set seq_stats "\n---   Sequence Statistics: Disabled   ---\n"
    	set struct_stats "\n---   Structure Statistics: Disabled   ---\n"
        set patt_stats "\n---   Pattern Statistics: Disabled   ---\n"

        set cs_options [CalcBaseOptionsStr]
        
        if {$show_seq_stats} {
            set seq_stats [SeqStatsStr seq bpp]
        }
        if {$show_struct_stats} {
            set struct_stats [CreateStructStats seq bpp]
        }
        if {$show_patt_stats} {
            pattern::Init [array get seq] [array get bpp]
            set patt_stats "\n\n---   Pattern Statistics   ---\n"
            append patt_stats [pattern::PatScanStats]
        }
        
        set stats "$cs_options\n"
        append stats "$seq_stats\n"
        append stats "$struct_stats\n"
        append stats "$patt_stats\n"
        
    	CreateStatsWin $w_name_stat "Statistics" $stats
    }
    # Run




    ###   SetRowheaderFontWidth
    #
    # sets Font(NameSpace) = width for seq-ids
    #
    proc SetRowheaderFontWidth {s} {
    	##########################
    	variable Font
    	variable CSTRUCT_HEADER
    	variable CSEQ_HEADER

    	upvar $s seq

    	set id_font_width 0
    	for {set si 1} {$si<=$seq(n_seq)} {incr si} {
    		set dummy [font measure $Font(Name) "$seq(id,$si) "]
    		if {$id_font_width < $dummy} {
    			set id_font_width $dummy
    		}
    	}
    	# overwrite if smaller then font measure for header
    	set dummy [font measure $Font(Name) "$CSTRUCT_HEADER "]
    	if {$id_font_width < $dummy} {
    		set id_font_width $dummy
    	}
    	set dummy [font measure $Font(Name) "$CSEQ_HEADER "]
    	if {$id_font_width < $dummy} {
    		set id_font_width $dummy
    	}

    	set Font(NameSpace) [expr {int((1*$id_font_width)/$Font(Width))}]
    }
    # SetRowheaderFontWidth





    ###   OpenStructAlignWindow
    #
    #
    proc OpenStructAlignWindow {win  title seq_an bpp_an} {
    	#################################################
    	variable EMPTY_LINES
    	variable Font
    	variable LINE_LENGTH

    	upvar $seq_an seq
    	upvar $bpp_an bpp

    	set outfile "${title}_msf.ps"
        set title "ConStruct structural alignment: $title"

    	set np3          [expr {$seq(n_seq) +1+ $EMPTY_LINES}]
    	set num2         [expr {2+$np3}]
    	set height       [expr {($num2>40) ? 40 : $num2}];
    	set height       [expr {$Font(Height)*$height}]
    	set scrollheight [expr {$Font(Height)*(2+((1+$seq(aln_len)/($LINE_LENGTH-1)))*$np3)}]
    	set       width  [expr {$Font(Width) *(2+$Font(NameSpace)+$LINE_LENGTH)}]
    	set scrollwidth  $width

    	if {[winfo exists $win]} { destroy $win }
    	toplevel	 $win
    	wm title	 $win $title
    	wm iconname  $win $title
    	# FIXME: wm geometry  $win $pos_w(cons)


    	frame $win.menu -relief raised -bd 2
    	pack  $win.menu -side top -fill x -anchor n


    	### File menu
    	#
    	menubutton	 $win.menu.file		  -text File -menu $win.menu.file.menu
    	set m		 $win.menu.file.menu
    	menu		 $m
    	$m add cascade			 -label "Save ..."	-menu $m.save
    	$m add cascade			 -label "Print ..." -menu $m.print
    	$m add separator
    	$m add command			 -label "Close" -command "destroy $win"


    	### File menu, submenu Save

    	menu  $m.save
    	$m.save add command	-label "Consensus" -command "[namespace code SaveCons] cs-consensus"
    	$m.save add cascade	-label "Connect..." -menu $m.ct
    	$m.save add command	-label "Rnaml" -command "[namespace code SaveCons] rnaml"
    	$m.save add command	-label "Stockholm" -command "[namespace code SaveCons] stockholm"
    	
    	menu  $m.ct
    	$m.ct add command -label "Zuker" -command "[namespace code SaveCons] connect zuker"
    	$m.ct add command -label "Gcg" -command "[namespace code SaveCons] connect gcg"
    	$m.ct add command -label "Rnaviz" -command "[namespace code SaveCons] connect rnaviz"

    	### File menu, submenu Print

    	menu  $m.print
    	$m.print add cascade -label "to printer"		-menu $m.printer
    	$m.print add cascade -label "to color printer"	-menu $m.printer_color
    	$m.print add cascade -label "to file"			-menu $m.file
    	$m.print add cascade -label "to screen"			-menu $m.screen


    	### File menu, submenu Print, subsubmenu printer
    	menu  $m.printer

    	$m.printer add command	-label "Full alignment" \
    		-command "[namespace code PrintAlnWin] $win.text full $height $scrollheight $scrollwidth $outfile print"
    	$m.printer add command	-label "View" \
    		-command "[namespace code PrintAlnWin] $win.text part $height $scrollheight $scrollwidth $outfile print"


    	### File menu, submenu Print, subsubmenu printer_color
    	menu  $m.printer_color

    	$m.printer_color add command	-label "Full alignment" \
    		-command "[namespace code PrintAlnWin] $win.text full $height $scrollheight $scrollwidth $outfile printcolor"
    	$m.printer_color add command	-label "View" \
    		-command "[namespace code PrintAlnWin] $win.text part $height $scrollheight $scrollwidth $outfile printcolor"


    	### File menu, submenu Print, subsubmenu file
    	menu  $m.file
    	$m.file add command	-label "Full alignment" \
    		-command "[namespace code PrintAlnWin] $win.text full $height $scrollheight $scrollwidth $outfile save"
    	$m.file add command	-label "View" \
    		-command "[namespace code PrintAlnWin] $win.text part $height $scrollheight $scrollwidth $outfile save"


    	### File menu, submenu Print, subsubmenu screen
    	menu  $m.screen
    	$m.screen add command	-label "Full alignment" \
    		-command "[namespace code PrintAlnWin] $win.text full $height $scrollheight $scrollwidth $outfile screen"
    	$m.screen add command	-label "View" \
    		-command "[namespace code PrintAlnWin] $win.text part $height $scrollheight $scrollwidth $outfile screen"


    	pack         $win.menu.file -side left



    	### Alignment Canvas incl. Scrollbars

    	frame	$win.text
    	pack         $win.text -expand yes -fill both -padx 1 -pady 1
    	canvas	$win.text.text -yscrollcommand	"$win.text.yscroll set" \
    		-xscrollcommand	"$win.text.xscroll set"     \
    		-background white  \
            -borderwidth	2         \
    		-width ${width}p   -height ${height}p     \
    		-scrollregion      "0p 0p ${scrollwidth}p ${scrollheight}p" \
    		-yscrollincrement  $Font(Height)p   \
    		-xscrollincrement  $Font(Width)p    \
    		-relief raised  -bd 2

    	scrollbar   $win.text.xscroll -command "$win.text.text xview" \
    		-highlightthickness 0 -orient horizontal
    	scrollbar   $win.text.yscroll -command "$win.text.text yview" \
    		-highlightthickness 0 -orient vertical

    	grid        $win.text.text -in $win.text -padx 1 -pady 1 \
    		-row 0 -column 0 -rowspan 1 -columnspan 1 -sticky news
    	grid        $win.text.yscroll -in $win.text -padx 1 -pady 1 \
    		-row 0 -column 1 -rowspan 1 -columnspan 1 -sticky news
    	grid        $win.text.xscroll -in $win.text -padx 1 -pady 1 \
    		-row 1 -column 0 -rowspan 1 -columnspan 1 -sticky news
    	grid rowconfig    $win.text 0 -weight 1 -minsize 0
    	grid columnconfig $win.text 0 -weight 1 -minsize 0




    	### setup mouse scrolling support if available
    	#
    	MouseScroll $win.text.text

    	return $win.text.text
    }
    # OpenStructAlignWindow




    ###   CreateNiceAlignment
    #
    #
    proc CreateNiceAlignment {win header seq_an bpp_an} {
    	###############################################
    	upvar $seq_an seq
    	upvar $bpp_an bpp

    	DrawHeader            $win $header

    	DrawNumbering         $win seq

    	DrawSeqIds            $win seq

        DrawColoredSeqs       $win seq bpp

    	DrawColoredConsSeq    $win seq

    	DrawColoredConsStruct $win seq bpp
    }
    # CreateNiceAlignment




    ###   DrawNumbering
    #
    # Write the numbering
    # add tags to the numbering lines to separate the blocks during printing
    #
    #
    proc DrawNumbering {win seq_an} {
    	##########################
    	variable LINE_LENGTH
    	variable EMPTY_LINES
    	variable Font
    	variable NUMBERING_BLOCK_TAG

    	upvar $seq_an seq


    	set np3  [expr {$seq(n_seq)+$EMPTY_LINES}]

    	for {set i 0} {$i<=[expr {($seq(aln_len)-1)/$LINE_LENGTH}]} {incr i} {
    		set line [expr {$Font(Height)*($EMPTY_LINES-1+$i*$np3)}]p
    		set text [expr {$i*$LINE_LENGTH+ 1}]
    		$win create text [expr {$Font(Width)*($Font(NameSpace)+ 0.5)}]p  $line \
    			-text $text -font $Font(Name) \
    			-anchor s -tags [list $NUMBERING_BLOCK_TAG $text]

    		for {set j 10} {$j<=$LINE_LENGTH && [expr {$i*$LINE_LENGTH+$j}]<=$seq(aln_len)} {incr j 10} {
    			set text [expr {$i*$LINE_LENGTH+$j}]
    			$win create text [expr {$Font(Width)*($Font(NameSpace)+$j-0.5)}]p  $line \
    				-text $text -font $Font(Name) -anchor s
    		}
    	}
    	$win create text [expr {$Font(Width) *($Font(NameSpace)+ 0.5)}]p \
    		[expr {$Font(Height)*($EMPTY_LINES-1+$i*$np3)}]p -text " "  \
    		-anchor s -tags [list $NUMBERING_BLOCK_TAG [expr {$i*$LINE_LENGTH+ 1}]]
    }
    # DrawNumbering





    ###   DrawHeader
    #
    #
    proc DrawHeader {win header} {
    	########################
    	variable Font

    	set now  [clock format [clock seconds] -format "%b %d, %Y  %R"]
    	set text [format "Project: %s; Date: %s" $header $now]
    	$win create text 0p [format "%dp" [expr {2*$Font(Height)}]] \
    		-text "$text" -font $Font(Name) -anchor sw
    }
    # DrawHeader




    ###   DrawSeqIds
    #
    #
    proc DrawSeqIds {win seq_an} {
    	########################
    	variable LINE_LENGTH
    	variable EMPTY_LINES
    	variable CSTRUCT_HEADER
    	variable Font
    	variable cs_seq

# GST
#puts "Font(Name):   $Font(Name)"
#puts "Font(Width):  $Font(Width)"
#puts "Font(Height): $Font(Height)"
#puts "tk scaling:   [tk scaling]"
# GST

    	upvar $seq_an seq

    	set pos 0p
    	set np3  [expr {$seq(n_seq)+$EMPTY_LINES}]
    	for {set i 0} {$i<=[expr {($seq(aln_len)-1)/$LINE_LENGTH}]} {incr i} {
    		set dummy [expr {$EMPTY_LINES+$i*$np3}]
    		for {set j 1} {$j<=$seq(n_seq)} {incr j} {
    			$win create text $pos [expr {$Font(Height)*($dummy+$j  )}]p \
    				-text $seq(id,$j) -font $Font(Name) -anchor sw
    		}
    		# consensus seq

    		$win create text $pos [expr {$Font(Height)*($dummy+$j  )}]p \
    			-text "$cs_seq(id)" -font $Font(Name) -anchor sw
    		incr j

    		$win create text $pos [expr {$Font(Height)*($dummy+$j  )}]p \
    			-text "$CSTRUCT_HEADER" -font $Font(Name) -anchor sw
    		incr j

    		$win create text $pos [expr {$Font(Height)*($dummy+$j)}]p \
    			-text "$CSTRUCT_HEADER" -font $Font(Name) -anchor sw
    	}
    }
    # DrawSeqIds




    ###   DrawColoredSeqs
    #
    #
    #
    proc DrawColoredSeqs {win seq_an bpp_an} {
    	###################################
    	variable LINE_LENGTH
    	variable EMPTY_LINES
    	variable cs_seq
    	variable Font
    	variable Offset
    	variable COLOR
    	variable MAP_VIOL_TAG
    	variable NUMBERING_BLOCK_TAG
    	variable displ_mapviols
    	variable nmapviols

    	upvar $seq_an seq
    	upvar $bpp_an bpp


    	set nmapviols 0
    	set tag_prefix "StructAln_Seq"

    	set np3  [expr {$seq(n_seq)+$EMPTY_LINES}]

    	for {set s 1} {$s <= $seq(n_seq)} {incr s} {
    		for {set i 1} {$i <=$seq(aln_len)} {incr i} {

    			set thisnttag "${tag_prefix}_${s}_${i}"

    			###  draw the nt char
    			#
    			#
    			set i_offset [expr {$i-1}]
    			set char [string index $seq(nt,$s) $i_offset]
    			if {$char==[string tolower [string index $cs_seq(nt) $i_offset]]} {
    				set char [string toupper $char]
    			}
    			set pos  [expr {$Font(Width) *($Font(NameSpace) +($i_offset%$LINE_LENGTH)) + $Offset($char)}]
    			set line [expr {$Font(Height)*($EMPTY_LINES+($i_offset/$LINE_LENGTH)*$np3+$s)}]


    			### colorize the basepair types
    			#
    			# chars should be centered, colored rectangle not -> -$Offset($char)
    			#
                set x1 [expr {$pos-$Offset($char)}]p
                set y1 ${line}p
                set x2 [expr {$Font(Width)+$pos}]p
                set y2 [expr {-$Font(Height)+$line}]p

                if {$bpp($i,prob)<=0.0} {
                        set color $COLOR(SS)
                } else {
                    set bpt [BpType $s $i seq bpp]

                    if {$bpt=="NONE"} {
                        set color $COLOR(NoBp)
                    } elseif {$bpt=="CS"} {
                        set color $COLOR(CsBp)
                    } elseif {$bpt=="NORM"} {
                        set color $COLOR(Bp)
                    } else {
                        puts "WARNING: unknown bp-type $bpt, assuming NONE"
                        set color $COLOR(NoBp)
                    }
                }
                $win create rectangle $x1 $y1 $x2 $y2 -fill $color -outline ""; # -tag $col_bp_tag

    			$win create text ${pos}p ${line}p -text $char \
    				-font $Font(Name) -anchor sw \
    				-tag $thisnttag;    # -fill $color

    			if {[IsGap $char]} {
    				continue
    			}

    			### balloon
    			#
    			set degapindex [MapNtIndex $seq(nt,$s) $i ALIGNED_TO_DEGAPPED]
    			if {$degapindex==-1} {
    				p::error "Upss...MapNtIndex failed for seq(nt,$s) ntidx $i"
    			}
    			BalloonHelp::set_balloon $win "Nt $degapindex/$i" $thisnttag

    			### mapping violation from project file with \
    			#   tag MAP_VIOL_TAG and color COLOR(MapViol)
    			#
    			#
    			if {[cs_proj::mapinfo::exists $s] && $degapindex!=-1} {
    				set paired [cs_proj::mapinfo::is_paired $s $degapindex]

    				# set to pair but calculated not to pair
    				if {$paired==1 && $bpp($i,prob)<=0.0} {
    					set mapviol 1
    				# set to pair and calculated to pair
    				} elseif  {$paired==1 && $bpp($i,prob)>0.0} {
    					set settopairwith  [cs_proj::mapinfo::is_paired_with $s $degapindex]
    					set calctopairwith [MapNtIndex $seq(nt,$s) $bpp($i,partner) ALIGNED_TO_DEGAPPED]
    					if {$settopairwith!=0 && $settopairwith!=$calctopairwith} {
    						set mapviol 1
    					} else {
    						set mapviol 0
    					}
    					# set not to pair but calculated to pair
    				} elseif {[cs_proj::mapinfo::is_not_paired $s $degapindex]==1 && $bpp($i,prob)>0.0} {
    					set mapviol 1
    				} else {
    					set mapviol 0
    				}

    				if {$mapviol} {
    					incr nmapviols
    					if {$displ_mapviols} {
    						$win create rectangle $x1 $y1 $x2 $y2 \
    							-fill $COLOR(MapViol) -outline "" \
    							-tag "$MAP_VIOL_TAG ${MAP_VIOL_TAG}_seq${s}_nt${i}"
    					}
    				}
    			}
    		}
    	}
    	$win raise $MAP_VIOL_TAG $NUMBERING_BLOCK_TAG

    }
    # DrawColoredSeqs




    ###   CalcSumOfPairs
    #
    #
    proc CalcSumOfPairs {seq_an} {
    	#######################
    	variable sop_score
    	variable sop_col_score
    	variable sop_col_max
    	variable sop_max

    	upvar $seq_an seq

    	set sop_col_max [expr {($seq(n_seq)*($seq(n_seq)-1))/2}]
    	set sop_max     [expr {$sop_col_max*$seq(aln_len)}]

    	set sop_score   [SumOfPairs  sop_col_score]
    }
    # CalcSumOfPairs



    ###   DrawColoredConsSeq
    #
    #
    proc DrawColoredConsSeq {win seq_an} {
    	################################
    	variable cs_seq
    	variable Font
    	variable LINE_LENGTH
    	variable Offset
    	variable EMPTY_LINES
    	variable sop_score
    	variable sop_col_score
    	variable sop_col_max

    	upvar $seq_an seq

    	set nt_tag_prefix "StructAln_CsSeq_Nt_"

    	set np3  [expr {$seq(n_seq)+$EMPTY_LINES}]



    	for {set i 0} {$i < $seq(aln_len)} {incr i} {
    		set char [string index $cs_seq(nt) $i]
    		set pos  [expr {$Font(Width) *($Font(NameSpace) +($i%$LINE_LENGTH)) + $Offset($char)}]
    		set line [expr {$Font(Height)*($EMPTY_LINES+($i/$LINE_LENGTH)*$np3+$seq(n_seq)+1)}]

    		###  colorize by column cost sum of pairs cost
    		#

    		set rel_col_score [expr {$sop_col_score([expr {$i+1}])/double($sop_col_max)}]

    		# puts "rel. column score([expr {$i+1}]) = $rel_col_score"

            set x1 [expr {$pos-$Offset($char)}]p
    		set coords [list ${x1}p ${line}p \
    							[expr {$Font(Width)+$pos}]p [expr {-$Font(Height)+$line}]p]

    		if {$rel_col_score!=0.0} {
    			$win create rectangle $coords -fill [ColorProb $rel_col_score] \
    				-outline "";    # -tag $col_tag
    		} else {
                $win create rectangle $coords -fill white \
    				-outline "";    # -tag $col_tag
            }

    		$win create text ${pos}p ${line}p -text $char \
    			-font $Font(Name) -anchor sw \
    			-tag "${nt_tag_prefix}[expr {$i+1}]"
    	}
    }
    # DrawColoredConsSeq




    ###   DrawColoredConsStruct
    #
    #
    proc DrawColoredConsStruct {win seq_an bpp_an} {
    	##########################################
    	variable cs_seq
    	variable Font
    	variable LINE_LENGTH
    	variable Offset
    	variable EMPTY_LINES
    	variable COLOR
    	variable struct_repr

    	upvar $seq_an seq
    	upvar $bpp_an bpp

    	set char_tag   "StructAln_CsStruct"

    	set np3  [expr {$seq(n_seq)+$EMPTY_LINES}]


    	### Colorize both ! consensus structure lines
    	#
    	for {set i 0} {$i < $seq(aln_len)} {incr i} {
   		    set x1 [expr {$Font(Width) *($Font(NameSpace) +($i%$LINE_LENGTH))}]
    		set y1 [expr {$Font(Height)*($EMPTY_LINES+($i/$LINE_LENGTH)*$np3+$seq(n_seq)+1)}]
    		set x2 [expr {$Font(Width) +$x1}]p
    		set y2 [expr {$Font(Height)*2+$y1}]p
    		set coords [list ${x1}p ${y1}p ${x2}p ${y2}p]
    		set prob $bpp([expr {$i+1}],prob)
    		if {$prob<=0.0} {
    			$win create rectangle $coords -fill $COLOR(SS) -outline ""
    		} else {
    			$win create rectangle $coords -fill [ColorProb $prob] -outline ""
    		}
    	}

    	# Create and colorize the consensus sequence from
    	# white (= no bp) to red (= always a bp)
    	# according to the probabilities of base pairing
    	#
    	for {set i 0} {$i < $seq(aln_len)} {incr i} {

    		# dotbracket
    		#
    		set char [string index $struct_repr(dotbracket) $i]
    		set pos  [expr {$Font(Width) *($Font(NameSpace) +($i%$LINE_LENGTH)) + $Offset($char)}]
    		set line [expr {$Font(Height)*($EMPTY_LINES+($i/$LINE_LENGTH)*$np3+$seq(n_seq)+2)}]

    		$win create text ${pos}p ${line}p -text $char \
    			-font $Font(Name) -anchor sw \
    			-tag "$char_tag"

    		# alnum
    		#
    		set char  [string index $struct_repr(al_num) $i]
    		set line  [expr {$line+$Font(Height)}]

    		$win create text ${pos}p ${line}p -text $char \
    			-font $Font(Name) -anchor sw \
    			-tag "$char_tag"

    	}
    }
    # DrawColoredConsStruct






    ###   PrintAlnWin
    #
    #
    #
    proc PrintAlnWin {w full height scrollheight scrollwidth outfile save_print} {
    	########################################################################
    	variable NUMBERING_BLOCK_TAG
    	variable PS_WIDTH

# GST
# puts "tk scaling 1x:   [tk scaling]"
set tk_scale [tk scaling]
tk scaling 1.0
# puts "tk scaling 2x:   [tk scaling]"
# GST

    	set length_A4_cm  25.0


    	if {$full=="full"} {
    		set oldview [$w.text yview];				# remember the viewport
    		$w.text yview moveto 0.;					# view the top part of the canvas
    		set totalheight   [expr {1.*$scrollheight*$PS_WIDTH/$scrollwidth}]


    		# canvas is heigher than DIN A4 page is long;
    		# thus print each block on a separate page
    		if {$totalheight>$length_A4_cm} {
    			set i 0
    			# get the canvas coords of the numbering lines, which start a block
    			foreach tag [$w.text find withtag $NUMBERING_BLOCK_TAG] {
    				set range($i) [lrange [$w.text coords $tag] end end]
    				incr i
    			}
    			# number of coord pairs
    			set last    [expr {$i-1}]
    			# first block starts with the header; thus begin printing with 0. but not $range(0)
    			set start   0.;
    			for {set j 1} {$j <= $last} {incr j} {
    				$w.text yview moveto $start
    				set end        [expr {$range($j)/$range($last)}]
    				set height_pnt [expr {(($totalheight*($end-$start))/$PS_WIDTH)*$scrollwidth}]

    				PrintPS  $w.text [format "%s_%d" $outfile $j] $save_print -height ${height_pnt}p -width ${scrollwidth}p

    				set start $end
    			}
    			# canvas fits on a single DIN A4 page
    		} else {

    			PrintPS  $w.text $outfile $save_print -height ${scrollheight}p -width ${scrollwidth}p
    		}
    		# restore the viewport
    		$w.text yview moveto [lindex $oldview 0];

    		# print only the visible part
    	} else {
    		PrintPS  $w.text $outfile $save_print;
    	}

# GST
#puts "tk scaling 3x:   [tk scaling]"
tk scaling $tk_scale
#puts "tk scaling 4x:   [tk scaling]"
# GST
    }
    # PrintAlnWin




    ###   SaveCons
    #
    # save consensus
    #
    proc SaveCons {structformat {ctsubformat "zuker"}} {
    	#############
    	variable cs_seq
    	variable bplist

    	SaveStructFrontend $cs_seq(id) $cs_seq(nt) $bplist \
    		$structformat $ctsubformat
    }
    # SaveCons




    ###   CreateStructStats
    #
    #
    proc CreateStructStats {seq_an bpp_an} {
    	##################################
    	variable prob_est
    	variable cs_seq
    	variable nmapviols
    	variable bplist
    	upvar $seq_an seq
    	upvar $bpp_an bpp

    	# stats for one helix
    	set helix_NoBp   0
    	set helix_Bp     0   ;# normal bps
    	set helix_CsBp   0   ;# consensus bps, e.g compensatory changes (used for gemoetric mean comp.)
    	set helix_csprob 1.0 ;# product of consensus probability (mic+td) (used for gemoetric mean comp.)
    	set helix_tdprob 1.0 ;# product of consensus td probability (used for gemoetric mean comp.)
    	set helix_mic    1.0 ;# product of mutual informatic content (used for gemoetric mean comp.)
    	set helix_length 0

    	set overall_struct(NoBp)   0 ;# as above
    	set overall_struct(Bp)     0
    	set overall_struct(CsBp)   0
    	set overall_struct(csprob) 0.0
    	set overall_struct(tdprob) 0.0
    	set overall_struct(mic)    0.0
    	set overall_struct(length) 0

    	set helixChar 97; # a

    	Upe::SetupChi2Distr


    	set TextAlignStatistics "---   Structure Statistics ---\n\n"

    	append TextAlignStatistics "\nNumber of Mapping Violations: $nmapviols\n\n"

    	append TextAlignStatistics [format "Helix %c:\n" $helixChar]
    	set helixheader "              BP     NoBP   BP CsBP"
    	append helixheader "   Prob(CS) Prob(TD) I(x,y)"
    	if {$prob_est(bit)} {
    		append helixheader "  R1    R2     Pairs\n"
    	} else {
    		append helixheader "     chi^2(df,p)       R1    R2     Pairs\n"
    	}
    	append TextAlignStatistics $helixheader

    	for {set n 1} {$n<$seq(aln_len)} {incr n} {

    		if {$bpp($n,prob)<=0.0} {
    			continue
    		}
    		# don't work twice
    		if {$bpp($n,partner)<$n} {
    			continue
    		}

    		if {$helixChar==123} { set helixChar 65 };	# A follows z
    		if {$helixChar== 91} { set helixChar 97 };	# a follows Z

    		### setup listNt
    		#
    		# will be shortened to {A U G C} if no gap present
    		set listNt [list . N A U G C]
    		# setup arrays
    		foreach a $listNt {;
    			set fi($a) 0
    			set fj($a) 0
    			foreach b $listNt {
    				set f($a,$b) 0
    			}
    		}

    		# count the numbers of nts at pos x, at pair pos
    		# and the numbers of nt combinations
    		for {set s 1} {$s<=$seq(n_seq)} {incr s} {
    			set idx [expr {$n-1}]
    			set typ_x [string toupper [string index $seq(nt,$s) $idx]]
    			# unknown character -> n
    			if {[lsearch -exact $listNt $typ_x]==-1} {set typ_x N}

    			set idx [expr {$bpp($n,partner)-1}]
    			set typ_y [string toupper [string index $seq(nt,$s) $idx]]
    			# unknown character -> n
    			if {[lsearch -exact $listNt $typ_y]==-1} {set typ_y N}

    			incr fi($typ_x)
    			incr fj($typ_y)
    			incr f($typ_x,$typ_y)
    		}


    		set listNt [UpdateNtList fi fj f $listNt]


    		set probstat [Upe::upe  $seq(n_seq)     \
    						  fi fj f $listNt \
    						  $prob_est(bit) $prob_est(unbiased)]

    		###  count bp types (double! see above)
    		#
    		set NoBp 0
    		set   Bp 0
    		set CsBp 0
    		#
    		for {set s 1} {$s<=$seq(n_seq)} {incr s} {

    			set bpt [BpType $s $n seq bpp]

    			if {$bpt=="NONE"} {
    				incr NoBp
    				incr helix_NoBp
    			} elseif {$bpt=="CS"} {
    				incr CsBp
    				incr helix_CsBp
    			} elseif {$bpt=="NORM"} {
    				incr Bp
    				incr helix_Bp
    			} else {
    				puts "WARNING: unknown bp-type $bpt, assuming NONE"
    				incr NoBp
    				incr helix_NoBp
    			}
    		}

    		set cs_nt_i [string index $cs_seq(nt) [expr {$n-1}]]
    		set cs_nt_j [string index $cs_seq(nt) [expr {$bpp($n,partner)-1}]]

    		#                                       nti ntj nt nt  No  BP  CS
    		append TextAlignStatistics \
    			[format "\t %3d:%3d=%s:%s %4d %4d %4d" \
    				 $n $bpp($n,partner) $cs_nt_i $cs_nt_j \
    				 $NoBp $Bp $CsBp]

    		if {$n>$bpp($n,partner)} {
    			set tdprob [Get_ConsBpProb $n $bpp($n,partner)]
    		} else {
    			set tdprob [Get_ConsBpProb $bpp($n,partner) $n]
    		}
    		if {$prob_est(bit)} {
    			# csprob tdprob   I     R1   R2  Pairs
    			append TextAlignStatistics \
    				[format "     %.3f    %.3f  %.3f   %.3f %.3f  %s\n" \
    					 $bpp($n,prob) $tdprob \
    					 [lindex $probstat 0] [lindex $probstat 2] [lindex $probstat 3] \
    					 [lindex $probstat 4] [lindex $probstat 4]]
    		} else {
    			#  csprob tdprob    I  chi   R1  R2  Pairs
    			append TextAlignStatistics \
    				[format "     %.3f    %.3f  %.3f%-23s %.3f %.3f  %s\n" \
    					 $bpp($n,prob) $tdprob \
    					 [lindex $probstat 0] [lindex $probstat 1] [lindex $probstat 2] \
    					 [lindex $probstat 3] [lindex $probstat 4]]
    		}


    		set  helix_csprob [expr {$helix_csprob * $bpp($n,prob)}]
    		set  helix_tdprob [expr {$helix_tdprob * $tdprob}]
    		set  helix_mic    [expr {$helix_mic    * [lindex $probstat 0]}]
    		incr helix_length

    		set  overall_struct(csprob) [expr {$overall_struct(csprob) + $bpp($n,prob)}]
    		set  overall_struct(tdprob) [expr {$overall_struct(tdprob) + $tdprob}]
    		set  overall_struct(mic)    [expr {$overall_struct(mic)    + [lindex $probstat 0]}]


    		if {[HelixEnd $n $seq(aln_len) $bplist]} {

    			if {$bpp($n,partner)>$n} {
    				incr helixChar
    			}
    			append TextAlignStatistics \
    				[format "       Helix len=%3d %4d %4d %4d     %.3f    %.3f  %.3f (geometric means)\n" \
    					 $helix_length $helix_NoBp $helix_Bp $helix_CsBp \
    					 [expr {pow($helix_csprob, 1./$helix_length)}] \
    					 [expr {pow($helix_tdprob, 1./$helix_length)}] \
    					 [expr {pow($helix_mic,    1./$helix_length)}] ]

    			append TextAlignStatistics [format "\nHelix %c:\n" $helixChar]

    			append TextAlignStatistics $helixheader

    			incr overall_struct(length) $helix_length
    			incr overall_struct(NoBp)   $helix_NoBp
    			incr overall_struct(Bp)     $helix_Bp
    			incr overall_struct(CsBp)   $helix_CsBp

    			set helix_NoBp   0
    			set helix_Bp     0
    			set helix_CsBp   0
    			set helix_csprob 1.
    			set helix_mic    1.
    			set helix_length 0
    		}
    	}

    	Upe::UnsetChi2Distr


    	# The last two lines so far are the table headers of the next (not present) helix
    	# -> delete
    	set TextAlignStatistics [string range $TextAlignStatistics 0 [expr {[string last "Helix" $TextAlignStatistics]-2}]]

    	if {$overall_struct(length)==0} {
    		set mean(csprob) 0.0
    		set mean(tdprob) 0.0
    		set mean(mic)    0.0
    	} else {
    		set mean(csprob) [expr {$overall_struct(csprob) / $overall_struct(length)}]
    		set mean(tdprob) [expr {$overall_struct(tdprob) / $overall_struct(length)}]
    		set mean(mic)    [expr {$overall_struct(mic)    / $overall_struct(length)}]
    	}

    	# append overall structure stat
    	append TextAlignStatistics "\nOverall:\n             "
    	append TextAlignStatistics [format "len=%3d %4d %4d %4d     %.3f    %.3f  %.3f (arithmetic means)\n" \
    									$overall_struct(length) $overall_struct(NoBp) \
    									$overall_struct(Bp)     $overall_struct(CsBp) \
    									$mean(csprob) $mean(tdprob) $mean(mic)]

    	return $TextAlignStatistics
    }
    # CreateStructStats





    ###   CreateStatsWin
    #
    #
    proc CreateStatsWin {w header stat_text} {
    	##########################################

#        set ft_error -*-courier-medium-r-*-*-12-120-*-*-*-*-*-*
        set ft_error -*-courier-medium-r-*-*-18-180-*-*-*-*-*-*
    	if {![winfo exists $w]} {
    		toplevel 	$w
    		wm title    $w $header
    		wm iconname $w $header

    		frame     $w.buttons
    		pack      $w.buttons           -side bottom  	-fill x
            button    $w.buttons.save      -text Save       -command  "[namespace code SaveStructStats] $w"
    		button    $w.buttons.print     -text Print      -command {puts "PLEASE IMPLEMENT ME"}
    		button    $w.buttons.dismiss   -text Dismiss    -command "destroy $w"
    		pack      $w.buttons.save     \
    		          $w.buttons.print    \
    		          $w.buttons.dismiss  -side left -expand 1
    		scrollbar $w.scrolly           -command "$w.text yview"
    		scrollbar $w.scrollx           -orient horiz    -command "$w.text xview"
    		pack      $w.scrolly           -side right		-fill y
    		pack      $w.scrollx           -side bottom  	-fill x
    		text      $w.text              -height 40 -width  90              \
    		                               -yscrollcommand "$w.scrolly set"   \
    		                               -xscrollcommand "$w.scrollx set"   \
    		                               -wrap none  -setgrid 1             \
    		                               -font $ft_error
    		pack      $w.text              -side left       -expand 1   -fill both
    	} else {
    		wm title     $w $header
    		wm iconname  $w $header
    		wm deiconify $w
    		raise        $w
    	}
    	$w.text config -state normal
    	$w.text delete 1.0 end


    	$w.text insert end $stat_text
    	$w.text config -state disabled
    }
    # CreateStatsWin



    ###   SaveStructStats
    #
    #
    proc SaveStructStats {win} {
    	######################
    	variable work_dir
    	variable project_name
    	variable filetypestxt

    	set init_file "${project_name}_structstat.txt"


    	set filename [tk_getSaveFile \
    					  -filetypes $filetypestxt       \
    					  -initialfile "$init_file"    \
    					  -initialdir  $work_dir \
    					  -title       "Save StructStats to file"]
    	update idletasks
    	if {$filename == ""} {
    		return
    	}
    	DumpTextWindowContent $win $filename
    }
    # SaveStructStats




    ###   UpdateNtList
    #
    #
    proc UpdateNtList {freq_i freq_j freq listNt} {
    	#########################################
        upvar $freq_i fi
        upvar $freq_j fj
        upvar $freq   f

        set noGap 1
        set noN 1

        # check for gaps; if none remove from list
        foreach a $listNt {
            if { $f(.,$a)!=0} {set noGap 0}
        }
        if {$fi(.)!=0} {set noGap 0}
        if {$fj(.)!=0} {set noGap 0}


        # check for N; if none present remove from list
        foreach a $listNt {
            if { $f(N,$a)!=0} {set noN 0}
        }
        if {$fi(N)!=0} {set noN 0}
        if {$fj(N)!=0} {set noN 0}

        # now remove
        if {$noGap && $noN} {
            set listNt [list A U G C]
        } elseif {$noGap  		} {
            set listNt [list N A U G C]
        } elseif {  		  $noN} {
            set listNt [list . A U G C]
        }

    	return $listNt
    }
    # UpdateNtList




    ###   SeqStatsStr
    #
    #
    proc SeqStatsStr {seq_an bpp_an} {
    	############################
    	variable struct_repr
    	variable sop_score
    	variable sop_col_score
    	variable sop_col_max
    	variable sop_max
    	variable cs_seq

    	upvar $seq_an seq
    	upvar $bpp_an bpp

    	set non_helix_c  "." ;# char in struct_repr representing single strand

    	set n_ssreg 0 ;# number of unpaired regions
    	set prev_paired -1
    	set this_sop(ss) 0
    	set this_sop(ds) 0
    	set this_len(ss) 0
    	set this_len(ds) 0

    	set    seq_stat "---   Sequence Statistics   ---\n\n"
    	append seq_stat [format "Avg. pairwise sequence identity = %.2f %%\n\n\n\n" [expr {100.0*[PwIdent]}]]

    	append seq_stat "Sum-Of-Pairs Score (SOP):\n"
    	append seq_stat "   absolute = $sop_score\n"
    	append seq_stat [format "   relative = %.3f\n\n" [expr {double($sop_score)/double($sop_max)}]]


    	for {set i 1} {$i<=$seq(aln_len)} {incr i} {

    		set helix_c [string index $struct_repr(al_num) [expr {$i-1}]]

    		if {$helix_c==$non_helix_c} {
    			set paired 0
    		} else {
    			set paired 1
    		}

    		if {$paired} {
    			incr this_sop(ds) $sop_col_score($i)
    			incr this_len(ds)
    		} else {
    			if {$prev_paired==1 || $i==1} {
    				incr n_ssreg

    				set ss_reg($n_ssreg,sop)   $sop_col_score($i)
    				set ss_reg($n_ssreg,len)   1
    				set ss_reg($n_ssreg,start) $i
    			} else {
    				incr ss_reg($n_ssreg,sop) $sop_col_score($i)
    				incr ss_reg($n_ssreg,len)

    			}
    			incr this_sop(ss) $sop_col_score($i)
    			incr this_len(ss)
    		}
    		set prev_paired $paired
    	}

    	### sop ss vs. ds
    	#
    	#
    	append seq_stat "Paired vs. unpaired regions\n"
    	append seq_stat "             SOP(abs.)   SOP(rel.)    len\n"

    	set max(ss) [expr {double($sop_col_max)*double($this_len(ss))}]
    	set max(ds) [expr {double($sop_col_max)*double($this_len(ds))}]

    	if {$max(ds)>0.0} {
    		set this_rel_sop(ds) [expr {double($this_sop(ds))/$max(ds)}]
    	} else {
    		set this_rel_sop(ds) 0.0
    	}
    	if {$max(ss)>0.0} {
    		set this_rel_sop(ss) [expr {double($this_sop(ss))/$max(ss)}]
    	} else {
    		set this_rel_sop(ss) 0.0
    	}

    	append seq_stat [format "   unpaired   %6d       " $this_sop(ss)]
    	append seq_stat [format "%.3f     %4d\n" $this_rel_sop(ss) $this_len(ss)]

    	append seq_stat [format "     paired   %6d       " $this_sop(ds)]
    	append seq_stat [format "%.3f     %4d\n" $this_rel_sop(ds) $this_len(ds)]


    	### sop per ss
    	#
    	append seq_stat "\nPer unpaired region:\n"
    	append seq_stat "   start:end   SOP(abs.)  SOP(rel.)   Cons.\n"
    	for {set i 1} {$i<=$n_ssreg} {incr i} {
    		set    txt "    "
    		append txt [format "%4d:" $ss_reg($i,start)]
    		set    end [expr {$ss_reg($i,start)+$ss_reg($i,len)-1}]
    		append txt [format "%4d" $end]
    		append txt [format "   %6d" $ss_reg($i,sop)]
    		set    this_max [expr {double($sop_col_max)*double($ss_reg($i,len))}]
    		if {$this_max>0.0} {
    			set this_rel [expr {double($ss_reg($i,sop))/$this_max}]
    		} else {
    			set this_rel 0.0
    		}
    		append txt [format "      %.3f" $this_rel]

    		append txt [format "     %s" [string range $cs_seq(nt) \
    										  [expr {$ss_reg($i,start)-1}] \
    										  [expr {$end-1}]]]
    		append seq_stat "$txt\n"
    	}

    	append seq_stat "\n\n"
    	return $seq_stat
    }
    # SeqStatsStr





    ###   BpType
    #
    # detects basepair type
    #
    # return:
    #   NONE:   no basepair
    #   CS:     consensus basepair
    #   NORM:   normal basepair
    #
    proc BpType {seq_idx nt_idx seq_an bpp_an} {
    	#####################################
    	variable cs_seq

    	upvar $seq_an seq
    	upvar $bpp_an bpp

    	set cs_nt(i) [string index $cs_seq(nt)       [expr {$nt_idx-1}]]
    	set cs_nt(j) [string index $cs_seq(nt)       [expr {$bpp($nt_idx,partner)-1}]]
    	set nt(i)    [string index $seq(nt,$seq_idx) [expr {$nt_idx-1}]]
    	set nt(j)    [string index $seq(nt,$seq_idx) [expr {$bpp($nt_idx,partner)-1}]]

    	set cs_nt(i) [string tolower $cs_nt(i)]
    	set cs_nt(j) [string tolower $cs_nt(j)]
    	set nt(i)    [string tolower $nt(i)]
    	set nt(j)    [string tolower $nt(j)]

    	if {[IsBasepair $nt(i) $nt(j)]==0} {
    		return NONE

    	} else {
    		# Count "consensus base pair"="conserved base pair change"=CsBp only
    		# if a WC is exchanged into another WC;
    		# A:U <-> U:A => CsBp
    		# G:U <-> G:C =>   Bp

    		if {(($nt(i)!=[string tolower $cs_nt(i)]) && ($cs_nt(j)!=[string tolower $nt(j)]))} {
    			return CS
    		} elseif {[IsBasepair $cs_nt(i) $cs_nt(j)]==0} {
    			return CS

    		} else {
    			return NORM
    		}
    	}
    }
    # BpType
}
# namespace eval StructAln
