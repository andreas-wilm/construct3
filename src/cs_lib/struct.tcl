##############################################################################
#
# struct.tcl - structure routines
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
#  CVS $Id: struct.tcl,v 1.18 2007-10-22 10:43:23 steger Exp $
#



###   OptimalConStruct
#
# Arguments:
#   display_method: drawstruct|circles|struct_aln
#
#
proc OptimalConStruct {display_method} {
    ##################################
    global opts      ;# r
    global seq       ;# r
    global selection ;# r
    global mic_is_computed_and_valid ;# r



    if {$display_method!="drawstruct" && \
    		$display_method!="circles" && \
    		$display_method!="struct_aln" &&\
			$display_method!="struct_logo"} {
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
    OptimalStructure  bpp   $opts(td,factor)    $opts(mic,factor) \
                            $opts(td,threshold) $opts(mic,Colormap,0) \
                            $opts(struct,allow_single_bp)
    #
    #####
    Debugger "struct"

    # debug
    # foreach an [lsort [array names bpp]] {
    #     puts "OptimalConStruct: bpp($an)=$bpp($an)"
    # }


    if {$display_method=="drawstruct"} {

        set h1 "Optimal Consensus Structure"
        set h2 "Project: $cs_proj::proj(name)/$seq_id"
        set h3 "[clock format [clock seconds] -format "%B %d, %Y"]"
        #
        set header [list $h1 $h2 $h3]

        DrawStructFrontend bpp $seq_nt $header $cs_proj::proj(name)

    } elseif {$display_method=="circles"} {

        set header_l1 "ConsensusCircles of $cs_proj::proj(name)/$seq_id"
        set header_l2 " "

        CirclesFrontend bpp  $seq_id  $seq_nt  "c"  $header_l1 $header_l2

    } elseif {$display_method=="struct_aln"} {
        StructAlnFrontend bpp

	} elseif {$display_method=="struct_logo"} {
        StructLogoFrontend bpp

    } else {
        p::error "Invalid display_method=$display_method"
    }

    SetCursor normal

}
# OptimalConStruct




###   SuboptimalConStruct
#
#
proc SuboptimalConStruct {display_method} {
########################
    global opts
    global FONT

    if {$display_method!="drawstruct" && \
    		$display_method!="circles" && \
    		$display_method!="struct_aln" && \
			$display_method!="struct_logo"} {
        p::error "Invalid display_method \"$display_method\""
        return
    }


    if {$opts(mic,factor)>0.0} {
        if {[Validate_MIC]=="abort"} {
            return
        }
    }
    SetCursor busy




    #####
    #
    set n_bp_maxs [SuboptimalPairlist pairlist \
    				   $opts(td,factor) $opts(mic,factor) \
    				   $opts(td,threshold) $opts(mic,Colormap,0) \
    				   $opts(suboptstruct,number)]
    #
    ###
    Debugger "struct"


    # debug
    # foreach an [lsort -dictionary [array names pairlist]] {
    #    p::fixme "pairlist($an)= $pairlist($an)"
    # }

    p::debug "got n_bp_maxs=$n_bp_maxs from SuboptimalPairlist"


    set w .subbuttons
    if {[winfo exists $w]} {
        destroy $w
    }
    toplevel     $w
    wm title     $w "Start-Bp's for Backtrack"
    wm iconname  $w "Start-Bp's for Backtrack"




    ###   label
    #
    label $w.label -wraplength 200 -justify left \
    	-font [format $FONT(dp_nm,format) 120]  \
    	-text "Select a basepair from desired maximum to start Backtrack"
    pack  $w.label -side top -padx 5 -pady 5


    ###   bp-buttons
    #
    #
    frame $w.base -borderwidth 2 -relief sunken
    set b $w.base
    pack  $b
    frame $b.body
    set t $b.body.text
    text  $t -yscrollcommand "$b.scrolly set" \
             -width 30  -height 24 -borderwidth 0 \
             -wrap word -highlightthickness 0 -setgrid true
    pack $t  -expand  yes -fill both
    scrollbar $b.scrolly -command "$t yview"
    pack $b.scrolly $b.body -side right  -fill y


    for {set maxctr 1} {$maxctr<=$n_bp_maxs} {incr maxctr} {
        set relprob [format "%0.3f" $pairlist($maxctr,rel_prob)]
        p::debug "creating max button for maxctr=$maxctr relprob=$relprob"
        menubutton $t.$maxctr -text "Maximum $maxctr (prob=[expr {$relprob*100}]%)" \
                              -direction right -menu $t.$maxctr.m -relief raised -width 30
        menu $t.$maxctr.m -tearoff 0 -disabledforeground black

        foreach bp [lsort -dictionary $pairlist($maxctr,pairs)] {
            p::debug "creating bp button for maxctr=$maxctr: bp=$bp"
            scan $bp "%d:%d" nti ntj
            $t.$maxctr.m add command -label "Bp ${nti}:${ntj}" \
    			-command "DrawSubStruct $display_method $nti $ntj $maxctr $relprob"
        }
    }

    for {set maxctr 1} {$maxctr<=$n_bp_maxs} {incr maxctr} {
        $t window create end -window $t.$maxctr -padx 3 -pady 2
        $t insert end "\n"
    }



    ## buttons

    frame $w.buttons
    pack $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text "Dismiss"   -command "destroy $w"
    pack $w.buttons.dismiss -side left -expand 1 -pady 5


    SetCursor normal
}
# SuboptimalConStruct




###   DrawSubStruct
#
#
proc DrawSubStruct {display_method nt_i nt_j maxbpctr rel_prob} {
##################
    global opts
    global selection
    global seq

    SetCursor busy

    #####
    #
    SubBacktrack  bpp  $nt_i $nt_j $opts(struct,allow_single_bp)
    #
    #####

    if {$opts(displ,use_cons_seq)==1} {
        set seq_nt  [Get_ConsSeq]
        set seq_id "Consensus Sequence"
    } else {
        set seq_nt  $seq(nt,$selection(seq_no))
        set seq_id  $seq(id,$selection(seq_no))
    }


    if {$display_method=="drawstruct"} {

        set    h1 "SubOptimal Consensus Structure"
        set    h2 "Start Bp.=$nt_i:$nt_j, Max.=$maxbpctr, Prob.=[expr {$rel_prob*100}]%"
        append h2 "\nProject: $cs_proj::proj(name)/$seq_id"
        set    h3 "[clock format [clock seconds] -format "%B %d, %Y"]"
        #
        set header [list $h1 $h2 $h3]

        DrawStructFrontend bpp $seq_nt $header $cs_proj::proj(name)

    } elseif {$display_method=="circles"} {

        set header_l1 "SubOptimal of $cs_proj::proj(name)/$seq_id"
        set header_l2 "Start Bp.=$nt_i:$nt_j, Max.=$maxbpctr, Prob.=[expr {$rel_prob*100}]%"

        CirclesFrontend bpp  $seq_id  $seq_nt  "c"  $header_l1 $header_l2

    } elseif {$display_method=="struct_aln"} {
        StructAlnFrontend bpp

	} elseif {$display_method=="struct_logo"} {
        StructLogoFrontend bpp

    } else {
        p::error "Invalid display_method=$display_method"
    }


    SetCursor normal
}
#  DrawSubStruct





###   remGaps
#
#
#
proc remGaps {seq1 l_basepairs1 l_prob1} {
############
    set seq ""
    set l_basepairs {}
    set l_prob {}
    set j 0
    for {set i 0} {$i<[llength $l_basepairs1]} {incr i} {
        incr j
        set a_seq($j)  [string index $seq1   $i]
        set a_bp($j)   [lindex $l_basepairs1 $i]
        set a_prob($j) [lindex $l_prob1      $i]
    }
    set j 0
    for {set i 1} {$i<=[llength $l_basepairs1]} {incr i} {
        if {$a_seq($i)!="-"} {;                             # aktuelle Base i
            incr j
            # puts -nonewline [format "%2d %2d" $i $j]
            # puts -nonewline [format " %4s" $a_seq($i)]
            set seq "${seq}$a_seq($i)";                     #   an Sequenz anhaengen (an Pos. j)
            if {$a_bp($i)>$i} {;                            #   Basenpaar i:k (k>i)
                if {$a_seq($a_bp($i))!="-"} {;              #       Partner k ist Base
                    set a_bp($a_bp($i)) $j;                 #           beim Partner k schon die Base i (als j) als Partner eintragen
                    # puts        [format " a_bp(%2d)=%2d" $a_bp($i) $j]
                } else {;                                   #       Partner k ist Gap
                    set a_bp($j) 0;                         #       also kein Basenpaar
                    # puts "     Complement is gap"
                }
            } elseif {$a_bp($i)==0} {;                      #   i ist nicht basengepaart
                set a_bp($j) 0;                             #       und ist's auch weiterhin nicht
                # puts            [format " a_bp(%2d)=%2d" $a_bp($j) 0]
            } else {;                                       #   Basenpaar i:k (k<i)
                set a_bp($j) $a_bp($i);                     #       bei sich selbst eintragen
                set a_bp($a_bp($j)) $j;                     #       beim Partner eintragen
                # puts -nonewline [format " a_bp(%2d)=%2d" $j        $a_bp($i)]
                # puts            [format " a_bp(%2d)=%2d" $a_bp($j) $j]
            }
            lappend l_prob $a_prob($i)
        } elseif {$a_bp($i)!=0 &&   \
                  $a_bp($i)>$i} {;                          # aktuelle Base ist Gap, aber basengepaart
            set a_bp($a_bp($i)) 0;                          #   beim Partner loeschen
            # puts                [format " a_bp(%2d)=%2d     complement is gap" $a_bp($i) 0]
        }
    }
    for {set i 1} {$i<=[string length $seq]} {incr i} {
        lappend l_basepairs $a_bp($i)
    }
# > test
#          1         2         3
# 123456789 123456789 123456789
# AA-AA-AAnAANNNNUU-UUnUU-UU
# (((((((((((....)))))))))))
# (( (( ((.((....)) )).)) ))
# ((((((.((....)))).))))
# 
#          1         2         3
# 123456789 123456789 123456789
# AAAAAAnAANNNNUUUUnUUUU
# ((((((.(.....)))).))))
# ((((((.((....)))).))))

# #     for {set i 0} {$i<[string length $seq]} {incr i} {
# #         puts -nonewline [format "%4s" [string index $seq $i]]
# #     }
# #     puts " "
# #     for {set i 0} {$i<[string length $seq]} {incr i} {
# #         puts -nonewline [format "%4d" [lindex $l_basepairs $i]]
# #     }
# #     puts " "
# #     for {set i 0} {$i<[string length $seq]} {incr i} {
# #         puts -nonewline [format "%3.1f" [lindex $l_prob $i]]
# #     }
# #     puts " "
#     for {set i 1} {$i<=[string length $seq]} {incr i +10} {
#         puts -nonewline [format "        %2d" [expr {1+$i/10}]]
#     }
#     puts " "
#     for {set i 1} {$i<=[string length $seq]} {incr i +10} {
#         puts -nonewline "123456789 "
#     }
#     puts " "
#     puts $seq
#     for {set i 1} {$i<=[string length $seq]} {incr i} {
#         if {$a_bp($i)==0} {
#             puts -nonewline "."
#         } elseif {$a_bp($i)>$i} {
#             puts -nonewline "("
#         } elseif {$a_bp($i)<$i} {
#             puts -nonewline ")"
#         } else {
#             puts "Dat jeht nicht"; exit
#         }
#     }
#     puts " "
    return [list $seq $l_basepairs $l_prob]
}
# remGaps
 
###   DrawStructFrontend
#
#
#
proc DrawStructFrontend {bpp seqnt header projname} {
    ###############################################
    global opts

    upvar $bpp bp_prob

    ### setup the structure for DrawStruct
    #
    set l_basepairs {}
    for {set i 1} {$i<=[string length $seqnt]} {incr i} {
        lappend l_basepairs $bp_prob($i,partner)
        lappend l_prob $bp_prob($i,prob)
    	#p::debug "$i: $bp_prob($i,partner) $bp_prob($i,prob)"
    }

    if {$opts(displ,remove_gaps)==1} {
        set result      [remGaps $seqnt $l_basepairs $l_prob]
        set seqnt       [lindex $result 0]
		set l_basepairs [lindex $result 1]
		set l_prob      [lindex $result 2]
    }

    set optimal     1
    #set projname    $cs_proj::proj(name)
    set win_pos    "+115+115"
    set DrawStruct::tmp_dir $opts(dir,work)
    DrawStruct::DrawStructure $seqnt $l_basepairs $optimal "$win_pos" \
    	-prob $l_prob -projname "$projname" -header $header

}
# DrawStructFrontend




###   CirclesFrontend
#
# setup variables prior to call of Circles
# Arguments:
#   bpp          : upvared basepair array
#                  index1: 1-alignmentlength, index2: "partner"|"prob"
#   seqid        ; sequence id
#   seqnt        : sequence nucleotides
#   match_method : "i"=imatch, "b"=bmatch,"c"=consensus-circles
#   header       : canvas header
#
#
proc CirclesFrontend {bpp  seqid  seqnt  match_method  header_l1 header_l2} {
    ######################################################################
    global opts
    upvar $bpp bp_prob


    if {$match_method!="i" && $match_method!="b" && $match_method!="c"} {
        p::error "Invalid match_method \"$match_method\""
        return
    }



    # convert to arrays used in cs2
    #
    if {$match_method=="c"} {
        for {set i 1} {$i<=[string length $seqnt]} {incr i} {
            set prob($i) $bp_prob($i,prob)
            set pair($i) $bp_prob($i,partner)
        }
    }


    ### remove rematched bps (neg probs) in imatch
    #   and construct the formerly (cs2) used array pair and prob
    #
    #
    if {$match_method=="i"} {
        for {set i 1} {$i<=[string length $seqnt]} {incr i} {

            # FIXME: formerly precircles checked 0 > prob > $opts(td,threshold)

            set prob($i) $bp_prob($i,prob)
            set pair($i) $bp_prob($i,partner)

            if {$prob($i) < 0} {
                # FIXME: use ?
                # set origpair($i) pair($i)
                # set origprob($i) [expr {abs(prob($i))}]
                # set rematched($i) "*"
    			# remove this and the above lines if it's message never appears
    			p::error "This shouldn't happen: neg probs should be ruled out in core"
                set pair($i) 0
                set prob($i) 0.0
            } else {
                # FIXME: use ?
                # set origpair($i) pair($i)
                # set origprob($i) prob($i)
                # set rematched($i) " "
            }
        }
    }


    ### remove rematched bps (neg probs) in  bmatch
    #   beware: bmatch returns a list of triple basepairs
    #
    #
    if {$match_method=="b"} {


        for {set i 1} {$i<=[string length $seqnt]} {incr i} {

            # FIXME: formerly precircles checked 0 > prob > $opts(td,treshold)

            set problist $bp_prob($i,prob)
            set pairlist $bp_prob($i,partner)

            # sanity check
            # we need basetriples -> 2 pairs and probs per nucleotide
            if {([llength $problist] != 2) || ([llength $pairlist] !=2 )} {
                p::error "Number of elements for problists and pairlists seems sick"
                return
            }

            # construct the formerly (cs2) used array pair and prob
            # whereby the indices range from 1 + 2*aln_len
            # in cs2: triple basepairs were appended to l_optsruct as
            # "second" pairlist from N+1...2N
            #
            #
            set idx  $i
            set prob($idx) [lindex $problist 0]
            set pair($idx) [lindex $pairlist 0]
            # remove rematched bps (neg probs) in bmatch
            if {$prob($idx) < 0} {
                set pair($idx) 0
                set prob($idx) 0.0
                # FIXME: use ?
                # set origpair($idx) pair($idx)
                # set origprob($idx) [expr {abs(prob($idx))}]
                # set rematched($idx) "*"
            } else {
                # FIXME: use ?
                # set origpair($idx) pair($idx)
                # set origprob($idx) prob($idx)
                # set rematched($idx) " "
            }
            #
            #####
            #
            set idx [expr {$i+[string length $seqnt]}]
            set prob($idx) [lindex $problist 1]
            set pair($idx) [lindex $pairlist 1]
            if {$prob($idx) < 0} {
                set pair($idx) 0
                set prob($idx) 0.0
                # FIXME: use ?
                # set origpair($idx) pair($idx)
                # set origprob($idx) [expr {abs(prob($idx))}]
                # set rematched($idx) "*"
            } else {
                # FIXME: use ?
                # set origpair($idx) pair($idx)
                # set origprob($idx) prob($idx)
                # set rematched($idx) " "
            }

        }
    }



##   FIXME: recode mapping
##
##    set file [format "%s/%s.map" $CS_DATDIR $seqname($SelSeq)]
##    if {[file exists $file]==1} {
##        set j 1;                      # position j in seq_w/o_gaps is position $nt($j) in aligned sequence
##        for {set i 1} {$i<=$alignlength} {incr i} {
##           if {$seq($i)!="-"} {
##              set nt($j) $i
##         incr j
##           }
##        }
##        vputs "Reading mapping data from $file ..." 5
##        set doMapping 1
##        set Map_ID [open $file r]
##        while {[gets $Map_ID zeile]>=0 && [string match ".." $zeile]!=1 } {vputs "1: $zeile" 3};    # Kommentare bis zum ".." ueberlesen
##        set line 0
##        while {[gets $Map_ID zeile]>=0} {
##           if {[string index $zeile 0]=="!"} {;                    # Zeile ueberlesen, die mit "!" beginnt
##              vputs "2: Kommentar: $zeile" 5
##           } else {
##              incr line
##              vputs "3: color,offset: $zeile" 5
##              scan $zeile "%s %f" mapData($line,color) mapData($line,offset);
##                vputs "$mapData($line,color) $mapData($line,offset)" 5;    # "color" #offset
##              gets $Map_ID zeile
##              while {[regsub -all "  " $zeile " " zeile]!=0} {};                # doppelte Leerzeichen entfernen
##              if {[string index $zeile 0]==" "} {set zeile [string range $zeile 1 end]};    # eventuelles 1. Leerzeichen entfernen
##              vputs "3: nts: $zeile" 5
##              set list [split $zeile]
##                foreach i $list {set mapData($line,$nt($i)) 0}
##           }
##        }
##        close $Map_ID
##    } else {
##        set doMapping 0
##    }

    global selection
    set seq_no $selection(seq_no)
    if {$opts(mapping,use_in_structaln) && \
        [cs_proj::mapinfo::exists $seq_no]} {
        set doMapping 1
        set line 2
        set mapData(1,color) blue;  # NotPaired
        set mapData(2,color) red;   # Paired
        set mapData(1,offset) 0.
        set mapData(2,offset) 1.
        for {set i 1} {$i<=[string length $seqnt]} {incr i} {
            if {[cs_proj::mapinfo::is_not_paired $seq_no $i]==1} {set mapData(1,$i) 0}
            if {[cs_proj::mapinfo::is_paired     $seq_no $i]==1} {set mapData(2,$i) 0}
        }
    } else {
        set doMapping 0
    }

    # FIXME hardcoded mapfile"
    set mapfile "$cs_proj::proj(name)"


    # FIXME hardcoded win_pos
    set win_pos "+413+108"

    if {$doMapping==1} {
        if {$match_method=="i"} {
            Circles::Circles $mapfile $header_l1 $header_l2 \
    			$seqid $seqnt pair prob $win_pos $doMapping $line mapData 1
        } elseif {$match_method=="b"} {
            Circles::Circles $mapfile $header_l1 $header_l2 \
    			$seqid $seqnt pair prob $win_pos $doMapping $line mapData 2
        } else {
            Circles::Circles $mapfile $header_l1 $header_l2 \
    			$seqid $seqnt pair prob $win_pos $doMapping $line mapData
        }
    } else {
        if {$match_method=="i"} {
            Circles::Circles $mapfile $header_l1 $header_l2 \
    			$seqid $seqnt pair prob $win_pos 0 0 0 1
        } elseif {$match_method=="b"} {
            Circles::Circles $mapfile $header_l1 $header_l2 \
    			$seqid $seqnt pair prob $win_pos 0 0 0 2
        } else {
            Circles::Circles $mapfile $header_l1 $header_l2 \
    			$seqid $seqnt pair prob $win_pos
        }
    }

    # FIXME: recode BmatchCircles
    #    if {$imatch==2} {
    #        GenBmatchBpWin origpair origprob rematched seq
    #    }
    SetCursor normal

}
#  CirclesFrontend




###   StructAlnFrontend
#
#
# tertiary: bool if structure may contain 3d basepairs
#
proc StructAlnFrontend {bpp_arrayname {tertiary 0}} {
######################
    global opts
    global seq
    global mic_is_computed_and_valid ;# r

    upvar $bpp_arrayname bpp

    if {$tertiary} {
    	set msg "The structure could contain tertiary interactions..."
    	append msg "These may corrupt some structure representations (e.g. dotbracket)"
    	set answer [tk_messageBox -parent . \
    					-title "Warning" -type ok \
    					-icon warning -message "$msg"]
    }

    SetCursor busy

    set csseq [Get_ConsSeq]

    set StructAln::work_dir $opts(dir,work)
    set StructAln::project_name ${cs_proj::proj(name)}

    StructAln::Run seq $csseq bpp \
        $opts(mic,unbiased) $opts(mic,bit) $cs_proj::proj(name) \
        $opts(structaln,show_seq_stat) $opts(structaln,show_struct_stat) $opts(structaln,show_patt_stat)

    SetCursor normal
}
# StructAlnFrontend




###   DotbracketToAlNum
#
#
proc DotbracketToAlNum {db_str} {
    ##########################

    set len [string length $db_str]
        for {set n 1} {$n<=$len} {incr n} {

            set c [string index $db_str [expr {$n-1}]]
            if {$c=="("} {
                lappend ob $n

            } elseif {$c==")"} {
                set struct($n) [lindex $ob end]
                set struct([lindex $ob end])  $n
                set ob [lreplace $ob end end]

            } else {
                set struct($n) "-1"
            }
        }

        set hc 97; # a
        for {set n 1} {$n<=$len} {incr n} {
            if {$hc==123} { set hc 65 };	# A follows z
            if {$hc== 91} { set hc 97 };	# a follows Z

            if {$struct($n)==-1} {
                set str_a($n) "."
                continue
            }
            if {$struct($n)<$n} {
                continue
            }

            set str_a($n)          [format "%c" $hc]
            set str_a($struct($n)) [format "%c" $hc]

            if {[HelixEnd_NoPk $n $len struct]} {
                incr hc
            }
        }

        set struct_str ""
        for {set i 1} {$i<=$len} {incr i} {
            append struct_str "$str_a($i)"
        }

        return $struct_str
}
# DotbracketToAlNum



###   HelixEnd_NoPk
#
# return 1 if nt at <nt_idx> or it's partner mark the end of helix
# struct_array should contain all partner indices and -1 for nonpaired
#
# !!!   NOTE   !!!
# Use this function only if you know it doesn't contain a pseudoknot !
# A correct version (you need minor info than just the dotbracket notation)
# is StructAln::HelixEnd
#
proc HelixEnd_NoPk {nt_idx aln_len struct_arrayname} {
    ###############################################
    upvar $struct_arrayname struct

    # helix end on either side?
    set next [expr {$nt_idx+1}]
    if {$next<=$aln_len} {
        if {$struct($next)==-1} {
            return 1
        }
    }
    set next [expr {$struct($nt_idx)-1}]
    if {$next<=$aln_len && $next>0} {
        if {$struct($next)==-1} {
            return 1
        }
    }
    return 0

}
# HelixEnd_NoPk



###   StructStringToBpList
#
#
proc StructStringToBpList {str} {
    ###########################

    #****f* StructStringToBpList
    # NAME
    #  StructStringToBpList
    # SYNOPSIS
    #
    # FUNCTION
    #
    # INPUTS:
    #  -str: either dotbracket notation or its alphanumeric counterpart
    #
    # RESULT
    #  Returns list where each element represents a corresponding basepair
    #  index for this element's index
    #  or <=0 if unpaired
    #  Returns empty string on error
    # EXAMPLE
    #  StructStringToBpList (((.....)))...
    #  or
    #  StructStringToBpList aaa..b..aaa..b
    # NOTES
    # BUGS
    #
    # SEE ALSO
    #
    #********

    set bracketstack {}
    array set alnumstack {}
    array set bparray {}
    array set locked {}
    set len [string length $str]
    set isdb [regexp -- {^[\(\.\)]*$} $str]

    for {set p 0} {$p<$len} {incr p} {
    	set char [string index $str $p]

    	if {[IsGap $char]} {
    		set bparray($p) -1
    		continue
    	}

    	# dotbracket makes things easy
    	if {$isdb} {
    		# push or pop
    		if {$char=="("} {
    			# push
    			lappend bracketstack "$p"
    		} else {
    			# pop
    			set partneridx [lindex $bracketstack end]
    			set bracketstack [lreplace $bracketstack end end]
    			set bparray($p) $partneridx
    			set bparray($partneridx) $p
    		}
    	# alphanumeric representation may store pseudoknots as well
    	} else {

    		if {[info exists alnumstack($char)]} {
    			if { ! [info exists locked($char)]} {
    				# push
    				lappend alnumstack($char) $p
    				#puts "  pushin pos=$p char=$char (now alnumstack($char)=$alnumstack($char))"
    			} else {
    				# pop
    				#puts "  poppin"
    				set partneridx [lindex $alnumstack($char) end]
    				set alnumstack($char) [lreplace $alnumstack($char) end end]
    				set bparray($p) $partneridx
    				set bparray($partneridx) $p
    			}
    		} else {
    			# push onto virgin stack
    			lappend alnumstack($char) $p
    			#puts "  pushin pos=$p char=$char (now alnumstack($char)=$alnumstack($char))"
    		}

    		# lock if gap or new helixchar behind this one
    		if {$p<[expr {$len-1}]} {
    			if {$char!=[string index $str [expr {$p+1}]]} {
    				set locked($char) 1
    			}
    		}

    	}
    }

    # construct return list
    set bplist {}
    for {set i 0} {$i<$len} {incr i} {
    	lappend bplist [expr {$bparray($i)+1}]
    }

    # check
    if {[llength $bplist]!=$len} {
    	p::error "internal error: created \#bplist != input string length"
    	return {}
    } else {
    	return $bplist
    }
}
# StructStringToBpList



###   StructTransform
#
# Convert a structure to a string representation
#
# OUT:
#  structure representation as string of requested format
#  or empty string on error
# IN:
#  len:           length of sequence/structure
#  bplist:        list of basepair indices or <=0 if unpaired
#                 array must have  indices:
#                 <i>,prob    basepair probability (where 1<=i<=len)
#                 <i>,partner basepair partner (where 1<=i<=len)
#  optional repr: alnum (default) or dotbracket
#
# FIXME: rewrite and use a bplist as arg instead to make things easier (don't need prob)
#
proc StructTransform {bplist len {outformat "alnum"}} {
    ##################################################

    if {$outformat!="alnum" && $outformat!="dotbracket" && $outformat!="stockholm"} {
    	p:error "wrong arg \"$outformat\""
    	return ""
    }

    set helixchar 97; # a

    set structrepr ""
    for {set i 1} {$i<=$len} {incr i} {
    	append structrepr "."
    }


    for {set i 1} {$i<=$len} {incr i} {

    	if {$outformat=="alnum"} {
    		if {$helixchar==123} {
    			set helixchar 65;# A follows z
    		}
    		if {$helixchar== 91} {
    			set helixchar 97;# a follows Z
    		}
    	}
    	set partner [lindex $bplist [expr {$i-1}]]
    	if {$partner<=0} {
    		continue
    	}
    	if {$partner<$i} {
    		continue
    	}

    	set ob1idx1 [expr {$i-1}]
    	set ob1idx2 [expr {$partner-1}]
    	if {$outformat=="dotbracket"} {
    		set char1 "("
    		set char2 ")"

    	} elseif {$outformat=="stockholm"} {
    		set char1 "<"
    		set char2 ">"
            
        } else {
    		set char1 [format "%c" $helixchar]
    		set char2 [format "%c" $helixchar]
    	}

    	set structrepr [string replace $structrepr $ob1idx1 $ob1idx1 $char1]
    	set structrepr [string replace $structrepr $ob1idx2 $ob1idx2 $char2]

    	if {[HelixEnd $i $len $bplist]} {
    			incr helixchar
    	}

    }
    return $structrepr
}
# StructTransform



###   HelixEnd
#
#
# Helix end _after_ given _paired_ nt
# -> process from 5 to 3'
#
# !!!   NOTE   !!!
# This function is (in contrast to HelixEnd_NoPk) PK-sensitiv , hahaha :)
#
proc HelixEnd {nt_idx len bplist} {
    ################################

    # already at end of sequence?
    if {[expr {$nt_idx+1}]>=$len} {
    	return 1
    }

    # partner of next nt is directly before partner of this one ?
    set nextpartner [lindex $bplist $nt_idx]
    set partneroffset [expr {[lindex $bplist [expr {$nt_idx-1}]]-1}]
    if {$nextpartner == $partneroffset} {
    	return 0
    } else {
    	return 1
    }
}
# HelixEnd



###   SaveStructFrontend
#
#
# bplist: list of basepair indices, or <=0 if unpaired
# format: csconsensus|connect
#
proc SaveStructFrontend {seqid seqnt bplist format {ctsubformat "zuker"}} {
    #####################################################################
    global opts
    global FILE_TYPES
    global seq;# rnaml only

    set filetypes $FILE_TYPES(struct)
    set defaultext [lindex [lsearch -glob -inline $FILE_TYPES(struct) "*${format}*"] end]

    regsub -all -- { } $seqid "" prefix
    regsub -- {\.$} $prefix "" prefix
    set defaultfile "${prefix}${defaultext}"

    if {[llength $bplist] != [string length $seqnt]} {
    	p::error "upps...basepair list and sequence length doesn't match"
    	p::fixme "bplist = $bplist"
    	p::fixme "seqnt = $seqnt"
    	return
    }

    set filename [tk_getSaveFile -title "Browse Structure Files" \
    				  -initialfile "$defaultfile" -initialdir $opts(dir,work) \
    				  -filetypes "$filetypes"]

    if {$filename==""} {
    	set fid stdout
    	set closefid 0
    } else {
    	set fid [open $filename w]
    	set closefid 1
    }


    if {$format=="connect"} {
    	SaveConnect $seqid $seqnt $bplist $ctsubformat $fid

    } elseif {$format=="cs-consensus"} {
    	SaveCsConsensus $seqid $seqnt $bplist $fid

    } elseif {$format=="rnaml"} {
    	rnaml::set_seq [array get seq]
    	rnaml::set_cs "consensus" [Get_ConsSeq] $bplist
    	rnaml::write $cs_proj::proj(name) $fid
    	
    } elseif {$format=="stockholm"} {
    	# FIXME: no pseudoknot support!
    	set sscons [StructTransform $bplist [string length $seqnt]  "stockholm"]
    	SaveStockholm [array get seq] $sscons  $fid
    	
    } else {
    	p::error "unknown requested structure format \"$format\""
    	# close and return
    }
    if {$closefid} {
    	close $fid
    }
}
# SaveStructFrontend




###   LoadConnect
#
#
# See also SaveConnect
#
proc LoadConnect {fname} {
    ####################

    set seqid "NOT_FOUND"
    set seqnt ""
    set blist {}

    set instruct 0
    set oldidx 0

    set fid [open $fname r]
    while {[gets $fid line]!=-1} {
    	# catch header with id
    	if {! $instruct} {
    		if {[regexp -- {\[initially[^\]]*\] *([^ ]*)} $line all seqid]} {
    			p::debug "got seqid $seqid"
    			continue
    		}
    	}

    	if {[scan $line "%d %s %d %d %d %d" realidx nt prev next partner realidx2]==6} {
    		# catch first nucleotide
    		if {! $instruct} {
    			if {$realidx==1} {
    				set instruct 1
    			} else {
    				continue
    			}
    		}
    		# new nt entry
    		if {$realidx!=[expr {$oldidx+1}]} {
    			p::warn "index mismatch ignoring line \"$line\""
    			continue
    		}
    		set oldidx $realidx
    		append seqnt $nt
    		lappend bplist $partner
    	} else {
    		p::warn "ignoring line \"$line\""
    		continue
    	}
    }
    close $fid
    return [list $seqid $seqnt $bplist]
}
# LoadConnect



###   SaveConnect
#
# save structure and sequence as connect file
# bplist must be a list (length=length seqnt) where each index
# lists the basepair partner the corresponding, or <=0 if unpaired
# See also LoadConnect
#
proc SaveConnect {seqid seqnt bplist {subformat "zuker"} {fid stdout}} {
    #################################################################

    if {$subformat!="zuker" && $subformat!="gcg" && $subformat!="rnaviz"} {
    	p::error "unknown arg \"$subformat\""
    	return
    }

    set seqlen [string length $seqnt]
    set title "internal error"

    # only first two lines change between the different formats
    #
    if {$subformat=="zuker"} {
    	set title [format "%5d   dG = -0.0  \[initially   -0.0\]    %-33s" \
    			 $seqlen $seqid]

    } elseif {$subformat=="gcg"} {
    	set title [format "ConStruct of:  \[initially -0.0\] %s Check: 0 from: 1 to:%4d" \
    			 $seqid $seqlen]
    	append title "\n"
    	append title [format "Length:%4d Energy:  -0.0 .." $seqlen ]

    } elseif {$subformat=="rnaviz"} {
    	set title [format "%5d   ENERGY = -0.0  \[initially   -0.0\]    %-43s" \
    			 $seqlen $seqid]
    }


    puts $fid $title

    for {set i 1} {$i<=$seqlen} {incr i} {
    	set previdx [expr {$i-1}]
    	set nextidx [expr {$i+1}]
    	set thisnt  [string toupper [string index $seqnt $previdx]]

    	set partner [lindex $bplist $previdx]
    	if {$partner<=0} {
    		set partner 0
    	} elseif {$partner>$seqlen} {
    		set msg "upps...corrupted structure:"
    		append msg "found basepair partner > sequence length!"
    		p::error  $msg
    		return
    	}

    	puts $fid [format "%5d %s %7d %4d %4d %4d" $i $thisnt \
    				   $previdx $nextidx $partner $i]
    }
    puts -nonewline $fid "\n"
}
# SaveConnect



###   LoadCsConsensus
#
# Loads consensus sequence  and stucture (as dotbracket and alphanumerical string)
# saved in a vienna-like format
# See also SaveCsConsensus
#
proc LoadCsConsensus {fname} {
    ################

    set seqid "NOT_FOUND"
    set seqnt "NOT_FOUND"
    set struct_db ""
    set struct_alnum ""

    set start 0
    set fid [open $fname r]
    while {[gets $fid line]!=-1} {
    	set line [string trim $line]
    	# get seqid==start
    	if {!$start && [regexp -- {^> (.*)$} $line all seqid]} {
    		set start 1
    	}
    	if {[gets $fid line]==-1} {
    		p::error "couldn't get seqnt"
    		return {}
    	}
    	set seqnt $line

    	if {[gets $fid line]!=-1} {
    		set struct_db $line
    	} else {
    		break
    	}

    	if {[gets $fid line]!=-1} {
    		set struct_alnum $line
    	} else {
    		break
    	}
    }
    close $fid

    # check
    # sequence length and structure have to be of same length
    if {[string length $seqnt]!=[string length $struct_db]} {
    	p::error "sequence length and structure length mismatch"
    	return {}
    }
    # we need at least dotbracket or alphanumeric structure representation
    if {$struct_alnum=="" && $struct_db==""} {
    	p::error "couldn't find a structure representation"
    	return {}
    }
    if {$struct_alnum!=""} {
    	set bplist [StructStringToBpList $struct_alnum]
    } else {
    	set bplist [StructStringToBpList $struct_db]
    }

    return [list $seqid $seqnt $bplist]
}
# LoadCsConsensus



###   SaveCsConsensus
#
# Saves consensus sequence  and stucture (as dotbracket and alphanumerical string)
# in a vienna-like format
# See also LoadCsConsensus
#
proc SaveCsConsensus {seqid seqnt bplist {fid stdout}} {
    ################

    set seqlen [string length $seqnt]
    set struct_alnum [StructTransform $bplist $seqlen]
    set struct_db [StructTransform $bplist $seqlen "dotbracket"]

    puts $fid "> $seqid"
    puts $fid "$seqnt"
    puts $fid "$struct_db"
    puts $fid "$struct_alnum"
}
# SaveCsConsensus



###   SaveStockholm
#
# Spec according to:
# http://www.cgr.ki.se/cgb/groups/sonnhammer/Stockholm.html
#
# IN:
#  seq: array list constructed with array get
#  ss_cons: consensus structure: [<>\.]*
#  fid: file descriptor, defaults to stdout
#
# OUT:
#  -1 on error
#  0  on success
#
proc SaveStockholm {seq_al sscons {fid stdout}} {
    #####################################################

    array set seq $seq_al
    set sscons_header "#=GC SS_cons"
    
    if { ! [regexp {^[\<\>\.]*$} $sscons]} {
    	p::error "sscons parse error...writing skipped"
    	return -1
    }
    if {[string length $sscons] != $seq(aln_len)} {
    	p::error "sscons length != alignment length...writing skipped"
    	return -1
    }
    
    set maxidlen [string length $sscons_header]
    for {set i 1} {$i<=$seq(n_seq)} {incr i} {
    	set this_len [string length $seq(id,$i)]
    	if {$this_len>$maxidlen} {
    		set maxidlen $this_len
    	}
    }
    if {$maxidlen>255} {
    	set maxidlen 255
    }

    
    puts $fid "# STOCKHOLM 1.0\n"
    for {set i 1} {$i<=$seq(n_seq)} {incr i} {
    	set this_nt  $seq(nt,$i)
        set this_id  $seq(id,$i)

    	set this_id [string range $seq(id,$i) 0 [expr {$maxidlen-1}]]
    	# "Sequence letters may include any characters except whitespace."
    	# "Gaps may be indicated by '.' or '-'."
    	# "Wrap-around alignments are allowed in principle,..."
    	# "Wrapped alignments are discouraged since they are much harder to parse."
    	set this_nt [string toupper $this_nt]
    	set this_nt [string map {- .} $this_nt]
    	
    	puts $fid [format "%-${maxidlen}s   %s" $this_id $this_nt]
    }
    puts $fid [format "%-${maxidlen}s   %s" $sscons_header $sscons]
    # "The "//" line indicates the end of the alignment."
    puts $fid "//\n"
    
    return 0
}
# SaveStockholm
