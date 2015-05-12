##############################################################################
#
# paatern.tcl - create some fancy statistics useful for pattern searches
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
#  CVS $Id: pattern.tcl,v 1.2 2004-12-07 09:27:31 steger Exp $
#


# -- pattern
#
#
# namespace for creating patterns from an alignment
#


package provide pattern 0.1


###   namespace pattern
#
#
namespace eval pattern {
##################

    variable seq
    # seq ;# array with indices:
    #    aln_len (alignment length == length of all seqs)
    #    n_seq  (number of sequences)
    #    nt,<i> (where 1<=i<=n_seq   sequence nucleotides)
    #    id,<i> (where 1<=i<=n_seq   sequence ids)
    
    variable bpp
    # bpp ;# array with indices:
    #    <i>,prob    (where 1<=i<=aln_len   basepair probability)
    #    <i>,partner (where 1<=i<=aln_len   basepair partner)

    variable init_done

    variable pairs

    
    # ---   Init
    #
    # <FIXME:shortdescription>
    #
    # IN:
    #  alignment array as array list (described above)
    # OUT:
    # SIDEEFFECTS:
    # NOTES:
    #
    proc Init {seq_al bpp_al} {
        #######################
        variable seq
        variable bpp
        variable pairs
        variable init_done

        array set seq {}
        array set seq $seq_al
        array set bpp {}
        array set bpp $bpp_al


        set pairs [list "A,A" "A,C" "A,G" "A,U" "A,-" "A,N"]
        lappend pairs "C,C" "C,G" "C,U" "C,-" "C,N"
        lappend pairs "G,G" "G,U" "G,-" "G,N"
        lappend pairs "U,U" "U,-" "U,N"
        lappend pairs "-,-" "-,N" "N,N"
        

        set init_done 1
    }
    # Init


    
    ####   PatScanStats
    #
    # Create PatScan statistics and return as text
    #
    # IN:
    # OUT:
    # SIDEEFFECTS:
    # NOTES:
    #
    proc PatScanStats {} {
        ################
        variable init_done
        variable seq
        
        if { ! $init_done} {
            return "ERROR: Call init before calling me!"
        }
        
        set ret_txt ""

        set sumweight 0.0
        for {set n 1} {$n <= $seq(n_seq)} {incr n} {
            set sumweight [expr {$sumweight + $cs_proj::weight($n)}]
        }
        for {set n 1} {$n <= $seq(aln_len)} {incr n} {
            append ret_txt [PatScanPerColumn $n $sumweight]
        }
        return $ret_txt
        
        #    for {set s 1} {$s <= $seq(n_seq)} {incr s} {
        #        for {set n 0} {$n < $seq(aln_len)} {incr n} {
        #            set nt [string index $seq(nt,$s) $n]
        #        }
        #    }

    }
    # PatScanStats

    
    
    ###   PatScanHeader
    #
    # return a header line for description of the table produced by outPatScan
    #
    proc PatScanHeader {} {
        ######################

        variable pairs
        
        set ret_txt " #nt    #A   #C   #G   #U     #gap            profile       i:j   prob  "
        #foreach element [list "A,A" "A,C" "A,G" "A,U" "A,." "A,N" "C,C" "C,G" "C,U" "C,." "C,N" "G,G" "G,U" "G,." "G,N" "U,U" "U,." "U,N" ".,." ".,N" "N,N"]
        foreach element $pairs {
            append ret_txt " $element"
        }
        
        return $ret_txt
    }
    # PatScanHeader




    ###   PatScanPerColumn
    #
    # output a table describing the structure in terms useful for creation of a 
    # PatScan pattern
    #
    # Output is produced per nt position, that is per function call one line is written.
    #	pos			:= position of nt for which the line is written [1-alnlen]
    #	sumweight	:= sum of sequence weights; this is num of sequences if all weights are 1.
    #
    proc PatScanPerColumn {pos sumweight} {
        #########################################
        variable seq
        variable pairs
        variable bpp

        
        set ret_txt ""
        set pos_off [expr {$pos-1}]
        
        # define array with elements that are the possible nucleotides;
        #	i.e, all IUPAC symbols are summarized as "N"
        set types [list A C G U - N]
        foreach element $types {
            set nuc($element) 0
        }

        # count the nucleotides at the position pos in all sequences;
        #	IUPAC symbols are summarized as "N"
        for {set s 1} {$s <= $seq(n_seq)} {incr s} {
            set typ [string toupper [string index $seq(nt,$s) $pos_off]]
            if {[lsearch -exact $types $typ]==-1} {
                set typ N
            }
            incr nuc($typ)
        }

        # create a descritpive table header in case a new helix starts
        set bi $pos
        set bj $bpp($pos,partner)
        if {$bpp($pos,prob)>0.0} {
            set wahr 1
        } else {
            set wahr 0
        }
        
        if {$pos==1} {
            append ret_txt " [PatScanHeader]\n"
        } else {
            if {$bj>0} {
                set bim1 $pos_off
                set bjm1 $bpp($pos_off,partner)
                if {$bpp($pos_off,prob)>0.0} {
                    set wahrm1 1
                } else {
                    set wahrm1 0
                }
                if {$bjm1==0 || ($bjm1>0 && $bjm1!=[expr {$bj+1}])} {
                    append ret_txt "[PatScanHeader]\n"
                }
            }
        }

        # Now output the table line

        append ret_txt [format "%4d (" $pos];	# output of #nt
        set patscanStr " (";	# for profile
        set sumnucs 0
        set sumgaps 0
        foreach element [list A C G U -] {
            set nucweight $nuc($element)
            if {$element!="-"} {
                append ret_txt [format "%4.0f" $nucweight]
                set patscanStr [format "%s%3d" $patscanStr [expr {int(0.5+$nucweight/$sumweight*100.)}]]
                set sumnucs [expr {$sumnucs+$nucweight}]
                if {$element!="U"} {;	#	don't add a comma for the last #nt
                    append ret_txt ","
                    set patscanStr [format "%s," $patscanStr]
                }
            } else {
                set sumgaps [expr {$sumgaps+$nucweight}]
            }
        }
        append ret_txt [format ")   (%4.0f)" $sumgaps]
        append ret_txt [format "%18s)" $patscanStr]

        # try to determine a IUPAC symbol for the consensus nt
        if {$sumgaps< $seq(n_seq)} {
            set anz_o_gaps [expr { $seq(n_seq)-$sumgaps}]
            set outNuc 0
            foreach element [list A C G U] {
                if {$nuc($element)==$anz_o_gaps} {
                    append ret_txt " $element"; set outNuc 1
                }
            }
            if {$outNuc==0} {;	# it's not a single nt
                if       {[expr {$anz_o_gaps-$nuc(A)-$nuc(G)}]==0} {
                    append ret_txt " R"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(C)-$nuc(U)}]==0} {
                    append ret_txt " Y"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(A)-$nuc(C)}]==0} {
                    append ret_txt " M"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(A)-$nuc(U)}]==0} {
                    append ret_txt " W"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(C)-$nuc(G)}]==0} {
                    append ret_txt " S"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(G)-$nuc(U)}]==0} {
                    append ret_txt " K"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(A)-$nuc(C)-$nuc(G)}]==0} {
                    append ret_txt " v"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(A)-$nuc(C)-$nuc(U)}]==0} {
                    append ret_txt " h"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(A)-$nuc(G)-$nuc(U)}]==0} {
                    append ret_txt " d"; set outNuc 1
                } elseif {[expr {$anz_o_gaps-$nuc(C)-$nuc(G)-$nuc(U)}]==0} {
                    append ret_txt " b"; set outNuc 1
                }
            }
            if {$outNuc==0} {;	# it's neither a IUPAC nor a single nt
                append ret_txt "  "
            }
        }

        # output a char for the pairing direction, if the nt at pos does pair;
        #	$bi, $bj, $wahr are already determined for table header
        if {$bj>0} {
            append ret_txt [format "%4d%4d%6.3f" $bi $bj $wahr]
            if {$bi<$bj} {
                append ret_txt " ("
            } else {
                append ret_txt " )"
            }
        } else {
            #              "12341234123456 ("
            append ret_txt "                "
        }


        # output the number of basepairs types, if the nt at pos does pair;
        #	do this even if $bi>$bj.
        #	all IUPAC symbols are summarized as "N"
        if {$bj>0} {
            foreach elementI $types {
                foreach elementJ $types {
                    set bp($elementI,$elementJ) 0
                }
            }
            incr bi -1
            incr bj -1
            for {set i 1} {$i <= $seq(n_seq)} {incr i} {;	# count the
                                                            # base pairs at the position pos
                
                set typI [string toupper [string index $seq(nt,$i) $bi]]
                if {[lsearch -exact [list A C G U -] $typI]==-1} {
                    set typI N
                }
                set typJ [string toupper [string index $seq(nt,$i)) $bj]]
                if {[lsearch -exact [list A C G U -] $typJ]==-1} {
                    set typJ N
                }
                incr bp($typI,$typJ)
            }
            foreach element $pairs  {
                set elementI [string index $element 0]
                set elementJ [string index $element 2]
                if {$elementI!=$elementJ} {
                    set bp($elementI,$elementJ) [expr {$bp($elementI,$elementJ)+$bp($elementJ,$elementI)}]
                }
                append ret_txt [format %4d $bp($elementI,$elementJ)]
            }
        } else {
            #               "         1         2         3         4         5         6         7         8         9"
            #               "123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 "
            append ret_txt "                                                                                    "
        }


        # if any sequence contains a IUPAC nt at pos, write a remark
        if {[expr {$seq(n_seq)-$sumgaps-$sumnucs}]!=0} {
            append ret_txt " $seq(n_seq)!=$sumgaps+$sumnucs <== IUPAC nts\n"
        } else {
            append ret_txt " \n"
        }
        return $ret_txt
    }
    # PatScanPerColumn
}
# namespace pattern
