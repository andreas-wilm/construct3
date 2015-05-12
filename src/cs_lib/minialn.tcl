##############################################################################
#
# minialn.tcl - realign a a subrange of the alignment by external tools
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

##### MiniAl
#
#
proc MiniAl {} {
################
    puts stderr "MiniAl feature disabled"; return

    global w
    global minial

    set minial(addseq1) "cgcg"
    set minial(addseq2) "gcgcuccggcgc"
    set minial(addseq3) "cgcg"
    
    set w(minialn) .spin 
    catch {destroy $w(minialn)}
    toplevel       $w(minialn)
    wm title       $w(minialn) "Minialignment"
    wm iconname    $w(minialn) "Minialignment"
  # wm geometry    $w(minialn) +300+300

    set minial(fontnormal) "Helvetica 12"
    set minial(fontsmall)  "Helvetica 10"
   
    frame $w(minialn).buttons
    pack  $w(minialn).buttons -side bottom -fill x -pady 2m
	 
    set tmpTopDir "/tmp"
    set rand [string range [expr rand()] 2 end]
    set tmpDir [file join $tmpTopDir "cs-3.2.4_$rand"]
    if {[file exists $tmpTopDir]} {
        if {![file exists $tmpDir]} {
            if {[catch {file mkdir $tmpDir}]} {
                puts stderr "ERROR: Minial couldn't create directory $tmpDir"
                return 0
            }
        }
    } else {
        puts stderr "ERROR: Minial needs directory $tmpTopDir"
        return 0
    }
    set minial(file)     [file join $tmpDir neueoutDatei.vie]
    set minial(newfile)  [file join $tmpDir neueinDatei.vie]
    set stralDir         [file join $tmpDir resultDIR]
    if {[catch {file mkdir $stralDir}]} {
        puts stderr "ERROR: Minial couldn't create directory $stralDir"
        return 0
    }
    set minial(stralfile)   [file join $stralDir neueoutDatei.vie]
    set minial(mafftfile)   [file join $tmpDir mafftinDatei.ginsi]
    set minial(newalnfile)  [file join [file dirname $cs_proj::proj(aln_path)] [file rootname $cs_proj::proj(file)]_${rand}.vie]
    set minial(newprojfile) [file rootname $minial(newalnfile)].proj

    label  $w(minialn).msg0 -font $minial(fontnormal) -wraplength 5i -justify left -text "Realign a stem-loop structure"
    label  $w(minialn).msg1 -font $minial(fontsmall)  -wraplength 5i -justify left -text "(Please insert nucleotide numbers or use right mouse click on a consensus basepair.)"
    pack   $w(minialn).msg0     \
           $w(minialn).msg1    -side top  -padx 1m -pady 2m

    frame  $w(minialn).w1      -bd 2
    label  $w(minialn).w1.ln1  -text "1st BP 5': "
    pack   $w(minialn).w1.ln1  -side left
    entry  $w(minialn).w1.nuc1 -width  9 -relief sunken -textvariable minial(nuc1)
    pack   $w(minialn).w1.nuc1 -side right
    pack   $w(minialn).w1      -side top -padx 1m -pady 2m

    frame  $w(minialn).w2      -bd 2
    label  $w(minialn).w2.ln2  -text "2nd BP 5': "
    pack   $w(minialn).w2.ln2  -side left
    entry  $w(minialn).w2.nuc2 -width  9 -relief sunken -textvariable minial(nuc2)
    pack   $w(minialn).w2.nuc2 -side right
    pack   $w(minialn).w2      -side top -padx 1m -pady 2m

    frame  $w(minialn).w3      -bd 2
    label  $w(minialn).w3.ln3  -text "1st BP 3': "
    pack   $w(minialn).w3.ln3  -side left
    entry  $w(minialn).w3.nuc3 -width  9 -relief sunken -textvariable minial(nuc3)
    pack   $w(minialn).w3.nuc3 -side right
    pack   $w(minialn).w3      -side top -padx 1m -pady 2m

    frame  $w(minialn).w4      -bd 2
    label  $w(minialn).w4.ln4  -text "2nd BP 3': "
    pack   $w(minialn).w4.ln4  -side left
    entry  $w(minialn).w4.nuc4 -width  9 -relief sunken -textvariable minial(nuc4)
    pack   $w(minialn).w4.nuc4 -side right
    pack   $w(minialn).w4      -side top -padx 1m -pady 2m

#    entry  $w(minialn).ent  -width 30 -relief sunken -textvariable minial(newalnfile)
    label  $w(minialn).msg3 -font $minial(fontnormal) -text "Using program:"
    pack   $w(minialn).msg3 -side top -padx 1m -pady 1m
    if {[check4pgm stral]} {
        button $w(minialn).stral -text StrAl -command "DoStral"
    } else {
        button $w(minialn).stral -text StrAl -state disabled
    }
    if {[check4pgm mafft]} {
        button $w(minialn).mafft -text Mafft -command "DoMafft"
    } else {
        button $w(minialn).mafft -text Mafft -state disabled
    }
    if {[check4pgm stemloc]} {
        button $w(minialn).stemloc -text stemloc -command "DoStemloc"
    } else {
        button $w(minialn).stemloc -text stemloc -state disabled
    }
    
#    button $w(minialn).files   -text Files     -command "Mini_In_Out"
    button $w(minialn).clear   -text "Clear"   -command "ClearMinialn"
    pack   $w(minialn).msg3     \
           $w(minialn).stral    \
           $w(minialn).mafft    \
           $w(minialn).stemloc  \
           $w(minialn).clear    -side top -padx 1m -pady 1m
}
### MiniAl

#### check4pgm
#
# if program is in path and executable
#   return 1
# else
#   return 0
#
proc check4pgm {program} {
##############
    global env

    foreach dir [split $env(PATH) ":"] {
        if {$dir == ""} {set dir .}
        set file [file join $dir $program]
        if {[catch {file exists $file} result] == 0 && $result}  {
            if {[file executable $file]} {
                # puts "check4pgm: $file"
                return 1
            } else {
                # puts "check4pgm: $file; $result ???"
                return 0
            }
        }
    }
    # puts "check4pgm: $program not found"
    return 0
}
#### check4pgm


#### Mini_In_Out
#
#
#
proc  Mini_In_Out {} {
#################
    global w
    global minial

    catch {destroy $w(minialn)}
    toplevel    $w(minialn)
    wm title    $w(minialn) "Minialignment"
    wm iconname $w(minialn) "Minialignment"
  # wm geometry $w(minialn) +300+300

    frame $w(minialn).buttons
    pack  $w(minialn).buttons -side bottom -fill x -pady 2m

	 
    label  $w(minialn).filen -font $minial(fontnormal) -wraplength 2i -justify left -text "Save Minialignment in: " 
    entry  $w(minialn).file -width 30 -relief sunken -textvariable minial(file)
    button $w(minialn).buttonsu -text "MiniAl" -command "ShowMiniAln  "    
    pack   $w(minialn).filen   \
           $w(minialn).file    \
           $w(minialn).buttonsu  -side top -padx 1m -pady 1m
	 
    label  $w(minialn).filem -font $minial(fontnormal) -wraplength 4i -justify left -text "Insert Minialignment: " 
    entry  $w(minialn).film -width 30 -relief sunken -textvariable minial(newfile)			  
    button $w(minialn).buttonst -text "ReplaceAln" -command " ReadnewMiniAlnFile selba"
    label  $w(minialn).lab -font $minial(fontnormal) -wraplength 2i -justify left -text "...and save as: "
    entry  $w(minialn).ent -width 30 -relief sunken -textvariable minial(newalnfile)
    button $w(minialn).clear -text "Clear" -command "ClearMinialn"
    pack   $w(minialn).filem    \
           $w(minialn).film     \
           $w(minialn).buttonst \
           $w(minialn).lab      \
           $w(minialn).ent      \
           $w(minialn).clear    -side top  -padx 1m -pady 1m
}  
### Mini_In_Out


##### ShowMiniAln  
#
#
#
proc ShowMiniAln  {} {
################
    global seq 
    global minial
    global c
    global c3
	 
   
	# faerbt gewaehltes Minialignment gruen
    
    ClearMinialn

    set coords_1  [$c(al) bbox AlnNt_seq_1_nt_$minial(nuc1)]
    set coords_2  [$c(al) bbox AlnNt_seq_${cs_proj::nseq}_nt_$minial(nuc2)]

    set coords_3  [$c(al) bbox AlnNt_seq_1_nt_$minial(nuc3)]
    set coords_4  [$c(al) bbox AlnNt_seq_${cs_proj::nseq}_nt_$minial(nuc4)]

    set x1_coord_list [list [lindex $coords_1 0] [lindex $coords_2 0]]
    set x2_coord_list [list [lindex $coords_1 2] [lindex $coords_2 2]]
    set y1_coord_list [list [lindex $coords_1 1] [lindex $coords_2 1]]
    set y2_coord_list [list [lindex $coords_1 3] [lindex $coords_2 3]]

    set x3_coord_list [list [lindex $coords_3 0] [lindex $coords_4 0]]
    set x4_coord_list [list [lindex $coords_3 2] [lindex $coords_4 2]]
    set y3_coord_list [list [lindex $coords_3 1] [lindex $coords_4 1]]
    set y4_coord_list [list [lindex $coords_3 3] [lindex $coords_4 3]]

    set x1 [lindex [lsort -real $x1_coord_list] 0]
    set x2 [lindex [lsort -real $x2_coord_list] 1]
    set y1 [lindex [lsort -real $y1_coord_list] 0]
    set y2 [lindex [lsort -real $y2_coord_list] 1]

    set x3 [lindex [lsort -real $x3_coord_list] 0]
    set x4 [lindex [lsort -real $x4_coord_list] 1]
    set y3 [lindex [lsort -real $y3_coord_list] 0]
    set y4 [lindex [lsort -real $y4_coord_list] 1]

    foreach ca [list $c(al) $c3(al)] {
        set box [$ca create rectangle $x1 $y1 $x2 $y2  -tags Mini1AlnSelBox \
                                  -fill green -outline "" ]
        $ca lower $box
        $ca lower REMatchBox					 
    										 
        set box [$ca create rectangle $x3 $y3 $x4 $y4  -tags Mini2AlnSelBox \
                                  -fill green -outline "" ]					 
        $ca lower $box
        $ca lower REMatchBox
    }
   
	
	
	#speichert Minialignment eingeschlossen von kuenstlichen Helices (ccc...ggg) und 
	#mit "uncg"???? zwischen den beiden Abschnitten
    set fid [open $minial(file) w]
        for {set s 1} {$s<=$seq(n_seq)} {incr s} {
            puts $fid "> $cs_proj::seq_id($s)"
            puts -nonewline $fid $minial(addseq1)
            for {set n $minial(nuc1)} {$n<=$minial(nuc2)} {incr n} {
                puts -nonewline $fid [string index $seq(nt,$s) [expr {$n-1}]]
            }		
            puts -nonewline $fid $minial(addseq2)
            for {set n $minial(nuc3)} {$n<=$minial(nuc4)} {incr n} {
                puts -nonewline $fid [string index $seq(nt,$s) [expr {$n-1}]]
            }
            puts -nonewline $fid  $minial(addseq3)
            puts $fid "\n"
        }
    close $fid

    # Minialignment Laengenberechnung:
    # set minial(len) [expr {16+$minial(nuc2)-$minial(nuc1)+$minial(nuc4)-$minial(nuc3)}]
}
## ShowMiniAln
 
##### DoStral
#
#
proc DoStral {} {
############
    global minial

    ShowMiniAln
    set pwd [pwd]
    cd [file dirname $minial(file)]
    exec stral $minial(file)
    cd $pwd
    ReadnewMiniAlnFile stral
} 

   
##### DoMafft
#
#
proc DoMafft {} {
############
    global minial

    ShowMiniAln
    set doit "exec ginsi $minial(file) > $minial(mafftfile)"
	 
    if {[catch $doit fehler]!=0} {
        if {file size $minial(mafftfile) != 0 } {
            puts "SUCCESS in DoMafft:\n$fehler\n"
        } else {
            puts "ERROR in DoMafft:\n$fehler\n"
        }
    } else {
        puts "ERROR in DoMafft:\n$fehler\n"
    }
    ReadnewMiniAlnFile mafft
}
### DoMafft
 
 
 
  
##### DoStemloc
#
#
proc DoStemloc {} {
##############
    puts "Stemloc: not implemented!"
}
#### DoStemloc
 
 
#### ClearMinialn
#
# 
proc ClearMinialn {} {
#################
    global c
    global c3

    $c(al)  delete Mini1AlnSelBox
    $c(al)  delete Mini2AlnSelBox
    $c3(al) delete Mini1AlnSelBox
    $c3(al) delete Mini2AlnSelBox
}	  
#### ClearMinialn


#### ReadnewMiniAlnFile
#
#liest Minialignmentfile mit zwei Teilabschnitten (fuer Helix), welche von kuenstlichen Helices umgeben sind.
#
#AENDERN:SEQ'IDs , ggg ,ccc und uncg entfernen!!
#
#
proc ReadnewMiniAlnFile {meth} {
#######################
    global minial
    global seq   

    if {$meth == "stral"} {
        set newf $minial(stralfile)
    } elseif {$meth=="mafft"} {
        set newf $minial(mafftfile)
    } elseif {$meth=="stemloc"} {
        set newf $minial(stemlocfile)
    } elseif {$meth=="selba"} {
        set newf $minial(newfile)
    } else {
        puts "Please give method!"
        ClearMinialn
        return
    }
    set d 1
    set fid [open $newf r]
        for {set s 0} {$s <= [expr 3*$seq(n_seq)]} {incr s} {
            gets $fid line
            if  {[regexp {^(>)} $line]} {
            } elseif {[regexp {(?n)^\s*$} $line]} {
            } else {
                set i [string first $minial(addseq2) $line]
                if {$i!=0   &&  \
                    [string equal $minial(addseq1) [string range $line 0 [expr [string length $minial(addseq1)]-1]]]    && \
                    [string equal $minial(addseq3) [string range $line end-[expr [string length $minial(addseq3)]-1] end]]} {
                    set seq(mta,$d) [string range $line [string length $minial(addseq1)] [expr $i-1]]
                    set length [string length $line]
                    set seq(mtb,$d) [string range $line [expr $i+[string length $minial(addseq2)]]  \
                                                        [expr $length-1-[string length $minial(addseq3)]]]

                    incr d
                } else {
                    puts "MiniAl: New alignment wasn't consistent; try something else ..."
                    close $fid
                    ClearMinialn
                    return
                }
            }		  
        }
    close $fid
    InsertnewMinialignment   
}
#### ReadnewMiniAlnFile


#### InsertnewMinialignment
#
# inserts selected Minialignment and writes it into newalnfile,
# invokes procs for creating the new Dotplot
#
proc InsertnewMinialignment {} {
###########################
    global seq
    global minial

    for {set s 1} {$s<=$seq(n_seq)} {incr s} {
        # puts "> $seq(id,$s)"
        # puts "$seq(nt,$s)"
        # puts "$seq(mta,$s)"
        # puts "$seq(mtb,$s)"
        set dummy    [string range $seq(nt,$s) 0                      [expr {$minial(nuc1)-2}]]
        append dummy $seq(mta,$s)
        append dummy [string range $seq(nt,$s) [expr {$minial(nuc2)}] [expr {$minial(nuc3)-2}]]
        append dummy $seq(mtb,$s)
        append dummy [string range $seq(nt,$s) [expr {$minial(nuc4)}] end]
        set seq(nt,$s) $dummy
        # puts "$seq(nt,$s)\n"
    }
		
    set fid [open $minial(newalnfile) w]
        for {set s 1} {$s<=$seq(n_seq)} {incr s} {
            puts $fid "> $seq(id,$s)"
            puts $fid "$seq(nt,$s)\n"
        }
    close $fid
    set seq(aln_len) [string length $seq(nt,1)]

    WriteProj_tcl
    LoadProject $minial(newprojfile)
    # puts "LoadProj done"
    ClearMinialn
}
#### InsertnewMinialignment



#### WriteProj_tcl   
#
# Write a project file (using cs_proj) in the same directory as the sequence
# file with the same name extension $FILE_EXT(proj)
#
#
proc WriteProj_tcl {} {
##################
    global minial

    set oldprojfile [file join [file dirname $cs_proj::proj(aln_path)] ${cs_proj::proj(name)}.proj]
    file copy -- $oldprojfile $minial(newprojfile)
    set doit    {eval exec sed -ie \"s/^Alignment: }
    append doit [file tail $cs_proj::proj(aln_path)]
    append doit {/Alignment: }
    append doit [file tail $minial(newalnfile)]
    append doit {/\" }
    append doit $minial(newprojfile)
    puts $doit
    if {[catch $doit fehler]!=0} {puts "ERROR during WriteProj_tcl:\n$fehler\n"; exit}
}
#### WriteProj_tcl
