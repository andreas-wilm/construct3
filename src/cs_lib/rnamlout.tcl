##############################################################################
#
# rnamlout.tcl - rnaml writing routines
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
#  CVS $Id: rnamlout.tcl,v 1.3 2004-08-11 08:14:22 wilm Exp $
#


package provide rnaml 0.1


namespace eval rnaml {

    variable seq
    variable cs
    
    variable ilevel
    variable fid
    variable cs_lit_ref "construct_lueck_2001"


    
    ###   set_seq
    #
    # seq_arraylist: aligned sequences as array list:
    #                with elements n_seq and id,nt for each sequence
    #
    proc set_seq {seq_arraylist} {
    	########################
    	variable seq
    	
    	array unset seq
    	array set seq $seq_arraylist
    }
    # set_seq



    ###   set_cs
    #
    #
    proc set_cs {id nt {bplist ""}} {
    	##########
    	variable cs
    	
    	set cs(id) $id
    	set cs(nt) $nt
    	if {[llength $bplist]} {
    		if {[llength $bplist]==[string length $nt]} {
    			set cs(bp) $bplist
    		} else {
    			set msg "bplist length and sequence length mismatch"
    			append msg " ignoring bplist"
    			p::error $msg
    		}
    	}
    }
    # set_cs


    
    ###   pindent
    #
    proc pindent {str} {
    	variable ilevel
    	variable fid

    	set out ""
    	for {set i 1} {$i<=$ilevel} {incr i} {
    		append out "\t"
    	}
    	append out $str

    	puts $fid $out
    }
    # pindent

    

    ###   write
    #
    #
    proc write {projname {thisfid stdout}} {
    	####################################
    	variable ilevel
    	variable fid
    	variable seq
    	variable cs
    	

    	if {![info exists cs] || ![info exists seq]} {
    		p::error "use set_cs and set_seq before calling me"
    		return -1
    	}
    	
    	set projname [string map {\  _} "$projname"]
    	set fid $thisfid
    	
    	# preamble
    	#
    	set ilevel 0
    	set dtduri "http://www-lbit.iro.umontreal.ca/rnaml/current/rnaml.dtd"
    	pindent "<?xml version=\"1.0\"?>"
    	pindent "<!DOCTYPE rnaml SYSTEM \"${dtduri}\">"
    	pindent "<rnaml comment=\"Created with ConStruct\" version=\"1.1\">"

    	set ilevel 1
    	pindent "<molecule-class>"

    	incr ilevel
    	pindent "<identity>"
    	incr ilevel
    	pindent "<name>Consensus Sequence for ${projname}</name>"
    	incr ilevel -1
    	pindent "</identity>"

    	
    	# all molecules
    	#
    	for {set i 1} {$i<=$seq(n_seq)} {incr i} {
    		set id $seq(id,$i)
    		set idnospace [string map {\  _} "$id"]
    		set name $seq(id,$i)
    		set nalseq [string tolower [Degap $seq(nt,$i)]]
    		set seqlen [string length $nalseq]

    		pindent "<molecule id=\"${idnospace}\">"
    		incr ilevel
    		pindent "<identity><name>${name}</name></identity>"
    		pindent "<sequence length=\"$seqlen\"><seq-data>${nalseq}</seq-data></sequence>"
    		
    		# FIXME support writing of unaligned matrices here as structuremodel

    		incr ilevel -1
    		pindent "</molecule>"
    	}


    	# alignment
    	#
    	set now [clock seconds]
    	set day [clock format $now -format "%d"]
    	set month [clock format $now -format "%m"]
    	set year [clock format $now -format "%y"]
    	set hour [clock format $now -format "%H%M"]
    	set datestr "${year}${month}${day}_${hour}"
    	#
    	set alnid "${projname}" ;# see below
    	set analid "${projname}_${datestr}"
    	#
    	pindent "<alignment reference-ids=\"${analid}\" id=\"${alnid}\">"
    	incr ilevel
    	for {set i 1} {$i<=$seq(n_seq)} {incr i} {
    		set id $seq(id,$i)
    		set alnseq [string tolower $seq(nt,$i)]
    		pindent "<ali-sequence>"
    		incr ilevel
    		pindent "<molecule-id ref=\"${id}\"/>"
    		pindent "<seq-data>${alnseq}</seq-data>"
    		incr ilevel -1
    		pindent "</ali-sequence>"
    	}
    	incr ilevel -1
    	pindent "</alignment>"


    	if {[info exists cs]} {
    		pindent "<consensus-molecule>"

    		incr ilevel
    		pindent "<alignment-id ref=\"$alnid\"/>"
    		set aln_id "${projname}" ;# see above
    		set modelid "model"
    		pindent "<molecule id=\"${cs(id)}\">"
    		incr ilevel
    		pindent "<sequence><seq-data>${cs(nt)}</seq-data></sequence>"

    		pindent "<structure>"
    		incr ilevel
    		pindent "<model id=\"${modelid}\">"
    		incr ilevel

    		write_strannotation $cs(bp)
    	
    		incr ilevel -1
    		pindent "</model>"
    		incr ilevel -1
    		pindent "</structure>"
    		incr ilevel -1
    	
    		pindent "</molecule>"
    		incr ilevel -1
    		pindent "</consensus-molecule>"
    		incr ilevel -1
    	}
    	
    	pindent "</molecule-class>"


    	write_analysis $analid $year $month $day

    	write_cs_ref

    	puts $fid "</rnaml>"
    }
    # write


    
    ###   write_strannotation
    #
    #
    proc write_strannotation {bplist} {
    	################
    	variable ilevel
    	variable fid

    	pindent "<str-annotation>"
    	incr ilevel
    	
    	for {set i 0} {$i<[llength $bplist]} {incr i} {
    		set j [lindex $bplist $i]
    		if {$j>0} {
    			if {$i<$j} {
    				set p5 [expr {$i+1}]
    				set p3 $j
    			} else {
    				continue
    			}
    			
    			pindent "<base-pair>"
    			incr ilevel

    			pindent "<base-id-5p>"
    			incr ilevel
    			pindent "<base-id><position>${p5}</position></base-id>"
    			incr ilevel -1
    			pindent "</base-id-5p>"

    			pindent "<base-id-3p>"
    			incr ilevel
    			pindent "<base-id><position>${p3}</position></base-id>"
    			incr ilevel -1
    			pindent "</base-id-3p>"

    			incr ilevel -1
    			pindent "</base-pair>"
    		}
    	}

    	incr ilevel -1
    	pindent "</str-annotation>"
    }
    # write_strannotation


    
    ###   write_cs_ref
    #
    #
    proc write_cs_ref {} {
    	###################
    	variable ilevel
    	variable fid
    	variable cs_lit_ref
    	
    	pindent "<!-- Info about this program -->"
    	pindent "<reference id=\"${cs_lit_ref}\">"
    	# FIXME add author/person/first-name/last-name
    	incr ilevel
    	pindent "<author><person><first-name>Rupert</first-name><last-name>Lück</last-name></person></author>"
    	pindent "<author><person><first-name>Stefan</first-name><last-name>Gräf</last-name></person></author>"
    	pindent "<author><person><first-name>Gerhard</first-name><last-name>Steger</last-name></person></author>"
        pindent "<title>"
    	pindent "ConStruct: a tool for thermodynamic controlled prediction of conserved secondary structure."
    	pindent "</title>"
        pindent "<journal>Nucleic Acids Res.</journal>"
        pindent "<pubmed-id>10518612</pubmed-id>"
        pindent "<volume>27</volume>"
        pindent "<issue>21</issue>"
        pindent "<pages>4208-17</pages>"
    	incr ilevel -1
    	pindent "</reference>"
    }
    # write_cs_ref


    ###   write_analysis
    #
    proc write_analysis {analid year month day} {
    	#####################
    	variable ilevel
    	variable fid
    	variable cs_lit_ref

    	set progname "ConStruct"
    	set progvers "3.0"

    	pindent "<!-- Info about this analysis -->"
    	pindent "<analysis reference-ids=\"${cs_lit_ref}\" id=\"${analid}\">"
    	incr ilevel
    	pindent "<program>"
    	incr ilevel
    	pindent "<prog-name>${progname}</prog-name>"
    	pindent "<prog-version>${progvers}</prog-version>"
    	incr ilevel -1
    	pindent "</program>"

    	pindent "<date><day>${day}</day><month>${month}</month><year>${year}</year></date>"
    	incr ilevel -1
    	pindent "</analysis>"
    }
    # write_analysis

}
# namespace eval rnaml
