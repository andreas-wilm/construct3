##############################################################################
#
# initcore.tcl - routines to initialize cs_wish
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
#  CVS $Id: initcore.tcl,v 1.19 2004-08-11 08:11:43 wilm Exp $
#

###   InitCore   
#
#
# Only called once, after cs_wish has been loaded
#
proc InitCore {} {
#############
    global opts         ;# w
    global FILE_TYPES   ;# w
    global DEBUGGER_SKIP_ALL;# w

    set opts(verbose)    0
    set opts(debug)      0
    set opts(do_timing)  0


    # seqio read supported formats
    set FILE_TYPES(seqin) {
        {{Vienna}     {.vie}}
        {{FASTA}      {.fa} }
        {{Clustal}    {.aln}}
        {{EMBL}       {.embl}}
        {{MSF/GCG}    {.msf}}
        {{PHYLIP}     {.phy}}
        {{PIR/CODATA} {.pir}}
        {{ASN.1}      {.asn}}
        {{All Files}     *}
    }
    # seqio write supported formats: should be the same as in cs_seqio
    set FILE_TYPES(seqout) {
        {{Vienna}     {.vie}}
        {{FASTA}      {.fa} }
        {{Clustal}    {.aln}}
        {{MSF/GCG}    {.msf}}
        {{PHYLIP}     {.phy}}
        {{ASN.1}      {.asn}}
        {{All Files}     *}
    }

    set FILE_TYPES(struct) {
        {{cs-consensus}   {.cs}}
        {{connect}     {.ct}}
    	{{rnaml}       {.xml}}
        {{stockholm}     {.stk}}
        {{All Files}     *}
    }

    set FILE_TYPES(proj) {
        {{Project Files} {.proj}}
        {{All Files}        *}
    }
    set FILE_TYPES(ps) {
        {{Postscript}              {.ps}}
        {{Encapsulated Postscript} {.eps}}
        {{All Files}                  *}
    }
    set FILE_TYPES(txt) {
        {{Text Files}              {.txt}}
        {{All Files}                  *}
    }

    set ::tcl_prompt1 "puts -nonewline \"ConStruct>\""


    ###   packages
    #
    package require GetOpts 1.1
    package require p 0.1


    # set to 0 for online debugging
    set DEBUGGER_SKIP_ALL 1

}
# InitCore



