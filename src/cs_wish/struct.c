/******************************************************************************
* 
* struct.c - some useful macros and functions for structure handling
*
* Copyright (C) 2001-2004 Institute of Physical Biology,
*                    Heinrich-Heine Universität Düsseldorf, Germany
*                    <construct@biophys.uni-duesseldorf.de>
*
* This file is part of ConStruct.
* 
* ConStruct is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
* 
* ConStruct is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with ConStruct; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
* 
*******************************************************************************/

/*
 *  CVS $Id: struct.c,v 1.13 2007-10-22 10:43:23 steger Exp $    
 */



#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "struct.h"



/***   RemoveSingleBps   ******************************************************
 *
 * Removes single basepairs from structure [0...alnlen-1]
 * Walk through all basepairs and
 *
 * if i is paired
 *      if i-1 is same brace or if i+1 is same brace
 *          ok
 *      else
 *          remove
 * FIXME: this is not correct: (((...).)) is left unchecked
 *
 *

void
RemoveSingleBps(pair_prob_t *structure, int aln_len)
{
    int  i;
    int this_brace, next_brace, prev_brace;
    int removed_one=1;   
    
    while (removed_one)
    {
        removed_one=0;

        * start *
        i=0;

        if (structure[i].pair!=-1 && structure[i].pair>i)
        {
            if (structure[i].pair > i)
                this_brace = -1;
            else
                this_brace = 1;   
            
            if     (structure[i+1].pair==-1)
                next_brace = this_brace * -1;
            else if (structure[i+1].pair > i+1)
                next_brace = -1;
            else
                next_brace = 1;
                
            if (this_brace!=next_brace)
            {
                VERBOSE_P("removing single basepair %d:%d\n", structure[i].pair,
                                             structure[structure[i].pair].pair);
                DEBUG_P("removing single basepair %d:%d\n", structure[i].pair,
                                           structure[structure[i].pair].pair);
                structure[structure[i].pair].pair=-1;
                structure[structure[i].pair].prob=-1.0;
                structure[i].pair=-1;
                structure[i].prob=-1.0;
                removed_one=1;
            }
        }

        for (i=1; i<aln_len-1; i++)
        {
            if (structure[i].pair==-1)
                continue;
            
            if (structure[i].pair > i)
                this_brace = -1;
            else
                this_brace = 1;
                
            if     (structure[i+1].pair==-1)
                next_brace = this_brace * -1;
            else if (structure[i+1].pair > i+1)
                next_brace = -1;
            else
                next_brace = 1;
                
            if     (structure[i-1].pair==-1)
                prev_brace = this_brace * -1;
            else if (structure[i-1].pair > i-1)
                prev_brace = -1;
            else
                prev_brace = 1;
        
            if (this_brace!=next_brace && this_brace!=prev_brace)
            {
                VERBOSE_P("removing single basepair %d:%d\n", structure[i].pair,
                                             structure[structure[i].pair].pair);
                DEBUG_P("removing single basepair %d:%d\n", structure[i].pair,
                                           structure[structure[i].pair].pair);
                structure[structure[i].pair].pair=-1;
                structure[structure[i].pair].prob=-1.0;
                structure[i].pair=-1;
                structure[i].prob=-1.0;
                removed_one=1;
            }
        }
        
        * end *
        i=aln_len-1;
        if (structure[i].pair!=-1)
        {
            if (structure[i].pair > i)
                this_brace = -1;
            else
                this_brace = 1;
                
            if     (structure[i-1].pair==-1)
                prev_brace = this_brace * -1;
            else if (structure[i-1].pair > i+1)
                prev_brace = -1;
            else
                prev_brace = 1;
            if (this_brace!=prev_brace)
            {
                VERBOSE_P("removing single basepair %d:%d\n", structure[i].pair,
                                             structure[structure[i].pair].pair);
                DEBUG_P("removing single basepair %d:%d\n", structure[i].pair,
                                           structure[structure[i].pair].pair);
                structure[structure[i].pair].pair=-1;
                structure[structure[i].pair].prob=-1.0;
                structure[i].pair=-1;
                structure[i].prob=-1.0;
                removed_one=1;
            }
        }
        
    }
}
***   RemoveSingleBps   ***/

void
RemoveSingleBps(pair_prob_t *structure, int aln_len)
{
    int   nt;
    int   nt_pair;
    float prob;

    for(nt=1; nt<=aln_len; nt++)
    {
        nt_pair = structure[nt].pair;
        if (nt_pair<0)
            nt_pair = 0;
        prob    = structure[nt].prob;
        /* printf("%3d:%3d => ", nt, nt_pair); */
        if (nt_pair > 0 &&     /* nt:nt_pair is a potential pair */
            nt_pair > nt) {     /* check only for upstream pair; the downstream part is handled simultaneously */
            if (nt > 1       && /* nt  might have a predecessor */
                nt < aln_len) { /* check for a neighboring, stacking basepair */
                if (!(nt_pair==structure[nt-1].pair-1 ||
                      nt_pair==structure[nt+1].pair+1)) { /* if no stack, then delete the base pair*/
                    structure[nt     ].prob = 0.0;
                    structure[nt_pair].prob = 0.0;
                    structure[nt     ].pair = -1;
                    structure[nt_pair].pair = -1;
                    /* printf("Deleted: %d:%d\n", nt, nt_pair); */
                }
            } else {
                if (nt == 1) { /* there can't be a predecessor */
                    if (!(nt_pair==structure[nt+1].pair+1)) {
                        structure[nt     ].prob = 0.0;
                        structure[nt_pair].prob = 0.0;
                        structure[nt     ].pair = -1;
                        structure[nt_pair].pair = -1;
                        /* printf("Deleted: %d:%d\n", nt, nt_pair); */
                    }
                }
            }
        }
        /* printf("%3d:%3d\n", nt, nt_pair); */
    }
}
