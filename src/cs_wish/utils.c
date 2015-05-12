/******************************************************************************
*
* utils.c - various routines
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
 *  CVS $Id: utils.c,v 1.31 2007-10-22 10:43:23 steger Exp $
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "public_datatypes.h"
#include "main.h"
#include "stopwatch.h"
#include "generic_utils.h"

#include "utils.h"




/***   TimerStart   ***********************************************************
 *
 * Start a timer if cs_opts.do_timing is 1 and return it
 * Return NULL otherwise
 */
Stopwatch_t *
TimerStart(void)
{
    Stopwatch_t *watch = NULL;
    if (cs_opts.do_timing)
    {
        watch = StopwatchCreate();
	    StopwatchZero(watch);
	    StopwatchStart(watch);
    }
    return watch;
}
/***   TimerStart   ***/



/***   TimerStop   ***********************************************************
 *
 */
void
TimerStop(const char *prefix, Stopwatch_t *watch)
{
    char buf[1024];
    if (cs_opts.do_timing && watch!=NULL)
    {
	    StopwatchStop(watch);
        StopwatchDisplay(buf, watch);
        fprintf(stdout, "%s: CPU time  %s", prefix, buf);
        StopwatchFree(watch);
    }
}
/***   TimerStop   ***/



/***   SetPair   ***********************************************************
 *
 * Sets a pair with the corresponding indices
 * All matrices can only be accessed by [x][y], where x>y
 *
 *
 */
void
SetPair(int nti, int ntj, pair_t *p)
{
    p->nti = MAX(nti, ntj);
    p->ntj = MIN(nti, ntj);
}
/***   SetPair   ***/



/***   CreateGapShift   *******************************************************
 *
 * Compute the gapshift from an aligned sequence and return it as int array
 * (zero-offset)
 *
 * Caller must free
 *
 */
int *
CreateGapShift(char *aln_seq)
{
    int   dealn_seq_len, aln_seq_len;
    char *dealn_seq;
    int  *gap_shift;
    int   introd_gaps;
    int   i, j;


    aln_seq_len   = strlen(aln_seq);
    dealn_seq = (char*) Scalloc(aln_seq_len+2, sizeof(char));
    Degap(aln_seq, dealn_seq);
    dealn_seq_len = strlen(dealn_seq);

    
    /* gap_shift = (int*) Scalloc(dealn_seq_len, sizeof(int)); */
    gap_shift = (int*) Scalloc(aln_seq_len, sizeof(int));

    j=0;
    introd_gaps=0;
    for (i=0; i<aln_seq_len; i++)
    {
        if ( IS_GAP(aln_seq[i]) )
        {
            introd_gaps++;
            /* new */
            gap_shift[i] = 0;
        }
        else
        {
            /* gap_shift[j] = introd_gaps;
            j++;
             */
            gap_shift[i] = introd_gaps;
        }
    }

    /*
    TMPDEBUG_P("gapshift for %s\n", aln_seq);
    for (i=0; i<aln_seq_len; i++)
        fprintf(stdout, "%3d ", gap_shift[i]);
    fprintf(stdout, "\n");
    */

    free(dealn_seq);

    return gap_shift;
}
/***   CreateGapShift   ***/



/***   GetBpFromAlnSeq   **************************************************
 *
 * nt3<nt5: aligned nt indices
 * return pointer to cooresponding bpair_t
 * or null if one of it's nts is a gap
 */
bpair_t *
GetBpFromAlnSeq(aln_t *aln, int seqidx, pair_t pair)
{
    int gs3, gs5;

    /* TMPDEBUG_P("seq %d: %d:%d=%c:%c\n", seqidx, pair.nti, pair.ntj, aln->seq[seqidx].nt[pair.nti], aln->seq[seqidx].nt[pair.ntj]);    */
    if (IS_GAP(aln->seq[seqidx].nt[pair.nti]) \
        || \
        IS_GAP(aln->seq[seqidx].nt[pair.ntj]))
        return NULL;

    gs5 = aln->seq[seqidx].gapshift[pair.ntj];
    gs3 = aln->seq[seqidx].gapshift[pair.nti];


    return aln->seq[seqidx].bp[pair.nti-gs3][pair.ntj-gs5];
}
/***   GetBpFromAlnSeq   ***/



/***   GetConsBp   **************************************************
 *
 * Same as GetBpFromAlnSeq, just for convinience
 *
 */
cons_bpair_t *
GetConsBp(aln_t *aln, pair_t pair)
{
    return &aln->cons_seq->cons_bp[pair.nti][pair.ntj];
}
/***   GetConsBp   ***/



/***   IsBasepair   ***********************************************************
 *
 * return 1 if two nucleotides, may form a basepair
 *
 * wobble-bp (if allow_wobble==1) and DNA allowed
 * case-insensitiv
 */
int
IsBasepair(char nt1, char nt2)
{
    const int allow_wobble = 1;
    int      is_bp = 0;

    /* case insensitiv */
    nt1=(char) tolower((int) nt1);
    nt2=(char) tolower((int) nt2);

    /* allow dna */
    if (nt1=='t')
        nt1='u';
    if (nt2=='t')
        nt2='u';

    /* a:u */
    if      ( (nt1=='a' && nt2=='u') || (nt1=='u' && nt2=='a') )
    {
        is_bp = 1;
    }
    /* g:c */
    else if ( (nt1=='c' && nt2=='g') || (nt1=='g' && nt2=='c') )
    {
        is_bp = 1;
    }
    /* g:u (wobble) */
    else if ( (nt1=='g' && nt2=='u') || (nt1=='u' && nt2=='g') )
    {
        if (allow_wobble)
            is_bp = 1;
    }
    else
    {
         is_bp = 0;
    }
    return is_bp;
}
/***   IsBasepair   ***/




/***   Degap   ****************************************************************
 *
 * Degaps provided sequence string and returns the string
 * out must be preallocated
 *
 */
void
Degap(char *aln_seq, char *out)
{

    int  i=0;
    int  counter=0;


    for (i=0; i<strlen(aln_seq); i++)
    {
        if (! IS_GAP(aln_seq[i]) )
        {
            out[counter] = aln_seq[i];
            counter++;
        }
    }
    out[counter] = '\0';
}
/***   Degap   ***/


