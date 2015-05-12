/******************************************************************************
* 
* utils.h - various routines
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
 *  CVS $Id: utils.h,v 1.23 2007-10-22 10:43:23 steger Exp $    
 */


#ifndef UTILS_H_INCL
#define UTILS_H_INCL


#define IS_GAP(c)           ((c) == '-' || (c) == '.')


/* frontends to sr eddy's stopwatch
 */
extern Stopwatch_t *
TimerStart(void);

extern void
TimerStop(const char *prefix, Stopwatch_t *s);

/* Return pointer to basepair from indices
 * of an aligned sequence
 */
extern bpair_t *
GetBpFromAlnSeq(aln_t *aln, int seqidx, pair_t pair);

/* Same as GetBpFromAlnSeq, just for convinience
 */
extern cons_bpair_t *
GetConsBp(aln_t *aln, pair_t pair);

/* Degaps provided sequence string and returns the string
 * Caller must free
 */
extern void
Degap(char *aln_seq, char *out);

/* Sets a pair with the corresponding indices
 */
void
SetPair(int nti, int ntj, pair_t *p);

/* Compute the gapshift from an aligned sequence and return it as int array
 * Caller must free
 */
extern int *
CreateGapShift(char *aln_seq);

/* return 1 if two nucleotides, may form a basepair */
extern int
IsBasepair(char nt1, char nt2);


#endif
