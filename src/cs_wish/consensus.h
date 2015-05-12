/******************************************************************************
* 
* consensus.h - routines for creating consensus sequences and consensus
*               basepair probability matrices
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
 *  CVS $Id: consensus.h,v 1.11 2004-05-25 13:29:13 wilm Exp $    
 */



#ifndef CONSENSUS_H_INCL
#define CONSENSUS_H_INCL




/* Calculate consensus sequence and it's bp probabilites
 * Caller must free
 */
extern cons_seq_t *
CalcConsensus(aln_t *aln, float *weight);

/* Compute the consensus sequence nt's of given range (inclusive)
 * Caller must free
 */
extern void
CalcConsSeq(aln_t *aln, int min_nt, int max_nt,  float *weight, char *cons_seq_nt);

/* Compute the number of gaps at each position of given range (inclusive)
 * and write them to last argument
 */
extern void
CalcConsNumGaps(aln_t *aln, int min_nt, int max_nt, int *ret_ngaps);

/* Returns consensus prob of a bp
 * or 0 if below prob_limit
 * nt_j>nt_i
 */
extern float
GetConsBpProb(aln_t *aln, pair_t pair, float prob_limit);

/* Pretty print out matrix to stream */
extern void
PrintCsTdProbMat(cons_seq_t *csseq, FILE *stream);

/* Free alle probs and their contrib lists */
extern void
FreeConsBpMx(cons_bpair_t **bpmx, int seqlen);

#endif

