/******************************************************************************
* 
* bp_prob_mat.h - routines for reading basepair probability matrices
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
 *  CVS $Id: bp_prob_mat.h,v 1.10 2004-05-25 13:29:13 wilm Exp $    
 */


#ifndef BP_PROB_MAT_H_INCL
#define BP_PROB_MAT_H_INCL


extern void
FreeBpProbMx(bpair_t ***bpmx, int seqlen);

/* Read given rnafold basepair matrix file (unaligned),
   cutoff probs and return dp bpair_t matrix
   seq must be raw and not aligned
  */
extern bpair_t ***
SetupBpProbMx (char *path2mat, char *seq, float prob_cutof);

/* print the td probmat of one sequence to stream
 */
extern void
PrintTdProbMat(char *nalseq, bpair_t ***bpmx, FILE *stream);



#endif
