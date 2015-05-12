/******************************************************************************
* 
* mic.h - routines for the mutual information content (single nucleotide dependencies)
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
 *  CVS $Id: mic.h,v 1.13 2007-10-22 10:43:23 steger Exp $    
 */



#ifndef MIC_H_INCL
#define MIC_H_INCL



/***   Calculate single nucleotide dependencies
 *
 *  using unbiased probability estimation
 *   (see Chiu & Kolodziejczak, CABIOS Vol.7, no.3, 1991 (347-352))
 *  using maximum likelyhood estimation
 *   (see Gutell, Power, Hertz, Putz & Stormo (1992))
 *
 ***   Prob estimation factor
 *
 *  do calculations with log_2
 *   (see Schneider, Stormo, Gold & Ehrenfeucht (1986))
 *  do calculations with log_e
 *   (see Chiu & Kolodziejczak (1991))
 *
 */

/* counted characters: -,a,c,g,u
   (all others will be proportionate converted)
*/

/* frontend
 * alifoldscore: 0 use mic instead of alifoldscore:
 *               1 use alifold scoring instead of mic
 *               2 use alifold score incl stacking
 */                   
extern void
ComputeMutualInfoContent(char **alnseq, int alnlen, int n_seq,
                         int unbiased, int bit, int pair_entropy_norm,
                         int alifoldscore);

/* get the consensus mic for a basepair-position, that is
 * a particular mic (GetMIC) normalized to 1
 * return 0.0 if <= mic_thresh
 * nt_i and nt_j are zero-offset nucleotide column indices
 */                   
extern float
GetConsMIC(pair_t pair, float mic_thresh);
   
   
/* get the mic for a basepair-position
 * nt_i and nt_j are zero-offset nucleotide column indices
 */                   
extern float
GetMIC(pair_t pair);


/* test if mutinfo has already been computed
 */
extern int
MIC_IsComputed(void);

extern float
GetMaxMIC(void);


extern void
PrintMic(char *csseq, FILE *stream);

#endif
