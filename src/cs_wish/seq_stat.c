/******************************************************************************
* 
* seq_stat.c - routines for sequence alignment statistics
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
 *  CVS $Id: seq_stat.c,v 1.9 2004-05-25 13:29:17 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE
  
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h> /* tolower */

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "seq_stat.h"


/***   private
 */

#if 0
    #define DEBUG
#endif

static int
WeightNt(char s_il, char s_jl);

static void
MakeIdentityMx(char **aseqs, int num, double **ret_imx);

static double
PairwiseIdentity(char *s1, char *s2);

/*
 ***/



/***   Sop   ******************************************************************
 *
 * Calculate the Sum-of-pairs-Cost using WeightNt as weighting scheme
 * sop is returned
 * all column costs/scores will be written to preallocated int *col_cost
 * returns -1 on error
 *
 */
int
Sop(char **aln_seq, int aln_len, int n_seq, int *col_cost)
{
    int l, i, j;
    int col_score, sum;
    char nt[2];
    

    /***** for all columns
     */
    sum = 0;
    for (l=0; l<aln_len; l++)
    {
        col_score = 0;        

        /* for all pairs i,j in one column
         */        
        for (i=0; i<n_seq-1; i++)
        {
            for (j=i+1; j<n_seq; j++)
            {
                nt[0] = (char)tolower(aln_seq[i][l]);
                nt[1] = (char)tolower(aln_seq[j][l]);
                
                col_score += WeightNt(nt[0], nt[1]);
            }
        }
        #ifdef DEBUG
            DEBUG_P("col_score[%d]=%d\n", l, col_score);
        #endif
        col_cost[l] = col_score;
        sum += col_score;
    }
    
    return sum;
}
/***   Sop   ***/




/***   WeightNt   *************************************************************
 *
 * a simple 'weight-function' for two nucleotides
 * valid nucleotide symbols are acgun-ry (lowercase!)
 *
 */
int
WeightNt(char s_il, char s_jl)
{
    /* gaps don't match */
    if        ( s_il=='-' || s_jl=='-')
        return 0;
    
    /* 'n' matches all */
    else if ( s_il=='n' || s_jl=='n')
        return 1;
    
    /* a real match */
    else if ( s_il == s_jl )
       return 1;
        
    /* s_il is a purine and s_jl matches */
    else if ( s_il=='r' && (s_jl=='a' || s_jl=='g'))
        return 1;
    
    /* s_jl is a purine and s_il matches */
    else if ( s_jl=='r' && (s_il=='a' || s_il=='g'))
        return 1;
    
    /* s_il is a pyrimidine and s_jl matches */
    else if ( s_il=='y' && (s_jl=='c' || s_jl=='u'))
        return 1;
        
    /* s_jl is a pyrimidine and s_il matches */
    else if ( s_jl=='y' && (s_il=='c' || s_il=='u'))
        return 1;
        
    else
        return 0;
}
/***   WeightNt   ***/




/***   PwIdent   ************************************************************
 *
 * Purpose:     Calculate similarity of aligned sequence set.
 *              Adopted from Sean Eddys squid-library (alistat_main.c)
 *              (see below)
 *
 * Args:        Alignment aseq and number of sequences nseq
 *
 * Return:      Relative average identity 
 *
 */
double
PwIdent(char **aseq, int alnlen, int nseq)
{
    double   **imx;          /* identity matrix  */
    double     sum = 0;
	double     avg_pwident;   /* average identity */
    int   i,j;

    imx = (double**) mx_new(nseq+1, nseq+1, sizeof(double), sizeof(double*), MX_DP);
    
    /* Calculate identity matrix where each sequence pair index
     * is set to (idents / MIN(len1, len2))
	 * (imx has zero offset and is symmetric)
     */
    MakeIdentityMx(aseq, nseq, imx);
    
	/*** Sum up identity and divide it by number of possible sequence pairs
	 *
     */
    for (i = 0; i < nseq-1; i++)
        for (j = i+1; j < nseq; j++)
            sum += imx[j][i];
   
    avg_pwident = sum / (double) (nseq * (nseq-1)/2.0);
    
    
    
    mx_free((void**)imx, nseq+1, MX_DP);
    
    return avg_pwident;
}
/***   PwIdent   ***/




/*****************************************************************
 *
 * MakeIdentityMx, PairwiseIdentity are from
 *                   Sean Eddy's HMMER (squid)
 *
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1999 Washington University School of Medicine
 * All Rights Reserved
 * 
 *
 *****************************************************************/


/***   MakeIdentityMx   *******************************************************
 *
 * Purpose:  Given a set of aligned sequences, construct
 *           an NxN fractional identity matrix. (i.e. 1.0 is
 *           completely identical, 0.0 is completely different).
 *           
 * Args:     aseqs        - flushed, aligned sequences
 *           num          - number of aseqs
 *           ret_imx      - must already be allocated
 *           
 */
void
MakeIdentityMx(char **aseqs, int num, double **ret_imx)
{
    int     i,j;               /*  counters over sequences  */


    /* Calculate distances, symmetric matrix
     */
    for (i = 0; i < num; i++)
        for (j = i+1; j < num; j++)
            ret_imx[j][i] = PairwiseIdentity(aseqs[i], aseqs[j]);

    return;
}
/***   MakeIdentityMx   ***/


/***   PairwiseIdentity   ******************************************************
 * 
 * Purpose:  Calculate the pairwise fractional identity between
 *           two aligned sequences s1 and s2. This is simply
 *           (idents / MIN(len1, len2)).
 *
 *           Note how many ways there are to calculate pairwise identity,
 *           because of the variety of choices for the denominator:
 *           idents/(idents+mismat) has the disadvantage that artifactual
 *             gappy alignments would have high "identities".
 *           idents/(AVG|MAX)(len1,len2) both have the disadvantage that 
 *             alignments of fragments to longer sequences would have
 *             artifactually low "identities".
 *           
 *           Case sensitive; also, watch out in nucleic acid alignments; 
 *           U/T RNA/DNA alignments will be counted as mismatches!
 */
double
PairwiseIdentity(char *s1, char *s2)
{
    int     idents;        /*  total identical positions  */
    int     len1, len2;    /*  lengths of seqs            */
    int     x;             /*  position in aligned seqs   */

    idents = len1 = len2 = 0;
    for (x = 0; s1[x] != '\0' && s2[x] != '\0'; x++)
    {
        if (!IS_GAP(s1[x]))
        {
            len1++;
            if (s1[x] == s2[x])
                idents++; 
      }
      if (!IS_GAP(s2[x]))
            len2++;
    }
    if (len2 < len1)
        len1 = len2;
    
    return (len1 == 0 ? 0.0 : (double) idents / (double) len1);
}
/***   PairwiseIdentity   ***/
