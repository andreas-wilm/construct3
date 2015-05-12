/******************************************************************************
* 
* consensus.c - routines for creating consensus sequences and consensus
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
 *  CVS $Id: consensus.c,v 1.33 2005-10-20 11:43:59 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>
#include <tk.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"
#include "proj.h"

#include "consensus.h"



/***   private
 */

#if 0
    #define DEBUG
#endif

#define NT_GAP 0
#define NT_A   1
#define NT_C   2
#define NT_G   3
#define NT_U   4


static char
ConsNt(float *nt_freq);
 
/*
 ***/



/***   FreeConsBpMx   ***
 * 
 * free all consensus dlists and the matrix itself
 */
void
FreeConsBpMx(cons_bpair_t **bpmx, int seqlen)
{
    int nt_i=0, nt_j=0;
    
    if (bpmx==NULL)
        return;
    for (nt_i=0; nt_i<seqlen; nt_i++)
        for (nt_j=nt_i+1; nt_j<seqlen; nt_j++)
            dlist_destroy(bpmx[nt_j][nt_i].contribs);

    mx_free((void**)bpmx, seqlen, MX_DP);
}
/***   FreeConsBpMx   ***/



/***   CalcConsensus   ********************************************************
 *
 * calculate consensus sequence and it's bp probabilites
 *
 * Caller must free
 */
cons_seq_t *
CalcConsensus (aln_t *aln, float *weight)
{
    
    int        seq_idx    = 0;    /* iteration index over all sequences    */
    float      csbp_prob  = 0.0;  /* consensus probability of one basepair */
    cons_seq_t *ret_cons  = NULL;  /* returned consensus sequence_t         */
    bpair_t    *bp;
    pair_t      pair;
    int nt_i=0, nt_j=0;
    Stopwatch_t *watch;
    int *contrib;
    
    
    watch = TimerStart();
            
    ret_cons = (cons_seq_t *) Scalloc(1, sizeof(cons_seq_t));

    ret_cons->id = strdup("Consensus");
    ret_cons->nt = (char *) Scalloc(aln->len+1, sizeof(char));
    CalcConsSeq(aln, 0, aln->len-1, weight, ret_cons->nt);
    ret_cons->nt[aln->len] = '\0';
    
    ret_cons->num_gaps = (int *) Scalloc(aln->len, sizeof(int));
    CalcConsNumGaps(aln, 0, aln->len-1, ret_cons->num_gaps);
    
    ret_cons->cons_bp = (cons_bpair_t **) \
                        mx_new (aln->len, aln->len, sizeof(cons_bpair_t), sizeof(cons_bpair_t*), MX_DP);
    DEBUG_P("Consensus Sequence = \"%s\"\n", ret_cons->nt);

        
        
    /*****   setup cons_basepair_t
     *
     */
    for (nt_i=0; nt_i<aln->len; nt_i++)
    {
        for (nt_j=nt_i+1; nt_j<aln->len; nt_j++)
        {
            ret_cons->cons_bp[nt_j][nt_i].contribs = dlist_new(sizeof(int*), free);
                        
            /* Set up prob later */
        }
    }

    /*****   setup prob and contrib list foreach bp
     *
     */
    
    /* foreach basepair
     */
    for (nt_i=0; nt_i<aln->len; nt_i++)
    {
        for (nt_j=nt_i+1; nt_j<aln->len; nt_j++)
        {
            /* reset */
            csbp_prob  = 0.0;
            SetPair(nt_i, nt_j, &pair);
            
            /* foreach sequence
             */
            for (seq_idx=0; seq_idx < aln->num_seq; seq_idx++)
            {   
                bp = GetBpFromAlnSeq(aln, seq_idx, pair);
                if (bp==NULL)
                    continue;
                if (bp->prob==0.0)
                    continue;

                csbp_prob += (pow(bp->prob, 1.0/(double)POW_A)* weight[seq_idx]);

#ifdef DEBUG
                DEBUG_P("appending seq_idx=%d with prob=%f and weight=%f to cons nt_i=%d nt_j=%d\n",
                                                                  seq_idx, bp->prob, weight[seq_idx], nt_i, nt_j);
#endif

                contrib = Scalloc(1, sizeof(int));
                *contrib=seq_idx;
                if (dlist_append(ret_cons->cons_bp[nt_j][nt_i].contribs,
                                 contrib) == DLIST_ERROR)
                    ERROR_P("couldn't dlist_append to ret_cons->cons_bp[%d][%d].contribs\n", nt_j, nt_i);
                
            
            } /* end : foreach sequence */
            
    
            ret_cons->cons_bp[nt_j][nt_i].prob = csbp_prob;
        }
    }
    /* end : foreach basepair */
    
   
    
    TimerStop(__FUNCTION__, watch);
    
    return ret_cons;
}
/***   CalcConsensus   ***/





/***   CalcConsSeq   **********************************************************
 *
 * Compute the consensus sequence nt's
 * Done by simple majority rule
 *
 * no trailing \0 is added !
 *
 * Sequences must be lowercase
 *
 * Range: min_nt, max_nt (zero-offset, inclusive)
 *
 * User must free
 *
 */
void
CalcConsSeq(aln_t *ali, int min_nt, int max_nt,  float *weight, char *cons_seq_nt)
{
    int   seq_idx = 0;       /* iteration index over all sequences    */
    int   nt_idx  = 0;       /* iteration index over all nts          */
    char  nt;
    float nt_freq[5];
    float weight_sum;
    
#ifdef DEBUG
        DEBUG_P("update from min_nt=%d  to max_nt=%d\n", min_nt, max_nt);
#endif

    weight_sum = GetSumOfWeights();
    
    /*** foreach nt == column
     */
    for (nt_idx=min_nt; nt_idx <= max_nt; nt_idx++)
    {
        nt_freq[NT_A]   = 0.0;
        nt_freq[NT_C]   = 0.0;
        nt_freq[NT_G]   = 0.0;
        nt_freq[NT_U]   = 0.0;
        nt_freq[NT_GAP] = 0.0;
        
        /*** foreach sequence
         */    
        for (seq_idx=0; seq_idx < ali->num_seq; seq_idx++)
        {   
            nt = (char) tolower(ali->seq[seq_idx].nt[nt_idx]);

            switch (nt)
            {
                case 'a':
                    nt_freq[NT_A]+=weight[seq_idx];
                    break;
                case 'c':
                    nt_freq[NT_C]+=weight[seq_idx];
                    break;
                case 'g':
                    nt_freq[NT_G]+=weight[seq_idx];
                    break;
                case 'u':
                    nt_freq[NT_U]+=weight[seq_idx];
                    break;
                case 'r':  /* -> y:purine     : a || g*/
                    nt_freq[NT_A]+=(weight[seq_idx]*0.5);
                    nt_freq[NT_G]+=(weight[seq_idx]*0.5);
                    break;
                case 'y':  /* -> r:pyrimidine : u || c*/
                    nt_freq[NT_U]+=(weight[seq_idx]*0.5);
                    nt_freq[NT_C]+=(weight[seq_idx]*0.5);
                    break;
                case 'n':
                    nt_freq[NT_A]+=(weight[seq_idx]*0.25);
                    nt_freq[NT_C]+=(weight[seq_idx]*0.25);
                    nt_freq[NT_G]+=(weight[seq_idx]*0.25);
                    nt_freq[NT_U]+=(weight[seq_idx]*0.25);
                    break;
                default:
                    if ( IS_GAP(nt) )
                    {
                        nt_freq[NT_GAP]+=weight[seq_idx];
                    }
                    else
                    {
                        WARN_P("ignoring unexpected nt character \"%c\"\n", nt);
                    }
                    break;
            }
            
        } /* foreach sequence */
        
        
        /* compute fractions
         */
        nt_freq[NT_A]   /= weight_sum;
        nt_freq[NT_C]   /= weight_sum;
        nt_freq[NT_G]   /= weight_sum;
        nt_freq[NT_U]   /= weight_sum;
        nt_freq[NT_GAP] /= weight_sum;
        
#ifdef DEBUG
            DEBUG_P("%03d: A=%0.5f C=%0.5f G=%0.5f U=%0.5f GAP=%0.5f\n",
                                           nt_idx, nt_freq[NT_A], nt_freq[NT_C],
                                  nt_freq[NT_G], nt_freq[NT_U], nt_freq[NT_GAP]);            
#endif
	    
        cons_seq_nt[nt_idx] = ConsNt(nt_freq);

#ifdef DEBUG
            DEBUG_P("cons_seq_nt[%d]=%c\n", nt_idx, cons_seq_nt[nt_idx]);
#endif


    } /* foreach nt == column */
}
/***   CalcConsSeq   ***/







/***   ConsNt   ***************************************************************
 * one nt is in majority:
 *  >50%: ACGU-
 *  <50%: acgun
 *
 * two nts are equally frequent:
 *  a/g -> r
 *  u/c -> y
 *  n otherwise
 *
 * otherwise:
 *  -
 */
static char
ConsNt(float *nt_freq)
{
    char  cons_nt  = '#';
    int   best_idx = 0;
    float max_freq = 0.0;
    int   nof_sims = 0;
    int   idx      = 0;
        
    for (idx=0; idx<5; idx++)   
    {
        if (nt_freq[idx]>max_freq)
        {
            max_freq = nt_freq[idx];
            best_idx = idx;
        }
    }

    /* shared first place?
     */
    nof_sims=-1;
    for (idx=0; idx<5; idx++)   
    {
        DEBUG_P("max_freq=%f nt_freq[%d]=%f\n", max_freq, idx, nt_freq[idx]);
        if (max_freq==nt_freq[idx])
            nof_sims++;
    }
    
#ifdef DEBUG
    for (idx=0; idx<5; idx++) {
        DEBUG_P("idx=%d nt_freq=%f  max_freq=%f  nof_sims=%d\n",
                idx, nt_freq[idx], max_freq, nof_sims);
    }
#endif
    
    /* one nt rules them all and got more then 50%
     */
    if (max_freq>=0.5)
    {
        if (best_idx==NT_A)
            cons_nt = 'A';
        else if (best_idx==NT_C)
            cons_nt = 'C';
        else if (best_idx==NT_G)
            cons_nt = 'G';
        else if (best_idx==NT_U)
            cons_nt = 'U';
        else if (best_idx==NT_GAP)
            cons_nt = '-';
    }
    /* one nt rules them all but got less then 50%
     */
    else if (nof_sims==0)
    {
        if (best_idx==NT_A)
            cons_nt = 'a';
        else if (best_idx==NT_C)
            cons_nt = 'c';
        else if (best_idx==NT_G)
            cons_nt = 'g';
        else if (best_idx==NT_U)
            cons_nt = 'u';
        else if (best_idx==NT_GAP)
            cons_nt = '-';
    }
    /* there are two equally prominent nts
     */
    else if (nof_sims)
    {
        if ((nt_freq[NT_A] == nt_freq[NT_G]) && (nt_freq[NT_U] == nt_freq[NT_C]) && (nt_freq[NT_A] == nt_freq[NT_U]))
            cons_nt = 'n';
        
        /* case 'r':   -> a || g
           case 'y':   -> u || c
        */      
        else if ((best_idx==NT_A || best_idx==NT_G)
                            &&
              nt_freq[NT_A] == nt_freq[NT_G])
            cons_nt = 'r';

        else if ((best_idx==NT_U || best_idx==NT_C)
                                 &&
                   nt_freq[NT_U] == nt_freq[NT_C])
            cons_nt = 'y';
        else
            cons_nt = 'n';          
    }
    else
    {
        cons_nt = '-';
    }

    return cons_nt;
}
/***   ConsNt   ***/





/***   CalcConsNumGaps   ******************************************************
 *
 * Compute the number of gaps at each position
 *
 * range: min_nt, max_nt (inclusive)
 *
 */
void
CalcConsNumGaps(aln_t *ali, int min_nt, int max_nt, int *ret_ngaps)
{
    int nt_idx     = 0;
    int seq_idx    = 0;
    int gap_count  = 0;

    
    /*** foreach nt == column
     */
    for (nt_idx=min_nt; nt_idx <= max_nt; nt_idx++)
    {   
        gap_count = 0;
        
        /*** foreach sequence
         */    
        for (seq_idx=0; seq_idx < ali->num_seq; seq_idx++)
        {    
            if (IS_GAP(ali->seq[seq_idx].nt[nt_idx]))
                gap_count++;
        }
        ret_ngaps[nt_idx] = gap_count;
    }
}
/***   CalcConsNumGaps   ***/





/***   PrintCsTdProbMat   *********************************************************
 *
 * Prints the Consensus Basepair-ProbMat to console
 *
 *
 */
void
PrintCsTdProbMat(cons_seq_t *csseq, FILE *stream)
{
    int seqlen = 0;
    int i=0, j=0;
    pair_t pair;
    
    seqlen = strlen(csseq->nt);

    
    /* print sequence nt
     */
    fprintf(stream, "   ");
    for (i=0; i<seqlen ; i++)
        fprintf(stream, "  %4c  ", csseq->nt[i]);
    fprintf(stream, "\n");



    for (i=0; i<seqlen ; i++)
    {
    
   	    fprintf(stream, "%c  ", csseq->nt[i]);
        for (j=0; j<seqlen ; j++)
        {
            if (i>j)
            {
                fprintf(stream, "  ----  ");
            }
            else if (j==i)
            {
                fprintf(stream, "  ++++  ");
            }
            else
            {
                SetPair(i, j, &pair);
                fprintf(stream, " %0.4f ", GetConsBpProb(aln, pair, 0.0));
            }
        }
        fprintf(stream, "\n");
    }    
    fflush(stream);
}
/***   PrintCsTdProbMat   ***/





/***   GetConsBpProb   ********************************************************
 *
 * nt_i<nt_j = zero offset nt indices
 * prob must be greater then prob_thresh, otherwise 0.0 will be returned
 */
float
GetConsBpProb (aln_t *aln, pair_t pair, float prob_limit)
{
    float cs_prob;
    float norm_factor = GetSumOfWeights();
    
    if (aln->cons_seq->cons_bp[pair.nti][pair.ntj].prob==0)
        return 0.0;
        
    cs_prob = pow(aln->cons_seq->cons_bp[pair.nti][pair.ntj].prob / norm_factor, (double)POW_B);
    if (cs_prob <= prob_limit)
        cs_prob = 0.0;

    return cs_prob;
}
 /***   GetConsBpProb   ***/           

