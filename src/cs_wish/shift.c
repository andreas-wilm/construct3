/******************************************************************************
* 
* shift.c - compute cs-shift between two alignments
*           that is the avg. number of move operations in ConStruct
*           to transform a predicted alignment into a trusted one
*          
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
 *  CVS $Id: shift.c,v 1.7 2004-11-11 13:39:55 wilm Exp $    
 */



#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>
#include <string.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"


#include "shift.h"



/***   private
 */
 
#if 0
    #define DEBUG
#endif

#define PRETTY_PRINT_SHIFT 1

/* return 0 if 0
 *       -1 if negative (left)
 *        1 if positive (right)
 *         
 */
#define GET_ORIENTATION(shift) ((shift)==0 ? 0 : ((shift)/abs(shift)));



/*** shift_seq
 *
 */
typedef struct
{
    short aln_pos;
    char   nt;
} residue_t;


typedef struct
{
    short    aln_len;
    short    nal_len;
    residue_t *residue; /* unaligned */
} shift_seq_t;



/*** block
 */
typedef struct
{
    int start;
    int end;
    int offset;
} block_move_t;





/*******************   Shift
 */

/* Construct a shift_seq_t from a aligned sequence */
static shift_seq_t *
Construct_ShiftSeqType(seq_t aln_seq);

/* Free the constructed shift_seq_t */
static void
Destruct_ShiftSeqType(shift_seq_t *shift_seq);

/* looks up a block of equally orientated non-0 shifts
   and store its start, its end and the minimum offset in block_move
*/
static int
FindBlock(int *rawshift, int start, int end, block_move_t *block_move);

/* looks up the new start index after a block has been moved */
static void
SetNewBlockStart(int *rawshift, int *start, int end);

/* moves a block via blockmove */
static void
MoveBlock(int *rawshift, block_move_t blockmove);

/* Calculates the score from a rawshift */
static int
ScoreFromRawShift(int *pwshift, int size);

/* Inserts Pairmask to rawshift */
static void
InsertPairmask(shift_seq_t *trust, int *pairmask, int *rawshift);


/* FIXME: Doku */
static int
SeqCmp(const void *seqt1, const void *seqt2);


/***   Scoring
 */

/* calculate the shift (int array of moves) of to shift_seqs */
static int
CalcRawShift(shift_seq_t *trust, shift_seq_t *pred, int *pairmask);






/***   CsShiftScore   *********************************************************
 *
 * Sequences must be sorted and of same case
 * NULL on error
 * FIXME: pairmask int array 1 for paired 0 for unpaired, idx: 0 to aln_len -1 
 */
int *
CsShiftScore(aln_t *aln_trust, aln_t *aln_pred, int *pairmask)
{
    int  *score;
    char *nal_seq_trust, *nal_seq_pred, *id;
    shift_seq_t *pred_seq, *trust_seq;
    int i;


    /***   check number of sequences
     */
    if (aln_pred->num_seq!=aln_trust->num_seq) {
        ERROR_P("%s\n", "Number of sequences differ!");
        return NULL;
    }
    
        
    /* sort alignments by nts's */
    qsort(aln_trust->seq, aln_trust->num_seq, sizeof(seq_t), SeqCmp);
    qsort(aln_pred->seq, aln_pred->num_seq, sizeof(seq_t), SeqCmp);


    /* check seqs
     */
    nal_seq_trust = (char*) Scalloc(aln_trust->len+1, sizeof(char));
    nal_seq_pred  = (char*) Scalloc( aln_pred->len+1, sizeof(char));
    for (i=0; i<aln_pred->num_seq; i++) {
        Degap(aln_trust->seq[i].nt, nal_seq_trust);
        Degap(aln_pred->seq[i].nt,  nal_seq_pred);
       
        if ( ! STR_EQ(nal_seq_trust, nal_seq_pred)) {
            ERROR_P("Raw sequences differ (at least for seq %s)\n", aln_trust->seq[i].id);
            return NULL;
        }
    }
    free(nal_seq_trust);
    free(nal_seq_pred);
    
    
    

    if (pairmask==NULL)
        VERBOSE_P("not using a pairmask");
    else
        VERBOSE_P("using a pairmask");
    
    pred_seq  = (shift_seq_t *) calloc(1, sizeof(shift_seq_t));
    trust_seq = (shift_seq_t *) calloc(1, sizeof(shift_seq_t));
    
    
    /*****   foreach sequence
     *
     */
    score = Scalloc(aln_pred->num_seq, sizeof(int));
    for (i=0; i<aln_pred->num_seq; i++)
    {
        id = aln_trust->seq[i].id;    /* or aln_pred */
        
        trust_seq = Construct_ShiftSeqType(aln_trust->seq[i]);
        pred_seq  = Construct_ShiftSeqType(aln_pred->seq[i]);
        
        VERBOSE_P("Computing Shift for %s\n", id);

        score[i] = CalcRawShift(trust_seq, pred_seq, pairmask);
        
        VERBOSE_P("Shift for %s = %d\n", id, score[i]);

        
        Destruct_ShiftSeqType(pred_seq);
        Destruct_ShiftSeqType(trust_seq);
            
    } /* foreach sequence */
    

    return score;
}
/***   CsShiftScore   ***/





/***   Construct_ShiftSeqType   *****************************************************
 *
 * Takes two aligned seq_t's as input and construct a new shift_seq_t
 *
 */
shift_seq_t *
Construct_ShiftSeqType(seq_t aln_seq)
{
    shift_seq_t *shiftseq  = NULL;  /* return */
    int       aln_idx   = 0;
    int       nal_idx   = 0;
    
    
    shiftseq = (shift_seq_t *) calloc(1, sizeof(shift_seq_t));

    shiftseq->aln_len  = strlen(aln_seq.nt);

    shiftseq->residue  = (residue_t *) calloc(shiftseq->aln_len, sizeof(residue_t));
    
    
    nal_idx = 0;
    for (aln_idx=0; aln_idx < shiftseq->aln_len; aln_idx++)
    {
        
        if ( ! IS_GAP(aln_seq.nt[aln_idx]) )
        {    
            shiftseq->residue[nal_idx].aln_pos = aln_idx;
            shiftseq->residue[nal_idx].nt = aln_seq.nt[aln_idx];

            nal_idx++;
        }
    }
    shiftseq->nal_len = nal_idx;
    shiftseq->residue  = (residue_t *) realloc(shiftseq->residue,
                                               shiftseq->nal_len * sizeof(residue_t));

        #ifdef DEBUG    
            DEBUG_P("shiftseq for id %s\n", aln_seq.id);
            DEBUG_P("\taligned seq = %s\n", aln_seq.nt);
            fprintf(stdout, "\tshift   seq = ");
            for (nal_idx=0; nal_idx < shiftseq->nal_len; nal_idx++)
                fprintf(stdout, " %c:%d ", shiftseq->residue[nal_idx].nt,
                                           shiftseq->residue[nal_idx].aln_pos);
            fprintf(stdout, "\n");
            fflush(stdout);
       #endif

    return shiftseq;
    
}
/***   Construct_ShiftSeqType   ***/



/***   Destruct_ShiftSeqType   ************************************************
 */
void
Destruct_ShiftSeqType(shift_seq_t *shiftseq)
{
    free(shiftseq->residue);
    free(shiftseq);
    
}
/***   Destruct_ShiftSeqType   ***/





/***   InsertPairmask   *******************************************************
 *
 * We only want basepaired nts to contribute to a shift
 * In other words: If a nt is not basepaired , it shouldn't contribute to a block
 * -> set this position to its first left/right neighbour which is 0 or paired
 *
 */
void
InsertPairmask(shift_seq_t *trust, int *pairmask, int *rawshift)
{
    int nal_idx = 0;    /* not aligned pos as loop counter */
    int aln_pos = 0;    /* corresponding aligned postition */
    int success = 0;    /* altering successfull            */
    int i = 0;

    for (nal_idx=0; nal_idx < trust->nal_len; nal_idx++)
    {
        aln_pos = trust->residue[nal_idx].aln_pos;

        /* only change if not basepaired
         */
        if ( ! pairmask[aln_pos])
        {     
            #ifdef DEBUG
                DEBUG_P("Mask alter at aln_pos %d==%c\n",
                         aln_pos, trust->residue[nal_idx].nt);
            #endif

            success = 0;

            /*** lookup to left
             */
            for (i=nal_idx; i>=0; i--)
            {
                if (pairmask[trust->residue[i].aln_pos])
                {
                    #ifdef DEBUG
                        DEBUG_P("Next Bp at aln_pos %d==%c (alter from %d to %d)\n",
                                 trust->residue[i].aln_pos, trust->residue[i].nt,
                                 rawshift[nal_idx], rawshift[i]);
                    #endif

                    rawshift[nal_idx] = rawshift[i];
                    success = 1;
                    break;
                }
            }
            
            /*** if not successfull lookup to right
             */
            if ( ! success )
            {
                for (i=nal_idx; i<trust->nal_len; i++)
                {
                    if (pairmask[trust->residue[i].aln_pos])
                    {
                        #ifdef DEBUG
                            
                            DEBUG_P("Next Bp at aln_pos %d==%c (alter from %d to %d)\n",
                                     trust->residue[i].aln_pos, trust->residue[i].nt,
                                     rawshift[nal_idx], rawshift[i]);
                        #endif
                        rawshift[nal_idx] = rawshift[i];
                        success = 1;
                        break;
                    }
                }
            }
            
            /* something went wrong
             */
            if ( ! success )
            {
                ERROR_P("%s\n", "Couldn't find next paired nt ! Maybe your pairmask file is corrupted ?");
            }
        }
    }
}
/***   InsertPairmask   ***/





/***   FindBlock   ************************************************************
 *
 * Creates a block_move_t
 * which defines the movement of a block inside start-end
 *
 * end is inclusive: rawshift(0 - end)
 *
 */
static int
FindBlock(int *rawshift, int start, int end, block_move_t *block_move)
{
    int     idx            = 0;
    int     orientation    = 0;
    int     old_orient     = 0;
    int     inside_block   = 0;
    int     min_shift      = 65535; 

    
    for (idx=start; idx<=end; idx++)
    {
        orientation = GET_ORIENTATION(rawshift[idx]);
        #ifdef DEBUG        
            DEBUG_P("idx=%d -> %d\n", idx, rawshift[idx]);
        #endif

        /* save minimum element
         */
        if (rawshift[idx] != 0.0)
        {
            min_shift = MIN(min_shift, abs(rawshift[idx]));
        }
        
            
        /* check if this is the beginning of a block
         */
        if (orientation != 0)
        {        
            if (inside_block == 0)
            {
                block_move->start = idx;
                #ifdef DEBUG
                    DEBUG_P("new block start = %d\n", idx);
                #endif
            }
            inside_block = 1;
        }
        
        
        /* first element needs special treatment
           since no old orientation is known
        */
        if (idx == start)
        {
            old_orient = orientation;
            continue;
        }
        
        
        /* check if this is the end of a block
         */
        if (  (orientation == 0 && inside_block == 1)
            ||
              (orientation != old_orient && old_orient != 0) )
        {
            block_move->end    = idx-1;
            block_move->offset = min_shift * old_orient;
            #ifdef DEBUG
                DEBUG_P("%s\n", "normal block end");
            #endif
            return 1;
        }
        
        old_orient = orientation;
    }
    
    if (inside_block)
    {
        block_move->end    = idx-1;
        block_move->offset = min_shift * old_orient;
        #ifdef DEBUG
            DEBUG_P("%s\n", "total block end");
        #endif
        return 1;
    }
    else
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "no block found");
        #endif
        return 0;    
    }

}
/***   FindBlock   ***/





/***   SetNewBlockStart   *****************************************************
 *
 * begins to read from start and stops at end or at first non-0
 *
 */
void
SetNewBlockStart(int *rawshift, int *start, int end)
{
    int i=0;
    for (i=*start; i<=end; i++)
    {
        if (rawshift[i]!=0)
        {
            *start = i;
            return;
        }
    }
    *start = end;
}
/***   SetNewBlockStart   ***/





/***   MoveBlock   ************************************************************
 *
 * Moves a block inside rawshift via provided blockmove
 *
 */
void
MoveBlock(int *rawshift, block_move_t blockmove)
{
    int new_shift = 0.0;
    int   i=0;
    for (i=blockmove.start; i<=blockmove.end; i++)
    {
        new_shift = rawshift[i] - blockmove.offset;
        memcpy(&rawshift[i], &new_shift, sizeof(int));
    }
}
/***   MoveBlock   ***/





/***   ScoreFromRawShift   *******************************************************
 *
 * Computes the score from a rawshift
 * This is done by looking up a block, then moving it by its minimal elements
 * and trying to find a new one
 *
 * size is exclusive : rawshift(0 - size-1)
 *
 */
int
ScoreFromRawShift(int *rawshift, int size)
{
    int          start     = 0;
    int          end       = size-1;
    int          score     = 0;
    block_move_t blockmove;

    
    while (FindBlock(rawshift, start, end, &blockmove))
    {
        if (blockmove.offset==0)
            break;
        
        DEBUG_P("blockmove->start=%d\n", blockmove.start);
        DEBUG_P("blockmove->end=%d\n", blockmove.end);
        DEBUG_P("blockmove->offset=%d\n", blockmove.offset);

        
        MoveBlock(rawshift, blockmove);
        score += abs(blockmove.offset);


        SetNewBlockStart(rawshift, &start, end);
    }
    
    return score;
}
/***   ScoreFromRawShift   ***/



/***   CalcRawShift   ************************************************************
 *
 * Calculates the rawshift out of a predicted and a trusted shift_seq
 *
 */
int
CalcRawShift(shift_seq_t *trust, shift_seq_t *pred, int *pairmask)
{
    int     *rawshift  = NULL;   /* shift of each nal nt */
    int      score     = 0;      /* return */
    int     nal_idx   = 0;
    
    
    rawshift = (int *) calloc(pred->nal_len, sizeof(int));
    
    /* creates a rawshift:
     * the shift of each nt to its predicted postition
     */
    for (nal_idx=0; nal_idx < pred->nal_len; nal_idx++) /* or trust->nal_len */
    {
        /* paranoia
         */
        if (trust->residue[nal_idx].nt != pred->residue[nal_idx].nt)
        {
            ERROR_P("Uppss...NTs at unaligned position %d differ\n", nal_idx);
            ERROR_P("%s\n", "I will ignore the fact but you should check the sequences!");
        }
    
        rawshift[nal_idx]  =  trust->residue[nal_idx].aln_pos;
        rawshift[nal_idx] -=  pred->residue[nal_idx].aln_pos;
        
    }
    
    if (pairmask!=NULL)
        InsertPairmask(trust, pairmask, rawshift);
    
        
        
    /*** pretty print the shift
     *   
     */
    if (PRETTY_PRINT_SHIFT && cs_opts.verbose)
    {
        int j=0;
        char nt = 0;
        fprintf(stdout, "Raw Shift:\n");
        for (j=0; j<trust->nal_len; j++)
        {
            nt = pred->residue[j].nt;
            if (rawshift[j]==0)
            {
                fprintf(stdout, " %c=%d | ",  nt, rawshift[j]);
            }
            else if (rawshift[j]<0)
            {
                fprintf(stdout, " %d<%c | ", rawshift[j], nt);
            }
            else
            {
                fprintf(stdout, " %c>%d | ", nt, rawshift[j]);
            }
        }
        fprintf(stdout, "\n");
    }
    
    
    score = ScoreFromRawShift(rawshift, pred->nal_len);
    
    free(rawshift);
    
    return score;
}
/***   CalcRawShift   ***/






/***   SeqCmp   ***
 *
 */
int
SeqCmp(const void *seqt1, const void *seqt2)
{
    seq_t *this_seqt1, *this_seqt2;
    char *nt1, *nt2;
    int score;
        
    this_seqt1 = (seq_t *) seqt1;
    this_seqt2 = (seq_t *) seqt2;      
    
    nt1 = (char*) Scalloc(strlen(this_seqt1->nt)+1, sizeof(char));
    nt2 = (char*) Scalloc(strlen(this_seqt1->nt)+1, sizeof(char));

    Degap(this_seqt1->nt, nt1);
    Degap(this_seqt2->nt, nt2);

#ifdef DEBUG
    DEBUG_P("comparing %05s : %05s = %05s... : %05s...\n",
            this_seqt1->id, this_seqt2->id, nt1, nt2);
#endif
    
    score=strcmp(nt1, nt2);
    
    free(nt1);
    free(nt2);
    
    return score;
}
/***   SeqCmp   ***/
