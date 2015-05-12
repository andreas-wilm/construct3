/******************************************************************************
*
* move_update.c - routines for handling cs_dp update after a move
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
 *  CVS $Id: move_update.c,v 1.39 2007-10-22 10:43:23 steger Exp $
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>
#include <tk.h>
#include <math.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"


#include "consensus.h"
#include "if.h"
#include "proj.h"
#include "gui.h"

#include "move_update.h"



/***   private
 */
#if 0
    #define DEBUG
#endif

#define BEWARE_THE_FLOAT

#define RIGHT 0
#define DOWN  1
#define LEFT  2
#define UP    3
#define RIGHTDOWN  4
#define RIGHTUP    5
#define LEFTDOWN   6
#define LEFTUP     7


typedef struct {
    int seq_idx;
    int nt_low;
    int nt_high;
    int offset;
} move_t;


static void
Update_SeqNt(aln_t *ali, move_t move);

static int
UpdateBp(Tcl_Interp *interp, aln_t *aln, move_t move);

static void
Update_GapShift(aln_t *aln, move_t move);

static int
DrawUpdate_Bp(Tcl_Interp *interp, aln_t *ali,
              int seq_idx, bpair_t *bp, pair_t old_pair, pair_t new_pair);

static void
Update_ConsSeq(proj_t *proj, move_t move);

static void
Update_ConsNumGaps(Tcl_Interp *interp, move_t move);

static int
DrawUpdate_ConsGap(Tcl_Interp *interp, int nt_idx, int new_num_gaps, int old_num_gaps);

static int
AddToConsBp(Tcl_Interp *interp, int seq_idx, pair_t dest_pair, pair_t src_pair);

static int
RemoveFromConsBp(Tcl_Interp *interp, int seq_idx, pair_t pair);

static int
DrawUpdate_ConsBp(Tcl_Interp *interp, pair_t pair, float old_raw_prob);

static void
PrintContribs(aln_t *aln, pair_t pair);

/*
 ***/






/***   MoveNtSelection_Cmd   ***************************************************
 *
 * Handles Aligned Nt movement
 * And returns modfied sequence
 *
 */
int
MoveNtSelection_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    move_t   move;
    Stopwatch_t     *watch;
    Tcl_Obj *obj_poi  = NULL;


 
    watch = TimerStart();

    /*****   check and get args
     *
     */
    if (objc!=5) {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq_no nt_low nt_high offset?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &move.seq_idx) != TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[2], &move.nt_low)  != TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[3], &move.nt_high) != TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[4], &move.offset)  != TCL_OK)
        return TCL_ERROR;
    if (move.nt_low>move.nt_high) {
        int tmp;
        WARN_P("%d\n", "nt_low>nt_high: swapping");
        tmp = move.nt_low;
        move.nt_low = move.nt_high;
        move.nt_high = tmp;
    }
    #ifdef DEBUG
        DEBUG_P("%s\n", "start");
        DEBUG_P("move.seq_idx=%d\n", move.seq_idx);
        DEBUG_P("move.nt_low =%d\n", move.nt_low);
        DEBUG_P("move.nt_high=%d\n", move.nt_high);
        DEBUG_P("move.offset =%d\n", move.offset);
    #endif



    /* c is zero-offset */
    move.seq_idx--;
    move.nt_low--;
    move.nt_high--;


    /*** Notes:
     *
     *   - graphical update of alignment canvas is
     *     done exlusively by tcl/tk
    */
    /* end and start are inclusive */
    /* left */
    if (move.offset==0) {
        char buf[1024];
        sprintf(buf, "Illegal Move Offset %d\n", move.offset);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }



    UpdateBp(interp, aln, move);

    Update_SeqNt(aln, move);

    /* this must follow after all functions accessing GetBpFromAlnSeq
       since this relies on the gapshift */
    Update_GapShift(aln, move);

    Update_ConsSeq(proj, move);

    Update_ConsNumGaps(interp, move);


    /* return updated sequence to interpreter
     */
    obj_poi = Tcl_NewStringObj(aln->seq[move.seq_idx].nt, -1);
    Tcl_SetObjResult(interp, obj_poi);

    TimerStop(__FUNCTION__, watch);

    return TCL_OK;
}
/***   MoveNtSelection_Cmd   ***/




/***   UpdateBp   *************************************************************
 *
 * This is tricky stuff, but it works
 * iterate overall possible basepairs wher nucleotides in moved
 * sequence slice may interact
 * If that position really forms a basepair
 *  - move and retag it
 *  - remove from consensus
 *  - add to new consensus
 *
 * Some pitfalls:
 * - Iteration over all possible basepairs must be done in a way so
 *   that retagging works
 * - If a basepair resides _inside_ the sequence slice then the
 *   move it directly to new diagonal offset position when
 *   encountering it the first time and skip it afterwards
 *
 */
int
UpdateBp(Tcl_Interp *interp, aln_t *aln, move_t move)
{
    int i, j, k, skip;
    pair_t dest_pair, src_pair;
    bpair_t *bp;
    pair_t  *exclude;
    int n_excl;
    int start, end;

    /* in case a seq slice is moved which contains a basepair inside:
       move this basepair only once and over diagonal
       moving it twice doesn't work for retagging
       already moved basepairs are stored in pair_t *exclude
    */
    exclude = (pair_t *) Scalloc(move.nt_high-move.nt_low+1, sizeof(pair_t));
    n_excl=0;
    skip=0;

    /* distinguish between move order so that retagging in
       DrawUpdate_Bp does not overwrite tags
    */
    if (move.offset<0) {
        start = move.nt_low;
        end   = move.nt_high;
    } else {
        start = move.nt_high;
        end   = move.nt_low;
    }

    for (i=start; i!=end-move.offset; i=i-move.offset) {
        for (j=0; j<aln->len; j++) {
            if (j==i || i+move.offset==j)
                continue;

            if (IS_GAP(aln->seq[move.seq_idx].nt[i])
                ||
                IS_GAP(aln->seq[move.seq_idx].nt[j]))
                continue;

            SetPair(i+move.offset, j, &dest_pair);
            SetPair(i, j, &src_pair);

            bp = GetBpFromAlnSeq(aln, move.seq_idx, src_pair);
            if (bp==NULL)
                continue;


            /* find in exclude list
             */
            for (k=0; k<n_excl; k++) {
                if (exclude[k].nti==src_pair.nti
                    &&
                    exclude[k].ntj==src_pair.ntj) {
                    skip=1;
                    break;
                }
            }
            if (skip) {
                skip=0;
                continue;
            }

            /* set exclude list */
            skip=0;
            if (j>=move.nt_low && j<=move.nt_high) {
                n_excl++;
                SetPair(j, i, &exclude[n_excl-1]);

                SetPair(i+move.offset, j+move.offset, &dest_pair);
            }


            if (DrawUpdate_Bp(interp, aln, move.seq_idx, bp, dest_pair, src_pair)
                ==
                TCL_ERROR)
                return TCL_ERROR;

            /* next two also call DrawUpdate_ConsBp */
            if (RemoveFromConsBp(interp, move.seq_idx, src_pair)
                ==
                TCL_ERROR)
                return TCL_ERROR;

            if (AddToConsBp(interp, move.seq_idx, dest_pair, src_pair)
                ==
                TCL_ERROR)
                return TCL_ERROR;
        }
    }

    /*
    TMPDEBUG_P("n_excl=%d\n", n_excl);
    for (k=0; k<n_excl; k++) {
        TMPDEBUG_P("exclude[k=%d].nti=%d\n", k, exclude[k].nti);
        TMPDEBUG_P("exclude[k=%d].ntj=%d\n", k, exclude[k].ntj);
    }
    */
    
    free(exclude);


    return TCL_OK;
}
/***   UpdateBp   ***/







/***   Update_GapShift   ******************************************************
 *
 * Update the Gapshift
 * done by recomputing
 * relies on the fact, that only one nt is moved
 *
 * gapshift before move start cannot have changed
 * gapshift after move end cannot have changed
 *
 * nts must be already have been updated
 */
void
Update_GapShift(aln_t *aln, move_t move)
{
    int i, shift;
    char *seq;
    int start, end;

    /*
    TMPDEBUG_P("%s\n", "old gapshift");
    for (i=0; i<aln->len; i++)
        fprintf(stdout, "%3d ", aln->seq[move.seq_idx].gapshift[i]);
    fprintf(stdout, "\n");
    */

    if (move.offset<0) {
        start = move.nt_low-1;
        end   = move.nt_high;
    } else {
        start = move.nt_low;
        end   = move.nt_high+1;
    }

    seq = aln->seq[move.seq_idx].nt;

    shift=0;
    for (i=0; i<start; i++)
        if (IS_GAP(seq[i]))
            shift++;


    for (i=start; i<=end; i++) {
        if (IS_GAP(seq[i])) {
            shift++;
            aln->seq[move.seq_idx].gapshift[i] = 0;
        } else {
            aln->seq[move.seq_idx].gapshift[i] = shift;
        }
    }
    /*
    TMPDEBUG_P("%s\n", "new gapshift");
    for (i=0; i<aln->len; i++)
        fprintf(stdout, "%3d ", aln->seq[move.seq_idx].gapshift[i]);
    fprintf(stdout, "\n");
    */
}
/***   Update_GapShift   ***/



/***   Update_SeqNt   **********************************************************
 *
 *
 * Modifies alignment by moving an nt range beginning with nt_low to nt_high
 * to left (if neg. offset) or right (if pos. offset)
 * the resulting backlog will be transformed to a gap
 *
 */
void
Update_SeqNt(aln_t *ali, move_t move)
{
    int i;
    seq_t *seq;

    seq = & ali->seq[move.seq_idx];

    /*** update seq
     */
    /* move selected nts */
    if (move.offset<0) {
        /* nt */
        for (i=move.nt_low-1; i<move.nt_high; i++)
            seq->nt[i] = seq->nt[i+1];
        /* flipped gap */
        seq->nt[move.nt_high] = '-';
    } else {
        /* nt */
        for (i=move.nt_high+1; i>move.nt_low; i--)
            seq->nt[i] = seq->nt[i-1];
        /* flipped gap */
        seq->nt[move.nt_low] = '-';
    }
}
/***   Update_SeqNt   ***/




/***   DrawUpdate_Bp   ********************************************************
 *
 * Moves and retags one Bp
 *
 */
int
DrawUpdate_Bp(Tcl_Interp *interp, aln_t *ali,
              int seq_idx, bpair_t *bp, pair_t new_pair, pair_t old_pair)
{
    double       dp_size, dp_virsize;
    float        cell_height, cell_width;
    char         cmd[1024];    /* used for tcl/tk command formatting */
    Tcl_Obj     *cmd_obj;      /* used for executing tk commands     */
    char         buf[1024];
    int          direction;
    char         oldtag[1024], newtag[1024];

    if (GetCanvasVirtualsize(interp, &dp_virsize)==TCL_ERROR)
        return TCL_ERROR;
    if (GetCanvasDotplotSize(interp, &dp_size)==TCL_ERROR)
        return TCL_ERROR;

    cell_height = (float) (dp_size / ali->len);
    cell_width  = (float) (dp_size / ali->len);

    cell_height *= (dp_virsize/dp_size);
    cell_width  *= (dp_virsize/dp_size);


    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);

    GetThisDpBpTag(seq_idx, old_pair, oldtag);
    GetThisDpBpTag(seq_idx, new_pair, newtag);

    sprintf(cmd, "%s move %s", CANVAS_DP_NAME, oldtag);
    Tcl_SetStringObj(cmd_obj, cmd, -1);


    /*** what direction ?
      * remember: pair indices are always set in a way,
      * so that j>i, although this contradicts
      * to normal dotplot indices, where x=i and y=j
      * and therefore i>j
     */

    /* we are moving y */
    if (old_pair.nti==new_pair.nti) {
        if (new_pair.ntj > old_pair.ntj)
             direction=DOWN;
        else
            direction=UP;
    /* moving x */
    } else if (old_pair.ntj==new_pair.ntj) {
        if (new_pair.nti > old_pair.nti)
             direction=RIGHT;
        else
            direction=LEFT;
    /* moving x and y */
    } else {
        if (new_pair.nti > old_pair.nti) {
            if (new_pair.ntj > old_pair.ntj)
                 direction=RIGHTDOWN;
            else /* && (new_pair.nti < old_pair.nti))*/
                direction=RIGHTUP;
        } else /* (new_pair.ntj < old_pair.ntj) */ {
            if (new_pair.ntj > old_pair.ntj)
                direction=LEFTDOWN;
            else
                direction=LEFTUP;
        }
    }
    /*
    if (direction==RIGHT)
        TMPDEBUG_P("moving %s %s\n", oldtag, "RIGHT");
    else if (direction==DOWN)
        TMPDEBUG_P("moving %s %s\n", oldtag, "DOWN");
    else if (direction==LEFT)
        TMPDEBUG_P("moving %s %s\n", oldtag, "LEFT");
    else if (direction==UP)
        TMPDEBUG_P("moving %s %s\n", oldtag, "UP");
    else if (direction==RIGHTDOWN)
        TMPDEBUG_P("moving %s %s\n", oldtag, "RIGHTDOWN");
    else if (direction==RIGHTUP)
        TMPDEBUG_P("moving %s %s\n", oldtag, "RIGHTUP");
    else if (direction==LEFTDOWN)
        TMPDEBUG_P("moving %s %s\n", oldtag, "LEFTDOWN");
    else /@ (direction==LEFTUP) @/
        TMPDEBUG_P("moving %s %s\n", oldtag, "LEFTUP");
    */

    if (direction==RIGHT)
        sprintf(cmd, " %fc %fc", cell_width, 0.0);
    else if (direction==DOWN)
        sprintf(cmd, " %fc %fc", 0.0, cell_height);
    else if (direction==LEFT)
        sprintf(cmd, " %fc %fc", -cell_width, 0.0);
    else if (direction==UP)
        sprintf(cmd, " %fc %fc", 0.0, -cell_height);
    else if (direction==RIGHTDOWN)
        sprintf(cmd, " %fc %fc", cell_width, cell_height);
    else if (direction==RIGHTUP)
        sprintf(cmd, " %fc %fc", cell_width, -cell_height);
    else if (direction==LEFTDOWN)
        sprintf(cmd, " %fc %fc", -cell_width, cell_height);
    else /* (direction==LEFTUP) */
        sprintf(cmd, " %fc %fc", -cell_width, -cell_height);

    Tcl_AppendToObj(cmd_obj, cmd, -1);


    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR) {
        sprintf(buf, "Couldn't eval %s : %s\n",
              Tcl_GetStringFromObj(cmd_obj, NULL),
                      Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /*****   retag
     */
    sprintf(cmd, "ReplaceTag %s %s %s", CANVAS_DP_NAME, oldtag, newtag);

    Tcl_SetStringObj(cmd_obj, cmd, -1);
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR) {
        sprintf(buf, "Couldn't eval %s : %s\n",
              Tcl_GetStringFromObj(cmd_obj, NULL),
                      Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
    }

    Tcl_DecrRefCount(cmd_obj);


    return TCL_OK;
}
/***   DrawUpdate_Bp   ***/






/***   Update_ConsSeq   ********************************************************
 *
 * Sequence must be already updated by move (UpdateSeqNt)!
 */
void
Update_ConsSeq(proj_t *proj, move_t move)
{
    float *weight=NULL;
    int nt_low = 0; int nt_high =0;

    if (move.offset<0) {
        nt_low  =  move.nt_low  + move.offset;
        nt_high =  move.nt_high;
    } else {
        nt_low  =  move.nt_low;
        nt_high =  move.nt_high + move.offset;
    }

    #ifdef DEBUG
        DEBUG_P("Old Cons Seq = \"%s\"\n", aln->cons_seq->nt);
    #endif

    weight = GetWeightsAsArray();
    CalcConsSeq(aln, nt_low, nt_high, weight, aln->cons_seq->nt);

    DEBUG_P("New Cons Seq = \"%s\"\n", aln->cons_seq->nt);

    free(weight);
}
/***   Update_ConsSeq   ***/





/***   Update_ConsNumGaps   ****************************************************
 *
 * Sequence must be already updated by move (UpdateSeqNt)!
 *
 * NumGaps can only change at those indices (and their neighbours)
 * where a gap inside the block resides
 *
 */
void
Update_ConsNumGaps(Tcl_Interp *interp, move_t move)
{
    seq_t      *seq;
    cons_seq_t *cs_seq;
    char old_nt, new_nt;
    int  old_ngaps = 0;
    int  i;

    seq    = & aln->seq[move.seq_idx];
    cs_seq = aln->cons_seq;

    /*****   Move to left
     *
     */
    if (move.offset<0) {
        /* Leader always looses a gap
         */
        old_ngaps = cs_seq->num_gaps[move.nt_low-1];
        aln->cons_seq->num_gaps[move.nt_low-1]--;
        DrawUpdate_ConsGap(interp, move.nt_low-1, cs_seq->num_gaps[move.nt_low-1], old_ngaps);

        for (i=move.nt_low; i<move.nt_high; i++) {
            old_nt = seq->nt[i-1];
            new_nt = seq->nt[i];
            old_ngaps = cs_seq->num_gaps[i];

            if (IS_GAP(old_nt)==IS_GAP(new_nt)) {
                continue;
            } else if ( (! IS_GAP(old_nt)) && (IS_GAP(new_nt))) {
                cs_seq->num_gaps[i]++;
            } else if ( (IS_GAP(old_nt)) && (! IS_GAP(new_nt))) {
                cs_seq->num_gaps[i]--;
            }

            DrawUpdate_ConsGap(interp, i, cs_seq->num_gaps[i], old_ngaps);
        }

        /* Trailing end always gains a gap
         */
        old_ngaps = cs_seq->num_gaps[i];
        cs_seq->num_gaps[i]++;
        DrawUpdate_ConsGap(interp, i, cs_seq->num_gaps[i], old_ngaps);

        RaiseGapLayer(interp, move.nt_low, move.nt_high+1);

    /*****   Move to right
     *
     */
    } else {
        /* Leader always gains a gap
         */
        old_ngaps = cs_seq->num_gaps[move.nt_high+1];
        cs_seq->num_gaps[move.nt_high+1]--;
        DrawUpdate_ConsGap(interp, move.nt_high+1, cs_seq->num_gaps[move.nt_high+1], old_ngaps);

        for (i=move.nt_high; i>move.nt_low; i--) {
            old_nt = seq->nt[i+1];
            new_nt = seq->nt[i];
            old_ngaps = cs_seq->num_gaps[i];

            if (IS_GAP(old_nt)==IS_GAP(new_nt)) {
                continue;
            } else if ( (! IS_GAP(old_nt)) && (IS_GAP(new_nt))) {
                cs_seq->num_gaps[i]++;
            } else if ( (IS_GAP(old_nt)) && (! IS_GAP(new_nt))) {
                cs_seq->num_gaps[i]--;
            }

            DrawUpdate_ConsGap(interp, i, cs_seq->num_gaps[i], old_ngaps);
        }

        /* Trailing end always gains a gap
         */
        old_ngaps = cs_seq->num_gaps[i];
        cs_seq->num_gaps[i]++;
        DrawUpdate_ConsGap(interp, i, cs_seq->num_gaps[i], old_ngaps);

        RaiseGapLayer(interp, move.nt_low+1, move.nt_high+2);
    }

}
/***   Update_ConsNumGaps   ***/





/***   DrawUpdate_ConsGap   ****************************************************
 *
 * n_gaps must have changed
 */
int
DrawUpdate_ConsGap(Tcl_Interp *interp, int nt_idx, int new_num_gaps, int old_num_gaps)
{
    Tcl_Obj    *cmd_obj;      /* used for executing tk commands     */
    char        cmd[1024];    /* used for tcl/tk command formatting */
    char        coords[1024];
    unsigned int red   = 0;
    unsigned int green = 0;
    unsigned int blue  = 0;


    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);


    /*****  Update coords and color (if needed)
     *
     */
    if (new_num_gaps == 0) {
        /* convert to dummy cons gap
         */
        sprintf(cmd, "RecoordConsGap %d", nt_idx+1);
        Tcl_SetStringObj(cmd_obj, cmd, -1);
    } else {
        red   = (1.0 - ( 0.5 * (double)new_num_gaps / (double)aln->num_seq )) * 255;
        green = red;
        blue  = 255;


        /* recoord
         */
        GetGapCoords(interp, nt_idx+1, coords);
        sprintf(cmd, "RecoordConsGap %d [list %s]", nt_idx+1, coords);
        Tcl_SetStringObj(cmd_obj, cmd, -1);

        #ifdef DEBUG
            DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
        #endif
        if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR) {
            char buf[1024];
            sprintf(buf, "Couldn't eval %s : %s\n",
                   Tcl_GetStringFromObj(cmd_obj, NULL),
                           Tcl_GetStringResult(interp));
            CsDpError(interp, buf);
            return TCL_ERROR;
        }


        /* recolor
         */
        sprintf(cmd, "RecolorConsGap %d #%06x", nt_idx+1, (red<<16) | (green<<8) | blue);
        Tcl_SetStringObj(cmd_obj, cmd, -1);
    }

    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR) {
        char buf[1024];
            sprintf(buf, "Couldn't eval %s : %s\n",
                   Tcl_GetStringFromObj(cmd_obj, NULL),
                           Tcl_GetStringResult(interp));
            CsDpError(interp, buf);
            return TCL_ERROR;
    }

    Tcl_DecrRefCount(cmd_obj);

    return TCL_OK;
}
/***   DrawUpdate_ConsGap   ***/





/***   RemoveFromConsBp   *********************************************************
 *
 * Remove from consensus if it is a contrib
 */
int
RemoveFromConsBp(Tcl_Interp *interp, int seq_idx, pair_t pair)
{
    cons_bpair_t *cs_bp;
    bpair_t      *contrib_bp;
    float            old_prob   = 0.0;
    dlist_elem      *dlist_elem;
    void          *freeme = NULL;
    int *contrib_seqidx;


    cs_bp = GetConsBp(aln, pair);


    /* requested bp participates in contrib? */
    if (dlist_size(cs_bp->contribs)==0)
        return TCL_OK;
    dlist_elem = dlist_head(cs_bp->contribs);
    while (1) {
        contrib_seqidx = (int*)dlist_data(dlist_elem);
        if (*contrib_seqidx==seq_idx)
            break;

        if (dlist_is_tail(dlist_elem))
            break;
        dlist_elem = dlist_next(dlist_elem);
    }
    if (*contrib_seqidx!=seq_idx)
        return TCL_OK;

    /*
    TMPDEBUG_P("%s\n", "contribs before removing");
    PrintContribs(aln, pair);
    */

    contrib_bp = GetBpFromAlnSeq(aln, seq_idx, pair);
    if (contrib_bp==NULL) {
        WARN_P("Uupss..the contrib (seq %d) I wanted to delete is NULL\n", seq_idx);
        return TCL_OK;
    }
    if (! contrib_bp->prob>0.0) {
        WARN_P("Uupss..the contrib (seq %d) I wanted to delete has prob 0\n", seq_idx);
        return TCL_OK;
    }


    /***   Update prob and the contribs
     *
     */
    old_prob     = cs_bp->prob;
    cs_bp->prob -= pow(contrib_bp->prob, 1.0/(double)POW_A) * proj->entry[seq_idx].weight;


    #ifdef BEWARE_THE_FLOAT
        /* rule out rounding problems */
        /* FIXME:  use 10exp-FLT_DIG * dlist_size ? */
        if (cs_bp->prob<-0.000001)
            cs_bp->prob=0.0;
    #else
        if (cs_bp->prob<-0.000001) {
            char buf[1024];
            sprintf(buf, "updated_prob for cs_bp_idx %d (=%f) < 0 ! (old = %f)\n",
                                                       bp_idx, cs_bp->prob, old_prob);
            Die(buf);
        } else if (cs_bp->prob==-0.000001) {
            fprintf(stdout, "FIXME(%s:%s): cs_bp->prob==-0.000001 setting to 0.0\n", __FILE__, __FUNCTION__);
            /* FIXME: should we treat this floating point precision error as an severe error ? */
            cs_bp->prob=0.0;
        }
    #endif

    /* remove the contrib_bp
     */
    if (dlist_remove(cs_bp->contribs, dlist_elem, &freeme) == DLIST_ERROR) {
        char buf[1024];
        sprintf(buf, "couldn't dlist_remove(cs_bp->contribs, dlist_elem, &freeme)\n");
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /*
    TMPDEBUG_P("%s\n", "contribs after removing");
    PrintContribs(aln, pair);
    */

    DrawUpdate_ConsBp(interp, pair, old_prob);

    free(freeme);

    return TCL_OK;
}
/***   RemoveFromConsBp   ***/




/***   AddToConsBp   *********************************************************
 *
 *
 */
int
AddToConsBp(Tcl_Interp *interp, int seq_idx, pair_t dest_pair, pair_t src_pair)
{
    cons_bpair_t *cs_bp;
    bpair_t      *contrib_bp;
    float         old_prob;
    int          *contrib_seqidx;

    contrib_bp = GetBpFromAlnSeq(aln, seq_idx, src_pair);
    if (contrib_bp==NULL)
        return TCL_OK;
    if (! contrib_bp->prob>0.0)
        return TCL_OK;
    /*
    TMPDEBUG_P("%s\n", "contribs before adding");
    PrintContribs(aln, dest_pair);
    */

    cs_bp = GetConsBp(aln, dest_pair);

    /***   Update prob and the contribs
     *
     */
    old_prob = cs_bp->prob;
    cs_bp->prob += pow(contrib_bp->prob, 1.0/(double)POW_A) * proj->entry[seq_idx].weight;

    contrib_seqidx = Scalloc(1, sizeof(int));
    *contrib_seqidx = seq_idx;


    if (dlist_append(cs_bp->contribs, contrib_seqidx) == DLIST_ERROR) {
        char buf[1024];
        sprintf(buf, "couldn't dlist_append(&cs_bp->contribs, contrib_seqidx)\n");
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /*
    TMPDEBUG_P("%s\n", "contribs after adding");
    PrintContribs(aln, dest_pair)
    */

    DrawUpdate_ConsBp(interp, dest_pair, old_prob);


    return TCL_OK;
}
/***   AddToConsBp   ***/




/***   DrawUpdate_ConsBp   ****************************************************
 *
 * Decides whether a new dpbp has to be drawn, and old one has to bee deleted
 * or resized
 *
 * old_raw_prob may not be powered by POW_B
 *
 */
int
DrawUpdate_ConsBp(Tcl_Interp *interp, pair_t pair, float old_raw_prob)
{
    cons_bpair_t *cs_bp;
    float  old_prob, new_prob;
    double dp_scale;
    float  norm_factor;

    norm_factor = GetSumOfWeights();
    cs_bp       = GetConsBp(aln, pair);

    old_prob = pow(old_raw_prob / norm_factor, (double)POW_B);
    new_prob = pow(cs_bp->prob  / norm_factor, (double)POW_B);


    #ifdef DEBUG
        DEBUG_P("Updating %d:%d %f->%f\n", pair.nti, pair.ntj, old_prob, new_prob);
    #endif
    if (GetDpScale(interp, &dp_scale)==TCL_ERROR)
        return TCL_ERROR;

    if (old_prob<CS_BP_PROB_DRAW_CUTOFF) {
        /* old and new prob are below cutoff means
           we don't have to draw anything
         */
        if (new_prob<CS_BP_PROB_DRAW_CUTOFF)
            return TCL_OK;

        /* old prob is below cutoff but new one not means
           we have to create a new dp item
         */
        if (NewConsDpBp(interp, pair, dp_scale)==TCL_ERROR)
            return TCL_ERROR;
    } else {
        /* if old prob is above cutoff but new is below
           means we have to delete the old one
         */
        if (new_prob<CS_BP_PROB_DRAW_CUTOFF) {
            if (DeleteConsDpBp(interp, pair)==TCL_ERROR)
                return TCL_ERROR;
        } else {
            if (ReColorAndResizeConsDpBp(interp, pair, dp_scale)==TCL_ERROR)
                return TCL_ERROR;
        }
    }

    return TCL_OK;
}
/***   DrawUpdate_ConsBp   ***/




/***   PrintContribs   *******************************************************
 *
 * call before gapshift update
 */
void
PrintContribs(aln_t *aln, pair_t pair)
{
    cons_bpair_t *cs_bp;
    dlist_elem      *dlist_elem;
    int *contrib_seqidx;

    TMPDEBUG_P("contribs for csbp %d:%d\n", pair.nti, pair.ntj);

    cs_bp = GetConsBp(aln, pair);

    if (dlist_size(cs_bp->contribs)==0) {
        TMPDEBUG_P("   %s\n", "none");
        return;
    }
    dlist_elem = dlist_head(cs_bp->contribs);
    while (1) {
        contrib_seqidx = (int*)dlist_data(dlist_elem);
        /* GetBpFromAlnSeq not allowed here, since gapshift not updated */

        TMPDEBUG_P("   seq %d\n", *contrib_seqidx);

        if (dlist_is_tail(dlist_elem))
            break;
        dlist_elem = dlist_next(dlist_elem);
    }
}
/***   PrintContribs   ***/
