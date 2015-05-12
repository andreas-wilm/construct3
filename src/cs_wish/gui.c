/******************************************************************************
*
* gui.c - enhanced tcl routines for the cs_dp gui
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
 *  CVS $Id: gui.c,v 1.51 2007-10-22 10:43:23 steger Exp $
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

#include "mic.h"
#include "if.h"
#include "consensus.h"

#include "gui.h"



/***   private
 */

#if 0
    #define DEBUG
#endif

static void
ConsDpBpColor(cons_bpair_t *cs_bp, int *r, int *g, int *b);

/*
 ***/







/***   SetupGui   *************************************************************
 *
 *
 *
 */
int
SetupGui(Tcl_Interp *interp)
{
    char cmd[8096];
    char buf[8096];


    /***   CleanGui
     */
    sprintf(cmd, "CleanGui");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /***   CreateAlWindow
     */
    sprintf(cmd, "CreateAlWindow %d %d", aln->num_seq, aln->len);
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }



    /***   CreateToolButtons
     */
    sprintf(cmd, "CreateToolButtons");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }



    /***   CreateAlnNtsAndLabels
     */
    sprintf(cmd, "CreateAlnNtsAndLabels");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /***   CreateLayers
     */
    sprintf(cmd, "CreateLayers");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }



    /*** Draw Basepairs
     */
#ifndef BATCH_NON_DP_GUI
    CreateDpBps(interp);
#endif


	/*** Use Mapping if requested
	 */
    sprintf(cmd, "DotplotMappingFrontend");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /***   Draw Gaps
     */
#ifndef BATCH_NON_DP_GUI
    CreateConsGaps(interp);
#endif


    /***   Draw Consensus Basepairs
     */
#ifndef BATCH_NON_DP_GUI
    CreateDpConsBps(interp);
#endif


    /***   CreateDpSeqLabels
     */
    sprintf(cmd, "CreateDpSeqLabels %s", aln->cons_seq->nt);
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /***   GenDpNumLabels
     */
    if (Tcl_EvalEx(interp, "GenDpNumLabels", -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"GenDpNumLabels\": %s", Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /***   CreateBindings
     */
    sprintf(cmd, "BindAlignment %d", aln->num_seq);
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    sprintf(cmd, "BindDotplot");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

#ifndef BATCH_NON_DP_GUI

    /* misc startup stuff
     */

    sprintf(cmd, "SeqSelectionChange 1");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /* show hide defaults
     */

    sprintf(cmd, "DpShowHide DpConsBp $opts(displ,cons_bp)");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    sprintf(cmd, "DpShowHide DpConsGap $opts(displ,gaps)");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    sprintf(cmd, "DpShowHide DpBp $opts(displ,bps)");
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", cmd);
    #endif
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
    {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
#endif



    return TCL_OK;
}
/***   SetupGui   ***/





/***   CreateDpBps   **********************************************************
 *
 * Creates basepairs in dotplot
 *
 */
int
CreateDpBps(Tcl_Interp *interp)
{
    int           seq_idx      = 0;       /* sequence index (zero offset)       */
    bpair_t      *bp_pointer   = NULL;
    double        dp_scale;
    Tcl_Obj      *cmd_obj      = NULL;    /* used for executing tk commands     */
    char          cmd[1024];       /* used for tcl/tk command formatting */
    int           nt_i, nt_j;
    rect_coord_t  coord;
    Stopwatch_t  *watch;
    pair_t        pair;

    watch = TimerStart();

    if (GetDpScale(interp, &dp_scale)==TCL_ERROR)
        return TCL_ERROR;


    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);


    /*** foreach sequence
     */
    for (seq_idx=0; seq_idx < aln->num_seq; seq_idx++)
    {
        DEBUG_P("Setting up basepair plot for %s (%d/%d)\n",
                          aln->seq[seq_idx].id, seq_idx+1, aln->num_seq);
        VERBOSE_P("Setting up basepair plot for %s (%d/%d)\n",
                          aln->seq[seq_idx].id, seq_idx+1, aln->num_seq);


        /*** foreach basepair
         */
        for (nt_j=0; nt_j<aln->len; nt_j++)
        {
            for (nt_i=nt_j+1; nt_i<aln->len; nt_i++)
            {
                SetPair(nt_i, nt_j, &pair);

                bp_pointer = GetBpFromAlnSeq(aln, seq_idx, pair);
                if (bp_pointer==NULL)
                    continue;
                if (bp_pointer->prob==0.0)
                    continue;

                GetProbCellCoords(bp_pointer->prob, pair, dp_scale, &coord);

                sprintf(cmd, "CreateNewDpBp %d %d %d %fc %fc %fc %fc",
                                            seq_idx+1, nt_i+1, nt_j+1,
                                                   coord.x1, coord.y1,
                                                   coord.x2, coord.y2);
                Tcl_SetStringObj(cmd_obj, cmd, -1);



                #ifdef DEBUG
                    DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
                #endif
                if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
                    return TCL_ERROR;
                }
        }
        /* end foreach basepair */

    } /* foreach sequence */



    /*** Raise the DpBps to appropiate position/layer
     * FIXME: necessary ?
     */
    Tcl_SetStringObj(cmd_obj, "RaiseDpBpsToLayer", -1);
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
        return TCL_ERROR;


    Tcl_DecrRefCount(cmd_obj);

    TimerStop(__FUNCTION__, watch);

    return TCL_OK;
}
/***   CreateDpBps   ***/





/***   CreateConsGaps   *********************************************************
 *
 *
 */
int
CreateConsGaps(Tcl_Interp *interp)
{
    Tcl_Obj    *cmd_obj     = NULL;    /* used for executing tk commands     */
    char        cmd[1024];      /* used for tcl/tk command formatting */
    char        coords[1024];
    unsigned int red   = 0;
    unsigned int green = 0;
    unsigned int blue  = 0;
    int nt_idx = 0;
    Stopwatch_t     *watch;


    watch = TimerStart();


    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);

    for (nt_idx=0; nt_idx < aln->len; nt_idx++)
    {
        int num_gaps = aln->cons_seq->num_gaps[nt_idx];


        /*** no gaps at all ?
         *   create a dummy
         *   so that a move update just needs to
         *   resize and move instead of deleting and creating a new stripe
         *  FIXME: SEVERE: how does a resize find this one (no tag!)
         */
        if (num_gaps == 0)
        {
            sprintf(cmd, "CreateDummyConsGap %d", nt_idx+1);
            Tcl_SetStringObj(cmd_obj, cmd, -1);
        }
        else
        {
            red   = (1.0 - ( 0.5 * (double)num_gaps / (double)aln->num_seq )) * 255;
            green = red;
            blue  = 255;

            GetGapCoords(interp, nt_idx+1, coords);

            sprintf(cmd, "CreateConsGap %d #%06x [list %s]", nt_idx+1,
                                        (red<<16) | (green<<8) | blue,
                                                              coords);
            Tcl_SetStringObj(cmd_obj, cmd, -1);
        }

        #ifdef DEBUG
            DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
        #endif

        if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
        {
            char buf[1024];
            sprintf(buf, "Couldn't eval %s : %s\n",
                   Tcl_GetStringFromObj(cmd_obj, NULL),
                           Tcl_GetStringResult(interp));
            CsDpError(interp, buf);
            return TCL_ERROR;
        }
    }


    /*** raise gaps sorted
     */
    RaiseGapLayer(interp, 1, aln->len);



    Tcl_DecrRefCount(cmd_obj);

    TimerStop(__FUNCTION__, watch);

    return TCL_OK;
}
/***   CreateConsGaps   ***/




/***   RaiseGapLayer   ********************************************************
 *
 * Raise all Gaps sorted to their priority
 *
 * min_nt / max_nt are unit-offset !
 */
int
RaiseGapLayer(Tcl_Interp *interp, int min_nt, int max_nt)
{
    Tcl_Obj    *cmd_obj;    /* used for executing tk commands     */
    char        cmd[1024];      /* used for tcl/tk command formatting */
    int i=0;



    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);



    for (i=min_nt; i<=max_nt; i++)
    {
        int num_gaps = aln->cons_seq->num_gaps[i-1];
        sprintf(cmd, "%s raise DpConsGap_nt_%d LAYER_ConsGaps_%d", CANVAS_DP_NAME, i, num_gaps);
        Tcl_SetStringObj(cmd_obj, cmd, -1);

        #ifdef DEBUG
            DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
        #endif

        if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
        {
            char buf[1024];
            sprintf(buf, "Couldn't eval %s : %s\n",
                   Tcl_GetStringFromObj(cmd_obj, NULL),
                          Tcl_GetStringResult(interp));
            CsDpError(interp, buf);
            return TCL_ERROR;
        }

    }

    Tcl_DecrRefCount(cmd_obj);

    return TCL_OK;
}
/***   RaiseGapLayer   ***/





/***   CreateDpConsBps   ******************************************************
 *
 */
int
CreateDpConsBps(Tcl_Interp *interp)
{
    int      nt_i, nt_j;
    double   dp_scale;
    char     buf[1024];
    pair_t   pair;
    Stopwatch_t     *watch;

    watch = TimerStart();

    if (GetDpScale(interp, &dp_scale)==TCL_ERROR)
        return TCL_ERROR;

    /*** foreach basepair
     */
    for (nt_i=0; nt_i<aln->len; nt_i++)
    {
        for (nt_j=nt_i+1; nt_j<aln->len; nt_j++)
        {
            SetPair(nt_i, nt_j, &pair);
            if (NewConsDpBp(interp, pair, dp_scale)==TCL_ERROR)
            {
                sprintf(buf, "NewConsDpBp failed (%s)\n", Tcl_GetStringResult(interp));
                CsDpError(interp, buf);
                return TCL_ERROR;
            }
        }
    }

    TimerStop(__FUNCTION__, watch);

    return TCL_OK;
}
/***   CreateDpConsBps   ***/





/***   NewConsDpBp   **********************************************************
 *
 * Creates a new consensus basepair in dotplot
 *
 */
int
NewConsDpBp(Tcl_Interp *interp, pair_t pair, float dp_scale)
{
    cons_bpair_t  *cs_bp;
    float          cs_prob;
    rect_coord_t   coords;
    char     cmd[1024];
    Tcl_Obj *cmd_obj;
    int red = 0; int green = 0; int blue = 0;


    cs_bp = GetConsBp(aln, pair);

    cs_prob = GetConsBpProb(aln, pair, CS_BP_PROB_DRAW_CUTOFF);

    if (cs_prob==0.0)
        return TCL_OK;

    #ifdef DEBUG
        DEBUG_P("nt_i=%d nt_j=%d cs_prob=%f\n",  pair.nti, pair.ntj,  cs_prob);
    #endif

    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);


    GetProbCellCoords(cs_prob, pair, dp_scale, &coords);

    ConsDpBpColor(cs_bp, &red, &green, &blue);

    sprintf(cmd, "CreateNewConsDpBp #%04x%04x%04x %d %d %fc %fc %fc %fc",
                                                       red, green,  blue,
                                                        pair.nti+1, pair.ntj+1,
                              coords.x1, coords.y1, coords.x2, coords.y2);
    Tcl_SetStringObj(cmd_obj, cmd, -1);

    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
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
/***   NewConsBp   ***/




/***   DeleteConsDpBp   **********************************************************
 *
 *
 * Deletes a cons_bp from dotplot
 * nti, ntj: zero offset nt indices
 */
int
DeleteConsDpBp(Tcl_Interp *interp, pair_t pair)
{
    Tcl_Obj *cmd_obj;     /* used for executing tk commands     */
    char     cmd[1024];  /* used for tcl/tk command formatting */

    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);


    sprintf(cmd, "DeleteConsBp %d %d", pair.nti+1, pair.ntj+1);
    Tcl_SetStringObj(cmd_obj, cmd, -1);

    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
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
/***   DeleteConsDpBp   ***/





/***   ReColorAndResizeConsDpBp   *********************************************
 *
 * Called if a cons bp in dp must me updated
 * e.g. prob has changed but is above cutoff
 *
 */
int
ReColorAndResizeConsDpBp(Tcl_Interp *interp, pair_t pair, float dp_scale)
{
    cons_bpair_t *cs_bp = NULL;
    float            cs_prob;
    rect_coord_t     coords;
    int red, green, blue;
    char     cmd[1024];      /* used for tcl/tk command formatting */
    Tcl_Obj *cmd_obj;        /* used for executing tk commands     */


    cs_bp = GetConsBp(aln, pair);

    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);


    /*** Resize
     */
    cs_prob = GetConsBpProb(aln, pair, CS_BP_PROB_DRAW_CUTOFF);
    GetProbCellCoords(cs_prob, pair, dp_scale, &coords);

    sprintf(cmd, "ResizeConsDp %d %d %fc %fc %fc %fc",
                                       pair.nti+1, pair.ntj+1,
                                 coords.x1, coords.y1,
                                 coords.x2, coords.y2);
    Tcl_SetStringObj(cmd_obj, cmd, -1);

    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
            char buf[1024];
            sprintf(buf, "Couldn't eval %s : %s\n",
                   Tcl_GetStringFromObj(cmd_obj, NULL),
                          Tcl_GetStringResult(interp));
            CsDpError(interp, buf);
            return TCL_ERROR;
    }


    /*** Recolor
     *
     */
    ConsDpBpColor(cs_bp, &red, &green, &blue);
    sprintf(cmd, "RecolorConsDp %d %d #%04x%04x%04x", pair.nti+1, pair.ntj+1,
                                                    red, green, blue);
    Tcl_SetStringObj(cmd_obj, cmd, -1);


    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
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
/***   ReColorAndResizeConsDpBp   ***/




/***   ConsDpBpColor   **********************************************************
 *
 * Sets RGB Value for one consensus-basepair
 * according to frequency of its contribs
 */
void
ConsDpBpColor(cons_bpair_t *cs_bp, int *r, int *g, int *b)
{
    double freq  = 0;    /* n_contribs / nseq  */

     /*** Setup CsDpColor
     *
     */
    freq = (double)dlist_size(cs_bp->contribs) / (double)aln->num_seq;

    if ( freq<=0.5)
    {
        *r = (int)(  86376*freq + 22347);
        *g = (int)(  39584*freq + 45743);
        *b = (int)(-131070*freq + 65535);
    }
    else
    {
        *r = 0xffff;
        *g = (int)(-131070*freq + 131070);
        *b = 0x0000;
    }
}
/***   ConsDpBpColor   ***/






/***   CreateInfoDp_Cmd   *****************************************************
 *
 */
int
CreateInfoDp_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    char     rgb[1024];
    char    *info_pair_tag; /* argument */
    char     cmd[1024];
    Tcl_Obj *cmd_obj;
    double   dp_scale  = 0.0;
    int      i=0, j=0;
    rect_coord_t coords;
    pair_t pair;
    float mic;

    double mi_colormap_0 = 0.0; /* cs2: InfoStat(Colormap,0) */
    double mi_colormap_1 = 1.0; /* cs2: InfoStat(Colormap,1)  */

    Stopwatch_t     *watch;


    /***   parse args
     *
     */
    if ( (objc!=4))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?info_pair_tag colormap_0 colormap_1?");
        return TCL_ERROR;
    }
    info_pair_tag = Tcl_GetStringFromObj(objv[1], NULL);
    if (info_pair_tag == NULL )
    {
        CsDpError(interp, "can´t get info_pair_tag");
        return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, objv[2], &mi_colormap_0)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &mi_colormap_1)==TCL_ERROR)
        return TCL_ERROR;



    /* check if info has been computed before
     * if not, throw an error and exit
     */
    if ( ! MIC_IsComputed())
    {
        CsDpError(interp, "No mutual information content found! Compute it before calling me.");
        return TCL_ERROR;
    }


    if (GetDpScale(interp, &dp_scale)==TCL_ERROR)
        return TCL_ERROR;

    watch = TimerStart();


    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);


#ifndef BATCH_NON_DP_GUI

    /*** Draw each info pair in dotplot
     */
    for (j=0; j<aln->len; j++)
    {
        for (i=0; i<j; i++)
        {
            SetPair(i, j, &pair);
            mic  = GetMIC(pair);

            sprintf(cmd, "%s create rectangle", CANVAS_DP_NAME);
            Tcl_SetStringObj(cmd_obj, cmd, -1);

            GetFullCellCoords(pair, dp_scale, &coords);
            sprintf(cmd, " %fc %fc %fc %fc", coords.x1, coords.y1,
                                             coords.x2, coords.y2);
            Tcl_AppendToObj(cmd_obj, cmd, -1);


            /* since the colormaps min and max ranges from 0 to 1
               we have to transform mic to a relative value
               before computing the rgb */
            GetRgb(rgb, mic/GetMaxMIC(), (float)mi_colormap_0,
                                         (float)mi_colormap_1);
            sprintf(cmd, " -fill %s", rgb);
            Tcl_AppendToObj(cmd_obj, cmd, -1);
#ifdef DEBUG
                DEBUG_P("%d:%d mic=%f (max=%f) rgb=%s\n", j+1, i+1, mic, GetMaxMIC(), rgb);
#endif

            sprintf(cmd, " -outline \"\" -tags \"%s\"", info_pair_tag);
            Tcl_AppendToObj(cmd_obj, cmd, -1);
#ifdef DEBUG
                DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
#endif
            if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
            {
                char buf[1024];
                sprintf(buf, "Couldn't eval %s : %s\n",
                       Tcl_GetStringFromObj(cmd_obj, NULL),
                              Tcl_GetStringResult(interp));
                CsDpError(interp, buf);
                return TCL_ERROR;
            }
        }
    }
    /*** end: draw each info pair */
#endif


    TimerStop(__FUNCTION__, watch);

    Tcl_DecrRefCount(cmd_obj);

    return TCL_OK;
}
/***   CreateInfoDp_Cmd   ***/




/***   GetGapCoords   *********************************************************
 *
 * writes gap coords for dp to coords on success
 * and an error message otherwise
 * nt_pos is unit-offset !
 *
 *
 *    seven points describe one gap:
 *    only four values must be computed since the
 *    coordinates can be swapped i.e.:
 *    1 = x1, y1
 *    5 = y1, x1
 *
 *        1 2
 *        | |
 *        7 |
 *         \3---4
 *          6---5
 *
 */
void
GetGapCoords(Tcl_Interp *interp, int nt_pos, char *coords)
{
    double dp_size, dp_virsize, dp_scale;
    double gapx1, gapy1, gapx2, gapy2;

    if (GetCanvasDotplotSize(interp, &dp_size)==TCL_ERROR)
    {
        sprintf(coords, "ERROR_P(%s) :GetCanvasDotplotSize failed", __FUNCTION__);
        return;
    }
    if (GetDpScale(interp, &dp_scale)==TCL_ERROR)
    {
        sprintf(coords, "ERROR_P(%s) :GetDpScale failed", __FUNCTION__);
        return;
    }
    if (GetCanvasVirtualsize(interp, &dp_virsize)==TCL_ERROR)
    {
        sprintf(coords, "ERROR_P(%s) :GetCanvasVirtualsize failed", __FUNCTION__);
        return;
    }

    gapx1  = (nt_pos-1)*dp_scale;
    gapy1  = 0.0;
    gapx2  = gapx1 + dp_scale;
    gapy2  = gapy1 + dp_virsize;

    sprintf(coords, "%fc %fc %fc %fc %fc %fc %fc %fc %fc %fc %fc %fc %fc %fc",
                                     gapx1, gapy1, gapx1, gapx1, gapx2, gapx2,
                                     gapy2, gapx2, gapy2, gapx1, gapx2, gapx1,
                                                                 gapx2, gapy1);
}
/***   GetGapCoords   ***/





/***   GetRgb   ***************************************************************
 *
 *  Writes RGB-Value to already allocated string <rgb>
 *  Argument 0.0<=hue<=1; returns RGB-string :
 *  0= yello -> green -> cyan -> blue -> magenta -> red =1
 *
 */
void
GetRgb(char *rgb, float hue, float lowlimit, float highlimit)
{
	int   i, v;
	float f, q, t;

	v         = 65535;

	if (lowlimit < 0.0)
		lowlimit = 0.0;
	if (highlimit > 1.0)
		highlimit = 1.0;


	if (hue < lowlimit)
    {
		sprintf(rgb, "#ffffffffffff");
		return;
	}

	if (hue > highlimit)
    {
		sprintf(rgb, "#000000000000");
		return;
	}

	hue = (hue - lowlimit)/(highlimit - lowlimit);

	hue = hue*5.0;
	if (hue >= 5.0)
		 hue = 4.999999;


	i = (int)hue;
	f = hue-i;
	q = 65535 * (1 - f);
	t = 65535 * f;

	switch(i)
    {
		case 0:
			sprintf(rgb, "#%04x%04x%04x", (int)q, v, 0);
			break;
		case 1:
			sprintf(rgb, "#%04x%04x%04x", 0 ,v ,(int)t);
			break;
		case 2:
			sprintf(rgb, "#%04x%04x%04x", 0, (int)q, v);
			break;
		case 3:
			sprintf(rgb, "#%04x%04x%04x", (int)t, 0, v);
			break;
		case 4:
			sprintf(rgb, "#%04x%04x%04x", v, 0, (int)q);
			break;
	}
}
/***   GetRgb   ***/




/***   GetProbCellCoords   *******************************************************
 *
 * Write Coords for a prob dependent basepair cell to
 * already allocated rect_coord_t *c
 *
 */
void
GetProbCellCoords(float prob, pair_t pair, double dp_scale, rect_coord_t *c)
{
    float mid_width, mid_height;
    float prob_area;


    prob_area  = sqrt(prob)/2.0;

    mid_width  = (pair.nti+0.5)*dp_scale;
    mid_height = (pair.ntj+0.5)*dp_scale;

    c->x1 = mid_width  - (dp_scale  * prob_area);
    c->y1 = mid_height - (dp_scale  * prob_area);
    c->x2 = mid_width  + (dp_scale  * prob_area);
    c->y2 = mid_height + (dp_scale  * prob_area);

}
/***   GetProbCellCoords   ***/




/***   GetFullCellCoords   *******************************************************
 *
 * Write Coords for a fully occupied basepair cell (mic!) to
 * already allocated rect_coord_t *c
 *
 */
void
GetFullCellCoords(pair_t pair, double dp_scale, rect_coord_t *c)
{
    c->x1 = pair.ntj *dp_scale;
    c->y1 = pair.nti *dp_scale;
    c->x2 = c->x1+dp_scale;
    c->y2 = c->y1+dp_scale;
}
/***   GetFullCellCoords   ***/

