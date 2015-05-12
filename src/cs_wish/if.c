/******************************************************************************
*
* if.c - main interface for c<->tcl vice versa
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
 *  CVS $Id: if.c,v 1.60 2007-10-22 10:43:23 steger Exp $
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
#include <math.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "proj.h"
#include "mx.h"
#include "cs_dpm.h"
#include "consensus.h"
#include "bp_prob_mat.h"
#include "cs_seqio.h"
#include "mic.h"
#include "gui.h"
#include "seq_stat.h"
#include "shift.h"


#include "if.h"



/***   private
 */

#if 0
    #define DEBUG
#endif


/* BASEPAIR ARRAY INDICES
*/
#define BP_IDX2_PARTNER "partner"
#define BP_IDX2_PROB    "prob"

/* PAIRLIST ARRAY INDICES
*/
#define PAIRLIST_IDX2_PROB    "rel_prob"
#define PAIRLIST_IDX2_PAIRS   "pairs"

/* SEQ_EXCHANGE
 */
#define DIRECTION_TCL2C  "TCL2C"
#define DIRECTION_C2TCL  "C2TCL"

#define SEQ_ID_IDX      "id"
#define SEQ_NT_IDX      "nt"
#define ALN_LEN_IDX     "aln_len"
#define ALN_SEQ_NO_IDX  "n_seq"

/* COLORS
 */
#define COLOR_ARRAYNAME  "COLOR"
#define COLOR_DP_BG      "CDpBg"
#define COLOR_SEL_SEQ    "SelSeq"
#define COLOR_UNSEL_SEQ  "UnselSeq"
#define COLOR_ACT_NT     "ActiveNt"
#define COLOR_UNACT_NT   "UnactiveNt"
#define COLOR_SEL_BP     "SelBp"
#define COLOR_UNSEL_BP   "UnselBp"

/* OPTS
 */
#define CS_OPTS_ARRAYNAME "opts"
#define VERBOSE_IDX       "verbose"
#define DEBUG_IDX         "debug"
#define TIMING_IDX        "do_timing"

/*
 ***/






/***   GetCanvasDotplotSize   *************************************************
 *
 * Returns canvas size as double
 */
int
GetCanvasDotplotSize(Tcl_Interp *interp, double *dp_size)
{
    Tcl_Obj *cs_size_obj = NULL;

    cs_size_obj = Tcl_GetVar2Ex(interp, CANVAS_ARRAYNAME, CANVAS_SIZE_IDX_NAME,
                                           TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (cs_size_obj==NULL)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, cs_size_obj, dp_size)==TCL_ERROR)
        return TCL_ERROR;

    return TCL_OK;
}
/***   GetCanvasDotplotSize   ***/





/***   GetCanvasVirtualsize   *************************************************
 *
 * Returns canvas virtual size as double
 */
int
GetCanvasVirtualsize(Tcl_Interp *interp, double *dp_virtsize)
{
    Tcl_Obj *cs_virtsize_obj = NULL;

    cs_virtsize_obj = Tcl_GetVar2Ex(interp, CANVAS_ARRAYNAME, CANVAS_VIRSIZE_IDX_NAME,
                                                  TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (cs_virtsize_obj==NULL)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, cs_virtsize_obj, dp_virtsize)==TCL_ERROR)
        return TCL_ERROR;

    return TCL_OK;
}
/***   GetCanvasVirtualsize   ***/




/***   GetDpScale   ***********************************************************
 *
 * Returns dp canvas scale as double
 */
int
GetDpScale(Tcl_Interp *interp, double *scale_dp)
{
    Tcl_Obj *cs_scale_obj = NULL;

    cs_scale_obj = Tcl_GetVar2Ex(interp, SCALE_ARRAYNAME, SCALE_DP_IDX_NAME,
                                               TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (cs_scale_obj==NULL)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, cs_scale_obj, scale_dp)==TCL_ERROR)
        return TCL_ERROR;

    return TCL_OK;
}
/***   GetDpScale   ***/




/***   Seq_Exchange_Cmd   ******************************************************
 *
 * C:   zero offset
 * Tcl: unit-offset
 *
 */
int
Seq_Exchange_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    char *seq_arrayname = NULL;    /* argument */
    char *direction     = NULL;    /* argument */

    int  seq_idx = 0;
    /* we reallocate these for each seq instead
       of using a static buffer to be on the safe side */
    char *tcl_idx = NULL;
    char *tcl_val = NULL;
    char buf[1024];
    int len;


    if ( (objc!=3))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq_array direction?");
        return TCL_ERROR;
    }


    seq_arrayname = Tcl_GetStringFromObj(objv[1], NULL);
    if (seq_arrayname == NULL )
        return TCL_ERROR;

    direction = Tcl_GetStringFromObj(objv[2], NULL);
    if (seq_arrayname == NULL )
        return TCL_ERROR;


    if ( (! STR_EQ(direction, DIRECTION_TCL2C)) && (! STR_EQ(direction, DIRECTION_C2TCL)))
    {
        sprintf(buf, "direction must be one of %s, %s\n",
                        DIRECTION_TCL2C, DIRECTION_C2TCL);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /*** FIXME: implement reverse
     */
    if ( STR_EQ(direction, DIRECTION_TCL2C))
    {
        sprintf(buf, "direction %s not implemented yet\n", DIRECTION_TCL2C);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /*** foreach sequence
     */
    for (seq_idx=0; seq_idx < aln->num_seq; seq_idx++)
    {

        /*** export id
         */
        len = strlen(SEQ_ID_IDX) + 1 + INT_TO_STRLEN(seq_idx+1) + 1;
        tcl_idx = realloc(tcl_idx, len * sizeof(char));
        sprintf(tcl_idx, "%s,%d", SEQ_ID_IDX, seq_idx+1);

        len = strlen(aln->seq[seq_idx].id) + 1;
        tcl_val = realloc(tcl_val, len * sizeof(char));
        sprintf(tcl_val, "%s", aln->seq[seq_idx].id);

        #ifdef DEBUG
            DEBUG_P("setting %s(%s)=%s\n", seq_arrayname, tcl_idx, tcl_val);
        #endif
        if (Tcl_SetVar2(interp, seq_arrayname,
                        tcl_idx, tcl_val, TCL_LEAVE_ERR_MSG) == NULL)
        {
            free(tcl_idx);
            free(tcl_val);
            return TCL_ERROR;
        }

        /*** export nt
         */
        len = strlen(SEQ_NT_IDX) + 1 + INT_TO_STRLEN(seq_idx+1) + 1;
        tcl_idx = realloc(tcl_idx, len * sizeof(char));
        sprintf(tcl_idx, "%s,%d", SEQ_NT_IDX, seq_idx+1);

        len = strlen(aln->seq[seq_idx].nt) + 1;
        tcl_val = realloc(tcl_val, len * sizeof(char));
        sprintf(tcl_val, "%s", aln->seq[seq_idx].nt);

        #ifdef DEBUG
            DEBUG_P("setting %s(%s)=%s\n", seq_arrayname, tcl_idx, tcl_val);
        #endif
        if (Tcl_SetVar2(interp, seq_arrayname,
                        tcl_idx, tcl_val, TCL_LEAVE_ERR_MSG) == NULL)
        {
            free(tcl_idx);
            free(tcl_val);
            return TCL_ERROR;
        }


    } /* foreach sequence */



    /*****   export alignment length
     *
     */
    len = strlen(ALN_LEN_IDX) + 1;
    tcl_idx = realloc(tcl_idx, len * sizeof(char));
    sprintf(tcl_idx, "%s", ALN_LEN_IDX);

    len = INT_TO_STRLEN(aln->len) + 1;
    tcl_val = realloc(tcl_val, len * sizeof(char));
    sprintf(tcl_val, "%d", aln->len);

    #ifdef DEBUG
        DEBUG_P("setting %s(%s)=%s\n", seq_arrayname, tcl_idx, tcl_val);
    #endif
    if (Tcl_SetVar2(interp, seq_arrayname,
                    tcl_idx, tcl_val, TCL_LEAVE_ERR_MSG) == NULL)
    {
        free(tcl_idx);
        free(tcl_val);
        return TCL_ERROR;
    }


    /*** export number of seqs
     */
    len = strlen(ALN_SEQ_NO_IDX) + 1;
    tcl_idx = realloc(tcl_idx, len * sizeof(char));
    sprintf(tcl_idx, "%s", ALN_SEQ_NO_IDX);

    len = INT_TO_STRLEN(aln->num_seq) + 1;
    tcl_val = realloc(tcl_val, len * sizeof(char));
    sprintf(tcl_val, "%d", aln->num_seq);
    #ifdef DEBUG
        DEBUG_P("setting %s(%s)=%s\n", seq_arrayname, tcl_idx, tcl_val);
    #endif
    if (Tcl_SetVar2(interp, seq_arrayname,
                    tcl_idx, tcl_val, TCL_LEAVE_ERR_MSG) == NULL)
    {
        free(tcl_idx);
        free(tcl_val);
        return TCL_ERROR;
    }

    free(tcl_idx);
    free(tcl_val);
    return TCL_OK;
}
/***   Seq_Exchange_Cmd   ***/




/***   ExportBasepair   *******************************************************
 *
 * nt_idx and nt_partner must be zero-offset
 *
 * Set two Tcl-Values:
 *   bp_arrayname(nt_idx1, BP_IDX2_PARTNER)  -> nt_idx1
 *   bp_arrayname(nt_idx1, BP_IDX2_PROB)     -> prob
 *
 * return array with name basepair_arrayname (arg)
 *  index1: 1-alignmentlength
 *  index2: BP_IDX2_PARTNER|BP_IDX2_PROB,
 *      where partner=nt-index and prob=probability
*
 */
void
ExportBasepair(Tcl_Interp *interp, char *bp_arrayname,
                   int nt_idx, int nt_partner, float prob)
{
    char tcl_idx[1024], tcl_val[1024];

    /* tcl-vars got unit-offset */
    nt_idx++;
    nt_partner++;



    /***  export pair partner
     */
    sprintf(tcl_idx, "%d,%s", nt_idx, BP_IDX2_PARTNER);
    sprintf(tcl_val, "%d", nt_partner);

    #ifdef DEBUG
        DEBUG_P("Setting %s(%s)=%s\n", bp_arrayname, tcl_idx, tcl_val);
    #endif
    if (Tcl_SetVar2(interp, bp_arrayname, tcl_idx, tcl_val,
                                         TCL_LEAVE_ERR_MSG) == NULL)
    {
        ERROR_P("Couldn't set %s(%s) to %s\n", bp_arrayname, tcl_idx, tcl_val);
        return;
    }

    /***  export pair prob
     */
    sprintf(tcl_idx, "%d,%s", nt_idx, BP_IDX2_PROB);
    sprintf(tcl_val, "%f", prob);

    #ifdef DEBUG
        DEBUG_P("Setting %s(%s)=%s\n", bp_arrayname, tcl_idx, tcl_val);
    #endif
    if (Tcl_SetVar2(interp, bp_arrayname, tcl_idx, tcl_val,
                                          TCL_LEAVE_ERR_MSG) == NULL)
    {
        ERROR_P("Couldn't set %s(%s) to %s", bp_arrayname, tcl_idx, tcl_val);
        return;
    }

}
/***   ExportBasepair   ***/




/***   ExportHigherBasepair   *************************************************
 *
 * Export a higher (i.e. triple basepair) to Tcl
 * Works the same way as ExportBasepair, but
 * the values may contain a list
 *   bp_arrayname(nt_idx, BP_IDX2_PARTNER)  -> nt_idx
 *   bp_arrayname(nt_idx, BP_IDX2_PROB)     -> prob
 *
 *
 * nof_pairs denotes the basepair kind, for example 2 means basetriples
 * (and therefore nt_partner and pair_prob must contain 2 values)
 *
 *
 */
void
ExportHigherBasepair(Tcl_Interp *interp, char *bp_arrayname, int nt_idx,
                        int nof_pairs, int *nt_partner, float *pair_prob)
{
    char tcl_idx[1024], tcl_val[1024];
    int i;
    char tmpbuf[1024];

    /* tcl-vars got unit-offset */
    nt_idx++;
    for (i=0; i<nof_pairs; i++)
        nt_partner[i]++;


    /***  export pair partners
     */
    sprintf(tcl_idx, "%d,%s", nt_idx, BP_IDX2_PARTNER);
    sprintf(tcl_val, "%d", nt_partner[0]);
    for (i=1; i<nof_pairs; i++)
    {
        sprintf(tmpbuf, " %d", nt_partner[i]);
        strcat(tcl_val, tmpbuf);
    }

    #ifdef DEBUG
        DEBUG_P("Setting %s(%s) to %s\n", bp_arrayname, tcl_idx, tcl_val);
    #endif
    if (Tcl_SetVar2(interp, bp_arrayname, tcl_idx, tcl_val,
                                         TCL_LEAVE_ERR_MSG) == NULL)
    {
        ERROR_P("Couldn't set %s(%s) to %s\n", bp_arrayname, tcl_idx, tcl_val);
        return;
    }



    /***  export pair probs
     */
    sprintf(tcl_idx, "%d,%s", nt_idx, BP_IDX2_PROB);
    sprintf(tcl_val, "%f", pair_prob[0]);
    for (i=1; i<nof_pairs; i++)
    {
        sprintf(tmpbuf, " %f", pair_prob[i]);
        strcat(tcl_val, tmpbuf);
    }

    #ifdef DEBUG
        DEBUG_P("Setting %s(%s) to %s", bp_arrayname, tcl_idx, tcl_val);
    #endif
    if (Tcl_SetVar2(interp, bp_arrayname, tcl_idx, tcl_val,
                                          TCL_LEAVE_ERR_MSG) == NULL)
    {
        ERROR_P("Couldn't set %s(%s) to %s", bp_arrayname, tcl_idx, tcl_val);
        return;
    }
}
/***   ExportHigherBasepair   ***/




/***   GetColors   ************************************************************
 *
 * Space for variables will be allocated
 * Variables will be freed if non-NULL (i.e. )
 *
 * C:   zero offset
 * Tcl: unit-offset
 *
 */
int
GetColors(Tcl_Interp *interp, tk_color_t *color)
{
    Tcl_Obj *get_obj  = NULL;
    char    *get_char = NULL;



    /*** color->dp_bg
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_DP_BG,
                             TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;


    if (color->dp_bg!=NULL)
        free(color->dp_bg);
    color->dp_bg = strdup(get_char);



    /*** color->sel_seq
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_SEL_SEQ,
                               TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;

    if (color->sel_seq!=NULL)
        free(color->sel_seq);
    color->sel_seq = strdup(get_char);



    /*** color->unsel_seq
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_UNSEL_SEQ,
                                 TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;

    if (color->unsel_seq!=NULL)
        free(color->unsel_seq);
    color->unsel_seq = strdup(get_char);



    /*** color->act_nt
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_ACT_NT,
                                 TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;


    if (color->act_nt!=NULL)
        free(color->act_nt);
    color->act_nt = strdup(get_char);



    /*** color->unact_nt
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_UNACT_NT,
                                 TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;

    if (color->unact_nt!=NULL)
        free(color->unact_nt);
    color->unact_nt = strdup(get_char);



    /*** color->sel_bp
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_SEL_BP,
                                 TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;

    if (color->sel_bp!=NULL)
        free(color->sel_bp);
    color->sel_bp = strdup(get_char);



    /*** color->unsel_bp
     */
    get_obj = Tcl_GetVar2Ex(interp, COLOR_ARRAYNAME, COLOR_UNSEL_BP,
                                 TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    get_char = Tcl_GetStringFromObj(get_obj, NULL);
    if (get_char==NULL)
        return TCL_ERROR;

    if (color->unsel_bp!=NULL)
        free(color->unsel_bp);
    color->unsel_bp = strdup(get_char);

    return TCL_OK;
}
/***   GetColors   ***/






/***   Get_BpProb_Cmd   *******************************************************
 *
 * Interface to fetch Basepair Probs from within Tcl
 * Arguments: seq_no nt_no1 nt_no2
 *
 *
 */
int
Get_BpProb_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    bpair_t *bp_poi   = NULL;
    Tcl_Obj *obj_poi  = NULL;
    float    bp_prob  = 0.0;
    int      seq_idx  = -1;
    int nt_i =-1; int nt_j=-1;
    char buf[1024];
    pair_t pair;

    if (objc!=4)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq_no nt_no1>nt_no2?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &seq_idx) != TCL_OK)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, objv[2], &nt_i)    != TCL_OK)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, objv[3], &nt_j)    != TCL_OK)
        return TCL_ERROR;



    /*** Check for valid ranges
     *
     */
    if (seq_idx<0 || seq_idx>aln->num_seq)
    {
        sprintf(buf, "Sequence index %d is out of range\n", seq_idx);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    if (nt_i==nt_j)
    {
        sprintf(buf, "Basepair \"%d\":\"%d\" is invalid\n\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    if ( (NtIndexIsValid(nt_i)==NO) || (NtIndexIsValid(nt_i)==NO) )
    {
        sprintf(buf, "Nt indices \"%d:%d\" is out of valid range\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /* c is zero-offset */
    seq_idx--;
    nt_j--;
    nt_i--;

    SetPair(nt_i, nt_j, &pair);
    bp_poi = GetBpFromAlnSeq(aln, seq_idx, pair);

    if (bp_poi==NULL)
        bp_prob = 0;
    else
        bp_prob = bp_poi->prob;

    /*** set result
     */
    obj_poi = Tcl_NewDoubleObj((double)bp_prob);
    Tcl_SetObjResult(interp, obj_poi);

    return TCL_OK;
}
/***   Get_BpProb_Cmd   ***/




/***   Get_ConsBpProb_Cmd   *******************************************************
 *
 * Interface to fetch Consensus-Basepair Probs from within Tcl
 * nt_no1 nt_no2 (unit-offset)
 *
 */
int
Get_ConsBpProb_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj         *obj_poi  = NULL;
    float           bp_prob   = 0.0;
    int nt_i =-1; int nt_j=-1;
    char buf[1024];
    pair_t pair;

    if ( (objc!=3)) {
        Tcl_WrongNumArgs(interp, 1, objv, "?nt_no1 nt_no2?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &nt_i)!=TCL_OK)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, objv[2], &nt_j)!=TCL_OK)
        return TCL_ERROR;



    /*** Check for valid ranges
     *
     */
    if (nt_i==nt_j)
    {
        sprintf(buf, "Basepair \"%d\":\"%d\" is invalid\n\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    if ( (NtIndexIsValid(nt_i)==NO) || (NtIndexIsValid(nt_i)==NO) )
    {
        sprintf(buf, "Nt indices \"%d:%d\" is out of valid range\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /* c is zero-offset */
    nt_j--;
    nt_i--;


    SetPair(nt_i, nt_j, &pair);
    bp_prob = GetConsBpProb(aln, pair, 0.0);

    /*** set result
     */
    obj_poi = Tcl_NewDoubleObj((double)bp_prob);
    Tcl_SetObjResult(interp, obj_poi);

    return TCL_OK;
}
/***   Get_ConsBpProb_Cmd   ***/




/***   Get_MutInfo_Cmd   ******************************************************
 *
 */
int
Get_MutInfo_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj    *obj_poi  = NULL;
    int nt_i, nt_j;
    float mic;
    char  buf[1024];
    pair_t pair;


    if (objc!=3)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?nt_no1 nt_no2?");
        return TCL_ERROR;
    }

    if (Tcl_GetIntFromObj(interp, objv[1], &nt_i)!=TCL_OK)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, objv[2], &nt_j)!=TCL_OK)
        return TCL_ERROR;

    if (nt_i==nt_j)
    {
        sprintf(buf, "Basepair \"%d\":\"%d\" is invalid\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    if ( (NtIndexIsValid(nt_i)==NO) || (NtIndexIsValid(nt_i)==NO) )
    {
        sprintf(buf, "Nt indices \"%d:%d\" is out of valid range\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    if (! MIC_IsComputed())
    {
        sprintf(buf, "No mutual information content found! Compute it before calling me.");
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    /* c is zero-offset */
    nt_j--;
    nt_i--;

    SetPair(nt_i, nt_j, &pair);

    /*** set result
     */
    mic = (double)GetMIC(pair);

    obj_poi = Tcl_NewDoubleObj(mic);
    Tcl_SetObjResult(interp, obj_poi);

    return TCL_OK;
}
/***   Get_MutInfo_Cmd   ***/




/***   NtIndexIsValid   *******************************************************
 *
 * Check if a given nt-index is > 0 and <=aln_len
 *
 *
 */
int
NtIndexIsValid(int nt_idx)
{
    if ( (nt_idx<=0) || (nt_idx>aln->len))
        return NO;
    else
        return YES;
}
/***   NtIndexIsValid   ***/




/***   Get_ConSeq_Cmd   *******************************************************
 *
 *
 */
int
Get_ConsSeq_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *obj_poi  = NULL;

    if (objc!=1)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?no args?");
        return TCL_ERROR;
    }

    obj_poi = Tcl_NewStringObj(aln->cons_seq->nt, -1);

    Tcl_SetObjResult(interp, obj_poi);

    return TCL_OK;
}
/***   Get_ConSeq_Cmd   ***/




/***   Version_Cmd   **********************************************************
 *
 *   return package name, version number ...
 */
int
Version_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *obj_poi;

    obj_poi = Tcl_NewStringObj(VERSION_STR, -1);

    Tcl_SetObjResult(interp, obj_poi);

    return TCL_OK;
}
/***   Version_Cmd   ***/




/***   LibZ_supported_Cmd   ***************************************************
 *
 *   return TRUE if libZ is supported, FALSE otherwise
 */
int
LibZ_supported_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *b=NULL;

	#ifdef HAVE_LIBZ
        b=Tcl_NewBooleanObj(1);
	#else
        b=Tcl_NewBooleanObj(0);
	#endif

    Tcl_SetObjResult(interp, b);

    return TCL_OK;
}
/***   LibZ_supported_Cmd   ***/




/***   GetCsOpts   ************************************************************
 *
 *
 */
int
GetCsOpts(Tcl_Interp *interp, cs_opts_t *cs_opts)
{
    Tcl_Obj *get_obj  = NULL;


    /*** cs_opts->verbose
     */
    get_obj = Tcl_GetVar2Ex(interp, CS_OPTS_ARRAYNAME, VERBOSE_IDX,
                             TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, get_obj, &cs_opts->verbose)!=TCL_OK)
        return TCL_ERROR;



    /*** cs_opts->debug
     */
    get_obj = Tcl_GetVar2Ex(interp, CS_OPTS_ARRAYNAME, DEBUG_IDX,
                             TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, get_obj, &cs_opts->debug)!=TCL_OK)
        return TCL_ERROR;



    /*** cs_opts->do_timing
     */
    get_obj = Tcl_GetVar2Ex(interp, CS_OPTS_ARRAYNAME, TIMING_IDX,
                             TCL_LEAVE_ERR_MSG | TCL_GLOBAL_ONLY);
    if (get_obj==NULL)
        return TCL_ERROR;

    if (Tcl_GetIntFromObj(interp, get_obj, &cs_opts->do_timing)!=TCL_OK)
        return TCL_ERROR;


    return TCL_OK;
}
/***   GetCsOpts   ***/




/***   LoadAlignment_Cmd   ****************************************************
 *
 * Load an alignment via cs_seqio and export to tcl via Seq_Exchange_Cmd
 * used sequence array is seq (see: Seq_Exchange_Cmd for indices)
 * ! Not for use in cs_dp (loads alignment via project) !
 *
 * args: alignment_filename
 *       sequence array name
 *
 *
 */
int
LoadAlignment_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{

    char *aln_file      = NULL;
    char *seq_arrayname = NULL;
    Tcl_Obj *op[3];


    if ( (objc!=3))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?alignment_file sequence_arrayname?");
        return TCL_ERROR;
    }


    /*** get arg: alignment filename
     *
     */
    aln_file = Tcl_GetStringFromObj(objv[1], NULL);
    if (aln_file == NULL)
        return TCL_ERROR;

    /*** get seq_arrayname
     *
     */
    seq_arrayname = Tcl_GetStringFromObj(objv[2], NULL);
    if (seq_arrayname == NULL )
        return TCL_ERROR;


    /* sanity check
     */
    if (aln!=NULL)
    {
        WARN_P("%s\n", "must kill existant alignment!");

        KillAln(aln);
    }

    if ((aln = ReadAlnFile(aln_file)) == NULL)
    {
        char buf[1024];
        sprintf(buf, "couldn't read sequence file %s!", aln_file);
        CsDpError(interp, buf);

        return TCL_ERROR;
    }

    op[0] = Tcl_NewStringObj("Seq_Exchange_Cmd", -1);
    op[1] = Tcl_NewStringObj(seq_arrayname, -1);
    op[2] = Tcl_NewStringObj(DIRECTION_C2TCL, -1);
    Tcl_IncrRefCount(op[0]);
    Tcl_IncrRefCount(op[1]);
    Tcl_IncrRefCount(op[2]);

    if (Seq_Exchange_Cmd(NULL, interp, 3, op)==TCL_ERROR)
    {
        Tcl_DecrRefCount(op[0]);
        Tcl_DecrRefCount(op[1]);
        Tcl_DecrRefCount(op[2]);
        KillAln(aln);
        /* Why this hack? */
        aln = NULL;
        return TCL_ERROR;
    }
    Tcl_DecrRefCount(op[0]);
    Tcl_DecrRefCount(op[1]);
    Tcl_DecrRefCount(op[2]);
    /* Why is this hack needed ? */
    KillAln(aln);
    aln = NULL;

    return TCL_OK;
}
/***   LoadAlignment_Cmd   ***/





/***   CsDpError   ************************************************************
 *
 */
void
CsDpError(Tcl_Interp *interp, char *msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    Tcl_AppendResult(interp, msg, "\n", (char *) NULL);
}
/***   CsDpError   ***/




/***   GetRgb_Cmd   ***********************************************************
 *
 * Tcl-Frontend to GetRgb
 */
int
GetRgb_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    char    rgb[1024];
    double  hue, lowlimit, highlimit;
    Tcl_Obj *obj_poi;


    if ( (objc!=4))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?hue lowlimit highlimit?");
        return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, objv[1], &hue)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[2], &lowlimit)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &highlimit)==TCL_ERROR)
        return TCL_ERROR;


    GetRgb(rgb, (float) hue, (float) lowlimit, (float) highlimit);

    obj_poi = Tcl_NewStringObj(rgb, -1);
    Tcl_SetObjResult(interp, obj_poi);

    #ifdef DEBUG
        DEBUG_P("returning %s (hue=%f, lowlimit=%f, highlimit=%f)\n",
                Tcl_GetStringFromObj(obj_poi, NULL), hue, lowlimit, highlimit);
    #endif

    return TCL_OK;
}
/***   GetRgb_Cmd   ***/





/***   ComputeMutualInfo_Cmd   ***************************************************
 *
 *
 */
int
ComputeMutualInfo_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    int bit ; /* boolean argument */
    int unbiased; /* boolean argument */
    int pair_entropy_norm; /* boolean argument */
    int use_alifoldscore;   /* boolean argument */
    int use_stacking;       /* boolean argument */
    int alifoldscore_type=0;
    char  **alnseq;
    int i;
    Stopwatch_t     *watch;

    if ( (objc!=6))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?bool:unbiased bool:bit  bool:pair_entropy_norm bool:use_alifoldscore bool:use_stacking?");
        return TCL_ERROR;
    }

    if (Tcl_GetBooleanFromObj(interp, objv[1], &unbiased)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[2], &bit)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[3], &pair_entropy_norm)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[4], &use_alifoldscore)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[5], &use_stacking)==TCL_ERROR)
        return TCL_ERROR;

    watch = TimerStart();


    alnseq = Scalloc(aln->num_seq, sizeof(char*));
    for (i=0; i<aln->num_seq; i++)
        alnseq[i] = aln->seq[i].nt;

    if (use_alifoldscore) {
        if (use_stacking==0) {
            alifoldscore_type=1;
        } else {
            alifoldscore_type=2;
        }
    }
    ComputeMutualInfoContent(alnseq, aln->len, aln->num_seq,
                             unbiased, bit, pair_entropy_norm,
                             alifoldscore_type);

    TimerStop(__FUNCTION__, watch);

    free(alnseq);

	return TCL_OK;
}
/***   ComputeMutualInfo_Cmd   ***/





/***   IsBasepair_Cmd   ***************************************************
 *
 *   return TRUE if two nucleotides
 *   may form a basepair
 *   Wobble-Bp (if allow_wobble==1) and DNA allowed
 */
int
IsBasepair_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *bool_obj;
    char    *c_nt1, *c_nt2;
    int     b=0;


    if ( (objc!=3)) {
        Tcl_WrongNumArgs(interp, 1, objv, "?nt nt?");
        return TCL_ERROR;
    }

    c_nt1 = Tcl_GetStringFromObj(objv[1], NULL);
    if (c_nt1 == NULL )
        return TCL_ERROR;

    c_nt2 = Tcl_GetStringFromObj(objv[2], NULL);
    if (c_nt2 == NULL )
        return TCL_ERROR;


    b = IsBasepair(c_nt1[0], c_nt2[0]);

    bool_obj = Tcl_NewBooleanObj(b);

    Tcl_SetObjResult(interp, bool_obj);

    return TCL_OK;
}
/***   IsBasepair_Cmd   ***/




/***   IsGap_Cmd   ************************************************************
 *
 *   return TRUE if a nucleotide symbol is a gap
 */
int
IsGap_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *bool_obj;
    char *c_nt;

    if ( (objc!=2))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?nt_smybol?");
        return TCL_ERROR;
    }

    c_nt = Tcl_GetStringFromObj(objv[1], NULL);
    if (c_nt == NULL )
        return TCL_ERROR;

    bool_obj = Tcl_NewBooleanObj(IS_GAP(c_nt[0]));

    Tcl_SetObjResult(interp, bool_obj);

    return TCL_OK;
}
/***   IsGap_Cmd   ***/





/***   Degap_Cmd   ************************************************************
 *
 *
 */
int
Degap_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *seq_poi;
    char    *str_aln;
    char    *str_dealn;

    if ( (objc!=2))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?aln_seq?");
        return TCL_ERROR;
    }

    str_aln = Tcl_GetStringFromObj(objv[1], NULL);
    if (str_aln == NULL )
        return TCL_ERROR;

    str_dealn = (char*) Scalloc(strlen(str_aln)+1, sizeof(char));
    Degap(str_aln, str_dealn);

    seq_poi = Tcl_NewStringObj(str_dealn, -1);
    Tcl_SetObjResult(interp, seq_poi);

    free(str_dealn);

    return TCL_OK;

}
/***   Degap_Cmd   ***/




/***   ConsBpInfo_Cmd   ***********************************************************
 *
 * Returns some info about one particular consensus-bp
 * (which must exist!)
 *
 */
int
ConsBpInfo_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    int   nt_i, nt_j;
    char  buf[1024];
    int          *contrib_seqidx;
    bpair_t      *contrib_bp;
    cons_bpair_t *csbp;
    dlist_elem   *contrib_list_elem = NULL;
    Tcl_Obj *op;
    float cprob;
    float cmic, mic;
    pair_t pair;


    /***  get and check args
     *
     */
    if (objc!=3)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?nt_no1>nt_no2?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &nt_i)    != TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[2], &nt_j)    != TCL_OK)
        return TCL_ERROR;

    if (nt_i==nt_j)
    {
        sprintf(buf, "Basepair \"%d\":\"%d\" is invalid\n\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    if ( (NtIndexIsValid(nt_i)==NO) || (NtIndexIsValid(nt_i)==NO) )
    {
        sprintf(buf, "Nt indices \"%d:%d\" is out of valid range\n", nt_i, nt_j);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    nt_j--;
    nt_i--;


    sprintf(buf, "Info for Consensus-Basepair %d:%d\n", nt_j+1, nt_i+1);
    op = Tcl_NewStringObj(buf, -1);

    SetPair(nt_i, nt_j, &pair);


    /*****   MIC and TD
     *
     */
    if (MIC_IsComputed())
    {
        cmic   = GetConsMIC(pair, 0.0);
        mic    = GetMIC(pair);
        sprintf(buf, "  MIC = %0.4f (normed: %0.4f)\n", mic, cmic);
        Tcl_AppendStringsToObj(op, buf, NULL);
    }
    else
    {
        sprintf(buf, "  No MIC computed\n");
        Tcl_AppendStringsToObj(op, buf, NULL);
    }
    #ifdef DEBUG
        DEBUG_P("%s\n", buf);
    #endif


    cprob  = GetConsBpProb(aln, pair, 0.0);
    sprintf(buf, "  Consensus-Prob. = %0.4f\n", cprob);
    Tcl_AppendStringsToObj(op, buf, NULL);
    #ifdef DEBUG
        DEBUG_P("%s\n", buf);
    #endif


    /***   Contribs
     */
    sprintf(buf, "  Contrib. seqs:\n");
    Tcl_AppendStringsToObj(op, buf, NULL);
    csbp = GetConsBp(aln, pair);

    contrib_list_elem = dlist_head(csbp->contribs);
    while (1)
    {
        contrib_seqidx   = (int*)dlist_data(contrib_list_elem);
        contrib_bp = GetBpFromAlnSeq(aln, *contrib_seqidx, pair);
        if (contrib_bp==NULL)
        {
            ERROR_P("%s\n", "Uups..Got NULL pointer as contrib");
            continue;
        }
        sprintf(buf, "    %-20s (prob %0.4f)\n", aln->seq[*contrib_seqidx].id, contrib_bp->prob);
        Tcl_AppendStringsToObj(op, buf, NULL);
        #ifdef DEBUG
            DEBUG_P("%s\n", buf);
        #endif

        if (dlist_is_tail(contrib_list_elem))
            break;
        contrib_list_elem = dlist_next(contrib_list_elem);
    }


    Tcl_SetObjResult(interp, op);

    return TCL_OK;
}
/***   ConsBpInfo_Cmd   ***/




/***   GetThisDpBpTag   ************************************************************
 *
 * Get Dotplot basepair tag for one particular basepair
 *
 * tag must be already allocated
 */
void
GetThisDpBpTag(int seq_idx, pair_t pair, char *tag)
{
    sprintf(tag, "DpBp_seq_%d_nti_%d_ntj_%d", seq_idx+1, pair.nti+1, pair.ntj+1);
}
/***   GetThisDpBpTag   ***/




/***   GetGeneralDpBpTags   ***************************************************
 *
 * Get Dotplot basepair tag for a group of basepairs
 *
 * tag must be already allocated
 */
void
GetGeneralDpBpTags(int seq_idx, char *tag)
{
    sprintf(tag, "DpBp DpBp_seq_%d", seq_idx+1);
}
/***   GetGeneralDpBpTags   ***/




/***   GetThisDpConsBpTag   **************************************************
 *
 * Get Dotplot consensus basepair tag for one particular consensus basepair
 *
 * tag must be already allocated
 */
void
GetThisDpConsBpTag(pair_t pair, char *tag)
{
    sprintf(tag, "DpConsBp_nti_%d_ntj_%d", pair.nti+1, pair.ntj+1);
}
/***   GetThisDpConsBpTag   ***/




/***   GetGeneralDpConsBpTags   ***********************************************
 *
 * Get the general dotplot consensus basepair tag
 * tag must be already allocated
 */
void
GetGeneralDpConsBpTags(char *tag)
{
    sprintf(tag, "%s", "DpConsBp");
}
/***   GetGeneralDpConsBpTags   ***/





/***   SumOfPairs_Cmd   *******************************************************
 *
 * Calculates the Sum-Of-Pairs Sequence Alignment Score
 * The sop will be returned to interpreter
 * and each column score will be written (unit-offset) to array named
 * column_cost_arrayname
 * max col score == number of all possible pairs (n*(n-1))/2
 */
int
SumOfPairs_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{

    int    sop;
    char **aln_seq;
    char  *colarray, colval[1024], colidx[1024];
    int   *col_cost;
    int i;
    Tcl_Obj *obj_ptr;

    if (objc!=2)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?[column_cost_arrayname]?");
        return TCL_ERROR;
    }

    colarray = Tcl_GetStringFromObj(objv[1], NULL);
    if (colarray == NULL )
        return TCL_ERROR;



    /*** setup vars prior to Sop call
     */
    aln_seq = Scalloc(aln->num_seq, sizeof(char*));
    for (i=0; i<aln->num_seq; i++)
        aln_seq[i] = aln->seq[i].nt;
	col_cost = (int *) Scalloc(aln->len, sizeof(int));


    /***
     */
    sop = Sop(aln_seq, aln->len, aln->num_seq, col_cost);
    /*
     ***/


    /* export column cost
     */
    for (i=0; i<aln->len; i++)
    {
        sprintf(colidx, "%d", i+1);
        sprintf(colval, "%d", col_cost[i]);
        #ifdef DEBUG
            DEBUG_P("Setting %s(%s) to %s\n", colarray, colidx, colval);
        #endif
        if (Tcl_SetVar2(interp, colarray, colidx, colval, TCL_LEAVE_ERR_MSG) == NULL)
        {
            char buf[1024];
            sprintf(buf, "Couldn't set %s(%s) to %s\n", colarray, colidx, colval);
            CsDpError(interp, buf);

            free(aln_seq);
            free(col_cost);

            return TCL_ERROR;
        }
    }

    free(col_cost);
    free(aln_seq);

    /***** Return Sum-of-pairs-Score
     */
    obj_ptr = Tcl_NewIntObj(sop);
    Tcl_SetObjResult(interp, obj_ptr);


    return TCL_OK;
}
/***   SumOfPairs_Cmd   ***/




/***   PwIdent_Cmd   *******************************************************
 *
 * return average pairwise identity between all seq pairs
 */
int
PwIdent_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{

    double   pwi;
    char   **aln_seq;
    Tcl_Obj *obj_ptr;
    int i;

    if (objc!=1)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?NONE?");
        return TCL_ERROR;
    }



    /*** setup aln_seq prior to PwIdent call
     */
    aln_seq = Scalloc(aln->num_seq, sizeof(char*));
    for (i=0; i<aln->num_seq; i++)
        aln_seq[i] = aln->seq[i].nt;


    /***
     */
    pwi = PwIdent(aln_seq, aln->len, aln->num_seq);
    /*
     ***/


    free(aln_seq);

    /***** Return Average Pairwise Identity
     */
    obj_ptr = Tcl_NewDoubleObj(pwi);
    Tcl_SetObjResult(interp, obj_ptr);


    return TCL_OK;
}
/***   PwIdent_Cmd   ***/




/***   PrintCsTdProbMat_Cmd   ***************************************************
 *
 * Print consensus td prob mat
 */
int
PrintCsTdProbMat_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[])
{
    char *filename;
    FILE *fstream;

    if (objc!=2)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?filename?");
        return TCL_ERROR;
    }

    filename = Tcl_GetStringFromObj(objv[1], NULL);
    if (filename == NULL )
        return TCL_ERROR;

    fstream=fopen(filename, "w");
    if (fstream==NULL)
    {
        char buf[1024];
        sprintf(buf, "Couldn't open %s", filename);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    PrintCsTdProbMat(aln->cons_seq, fstream);

    fclose(fstream);

    return TCL_OK;
}
/***   PrintCsTdProbMat_Cmd   ***/




/***   PrintTdProbMat_Cmd   *****************************************************
 *
 * Print td-prob-mat of one sequence
 */
int
PrintTdProbMat_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[])
{
    char *filename;
    FILE *fstream;
    int   seq_idx;
    char *nalseq;

    if (objc!=3)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq_idx filename?");
        return TCL_ERROR;
    }

    if (Tcl_GetIntFromObj(interp, objv[1], &seq_idx) != TCL_OK)
        return TCL_ERROR;
    filename = Tcl_GetStringFromObj(objv[2], NULL);
    if (filename == NULL)
        return TCL_ERROR;

    fstream=fopen(filename, "w");
    if (fstream==NULL)
    {
        char buf[1024];
        sprintf(buf, "Couldn't open %s", filename);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    nalseq = (char*) Scalloc(strlen(aln->seq[seq_idx-1].nt)+1, sizeof(char));
    Degap(aln->seq[seq_idx-1].nt, nalseq);
    PrintTdProbMat(nalseq, aln->seq[seq_idx-1].bp, fstream);

    fclose(fstream);

    return TCL_OK;
}
/***   PrintTdProbMat_Cmd   ***/




/***   PrintCsDpm_Cmd   ***********************************************************
 *
 * print merged mic cs-td.prob matrix
 */
int
PrintCsDpm_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[])
{
    char *filename;
    FILE *fstream;
    double  prob_fac, mic_fac;
    double  prob_thr, mic_thr;
    float **mpm;    /* the consensus matrix, merged out of mic and td */

    if (objc!=6)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?filename prob_fac mic_fac prob_thr mic_thr?");
        return TCL_ERROR;
    }

    filename = Tcl_GetStringFromObj(objv[1], NULL);
    if (filename == NULL)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[2], &prob_fac)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &mic_fac) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[4], &prob_thr)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[5], &mic_thr) ==TCL_ERROR)
        return TCL_ERROR;

    #ifdef DEBUG
        DEBUG_P("filename=%s\n", filename);
        DEBUG_P("prob_fac=%f\n", prob_fac);
        DEBUG_P("mic_fac=%f\n", mic_fac);
        DEBUG_P("prob_thr=%f\n", prob_thr);
        DEBUG_P("mic_thr=%f\n", mic_thr);
    #endif

    fstream=fopen(filename, "w");
    if (fstream==NULL)
    {
        char buf[1024];
        sprintf(buf, "Couldn't open %s", filename);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    mpm = CsDpm((float)prob_fac, (float)mic_fac, (float)prob_thr, (float)mic_thr);

    if (mpm==NULL)
    {
        CsDpError(interp, "CsDpm failed");
        return TCL_ERROR;
    }

    PrintCsDpm(aln->cons_seq->nt, mpm, fstream);

    FreeCsDpm(mpm);

    return TCL_OK;
}
/***   PrintCsDpm_Cmd   ***/




/***   PrintMic_Cmd   *********************************************************
 *
 */
int
PrintMic_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[])
{
    char *filename;
    FILE *fstream;

    if (objc!=2)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?filename?");
        return TCL_ERROR;
    }

    filename = Tcl_GetStringFromObj(objv[1], NULL);
    if (filename == NULL )
        return TCL_ERROR;

    fstream=fopen(filename, "w");
    if (fstream==NULL)
    {
        char buf[1024];
        sprintf(buf, "Couldn't open %s", filename);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }

    PrintMic(aln->cons_seq->nt, fstream);

    fclose(fstream);

    return TCL_OK;
}
/***   PrintMic_Cmd   ***/





/***   ExportSuboptimalPairlist   *********************************************
 *
 * Exports a suboptimal pairlist to tcl to selecting one pair for backtracking
 *
 * Sets the following Tcl-Values:
 *  pairlist_arrayname(max_no, PAIRLIST_IDX2_PROB)  -> prob for a list of basepairs
 *  pairlist_arrayname(max_no, PAIRLIST_IDX2_PAIRS) -> list of basepairs for this maximum
 *
 * the number of maximum basepairs is returned to the interpreter
 *
 */
int
ExportSuboptimalPairlist(Tcl_Interp *interp, char *pairlist_arrayname, int *pairlist)
{
    int       i;
    int       max_prob;
    char      tcl_idx[1024], tcl_val[1024];
    int       list_len;
    Tcl_Obj  *obj_ptr;
    int       n_maxbps;
    float    *probs;
    char    **bps;

    /* old code :
     *   /@ return to tcl with bp i j, probability and number of maximum @/
     *   sprintf(dummy_string, "%d 0 0 0", pairlist[0]);
     *   Tcl_AppendElement(interp, dummy_string);
     *   for (i=0; i<pairlist[0]; i++)
     *   {
     *       sprintf(dummy_string, "%d %d %d %d", pairlist[i*4+1], pairlist[i*4+2], pairlist[i*4+3], pairlist[i*4+4]+1);
     *       Tcl_AppendElement(interp, dummy_string);
     *   }
     *
     * bi             = pairlist[i*4+1]
     * bj             = pairlist[i*4+2]
     * prob           = pairlist[i*4+3]
     * maxbpctr       = pairlist[i*4+4]+1
     * pairlistlength = pairlist[0]
     *
     */


    /* no structure found ? */
    if (pairlist==NULL)
        return TCL_OK;
    
    list_len = pairlist[0];
    max_prob = pairlist[1*4+3]; /* list is sorted -> first prob is the max*/


    /* find the number of maximums
     */
    n_maxbps = 0;
    for (i=0; i<list_len; i++)
    {
        int maxbpctr = pairlist[i*4+4]+1;
        int prob     = pairlist[i*4+3];
        if (maxbpctr > n_maxbps)
            n_maxbps = maxbpctr;
        if (prob > max_prob)
            max_prob = prob;
    }
    #ifdef DEBUG
        DEBUG_P("n_maxbps=%d\n", n_maxbps);
        DEBUG_P("max_prob=%d\n", max_prob);
        DEBUG_P("list_len=%d\n", list_len);
    #endif

    probs = (float*) Scalloc(n_maxbps+1, sizeof(float));
    bps   = (char**) Scalloc(n_maxbps+1, sizeof(char*));
    for (i=1; i<=n_maxbps; i++)
    {
        bps[i]    = (char*)Scalloc(1024, sizeof(char));
        bps[i][0] = '\0';
    }

    /*   there is one prob for each maxbpctr
     *   and
     *   a list of basepairs (with the same probs)
     *   unfortunately we have to iterate over complete list
     */
    for (i=0; i<list_len; i++)
    {
        int   maxbpctr, bi, bj, prob;
        float rel_prob;
        char  cpybuf[1024];

        bi = pairlist[i*4+1];
        bj = pairlist[i*4+2];
        maxbpctr = pairlist[i*4+4]+1;
        prob = pairlist[i*4+3];

        /* this is redundant, probs[maxbpctr] will be overwritten
           several times with the same value
         */
        rel_prob = (float)prob/(float)max_prob;
        probs[maxbpctr] = rel_prob;

        /* write or append to bplist */
        if (bps[maxbpctr][0]=='\0')
        {
            sprintf(bps[maxbpctr], "%d:%d",bi, bj);
        }
        else
        {
            sprintf(cpybuf, "%s", bps[maxbpctr]);
            sprintf(bps[maxbpctr], "%s %d:%d", cpybuf, bi, bj);
        }
        #ifdef DEBUG
            DEBUG_P("probs[%d]=%f\n", maxbpctr, probs[maxbpctr]);
            DEBUG_P("bps[%d]=%s\n",   maxbpctr, bps[maxbpctr]);
        #endif
    }


    /* now we got the probs and bplist for each maximum
     * -> export
     */
    for (i=1; i<=n_maxbps; i++)
    {

        /***  export prob for this maximum
         */
        sprintf(tcl_idx, "%d,%s", i, PAIRLIST_IDX2_PROB);
        sprintf(tcl_val, "%f", probs[i]);
        #ifdef DEBUG
            DEBUG_P("Setting %s(%s)=%s\n", pairlist_arrayname, tcl_idx, tcl_val);
        #endif
        if (Tcl_SetVar2(interp, pairlist_arrayname, tcl_idx, tcl_val,
                                              TCL_LEAVE_ERR_MSG) == NULL)
        {
            ERROR_P("Couldn't set %s(%s) to %s", pairlist_arrayname, tcl_idx, tcl_val);

            free(probs);
            for (i=0; i<n_maxbps; i++)
              free(bps[i]);
            free(bps);

            return TCL_ERROR;
        }


        /***  export bplist for this maximum
         */
        sprintf(tcl_idx, "%d,%s", i, PAIRLIST_IDX2_PAIRS);
        sprintf(tcl_val, "%s", bps[i]);
        #ifdef DEBUG
            DEBUG_P("Setting %s(%s)=%s\n", pairlist_arrayname, tcl_idx, tcl_val);
        #endif
        if (Tcl_SetVar2(interp, pairlist_arrayname, tcl_idx, tcl_val,
                                              TCL_LEAVE_ERR_MSG) == NULL)
        {
            ERROR_P("Couldn't set %s(%s) to %s", pairlist_arrayname, tcl_idx, tcl_val);

            free(probs);
            for (i=0; i<n_maxbps; i++)
              free(bps[i]);
            free(bps);

            return TCL_ERROR;
        }
    }

    /* cleanup
     */
    free(probs);
    for (i=0; i<n_maxbps; i++)
        free(bps[i]);
    free(bps);


    obj_ptr = Tcl_NewIntObj(n_maxbps);
    Tcl_SetObjResult(interp, obj_ptr);

    return TCL_OK;
}
/***   ExportSuboptimalPairlist   ***/





/***   ConvertAln_Cmd   *******************************************************
 *
 * args:
 *  1: old seq file
 *  2: new seq file
 *  3: format
 */
int
ConvertAln_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    char *fname_in, *fname_out, *ext;


    if (objc!=4)
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                          "?<aln_in> <aln_out> <extension>?");
        return TCL_ERROR;
    }

    fname_in = Tcl_GetStringFromObj(objv[1], NULL);
    if (fname_in == NULL )
        return TCL_ERROR;
    fname_out = Tcl_GetStringFromObj(objv[2], NULL);
    if (fname_out == NULL )
        return TCL_ERROR;
    ext = Tcl_GetStringFromObj(objv[3], NULL);
    if (ext == NULL )
        return TCL_ERROR;

    #ifdef DEBUG
        DEBUG_P("fname_in = %s\n",  fname_in);
        DEBUG_P("fname_out = %s\n", fname_out);
        DEBUG_P("ext       = %s\n", ext);
    #endif


    if (FileExists(fname_in)==FALSE)
    {
        CsDpError(interp, "Input file does not exist\n");
        return TCL_ERROR;
    }

    if (ConvertAndUpdateAln(aln, fname_in, fname_out, ext)==ERROR)
    {
        CsDpError(interp, "Conversion failed\n");
        return TCL_ERROR;
    }


    return TCL_OK;
}
/***   ConvertAln_Cmd   ***/





/***   CsShiftScore_Cmd   *******************************************************
 *
 * return average pairwise identity between all seq pairs
 */
int
CsShiftScore_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    Tcl_Obj *obj_ptr;
    char    *faln_pred, *faln_trust, *f_mask;
    aln_t   *aln_trust, *aln_pred;
    int     *shift;
    int ovl_shift, i;


    if (objc<3 || objc>4 )
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                          "?aln_file_trusted aln_file_predicted [pairmask]?");
        return TCL_ERROR;
    }

    faln_trust = Tcl_GetStringFromObj(objv[1], NULL);
    if (faln_trust == NULL)
        return TCL_ERROR;
    faln_pred = Tcl_GetStringFromObj(objv[2], NULL);
    if (faln_pred == NULL)
        return TCL_ERROR;
    if (objc==4)
    {
        f_mask = Tcl_GetStringFromObj(objv[3], NULL);
        if (f_mask == NULL)
            return TCL_ERROR;
    }
    else
    {
        f_mask = NULL;
    }
    TMPDEBUG_P("faln_trust = %s\n", faln_trust);
    TMPDEBUG_P("faln_pred  = %s\n", faln_pred);
    TMPDEBUG_P("f_mask     = %s\n", f_mask);

    TMPDEBUG_P("%s\n", "Help: implement pairmask as int array");

    if (GetCsOpts(interp, &cs_opts)  == TCL_ERROR)
        return TCL_ERROR;

    aln_trust = ReadAlnFile(faln_trust);
    aln_pred  = ReadAlnFile(faln_pred);

    shift = CsShiftScore(aln_trust, aln_pred, NULL);
    if (shift==NULL)
        return TCL_ERROR;

    ovl_shift = 0;
    for (i=0; i<aln_pred->num_seq; i++)
        ovl_shift += shift[i];

    /***** Return Average Pairwise Identity
     */
    obj_ptr = Tcl_NewDoubleObj((double)ovl_shift/(double)aln_pred->num_seq);
    Tcl_SetObjResult(interp, obj_ptr);


    KillAln(aln_trust);
    KillAln(aln_pred);

    return TCL_OK;
}
/***   CsShiftScore_Cmd   ***/
