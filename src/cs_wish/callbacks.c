/******************************************************************************
* 
* callbacks.c - routines for handling gui callbacks from cs_dp
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
 *  CVS $Id: callbacks.c,v 1.26 2007-10-22 10:43:23 steger Exp $    
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
#include <string.h>

#include "public_datatypes.h"
#include "bp_prob_mat.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "if.h"

#include "callbacks.h"



/***   private
 */
#if 0
    #define DEBUG
#endif 

static int
Colorize_AlnNts(Tcl_Interp *interp, int seq_idx, pair_t pair, char *color);

/* Returns all pair indices for a particular nt 
 * where prob>0
 */
static int *
GetBpPartnersForNt(aln_t *aln, int seq_no, int nt_no);

/*
 ***/

/***   Highlight_DpBp_Cmd   ***************************************************
 *
 *
 */
int
Highlight_DpBp_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    int         seq_no;
    int         nt_no;
    int        *bp_indices;
    int         state       = 0;       /* 1 if state_str = on                */
    char       *state_str;             /* 0 if state_str = off               */
    Tcl_Obj    *cmd_obj     = NULL;    /* used for executing tk commands     */
    char        cmd[1024];      /* used for tcl/tk command formatting */
    char        tag[1024];
    pair_t      pair;
    int i=0;
   
   
    /*** check and get args
     *
     */    
    if ( (objc!=4))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq_no nt_no state?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &seq_no)!=TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[2], &nt_no)!=TCL_OK)
        return TCL_ERROR;
    state_str = Tcl_GetStringFromObj(objv[3], NULL);
    if (state_str == NULL)
    {
        CsDpError(interp, "Can't get state");
        return TCL_ERROR;
    }

    if (STR_EQ(state_str, "off"))
    {
        state = 0;
    }
    else if (STR_EQ(state_str, "on"))
    {
        state = 1;
    }
    else
    {
        CsDpError(interp, "State must be one of (on|off)");
        return TCL_ERROR;
    }
    
    
    /* c is zero-offset ! */
    seq_no -= 1;
    nt_no  -= 1;
    
    
    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);
    
    /* FIXME: check range */
    
    /*** get all bp indices and highlight them
     */    
    bp_indices = GetBpPartnersForNt(aln, seq_no, nt_no);
    i=0;
    while (bp_indices[i] != -1)
    {
        SetPair(bp_indices[i], nt_no, &pair);
        
        GetThisDpBpTag(seq_no, pair, tag);
        
        if (state)
        {
           sprintf(cmd, "%s itemconfigure %s -fill %s",
                        CANVAS_DP_NAME, tag, tk_color.sel_bp);
        }
        else
        {
            /* more correct is to check if this sequence is selected
             * if so, should use  COLOR_SEL_SEQ instead of unsel_bp
             * tcl keeps track of this
             */
            
            sprintf(cmd, "%s itemconfigure %s -fill %s",
                         CANVAS_DP_NAME, tag, tk_color.unsel_bp);
        }
        Tcl_SetStringObj(cmd_obj, cmd, -1);   
        
        #ifdef DEBUG    
            DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
        #endif
        
        if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
        {
            ERROR_P("Couldn't eval %s : %s\n",
                                           Tcl_GetStringFromObj(cmd_obj, NULL),
                                                   Tcl_GetStringResult(interp));
            return TCL_ERROR;
        }
        
        i++;
    }
    
     
    Tcl_DecrRefCount(cmd_obj);
    free(bp_indices);
    
    
    return TCL_OK;
}
/***   Highlight_DpBp_Cmd   ***/




/***   Highlight_AlnNt_From_DpConsBp_Cmd   *************************************
 *
 *
 */
int
Highlight_AlnNt_From_DpConsBp_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    int           state = 0;
    char         *state_str;
    cons_bpair_t *csbp;   
    dlist_elem   *contrib_list_elem ;
    int          *contrib_seqidx;
    bpair_t      *contrib_bp;
    int cs_nt_i,  cs_nt_j;
    Stopwatch_t  *watch;
    pair_t        pair;
    
    /*** check and get args
     *
     */    
    if ( (objc!=4))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?cs_nt_i cs_nt_j state?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &cs_nt_i)!=TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[2], &cs_nt_j)!=TCL_OK)
        return TCL_ERROR;
    state_str = Tcl_GetStringFromObj(objv[3], NULL);
    if (state_str == NULL )
    {
        CsDpError(interp, "can´t get state");
        return TCL_ERROR;
    }
    if (STR_EQ(state_str, "off"))
        state = 0;
    else if (STR_EQ(state_str, "on"))
        state = 1;

    #ifdef DEBUG
        DEBUG_P("%s\n", "start");
        DEBUG_P("cs_nt_i=%d\n",   cs_nt_i);
        DEBUG_P("cs_nt_j=%d\n",   cs_nt_j);
        DEBUG_P("state_str=%s\n", state_str);
    #endif
    

    /* c: zero-offset */
    cs_nt_i--;
    cs_nt_j--;
    
    watch = TimerStart();

    

    /*** now loop over all contributing basepairs
     *
     */
    SetPair(cs_nt_i, cs_nt_j, &pair);
    csbp = GetConsBp(aln, pair);

    /* PARANOIA */
    if ( ! dlist_size(csbp->contribs))
    {
        char buf[1024];
        sprintf(buf, "Upps...this consensus bp (%d:%d) got no contribs\n", cs_nt_i+1, cs_nt_j+1);
        CsDpError(interp, buf);
TMPDEBUG_P("%s\n", "simply delete it and don't bother?");
        return TCL_ERROR;
    }
    
    contrib_list_elem = dlist_head(csbp->contribs);
    while (1)
    {
        contrib_seqidx = (int*)dlist_data(contrib_list_elem);


        /*** if state==on then colorize via prob
         */      
        if (state)
        {
            char color[1024];
            int red   = 0xff;
            int green;
            int blue  = 0x00;
            
            contrib_bp = GetBpFromAlnSeq(aln, *contrib_seqidx, pair);
            if (contrib_bp==NULL)
            {
                ERROR_P("%s\n", "Uups..contrib is NULL");
                break;
            }
            green = (int)( (1.0 - contrib_bp->prob) * 255);
            

            sprintf(color, "#%06x", (red<<16) | (green<<8) | blue);

            Colorize_AlnNts(interp, *contrib_seqidx, pair, color);
        }
        /*** if state==off then decolorize
         */      
        else
        {
            Colorize_AlnNts(interp, *contrib_seqidx, pair, tk_color.unact_nt);
        }
        
        if (dlist_is_tail(contrib_list_elem))
            break;
             
        contrib_list_elem = dlist_next(contrib_list_elem);
    }
    
    
    TimerStop(__FUNCTION__, watch);
        
    return TCL_OK;
}
/***   Highlight_AlnNt_From_DpConsBp_Cmd   ***/





/***   Highlight_AlnNt_From_DpBp_Cmd   *****************************************
 *
 *
 * FIXME: Highligh only if visible
 * FIXME: This and Colorize_AlnNts should be Tcl coded
 */
int
Highlight_AlnNt_From_DpBp_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    int     nt_i, nt_j;
    int     seq_idx;
    int     state        = 0;
    char   *state_str;
    pair_t  pair;
    
    if ( (objc!=5))
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq_no nt_i nt_j state?");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[1], &seq_idx)!=TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[2], &nt_i)!=TCL_OK)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[3], &nt_j)!=TCL_OK)
        return TCL_ERROR;
    state_str = Tcl_GetStringFromObj(objv[4], NULL);
    if (state_str == NULL )
    {
        CsDpError(interp, "Can't get state");
        return TCL_ERROR;
    }
    if (STR_EQ(state_str, "off"))
    {
        state = 0;
    }
    else if (STR_EQ(state_str, "on"))
    {
        state = 1;
    }
    else
    {
        CsDpError(interp, "State must be one of (on|off)");
        return TCL_ERROR;
    }
    #ifdef DEBUG
        DEBUG_P("%s\n", "start");
        DEBUG_P("nt_i=%d\n",   nt_i);
        DEBUG_P("nt_j=%d\n",   nt_j);
        DEBUG_P("state_str=%s\n", state_str);
    #endif


   
    /* c is zero offset */
    seq_idx--;
    nt_i--;
    nt_j--;
    
    SetPair(nt_i, nt_j, &pair);
    if (state)
        Colorize_AlnNts(interp, seq_idx, pair, tk_color.act_nt);
    else
        Colorize_AlnNts(interp, seq_idx, pair, tk_color.unact_nt);
    
    return TCL_OK;
}
/***   Highlight_AlnNt_From_DpBp   ***/  




    
/***   Colorize_AlnNts   ***********************************************
 *
 * Highlights Alignments Nts nt_i and nt_j (zero-offset)
 * from seq seq_idx (zero-offset)
 * (where nt_i should be nt_j)
 */
int
Colorize_AlnNts(Tcl_Interp *interp, int seq_idx,
                        pair_t pair, char *color)
{
    char        cmd[1024];      /* used for tcl/tk command formatting */
    Tcl_Obj    *cmd_obj     = NULL;    /* used for executing tk commands     */
    
    cmd_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj);
    
    
    /*** colorize low nt = 5' nt left alnwindow
     *
     */
    sprintf(cmd, "%s itemconfigure", CANVAS_ALN_NAME );
    Tcl_SetStringObj(cmd_obj, cmd, -1);
	 
	 sprintf(cmd, " AlnNt_seq_%d_nt_%d -fill %s", seq_idx+1, pair.nti+1, color);
    Tcl_AppendToObj(cmd_obj, cmd, -1);
    
	 
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
            ERROR_P("Couldn't eval %s : %s\n",
                                           Tcl_GetStringFromObj(cmd_obj, NULL),
                                                   Tcl_GetStringResult(interp));
            return TCL_ERROR;
    }
    
    
    
    /*** colorize high nt = 3'nt left alnwindow
     *
     */
    sprintf(cmd, "%s itemconfigure", CANVAS_ALN_NAME  );
    Tcl_SetStringObj(cmd_obj, cmd, -1);
	 
	 sprintf(cmd, " AlnNt_seq_%d_nt_%d -fill %s", seq_idx+1, pair.ntj+1, color);
    Tcl_AppendToObj(cmd_obj, cmd, -1);

    
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj, NULL));
    #endif
    
    if (Tcl_EvalObjEx(interp, cmd_obj, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
        ERROR_P("Couldn't eval %s : %s\n",
                                       Tcl_GetStringFromObj(cmd_obj, NULL),
                                       Tcl_GetStringResult(interp));
        return TCL_ERROR;
      }
    
    
    Tcl_DecrRefCount(cmd_obj);

  
	 
	
 
	 
	 	  char        cmd3[1024];      /* used for tcl/tk command formatting */
    Tcl_Obj    *cmd_obj3     = NULL;    /* used for executing tk commands     */
    cmd_obj3 = Tcl_NewStringObj("THIS IS JUST A DUMMY\0", -1);
    Tcl_IncrRefCount(cmd_obj3);
    
    
    /*** colorize low nt = 5' nt right alnwindow
     *
     */
    sprintf(cmd3, "%s itemconfigure", CANVAS_ALN_NAME3);
    Tcl_SetStringObj(cmd_obj3, cmd3, -1);
    
    sprintf(cmd3, " AlnNt_seq_%d_nt_%d -fill %s", seq_idx+1, pair.nti+1, color);
    Tcl_AppendToObj(cmd_obj3, cmd3, -1);
    
    
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj3, NULL));
    #endif
    
    if (Tcl_EvalObjEx(interp, cmd_obj3, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
            ERROR_P("Couldn't eval %s : %s\n",
                                           Tcl_GetStringFromObj(cmd_obj3, NULL),
                                                   Tcl_GetStringResult(interp));
            return TCL_ERROR;
    }
    
    
    
    /*** colorize high nt = 3'nt right alnwindow
     *
     */
    sprintf(cmd3, "%s itemconfigure", CANVAS_ALN_NAME3);
    Tcl_SetStringObj(cmd_obj3, cmd3, -1);
    
    sprintf(cmd3, " AlnNt_seq_%d_nt_%d -fill %s", seq_idx+1, pair.ntj+1, color);
    Tcl_AppendToObj(cmd_obj3, cmd3, -1);
    
    
    #ifdef DEBUG
        DEBUG_P("Eval %s\n", Tcl_GetStringFromObj(cmd_obj3, NULL));
    #endif
    
    if (Tcl_EvalObjEx(interp, cmd_obj3, TCL_EVAL_GLOBAL) == TCL_ERROR)
    {
        ERROR_P("Couldn't eval %s : %s\n",
                                       Tcl_GetStringFromObj(cmd_obj3, NULL),
                                       Tcl_GetStringResult(interp));
        return TCL_ERROR;
      }
    
    
    Tcl_DecrRefCount(cmd_obj3);

    return TCL_OK;	 
	 
  }	 
	 

/***   Colorize_AlnNts   ***/






/***   GetBpPartnersForNt   **********************************************************
 *
 * Looks up bp partners for particular nt and returns
 * the ones, which form a bp with prob > 0
 * end element is -1
 *
 * Caller must free
 *
 */
int *
GetBpPartnersForNt(aln_t *ali, int seq_no, int nt_no)
{
    
    int         bp_partner = 0;        /* a possbile partner */
    int        *ret        = NULL;
    int         ret_idx;
    bpair_t    *bp;
    pair_t      pair;
       
    
    ret_idx=0;
    
    /*** lookup nt_no at 3'
     */
    for (bp_partner=0; bp_partner<nt_no; bp_partner++)
    {
        SetPair(bp_partner, nt_no, &pair);
        bp = GetBpFromAlnSeq(aln, seq_no, pair);
        if (bp==NULL)
            continue;
            
        if (bp->prob > 0.0)
        {
            ret_idx++;
            ret = realloc(ret, ret_idx * sizeof(int *));
            ret[ret_idx-1] = bp_partner;
        }
    }
    
    /*** lookup nt_no at 5'
     */
    for (bp_partner=nt_no+1; bp_partner<ali->len; bp_partner++)
    {
        SetPair(bp_partner, nt_no, &pair);
        bp = GetBpFromAlnSeq(aln, seq_no, pair);
        if (bp==NULL)
            continue;
        
        if (bp->prob > 0.0)
        {
            ret_idx++;
            ret = realloc(ret, ret_idx * sizeof(int *));
            ret[ret_idx-1] = bp_partner;
        }
    }

    /*** mark end
     */
    ret_idx++;
    ret = realloc(ret, ret_idx * sizeof(int));
    ret[ret_idx-1] = -1;

    return ret;
}
/***   GetBpPartnersForNt   ***/
