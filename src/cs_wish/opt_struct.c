/******************************************************************************
* 
* opt_struct.c - dynamic programming function to compute an optimal structure
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
 *  CVS $Id: opt_struct.c,v 1.26 2004-05-25 13:29:17 wilm Exp $    
 */


/**********
 * 
 *  Algorithm is based upon Nussinov et al 1978 and Waterman & Smith 1978
 *  
 *  See "Biological Sequence Analysis" Durbin, Eddy, Krogh, Mitchison
 *  for example
 *
**********/


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE
  
#include <tcl.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "if.h"
#include "mx.h"
#include "cs_dpm.h"
#include "struct.h"
#include "stack.h"

#include "opt_struct.h"



/***   private
 */

#if 0
    #define DEBUG
#endif


static void
Optimize(int seqlen, int **dp, int **opt);

static int *
Backtrack(int seqlen, int **dp, int **opt);
 
/*
 ***/






/***   Backtrack   ************************************************************
 *
 * Implementation of algorithm from
 * "Biological Sequence Analysis";  Durbin & Eddy & Krogh & Mitchison
 * Cambridge University Press; 1998
 *
 * Returns an int array [0..seqlen-1] with corresponding basepair partners
 * for each indexed nt or -1 for single stranded / unpaired
 *
 * Note:
 *  push: first nt_i then nt_j
 * 	pop : first nt_j then nt_i
 *
 */
int *
Backtrack(int seqlen, int **dp, int **opt)
{
    int *bps;
    int k;
    stack_t *s;
    int *push_data;
    int *nt_i, *nt_j;
    
    
    if ((bps = (int*) calloc(seqlen, sizeof(int)))==NULL)
    {
        ERROR_P("%s\n", "no more memory");
        exit(-1);
    }

    for (k=0; k<seqlen; k++)
        bps[k] = -1;
    
    s = stack_new(sizeof(int*), free);
	

    push_data=(int*)calloc(1, sizeof(int));
    *push_data=0;
    stack_push(s, push_data);

    
    push_data=(int*)calloc(1, sizeof(int));
    *push_data=seqlen-1;
    stack_push(s, push_data);


	while (stack_size(s))
    {	
		nt_j=(int*)stack_pop(s, 0);
	    nt_i=(int*)stack_pop(s, 0);
        
		#ifdef DEBUG
            DEBUG_P("invastigating %d:%d\n", *nt_i, *nt_j);
        #endif
        
        /*
		  order of checks arbitrary
		  downleft (incl. check for bp) + left + down  = steger's version
		 */

		if (*nt_i>=*nt_j)
		{
        	continue;    
		
		/* down */
		}
        else if (opt[*nt_j][*nt_i+1]==opt[*nt_j][*nt_i])
        {
            push_data=(int*)calloc(1, sizeof(int));
            *push_data=*nt_i+1;
            stack_push(s, push_data);

			
            push_data=(int*)calloc(1, sizeof(int));
            *push_data=*nt_j;
            stack_push(s, push_data);

			
		/* left */
		}
        else if (opt[*nt_j-1][*nt_i] == opt[*nt_j][*nt_i])
        {
            push_data=(int*)calloc(1, sizeof(int));
            *push_data=*nt_i;
            stack_push(s, push_data);


            push_data=(int*)calloc(1, sizeof(int));
            *push_data=*nt_j-1;
            stack_push(s, push_data);

		
        /* downleft */
		}
        else if ( (opt[*nt_j-1][*nt_i+1] + dp[*nt_j][*nt_i]) == opt[*nt_j][*nt_i])
        {
            #ifdef DEBUG
                DEBUG_P("basepair: %d:%d\n", *nt_i, *nt_j);
            #endif
            
			bps[*nt_i]=*nt_j;
			bps[*nt_j]=*nt_i;
			
            push_data=(int*)calloc(1, sizeof(int));
            *push_data=*nt_i+1;
            stack_push(s, push_data);


            push_data=(int*)calloc(1, sizeof(int));
            *push_data=*nt_j-1;
            stack_push(s, push_data);


		/* bifurcation */
		}
        else
        {
			for (k=*nt_i+1; k<=*nt_j-1; k++)
            {
				if ( (opt[k][*nt_i] + opt[*nt_j][k+1]) == opt[*nt_j][*nt_i])
                {
                    push_data=(int*)calloc(1, sizeof(int));
                    *push_data=k+1;
                    stack_push(s, push_data);


                    push_data=(int*)calloc(1, sizeof(int));
                    *push_data=*nt_j;
                    stack_push(s, push_data);


                    push_data=(int*)calloc(1, sizeof(int));
                    *push_data=*nt_i;
                    stack_push(s, push_data);


                    push_data=(int*)calloc(1, sizeof(int));
                    *push_data=k;
                    stack_push(s, push_data);
			
					k=*nt_j;/* =break */
				}
			
			}
		}
        free(nt_j);
        free(nt_i);
	} /* while (stack_size(s)) */
    
    stack_free(s);
	
    return bps;
}
/***   Backtrack   ***/




/***   Optimize   *************************************************************
 *
 * maximizes number of basepairs in matrix opt
 * based upon Nussinov et al 1978 and Waterman & Smith 1978
 */
void
Optimize(int seqlen, int **dp, int **opt)
{
    int i, j, k, max;
    
    
	for (j=0; j<seqlen; j++)
    {
		for (i=0; i<j-MIN_HP_SIZE; i++)
        {
			max=0;
			for (k=i; k<j-MIN_HP_SIZE; k++)
            {

				if (dp[j][k] > 0)
                {
                    if (k-1<0)
                        max = MAX(max, 0 + opt[j-1][k+1] + dp[j][k]);
                    else
                        max = MAX(max, opt[k-1][i] + opt[j-1][k+1] + dp[j][k]);                    
                }                   
            }
			opt[j][i] = MAX(opt[j-1][i], max);
            #ifdef DEBUG
                DEBUG_P("set opt[%d][%d] to :%d\n", j, i, opt[j][i]);
	        #endif
    	}
	}
}
/***   Optimize   ***/




/***   GenOptConStructure   ***************************************************
 *
 * takes float matrix mcm and computes the best structure via
 * dynamic programming
 * the matrix will be converted to ints before to vaoid
 * rounding problems while backtracking
 *
 */
pair_prob_t *
GenOptConStructure(float **mcm, int alnlen)
{
    int         **dyntable;
    float         to_int_factor;
    int         **intdp;
    float         tmp;
    int           *bps;
    pair_prob_t   *pairs;    
    int i, j;    
    
    

    intdp   = (int**) mx_new(alnlen, 0, sizeof(int), sizeof(int*), MX_DP);
    /* some elements of dyntable accessed while backtracking
       are not within a dotplot, therefore a standard matrix
       (MX_DEFAULT) is needed */
	dyntable  = (int**) mx_new(alnlen, alnlen, sizeof(int), sizeof(int*), MX_DEFAULT);



   
    
    /* old norm factor
     *    Convert probabilities from float to integers
     *    max. number of bp is aligned_len/2;
     *    each bp has a max. probability of 1.;
     * norm = 2.*(float)INT_MAX/(float)aln_len;
     */
    to_int_factor = FLOAT2INT;
    /*
    to_int_factor = 2.*(float)INT_MAX/(float)alnlen;
    */

    for (i=0; i<aln->len; i++)
		for (j=i+1; j<aln->len; j++)
            intdp[j][i] = (int) (to_int_factor * mcm[j][i]);
    


    /***   filling stage
     */
	Optimize(alnlen, intdp, dyntable);


    /***   backtrack
     */
    bps = Backtrack(alnlen, intdp, dyntable);

    
    pairs = (pair_prob_t *) Scalloc(alnlen, sizeof(pair_prob_t));
     

    for(i=0; i<alnlen; i++)
    {
        if (bps[i]!=-1)
        {
            if (i>bps[i])
                tmp = (float) (intdp[i][bps[i]] / to_int_factor);
            else
                tmp = (float) (intdp[bps[i]][i] / to_int_factor);

                
            pairs[i].pair      = bps[i];
            pairs[bps[i]].pair = i;
            
            pairs[i].prob      = tmp;
            pairs[bps[i]].prob = tmp;
        }
        else
        {
            tmp = -1.0;
            
            pairs[i].pair      = -1;
            pairs[i].prob      = tmp;
        }
/*        #ifdef DEBUG*/
            DEBUG_P("setting %d:%d -> %f\n", i, bps[i], tmp);
/*        #endif*/
    }
    
    mx_free((void**)intdp, alnlen, MX_DP);    
    mx_free((void**)dyntable, alnlen, MX_DEFAULT);
    free(bps);

    return (pairs);
}
/***   GenOptConStructure   ***/




/***   OptimalStructure_Cmd   *************************************************
 *
 * calculate only best structure, incl. backtrack
 *
 * return array with name basepair_arrayname (arg)
 *  index1: 1-alignmentlength
 *  index2: BP_IDX2_PARTNER|BP_IDX2_PROB,
 *      where partner=nt-index and prob=probability
 *
 */
int
OptimalStructure_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    pair_prob_t *pair;
    int          nt_idx=0;
    float      **mcm; /* the new merged consensus matrix, merged out of mic and td.prob */
    /* arguments */
    char        *basepair_arrayname;
    double       prob_fac, prob_thr; 
    double       mic_fac, mic_thr;   
    int          allow_single_bp;
    
    /***   parse args
     *
     */
    if (objc!=7)
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                         "?<pair_array_name> <prob_fac> <mic_fac> <prob_thr> <mic_thr> <allow_single_bp>?");
        return TCL_ERROR;
    }
    basepair_arrayname = Tcl_GetStringFromObj(objv[1], NULL);
    if (basepair_arrayname == NULL)
    {
        CsDpError(interp, "Can´t get basepair_arrayname");
        return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, objv[2], &prob_fac)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &mic_fac) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[4], &prob_thr)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[5], &mic_thr) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[6], &allow_single_bp)==TCL_ERROR)
        return TCL_ERROR;

    #ifdef DEBUG 
        DEBUG_P("pair_array_name=%s\n", basepair_arrayname);
        DEBUG_P("prob_fac=%f\n", prob_fac);
        DEBUG_P("mic_fac=%f\n",  mic_fac);
        DEBUG_P("prob_thr=%f\n", prob_thr);
        DEBUG_P("mic_thr=%f\n",  mic_thr);
        DEBUG_P("allow_single_bp=%d\n", allow_single_bp);
    #endif
    
    
    mcm = CsDpm((float) prob_fac, (float) mic_fac,
                (float) prob_thr, (float) mic_thr);
    
    if (mcm==NULL)
    {
        CsDpError(interp, "CsDpm failed");
        return TCL_ERROR;
    }
    #ifdef DEBUG
        DEBUG_P("Merged matrix (prob_fac=%f, mic_fac=%f)\n", (float)prob_fac, (float)mic_fac);
        PrintCsDpm(aln->cons_seq->nt, mcm, stdout);
    #endif
    

    pair = GenOptConStructure(mcm, aln->len);

    /*****   remove single basepairs if requested
     */
    if (!allow_single_bp)
        RemoveSingleBps(pair, aln->len);

    
    /*****   export calculated structure
     */
    for(nt_idx=0; nt_idx<aln->len; nt_idx++)
        ExportBasepair(interp, basepair_arrayname, nt_idx, pair[nt_idx].pair, pair[nt_idx].prob);
     
    free(pair);

    FreeCsDpm(mcm);
    
    return TCL_OK;
}
/***   OptimalStructure_Cmd   ***/


