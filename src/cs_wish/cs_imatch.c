/******************************************************************************
* 
* cs_imatch.c - ConStruct interface to imatch.c
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
 *  CVS $Id: cs_imatch.c,v 1.23 2007-10-22 10:43:23 steger Exp $    
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


#include "if.h"
#include "mwmatch.h"
#include "imatch.h"

#include "cs_imatch.h"



/***   private
 */
#if 0
    #define DEBUG
#endif

typedef struct  {
    int   pair;
    float prob;
} pair_prob_t;


static int
CountPositiveEdges(TempGraph_t *g);

static int
ReadGraph(TempGraph_t *g, float prob_fac, float mic_fac,
                          float prob_thr, float mic_thr);

static char *
CreateImatchResult(char *aln_seq, pair_prob_t *ipair, int iflag, int *Mate);

/*
 ***/

extern void
RemoveSingleBps(pair_prob_t *structure, int aln_len);



/***   CS_Imatch_Cmd   *********************************************************
 *
 * Runs imatch algorithm
 *
 * Argmuents:
 *  aln_seq:               aligned sequence string
 *  basepair_arrayname:    name of basepair array where data will be stored
 *  add_bpprob:            boolean: add csprob
 */
int
CS_Imatch_Cmd(ClientData clientData, Tcl_Interp *interp,
                              int objc, Tcl_Obj *CONST objv[])
{

	int     *Mate;
	int      j, nvert;
	float    v;
    char    *aln_seq, *basepair_arrayname;  /* argument        */
    char    *imatch_str_result;
    Tcl_Obj *obj_result;
	int      iflag = 1;   /* boolean: report rematched vertices for i-matching ? */
    int      k = 0;       /* Print maximum k-edge matching */
    double   prob_fac, prob_thr; 
    double   mic_fac,  mic_thr;   
	pair_prob_t *ipair;
	TempGraph_t  g;
    int      allow_single_bp;
    int      nt;
    
    
    if (objc!=8)
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                         "?<aln_seq> <basepair_arrayname> <prob_fac> <mic_fac> <prob_thr> <mic_thr> <allow_single>?");
        return TCL_ERROR;
    }

    aln_seq = Tcl_GetStringFromObj(objv[1], NULL);
    if (aln_seq == NULL)
    {   
        char buf[1024];
        sprintf(buf, "can't get aligned sequence string.");
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    basepair_arrayname = Tcl_GetStringFromObj(objv[2], NULL);
    if (basepair_arrayname == NULL)
    {
        char buf[1024];
        sprintf(buf, "can´t get basepair_arrayname.");
        CsDpError(interp, buf);
        return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, objv[3], &prob_fac)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[4], &mic_fac) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[5], &prob_thr)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[6], &mic_thr) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[7], &allow_single_bp)==TCL_ERROR)
        return TCL_ERROR;


    #ifdef DEBUG
        DEBUG_P("aln_seq            = %s\n", aln_seq);
        DEBUG_P("basepair_arrayname = %s\n", basepair_arrayname);
        DEBUG_P("prob_fac           = %f\n", prob_fac);
        DEBUG_P("mic_fac            = %f\n", mic_fac);
        DEBUG_P("prob_thr           = %f\n", prob_thr);
        DEBUG_P("mic_thr            = %f\n", mic_thr);
        DEBUG_P("allow_single_bp    = %d\n", allow_single_bp);
    #endif


    nvert = strlen(aln_seq);
	g.n_vert = nvert;
	g.n_edge = nvert*(nvert-1);


	MakeTempGraph(&g, nvert);
    
	if (!ReadGraph(&g, (float)prob_fac, (float)mic_fac,
                       (float)prob_thr, (float)mic_thr))
    {
        CsDpError(interp, "ReadGraph failed");
        return TCL_ERROR;
	}


	Mate = Weighted_Match(1,1,k);

    #ifdef DEBUG
        DEBUG_P("%s\n", "Mating done");
        for (j=1; j<=(imvar.U); j++)
            DEBUG_P("Mate[%d]=%d\n", j, Mate[j]);
    #endif

	ipair = (pair_prob_t *) Scalloc((1+nvert), sizeof(pair_prob_t));
	ipair[0].pair=0;
	ipair[0].prob=0;
    
	for (j=1; j<=(imvar.U); j++)
    {
        
		ipair[j].pair=Mate[j];
        /* Set prob to negative value if bp is rematched */
		if ((iflag) && ((imvar.REMATCH_FLAGS)[j] || (imvar.REMATCH_FLAGS)[Mate[j]]))
            v=-1.0;
        else
			v=1.0;
        
        #ifdef PREVENT_MATCH_COREDUMP_BY_MEM_ACCESS_TEST
            if (j==1 && Mate[j]==0)
            {
                ipair[j].prob=0.0;             
                DEBUG_P("%s\n", "j==1 && Mate[j]==0 ! -> setting ipair[j].prob=0.0 and cont");
                continue;
            }
        #endif
        
		if (j > Mate[j])
        {
			ipair[j].prob=(float)g.edges[j][Mate[j]] / FLOAT2INT * v;
            #ifdef DEBUG                
			    DEBUG_P("ipair[j].prob=%f, g.edges[%d][%d]=%d\n",
                                 ipair[j].prob, j, Mate[j], g.edges[j][Mate[j]]);
            #endif
		}
        else
        {
			ipair[j].prob=(float)g.edges[Mate[j]][j] / FLOAT2INT *v ;
            #ifdef DEBUG
                DEBUG_P("ipair[j]prob=%f, g.edges[%d][%d]=%d\n",
                                 ipair[j].prob, Mate[j],j, g.edges[Mate[j]][j]);
            #endif
		}
	}
    

	/* write values to given basepair_arrayname
     */
    for(nt=1; nt<=imvar.U; nt++) /* == g.n_vert == strlen(aln_seq) */
    {
        float prob;
        int   nt_pair;
        
		if ( (iflag) && ((imvar.REMATCH_FLAGS)[nt] || (imvar.REMATCH_FLAGS)[Mate[nt]]) )
        {
	   	    nt_pair = 0;
            prob    = 0.0;
		}
        else
        {
            nt_pair = ipair[nt].pair;
            prob    = ipair[nt].prob;
		}
    }

    /*****   remove single basepairs if requested */
    if (!allow_single_bp)
        RemoveSingleBps(ipair, imvar.U);

    for(nt=1; nt<=imvar.U; nt++) /* == g.n_vert == strlen(aln_seq) */
    {
        float prob;
        int   nt_pair;

        nt_pair = ipair[nt].pair;
        prob    = ipair[nt].prob;
        ExportBasepair(interp, basepair_arrayname, nt-1, nt_pair-1, prob);
	}

    imatch_str_result = CreateImatchResult(aln_seq, ipair, iflag, Mate);
    
    
    
    obj_result = Tcl_NewStringObj(imatch_str_result, -1);
    Tcl_SetObjResult(interp, obj_result);
    

/*	AW FIXME crashed sometimes: FreeTempGraph(&g);*/
    FreeUpImvar();
    free(imatch_str_result);

	return TCL_OK;
}
/***   CS_Imatch_Cmd   ***/




/***   ReadGraph   *************************************************************
 *
 * return 1 on success, 0 otherwise
 * data below threshold will be set to 0
 */
int
ReadGraph(TempGraph_t *g, float prob_fac, float mic_fac,
                          float prob_thr, float mic_thr)
{
	int *row;
	int alloc_size;
	int curr_edge;
	int j, k;

	if (!GetCsGraph(g, prob_fac, mic_fac, prob_thr, mic_thr))
		return 0;

	
	imvar.U = g->n_vert;
	imvar.V = CountPositiveEdges(g);
	alloc_size = imvar.U + 2*imvar.V + 2;

	imvar.A             = (int *)  Scalloc(alloc_size, sizeof(int));
	imvar.END           = (int *)  Scalloc(alloc_size, sizeof(int));
	imvar.WEIGHT        = (int *)  Scalloc(alloc_size, sizeof(int));
	imvar.REMATCH_FLAGS = (char *) Scalloc((imvar.U) + 1, sizeof(char));

	curr_edge = (imvar.U) + 2;

	for (j = (imvar.U); j > 1; j--)
    {
		row = g->edges[j];
		for (k = j - 1; k >= 1; k--)
        {
			if (row[k] > 0) {
				(imvar.WEIGHT)[curr_edge - 1] = (imvar.WEIGHT)[curr_edge] = 4*row[k];
				(imvar.END)[curr_edge - 1]    = j;
				(imvar.END)[curr_edge]        = k;
				(imvar.A)[curr_edge]          = (imvar.A)[j];
				(imvar.A)[j]                  = curr_edge;
				(imvar.A)[curr_edge - 1]      = (imvar.A)[k];
				(imvar.A)[k]                  = curr_edge - 1;
				curr_edge += 2;
			}
	    }
    }
    
    return 1;
}
/***   ReadGraph   ***/




/***   CountPositiveEdges   ***************************************************
 *
 */
int
CountPositiveEdges(TempGraph_t *g)
{
	int j, k, count = 0;

	for (j = 2; j <= g->n_vert; j++)
		for (k = 1; k < j; k++)
			if (Edge(*g, j, k) > 0)
		count++;

	return count;
}
/***   CountPositiveEdges   ***/




/***   CreateImatchResult   ***************************************************
 *
 * Create a string which represents the imatch output.
 * Caller must free
 *
 */
char *
CreateImatchResult(char *aln_seq,      /* sequence */
                   pair_prob_t *ipair, /* base pairs */
                   int iflag,
                   int *Mate)
{
	int    i, j1, j2=0, aln_len;
	char   im1, im2=' ', nt1, nt2=' ', line[66];
    char   *ret_imatch_result;
    
	aln_len = strlen(aln_seq);    
    ret_imatch_result = (char*) Scalloc(aln_len*66+1, sizeof(char));
    ret_imatch_result[0] = '\0'; /* strcat paranoia :) */
    

	for(i=1; i<=aln_len; i+=2)
    {
		j1 = ipair[i].pair;
		if ( (iflag)  &&  ((imvar.REMATCH_FLAGS)[i] || (imvar.REMATCH_FLAGS)[Mate[i]]) )
			im1='*';
        else
			im1=' ';

		if (j1 != 0)
			nt1=*(aln_seq+j1-1);
		else
			nt1=' ';


		if (i < aln_len)
        {
			j2 = ipair[i+1].pair;
			if ( (iflag) && ((imvar.REMATCH_FLAGS)[i+1] || (imvar.REMATCH_FLAGS)[Mate[i+1]]) )
				im2='*';
			else
				im2=' ';

			if (j2 != 0)
				nt2=*(aln_seq+j2-1);
			else
				nt2=' ';
		}
		
		if (i < aln_len)
        {
			/*                      1   23  456789   01  23  4     56789   01  234567   89  01  2*/
            /* 32+(4+1+4+1+1+6)*2=66 */
			sprintf(line,     "\n%4d (%c) <-> %4d (%c) %c %1.4f   | %4d (%c) <-> %4d (%c) %c %1.4f",
		                                   i,   *(aln_seq+i-1), j1, nt1, im1, ABS(ipair[i].prob),
                                           i+1, *(aln_seq+i),   j2, nt2, im2, ABS(ipair[i+1].prob));

		} else {
			sprintf(line,     "\n%4d (%c) <-> %4d (%c) %c %1.4f   |", 
		         i, *(aln_seq+i-1), j1, nt1, im1, ABS(ipair[i].prob));
		}
		
				strcat(ret_imatch_result, line);
	}

    return ret_imatch_result;
}
/***   CreateImatchResult   ***/
