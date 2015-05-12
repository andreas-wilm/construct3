/******************************************************************************
* 
* cs_bmatch.c - ConStruct interface to bmatch.c
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
 *  CVS $Id: cs_bmatch.c,v 1.18 2004-05-25 13:29:13 wilm Exp $    
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

#include "if.h"
#include "mwmatch.h"
#include "bmatch.h"

#include "cs_bmatch.h"



/***   private
 */

#if 0
    #define DEBUG
#endif


/* matching degree : don't change !*/
#define DEFAULT_B 2



static int
ReadDerivedGraph(TempGraph_t *g, float prob_fac, float mic_fac,
                                 float prob_thr, float mic_thr);

static char *
CreateBmatchResult(int **bm_out, TempGraph_t g, char *aln_seq);
 
/*
 ***/


/***   CS_Bmatch_Cmd   ******************************************************** 
 *
 */
int
CS_Bmatch_Cmd(ClientData clientData, Tcl_Interp *interp,
                              int objc, Tcl_Obj *CONST objv[])
{
    char *aln_seq, *basepair_arrayname;                 /* argument        */
	int  *Mate;
	int   k = 0;
	int   m = 0;
	int   j;
	int   pair = 0;
    int   pair_list[DEFAULT_B];
    float prob_list[DEFAULT_B];
	int **output;
	int   nvert;
	TempGraph_t g;
    char    *bmatch_str_result;
	Tcl_Obj *obj_result;
    double   prob_fac, prob_thr; 
    double   mic_fac,  mic_thr;   
    
    
	bmvar.B     = DEFAULT_B; /* matching degree */
	bmvar.IFlag = 1;         /* don't change !  */
 
    if (objc!=7)
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                         "?<aln_seq> <basepair_arrayname> <prob_fac> <mic_fac> <prob_thr> <mic_thr>?");
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



    #ifdef DEBUG
        DEBUG_P("aln_seq            = %s\n", aln_seq);
        DEBUG_P("basepair_arrayname = %s\n", basepair_arrayname);
        DEBUG_P("prob_fac           = %f\n", prob_fac);
        DEBUG_P("mic_fac            = %f\n", mic_fac);
        DEBUG_P("prob_thr           = %f\n", prob_thr);
        DEBUG_P("mic_thr            = %f\n", mic_thr);
    #endif




    nvert    = strlen(aln_seq);
	g.n_vert = nvert;
	g.n_edge = nvert*(nvert-1);


	MakeTempGraph(&g, nvert);

	if ( ! ReadDerivedGraph(&g, (float) prob_fac, (float) mic_fac,
                                (float) prob_thr, (float) mic_thr))
    {
        CsDpError(interp, "ReadDerivedGraph failed");
        return TCL_ERROR;
	}


	Mate = Weighted_Bmatch(1,1,k);


	output    = (int **) Scalloc((bmvar.OriginalVert) + 1, sizeof(int *));
	output[0] = (int *)  Scalloc((bmvar.B)*((bmvar.OriginalVert) + 1), sizeof(int));

	for ((bmvar.i) = 1; (bmvar.i) <= (bmvar.OriginalVert); (bmvar.i)++)
		output[(bmvar.i)] = output[(bmvar.i) - 1] + (bmvar.B);

	for ((bmvar.i) = (bmvar.FirstOuter) - 1; (bmvar.i) >= 1; (bmvar.i)--)
    {
		if (Mate[(bmvar.i)] >= (bmvar.FirstOuter))
        {
			for (j = 0; j < (bmvar.B); j++)
				if (output[(bmvar.FromVtx)[(bmvar.i)]][j] == 0)
					break;
            output[(bmvar.FromVtx)[(bmvar.i)]][j] = (bmvar.ToVtx)[(bmvar.i)];
		}
	}


	if (bmvar.IFlag)
    {
		for ((bmvar.i) = 1; (bmvar.i) <= (bmvar.OriginalVert); (bmvar.i)++)
        {
			for (j = 0; j < (bmvar.B); j++)
            {
				if ((pair = output[(bmvar.i)][j]) > 0)
					for (m = 0; m < (bmvar.B); m++)
						if ((bmvar.MatchStore)[(bmvar.i)][m] == pair)
							break;
				if (m == (bmvar.B))
                {
					output[(bmvar.i)][j] = -pair;
					for (m = 0; m < (bmvar.B); m++)
						if (output[pair][m] == (bmvar.i))
                        {
							output[pair][m] = -(bmvar.i);
							break;
						}
				}
			}
		}
	}
    #ifdef DEBUG
        DEBUG_P("%s", "Calculation finished exporting values to tcl\n");
    #endif


    
	/* write values to given basepair_arrayname
     */
    for(bmvar.i = 1; bmvar.i <= bmvar.OriginalVert; (bmvar.i)++)
    {       
	    for (j = 0; j < (bmvar.B); j++)   /* bmvar.B == matching degree DEFAULT_B */
        {

 			if ( bmvar.IFlag ) /* ? */
            {
				pair = (output[(bmvar.i)][j] < 0 ? -output[(bmvar.i)][j] : output[(bmvar.i)][j]);
              
                #ifdef PREVENT_MATCH_COREDUMP_BY_MEM_ACCESS_TEST
                    if (MAX(bmvar.i, pair)==1 && MIN(bmvar.i, pair)==0)
                    {
                        prob_list[j] = 0.0;
                        pair_list[j] = 0;           
                        DEBUG_P("%s\n", "MAX(bmvar.i, pair)==1 && MIN(bmvar.i, pair)==0 -> setting prob_list[j]=0.0 and pair_list[j] and cont");
                        continue;
                    }
                #endif                

				if (output[(bmvar.i)][j] < 0)
                {
					/* negative prob defines rematched vertex */
                    /* prob_list[j] = -1. * Edge(g, bmvar.i, pair)/FLOAT2INT; */
                    prob_list[j] = 0.0;
                    pair_list[j] = 0;
                }
                else
                {
					prob_list[j] = Edge(g, bmvar.i, pair)/FLOAT2INT;
                    pair_list[j] = pair;
				}

                /* ExportHigherBasepair uses zero-offset vars */
                pair_list[j]--;

            }
            
		} /* for (j = 0; j < (bmvar.B) */
        


        ExportHigherBasepair(interp, basepair_arrayname, (bmvar.i)-1,
                                     DEFAULT_B, pair_list, prob_list);
                        
	} /* for(bmvar.i = 1; bmvar.i <= bmvar.OriginalVert */



    /* write text output to tcl interpreter */
    bmatch_str_result = CreateBmatchResult(output, g, aln_seq);
    obj_result        = Tcl_NewStringObj(bmatch_str_result, -1);
    Tcl_SetObjResult(interp, obj_result);



    #ifdef DEBUG
        DEBUG_P("%s\n", "output");
        for (bmvar.i = 1; bmvar.i <= bmvar.OriginalVert; (bmvar.i)++)
        {
    		printf("%d: ", (bmvar.i));
    		for (j = 0; j < (bmvar.B); j++)
            {
    			printf("%d", (output[(bmvar.i)][j] < 0 ? -output[(bmvar.i)][j] : output[(bmvar.i)][j]));
    			if (output[(bmvar.i)][j] < 0)
    				printf("* ");
    			else
    				printf("  ");
    		}
    		printf("\n");
    	}
    
    	if ((bmvar.IFlag))
    		printf("\n* = Pairing involving a rematched vertex\n\n");
    #endif
    
    
    free(output[0]);
    free(output);
	FreeTempGraph(&g);
    FreeUpBmvar();
    free(bmatch_str_result);
    free(Mate);
    
	return TCL_OK;
}
/***   CS_Bmatch_Cmd   ***/




/***   CreateBmatchResult   ***************************************************
 *
 * Create a string which represents the bmatch output.
 * Caller must free
 * accesses bmvar
 *
 */
char *
CreateBmatchResult(int **bm_out, TempGraph_t g, char *aln_seq)
{
    const int chars_per_line = 54;    /* see sprintf below */    
    int   nt_idx;
    float prob[2];
    int   pair[2];
    char  rematch[2];
    char  pair_nt[2];
    int   j;
    char  *ret_buf;
    char  *line;
    
    ret_buf = (char *) Scalloc((strlen(aln_seq)+1)*chars_per_line, sizeof(char));
    ret_buf[0] = '\0'; /* strcat paranoia :) */
        
    if (bmvar.B!=2)
    {
        sprintf(ret_buf, "CreateBmatchResult only works for base triples (DEFAULT_B=2)");
        return ret_buf;
    }
   

    #ifdef DEBUG
        DEBUG_P("%s\n", "Creating string result");
    #endif

    line = (char *) Scalloc(chars_per_line, sizeof(char));
    
    for(bmvar.i = 1; bmvar.i <= bmvar.OriginalVert; (bmvar.i)++)
    {       
    
	    /* for (j = 0; j < (bmvar.B); j++) */
        
        if ( bmvar.IFlag ) /* ? */
        {
            nt_idx   = bmvar.i;
            
            j = 0;
            
            pair[0] = (bm_out[(bmvar.i)][j] < 0 ? -bm_out[(bmvar.i)][j] : bm_out[(bmvar.i)][j]);
            /* TMPDEBUG_P("bmvar.i=%d   pair[0]=%d\n", bmvar.i, pair[0]); */
            if (bm_out[(bmvar.i)][j] < 0)
                rematch[0] = '*';
            else
                rematch[0] = ' ';

              
            #ifdef PREVENT_MATCH_COREDUMP_BY_MEM_ACCESS_TEST
                if (MAX(bmvar.i, pair[0])==1 && MIN(bmvar.i, pair[0])==0)
                {
                    prob[0] = 0.0;
                    DEBUG_P("%s\n", "MAX(bmvar.i, pair[0])==1 && MIN(bmvar.i, pair[0])==0 -> setting prob[0]=0");
                }
                else
                {
                    prob[0] = ABS(Edge(g, bmvar.i, pair[0])/FLOAT2INT);
                }
            #else
                prob[0] = ABS(Edge(g, bmvar.i, pair[0])/FLOAT2INT);            
            #endif



            /* TMPDEBUG_P("bmvar.i=%d   prob[0]=%f\n", bmvar.i, prob[0]); */


            j = 1;
            
            pair[1] = (bm_out[(bmvar.i)][j] < 0 ? -bm_out[(bmvar.i)][j] : bm_out[(bmvar.i)][j]);
            
            /* TMPDEBUG_P("bmvar.i=%d   pair[1]=%d\n", bmvar.i, pair[1]); */
            if (bm_out[(bmvar.i)][j] < 0)
                rematch[1] = '*';
            else
                rematch[1] = ' ';
                
            #ifdef PREVENT_MATCH_COREDUMP_BY_MEM_ACCESS_TEST
                if (MAX(bmvar.i, pair[1])==1 && MIN(bmvar.i, pair[1])==0)
                {
                    prob[1] = 0.0;
                    VERBOSE_P("%s\n", "Circumvating bmatch error");
                    VERBOSE_P("%s\n", "Please increase the noise, i.e. reduce colormapping and increase mic-factor");
                }
                else
                {
                    prob[1] = -ABS(Edge(g, bmvar.i, pair[1])/FLOAT2INT);    
                }
            #else
                prob[1] = -ABS(Edge(g, bmvar.i, pair[1])/FLOAT2INT);
            #endif

            
            if (pair[0]!=0)
                pair_nt[0]=aln_seq[pair[0]-1];
            else
                pair_nt[0]=' ';

            if (pair[1]!=0)
                pair_nt[1]=aln_seq[pair[1]-1];
            else
                pair_nt[1]=' ';

            /* TMPDEBUG_P("bmvar.i=%d   pair_nt[0]=%c pair_nt[1]=%c\n", bmvar.i, pair_nt[0], pair_nt[1]); */
            
                           /* 2x (5 + 4 + 2 + 7) + 4 + 5 + 4 +1+1= 53*/            
            sprintf(line, "%4d (%c) %c %.4f <-> %4d (%c) <-> %4d (%c) %c %.4f\n",
                                pair[0], pair_nt[0], rematch[0], ABS(prob[0]),
                                nt_idx,  aln_seq[nt_idx-1],
                                pair[1], pair_nt[1], rematch[1], ABS(prob[1]));
            
            strcat(ret_buf, line);
        }
    }
    
    sprintf(line, "\n* = Pairing involving a rematched vertex\n\n");
    strcat(ret_buf, line);
    
    
    free(line);
    
    return ret_buf;
}
/***   CreateBmatchResult   ***/




/***   ReadDerivedGraph   *****************************************************
 *
 * data below threshold will be set to 0
 */
int
ReadDerivedGraph(TempGraph_t *g, float prob_fac, float mic_fac,
                                 float prob_thr, float mic_thr)
{

	int *row;
	int *degree;
	int *next_available;
	int alloc_size;
	int curr_edge;
	int curr_outer;
	int inner_needed = 0;
	int outer_needed = 0;
	int edges_needed = 0;
	int edge_mult;
	int this_weight;
	int v1, v2;
	int j, k, m;
    
    

	if (!GetCsGraph(g, (float)prob_fac, (float)mic_fac,
                       (float)prob_thr, (float)mic_thr))
    {
        ERROR_P("%s", "GetCsGraph failed\n");
		return 0;
	}

  

	bmvar.OriginalVert = g->n_vert;

	degree         = (int *) Scalloc(g->n_vert + 1, sizeof(int));
	next_available = (int *) Scalloc(g->n_vert + 1, sizeof(int));



	degree[0] = 0;
	next_available[0] = 1;
	edge_mult = 1 + 2*bmvar.B;
	bmvar.W = 0;

	for (j = 1; j <= g->n_vert; j++)
    {
		for (k = 1; k <= g->n_vert; k++)
			if ( (j != k) && ( (this_weight = Edge(*g, j, k)) > 0 ) )
            {
				bmvar.W = this_weight > bmvar.W ? this_weight : bmvar.W;
				degree[j]++;
				if (j < k)
					edges_needed++;
			}

		if (degree[j] > 0)
			outer_needed += bmvar.B;
		inner_needed += degree[j];
		next_available[j] = next_available[j - 1] + degree[j - 1];
	}

	edges_needed *= edge_mult;
	bmvar.W *= 4;
	bmvar.LAST_D = bmvar.W;

	bmvar.U = outer_needed + inner_needed;
	bmvar.FirstOuter = inner_needed + 1;
	bmvar.V = edges_needed;

	alloc_size = bmvar.U + 2*bmvar.V + 2;

	bmvar.A  	  = (int *) Scalloc(alloc_size, sizeof(int));
	bmvar.END	  = (int *) Scalloc(alloc_size, sizeof(int));
	bmvar.WEIGHT  = (int *) Scalloc(alloc_size, sizeof(int));
	bmvar.FromVtx = (int *) Scalloc(inner_needed, sizeof(int));
	bmvar.ToVtx   = (int *) Scalloc(inner_needed, sizeof(int));

	bmvar.MatchStore = (int **) Scalloc(bmvar.OriginalVert + 1, sizeof(int *));
	bmvar.MatchFlags = (char *) Scalloc(bmvar.B, sizeof(char));


	(bmvar.MatchStore)[0] = Scalloc((bmvar.OriginalVert + 1)*bmvar.B, sizeof(int));

	for (j = 1; j <= bmvar.OriginalVert; j++)
		(bmvar.MatchStore)[j] = (bmvar.MatchStore)[j - 1] + bmvar.B;

	bmvar.DUMMYVERTEX = bmvar.U + 1;
	bmvar.DUMMYEDGE	  = bmvar.U + 2*(bmvar.V) + 1;
	(bmvar.END)[bmvar.DUMMYEDGE] = bmvar.DUMMYVERTEX;

	bmvar.MATE     = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.LINK     = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.BASE     = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.NEXTVTX  = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.LASTVTX  = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.Y  	   = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.NEXT_D   = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.NEXTEDGE = (int *) Scalloc(bmvar.U + 2, sizeof(int));
	bmvar.NEXTPAIR = (int *) Scalloc(bmvar.U + 2 * bmvar.V + 2, sizeof(int));


	curr_edge = bmvar.U + 2;
	curr_outer = bmvar.FirstOuter;

	for (j = g->n_vert; j >= 1; j--)
    {
		if (degree[j] > 0) {
			row = g->edges[j];
			for (k = j - 1; k >= 1; k--)
			if (row[k] > 0)
            {
				v1 = (next_available[j])++;
				v2 = (next_available[k])++;
				bmvar.FromVtx[v1] = bmvar.ToVtx[v2] = j;
				bmvar.ToVtx[v1] = bmvar.FromVtx[v2] = k;

				this_weight = 4*row[k];
				(bmvar.Y)[v1] = (bmvar.Y)[v2] = -this_weight;

				(bmvar.WEIGHT)[curr_edge - 1] = (bmvar.WEIGHT)[curr_edge] = -2*this_weight;
				(bmvar.END)[curr_edge - 1]    = v1;
				(bmvar.END)[curr_edge]        = v2;
				(bmvar.A)[v1]                 = curr_edge;
				(bmvar.A)[v2]                 = curr_edge - 1;
				(bmvar.MATE)[v1]              = curr_edge;
				(bmvar.MATE)[v2]              = curr_edge - 1;

				curr_edge += 2;
			}

			for (k = 0; k < bmvar.B; k++)
            {
				for (m = 1; m <= degree[j]; m++)
                {
					v1 = next_available[j] - m;
					(bmvar.WEIGHT)[curr_edge - 1] = (bmvar.WEIGHT)[curr_edge] = 0;
					(bmvar.END)[curr_edge - 1]    = curr_outer;
					(bmvar.END)[curr_edge]        = v1;
					(bmvar.A)[curr_edge]          = (bmvar.A)[curr_outer];
					(bmvar.A)[curr_outer]         = curr_edge;
					(bmvar.A)[curr_edge - 1]      = (bmvar.A)[v1];
					(bmvar.A)[v1]                 = curr_edge - 1;
					curr_edge+= 2;
	 			}
				(bmvar.MATE)[curr_outer] = bmvar.DUMMYEDGE;
				(bmvar.Y)[curr_outer++]  = bmvar.W;
			}
		}
	}

	(bmvar.Y)[bmvar.DUMMYVERTEX]    = bmvar.W;
	(bmvar.MATE)[bmvar.DUMMYVERTEX] = bmvar.DUMMYEDGE;

	for (j = 1; j <= bmvar.U + 1; j++)
    {
		(bmvar.NEXTEDGE)[j] = bmvar.DUMMYEDGE;
		(bmvar.NEXTVTX)[j]  = 0;
		(bmvar.LINK)[j]     = -bmvar.DUMMYEDGE;
		(bmvar.BASE)[j]     = j;
		(bmvar.LASTVTX)[j]  = j;
		(bmvar.NEXT_D)[j]   = bmvar.LAST_D;
	}

    free(degree);
    free(next_available);
    
    return 1;
}
/***   ReadDerivedGraph   ***/
