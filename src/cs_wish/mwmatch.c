/******************************************************************************
* 
* mwmatch.c - Routines for weighted matching
*            (stuff shared between imatch and bmatch)
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
 *  CVS $Id: mwmatch.c,v 1.15 2004-05-25 13:29:16 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>
#include <float.h>
#include <math.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "mx.h"
#include "cs_dpm.h"
#include "mwmatch.h"



/***   private
 */
#if 0
    #define DEBUG
#endif

#define TEST_HP_LOOP
#define MIN_HP_SIZE 4

static int
CountMaxEdges(TempGraph_t *g);
 
/*
 ***/






/***   GetCsGraph   ***********************************************************
 *
 * Read the cs_wish internally computed info content
 * and convert it to a graph
 *
 *
 * Return 1 on success, 0 otherwise
 */
int
GetCsGraph(TempGraph_t *g,  float prob_fac, float mic_fac,
                            float prob_thr, float mic_thr)
{
	int        i,j;
	float      idummy;
    float    **mpm_mx;
 
    mpm_mx    = CsDpm(prob_fac, mic_fac, prob_thr, mic_thr);
    
    if (mpm_mx==NULL)
    {
        ERROR_P("%s\n", "CsDpm failed");
        return 0;
    }
    #ifdef DEBUG
        DEBUG_P("Merged matrix (prob_fac=%f, mic_fac=%f)\n", (float)prob_fac, (float)mic_fac);
        PrintCsDpm(aln->cons_seq->nt, mpm_mx, stdout);
    #endif


	for (i=1; i<=g->n_vert; i++)
    {
		for (j=1; j<i; j++)
        {
            /* access new core */
			idummy = mpm_mx[i-1][j-1];
            
                       
            #ifdef TEST_HP_LOOP
                if (i-j < MIN_HP_SIZE)
                    idummy=0.0;
            #endif
            
            
            
            g->edges[i][j] = (int)(FLOAT2INT*idummy);
            #ifdef DEBUG
                DEBUG_P("g->edges[%d][%d]=%d (=%f * %d)\n", i, j, g->edges[i][j], idummy, (int)FLOAT2INT);
            #endif

            #ifdef PREVENT_MATCH_COREDUMP_BY_EDGE_INIT
                DEBUG_P("preventing coredump by setting g->edges[%d][%d]=1\n", i, j);
                if (g->edges[i][j]==0)
                    g->edges[i][j]=1;
            #endif
            
		}
	}
    
	return 1;
}
/***   GetCsGraph   ***/




/***   MakeTempGraph   *********************************************************
 *
 */
void
MakeTempGraph(TempGraph_t *g, int nvert)
{
	int size;
	int j;

	size = CountMaxEdges(g);

	g->edges    = (int **) Scalloc((nvert + 1), sizeof(int *));
    g->edges[2] = (int *)  Scalloc(size + 1, sizeof(int));    /* AW: ??? */

	for (j = 3; j <= nvert; j++)
		(g->edges)[j] = &((g->edges)[j - 1][j - 2]);

}
/***   MakeTempGraph   ***/




/***   CountMaxEdges   ********************************************************
 *
 */
int
CountMaxEdges(TempGraph_t *g)
{
	int x;

	x = g->n_vert;
	return (x*(x - 1))/2;
}
/***   CountMaxEdges   ***/




/***   FreeTempGraph   *********************************************************
 *
 */
void
FreeTempGraph(TempGraph_t *g)
{
	free(g->edges[2]);
	free(g->edges);
}
/***   FreeTempGraph   ***/

