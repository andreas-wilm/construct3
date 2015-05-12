/******************************************************************************
* 
* bmatch.c - routines for maximum weigthed matching
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
 *  CVS $Id: bmatch.c,v 1.10 2004-05-25 13:29:12 wilm Exp $    
 */



/***   Routines for Maximum Weigthed Matching    ******************************
 *
 * N-cubed weighted matching
 * Implementation of H. Gabow's Ph.D. thesis, Stanford Univ. 1973
 * Written by Edward Rothberg  7/85
 * For complete details, please refer to the original paper
 *
 * Modified by J. Riks for use in construct. See Tabaska, Bioinformatics
 * Vol. 14 no. 8 1998 (691-699) for further details
 *
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"


#include "bmatch.h"



/***   private
 */

#if 0
    #define DEBUG
#endif

static void
CheckRematchB();
static void
BUNPAIR();
static void
BREMATCH ();
static void
BUNLINK ();
static void
BSET_BOUNDS ();
static void
BUNPAIR_ALL ();
static void
BPOINTER ();
static void
BSCAN ();
static void
BPAIR ();
static void
BMERGE_PAIRS();
static void
BLINK_PATH();
static void
BINSERT_PAIR();
 	  
/*
 ***/



/***   Weighted_Bmatch   ******************************************************
 *
 */
int *Weighted_Bmatch (int type, int maximize, int k)
{

	int g, j, w, outcome, loop=0;



	while (1)
    {
		loop++;
		if (k && (loop > k))
			break;
	
		bmvar.DELTA = 0;
		for (bmvar.v=1; bmvar.v<=bmvar.U; ++bmvar.v)
			if ((bmvar.MATE)[bmvar.v] == bmvar.DUMMYEDGE)
				BPOINTER (bmvar.DUMMYVERTEX, bmvar.v, bmvar.DUMMYEDGE);

		
        while (1)
        {
			bmvar.i = 1;
			for (j=2; j<=bmvar.U; ++j)
				if ((bmvar.NEXT_D)[bmvar.i] > (bmvar.NEXT_D)[j])
					bmvar.i = j;
			bmvar.DELTA = (bmvar.NEXT_D)[bmvar.i];
			if (bmvar.DELTA == bmvar.LAST_D)
				goto done;
			bmvar.v = (bmvar.BASE)[bmvar.i];
			if ((bmvar.LINK)[bmvar.v] >= 0)
            {
				BPAIR (&outcome);
				if (outcome == 1)
					break;
			}
            else
            {
				w = BMATE (bmvar.v);
				if ((bmvar.LINK)[w] < 0)
					BPOINTER (bmvar.v, w, OPPEDGE((bmvar.NEXTEDGE)[bmvar.i]));
				else
                    BUNPAIR (bmvar.v, w);
			}
		}
	
		bmvar.LAST_D -=bmvar.DELTA;
		BSET_BOUNDS();
		g = OPPEDGE(bmvar.e);
		BREMATCH (BEND(bmvar.e), g);
		BREMATCH (BEND(g), bmvar.e);

		if (bmvar.IFlag)
			CheckRematchB();
	}


	done:
    
    #ifdef DEBUG
        DEBUG_P("Final matching: %d Edges\n", loop - 1);
    #endif
    
	BSET_BOUNDS();
	BUNPAIR_ALL();

	if (bmvar.IFlag)
		CheckRematchB();

	for (bmvar.i=1; bmvar.i<=bmvar.U;++(bmvar.i))
    {
		(bmvar.MATE)[bmvar.i] = (bmvar.END)[(bmvar.MATE)[bmvar.i]];
		if ((bmvar.MATE)[bmvar.i]==bmvar.DUMMYVERTEX)
			(bmvar.MATE)[bmvar.i]=0;
	}

	return(bmvar.MATE);
}
/***   Weighted_Bmatch   ***/




/***   FreeUpBmvar   **********************************************************
 *
 */
void FreeUpBmvar()
{
	free(bmvar.LINK);
	free(bmvar.BASE);
	free(bmvar.NEXTVTX);
	free(bmvar.LASTVTX);
	free(bmvar.Y);
	free(bmvar.NEXT_D);
	free(bmvar.NEXTEDGE);
	free(bmvar.NEXTPAIR);

	free(bmvar.A);
	free(bmvar.END);
	free(bmvar.WEIGHT);


	free(bmvar.FromVtx);
	free(bmvar.ToVtx);
	free((bmvar.MatchStore)[0]);
	free(bmvar.MatchStore);
	free(bmvar.MatchFlags);

}
/***   FreeUpBmvar   ***/




/***   CheckRematchB   ********************************************************
 * 
 * Clunky function to detect rematches
 *
 */
void CheckRematchB()
{
	int from;
	int to;
	int count;
	int this_vert;
	int j, k, m;

	for (j = bmvar.FirstOuter, from = bmvar.OriginalVert; j <= bmvar.U; j += bmvar.B, from--)
    {
		for (k = 0; k < bmvar.B; k++)
			(bmvar.MatchFlags)[k] = 0;
		count = 0;

		for (k = 0; k < bmvar.B; k++)
        {
			this_vert = j + k;

		 	if (((bmvar.MATE)[this_vert] != bmvar.DUMMYEDGE) &&
		 	    ((bmvar.FromVtx)[(bmvar.END)[(bmvar.MATE)[this_vert]]] == from))
            {
				to = (bmvar.ToVtx)[(bmvar.END)[(bmvar.MATE)[this_vert]]];
				for (m = 0; m < bmvar.B; m++)
                {
					if ((bmvar.MatchStore)[from][m] == 0)
						(bmvar.MatchStore)[from][m] = to;

					if ((bmvar.MatchStore)[from][m] == to)
                    {
						(bmvar.MatchFlags)[m] = 1;
						count++;
						break;
					}
				}
			}
            else
            {
				count++;
		    }
        }

		if (count < bmvar.B) {
			for (k = 0; k < bmvar.B; k++)
				if ((bmvar.MatchFlags)[k] == 0)
                {
					(bmvar.MatchStore)[from][k] = -1;
					break;
				}
		}
	}
}
/***   CheckRematchB   ***/




/***   BUNPAIR   **************************************************************
 *
 * Expands a blossom.  Fixes up LINK and MATE.
 *
 */
void BUNPAIR (int oldbase, int oldmate)
{
	int e, newbase, u;

    #ifdef DEBUG
	    printf("Unpair oldbase, oldmate=%d %d\n",oldbase, oldmate);
    #endif

	BUNLINK (oldbase);
	newbase = BMATE (oldmate);
	if (newbase != oldbase)
    {
		(bmvar.LINK)[oldbase] = -bmvar.DUMMYEDGE;
		BREMATCH (newbase, (bmvar.MATE)[oldbase]);
		if (bmvar.f == (bmvar.LASTEDGE)[1])
			(bmvar.LINK)[bmvar.secondmate] = -(bmvar.LASTEDGE)[2];
		else (bmvar.LINK)[bmvar.secondmate] = -(bmvar.LASTEDGE)[1];
	}
	e = (bmvar.LINK)[oldmate];
	u = BEND (OPPEDGE (e));
	if (u == newbase)
    {
		BPOINTER (newbase, oldmate, e);
		return;
	}
	(bmvar.LINK)[BMATE (u)] = -e;
	do
    {
		e = -(bmvar.LINK)[u];
		bmvar.v = BMATE (u);
		BPOINTER (u, bmvar.v, -(bmvar.LINK)[bmvar.v]);
		u = BEND (e);
	} while (u != newbase);
	e = OPPEDGE (e);
	BPOINTER (newbase, oldmate, e);
}
/***   BUNPAIR   ***/




/***   BREMATCH   *************************************************************
 *
 * changes the matching along an alternating path
 * firstmate is the first base vertex on the path
 * edge e is the new matched edge for firstmate
 *
 */
void BREMATCH (int firstmate, int e)
{

    #ifdef DEBUG
    	printf("Rematch firstmate=%d e=%d-%d\n", firstmate,
                                   (bmvar.END)[OPPEDGE(e)],
                                            (bmvar.END)[e]);
    #endif
  
	(bmvar.MATE)[firstmate] = e;
    bmvar.nexte = -(bmvar.LINK)[firstmate];

	while (bmvar.nexte != bmvar.DUMMYEDGE)
    {
		e                = bmvar.nexte;
		bmvar.f          = OPPEDGE (e);
		firstmate        = BEND (e);
		bmvar.secondmate = BEND (bmvar.f);

		bmvar.nexte                    = -(bmvar.LINK)[firstmate];
		(bmvar.LINK)[firstmate]        = -(bmvar.MATE)[bmvar.secondmate];
		(bmvar.LINK)[bmvar.secondmate] = -(bmvar.MATE)[firstmate];
		(bmvar.MATE)[firstmate]        = bmvar.f;
		(bmvar.MATE)[bmvar.secondmate] = e;
	}
}
/***   BREMATCH   ***/




/***   BUNLINK   **************************************************************
 *
 * unlinks subblossoms in a blossom.  oldbase is the base of the blossom to
 * be unlinked.
 *
 */
void BUNLINK (int oldbase)
{
	int k, j=1;

    #ifdef DEBUG
    	printf("Unlink oldbase=%d\n",oldbase);
    #endif
	
	bmvar.i = (bmvar.NEXTVTX)[oldbase];
	bmvar.newbase = (bmvar.NEXTVTX)[oldbase];
	bmvar.nextbase = (bmvar.NEXTVTX)[(bmvar.LASTVTX)[bmvar.newbase]];
	bmvar.e = (bmvar.LINK)[bmvar.nextbase];
    
UL2:
	do
    {
	    bmvar.nextedge = OPPEDGE ((bmvar.LINK)[bmvar.newbase]);
    	for (k=1; k <= 2; ++k)
        {
    		(bmvar.LINK)[bmvar.newbase] = -(bmvar.LINK)[bmvar.newbase];
    		(bmvar.BASE)[bmvar.i] = bmvar.newbase;
    		bmvar.i = (bmvar.NEXTVTX)[bmvar.i];
    		while (bmvar.i != bmvar.nextbase)
            {
    		    (bmvar.BASE)[bmvar.i] = bmvar.newbase;
    		    bmvar.i = (bmvar.NEXTVTX)[bmvar.i];
    		}
    		bmvar.newbase = bmvar.nextbase;
    		bmvar.nextbase = (bmvar.NEXTVTX)[(bmvar.LASTVTX)[bmvar.newbase]];
        }
	} while ((bmvar.LINK)[bmvar.nextbase] == bmvar.nextedge);
    
	if (j==1)
    {
		(bmvar.LASTEDGE)[1] = bmvar.nextedge;
		j++;
		bmvar.nextedge = OPPEDGE (bmvar.e);
		if ((bmvar.LINK)[bmvar.nextbase] == bmvar.nextedge)
			goto UL2;
	}
	(bmvar.LASTEDGE)[2] = bmvar.nextedge;

	if ((bmvar.BASE)[(bmvar.LASTVTX)[oldbase]] == oldbase)
    {
		(bmvar.NEXTVTX)[oldbase] = bmvar.newbase;
	}
    else
    {
		(bmvar.NEXTVTX)[oldbase] = bmvar.DUMMYVERTEX;
		(bmvar.LASTVTX)[oldbase] = oldbase;
	}
}
/***   BUNLINK   ***/




/***   BSET_BOUNDS   **********************************************************
 *
 *
 * updates numerical bounds for linking paths.
 * called with LAST_D set to the bound on DELTA for the next search
 *
 */
void BSET_BOUNDS ()
{
	int del;

	for (bmvar.v=1; bmvar.v <= bmvar.U; ++(bmvar.v))
    {
		if ((bmvar.LINK)[bmvar.v] < 0 || (bmvar.BASE)[bmvar.v] != bmvar.v)
        {
			(bmvar.NEXT_D)[bmvar.v] = bmvar.LAST_D;
			continue;
        }
		(bmvar.LINK)[bmvar.v] = -(bmvar.LINK)[bmvar.v];
		bmvar.i = bmvar.v;
		while (bmvar.i != bmvar.DUMMYVERTEX)
        {
			(bmvar.Y)[bmvar.i] -= bmvar.DELTA;
			bmvar.i = (bmvar.NEXTVTX)[bmvar.i];
		}
		bmvar.f = (bmvar.MATE)[bmvar.v];
		if (bmvar.f != bmvar.DUMMYEDGE)
        {
			bmvar.i = BEND(bmvar.f);
			del = SLACK(bmvar.f);
			while (bmvar.i != bmvar.DUMMYVERTEX)
            {
				(bmvar.Y)[bmvar.i] -= del;
				bmvar.i = (bmvar.NEXTVTX)[bmvar.i];
			}
		}
		(bmvar.NEXT_D)[bmvar.v] = bmvar.LAST_D;
	}
}
/***   BSET_BOUNDS   ***/




/***   BUNPAIR_ALL   **********************************************************
 *
 * undoes all blossoms to get the final matching
 *
 */
void BUNPAIR_ALL ()
{
	int u;

	for (bmvar.v=1; bmvar.v <= bmvar.U; ++(bmvar.v))
    {
		if ((bmvar.BASE)[bmvar.v] != bmvar.v || (bmvar.LASTVTX)[bmvar.v] == bmvar.v)
			continue;
		bmvar.nextu = bmvar.v;
		(bmvar.NEXTVTX)[(bmvar.LASTVTX)[bmvar.nextu]] = bmvar.DUMMYVERTEX;
		while (1)
        {
			u = bmvar.nextu;
			bmvar.nextu = (bmvar.NEXTVTX)[bmvar.nextu];
			BUNLINK (u);
			if ((bmvar.LASTVTX)[u] != u)
            {
				bmvar.f = ((bmvar.LASTEDGE)[2] == OPPEDGE(bmvar.e)) ? (bmvar.LASTEDGE)[1] : (bmvar.LASTEDGE)[2];
				(bmvar.NEXTVTX)[(bmvar.LASTVTX)[BEND(bmvar.f)]] = u;
			}
			bmvar.newbase = BMATE (BMATE(u));
			if (bmvar.newbase != bmvar.DUMMYVERTEX && bmvar.newbase != u)
            {
				(bmvar.LINK)[u] = -bmvar.DUMMYEDGE;
				BREMATCH (bmvar.newbase, (bmvar.MATE)[u]);
			}
			while ((bmvar.LASTVTX)[bmvar.nextu] == bmvar.nextu && bmvar.nextu != bmvar.DUMMYVERTEX)
				bmvar.nextu = (bmvar.NEXTVTX)[bmvar.nextu];
			if ((bmvar.LASTVTX)[bmvar.nextu] == bmvar.nextu && bmvar.nextu == bmvar.DUMMYVERTEX)
				break;
		}
	}
}
/***   BUNPAIR_ALL   ***/




/***   BPOINTER   *************************************************************
 *
 * Assign a pointer link to a vertex.  Edge e joins a vertex in blossom
 * u to a linked vertex.
 *
 */
void BPOINTER (int u, int v, int e)
{
	int i, del;

    #ifdef DEBUG
    	printf("Pointer u,v,e=%d %d %d, %d-%d\n",u,v,e,(bmvar.END)[OPPEDGE(e)],(bmvar.END)[e]);
    	fflush(NULL);
    #endif

	(bmvar.LINK)[u] = -bmvar.DUMMYEDGE;
	(bmvar.NEXTVTX)[(bmvar.LASTVTX)[u]] = bmvar.DUMMYVERTEX;
	(bmvar.NEXTVTX)[bmvar.LASTVTX[v]] = bmvar.DUMMYVERTEX;
	
	if ((bmvar.LASTVTX)[u] != u)
    {
		i = (bmvar.MATE)[(bmvar.NEXTVTX)[u]];
		del = -SLACK(i) / 2;
	}
    else
    {
        del = bmvar.LAST_D;
    }

	i = u;
	while (i != bmvar.DUMMYVERTEX)
    {
		(bmvar.Y)[i] += del;
		(bmvar.NEXT_D)[i] += del;
		i = (bmvar.NEXTVTX)[i];
	}
	if ((bmvar.LINK)[v] < 0)
    {
		(bmvar.LINK)[v] = e;
		(bmvar.NEXTPAIR)[bmvar.DUMMYEDGE] = bmvar.DUMMYEDGE;
		BSCAN (v, bmvar.DELTA);
		return;
	}
    else
    {
		(bmvar.LINK)[v] = e;
		return;
	}
}
/***   BPOINTER   ***/




/***   BSCAN   ****************************************************************
 *
 * Scan each vertex in the blossom whose base is x
 *
 */
void BSCAN (int x, int del)
{
	int u, del_e;

    #ifdef DEBUG
    	printf("Scan del=%d x=%d\n",del,x);
    #endif

	bmvar.newbase = (bmvar.BASE)[x];
	bmvar.stopscan = (bmvar.NEXTVTX)[(bmvar.LASTVTX)[x]];
	while (x != bmvar.stopscan)
    {
		(bmvar.Y)[x] += del;
		(bmvar.NEXT_D)[x] = (bmvar.LAST_D);
		bmvar.pairpoint = bmvar.DUMMYEDGE;
		bmvar.e = (bmvar.A)[x];
		while (bmvar.e != 0)
        {
			bmvar.neighbor = (bmvar.END)[bmvar.e];
			u = (bmvar.BASE)[bmvar.neighbor];
			if ((bmvar.LINK)[u] < 0)
            {
				if ((bmvar.LINK)[BMATE (u)] < 0 || (bmvar.LASTVTX)[u] != u)
                {
					del_e = SLACK (bmvar.e);
					if ((bmvar.NEXT_D)[bmvar.neighbor] > del_e)
                    {
						(bmvar.NEXT_D)[bmvar.neighbor] = del_e;
						(bmvar.NEXTEDGE)[bmvar.neighbor] = bmvar.e;
					}
				}
			}
            else if (u != bmvar.newbase)
            {
				BINSERT_PAIR();
			}
			bmvar.e = (bmvar.A)[bmvar.e];
		}
		x = (bmvar.NEXTVTX)[x];
	}
	(bmvar.NEXTEDGE)[bmvar.newbase] = (bmvar.NEXTPAIR)[bmvar.DUMMYEDGE];
}
/***   BSCAN   ***/




/***   BPAIR   ****************************************************************
 *
 * Process an edge linking two linked vertices
 * Note: global variable v set to the base of one end of the linking edge
 *
 */
void BPAIR (int *outcome)
{
	int u, w, temp;

    #ifdef DEBUG
        printf("Pair v=%d\n",(bmvar.v));
    #endif

	(bmvar.e) = (bmvar.NEXTEDGE)[bmvar.v];
	while (SLACK((bmvar.e)) != 2*(bmvar.DELTA))
		(bmvar.e) = (bmvar.NEXTPAIR)[(bmvar.e)];
	w = BEND ((bmvar.e));
	(bmvar.LINK)[BMATE (w)] = -(bmvar.e);
	u = BMATE ((bmvar.v));
	while ((bmvar.LINK)[u] != -(bmvar.e))
    {
		(bmvar.LINK)[u] = -(bmvar.e);
		if ((bmvar.MATE)[w] != (bmvar.DUMMYEDGE))
        {
			temp = (bmvar.v);
			(bmvar.v) = w;
			w = temp;
		}
		(bmvar.v) = BLINK ((bmvar.v));
		u = BMATE ((bmvar.v));
	}
	if (u == (bmvar.DUMMYVERTEX) && (bmvar.v) != w)
    {
		*outcome = 1;
		return;
	}
	bmvar.newlast = bmvar.v;
	bmvar.newbase = bmvar.v;
	bmvar.oldfirst = (bmvar.NEXTVTX)[bmvar.v];
	BLINK_PATH (bmvar.e);
	BLINK_PATH (OPPEDGE (bmvar.e));
	(bmvar.NEXTVTX)[bmvar.newlast] = bmvar.oldfirst;
	if ((bmvar.LASTVTX)[bmvar.newbase] == bmvar.newbase)
		(bmvar.LASTVTX)[bmvar.newbase] = bmvar.newlast;
	(bmvar.NEXTPAIR)[bmvar.DUMMYEDGE] = bmvar.DUMMYEDGE;
	BMERGE_PAIRS (bmvar.newbase);
	bmvar.i = (bmvar.NEXTVTX)[bmvar.newbase];
	do
    {
		BMERGE_PAIRS (bmvar.i);
		bmvar.i = (bmvar.NEXTVTX)[(bmvar.LASTVTX)[bmvar.i]];
		BSCAN (bmvar.i, 2*bmvar.DELTA - SLACK((bmvar.MATE)[bmvar.i]));
		bmvar.i = (bmvar.NEXTVTX)[(bmvar.LASTVTX)[bmvar.i]];
	} while (bmvar.i != bmvar.oldfirst);
	*outcome = 0;
	return;
}
/***   BPAIR   ***/




/***   BMERGE_PAIRS   *********************************************************
 *
 * merges a subblossom's pair list into a new blossom's pair list
 * v is the base of the previously unlinked subblossom
 * Note: global variable newbase set to the base of the new blossom
 * 	called with NEXTPAIR[DUMMYEDGE] pointing to the first edge
 *		on newbase's pair list
 */
void BMERGE_PAIRS (int v)
{
    #ifdef DEBUG
    	printf("Merge Pairs v=%d\n",v);
    #endif

	(bmvar.NEXT_D)[v] = bmvar.LAST_D;
	bmvar.pairpoint = bmvar.DUMMYEDGE;
	bmvar.f = (bmvar.NEXTEDGE)[v];
	while (bmvar.f != bmvar.DUMMYEDGE)
    {
		bmvar.e = bmvar.f;
		bmvar.neighbor = (bmvar.END)[bmvar.e];
		bmvar.f = (bmvar.NEXTPAIR)[bmvar.f];
		if ((bmvar.BASE)[bmvar.neighbor] != bmvar.newbase)
			BINSERT_PAIR();
	}
}
/***   BMERGE_PAIRS   ***/



/***   BLINK_PATH   ***********************************************************
 *
 * links the unlinked vertices in the path P(END[e],newbase)
 * Note: global variable newbase is set to the base vertex of the new blossom
 *		newlast is set to the last vertex in newbase's current blossom
 */
void BLINK_PATH (int e)
{
	int u;

    #ifdef DEBUG
    	printf("Link Path e=%d-%d\n", (bmvar.END)[OPPEDGE(e)], (bmvar.END)[e]);
    #endif

	bmvar.v = BEND (e);
	while (bmvar.v != bmvar.newbase)
    {
		u = BMATE (bmvar.v);
		(bmvar.LINK)[u] = OPPEDGE (e);
		(bmvar.NEXTVTX)[bmvar.newlast] = bmvar.v;
		(bmvar.NEXTVTX)[(bmvar.LASTVTX)[bmvar.v]] = u;
		(bmvar.newlast) = (bmvar.LASTVTX)[u];
		bmvar.i = bmvar.v;
		(bmvar.BASE)[bmvar.i] = bmvar.newbase;
		bmvar.i = (bmvar.NEXTVTX)[bmvar.i];
		while (bmvar.i != bmvar.DUMMYVERTEX)
        {
			(bmvar.BASE)[bmvar.i] = bmvar.newbase;
			bmvar.i = (bmvar.NEXTVTX)[bmvar.i];
		}
		e = (bmvar.LINK)[bmvar.v];
		bmvar.v = BEND (e);
	}
}
/***   BLINK_PATH   ***/



/***   BINSERT_PAIR   *********************************************************
 *
 * Update a blossom's pair list.
 * Note: called with global variable e set to the edge to be inserted.
 *			neighbor set to the vertex at the end of e
 *			pairpoint set to the next pair on the pair list
 */
void BINSERT_PAIR ()
{
	int del_e;

    #ifdef DEBUG
    	printf("Insert Pair e=%d-%d\n",(bmvar.END)[OPPEDGE(bmvar.e)],(bmvar.END)[bmvar.e]);
    #endif

	del_e = SLACK(bmvar.e)/2;
	bmvar.nextpoint = (bmvar.NEXTPAIR)[bmvar.pairpoint];

	while ((bmvar.END)[bmvar.nextpoint] < bmvar.neighbor)
    {
		bmvar.pairpoint = bmvar.nextpoint;
		bmvar.nextpoint = (bmvar.NEXTPAIR)[bmvar.nextpoint];
	}
	if ((bmvar.END)[bmvar.nextpoint] == bmvar.neighbor)
    {
		if (del_e >= SLACK (bmvar.nextpoint)/2)
			return;
		bmvar.nextpoint = (bmvar.NEXTPAIR)[bmvar.nextpoint];
	}
	(bmvar.NEXTPAIR)[bmvar.pairpoint] = bmvar.e;
	bmvar.pairpoint = bmvar.e;
	(bmvar.NEXTPAIR)[bmvar.e] = bmvar.nextpoint;
	if ((bmvar.NEXT_D)[bmvar.newbase] > del_e)
		(bmvar.NEXT_D)[bmvar.newbase] = del_e;
}
/***   BINSERT_PAIR   ***/
