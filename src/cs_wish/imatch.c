/******************************************************************************
* 
* imatch.c - procedures for maximum weighted matching
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
 *  CVS $Id: imatch.c,v 1.11 2004-05-25 13:29:15 wilm Exp $    
 */



/***   Routines for Maximum Weigthed Matching    ******************************
 *
 * N-cubed weighted matching
 * Implementation of H. Gabow's Ph.D. thesis, Stanford Univ. 1973
 * Written by Edward Rothberg  7/85
 * For complete details, please refer to the original paper
 * Modified by Jack Tabaska 12/96 to allow k-matching to be performed
 *
 * Modified by J. Riks for use in construct. See Tabaska, Bioinformatics
 * Vol. 14 no. 8 1998 (691-699) for further details
 *
 *****************************************************************************/


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"

#include "imatch.h"



/***   private
 */
#if 0
    #define DEBUG
#endif


#define	MAXWT  1000000


/* the number of the blossom entered by edge e */
#define BEND(e) ((imvar.BASE)[imvar.END[e]])

/* the blossom matched with v's blossom */
#define BMATE(v) ((imvar.BASE)[(imvar.END)[(imvar.MATE)[v]]])

/* the blossom entered by the edge that links v's blossom */
#define BLINK(v) ((imvar.BASE)[(imvar.END)[(imvar.LINK)[v]]])

/* the edge e with it's direction reversed */
#define OPPEDGE(e) (((e - (imvar.U)) % 2 == 0) ? (e - 1) : (e + 1))

/* the slack of edge e */
#define SLACK(e) ((imvar.Y)[(imvar.END)[e]] + (imvar.Y)[(imvar.END)[OPPEDGE(e)]] - (imvar.WEIGHT)[e])


static void
Initialize(int maximize);

static void
PAIR (int *outcome);

static void
MERGE_PAIRS (int v);

static void
LINK_PATH (int e);

static void
INSERT_PAIR ();

static void
POINTER (int u, int v, int e);

static void
SCAN (int x, int del);

static void
SET_BOUNDS ();

static void
UNPAIR_ALL ();

static void UNPAIR
(int oldbase, int oldmate);

static void
REMATCH (int firstmate, int e);

static void
UNLINK (int oldbase);

static void 
Initialize(int maximize);
 
/*
 ***/






/***   Weighted_Match   *******************************************************
 *
 */
int
*Weighted_Match (int type, int maximize, int k)
{

	int g, j, w, outcome, loop=0;
	
	/* set up internal data structure */
	Initialize(maximize);

	for(;;) {
		loop++;
	
		if (k && (loop > k))
			break;
	
		imvar.DELTA = 0;
		for (imvar.v=1; imvar.v<=imvar.U; ++imvar.v)
			if ((imvar.MATE)[(imvar.v)] == (imvar.DUMMYEDGE))
				POINTER ((imvar.DUMMYVERTEX), (imvar.v), (imvar.DUMMYEDGE));
		for(;;) {
			imvar.i = 1;
			for (j=2; j<=(imvar.U); ++j)
				if ((imvar.NEXT_D)[imvar.i] > (imvar.NEXT_D)[j])
					imvar.i = j;
			imvar.DELTA = (imvar.NEXT_D)[imvar.i];
			if (imvar.DELTA == imvar.LAST_D)
				goto done;
			(imvar.v) = (imvar.BASE)[imvar.i];
			if ((imvar.LINK)[imvar.v] >= 0) {
				PAIR (&outcome);
				if (outcome == 1)
					break;
			} else {
				w = BMATE (imvar.v);
				if ((imvar.LINK)[w] < 0) {
					POINTER (imvar.v, w, OPPEDGE((imvar.NEXTEDGE[imvar.i])));
				} else UNPAIR (imvar.v, w);
			}
		}
	
		imvar.LAST_D -=imvar.DELTA;
		SET_BOUNDS();
		g = OPPEDGE(imvar.e);
		REMATCH (BEND(imvar.e), g);
		REMATCH (BEND(g), imvar.e);
	}
	done:
    
    #ifdef DEBUG
    	printf("\nResult of i-matching: % 3d Edges\n-------------------------------\n", loop - 1);
    #endif
    
	SET_BOUNDS();
	UNPAIR_ALL();
	for (imvar.i=1; imvar.i<=imvar.U; ++(imvar.i))
    {
		(imvar.MATE)[imvar.i] = (imvar.END)[(imvar.MATE)[imvar.i]];
		if ((imvar.MATE)[imvar.i]==(imvar.DUMMYVERTEX))
			(imvar.MATE)[imvar.i]=0;
	}

    /*	FreeUp(); */
	return(imvar.MATE);
}
/***   Weighted_Match   ***/




/***   PAIR   *****************************************************************
 *
 * Process an edge linking two linked vertices
 * Note: global variable imvar.v set to the base of
 * one end of the linking edge
 */
void
PAIR (int *outcome)
{

	int u, w, temp;

    #ifdef DEBUG
        printf("Pair imvar.v=%d\n",imvar.v);
    #endif

	imvar.e = imvar.NEXTEDGE[imvar.v];
	while (SLACK(imvar.e) != 2*imvar.DELTA)
		imvar.e = (imvar.NEXTPAIR)[imvar.e];
	
    w = BEND (imvar.e);
	(imvar.LINK)[BMATE (w)] = -imvar.e;
	u = BMATE (imvar.v);
	
    while ((imvar.LINK)[u] != -imvar.e)
    {
		(imvar.LINK)[u] = -(imvar.e);
		if ((imvar.MATE)[w] != imvar.DUMMYEDGE)
        {
			temp = imvar.v;
			imvar.v = w;
			w = temp;
		}
		imvar.v = BLINK (imvar.v);
		u = BMATE (imvar.v);
	}
    
	if (u == imvar.DUMMYVERTEX && imvar.v != w)
    {
		*outcome = 1;
		return;
	}
    
	imvar.newlast  = imvar.v;
	imvar.newbase  = imvar.v;
	imvar.oldfirst = (imvar.NEXTVTX)[imvar.v];
	LINK_PATH (imvar.e);
	LINK_PATH (OPPEDGE (imvar.e));
	(imvar.NEXTVTX)[imvar.newlast] = imvar.oldfirst;
	
    if ((imvar.LASTVTX)[imvar.newbase] == imvar.newbase)
		(imvar.LASTVTX)[imvar.newbase] = imvar.newlast;
	
    (imvar.NEXTPAIR)[imvar.DUMMYEDGE] = imvar.DUMMYEDGE;
	MERGE_PAIRS (imvar.newbase);
	imvar.i = (imvar.NEXTVTX)[imvar.newbase];
	do {
		MERGE_PAIRS (imvar.i);
		imvar.i = (imvar.NEXTVTX)[(imvar.LASTVTX)[imvar.i]];
		SCAN (imvar.i, 2*(imvar.DELTA) - SLACK((imvar.MATE)[imvar.i]));
		imvar.i = (imvar.NEXTVTX)[(imvar.LASTVTX)[(imvar.i)]];
	} while (imvar.i != imvar.oldfirst);
	
    *outcome = 0;
	return;
}
/***   PAIR   ***/



/***   MERGE_PAIRS   **********************************************************
 *
 * merges a subblossom's pair list into a new blossom's pair list
 * v is the base of the previously unlinked subblossom
 * Note: global variable (imvar.newbase) set to the base of the new blossom
 *       called with (imvar.NEXTPAIR)[(imvar.DUMMYEDGE)] pointing to the first edge
 *       on imvar.newbase's pair list
 *
 */
void
MERGE_PAIRS (int v)
{

    #ifdef DEBUG
        printf("Merge Pairs v=%d\n",v);
    #endif

	(imvar.NEXT_D)[v] = imvar.LAST_D;
	imvar.pairpoint = imvar.DUMMYEDGE;
	imvar.f = (imvar.NEXTEDGE)[v];
	
    while (imvar.f != imvar.DUMMYEDGE)
    {
		imvar.e = imvar.f;
		imvar.neighbor = (imvar.END)[imvar.e];
		imvar.f = (imvar.NEXTPAIR)[imvar.f];
		if ((imvar.BASE)[imvar.neighbor] != imvar.newbase)
			INSERT_PAIR();
	}
}
/***   MERGE_PAIRS   ***/



/***   LINK_PATH   ************************************************************
 *
 * links the unlinked vertices in the path P((imvar.END)[e],(imvar.newbase))
 * Note: global variable imvar.newbase is set to the base vertex of the new blossom
 *		imvar.newlast is set to the last vertex in imvar.newbase's current blossom
 *
 */
void
LINK_PATH (int e)
{

	int u;

    #ifdef DEBUG
       printf("Link Path e=%d-%d\n", (imvar.END)[OPPEDGE(e)], (imvar.END)[e]);
    #endif

	imvar.v = BEND (e);
	while (imvar.v != imvar.newbase)
    {
		u = BMATE (imvar.v);
		(imvar.LINK)[u] = OPPEDGE (e);
		(imvar.NEXTVTX)[imvar.newlast] = imvar.v;
		(imvar.NEXTVTX)[(imvar.LASTVTX)[imvar.v]] = u;
		imvar.newlast = (imvar.LASTVTX)[u];
		imvar.i = imvar.v;
		(imvar.BASE)[imvar.i] = imvar.newbase;
		imvar.i = (imvar.NEXTVTX)[imvar.i];
		
        while (imvar.i != imvar.DUMMYVERTEX)
        {
			(imvar.BASE)[imvar.i] = imvar.newbase;
			imvar.i = (imvar.NEXTVTX)[imvar.i];
		}
		
        e = (imvar.LINK)[imvar.v];
		imvar.v = BEND (e);
	}
}
/***   LINK_PATH   ***/



/***   INSERT_PAIR   **********************************************************
 *
 *
 * Update a blossom's pair list.
 * Note: called with global variable imvar.e set to the edge to be inserted.
 *			imvar.neighbor set to the vertex at the end of imvar.e
 *			imvar.pairpoint set to the next pair on the pair list
 *
 */
void INSERT_PAIR ()
{

	int del_e;

    #ifdef DEBUG
    	printf("Insert Pair imvar.e=%d-%d\n",(imvar.END)[OPPEDGE(imvar.e)],(imvar.END)[imvar.e]);
    #endif

	del_e = SLACK(imvar.e)/2;
	imvar.nextpoint = (imvar.NEXTPAIR)[imvar.pairpoint];

	while ((imvar.END)[(imvar.nextpoint)] < (imvar.neighbor)) {
		imvar.pairpoint = imvar.nextpoint;
		imvar.nextpoint = (imvar.NEXTPAIR)[imvar.nextpoint];
	}
	if ((imvar.END)[imvar.nextpoint] == imvar.neighbor) {
		if (del_e >= SLACK (imvar.nextpoint)/2)
			return;
		imvar.nextpoint = (imvar.NEXTPAIR)[imvar.nextpoint];
	}
	(imvar.NEXTPAIR)[(imvar.pairpoint)] = imvar.e;
	imvar.pairpoint = imvar.e;
	(imvar.NEXTPAIR)[imvar.e] = imvar.nextpoint;
	if ((imvar.NEXT_D)[imvar.newbase] > del_e)
		(imvar.NEXT_D)[imvar.newbase] = del_e;
}
/***   INSERT_PAIR   ***/



/***   POINTER   **************************************************************
 *
 *
 * Assign a pointer link to a vertex.  Edge e joins a vertex in blossom
 * u to a linked vertex.
 */
void
POINTER (int u, int v, int e)
{

	int i, del;

    #ifdef DEBUG
    	printf("Pointer u,v,e=%d %d %d-%d\n",u,v,(imvar.END)[OPPEDGE(e)],(imvar.END)[e]);
    #endif

	(imvar.LINK)[u] = -imvar.DUMMYEDGE;
	(imvar.NEXTVTX)[(imvar.LASTVTX)[u]] = imvar.DUMMYVERTEX;
	(imvar.NEXTVTX)[(imvar.LASTVTX)[v]] = imvar.DUMMYVERTEX;
	
	if ((imvar.LASTVTX)[u] != u) {
		i = (imvar.MATE)[(imvar.NEXTVTX)[u]];
		del = -SLACK(i) / 2;
	}
	else del = imvar.LAST_D;

	i = u;
	while (i != imvar.DUMMYVERTEX)
    {
		(imvar.Y)[i] += del;
		(imvar.NEXT_D)[i] += del;
		i = (imvar.NEXTVTX)[i];
	}
	if ((imvar.LINK)[v] < 0)
    {
		(imvar.LINK)[v] = e;
		(imvar.NEXTPAIR)[imvar.DUMMYEDGE] = imvar.DUMMYEDGE;
		SCAN (v, imvar.DELTA);
		return;
	}
    else
    {
		(imvar.LINK)[v] = e;
		return;
	}
}
/***   POINTER   ***/



/***   SCAN   *****************************************************************
 *
 * Scan each vertex in the blossom whose base is x
 *
 */
void
SCAN (int x, int del)
{

	int u, del_e;

    #ifdef DEBUG
    	printf("Scan del=%d x=%d\n",del,x);
    #endif
    
	imvar.newbase = (imvar.BASE)[x];
	imvar.stopscan = (imvar.NEXTVTX)[(imvar.LASTVTX)[x]];
	while (x != imvar.stopscan)
    {
		(imvar.Y)[x] += del;
		(imvar.NEXT_D)[x] = imvar.LAST_D;
		imvar.pairpoint = imvar.DUMMYEDGE;
		imvar.e = (imvar.A)[x];
		while (imvar.e != 0)
        {
			imvar.neighbor = (imvar.END)[imvar.e];
			u = (imvar.BASE)[imvar.neighbor];
			if ((imvar.LINK)[u] < 0)
            {
				if ((imvar.LINK)[BMATE (u)] < 0 || (imvar.LASTVTX)[u] != u)
                {
					del_e = SLACK (imvar.e);
					if ((imvar.NEXT_D)[imvar.neighbor] > del_e)
                    {
						(imvar.NEXT_D)[imvar.neighbor] = del_e;
						(imvar.NEXTEDGE)[imvar.neighbor] = imvar.e;
					}
				}
			}
            else if (u != imvar.newbase)
            {
				INSERT_PAIR();
			}
			imvar.e = (imvar.A)[imvar.e];
		}
		x = (imvar.NEXTVTX)[x];
	}
	(imvar.NEXTEDGE)[imvar.newbase] = (imvar.NEXTPAIR)[imvar.DUMMYEDGE];
}
/***   SCAN   ***/



/***   SET_BOUNDS   ***********************************************************
 *
 *
 * updates numerical bounds for linking paths.
 * called with imvar.LAST_D set to the bound on imvar.DELTA for the next search
 */
void
SET_BOUNDS ()
{
	int del;

	for ((imvar.v)=1; (imvar.v) <= (imvar.U); ++(imvar.v))
    {
		if ((imvar.LINK)[(imvar.v)] < 0 || (imvar.BASE)[(imvar.v)] != (imvar.v))
        {
			(imvar.NEXT_D)[(imvar.v)] = (imvar.LAST_D);
			continue;
		}
		(imvar.LINK)[(imvar.v)] = -(imvar.LINK)[(imvar.v)];
		(imvar.i) = (imvar.v);
		while ((imvar.i) != (imvar.DUMMYVERTEX))
        {
			(imvar.Y)[(imvar.i)] -= (imvar.DELTA);
			(imvar.i) = (imvar.NEXTVTX)[(imvar.i)];
		}
		(imvar.f) = (imvar.MATE)[(imvar.v)];
		if ((imvar.f) != (imvar.DUMMYEDGE))
        {
			(imvar.i) = BEND((imvar.f));
			del = SLACK((imvar.f));
			while ((imvar.i) != (imvar.DUMMYVERTEX))
            {
				(imvar.Y)[(imvar.i)] -= del;
				(imvar.i) = (imvar.NEXTVTX)[(imvar.i)];
			}
		}
	(imvar.NEXT_D)[(imvar.v)] = (imvar.LAST_D);
	}
}
/***   SET_BOUNDS   ***/



/***   UNPAIR_ALL   ***********************************************************
 *
 * undoes all blossoms to get the final matching
 *
 */
void
UNPAIR_ALL ()
{

	int u;

	for ((imvar.v)=1; (imvar.v) <= (imvar.U); ++(imvar.v))
    {
		if ((imvar.BASE)[(imvar.v)] != (imvar.v) || (imvar.LASTVTX)[(imvar.v)] == (imvar.v))
			continue;
		(imvar.nextu) = (imvar.v);
		(imvar.NEXTVTX)[(imvar.LASTVTX)[(imvar.nextu)]] = (imvar.DUMMYVERTEX);
		while (1)
        {
			u = (imvar.nextu);
			(imvar.nextu) = (imvar.NEXTVTX)[(imvar.nextu)];
			UNLINK (u);
			if ((imvar.LASTVTX)[u] != u)
            {
				(imvar.f) = ((imvar.LASTEDGE)[2] == OPPEDGE((imvar.e))) ? (imvar.LASTEDGE)[1] : (imvar.LASTEDGE)[2];
				(imvar.NEXTVTX)[(imvar.LASTVTX)[BEND((imvar.f))]] = u;
			}
			(imvar.newbase) = BMATE (BMATE(u));
			if ((imvar.newbase) != (imvar.DUMMYVERTEX) && (imvar.newbase) != u)
            {
				(imvar.LINK)[u] = -(imvar.DUMMYEDGE);
				REMATCH ((imvar.newbase), (imvar.MATE)[u]);
			}
			while ((imvar.LASTVTX)[(imvar.nextu)] == (imvar.nextu) && (imvar.nextu) != (imvar.DUMMYVERTEX))
				(imvar.nextu) = (imvar.NEXTVTX)[(imvar.nextu)];
			if ((imvar.LASTVTX)[(imvar.nextu)] == (imvar.nextu) && (imvar.nextu) == (imvar.DUMMYVERTEX))
				break;
		}
	}
}
/***   UNPAIR_ALL   ***/




/***   UNPAIR   ***************************************************************
 *
 * Expands a blossom.  Fixes up imvar.LINK and imvar.MATE.
 *
 */
void
UNPAIR (int oldbase, int oldmate)
{

   int e, newbase, u;

    #ifdef DEBUG
        printf("Unpair oldbase, oldmate=%d %d\n",oldbase, oldmate);
    #endif

	UNLINK (oldbase);
	newbase = BMATE (oldmate);
	if (newbase != oldbase)
    {
		(imvar.LINK)[oldbase] = -imvar.DUMMYEDGE;
		REMATCH (newbase, imvar.MATE[oldbase]);
		if (imvar.f == imvar.LASTEDGE[1])
			(imvar.LINK)[imvar.secondmate] = -imvar.LASTEDGE[2];
		else (imvar.LINK)[imvar.secondmate] = -imvar.LASTEDGE[1];
	}
	e = (imvar.LINK)[oldmate];
	u = BEND (OPPEDGE (e));
	if (u == newbase)
    {
		POINTER (newbase, oldmate, e);
		return;
	}
	(imvar.LINK)[BMATE (u)] = -e;
	do {
		e = -(imvar.LINK)[u];
		imvar.v = BMATE (u);
		POINTER (u, imvar.v, -(imvar.LINK)[imvar.v]);
		u = BEND (e);
	} while (u != newbase);
	e = OPPEDGE (e);
	POINTER (newbase, oldmate, e);
}
/***   UNPAIR   ***/



/***   REMATCH   **************************************************************
 *
 * changes the matching along an alternating path
 * firstmate is the first base vertex on the path
 * edge e is the new matched edge for firstmate
 *
 */
void
REMATCH (int firstmate, int e)
{

    #ifdef DEBUG
         printf("Rematch firstmate=%d e=%d-%d\n",firstmate, (imvar.END)[OPPEDGE(e)], (imvar.END)[e]);
    #endif

	(imvar.MATE)[firstmate] = e;
	(imvar.nexte) = -(imvar.LINK)[firstmate];
	if ((imvar.nexte) != (imvar.DUMMYEDGE))
    {
		(imvar.REMATCH_FLAGS)[firstmate] = (imvar.REMATCH_FLAGS)[(imvar.END)[e]] = 1;
        /* printf("Vertices rematched:\n");
		   printf("  %d\n  %d\n", firstmate, (imvar.END)[e]); */
	}
	while ((imvar.nexte) != (imvar.DUMMYEDGE))
    {
		e = (imvar.nexte);
		(imvar.f) = OPPEDGE (e);
		firstmate = BEND (e);
		(imvar.secondmate) = BEND (imvar.f);
		(imvar.REMATCH_FLAGS)[firstmate] = (imvar.REMATCH_FLAGS)[(imvar.secondmate)] = 1;
        /*	printf("  %d\n  %d\n", firstmate, (imvar.secondmate)); */
		(imvar.nexte) = -(imvar.LINK)[firstmate];
		(imvar.LINK)[firstmate] = -(imvar.MATE)[(imvar.secondmate)];
		(imvar.LINK)[(imvar.secondmate)] = -(imvar.MATE)[firstmate];
		(imvar.MATE)[firstmate] = (imvar.f);
		(imvar.MATE)[(imvar.secondmate)] = e;
	}
}
/***   REMATCH   ***/



/***   UNLINK   ***************************************************************
 *
 * unlinks subblossoms in a blossom.  oldbase is the base of the blossom to
 * be unlinked.
 *
 */
void
UNLINK (int oldbase)
{
	int k, j=1;

    #ifdef DEBUG
        printf("Unlink oldbase=%d\n",oldbase);
    #endif

	(imvar.i) = (imvar.NEXTVTX)[oldbase];
	(imvar.newbase) = (imvar.NEXTVTX)[oldbase];
	(imvar.nextbase) = (imvar.NEXTVTX)[(imvar.LASTVTX)[(imvar.newbase)]];
	(imvar.e) = (imvar.LINK)[(imvar.nextbase)];
UL2:
	do {
		(imvar.nextedge) = OPPEDGE ((imvar.LINK)[(imvar.newbase)]);
		for (k=1; k <= 2; ++k) {
			(imvar.LINK)[(imvar.newbase)] = -(imvar.LINK)[(imvar.newbase)];
			(imvar.BASE)[(imvar.i)] = (imvar.newbase);
			(imvar.i) = (imvar.NEXTVTX)[(imvar.i)];
			while ((imvar.i) != (imvar.nextbase))
            {
				(imvar.BASE)[(imvar.i)] = (imvar.newbase);
				(imvar.i) = (imvar.NEXTVTX)[(imvar.i)];
			}
			(imvar.newbase) = (imvar.nextbase);
			(imvar.nextbase) = (imvar.NEXTVTX)[(imvar.LASTVTX)[(imvar.newbase)]];
		}
	} while ((imvar.LINK)[(imvar.nextbase)] == (imvar.nextedge));
	if (j==1)
    {
		(imvar.LASTEDGE[1]) = (imvar.nextedge);
		j++;
		(imvar.nextedge) = OPPEDGE ((imvar.e));
		if ((imvar.LINK)[(imvar.nextbase)] == (imvar.nextedge))
			goto UL2;
	}
	(imvar.LASTEDGE[2]) = (imvar.nextedge);

	if ((imvar.BASE)[(imvar.LASTVTX)[oldbase]] == oldbase)
    {
		(imvar.NEXTVTX)[oldbase] = (imvar.newbase);
	} else {
		(imvar.NEXTVTX)[oldbase] = (imvar.DUMMYVERTEX);
		(imvar.LASTVTX)[oldbase] = oldbase;
	}
}




/***   Initialize   ***********************************************************
 *
 */
void
Initialize(int maximize)
{
	int i, max_wt= -MAXWT, min_wt=MAXWT;

	(imvar.DUMMYVERTEX) = imvar.U+1;
	(imvar.DUMMYEDGE)   = imvar.U+2*imvar.V+1;
	(imvar.END)[(imvar.DUMMYEDGE)] = imvar.DUMMYVERTEX;

	for (i=imvar.U+2; i<=imvar.U+2*imvar.V; i+=2)
    {
		if ((imvar.WEIGHT)[i] > max_wt)
			max_wt = (imvar.WEIGHT)[i];
		if ((imvar.WEIGHT)[i] < min_wt)
			min_wt = (imvar.WEIGHT)[i];
	}
	if (!maximize)
    {
		if (imvar.U % 2 != 0)
        {
			printf("Must have an even number of vertices to do a\n");
			printf("minimum complete matching.\n");
			exit(0);
		}
		max_wt += 2; /* Don't want all zero weight */
		for (i=imvar.U+1; i<=imvar.U+2*imvar.V; i++)
			(imvar.WEIGHT)[i] = max_wt-(imvar.WEIGHT)[i];
		max_wt = max_wt-min_wt;
	}
	imvar.LAST_D = max_wt/2;

	imvar.MATE  	= (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.LINK  	= (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.BASE  	= (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.NEXTVTX   = (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.LASTVTX   = (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.Y  		= (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.NEXT_D	= (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.NEXTEDGE  = (int *) Scalloc(imvar.U+2, sizeof(int));
	imvar.NEXTPAIR = (int *) Scalloc(imvar.U+2*imvar.V+2, sizeof(int));

	for (i = 1; i <= imvar.U+1; ++i)
    {
		(imvar.MATE)[i] = imvar.DUMMYEDGE;
		(imvar.NEXTEDGE[i]) = imvar.DUMMYEDGE;
		(imvar.NEXTVTX)[i] = 0;
		(imvar.LINK)[i] = -imvar.DUMMYEDGE;
		(imvar.BASE)[i] = i;
		(imvar.LASTVTX)[i] = i;
		(imvar.Y)[i] = imvar.LAST_D;
		(imvar.NEXT_D)[i] = imvar.LAST_D;
	}
}
/***   Initialize   ***/



/***   FreeUpImvar   **********************************************************
 *
 */
void
FreeUpImvar()
{
   free(imvar.MATE); 
   free(imvar.LINK);    
   free(imvar.BASE);   
   free(imvar.NEXTVTX); 
   free(imvar.LASTVTX);
   free(imvar.Y);       
   free(imvar.NEXT_D);  
   free(imvar.NEXTEDGE);
   free(imvar.NEXTPAIR);
}
/***   FreeUpImvar   ***/
