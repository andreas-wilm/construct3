/******************************************************************************
* 
* bmatch.h - routines for maximum weigthed matching
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
 *  CVS $Id: bmatch.h,v 1.6 2004-05-25 13:29:12 wilm Exp $    
 */


#ifndef BMATCH_H_INCL
#define BMATCH_H_INCL

extern int *
Weighted_Bmatch();

extern void
FreeUpBmvar();



#define	MAXWT  1000000

/* the number of the blossom entered by edge e */
#define BEND(e) ((bmvar.BASE)[bmvar.END[e]])

/* the blossom matched with v's blossom */
#define BMATE(v) ((bmvar.BASE)[(bmvar.END)[(bmvar.MATE)[v]]])

/* the blossom entered by the edge that links v's blossom */
#define BLINK(v) ((bmvar.BASE)[(bmvar.END)[(bmvar.LINK)[v]]])

/* the edge e with it's direction reversed */
#define OPPEDGE(e) (((e - (bmvar.U)) % 2 == 0) ? (e - 1) : (e + 1))

/* the slack of edge e */
#define SLACK(e) ((bmvar.Y)[(bmvar.END)[e]] + (bmvar.Y)[(bmvar.END)[OPPEDGE(e)]] - (bmvar.WEIGHT)[e])



/* global structure Bmatch_Var bmvar */
struct Bmatch_Var {
	int *A;
	int *END;
	int *WEIGHT;
	int *NEXTPAIR;
	int *MATE;
	int *LINK;
	int *BASE;
	int *NEXTVTX;
	int *LASTVTX;
	int *Y;
	int *NEXT_D;
	int *NEXTEDGE;
	int LAST_D;
	int DELTA;
	int LASTEDGE[3];
	int DUMMYVERTEX;
	int DUMMYEDGE;
	int U;
	int V;
	int newbase;
	int oldbase;
	int nextbase;
	int stopscan;
	int pairpoint;
	int neighbor;
	int nextpoint;
	int newlast;
	int newmate;
	int oldmate;
	int oldfirst;
	int firstmate;
	int secondmate;
	int f;
	int nextedge;
	int nexte;
	int nextu;
	int v;
	int i;
	int e;
	int **MatchStore;
	int B;
	int W;
	int OriginalVert;
	int FirstOuter;
	int *FromVtx;
	int *ToVtx;
	char IFlag;
	char *MatchFlags;
};


struct Bmatch_Var bmvar;

#endif

