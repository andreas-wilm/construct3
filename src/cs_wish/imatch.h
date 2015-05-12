/******************************************************************************
* 
* imatch.h - procedures for maximum weighted matching
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
 *  CVS $Id: imatch.h,v 1.7 2004-05-25 13:29:15 wilm Exp $    
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
 */


#ifndef IMATCH_H_INCL
#define IMATCH_H_INCL

extern int *
Weighted_Match (int type, int maximize, int k);

extern void
FreeUpImvar(void);


/* Gglobal structure Imatch_Var imvar
 */
struct Imatch_Var {
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
	char *REMATCH_FLAGS;
};



struct Imatch_Var imvar;


#endif
