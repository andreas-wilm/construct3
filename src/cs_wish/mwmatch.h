/******************************************************************************
* 
* mwmatch.h - Routines for weighted matching
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
 *  CVS $Id: mwmatch.h,v 1.10 2004-05-25 13:29:16 wilm Exp $    
 */



/******************************************************************************
 ******************************************************************************
 *                      bmatch and imatch are buggy
 ******************************************************************************
 *   
 *   -> the routines access an invalid memory range: g.edges[1][0]
 *   (given very little data stored in the graph)
 *
 *   there are three ways out:
 *   1. Feed the matching routines with more data
 *      (low colormapping, high mic factor)
 *   2. Init all g.edges with 1 rather then 0
 *                      (PREVENT_MATCH_COREDUMP_BY_EDGE_INIT)
 *   3. Simply forbid access to g.edges[1][0]
 *                      (PREVENT_MATCH_COREDUMP_BY_MEM_ACCESS_TEST)
 *   The latter (in combination with the first) is prefered since that
 *   doesn't lead to the prediction of false basepairs
 *
 *
 ******************************************************************************
 ******************************************************************************/
 
#if 1
    #define PREVENT_MATCH_COREDUMP_BY_MEM_ACCESS_TEST
#else
    #define PREVENT_MATCH_COREDUMP_BY_EDGE_INIT
#endif



#ifndef MWMATCH_H_INCL
#define MWMATCH_H_INCL



#define Edge(g, x, y) ((g).edges)[MAX((x), (y))][MIN((x), (y))]


typedef struct  {
  int **edges;
  int n_vert;
  int n_edge;
} TempGraph_t;


extern void
MakeTempGraph(TempGraph_t *g, int nvert);

extern void
FreeTempGraph(TempGraph_t *g);

extern int
GetCsGraph(TempGraph_t *g, float prob_fac, float mic_fac,
                           float prob_thr, float mic_thr);


#endif
