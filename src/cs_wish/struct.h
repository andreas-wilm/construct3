/******************************************************************************
* 
* struct.h - some useful macros and functions for structure handling
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
 *  CVS $Id: struct.h,v 1.5 2004-05-25 13:29:18 wilm Exp $    
 */


#ifndef STRUCT_H_INCL
#define STRUCT_H_INCL


/* FIXME: should be identical in size to BRANCH in CS_DrawStruct.c */
#define MAX_NUMBIF 100



/* minimal hairpin size */
#define MIN_HP_SIZE 4

/* FIXME: maximum number of bifurcations
          should be identical in size to BRANCH in CS_DrawStruct.c 
 */
#define MAX_NUMBIF 100


#ifndef MAX
    #define MAX(x,y) (((x)>(y)) ? (x) : (y))
#endif

#ifndef INT_MAX
    #define INT_MAX			0x7fffffff
#endif


typedef struct  {
    int   pair;
    float prob;
} pair_prob_t;


extern void
RemoveSingleBps(pair_prob_t *structure, int aln_len);

#endif
