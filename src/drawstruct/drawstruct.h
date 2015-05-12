/******************************************************************************
* 
* drawstruct.h - main procedures for ConStruct's wish (cs_wish)
*
* Copyright (C) 2001-2004 Institute of Physical Biology,
#                    Heinrich-Heine Universität Düsseldorf, Germany
#                    <construct@biophys.uni-duesseldorf.de>
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
 *  CVS $Id: drawstruct.h,v 1.6 2004-05-25 13:29:22 wilm Exp $    
 */





/*****************************************************************************\
 **********************     par.h     ****************************************
\*****************************************************************************/

#define MAXNUCS   100
#define MAXANGLES	100

typedef struct {
	int	numbers, backbone, quiet, text, bonds;
	int	nuc,   numnuc,     nuc0[MAXNUCS],   nuc9[MAXNUCS];
	int	range, rangestep,  range0,     range1;
	int	pivot, numangle, angle0[MAXANGLES], angle9[MAXANGLES];
	double angle[MAXANGLES];
	float  scf;
	int	init;
	int	single;
	int	straight;
} PAR;




/*****************************************************************************\
 ***********************     cs.h     ****************************************
\*****************************************************************************/


#define MAX(x,y) (((x)>(y)) ? (x) : (y))
 
struct pair_prob {
          int   pair;
          float prob;
       };

struct xyplot {
          int	pair;
          char	nuc;
          int	tag;
          float nuc_x;
          float nuc_y;
          float back_x1;
          float back_y1;
          float back_x2;
          float back_y2;
          float bond_x1;
          float bond_y1;
          float bond_x2;
          float bond_y2;
          float line_x1;
          float line_y1;
          float line_x2;
          float line_y2;
          float lpos_x;
          float lpos_y;
          int	label;
       };
 
extern struct pair_prob *pair;
extern struct xyplot *my_xyplot;




/*****************************************************************************\
\*****************************************************************************/


#define MAXLENGTH  2000
