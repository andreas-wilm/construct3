/******************************************************************************
* 
* tinoco.h - Create tinoco dotplots
*
* Copyright (C) 2002 Gerhard Steger <steger@biophys.uni-duesseldorf.de>
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
 *  CVS $Id: tinoco.h,v 1.5 2007-05-10 16:24:25 steger Exp $    
 */


#ifndef TINOCO_H_INCL
#define TINOCO_H_INCL

#define MAXLENHELIX  50
#define TRUE          1
#define FALSE         0
#define LEN_FNAME    21 
#define LEN_DP_FNAME LEN_FNAME+12	/* max 20 chars + "_ti_pr.dat.Z" */

void usage(void);
void nrerror(char *message);

char scale1[]     = ".........1.........2.........3.........4";
char scale2[]     = ".........5.........6.........7.........8";
char argstr[200]  = "";

float *matrix;

#endif
