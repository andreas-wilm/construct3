/******************************************************************************
* 
* main.h - main procedures for ConStruct's wish (cs_wish)
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
 *  CVS $Id: main.h.in,v 1.40 2004-05-25 13:29:15 wilm Exp $
 */

#ifndef MAIN_H_INCL
#define MAIN_H_INCL



#if 0
    #define DEBUG
#endif


#define VERSION_STR "construct 3.2.5"

#define RC_FILE           "~/.cs_wish_v3.rc"




/* if bpprob is below this value
   it will be set to NULL */
#define BP_PROB_CUTOFF 0.01
/* if a consensus basepair prob is below this
   value it won't be drawn */
#define CS_BP_PROB_DRAW_CUTOFF 0.01



/* MAXs and floating point precision */
#include <float.h>
#include <math.h>
#include <limits.h>

#ifndef FLT_DIG
    #define FLT_DIG 6
#endif
#define FLOAT2INT  (pow(10, FLT_DIG))

#ifndef INT_MAX
    #define INT_MAX			0x7fffffff
#endif



/* TD-Consensus Factors
*/
#define POW_A  3
#define POW_B  3





/*********************************   common global variables
 *
 */
extern  proj_t      *proj;
extern  aln_t       *aln;
extern  cs_opts_t    cs_opts;
extern  tk_color_t   tk_color;



                       
                       
#endif
