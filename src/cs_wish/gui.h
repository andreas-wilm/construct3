/******************************************************************************
* 
* gui.h - enhanced tcl routines for the cs_dp gui
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
 *  CVS $Id: gui.h,v 1.16 2007-09-30 12:05:08 wilm Exp $    
 */


#ifndef GUI_H_INCL
#define GUI_H_INCL


/* coords for a tk rectangle */
typedef struct {
    float x1;
    float y1;
    float x2;
    float y2;
} rect_coord_t;


extern int
SetupGui(Tcl_Interp *interp);

extern int
CreateDpBps(Tcl_Interp *interp);

extern int
CreateConsGaps(Tcl_Interp *interp);

extern int
CreateDpConsBps(Tcl_Interp *interp);

extern int 
RaiseGapLayer(Tcl_Interp *interp, int min_nt, int max_nt);

extern int
NewConsDpBp(Tcl_Interp *interp, pair_t pair, float dp_scale);

extern int
DeleteConsDpBp(Tcl_Interp *interp, pair_t pair);


extern int
ReColorAndResizeConsDpBp(Tcl_Interp *interp, pair_t pair, float dp_scale);

extern int
CreateInfoDp_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);
extern void
GetGapCoords(Tcl_Interp *interp, int nt_pos, char *coords);

extern void
GetRgb(char *rgb, float hue, float lowlimit, float highlimit);

extern void
GetFullCellCoords(pair_t pair, double dp_scale, rect_coord_t *c);

extern void
GetProbCellCoords(float prob, pair_t pair, double dp_scale, rect_coord_t *c);



/* for batch structure prediction  of many files via remote tcl interpreter:
   don't draw anything in dotplot */
#if 0
# define BATCH_NON_DP_GUI
#endif


#endif

