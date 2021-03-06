/******************************************************************************
* 
* proj.h - interface to cs_proj.tcl
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
 *  CVS $Id: proj.h,v 1.8 2004-05-25 13:29:17 wilm Exp $
 */


#ifndef PROJ_H_INCL
#define PROJ_H_INCL

/* Converts proj weights to simple float array for easier handling
 * Caller must free
 */
extern float *
GetWeightsAsArray(void);

/* return sum of all project weights
 */
extern float
GetSumOfWeights();

extern proj_t *
NewProj(char *name, char *file, char *version, char *alnfile, int nseq);

extern void
KillProj(proj_t *p);

extern int
Proj_Exchange_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

extern void
ProjPrint(proj_t *p);

#endif
