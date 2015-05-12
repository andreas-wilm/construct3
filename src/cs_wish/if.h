/******************************************************************************
* 
* if.h - main interface for c<->tcl vice versa
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
 *  CVS $Id: if.h,v 1.37 2007-10-22 10:43:23 steger Exp $    
 */



#ifndef IF_H_INCL
#define IF_H_INCL


/* CANVAS
 */
#define CANVAS_DP_NAME          ".dotplot.f_dp.c_dp"
#define CANVAS_ALN_NAME         ".alignment.paned.links.c_al"
#define CANVAS_ALN_NAME3        ".alignment.paned.rechts.c3_al"
#define CANVAS_ARRAYNAME        "c"
#define CANVAS_SIZE_IDX_NAME    "size"
#define CANVAS_VIRSIZE_IDX_NAME "virtual_size"

/* SCALE
 */
#define SCALE_ARRAYNAME        "scale"
#define SCALE_DP_IDX_NAME      "dp"



/*** Write basepair to Tcl
 *   Indices here zero-offset !
 */
extern void
ExportBasepair(Tcl_Interp *interp, char *bp_arrayname,
                   int nt_idx, int nt_partner, float prob);

/*** Write higher basepairs to Tcl
 *   Indices here zero-offset !
 */
void
ExportHigherBasepair(Tcl_Interp *interp, char *bp_arrayname, int nt_idx,
                        int nof_pairs, int *nt_partner, float *pair_prob);

/***   Exports/Imports Seqs from Tcl to C vice versa
 */
extern int
Seq_Exchange_Cmd(ClientData clientData, Tcl_Interp *interp,
                        int objc, Tcl_Obj *CONST objv[]);

extern int
Get_BpProb_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

extern int
Get_MutInfo_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

extern int
Get_ConsBpProb_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

extern int
Get_ConsSeq_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

/***  Get canvas size as double
 */
extern int
GetCanvasDotplotSize(Tcl_Interp *interp, double *dp_size);


/*** Get canvas virtual size as double
 */
extern int
GetCanvasVirtualsize(Tcl_Interp *interp, double *dp_virtsize);

/*** Get dp canvas scale as double
 */
extern int
GetDpScale(Tcl_Interp *interp, double *scale_dp);

extern int
GetColors(Tcl_Interp *interp, tk_color_t *color);


/***   return package name and version number
 */
extern int
Version_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

/***   Import options to c-core
 */
extern int
GetCsOpts(Tcl_Interp *interp, cs_opts_t *cs_opts);

/***   Load alignment and export it to tcl via Seq_Exchange_Cmd
 *     (replaced in cs_dp via load_project)
 */
extern int
LoadAlignment_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);


/***   Return 1 if LibZ is compiled in, 0 otherwise
 */
extern int
LibZ_supported_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);
                           

/*** Print error message to console
 *   and append it to interpResult
 *   calling function must return TCL_ERROR
 */
extern void
CsDpError(Tcl_Interp *interp, char *msg);


/***
 */
extern int
ComputeMutualInfo_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

/* FIXME: implement: PrintBpProbMat_Cmd
*/

extern int
GetRgb_Cmd(ClientData clientData, Tcl_Interp *interp,
           int objc, Tcl_Obj *CONST objv[]);


extern int
NtIndexIsValid(int nt_idx);

extern int
IsGap_Cmd(ClientData clientData, Tcl_Interp *interp,
          int objc, Tcl_Obj *CONST objv[]);
                           
extern int
IsBasepair_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);


extern int
Degap_Cmd(ClientData clientData, Tcl_Interp *interp,
          int objc, Tcl_Obj *CONST objv[]);

extern int
ConsBpInfo_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);

/* FIXME: replace the functions  calling the following with tcl-procs */
extern void
GetThisDpBpTag(int seq_idx, pair_t pair, char *tag);
extern void
GetGeneralDpBpTags(int seq_idx, char *tag);
extern void
GetThisDpConsBpTag(pair_t pair, char *tag);
extern void
GetGeneralDpConsBpTags(char *tag);

/* Sum of pairs cost calculation */
extern int
SumOfPairs_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);

extern int 
PrintTdProbMat_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);
extern int
PrintCsTdProbMat_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);
extern int 
PrintCsDpm_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);
extern int 
PrintMic_Cmd(ClientData clientData, Tcl_Interp *interp,
               int objc, Tcl_Obj *CONST objv[]);
extern int 
ExportSuboptimalPairlist(Tcl_Interp *interp,
                         char *pairlist_arrayname, int *pairlist);

extern int
PwIdent_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);
extern int
CsShiftScore_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);
extern int
ConvertAln_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[]);

#endif
