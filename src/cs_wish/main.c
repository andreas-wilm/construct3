/******************************************************************************
*
* main.c - main procedures for ConStruct's wish (cs_wish)
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
 *  CVS $Id: main.c,v 1.63 2008-01-13 22:50:54 wilm Exp $
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <tk.h>
#include <tcl.h>

#include <string.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "proj.h"
#include "bp_prob_mat.h"
#include "cs_seqio.h"
#include "utils.h"
#include "consensus.h"
#include "gui.h"
#include "callbacks.h"
#include "move_update.h"
#include "if.h"
#include "opt_struct.h"
#include "cs_imatch.h"
#include "cs_bmatch.h"
#include "subopt_struct.h"


/***   private
 */
#if 0
    #define DEBUG
#endif

static void
SetResourceFile(Tcl_Interp *interp);

static int
LoadProject(Tcl_Interp *interp, char *projfile);

static void
SetAutoPath(Tcl_Interp *interp);

static void
RegisterNewCommands(Tcl_Interp *interp);

static int
LoadProject_Cmd(ClientData clientData, Tcl_Interp *interp,
                          int objc, Tcl_Obj *CONST objv[]);

extern int
Tcl_AppInit(Tcl_Interp *interp);

proj_t      *proj     = NULL;
aln_t       *aln      = NULL;
cs_opts_t    cs_opts;
tk_color_t   tk_color;
char        *cs_install_prefix = NULL;

/*
 ***/


/***   RegisterNewCommands   **************************************************
 *
 * Registers Extensions as Tcl-Commands
 *
 */
void
RegisterNewCommands(Tcl_Interp *interp)
{
    /* FIXME: Doku
     */
    Tcl_CreateObjCommand(interp,
                         "LoadProject", LoadProject_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Proj_Exchange", Proj_Exchange_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Seq_Exchange", Seq_Exchange_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Highlight_DpBp", Highlight_DpBp_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Highlight_AlnNt_From_DpConsBp", Highlight_AlnNt_From_DpConsBp_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Highlight_AlnNt_From_DpBp", Highlight_AlnNt_From_DpBp_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "MoveNtSelection", MoveNtSelection_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Get_BpProb", Get_BpProb_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Get_ConsBpProb", Get_ConsBpProb_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Get_ConsSeq", Get_ConsSeq_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Version", Version_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "LoadAlignment", LoadAlignment_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "LibZ_supported", LibZ_supported_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "OptimalStructure", OptimalStructure_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "ComputeMutualInfo", ComputeMutualInfo_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "CreateInfoDp", CreateInfoDp_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Get_MutInfo", Get_MutInfo_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "CS_Imatch", CS_Imatch_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "CS_Bmatch", CS_Bmatch_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "GetRgb", GetRgb_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "IsGap", IsGap_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "IsBasepair", IsBasepair_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "Degap", Degap_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "ConsBpInfo", ConsBpInfo_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "SumOfPairs", SumOfPairs_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "PrintCsTdProbMat", PrintCsTdProbMat_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "PrintTdProbMat", PrintTdProbMat_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "PrintCsDpm", PrintCsDpm_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "PrintMic", PrintMic_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "SuboptimalPairlist", SuboptimalPairlist_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "SubBacktrack", SubBacktrack_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "PwIdent", PwIdent_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "CsShiftScore", CsShiftScore_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "ConvertAln", ConvertAln_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp,
                         "LibZ_supported", LibZ_supported_Cmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);

}
/***   RegisterNewCommands   ***/




/***   SetAutoPath   **********************************************************
 *
 *
 * Appends CSLIBDIR and DSLIBDIR to auto_path, so that
 * packages can be loaded automatically
 * (compiler macros!)
 */
void
SetAutoPath(Tcl_Interp *interp)
{
    char cmd[1024];
    /*sprintf(cmd, "lappend auto_path %s %s", CSLIBDIR, DSLIBDIR);
      find internal stuff first
    */
    sprintf(cmd, "set auto_path [linsert $auto_path 0 %s %s]", CSLIBDIR, DSLIBDIR);
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR) {
        char buf[1024];
        sprintf(buf, "Couldn't execute \"%s\" :%s\n",
                       cmd, Tcl_GetStringResult(interp));
        Die(buf);
    }
}
/***   SetAutoPath   ***/





/***   LoadProject_Cmd   ******************************************************
 *
 * Frontend to LoadProject (see below)
 * Registered as LoadProject in wish
 * Write \"success\" to interp if LoadProject succeded, \"failed\" otherwise
 */
int
LoadProject_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
    char    *projfile=NULL;    /* argument */
    Tcl_Obj *op;

    /* check arg
     */
    if ( (objc!=2)) {
        Tcl_WrongNumArgs(interp, 1, objv, "?projfile?");
        return TCL_ERROR;
    }


    projfile = Tcl_GetStringFromObj(objv[1], NULL);
    if (projfile == NULL)
        return TCL_ERROR;



    /*** Init Stuff
     */
    if (GetCsOpts(interp, &cs_opts)  == TCL_ERROR)
        return TCL_ERROR;
    if (GetColors(interp, &tk_color) == TCL_ERROR)
        return TCL_ERROR;


    if (LoadProject(interp, projfile)==TCL_ERROR) {
        /* LoadProject is allowed to fail */
        op = Tcl_NewStringObj("failed", -1);
        Tcl_SetObjResult(interp, op);
        return TCL_OK;
    }


    op = Tcl_NewStringObj("success", -1);
    Tcl_SetObjResult(interp, op);

    return TCL_OK;
}
/***   LoadProject_Cmd   ***/





/***   LoadProject   **********************************************************
 *
 * Loads a provided project file
 *
 */
int
LoadProject(Tcl_Interp *interp, char *projfile)
{
    int    i;
    char   cmd[1024];
    float *weight;
    char   buf[1024];
    char  *degapped_seq;
    
    VERBOSE_P("Loading %s\n", projfile);


    /***   load project in namespace
     *    (also calls Proj_Exchange !)
     */
    sprintf(cmd, "cs_proj::read %s", projfile);
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR) {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }



    if ( ! STR_EQ(Tcl_GetStringResult(interp), "OK")) {
        sprintf(buf, "Loading of %s aborted (%s)\n", projfile, cmd);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    if (cs_opts.debug)
        ProjPrint(proj);


    /***   read alignment
     */
    if (aln!=NULL)
        KillAln(aln);


	if ( (aln = ReadAlnFile(proj->aln_file)) == NULL) {
        sprintf(buf, "Couldn't read sequence file %s!", proj->aln_file);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }



    /***  check consistence of read sequences and bp-prob-matrices
     *
     */
    if (proj->num_entries != aln->num_seq) {
        sprintf(buf, "Number of sequences in \"%s\"(=%d) and \
                      number of matrix-files in \"%s\"(=%d) differ!",
                      proj->aln_file, aln->num_seq, proj->file, proj->num_entries);
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /*** setup the gapshift
     */
	for (i=0; i<proj->num_entries; i++)
        aln->seq[i].gapshift = CreateGapShift(aln->seq[i].nt);

	/*** read matrices and setup basepair probabilities
	 */
	degapped_seq = (char*) Scalloc(aln->len+1, sizeof(char));
	for (i=0; i<proj->num_entries; i++) {
        Degap(aln->seq[i].nt, degapped_seq);
        aln->seq[i].bp = SetupBpProbMx(proj->entry[i].f_bp_prob_mat, degapped_seq, BP_PROB_CUTOFF);
        if (aln->seq[i].bp==NULL) {
            /* FIXME free */
            sprintf(buf, "Couldn't read %s\n", proj->entry[i].f_bp_prob_mat);
            CsDpError(interp, buf);
            return TCL_ERROR;
        }
    }


    /*** now it's save to setup the consensus sequence
     */
    weight = GetWeightsAsArray();
    aln->cons_seq = CalcConsensus(aln, weight);


    /***   InitProj
     */
    sprintf(cmd, "InitProj %s", projfile);
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR) {
        sprintf(buf, "Couldn't execute \"%s\": %s", cmd, Tcl_GetStringResult(interp));
        CsDpError(interp, buf);
        return TCL_ERROR;
    }


    /*** Setup the cs_dp gui
     */
    if (SetupGui(interp)==TCL_ERROR) {
        fflush(stdout);
        return TCL_ERROR;
    }

    free(weight);
    free(degapped_seq);

    return TCL_OK;
}
/***   LoadProject   ***/




/***   SetResourceFile   ******************************************************
 *
 */
void
SetResourceFile(Tcl_Interp *interp)
{
    Tcl_SetVar(interp, "tcl_rcFileName", RC_FILE, TCL_GLOBAL_ONLY);
}
/***   SetResourceFile   ***/







/***   Tcl_AppInit   **********************************************************
 *
 * Init of ConStruct Wish
 *
 */
int
Tcl_AppInit(Tcl_Interp *interp)
{
    char  cmd[1024];


    /***   initialize the Tcl and Tk libraries
     */
    if (Tcl_Init(interp)== TCL_ERROR)
        return TCL_ERROR;
    if (Tk_Init(interp) == TCL_ERROR)
        return TCL_ERROR;

    SetResourceFile(interp);

    SetAutoPath(interp);

    RegisterNewCommands(interp);


    /***   execute InitCore
     */
    sprintf(cmd, "InitCore");
    if (Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL)==TCL_ERROR)
        return TCL_ERROR;

    return TCL_OK;
}
/***   Tcl_AppInit   ***/




/***   main   *****************************************************************
 *
 */
int
main(int argc, char **argv)
{
    Tk_Main(argc, argv, Tcl_AppInit);

    return 0;
}
/***   main   ***/

