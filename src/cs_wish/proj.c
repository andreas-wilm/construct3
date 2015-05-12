/******************************************************************************
*
* proj.c - interface to cs_proj.tcl
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
 *  CVS $Id: proj.c,v 1.17 2004-05-25 13:29:17 wilm Exp $
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <tcl.h>
#include <string.h>


#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "cs_seqio.h"
#include "if.h"
#include "proj.h"



/***   GetWeightsAsArray   ****************************************************
 *
 * Converts proj weights to simple float array for easier handling
 * Caller must free
 */
float *
GetWeightsAsArray(void)
{
    int    i   = 0;
    float *w_a = NULL;

    w_a = (float*) calloc(proj->num_entries, sizeof(float));
    for (i=0; i<proj->num_entries; i++)
        w_a[i]=proj->entry[i].weight;

    return w_a;
}
/***   GetWeightsAsArray   ***/




/***   GetSumOfWeights   ******************************************************
 *
 */
float
GetSumOfWeights(void)
{
    float w_s = 0.0;
    int i;

    for (i=0; i<proj->num_entries; i++)
        w_s += proj->entry[i].weight;

    return w_s;
}
/***   GetSumOfWeights   ***/




/***   ProjPrint   ************************************************************
 *
 * Print out project
 * Usefull for debugging only
 */
void
ProjPrint(proj_t *p)
{
    int i=0;

    fprintf(stdout, "Project %s\n", p->name);
    fprintf(stdout, "\tFile=%s\n", p->file);
    fprintf(stdout, "\tVersion=%s\n", p->version);
    fprintf(stdout, "\tnseq=%d\n", p->num_entries);
    fprintf(stdout, "\talignment=%s\n", p->aln_file);

    for (i=0; i<p->num_entries; i++) {
        fprintf(stdout, "\tentry %d\n", i);
        fprintf(stdout, "\t\tseq_id=%s\n", p->entry[i].seq_id);
        fprintf(stdout, "\t\tf_bp_prob_mat=%s\n", p->entry[i].f_bp_prob_mat);
        fprintf(stdout, "\t\tweight=%f\n", p->entry[i].weight);
        fprintf(stdout, "\t\tfold_cmd=%s\n", p->entry[i].fold_cmd);
    }
    fprintf(stdout, "\n");
    fflush(stdout);

}
/***   ProjPrint   ***/




/***   NewProj   **************************************************************
 *
 * Construct a new proj_t
 * Caller must free with KillProj
 */
proj_t *
NewProj(char *name, char *file, char *version, char *alnfile, int n_entries)
{

	proj_t *p = NULL;

    p = (proj_t*) Scalloc(1, sizeof(proj_t));

    p->name         = strdup(name);
    p->file         = strdup(file);
    p->version      = strdup(version);

    p->num_entries  = n_entries;
    p->entry        = (proj_entry_t*) Scalloc(n_entries, sizeof(proj_entry_t));

    p->aln_file     = strdup(alnfile);

    return p;
}
/***   NewProj   ***/




/***   KillProj   *************************************************************
 *
 * Free a proj_t
 */
void
KillProj(proj_t *p)
{
    int i=0;

    free(p->aln_file);

    for (i=0; i<p->num_entries; i++) {
        free(p->entry[i].seq_id);
        free(p->entry[i].f_bp_prob_mat);
        free(p->entry[i].fold_cmd);
    }
    free(p->entry);

    free(p->version);
    free(p->file);
    free(p->name);
    free(p);
}
/***   KillProj   ***/




/***   Proj_Exchange_Cmd   ******************************************************
 *
 *
 * Must be called from inside namespace ! (TCL_NAMESPACE_ONLY)
 *
 *
 * FIXME: allow import!
 */
int
Proj_Exchange_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{

    char *projfile, *projname, *projversion;
    char *alnpath=NULL;
    int  nseq=0;

    char *seq_id_arrayname, *bpmat_arrayname;
    char *weight_arrayname, *foldcmd_arrayname;

	char    *seq_id=NULL;
	char    *f_bp_prob_mat=NULL;
	double   weight;
	char    *fold_cmd=NULL;
    char    msg[1024];

    Tcl_Obj *get_obj = NULL;
    int     i=0;
    char    idx[1024];


    if ( (objc!=10)) {
        Tcl_WrongNumArgs(interp, 1, objv, "?projfile projname projversion alignment nseq \
                                            seq_id_arrayname bpmat_arrayname \
                                            weight_arrayname foldcmd_arrayname?");
        return TCL_ERROR;
    }


    projfile = Tcl_GetStringFromObj(objv[1], NULL);
    if (projfile == NULL) {
        CsDpError(interp, "can´t get projfile");
        return TCL_ERROR;
    }
    projname = Tcl_GetStringFromObj(objv[2], NULL);
    if (projname == NULL) {
        CsDpError(interp, "can´t get projname");
        return TCL_ERROR;
    }
    projversion = Tcl_GetStringFromObj(objv[3], NULL);
    if (projversion == NULL) {
        CsDpError(interp, "can´t get projversion");
        return TCL_ERROR;
    }
    alnpath = Tcl_GetStringFromObj(objv[4], NULL);
    if (alnpath == NULL) {
        CsDpError(interp, "can´t get alnpath");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[5], &nseq)!=TCL_OK) {
        CsDpError(interp, "can´t get nseq");
        return TCL_ERROR;
    }

    seq_id_arrayname = Tcl_GetStringFromObj(objv[6], NULL);
    if (seq_id_arrayname == NULL) {
        CsDpError(interp, "can´t get seq_id_arrayname");
        return TCL_ERROR;
    }
    bpmat_arrayname = Tcl_GetStringFromObj(objv[7], NULL);
    if (bpmat_arrayname == NULL) {
        CsDpError(interp, "can´t get bpmat_arrayname");
        return TCL_ERROR;
    }
    weight_arrayname = Tcl_GetStringFromObj(objv[8], NULL);
    if (weight_arrayname == NULL) {
        CsDpError(interp, "can´t get weight_arrayname");
        return TCL_ERROR;
    }
    foldcmd_arrayname = Tcl_GetStringFromObj(objv[9], NULL);
    if (foldcmd_arrayname == NULL) {
        CsDpError(interp, "can´t get foldcmd_arrayname");
        return TCL_ERROR;
    }

    /* reentrance */
    if (proj!=NULL)
        KillProj(proj);

    proj = NewProj(projname, projfile, projversion, alnpath, nseq);


    for (i=0; i<nseq; i++) {
        sprintf(idx, "%d", i+1);

        /***   get seq_id
         */
        get_obj = Tcl_GetVar2Ex(interp, seq_id_arrayname, idx,
                       TCL_LEAVE_ERR_MSG | TCL_NAMESPACE_ONLY);
        if (get_obj==NULL) {
            sprintf(msg, "can't get %s(%s)", seq_id_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }
        seq_id = Tcl_GetStringFromObj(get_obj, NULL);
        if (seq_id == NULL) {
            sprintf(msg, "can´t read %s(%s)\n", seq_id_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }

        /***   get f_bp_prob_mat
         */
        get_obj = Tcl_GetVar2Ex(interp, bpmat_arrayname, idx,
                       TCL_LEAVE_ERR_MSG | TCL_NAMESPACE_ONLY);
        if (get_obj==NULL) {
            sprintf(msg, "can't get %s(%s)\n", bpmat_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }
        f_bp_prob_mat = Tcl_GetStringFromObj(get_obj, NULL);
        if (seq_id == NULL) {
            sprintf(msg, "can´t read %s(%s)\n", bpmat_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }

        /***   get weight
         */
        get_obj = Tcl_GetVar2Ex(interp, weight_arrayname, idx,
                       TCL_LEAVE_ERR_MSG | TCL_NAMESPACE_ONLY);
        if (get_obj==NULL) {

            sprintf(msg, "can't get %s(%s)\n", weight_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }
        if (Tcl_GetDoubleFromObj(interp, get_obj, &weight)==TCL_ERROR) {
            sprintf(msg, "can't read %s(%s)\n",  weight_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }

        /***   get fold_cmd
         */
        get_obj = Tcl_GetVar2Ex(interp, foldcmd_arrayname, idx,
                       TCL_LEAVE_ERR_MSG | TCL_NAMESPACE_ONLY);
        if (get_obj==NULL) {
            sprintf(msg, "can't get %s(%s)\n", foldcmd_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }
        fold_cmd = Tcl_GetStringFromObj(get_obj, NULL);
        if (fold_cmd == NULL) {
            sprintf(msg, "can´t read %s(%s)\n", foldcmd_arrayname, idx);
            CsDpError(interp, msg);
            return TCL_ERROR;
        }


        /***   and store
         */
        proj->entry[i].seq_id        = strdup(seq_id);
        proj->entry[i].f_bp_prob_mat = strdup(f_bp_prob_mat);
        proj->entry[i].weight        = (float) weight;
        proj->entry[i].fold_cmd      = strdup(seq_id);
    }

    return TCL_OK;
}
/***   Proj_Exchange_Cmd   ***/



