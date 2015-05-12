/******************************************************************************
* 
* TkCS_DrawStruct.c - 
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
 *  CVS $Id: TkCS_DrawStruct.c,v 1.11 2004-05-25 13:29:21 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <tcl.h>
#include <tk.h>
#include <stdlib.h>
#include <string.h>

#include "drawstruct.h"
#include "TkCS_DrawStruct.h"
#include "CS_DrawStruct.h"
#include "generic_utils.h"


#if 0
    #define DEBUG
#endif




struct pair_prob *pair;
struct xyplot    *my_xyplot;

#define TCLNAME_PIVOT       "ds_core_pivot"
#define TCLNAME_NUCRANGE    "ds_core_nucrange"
#define TCLNAME_SINGLE      "ds_core_single"
#define TCLNAME_TEXT        "ds_core_text"
#define TCLNAME_BACKBONE    "ds_core_backbone"


/***   Drawstructcore_Init   **************************************************
 *
 * Init is called when the package is loaded
 *
 */
int
Drawstructcore_Init(Tcl_Interp *interp)
{
    if (Tcl_InitStubs(interp, "8.3", 0) == NULL)
        return TCL_ERROR;

    Tcl_CreateObjCommand(interp,
                         "cs_drawstructure", TkCS_DrawStructCmd,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_PkgProvide(interp, "drawstructcore", "0.9");
    
    return TCL_OK;
}
/***   Drawstructcore_Init   ***/





 
/***   TkCS_DrawStructCmd   ***************************************************
 *
 */
int
TkCS_DrawStructCmd(ClientData  clientData,  Tcl_Interp  *interp,	   
                   int objc, Tcl_Obj *CONST objv[])
{

   int    i, l, seq_len;
   char   ret_string[200];
   float  coord[11];


    char     *seq, *loopmode;
    Tcl_Obj **listPtr;       /* Tcl_List pointer */
    int       listCtr;       /* Tcl_List counter */
    Tcl_Obj  *dummy_obj, *varname_obj;
    PAR       par;
    char     *dummy;
    struct pair_prob *pair;

   
    if (objc!=4)
    {
        Tcl_WrongNumArgs(interp, 1, objv, "?seq bp_list loopmode?");
        return TCL_ERROR;
    }

   par.quiet	= TRUE;        /* No output */
   par.pivot	= FALSE;       /* No pivot */

   par.single   = TRUE;	       /* remove single base pairs  */
   par.text	    = TRUE;        /* Show the labels			*/
   par.backbone = TRUE;        /* Plot the backbone			*/

   par.nuc	= FALSE;       /* Don't show the nucleotides		*/
   par.numnuc	= 0;       /* No ranges for showing nucleotides	*/
   par.numbers  = 10;      /* Do numbering in steps of 10		*/
   par.bonds	= TRUE;    /* Show the base pair bonds		*/
   par.numangle = 0;       /* No angles for pivot			*/
   par.scf	= 0.0;         /* scaling factor			*/
   par.init	= FALSE;       /* Plotter Init-Sequence			*/
   par.range	= FALSE;
   par.range0	=  1;
   par.range1	=  1;
   par.rangestep= 10;
   par.straight = TRUE;	       /* helices surrounding unsymmetric loops are plotted on a straight line */


    /***   get sequence
     */
    dummy = Tcl_GetStringFromObj(objv[1], NULL);
    
    
    if (dummy == NULL)
        return TCL_ERROR;
    seq = strdup(dummy);
    seq_len = strlen(seq);

    #ifdef DEBUG 
        DEBUG_P("got seq \"%s\" length %d\n", seq, seq_len);
    #endif

    /***   get bp list
     */
    if (Tcl_ListObjGetElements(interp, objv[2], &listCtr, &listPtr) != TCL_OK)
    {
        free(seq);
        return TCL_ERROR;
    }
    
    if (listCtr != seq_len)
    {
        Tcl_AppendResult(interp,
                         "sequence length and bp list length mismatch",
                         (char *) NULL);
        free(seq);
        return TCL_ERROR;
    }
    
    pair = (struct pair_prob *) Scalloc(seq_len+1, sizeof(struct pair_prob));

    for (i=0; i<listCtr; i++)
    {
        if (Tcl_GetIntFromObj(interp, listPtr[i], &pair[i+1].pair) != TCL_OK)
        goto ErrorExit;
        
        if (pair[i+1].pair<1)
            pair[i+1].pair=0;
        #ifdef DEBUG 
            DEBUG_P("got seq \"%s\" length %d\n", seq, seq_len);
        #endif
    }
    
    /***   get loopmode
     */
    dummy = Tcl_GetStringFromObj(objv[3], NULL);
    if (dummy == NULL )
        goto ErrorExit;
    
    loopmode = strdup(dummy);
    if (STR_EQ(loopmode, "straight"))
    {
        par.straight = TRUE;
    }
    else if (STR_EQ(loopmode, "symmetric"))
    {
        par.straight = FALSE;
    }
    else
    {
        WARN_P("%s\n", "loopmode must me straight or symmetric...using symmetric");
        par.straight = FALSE;
    }
    #ifdef DEBUG
        DEBUG_P("got loopmode \"%s\"\n", loopmode);
    #endif
    free(loopmode);
    
    
    /**********   other args
     *
     * not passed to function call: local vars
     *
     */
    varname_obj = Tcl_NewStringObj("THIS IS JUST A DUMMY", -1);
    
    
    /***   pivot = TCLNAME_PIVOT
     */
    Tcl_SetStringObj(varname_obj, TCLNAME_PIVOT, -1);
    dummy_obj = Tcl_ObjGetVar2(interp, varname_obj, NULL, TCL_NAMESPACE_ONLY);    
    
    /* no TCL_LEAVE_ERR_MSG since var may not exist */
    if (dummy_obj!=NULL)
    {
        par.pivot = TRUE;
        #ifdef DEBUG
            DEBUG_P("%s\n", "got pivot");
        #endif
        if (Tcl_ListObjGetElements(interp, dummy_obj, &listCtr, &listPtr) != TCL_OK)
                goto ErrorExit;
        
        for ( l=0; l<listCtr && l<MAXANGLES; l=l+3 )
        {
            if (Tcl_GetIntFromObj(interp, listPtr[l], &par.angle0[l/3]) != TCL_OK)
                goto ErrorExit;
            if (Tcl_GetIntFromObj(interp, listPtr[l+1], &par.angle9[l/3]) != TCL_OK)
                goto ErrorExit;
            if (Tcl_GetDoubleFromObj(interp, listPtr[l+2], &par.angle[ l/3]) != TCL_OK)
                goto ErrorExit;

            par.numangle++;
        }
    }
    else
    {
        #ifdef DEBUG
            DEBUG_P("no pivot from %s\n", Tcl_GetStringFromObj(varname_obj, NULL));
        #endif
    }



    /***   nucrange = TCLNAME_NUCRANGE
     */
    Tcl_SetStringObj(varname_obj, TCLNAME_NUCRANGE, -1);
    dummy_obj = Tcl_ObjGetVar2(interp, varname_obj, NULL, TCL_NAMESPACE_ONLY);
    
    /* no TCL_LEAVE_ERR_MSG since var may not exist */
    if (dummy_obj!=NULL)
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "got nucrange");
        #endif
        par.nuc = TRUE;
        
        if (Tcl_ListObjGetElements(interp, dummy_obj, &listCtr, &listPtr) != TCL_OK)
                goto ErrorExit;
        
        for ( l=0; l<listCtr && l<MAXNUCS; l=l+2 )
        {
            if (Tcl_GetIntFromObj(interp, listPtr[l], &par.nuc0[l/2]) != TCL_OK)
                goto ErrorExit;
            if (Tcl_GetIntFromObj(interp, listPtr[l+1], &par.nuc9[l/2]) != TCL_OK)
                goto ErrorExit;

            par.numnuc++;
        }
    }
    else
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "no nucrange");
        #endif
    }

    /***   single = TCLNAME_SINGLE (bool)
     */
    Tcl_SetStringObj(varname_obj, TCLNAME_SINGLE, -1);
    dummy_obj = Tcl_ObjGetVar2(interp, varname_obj, NULL, TCL_NAMESPACE_ONLY);
    /* no TCL_LEAVE_ERR_MSG since var may not exist */
    if (dummy_obj!=NULL)
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "got single");
        #endif
        if (Tcl_GetBooleanFromObj(interp, dummy_obj, &par.single) != TCL_OK)
            goto ErrorExit;
    }
    else
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "no single");
        #endif
    }
        
    

    /***   text = TCLNAME_TEXT (bool)
     */
    Tcl_SetStringObj(varname_obj, TCLNAME_TEXT, -1);
    dummy_obj = Tcl_ObjGetVar2(interp, varname_obj, NULL, TCL_NAMESPACE_ONLY);
    /* no TCL_LEAVE_ERR_MSG since var may not exist */
    if (dummy_obj!=NULL)
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "got text");
        #endif
        if (Tcl_GetBooleanFromObj(interp, dummy_obj, &par.text) != TCL_OK)
            goto ErrorExit;
    }
    else
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "no text");
        #endif
    }


    /***   backbone = TCLNAME_BACKBONE (bool)
     */
    Tcl_SetStringObj(varname_obj, TCLNAME_BACKBONE, -1);
    dummy_obj = Tcl_ObjGetVar2(interp, varname_obj, NULL, TCL_NAMESPACE_ONLY);
    /* no TCL_LEAVE_ERR_MSG since var may not exist */
    if (dummy_obj!=NULL)
    {
    
        #ifdef DEBUG
             DEBUG_P("%s\n", "got backbone");
        #endif
        if (Tcl_GetBooleanFromObj(interp, dummy_obj, &par.backbone) != TCL_OK)
            goto ErrorExit;
    }
    else
    {
        #ifdef DEBUG
            DEBUG_P("%s\n", "no backbone");
        #endif
    }




/*--------------------------------------------------------------------------*/


/* if ( par.nuc==TRUE && par.numnuc==0 ) par.backbone==FALSE; */
   if ( par.nuc==TRUE && par.numnuc>0  )
   {
       for (i=0; i<par.numnuc; i++)
           if ( par.nuc0[i]>par.nuc9[i] )
           {
               l	   = par.nuc0[i];
               par.nuc0[i] = par.nuc9[i];
               par.nuc9[i] = l;
           }
   }
   if ( par.backbone==TRUE && par.nuc==TRUE && par.numnuc==0 )
   {
       par.numnuc  = 1;
       par.nuc0[0] = 0;
       par.nuc9[0] = 0;
   }

   if ( par.pivot==TRUE && par.numangle>0 )
   {
       for (i=0; i<par.numangle; i++)
           if ( par.angle0[i]>par.angle9[i] )
           {
               l	     = par.angle0[i];
               par.angle0[i] = par.angle9[i];
               par.angle9[i] = l;
           }
   }

    if ( par.range==TRUE )
        par.numbers = FALSE;



/*--------------------------------------------------------------------------*/
  if (!par.quiet)
  {
      if ( par.nuc==TRUE ) {     printf("Show nucleotides (nucrange)\n");
         if ( par.numnuc>0 )
      	 for (i=0; i<par.numnuc; i++)
      			         printf("\tfrom nt %d to nt %d\n", par.nuc0[i],par.nuc9[i]);
      }
      else 		         printf("Don't show nucleotides (nucrange)\n");

      if ( par.pivot==TRUE ) {   printf("Pivot substructure (pivot)\n");
          if ( par.numangle>0 )
      	   for (i=0; i<par.numangle; i++)
      			         printf("\tfrom nt %d to nt %d by %.1f degrees\n", par.angle0[i],par.angle9[i],par.angle[i]);
      }
      else 		         printf("Don't pivot a substructure (pivot)\n");

      if ( par.text==FALSE )     printf("Don't show labels (text)\n");
      else 		         printf("Show labels (text)\n");

      if ( par.numbers==FALSE )  printf("Don't show numbers (label)\n");
      else 		         printf("Show numbers in steps of %d (label)\n", par.numbers);

      if ( par.range==FALSE )    printf("Don't show modified numbers (range)\n");
      else {
      			         printf("Show numbers in steps of %d (range)\n", par.rangestep);
      			         printf("\twith true nt. %d as nt. %d\n",  par.range0, par.range1);
      }

      if ( par.backbone==FALSE ) printf("Don't show backbone (backbone)\n");
      else 		         printf("Show backbone (backbone)\n");

      if ( par.bonds==FALSE )    printf("Don't show bonds (hbonds)\n");
      else 		         printf("Show bonds (hbonds)\n");

      if ( par.single==FALSE )   printf("Remove lonely base pairs (single)\n");
      else 		         printf("Lonely base pairs are ok (single)\n");

      if ( par.straight==TRUE )  printf("Loops are drawn straight (straight)\n");
      else 		         printf("Loops are drawn symmetrically (straight)\n");


      printf("\n");
   }
/*--------------------------------------------------------------------------*/


    #ifdef DEBUG
        DEBUG_P("%s\n", "calling CS_DrawStruct");
    #endif
    my_xyplot = (struct xyplot *)  Scalloc(seq_len, sizeof(struct xyplot));
      
    if (!CS_DrawStruct(pair, seq, &par, coord, my_xyplot))
    {
        Tcl_AppendResult(interp,
                         "ERROR in CS_DrawStruct",
                         (char *) NULL);
        free(my_xyplot);
        goto ErrorExit;
    }
      
    /* a tk-canvas has its origin at the top left edge, so all y-coords are mirrored horizontally */
    sprintf(ret_string, "%f %f %f %f %f %f %f %f %f %f %f %f", coord[ 0],	/* min x of canvas	*/
                                                                      coord[ 3]*(-1.),	/* min y		*/
                                                                      coord[ 2],	/* max x		*/
                                                                      coord[ 1]*(-1.),	/* max y		*/
                                                                      coord[ 4],	/* x of headerline1	*/
                                                                      coord[ 5]*(-1.),	/* y 			*/
                                                                      coord[ 6],	/* x of headerline2	*/
                                                                      coord[ 7]*(-1.),	/* y 			*/
                                                                      coord[ 8],	/* x of headerline3	*/
                                                                      coord[ 9]*(-1.),	/* y 			*/
                                                                      coord[10],	/* character size	*/
                                                                      coord[11]);	/* character width	*/
    Tcl_AppendElement(interp, ret_string);

    for (i=0; i<seq_len; i++)
    {
        sprintf(ret_string, "%4d %4d %c %4d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %4d", 
                                   i+1,
                                       my_xyplot[i].pair,
                                           my_xyplot[i].nuc,
                                              my_xyplot[i].tag,
                                                  my_xyplot[i].nuc_x,
                                                       my_xyplot[i].nuc_y*(-1.),
                                                            my_xyplot[i].back_x1,
                                                                 my_xyplot[i].back_y1*(-1.),
                                                                      my_xyplot[i].back_x2,
                                                                           my_xyplot[i].back_y2*(-1.),
                                                                                my_xyplot[i].bond_x1,
                                                                                     my_xyplot[i].bond_y1*(-1.),
                                                                                          my_xyplot[i].bond_x2,
                                                                                               my_xyplot[i].bond_y2*(-1.),
                                                                                                    my_xyplot[i].line_x1,
                                                                                                         my_xyplot[i].line_y1*(-1.),
                                                                                                              my_xyplot[i].line_x2,
                                                                                                                   my_xyplot[i].line_y2*(-1.),
                                                                                                                        my_xyplot[i].lpos_x,
                                                                                                                             my_xyplot[i].lpos_y*(-1.),
                                                                                                                                  my_xyplot[i].label);
        Tcl_AppendElement(interp, ret_string);
    }

    free(seq); free(pair);
    free(my_xyplot);

    return TCL_OK;
    
    
    ErrorExit:
        free(seq); free(pair);
        return TCL_ERROR;
    
}
/***   TkCS_DrawStructCmd   ***/
