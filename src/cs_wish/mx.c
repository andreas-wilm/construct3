/******************************************************************************
* 
* mx.c - generic (dotplot) matrix routines
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
 *  CVS $Id: mx.c,v 1.5 2004-05-25 13:29:16 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include "mx.h"




#if 0
    #define MX_DEBUG
#endif
 
#define ALLOC_ERROR       fprintf(stderr, "ERROR(%s|%s:%d): no more memory\n", __FILE__, __FUNCTION__, __LINE__); 
#define INVALID_VAL(arg)  fprintf(stderr, "ERROR(%s|%s): invalid value for arg %s\n", __FILE__, __FUNCTION__, (arg)); 


/***   mx_new   ***************************************************************
 *
 *
 */
void **
mx_new(int nr, int nc, size_t es, size_t pes, int type)
{
    void **mx;
    int i, is_a_dp;
    
    /*** paranoia
     *
     */
    if (es<=0)
    {
        INVALID_VAL("es");
        return NULL;
    }
    else if (pes<=0)
    {
        INVALID_VAL("pes");
        return NULL;
    }
    else if (nr<=0)
    {
        INVALID_VAL("nr");
        return NULL;
    }
    else if (type==MX_DEFAULT)
    {
        if (nc<=0)
        {
            INVALID_VAL("nc");
            return NULL;
        }
    }
    
    else if ((type!=MX_DEFAULT)&& (type!=MX_DP) && (type!=MX_DPDIAG))
    {
        INVALID_VAL("type");
        return NULL;
    }
    
    if ((type==MX_DP) || (type==MX_DPDIAG))
        is_a_dp = 1;
    else
        is_a_dp = 0;
    
    
        
    /* Allocate pointers to rows
     */
    /* FIXME: 0 is nonsense for  MX_DP */
    if ((mx=(void **) calloc(nr, pes))==NULL)
    {
        ALLOC_ERROR;
        free(mx);
        return NULL;
    }
        

   /* Allocate rows and set pointers to them
    */
    for (i=0; i<nr; i++)
    {
        /* 0:0 invalid for MX_DP */
        if (type==MX_DP && i==0) 
            continue;
            
        if (is_a_dp)
            mx[i]=(void *) calloc(i+1, es);
        else
            mx[i]=(void *) calloc(nc,  es);
            
        if (mx[i]==NULL)
        {
            ALLOC_ERROR;
            free(mx); /* FIXME: some memory is lost */
            return NULL;
        }
    }

   /* Return pointer to array of pointers to rows
    */
   return mx;
}
/***   mx_new   ***/





/***   mx_free   ****************************************************
 *
 *  Frees a matrix allocated with mx_new
 *
 */
void
mx_free (void **mx, int nr, int type)
{
    int i;
    
    if ((type!=MX_DEFAULT)&& (type!=MX_DP) && (type!=MX_DPDIAG))
        INVALID_VAL("type");
    if (nr<=0)
        INVALID_VAL("nr");

    
    for (i=nr-1; i>=0; i--)
    {
        /* 0:0 invalid for MX_DP */
        if (type==MX_DP && i==0) 
            continue;
        free(mx[i]);
    }
    free(mx);
}
/***   mx_free_matrix   ***/

