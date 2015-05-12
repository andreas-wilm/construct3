/******************************************************************************
* 
* cs_dpm.c - Procedures to construct a merged consensus matrix
*            out of weighted thermodynamic and info matrices
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
 *  CVS $Id: cs_dpm.c,v 1.14 2004-05-25 13:29:13 wilm Exp $    
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "consensus.h"
#include "mic.h"

#include "cs_dpm.h"




/***   CsDpm   **************************************************
 *
 * Returns float matrix with merged mic and td.prob values, or NULL on error
 * Accesses aln->mutinfo->mic directly and fetches td.prob via GetConsBpProb
 * 
 * if mic_fac>0 the mutinfo has to be precomputed !
 *
 * data = (prob_fac * prob, if prob > prob_thr)
 *                        +
 *        (mic_fac + (mic/mic_max),   if (mic/mic_max) > mic_thr)
 * 
 *
 * mic_thr is relative, i.e. range 0-1
 *
 * Caller must free (use FreeCsDpm)
 *
 */
float **
CsDpm(float prob_fac, float mic_fac, float prob_thr, float mic_thr)
{
    int     nt_i, nt_j;
    float **mmx;       /* merged matrix (mic and prob) */
    float   mic, prob;
    pair_t pair;
    

    /* sum of factors should be 1
     * THIS SHITTY STUFF DOESN'T WORK ALWAYS, SO TRUST
     */
    /*if ( (prob_fac+mic_fac) < 0.999999 || (prob_fac+mic_fac) > 1.000001)*/
    /* if ( ((int)(prob_fac*FLOAT2INT) + (int)(mic_fac*FLOAT2INT)) != (int)(FLOAT2INT)) 
    {
        ERROR_P("prob-factor (%f -> %d) + mic-factor (%f -> %d) = %f (->%d) != 1.0\n",
                 prob_fac,(int)(prob_fac*FLOAT2INT),
                 mic_fac,(int)(mic_fac*FLOAT2INT),
                 prob_fac+mic_fac,(int)(FLOAT2INT));
        return NULL;
    }
    */

    /* check thresholds
     */
    if (prob_thr<0.0)
    {
        ERROR_P("prob_thr=%f < 0.0\n", prob_thr);
        return NULL;
    }
    if (mic_thr<0.0 || mic_thr>1.0)
    {
        ERROR_P("mic_thr=%f <0.0 or >1.0\n", mic_thr);
        return NULL;
    }


    /* if mic should be added, it already has to be computed
     */
    if (mic_fac>0.0 && (! MIC_IsComputed()) )
    {
        ERROR_P("%s\n", "mic must be precomputed\n");
        return NULL;
    }



    mmx = (float**) mx_new(aln->len, 0, sizeof(float), sizeof(float*), MX_DP);

    
    /***** foreach basepair
     */
    for (nt_i=0; nt_i<aln->len; nt_i++)
    {
        for (nt_j=nt_i+1; nt_j<aln->len; nt_j++)
        {
            SetPair(nt_i, nt_j, &pair);
            mmx[nt_j][nt_i] = 0.0;

            if (mic_fac)
            {
                mic  = GetConsMIC(pair, mic_thr);
                mmx[nt_j][nt_i] +=  (mic_fac*mic);
                /* if (mic > mic_thr): already done by GetConsMIC */
            }

            if (prob_fac)
            {
                prob = GetConsBpProb(aln, pair, prob_thr);
                mmx[nt_j][nt_i] += (prob_fac*prob);
                /* if (prob > prob_thr) : already done by GetConsBpProb */
            }
        }       
    }
    
    
    return mmx;
}
/***   CsDpm   ***/



/***   FreeCsDpm   ************************************************************
 *
 */
void
FreeCsDpm(float **cs_dpm)
{
    mx_free((void**)cs_dpm, aln->len, MX_DP);
}
/***   FreeCsDpm   ***/






/***   PrintCsDpm   *********************************************************
 *
 *
 *
 */
void
PrintCsDpm(char *csseq, float **csdpm, FILE *stream)
{
    int i=0, j=0;
    int seqlen  = 0;


    fprintf(stream, "merged consensus dotplot matrix\n");
    fprintf(stream, "data = (prob_fac * prob, if prob > prob_thr) + (mic_fac + (mic/mic_max), if (mic/mic_max) > mic_thr)\n");

    
    seqlen = strlen(csseq);
    
    /* print sequence nt
     */
    fprintf(stream, "   ");
    for (i=0; i<seqlen ; i++)
        fprintf(stream, "  %3c  ", csseq[i]);
    fprintf(stream, "\n");



    for (i=0; i<seqlen ; i++)
    {
    
   	    fprintf(stream, "%c  ", csseq[i]);
        for (j=0; j<seqlen ; j++)
        {
            if (i>j)
                fprintf(stream, "  ---  ");
            else if (j==i)
                fprintf(stream, "   \\   ");
            else
                fprintf(stream, " %0.3f ", csdpm[j][i]);
        }
        fprintf(stream, "   %d\n", i+1);
    }
    
    fprintf(stream, "   ");
    for (i=0; i<seqlen ; i++)
        fprintf(stream, "  %3d  ", i+1);
    fprintf(stream, "\n");
    
    fflush(stream);
}
/***   PrintCsDpm   ***/

