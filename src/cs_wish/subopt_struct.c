/******************************************************************************
* 
* subopt_struct.c - dynamic programming routines for suboptimal structure prediction
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
 *  CVS $Id: subopt_struct.c,v 1.11 2007-10-22 10:43:23 steger Exp $    
 */


/**********
 * 
 *  Most functions are taken completely from ConStruct 2.1.
 *  A reimplementation would be nice
 * 
 *
**********/


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>


#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "if.h"
#include "cs_dpm.h"
#include "struct.h"

#include "subopt_struct.h"


/***   private
 */

#if 0
    #define DEBUG
#endif


#define AREA(i,j) s_area[((i)-1)*L - (i)*((i)-1)/2 + j-1]
#define V(i,j)    v[((i)-1)*L - (i)*((i)-1)/2 + j-1]
#define W(i,j)    w[((i)-1)*L - (i)*((i)-1)/2 + j-1]
#define VL(i,j)   vl[((i)-1)*L - (i)*((i)-1)/2 + j-1]
#define A(i,j)    a[((i)-1)*L - (i)*((i)-1)/2 + j-1]


int *s_area, *v, *w;
int L;      /* Dimension of arrays area, v, vl, w */
float norm; /*2.*(float)INT_MAX/(float)aligned_len; */

static int *
SearchStruct(int numStruct, int aligned_len);
static int *
SubOptPairList (float *sum_aligned_mat, int subdepth, int aligned_len);
static void
ReadArea(float *sum_aligned_mat, int aligned_len);
static void
SubOptimize (int aligned_len);
/*
static void
PrintOutVars (int *a, int length, char *name);
*/
static void
SetPairProbProb (pair_prob_t *pair, int i, int j, int N);
static void
SubBacktrack (pair_prob_t *pair, int i, int j, int aligned_len);
 
/*
 ***/




/***   PrintOutVars   *********************************************************
 *
 *   Print upper half contents of an integer array to StdOut
 *   Assumption: each value is printable in %2d.
 *   The diagonal i:i is marked by asterisks.
 *
 *   Input:
 *        a:        pointer to two-dimensional array a
 *        length:   size of array
 */
void
PrintOutVars (int *a, int length, char *name)
{
	int i, j, matdim, dummy;
    
	matdim = 2 * length;
	
    printf("SUBOPT ARRAY\n");
	
	printf("%s:\n          ", name);
	for (j=1;j<=matdim; j++)
		printf("%10d", j);
	printf (" j\n");
	for (i=1; i<=matdim; i++)
    {
		printf("%10d", i);
		for (j=1; j<=matdim; j++)
         {
			if (i==j)
            {
				printf("%10s", "*");
			}
            else
            {
				if (i>=j)
                {
					printf("%10s", " ");
				}
                else
                {
					dummy = A(i,j);
					printf("%10d", dummy);
				}
			}
		}
		printf("\n");
	}
	printf("  i\n");
}
/***   PrintOutVars   ***/




/***  SubOptimize   ***********************************************************
 *
 *  Fill in max pairings
 *
 *   Input:
 *   		s_area: 		array of consensus pairing probability
 *   		v: 			array of max pairing with bp(i:j)
 *   		w: 			array of max pairing from i to j
 *   		N:				length of sequence
 *			MIN_HP_SIZE:	minimum hairpin-size
 *    Output:
 *   		v
 *   		w
 */
void
SubOptimize (int aligned_len)
{
	int i, ip1, ip2, j, jm1, jm2, imax, imin, k;
	int M, N, loop, dummy, dummy1, dummy2, wdummy, bd;
	int *il, *bo;
    int **ww_mx;

	N  = aligned_len;
	il = (int *) Scalloc(N+1, sizeof(int));
	bo = (int *) Scalloc(N+1, sizeof(int));
	ww_mx = (int**) mx_new(N+1+1, 1, sizeof(int), sizeof(int*), MX_DEFAULT);
    
	for (j=0; j<N+1; j++)
    {
		il[j]=0;
		bo[j]=0;
	}
	for (j=0; j<=1; j++)
    {
		for (i=1; i<=N; i++)
			ww_mx[i][j]=0;
	}

	M = 2*N;
    for (j=MIN_HP_SIZE+2; j<=M; j++)
    {
      jm1 = j-1;
      jm2 = j-2;
		if (j<=N)		/* region I */
        {
			imax = j - MIN_HP_SIZE - 1;
			imin = 0;
		}
        else			/* region E */
        {
			imax = N;
			imin = j - N + MIN_HP_SIZE;
		}
		bd=0;
		for (i=imax; i>imin; i--)    /* ww_mx(i,j) == ww_mx[j*(N+1)+i] */
			ww_mx[i][0]=ww_mx[i][1];
		
		ww_mx[imax+1][0] = 0;

        for (i=imax; i>imin; i--)
        {
			ip1 = i+1;
			ip2 = i+2;
         if (AREA(i,j)==0)
         {
            V(i,j) = 0;
         }
        /* assume i:j form a pair;
           search for an optimal neighboring "loop", which might have size 0.
         */
         else
         {
				loop = 0;
				if (ip2 <= N)
                {
					dummy1 = il[i];
					dummy2 = V(ip2,jm2);

					il[i] = MAX(dummy1,dummy2);			/* internal loop */

					dummy1 = il[i];
					dummy2 = il[ip1];
					il[i]  = MAX(dummy1,dummy2);			/* il is stored for the next i (and j) */

					dummy2 = V(ip2,jm1);
					bd = MAX(bd,dummy2); 					/* bulge (down column (j-1));
																	bd is stored for the next i */
					dummy1 = il[i];
					loop = MAX(dummy1,bd);
				}
				if (ip1 <= N)
                {
					dummy1 = bo[i];
					dummy2 = V(ip1,jm2);
					bo[i]  = MAX(dummy1,dummy2);			/* bulge (over row (i+1)) */
				
					dummy1 = bo[i];
					loop   = MAX(loop,dummy1);

					dummy2 = V(ip1,jm1);
					loop   = MAX(loop,dummy2);				/* stack */

					if (ip2<=jm1)
                    {
						dummy2 = W(ip1,ip1)+W(ip2,jm1);
						loop   = MAX(loop,dummy2);			/* closed bifurcation */
					}

					dummy2 = ww_mx[ip1][0];
					loop   = MAX(loop,dummy2);
				}
				V(i,j) = loop + AREA(i,j);				/* the loop has to be closed by the current pair i:j */

			}

			wdummy = 0;
			for (k=ip1; k<jm1; k++)					/* open bifurcation */
            {
				dummy  = W(i,k) + W(k+1,j);
				wdummy = MAX(wdummy,dummy);
			}
			ww_mx[i][1] = wdummy;
			dummy2 = V(i,j);
			wdummy = MAX(wdummy,dummy2);					/* bp i:j */
			
			if (ip1 > N)
				dummy2 = 0;
		    else
				dummy2 = W(ip1,j);
			
			wdummy  = MAX(wdummy,dummy2);					/* dangling i */
			
			dummy2  = W(i,jm1);
			W(i,j) = MAX(wdummy,dummy2);					/* dangling j */
		}
   }

	/* Copy array w from included region I to I' */
	for (j=1; j<=N; j++)
    {
		for (i=j-MIN_HP_SIZE-1; i>0; i--)
			W(i+N,j+N) = W(i,j);
	}
    mx_free((void**)ww_mx, N+1+1, MX_DEFAULT);

	free(il);
	free(bo);
}
/***  SubOptimize   ***/




/***   ReadArea   *************************************************************
 *
 */
void
ReadArea(float *sum_aligned_mat, int aligned_len)
{
	int i,j, matdim, dummy;	
    
	matdim = aligned_len+1;
    
	for (i=1; i<=aligned_len; i++)
    {
		for (j=1; j<=i; j++)
        {
			dummy = 	(int)(norm*sum_aligned_mat[i*matdim+(matdim-j)]);
			if (dummy!=0)
				AREA(j,i) = dummy;
		}
	}
	return;
}
/***   ReadArea   ***/




/***   SearchStruct   ********************************************************
 *
 *  Find Vmax = MAX( V( I,J )+V( J,I ))
 *    and all base pairs which belong to "the" structure with
 *    this "energy". Mark all 'used' base pairs => these base
 *    pairs may not be used for construction of a structure with
 *    lower energy.
 *
 *  Input:
 *     s_area:	  array of consensus pairing probability
 *     v:         array of max pairing with bp(i:j)
 *     N:         length of sequence
 *     NumStruct: number of structures to be found
 *            If less optimal structures are found than 
 *            requested for, an error is written but the pairing 
 *            list is returned
 *  Output:
 *        returns a list of lists with starting pairs i, j and 
 *        their max pairing
 */
int *
SearchStruct(int numStruct, int aligned_len)
{
	int N, nMax, m, k, i, j;
	int vTotal, vMax, numBp;
	int *vl; 
	int *iMax, *jMax, *pairList;
	int pairListLength = 0;	/* Number of pairs found */

	N = aligned_len;
	vl = (int *) Scalloc((L*L+L)/2, sizeof(int));

	/* copy v -> vl */
	for (i=1; i<=L; i++)
    {
		for (j=i; j<=L; j++)
			VL(i,j) = V(i,j);
	}
	
	nMax =     N/2;		/* dimension of arrays iMax and jMax */
	iMax =     (int *)Scalloc(nMax,sizeof(int));
	jMax =     (int *)Scalloc(nMax,sizeof(int));
	pairList = (int *)Scalloc(4*(N*N/2)+1,sizeof(int));

	for (m=0; m<numStruct; m++)
    {
		vTotal = -1;
		vMax   =  0;
		numBp  =  0;

		for (j=1;j<=N; j++)
        {
			for (i=1;i<=j-MIN_HP_SIZE; i++)
            {
				if ( VL(i,j)>0 && ( vTotal=(VL(i,j) + VL(j,i+N)-AREA(i,j)) )>0 )
                {
					if (vTotal > vMax)
                    {
						numBp 		 = 1;
						vMax  		 = vTotal;
						iMax[numBp]  = i;
						jMax[numBp]  = j;
					}
                    else if (vTotal == vMax)
                    {
						if (numBp < nMax)
                        {
							numBp++;
							iMax[numBp] = i;
							jMax[numBp] = j;
						}
                        else
                        {
							ERROR_P("Overflow of Search matrix\nMaximum(%d) = %d bp\n", m, vMax);
						}
					}
				}
			}			/* for j */
		}				/* for i */

		if (numBp==0)
        {
			WARN_P("%s\n", "Found less than requested suboptimal structures");
			if (pairListLength>0)
            {
                free(vl);
                free(iMax);
	            free(jMax);
				return pairList;
			}
            else
			{
                free(vl);
                free(iMax);
	            free(jMax);
                free(pairList);
            	return NULL;
		    }
        }
        else
        {
			for (k=1; k<=numBp; k++)
            {
				/*	lappend pairList [list $iMax(1) $jMax(1) $sum] */
				pairList[pairListLength*4+1] = iMax[k];
				pairList[pairListLength*4+2] = jMax[k];
				pairList[pairListLength*4+3] = vMax;
				pairList[pairListLength*4+4] = m;
				pairListLength++;
				pairList[0] = pairListLength;

				/* mark energies (= base pairs) with energy vMax as used */			

				VL(iMax[k],jMax[k])   = 0;
				VL(jMax[k],iMax[k]+N) = 0;
			}
		}
	}
    free(vl);
    free(iMax);
	free(jMax);

	return pairList;
}
/***   SearchStruct   ***/




/***   SubOptPairList   *******************************************************
 *
 * Produce a list of pairs for suboptimal structure
 *
 */
int *
SubOptPairList (float *sum_aligned_mat, int subdepth, int aligned_len)
{
	int i, j;
	int *pairlist;
    
    /* init globals
     */
	if (s_area == NULL)
    {
		s_area = (int *) Scalloc((L*L+L)/2, sizeof(int));
		v      = (int *) Scalloc((L*L+L)/2, sizeof(int));
		w      = (int *) Scalloc((L*L+L)/2, sizeof(int));
	}
    /* re-init
     */
    else
    {
		free(s_area);
		free(v);
		free(w);
		s_area = (int *) Scalloc((L*L+L)/2, sizeof(int));
		v      = (int *) Scalloc((L*L+L)/2, sizeof(int));
		w      = (int *) Scalloc((L*L+L)/2, sizeof(int));
	}


	for (i=1; i<=L; i++)
    {
		for (j=i; j<=L; j++)
        {
			AREA(i,j) = 0;
			V(i,j)    = 0;
			W(i,j)    = 0;
		}
	}

	ReadArea(sum_aligned_mat, aligned_len);

    /* s_area is doubled to allow for an easy access
	   of x(i,j) AND x(j,i) with 1<=i<j<=N
     */

	for (i=1; i<=aligned_len; i++)
    {
		for (j=i; j<=aligned_len; j++)
			AREA(j,aligned_len+i) = AREA(i,j);
	}


    /********************
     */
	SubOptimize(aligned_len);
    /*
    ********************/

    #ifdef DEBUG
    	PrintOutVars(w, aligned_len, "w");
    #endif
    
    /********************************
     */
	pairlist = SearchStruct(subdepth, aligned_len);
    /*
     ********************************/
    
    #ifdef DEBUG
    	DEBUG_P("\n\nPairlist (aligned_len=%d):\n", pairlist[0]);
    	for (i=0; i<pairlist[0]; i++)
        {
		    fprintf(stdout, "%2d: %4d,%4d, %d, %d. best struct.\n",
                             i+1, pairlist[i*4+1], pairlist[i*4+2],
                                  pairlist[i*4+3], pairlist[i*4+4]);
	    }
    #endif

	return pairlist;
}
/***   SubOptPairList   ***/




/***   SuboptimalPairlist_Cmd   *************************************************
 *
 *
 * return FIXME
 *
 */
int
SuboptimalPairlist_Cmd(ClientData clientData, Tcl_Interp *interp,
                           int objc, Tcl_Obj *CONST objv[])
{
	float     *sum_aligned_mat;
	int        i=0, j=0;
    int        mat_dim = 0;
    int       *pairlist;
    int        export_success;
    float    **mcm_mx; /* the merged consensus matrix, out of mic and td.prob */
    /* arguments */
    char        *pairlist_arrayname;
    double      prob_fac, prob_thr; 
    double      mic_fac, mic_thr;   

    int         sub_depth;

    /***   parse args
     *
     */
    if (objc!=7)
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                         "?<pairlist_arrayname> <prob_fac> <mic_fac> <prob_thr> <mic_thr> <sub_depth>?");
        return TCL_ERROR;
    }
    pairlist_arrayname = Tcl_GetStringFromObj(objv[1], NULL);
    if (pairlist_arrayname == NULL)
    {
        CsDpError(interp, "Can´t get pairlist_arrayname");
        return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, objv[2], &prob_fac)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &mic_fac) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[4], &prob_thr)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[5], &mic_thr) ==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[6], &sub_depth)==TCL_ERROR)
        return TCL_ERROR;
        
    #ifdef DEBUG 
        DEBUG_P("pairlist_arrayname=%s\n", pairlist_arrayname);
        DEBUG_P("prob_fac=%f\n",           prob_fac);
        DEBUG_P("mic_fac=%f\n",            mic_fac);
        DEBUG_P("prob_thr=%f\n",           prob_thr);
        DEBUG_P("mic_thr=%f\n",            mic_thr);
        DEBUG_P("sub_depth=%d\n",          sub_depth);
    #endif
    
    
    /* Get the merged consensus dotplot matrix
     */
    mcm_mx = CsDpm((float) prob_fac, (float) mic_fac,
                (float) prob_thr, (float) mic_thr);
    if (mcm_mx==NULL)
    {
        CsDpError(interp, "CsDpm failed");
        return TCL_ERROR;
    }
    #ifdef DEBUG
        DEBUG_P("Merged matrix (prob_fac=%f, mic_fac=%f)\n", (float)prob_fac, (float)mic_fac);
        PrintCsDpm(aln->cons_seq->nt, mcm_mx, stdout);
    #endif
    
    
    
    /* setup globals
     */
	mat_dim = aln->len+1;
    L = 2*aln->len;
    sum_aligned_mat = (float *) Scalloc(pow(mat_dim,2), sizeof(float));    
    /* FIXME: optstruct uses FLOAT_2_INT */
    norm   = 2.*(float)INT_MAX/(float)aln->len;
    
    
    /* we use cs2 code, so we have to map the matrices/indices to sum_aligned_mat */
    
    for (i=1; i<=aln->len; i++)
    {
		for (j=1; j<=i; j++)
        {
            float cs_prob = 0.0;
            
            if (i==j)
            {
                sum_aligned_mat[i*mat_dim+(mat_dim-j)] = 0.0;                
                continue;
            }

            cs_prob = mcm_mx[i-1][j-1];

            sum_aligned_mat[i*mat_dim+(mat_dim-j)] = cs_prob;
        }
    }
    

    pairlist = SubOptPairList(sum_aligned_mat, sub_depth, aln->len);

    /*
     * bi             = pairlist[i*4+1]
     * bj             = pairlist[i*4+2]
     * prob           = pairlist[i*4+3]
     * maxbpctr       = pairlist[i*4+4]+1
     * pairlistlength = pairlist[0]
     */
    
    /* export (also returns the max bp count */
    export_success = ExportSuboptimalPairlist(interp, pairlist_arrayname, pairlist);
    
    
    free(sum_aligned_mat);
    FreeCsDpm(mcm_mx);
    
    if (export_success==TCL_ERROR)
    {
        CsDpError(interp, "ExportSuboptimalPairlist failed");
        return TCL_ERROR;
    }

    
    return TCL_OK;    
}
/***   SuboptimalPairlist_Cmd   ***/






/***   SubBacktrack_Cmd   ***************************************************
 *
 * Backtrack from selected basepair
 * call after SuboptimalPairlist_Cmd etc.
 *
 */
int
SubBacktrack_Cmd(ClientData clientData, Tcl_Interp *interp,
                              int objc, Tcl_Obj *CONST objv[])
{
	pair_prob_t *pair;
	int          nt_idx;
    /* Arguments */
    int nt_i, nt_j;
    char *basepair_arrayname;
    int allow_single_bp;

    

    /***   parse args
     *
     */
    if (objc!=5)
    {
        Tcl_WrongNumArgs(interp, 1, objv,
                         "?<bpp_array_name> <nt_i> <nt_j> <allow_single_bp>?");
        return TCL_ERROR;
    }

    basepair_arrayname = Tcl_GetStringFromObj(objv[1], NULL);
    if (basepair_arrayname == NULL)
    {
        CsDpError(interp, "Can´t get basepair_arrayname");
        return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[2], &nt_i)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetIntFromObj(interp, objv[3], &nt_j)==TCL_ERROR)
        return TCL_ERROR;
    if (Tcl_GetBooleanFromObj(interp, objv[4], &allow_single_bp)==TCL_ERROR)
        return TCL_ERROR;

    
    pair = (pair_prob_t *) Scalloc(aln->len+1, sizeof(pair_prob_t));


    #ifdef DEBUG_0
	    P_DEBUG("basepair_arrayname\n", basepair_arrayname);
	    P_DEBUG("start bp = %d:%d\n", nt_i, nt_j);
    #endif


    
    /* Since we start from an arbitrary basepairs we need to
     *  double the sequence/matrix/backtracking
     *  otherwise we would stop at this particular basepair
     *  and leave out the rest
     *
     * Und unsere deutschsprachigen Leser schlagen bitte in
     *  Gerhard Steger     
     *  Bioinformatik. Methoden zur Vorhersage von RNA- und Proteinstruktur.
     *  ISBN: 3764369515
     * nach!
     */

    /******************************
     * 1. round from nt_i to nt_j
     */
	SubBacktrack(pair, nt_i, nt_j, aln->len);
     /*
      * 2. round from nt_j to nt_i
      */
	SubBacktrack(pair, nt_j, nt_i, aln->len);
    /*
     ******************************/

    
    /* old routines use unit-offset
       -> transform to zero-offset
    */
    for(nt_idx=1; nt_idx<=aln->len; nt_idx++)
    {
        pair[nt_idx-1]      = pair[nt_idx];
        pair[nt_idx-1].pair = pair[nt_idx-1].pair-1;
    }


    
    /*****   remove single basepairs if requested
     */
    if (!allow_single_bp)
        RemoveSingleBps(pair, aln->len);

    
    #ifdef DEBUG
	    DEBUG_P("Pairs for suboptimal structure starting at %d:%d\n", nt_i, nt_j);
	    for (nt_idx=0; nt_idx<aln->len; nt_idx++)
			DEBUG_P("%d <-> %d\n", nt_idx, pair[nt_idx].pair);
    #endif


    /*****   export calculated structure
     */
    for(nt_idx=0; nt_idx<aln->len; nt_idx++)
        ExportBasepair(interp, basepair_arrayname, nt_idx, pair[nt_idx].pair, pair[nt_idx].prob);
    
    
	/*
	    FIXME tcl_SetSquigglesFile("test.squ" , aligned_seq, pair, interp, 0);
    */
    
    
	free(pair);
	
	return TCL_OK;
}
/***   SubBacktrack_Cmd   ***/




/***   SubBacktrack   *********************************************************
 *
 *  Find optimal structure with bp i:j
 *
 *   Input:
 *        ptrArea:	array of consensus pairing probability
 *        ptrV:		array of max pairing with bp(i:j)
 *        ptrW:		array of max pairing from i to j
 *        N:			length of sequence
 *        i:		   1<=i<=N
 *        j:		   1<=j<=N
 *                        
 *   Output:
 *        pair:		the pairing array is written
 */
void
SubBacktrack (pair_prob_t *pair, int i, int j, int aligned_len)
{
	int N, pos, ip1, jm1, jm2, l, k, isPair;
	int *ipos, *jpos;
	int max_i0, max_j0,
		max_i1, max_j1,
		max_i2, max_j2,
		max_i3, max_j3,
		max_v,  max_w;		

	ipos   = (int *) Scalloc(MAX_NUMBIF,sizeof(int));
	jpos   = (int *) Scalloc(MAX_NUMBIF,sizeof(int));

	N = aligned_len;

	/*	i:j is a base pair; we don't search for open regions	*/
	SetPairProbProb(pair,i,j,N);

	/* remember that this is the start */
	isPair = 1;

	/* backtrack from the optimal base pair i:j
		either in upstream direction (1<=i<=j<=N)
		in downstream direction (i<=N,1<=j or 1<=j<=i<=N) */
	if (i>j) {
	   j = j+N;
	}

	/* push the pair to the stack */
	pos			= 1;
	ipos[pos-1] = i;
	jpos[pos-1] = j;

	while (pos>0) {
		/* get next closed (=> find a bp from v) or
          open region (=> find bp in w) from stack
         */
		pos--;
		i = ipos[pos];
		j = jpos[pos];
		ip1 = i+1;
		jm1 = j-1;
		jm2 = j-2;

		max_i0 = 0;
		max_j0 = 0;
		max_i1 = 0;
		max_j1 = 0;
		max_v  = 0;
		
		/* allow for a bp i:j=N:(N+1)=N:1 or a hp loop */
		if ((j>N) || ((abs(j-i)-1)>MIN_HP_SIZE))
        {
			/* if a bp is possible */
			if (AREA(i,j)>0)
            {
            	/* hairpin or internal or bulge loop or bp: i:j is a bp */
				for (l=i+2; l<j; l++)
                {
					for (k=ip1; k<l; k++)
                    {
						if ((V(i,j)-AREA(i,j)-V(k,l)==0) && V(k,l)>max_v)
                        {
							max_i0 = k;
							max_j0 = l;
							max_i1 = 0;
							max_j1 = 0;
							max_v  = V(k,l);
						}
					}
				}
               	/* closed bifurcation with bp i:j */
				for (k=ip1; k<jm2; k++)
                {
					if ((V(i,j)-AREA(i,j)-W(ip1,k)-W(k+1,jm1)==0) \
                        && \
						W(ip1,k)+W(k+1,jm1)>max_v)
                    {
						max_i0 = ip1;
						max_j0 = k;
						max_i1 = k+1;
						max_j1 = jm1;
						max_v  = W(ip1,k)+W(k+1,jm1);
					}
				}
			}
			/* if this is the first i:j it has to be a pair */
			if (isPair==1)
            {
				ipos[pos] = max_i0;
				jpos[pos] = max_j0;
				pos++;
				if (max_i1>0)
                {
					ipos[pos] = max_i1;
					jpos[pos] = max_j1;
					pos++;
				}
                else
                {
					if (max_i0>0 && max_j0>0 && AREA(max_i0,max_j0)>0)
						SetPairProbProb(pair,i,j,N);
				}
				isPair = 0;
			}
            else
            {
				/* search for an open region */
				max_i2 = 0;
				max_j2 = 0;
				max_i3 = 0;
				max_j3 = 0;
				max_w  = max_v;
				if (W(i,j)>0)
                {
               		/* open bifurcation */
					for (k=jm2; k>i; k--)
                    {
						if ((W(i,k)+W(k+1,j)-W(i,j)==0) && (W(i,j)>max_w))
                        {
							max_i2 = i;
							max_j2 = k;
							max_i3 = k+1;
							max_j3 = j;
							max_w  = W(i,j);
						}
					}
					if ((W(i,j)==W(ip1,j)) && (W(ip1,j)>max_w))
                    {
						max_i2 = ip1;
						max_j2 = j;
						max_i3 = 0;
						max_j3 = 0;
						max_w  = W(ip1,j);
					}
					if ((W(i,j)==W(i,jm1)) && (W(i,jm1)>max_w))
                    {
						max_i2 = i;
						max_j2 = jm1;
						max_i3 = 0;
						max_j3 = 0;
						max_w  = W(i,jm1);
					}
				}
        		/* it's an open region */
				if (max_w>max_v)
                {
					ipos[pos] = max_i2;
					jpos[pos] = max_j2;
					pos++;
					if (max_i3>0)
                    {
						ipos[pos] = max_i3;
						jpos[pos] = max_j3;
						pos++;
					}
				}
                /* it's a closed region */
                else
                {
					ipos[pos] = max_i0;
					jpos[pos] = max_j0;
					pos++;
					if (max_i1>0)
                    {
						ipos[pos] = max_i1;
						jpos[pos] = max_j1;
						pos++;
						if (    i >0 &&     j >0 && AREA(    i,     j )>0)
                            SetPairProbProb(pair,    i,     j ,N);
					}
                    else
                    {
						if (max_i0>0 && max_j0>0 && AREA(max_i0,max_j0)>0)
                            SetPairProbProb(pair,max_i0,max_j0,N);
						if (    i >0 &&     j >0 && AREA(    i,     j )>0)
                            SetPairProbProb(pair,    i,     j ,N);
					}
				}
			}
		}
	}
	free(ipos);
	free(jpos);
    
   return;
}
/***   SubBacktrack   ***/



/***   SetPairProbProb   **************************************************************
 *
 */
void
SetPairProbProb (pair_prob_t *pair, int i, int j, int N)
{
	pair[i].pair  = j%N;
	pair[j%N].pair = i;
	pair[i].prob  = (float)AREA(i,j)/norm;
	pair[j%N].prob = (float)AREA(i,j)/norm;
}
/***   SetPairProbProb   ***/

