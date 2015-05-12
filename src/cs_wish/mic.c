/******************************************************************************
* 
* mic.c - routines for the mutual information content (single nucleotide dependencies)
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
 *  CVS $Id: mic.c,v 1.20 2007-10-22 10:43:23 steger Exp $    
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

  
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h> /* tolower */
#include <string.h> /* tolower */


#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"
#include "struct.h"

#include "mic.h"



#if 0
    #define DEBUG
#endif


/***   mutinfo_t
 */
typedef struct
{
    int        unbiased; /* bool */
    int        bit;      /* bool */
    float      max_mic;  /* maximum mic using the above factors */
    float    **mic;      /* mutual info content dp matrix */
    int        dim;      /* dimension = aln_len */
    int        use_alifoldscore;  /* bool; instead of MIC use the RNAalifold score*/
} mutinfo_t;

static mutinfo_t *mutinfo;


/* Backend to ComputeMutualInfo
 */
float **
MutualInfo(int **seq_matrix, int n_seq, int aln_len,
           int unbiased, int bit, int pair_entropy_norm);
float **
AliFoldScore(int **seq_mx, int n_seq, int aln_len, int use_stacking);

/* Constructor
 */
static mutinfo_t *
New_MutInfo(int aln_len, int unbiased, int bit, int use_alifoldscore);

/* Destructor
 */
static void
Kill_MutInfo(mutinfo_t *m);

/* Convert alignment to int matrix
 */
int **
AlnSeq2IntMat(char **alnseq, int alnlen, int n_seq);

/* don't change order
   or the loops use false indices ! */
enum {NT_IDX_GAP = 0, NT_IDX_A, NT_IDX_C, NT_IDX_G, NT_IDX_U, NT_IDX_Y, NT_IDX_R, NT_IDX_N};
 
/*
 ***/



/*     cs2     cs3
    gap -> 0 NT_IDX_GAP
    a   -> 1 NT_IDX_A
    c   -> 2 NT_IDX_C
    g   -> 3 NT_IDX_G
    u/t -> 4 NT_IDX_U
    y   -> 5 NT_IDX_Y
    r   -> 6 NT_IDX_R
    n,x -> 7 NT_IDX_N
*/





/***   New_MutInfo   **********************************************************
 *
 */
mutinfo_t *
New_MutInfo(int aln_len, int unbiased, int bit, int use_alifoldscore)
{
    mutinfo_t *info_mat = (mutinfo_t *) Scalloc(1, sizeof(mutinfo_t));

    info_mat->unbiased = unbiased;
    info_mat->bit      = bit;
    
    
    if (use_alifoldscore) {
        info_mat->max_mic = 1.;
    } else {
	    /* bit->  do calculations with log_2
	              (see Schneider, Stormo, Gold & Ehrenfeucht (1986))
	       else-> do calculations with log_e
	              (see Chiu & Kolodziejczak (1991))
	       maximum value is the log of options (=number of nucleotide symbols (=5))
	    */
		if (bit==1) /*base = 2*/
			info_mat->max_mic = log(5.0)/log(2.0); /* log_a(x)= log_b(x) / log_b(a) */
	    else        /*base = e*/
			info_mat->max_mic = log(5.0);
    }

    info_mat->mic      = NULL;
    info_mat->dim      = aln_len;
    
    return info_mat;
}
/***   New_MutInfo   ***/




/***   Kill_MutInfo   *********************************************************
 *
 */
void
Kill_MutInfo(mutinfo_t *m)
{
    if (m!=NULL)
    {
        if (m->mic!=NULL)
            mx_free((void**)m->mic, m->dim, MX_DP);
        free(m);
    }
    m = NULL;
}
/***   Kill_MutInfo   ***/




/***   ComputeMutualInfoContent   *********************************************
 *
 * Compute mutinfo_t *mutinfo
 *
 *
 */
void
ComputeMutualInfoContent(char **alnseq, int alnlen, int n_seq,
                                        int unbiased, int bit, int pair_entropy_norm,
                                        int alifoldscore)
{    
    int **seq_matrix;
    
    seq_matrix = AlnSeq2IntMat(alnseq, alnlen, n_seq);
    
    if (mutinfo!=NULL)
        Kill_MutInfo(mutinfo);
    mutinfo = New_MutInfo(alnlen, unbiased, bit, alifoldscore);

    if (alifoldscore == 0) {
    mutinfo->mic = MutualInfo(seq_matrix, n_seq, alnlen,
                              unbiased, bit , pair_entropy_norm);
    } else {
        int use_stacking;
        if (alifoldscore == 1) {
            use_stacking = 0;
        } else if (alifoldscore == 2) {
            use_stacking = 1;
        } else {
            ERROR_P("Unknown alifold/stacking option %d. Not using stacking\n", alifoldscore);
            use_stacking = 0;
        }
        mutinfo->mic = AliFoldScore(seq_matrix, n_seq, alnlen, use_stacking);
    }

	mx_free((void**)seq_matrix, n_seq+1, MX_DEFAULT);
}
/***   ComputeMutualInfoContent   ***/




/***   AlnSeq2IntMat   ********************************************************
 *
 * Transform nucleotide symbols to integer for easier handling
 * Unit-Offset
 * Caller must free
 */
int **
AlnSeq2IntMat(char **alnseq, int alnlen, int n_seq)
{
  	int  **ret;
    int    seq_idx;
    int    nt_idx;
    char   nt;
    char  *as;


    ret  = (int**) mx_new(n_seq+1, alnlen+1, sizeof(int), sizeof(int*), MX_DEFAULT);
	
    for (seq_idx=0; seq_idx<n_seq; seq_idx++)
    {
        as = alnseq[seq_idx];
        #ifdef DEBUG
            DEBUG_P("transforming seq (%s) int matrix\n", as);
        #endif
        for (nt_idx=0; nt_idx<alnlen; nt_idx++)
        {
            nt = (char)  tolower(as[nt_idx]);
            switch (nt)
            {
				case 'a':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_A;
					break;
				case 'c':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_C;
					break;
				case 'g':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_G;
					break;
				case 'u':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_U;
					break;
				case 'y':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_Y;
					break;
				case 'r':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_R;
					break;
				case 'n':
					ret[seq_idx+1][nt_idx+1] = NT_IDX_N;
					break;
				default:
                    if ( IS_GAP(nt) )
                    {
                        ret[seq_idx+1][nt_idx+1] = NT_IDX_GAP;
                    }
                    else
                    {
                        WARN_P("converting unexpected nt character \"%c\" to n\n", nt);
                        ret[seq_idx+1][nt_idx+1] = NT_IDX_N;
                    }
                    break;
			}
		}
	}
    return ret;
}
/***   AlnSeq2IntMat   ***/




/***    MutualInfo   **********************************************************
 *
 * Returns mic (unit-offset)
 *
 * unbiased: (bool)
 * bit: (bool)
 * pair_entropy_norm: (bool) normalize MIC by dividing with pair-entropy
 *                    (Martin et.al Bioformatics 2005)
 *
 */
float **
MutualInfo(int **seq_mx, int n_seq, int aln_len, int unbiased, int bit, int pair_entropy_norm)
{
	int        i,j,k,m,n,u,v;
    float    **icntr;
    float    **ijcntr;
	float      uprob, vprob, uvprob;  
	float      n16, n4, n1, loge2;
    float    **info_mx;

    float    **info_hxy;
    /* min_mx, max_mx;     mean_hxy;   */

    #ifdef DEBUG
        DEBUG_P("%s\n", "start");
    #endif

    info_mx = (float**)  mx_new(aln_len+1, 0, sizeof(float), sizeof(float*), MX_DP);
    icntr   = (float**)  mx_new(aln_len+1, 8, sizeof(float), sizeof(float*), MX_DEFAULT);   /* icntr[u]: number of X_i=u observed              */
	ijcntr  = (float**)  mx_new(8,         8, sizeof(float), sizeof(float*), MX_DEFAULT);   /* ijcntr[u][v]: number of X_i=u && X_j=v observed */
    if (pair_entropy_norm) {
        info_hxy= (float**)  mx_new(aln_len+1, 0, sizeof(float), sizeof(float*), MX_DP);
    } else {
        info_hxy = NULL;
    }
    
    
	/***   Calculate single nucleotide dependencies
     */
    /* using unbiased probability estimation
       (see Chiu & Kolodziejczak (1991)) */
	if (unbiased==1)
    { 
		n16 = (float)(n_seq + 25.0);
		n4  = (float)(n_seq +  5.0);
		n1  =  1.0;
	}
    /* using maximum likelyhood estimation
       (see Gutell, Power, Hertz, Putz & Stormo (1992)) */
    else
    {
		n16 = (float)(n_seq);
		n4  = (float)(n_seq);
		n1  =  0.0;
	}

    /* do calculations with log_2
      (see Schneider, Stormo, Gold & Ehrenfeucht (1986)) */
	if (bit==1)
		loge2 = log(2.0);
    /* do calculations with log_e
      (see Chiu & Kolodziejczak (1991)) */
    else
		loge2 = 1.0;


    #ifdef DEBUG
        fprintf(stdout, "seq_mx (%s|%s):", __FILE__, __FUNCTION__);
    	for (k=1; k<=n_seq; k++)
        {
            fprintf(stdout, "\n");
            fprintf(stdout, "seq %03d : ", k);
        	for (i=1; i<=aln_len; i++)
            {
                fprintf(stdout, "%d", seq_mx[k][i]);
            }
        }
        fprintf(stdout, "\n");
        fflush(stdout);
    #endif    



	
    /* calculate entropy of columns */
    /* These values will be used for calculation of MIC */
    for (i=1; i<=aln_len; i++)
    { 
        for (u=0; u<=7; u++)
        {
			icntr[i][u]=0.0;
        }

		for (k=1; k<=n_seq; k++)
        {
            /* icntr[u]: number of X_i=u observed */
			icntr[i][seq_mx[k][i]] += 1.0;
            
			switch(seq_mx[k][i])
            {               
				case NT_IDX_GAP:
				case NT_IDX_A:
				case NT_IDX_C:
				case NT_IDX_G:
				case NT_IDX_U:
					break;
				case NT_IDX_R: /* -> r:purine : a || g */
					icntr[i][NT_IDX_A] += 0.5;
					icntr[i][NT_IDX_G] += 0.5;
				case NT_IDX_Y: /* -> y:pyrimidine : u || c */
					icntr[i][NT_IDX_C] += 0.5;
					icntr[i][NT_IDX_U] += 0.5;
				case NT_IDX_N:
					icntr[i][NT_IDX_A] += 0.25;
					icntr[i][NT_IDX_C] += 0.25;
					icntr[i][NT_IDX_G] += 0.25;
					icntr[i][NT_IDX_U] += 0.25;
			}/* switch(seq_mx[k][i]) */

		} /* foreach seq */

        for (u=0; u<=4; u++)
        {
		    icntr[i][u]=(icntr[i][u]+n1)/n4;
        }

	} /* all X_i */







    /* calculate MIC */

    /* foreach possible nt combination */
    
    /* all X_i */
	for (i=1; i<=aln_len; i++)
    { 
        /* all X_j, j<=i */
		for (j=1; j<=i; j++)
        {
        
        
            #ifdef DEBUG
                DEBUG_P("X_i=%d X_j=%d\n", i, j);
            #endif
            /* zero all counters */
			for (m=0; m<=7; m++)
            {
				for (n=0; n<=7; n++)
					ijcntr[m][n]=0.0;
            }
            
			for (k=1; k<=n_seq; k++)
            {
                #ifdef DEBUG
                    DEBUG_P("next seq k=%d\n", k);
				#endif	

                /* ijcntr[u][v]: number of X_i=u && X_j=v observed */
				ijcntr[seq_mx[k][i]][seq_mx[k][j]] += 1.0;

                #ifdef DEBUG
                    DEBUG_P("seq_mx[k=%d][i=%d]=%d\n", k, i, seq_mx[k][i]);
				#endif	
				switch(seq_mx[k][i])
                {               
					case NT_IDX_GAP:
					case NT_IDX_A:
					case NT_IDX_C:
					case NT_IDX_G:
					case NT_IDX_U:
						break;
					case NT_IDX_R:  /* -> r:purine : a || g*/
						if   (seq_mx[k][j]==NT_IDX_GAP
                           || seq_mx[k][j]==NT_IDX_A || seq_mx[k][j]==NT_IDX_C
                           || seq_mx[k][j]==NT_IDX_G || seq_mx[k][j]==NT_IDX_U
                            )
                        {

							ijcntr[NT_IDX_A][seq_mx[k][j]] += 0.5;
							ijcntr[NT_IDX_G][seq_mx[k][j]] += 0.5;
						}
                        else if (seq_mx[k][j] == NT_IDX_R) /* -> r:purine : a || g*/
                        {
							ijcntr[NT_IDX_A][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_A][NT_IDX_G] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_G] += 0.25;
						}
                        else if (seq_mx[k][j] == NT_IDX_Y) /* -> y:pyrimidine : u || c*/
                        {
							ijcntr[NT_IDX_A][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_A][NT_IDX_U] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_U] += 0.25;
						}
                        else if (seq_mx[k][j] == NT_IDX_N)
                        {
							ijcntr[NT_IDX_A][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_A][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_A][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_A][NT_IDX_U] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_U] += 0.125;
						}
						break;
					case NT_IDX_Y: /* -> y:pyrimidine : u || c*/
						if   (seq_mx[k][j]==NT_IDX_GAP
                           || seq_mx[k][j]==NT_IDX_A || seq_mx[k][j]==NT_IDX_C
                           || seq_mx[k][j]==NT_IDX_G || seq_mx[k][j]==NT_IDX_U
                            )
                        {
							ijcntr[NT_IDX_C][seq_mx[k][j]] += 0.5;
							ijcntr[NT_IDX_U][seq_mx[k][j]] += 0.5;
						}
                        else if (seq_mx[k][j] == NT_IDX_R)
                        {
							ijcntr[NT_IDX_C][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_C][NT_IDX_G] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_G] += 0.25;
						}
                        else if (seq_mx[k][j] == NT_IDX_Y)
                        {
							ijcntr[NT_IDX_C][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_C][NT_IDX_U] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_U] += 0.25;
						}
                        else if (seq_mx[k][j] == NT_IDX_N)
                        {
							ijcntr[NT_IDX_C][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_U] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_U] += 0.125;
						}
						break;
					case NT_IDX_N:
						if   (seq_mx[k][j]==NT_IDX_GAP
                           || seq_mx[k][j]==NT_IDX_A || seq_mx[k][j]==NT_IDX_C
                           || seq_mx[k][j]==NT_IDX_G || seq_mx[k][j]==NT_IDX_U
                            )
                        {
							ijcntr[NT_IDX_A][seq_mx[k][j]] += 0.25;
							ijcntr[NT_IDX_C][seq_mx[k][j]] += 0.25;
							ijcntr[NT_IDX_G][seq_mx[k][j]] += 0.25;
							ijcntr[NT_IDX_U][seq_mx[k][j]] += 0.25;
						}
                        else if (seq_mx[k][j] == NT_IDX_R)
                        {
							for (m=1; m<=4; m++)
                            {
								ijcntr[m][NT_IDX_A] += 0.125;
								ijcntr[m][NT_IDX_G] += 0.125;
							}
						}
                        else if (seq_mx[k][j] == NT_IDX_Y)
                        {
							for (m=1; m<=4; m++)
                            {
								ijcntr[m][NT_IDX_C] += 0.125;
								ijcntr[m][NT_IDX_U] += 0.125;
							}
						}
                        else if (seq_mx[k][j] == NT_IDX_N)
                        {
							for (n=1; n<=4; n++)
								for (m=1; m<=4; m++)
									ijcntr[m][n] += 0.0625;
						}
						break;
				}
                /* switch(seq_mx[k][i]) */


                #ifdef DEBUG
                    DEBUG_P("seq_mx[k=%d][j=%d]=%d\n", k, j, seq_mx[k][j]);
				#endif	
				switch(seq_mx[k][j])
                {
					case NT_IDX_GAP:
					case NT_IDX_A:
					case NT_IDX_C:
					case NT_IDX_G:
					case NT_IDX_U:
						break;
					case NT_IDX_R: /* -> r:purine : a || g*/
						if   (seq_mx[k][i]==NT_IDX_GAP
                           || seq_mx[k][i]==NT_IDX_A || seq_mx[k][i]==NT_IDX_C
                           || seq_mx[k][i]==NT_IDX_G || seq_mx[k][i]==NT_IDX_U
                            )
                        {
							ijcntr[seq_mx[k][i]][NT_IDX_A] += 0.5;
							ijcntr[seq_mx[k][i]][NT_IDX_G] += 0.5;
						}
                        else if (seq_mx[k][i] == NT_IDX_R) /* -> r:purine : a || g*/
                        {
							ijcntr[NT_IDX_A][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_A][NT_IDX_G] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_G] += 0.25;
						}
                        else if (seq_mx[k][i] == NT_IDX_Y) /* -> y:pyrimidine : u || c*/
                        {
							ijcntr[NT_IDX_C][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_A] += 0.25;
							ijcntr[NT_IDX_C][NT_IDX_G] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_G] += 0.25;
						}
                        else if (seq_mx[k][i] == NT_IDX_N)
                        {
							ijcntr[NT_IDX_A][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_A] += 0.125;
							ijcntr[NT_IDX_A][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_G] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_G] += 0.125;
						}
						break;
					case NT_IDX_Y:
						if   (seq_mx[k][i]==NT_IDX_GAP
                           || seq_mx[k][i]==NT_IDX_A || seq_mx[k][i]==NT_IDX_C
                           || seq_mx[k][i]==NT_IDX_G || seq_mx[k][i]==NT_IDX_U
                            )
                        {
							ijcntr[NT_IDX_C][seq_mx[k][i]] += 0.5;
							ijcntr[NT_IDX_U][seq_mx[k][i]] += 0.5;
						}
                        else if (seq_mx[k][i] == NT_IDX_R)
                        {
							ijcntr[NT_IDX_A][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_A][NT_IDX_U] += 0.25;
							ijcntr[NT_IDX_G][NT_IDX_U] += 0.25;
						}
                        else if (seq_mx[k][i] == NT_IDX_Y)
                        {
							ijcntr[NT_IDX_C][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_C] += 0.25;
							ijcntr[NT_IDX_C][NT_IDX_U] += 0.25;
							ijcntr[NT_IDX_U][NT_IDX_U] += 0.25;
						}
                        else if (seq_mx[k][i] == NT_IDX_N)
                        {
							ijcntr[NT_IDX_A][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_C] += 0.125;
							ijcntr[NT_IDX_A][NT_IDX_U] += 0.125;
							ijcntr[NT_IDX_C][NT_IDX_U] += 0.125;
							ijcntr[NT_IDX_G][NT_IDX_U] += 0.125;
							ijcntr[NT_IDX_U][NT_IDX_U] += 0.125;
						}
						break;
					case NT_IDX_N:
						if   (seq_mx[k][i]==NT_IDX_GAP
                           || seq_mx[k][i]==NT_IDX_A || seq_mx[k][i]==NT_IDX_C
                           || seq_mx[k][i]==NT_IDX_G || seq_mx[k][i]==NT_IDX_U
                            )
                        {
							ijcntr[seq_mx[k][i]][NT_IDX_A] += 0.25;
							ijcntr[seq_mx[k][i]][NT_IDX_C] += 0.25;
							ijcntr[seq_mx[k][i]][NT_IDX_G] += 0.25;
							ijcntr[seq_mx[k][i]][NT_IDX_U] += 0.25;
						}
                        else if (seq_mx[k][i] == NT_IDX_R)
                        {
							for (m=1; m<=4; m++)
                            {
								ijcntr[NT_IDX_A][m] += 0.125;
								ijcntr[NT_IDX_G][m] += 0.125;
							}
						}
                        else if (seq_mx[k][i] == NT_IDX_Y)
                        {
							for (m=1; m<=4; m++)
                            {
								ijcntr[NT_IDX_C][m] += 0.125;
								ijcntr[NT_IDX_U][m] += 0.125;
							}
						}
                        else if (seq_mx[k][i] == NT_IDX_N)
                        {
							for (n=1; n<=4; n++)
								for (m=1; m<=4; m++)
									ijcntr[n][m] += 0.0625;
                        }
						break;
				}
                /* switch(seq_mx[k][j]) */
                
			} /* foreach seq */


            #ifdef DEBUG
    			DEBUG_P("i=%d, j=%d :", i, j);
                for (m=0; m<=4; m++)
                {
    			    DEBUG_P("ic[%d]=%f, jc[%d]=%f", m, icntr[i][m], m, icntr[j][m]);
    				fprintf(stdout, "ic[%d]=%f, ",  m, icntr[i][m]);
    				fprintf(stdout, "jc[%d]=%f\n",  m, icntr[j][m]);
    				for (n=0; n<=4; n++)
    					fprintf(stdout, "ijc[%d][%d]=%f ,  ", m, n, ijcntr[m][n]);
    				fprintf(stdout, "\n");
    			}   
    			fprintf(stdout, "\n\n");
    			fflush(stdout);
            #endif			

            if (i==j)
                continue;

			info_mx[i][j] = 0.0;
            if (pair_entropy_norm) {
                info_hxy[i][j] = 0.0;
            }
			if (icntr[i][0]<0.2 && icntr[j][0]<0.2) {
	            for (u=0; u<=4; u++)
	            {
			uprob =   icntr[i][u];
	
        		for (v=0; v<=4; v++)
	                {
	                    
	                    vprob =   icntr[j][v];
	                    uvprob= (ijcntr[u][v]+n1)/n16;
	
	                    if (uvprob>0.0)
	                    {
	                        if (pair_entropy_norm) {
	                            info_hxy[i][j] -= uvprob * log(uvprob);
	                        }
	                        if (uprob>0.0 && vprob>0.0)
							    info_mx[i][j] += uvprob * log(uvprob/(uprob * vprob));
	                    }
					}
				}
	            if (pair_entropy_norm) {
	                if (info_hxy[i][j]>0.0) {
	                    info_mx[i][j] = (log(5)/loge2)*info_mx[i][j]/info_hxy[i][j];
	                } else {
	                    info_mx[i][j] = 0.0;
	                }
	            } else {
	                info_mx[i][j] /= loge2;
	            }
            }

#ifdef DEBUG
            DEBUG_P("info_mx(i=%d,j=%d)=%f\n", i-1, j-1, info_mx[i][j]);
#endif

		} /* all X_j, j<=i */
	} /* all X_i */


/*
    min_mx = 1000.;
    max_mx =    0.;
    for (i=1; i<=aln_len; i++) {
		for (j=1; j<=i; j++) {
            min_mx = (min_mx > info_mx[i][j]) ? info_mx[i][j] : min_mx;
            max_mx = (max_mx < info_mx[i][j]) ? info_mx[i][j] : max_mx;
        }
    }
    printf("min_mx = %5.2f\n", min_mx);
    printf("max_mx = %5.2f\n", max_mx); fflush(stdout);
*/
	mx_free((void**)icntr,    aln_len+1, MX_DEFAULT);
	mx_free((void**)ijcntr,   8,         MX_DEFAULT);
    if (pair_entropy_norm) {
        mx_free((void**)info_hxy, aln_len+1, MX_DP);
    }
    
	return info_mx;
}
/***    MutualInfo   ***/


/***   AliFoldScore   ********************************************************
 *
 * Returns AliFoldScore (compare with routine mic)
 * code adapted from ViennaRNA-1.6.1/lib/alifold.c: make_pscores()
 * for equations see 
 *      Hofacker, IL, Fekete, M & Stadler, PF (2002) JMB 319, 1059-1066.
 *      Secondary structure prediction for aligned RNA sequences
 *
 * FIXME: add short doc. regarding stacking (Lindgreen 2006)
 *
 */
float **
AliFoldScore(int **seq_mx, int n_seq, int aln_len, int use_stacking)
{
    #define TURN MIN_HP_SIZE-1
    #define NONE 0. /* score for forbidden pairs */
    const double cv_fact = 1.;
    const double nc_fact = 1.;
    float norm = 3.*n_seq/4.;
    float dummy;

    int i, j, k, l, s;
    int score;
    static int pair[8][8]=
                /* -  A  C  G  U  Y  R  N */
                {{ 0, 0, 0, 0, 0, 0, 0, 0} /* - */,
                 { 0, 0, 0, 0, 5, 0, 0, 0} /* A */,
                 { 0, 0, 0, 1, 0, 0, 0, 0} /* C */,
                 { 0, 0, 2, 0, 3, 0, 0, 0} /* G */,
                 { 0, 6, 0, 4, 0, 0, 0, 0} /* U */,
                 { 0, 0, 0, 0, 0, 0, 0, 0} /* Y */,
                 { 0, 0, 0, 0, 0, 0, 0, 0} /* R */,
                 { 0, 0, 0, 0, 0, 0, 0, 0} /* N */};
    /* A:U = 6; 
       C:G = 2; 
       G:C = 1; G:U = 4;
       U:A = 5; U:G = 3 */
    /* hamming distance between pairs */
    /* eq. 3 */
    int dm[7][7]=
                /*    CG GC GU UG AU UA*/
                 {{0, 0, 0, 0, 0, 0, 0}, 
                  {0, 0, 2, 2, 1, 2, 2} /* CG */,
                  {0, 2, 0, 1, 2, 2, 2} /* GC */,
                  {0, 2, 1, 0, 2, 1, 2} /* GU */,
                  {0, 1, 2, 2, 0, 2, 1} /* UG */,
                  {0, 2, 2, 1, 2, 0, 2} /* AU */,
                  {0, 2, 2, 2, 1, 2, 0} /* UA */};

    float **pscore = (float**) mx_new(aln_len+1, 0, sizeof(float), sizeof(float*), MX_DP); /* precomputed array of pair types */
    float **bscore;

    if (use_stacking) {
        bscore = (float**) mx_new(aln_len+1, 0, sizeof(float), sizeof(float*), MX_DP);
    } else {
        bscore = NULL;
    }

    for (i=1; i<=aln_len; i++) {
        for (j=1; j<i; j++) {
            /* calculate number of base pair types at each dotplot position */
            int pfreq[8]={0,0,0,0,0,0,0,0};
            for (s=1; s<=n_seq; s++) {
                int type;
                if (seq_mx[s][i]==0 && seq_mx[s][j]==0) 
                    type = 7; /* gap-gap  */
                else 
                    type = pair[seq_mx[s][i]][seq_mx[s][j]];

                pfreq[type]++;
            }
            /* skip calculation if only non-WC */
            if (pfreq[0]*2>=n_seq) { 
                pscore[i][j] = NONE; 
                continue;
            }
            /* eq. 4.2 */
            for (k=1,score=0; k<=6; k++) { /* ignore pairtype 7 (gap-gap) */
                for (l=k+1; l<=6; l++) {
                    /* scores for replacements between pairtypes    */
                    /* consistent or compensatory mutations score 1 or 2  */
                    score += pfreq[k]*pfreq[l]*dm[k][l];
                }
            }
            /* counter examples score -1, gap-gap scores -0.25   */
            /* ??? pscore[i][j] = cv_fact * ((1.*score)/n_seq - nc_fact*(pfreq[0] + pfreq[7]*0.25)); */
            /* max pscore = (n_seq/4)^2 * 6 * max(dm) = n_seq^2 * 3/4
               Thus scale the SCORE by NORM=3.*N_SEQ/4., which results in 
               SCORE <= N_SEQ; than subtract number of nonWC and gap pairs.
               Why are gap pairs scaled by .25 ???
            */
            dummy = cv_fact * ((1.*score)/norm - nc_fact*(pfreq[0] + pfreq[7]*0.25));
            dummy = dummy/n_seq;
            /* fprintf(stdout, "pscore[%2d][%2d] = %f\n", i,j,dummy); */
            pscore[i][j] = (dummy>0.) ? dummy : 0.;
        }
    }

    if (use_stacking) {
        for (i=aln_len-1; i>1; i--) {
            for (j=2; j<i; j++) {
                bscore[i][j] = (2.*pscore[i][j]+pscore[i-1][j+1]+pscore[i+1][j-1])/4.;
            }
        }
        for (i=aln_len; i>1; i--) {
            bscore[i][1] = (2.*pscore[i][1]+pscore[i-1][2])/3.;
        }
        for (j=2; j<aln_len; j++) {
            bscore[aln_len][j] = (2.*pscore[aln_len][j]+pscore[aln_len-1][j+1])/3.;
        }
        /*
        aln_len=9
        i
        1 2 3 4 5 6 7 8 9 
          2 2 2 2 2 2 2 2 1 j
          1 1 1 1 1 1 1 3 2
            1 1 1 1 1 1 3 3
              1 1 1 1 1 3 4
                  1 1 1 3 5
                    1 1 3 6
                      1 3 7
                        3 8
                          9
        */                      

        mx_free((void**)pscore, aln_len+1, MX_DP);
        return bscore;
        
    } else {
        
        return pscore;
    }
}
/***   AliFoldScore   ***/


/***   GetConsMIC   ***********************************************************
 *
 * pendant to GetConsBpProb
 *
 * mic_thresh: range 0-1
 * returned mic is normalized to 1 (= mic') by dividing with max_mic !
 *

 *
 * mic must be greater then mic_thresh, otherwise 0.0 will be returned
 *
 * nt_i, nt_j = zero offset nt indices / columns
 * 
 */
float
GetConsMIC(pair_t pair, float mic_thresh)
{
    float cmic;

    cmic  = GetMIC(pair) / mutinfo->max_mic;
    if (cmic <= mic_thresh)
        cmic = 0.0;
    return cmic;
}
/***   GetConsMIC   ***/



/***   GetMIC   ***********************************************************
 *
 * return mic of two columns
 *
 * nt_i, nt_j = zero offset nt indices / columns
 * 
 */
float
GetMIC(pair_t pair)
{   
    /* internal computation was unit-offset */
    return  mutinfo->mic[pair.nti+1][pair.ntj+1];
}
/***   GetMIC   ***/



/***   MIC_IsComputed   *******************************************************
 *
 */
int
MIC_IsComputed(void)
{
    if (mutinfo==NULL)
        return 0;
    else
        return 1;
}
/***   MIC_IsComputed   ***/



/***   GetMaxMIC   ************************************************************
 *
 */
float
GetMaxMIC(void)
{
    if (! MIC_IsComputed())
    {
        ERROR_P("%s\n", "Compute MIC before calling me");
        return 0.0;
    }

    return mutinfo->max_mic;
}
/***   GetMaxMIC   ***/



/***   PrintCsDpm   *********************************************************
 *
 *
 *
 */
void
PrintMic(char *csseq, FILE *stream)
{
    int i=0, j=0;
    int seqlen  = 0;
    pair_t pair;

    fprintf(stream, "FIXME(%s:%s): UNTESTED\n", __FILE__, __FUNCTION__);
    
    if (! MIC_IsComputed())
    {
        fprintf(stream, "Compute Mutual Info Content before calling me\n");
        return;
    }
    fprintf(stream, "Mutual Info Content Matrix\n");
    fprintf(stream, "unbiased=%d bit=%d max_mic=%f\n\n",
                                      mutinfo->unbiased,
                                      mutinfo->bit,
                                      mutinfo->max_mic);

    
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
            SetPair(i, j, &pair);
            if (i>j)
                fprintf(stream, "  ---  ");
            else if (j==i)
                fprintf(stream, "   \\   ");
            else
                fprintf(stream, " %0.3f ", GetMIC(pair));
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
