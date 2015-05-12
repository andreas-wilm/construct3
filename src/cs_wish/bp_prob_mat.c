/******************************************************************************
 *
 * bp_prob_mat.c - routines for reading basepair probability matrices
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
 *  CVS $Id: bp_prob_mat.c,v 1.34 2007-05-10 16:30:44 steger Exp $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h> /* read */


#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "bp_prob_mat.h"




/***   private
 */

#if 0
#define DEBUG
#endif


#ifdef HAVE_LIBZ
#include <zlib.h>
#define FILE_TYPE gzFile
#define FILE_OPEN(x,y) gzopen(x, y)
#define FILE_CLOSE(x) gzclose(x)
#define FILE_GETS(line, len, fid) gzgets(fid, line, len)
#define FILE_PRINTF(fid, fmt, args...) gzprintf(fid, fmt, ## args);
#else
#define FILE_TYPE FILE
#define FILE_OPEN(x,y) fopen(x, y)
#define FILE_CLOSE(x) fclose(x)
#define FILE_GETS(line, len, fid) fgets(line, len, fid)
#define FILE_PRINTF(fid, fmt, args...) fprintf(fid, fmt, ## args);
#endif


static float *
TranformRnaFoldMatTo53(float *rnafold_mat, char *degapped_seq, float prob_cutoff);
static float *
ReadRnaFoldMat(char *path2mat, char *degapped_seq);
static int
Idx_2D_to_1D_inclDiag(int idx2, int idx1, int mat_len);
static int
RnaFoldPs_CheckHeader(FILE_TYPE *fd);
static char *
RnaFoldPs_GetSequence(FILE_TYPE *fd);
static bpair_t ***
RnaFoldPs_GetDpMatrix(FILE_TYPE *fd, int seqlen, float cutoff);
static int
FnameIsPs(char *fname);
/*
***/




/***   SetupBpProbMx   **********************************************************
 *
 * Frontend to most procedures defined here
 * Invokes reading of a basepair matrix (path2mat)
 * which may be either a binary construct version
 * or an rnafold postscript version
 *
 * and handles necessary transformations
 * returned bpair_t doesn't contain the unneccessary diagonal (i=j)
 * some intermediates do !
 * some bpair_t's maybe NULL for memory saving, if prob=0.0
 * seq must not contain any gaps !
 *
 * returns NULL on error
 */
bpair_t ***
SetupBpProbMx(char *path2mat, char *seq, float prob_cutoff)
{
    float *rnafold_mat;
    float *intermed_mat;
    int intermed_idx;
    int seqlen;
    bpair_t ***bpmx;
    int i, j;
    Stopwatch_t *watch;

    watch = TimerStart();


    DEBUG_P("Reading Basepair-Matrix-File \"%s\"\n", path2mat);
    VERBOSE_P("Reading Basepair-Matrix-File \"%s\"\n", path2mat);

    
    /* read from postscript */
    if (FnameIsPs(path2mat)) {
        FILE_TYPE *fd;
        char *seqbuf;

        DEBUG_P("%s", "Reading Postscript Matrix\n");

        if ((fd =  FILE_OPEN(path2mat, "r"))==NULL) {
            ERROR_P("file open failed for %s\n", path2mat);
            return NULL;
        }
        if ( ! RnaFoldPs_CheckHeader(fd)) {
            ERROR_P("CheckHeader failed for %s: incompatible RNAfold version?\n", path2mat);
            FILE_CLOSE(fd);
            return NULL;
        }

        if ((seqbuf=RnaFoldPs_GetSequence(fd))==NULL) {
            ERROR_P("GetSequence failed for %s\n", path2mat);
            return NULL;
        } else {
            if ( ! STR_NC_EQ(seqbuf, seq)) {
                WARN_P("sequences differ for %s -> \"%s\" != \"%s\"\n",
                       path2mat, seqbuf, seq);
            }
            free(seqbuf);
        }
        
        seqlen = strlen(seq);
        if ((bpmx=RnaFoldPs_GetDpMatrix(fd, seqlen, prob_cutoff))==NULL) {
            FILE_CLOSE(fd);
            return NULL;
        }

        FILE_CLOSE(fd);


        
    /* read from binary */
    } else {
        DEBUG_P("%s", "Reading Binary Matrix\n");
        rnafold_mat = ReadRnaFoldMat(path2mat, seq);
        
        DEBUG_P("%s", "Transforming Matrix\n");
        intermed_mat = TranformRnaFoldMatTo53(rnafold_mat, seq, prob_cutoff);
        
        
        /* clean returned bpair_t mx
         */
        DEBUG_P("%s", "Cleaning Matrix\n");
        seqlen = strlen(seq);
        bpmx   = (bpair_t ***) mx_new (seqlen, seqlen, sizeof(bpair_t*), sizeof(bpair_t**), MX_DP);

        /* create clean from intermediate
         */
        for (i=0; i<seqlen; i++) {
            for (j=i+1; j<seqlen; j++) {
                intermed_idx  = Idx_2D_to_1D_inclDiag(i, j, seqlen);
                /* isn't this nonsense ?
                bpmx[j][i] = (bpair_t*) Scalloc(1, sizeof(bpair_t));
                if (intermed_mat[intermed_idx] > 0.0)
                    bpmx[j][i]->prob = intermed_mat[intermed_idx];
                else
                    bpmx[j][i] = NULL;
                */
                if (intermed_mat[intermed_idx]>0.0) {
                    bpmx[j][i] = (bpair_t*) Scalloc(1, sizeof(bpair_t));
                    bpmx[j][i]->prob = intermed_mat[intermed_idx];
                }
#ifdef DEBUG
                DEBUG_P("bpmx[%d][%d] = %f\n", j, i, bpmx[j][i]->prob);
#endif
            }
        }

        free(rnafold_mat);
        free(intermed_mat);
    }
    
    TimerStop(__FUNCTION__, watch);

    /* return ret_basepair; */
    return bpmx;
}
/***   SetupBpProbMx   ***/



/***   FreeBpProbMx   *********************************************************
 *
 * seqlen must be that of unaligned seq
 *
 */
void
FreeBpProbMx(bpair_t ***bpmx, int seqlen)
{
    mx_free((void**)bpmx, seqlen, MX_DP);
}
/***   FreeBpProbMx   ***/



/***   ReadRnaFoldMat   *********************************************************
 *
 * Actually reads a rnafold basepair matrix
 * and returns it as a 1D untransformed array
 * Caller must free
 *
 */
float *
ReadRnaFoldMat(char *path2mat, char *degapped_orig_seq)
{
    
    /* read matrices in 16k-buffers */
    const int buflen = 16384;

    FILE_TYPE *FP_DP_DAT = NULL;
	char     buf[buflen];
	int      length=0, rest = 0;
	int      dummy = 0;
	char     err[buflen];


	char     *seq = NULL;
	int      seqlen = 0;
	char     cmd_str[200];

    /* RNAfold's basepair probability matrix
	   in 1D from 3 to 5 , 3 to 5 !!! (incl. diag.)
    */
	float    *rnafold_bp_prob_mat     = NULL;
    float    *rnafold_ptr             = NULL;
    int      rnafold_bp_prob_mat_dim  = 0;


    /*****   open matrix file
	 */
	FP_DP_DAT = FILE_OPEN(path2mat, "rb");
    if (FP_DP_DAT == NULL)
        {
            sprintf(err, "Can't open prob matrix file \"%s\"!\n", path2mat);
            Die(err);
        }


	/*****   read sequence-length from file
	 */
#ifdef HAVE_LIBZ
    if (gzread(FP_DP_DAT, &seqlen, sizeof(int)) != sizeof(int)) {
        sprintf(err, "Can't read sequence length from prob matrix file\"%s\"!\n", path2mat);
        Die(err);
    }
#else
    if (read(fileno(FP_DP_DAT), &seqlen, sizeof(int)) != sizeof(int)) {
        sprintf(err, "Can't read sequence length from prob matrix file\"%s\"!\n", path2mat);
        Die(err);
    }
#endif

#ifdef DEBUG
    DEBUG_P("sequence length = %d\n", seqlen);
#endif


	/*****   read matrix
	 *
	 */
    rnafold_bp_prob_mat_dim = (seqlen+1)*(seqlen+2)/2;
	length  = sizeof(float) * rnafold_bp_prob_mat_dim;
	rnafold_bp_prob_mat = (float *) Scalloc(rnafold_bp_prob_mat_dim, sizeof(float));
	rnafold_ptr = rnafold_bp_prob_mat;


	while (1) {
        rest = MIN(sizeof(buf), length);

#ifdef HAVE_LIBZ
        if (gzread(FP_DP_DAT, rnafold_ptr, (unsigned) rest) != rest)
            Die("Can't read prob matrix!");
#else
        if (read(fileno(FP_DP_DAT), rnafold_ptr, (unsigned) rest) != rest)
            Die("Can't read prob matrix!");
#endif

        length -= rest;
        if (length == 0)
            break;
        if (length <  0)
            Die("Can't find end of prob matrix!");
        
        rnafold_ptr += rest/sizeof(float);
    }

#ifdef DEBUG
    DEBUG_P("%s\n", "matrix successfully read");
#endif


	/*****   read seq
	 *
	 */
	seq  = (char *) Scalloc(1, (size_t)(sizeof(char)*(seqlen+1)));
#ifdef HAVE_LIBZ
    dummy = gzread(FP_DP_DAT, seq, sizeof(char) * (seqlen+1));
    if ( dummy != sizeof(char) * (seqlen+1) )
        Die("ERROR: can't read sequence from prob matrix file!");
#else
    dummy = read(fileno(FP_DP_DAT), seq, sizeof(char) * (seqlen+1));
    if ( dummy != sizeof(char) * (seqlen+1) )
        Die("ERROR: can't read sequence from prob matrix file!");
#endif
    StrTolower(seq);

#ifdef DEBUG
    DEBUG_P("read sequence = %s\n", seq);
#endif



    /* check sequence
     *
     */
    if ( ! STR_EQ(degapped_orig_seq, seq) ) {
            WARN_P("sequences differ for %s -> \"%s\" != \"%s\"\n",
                   path2mat, seq, degapped_orig_seq);
    }


	/*****   read cmd_str
	 *
	 */
#ifdef HAVE_LIBZ
    dummy = gzread(FP_DP_DAT, cmd_str, sizeof(char) * 200);
    if ( dummy != sizeof(char) * 200)
        Die("Can't read command from prob matrix file!");
#else
    dummy = read(fileno(FP_DP_DAT), cmd_str, sizeof(char) * 200);
    if ( dummy != sizeof(char) * 200)
        Die("Can't read command from prob matrix file!");
#endif

	/* FIXME: check cmd_str */


    FILE_CLOSE(FP_DP_DAT);

    free(seq);

    return rnafold_bp_prob_mat;
}
/*** ReadRnaFoldMat   ***/




/***   TranformRnaFoldMatTo53   ***********************************************
 *
 * Transforms the standard RNAfold Matrix which is 3'-5' to 5'-3'
 * and cuts all probs below prob_cutoff
 *
 * Returned mat includes diagonal !
 *
 * Caller must free
 *
 */
float *
TranformRnaFoldMatTo53(float *rnafold_mat, char *seq, float prob_cutoff)
{
    int        bp_prob_mat_dim;
    float     *bp_prob_mat_incldiag;
    int        mirror;
    int        vienna_idx;


    int        seqlen = strlen(seq);


    /* intermediate incl diag */
    bp_prob_mat_dim       = (seqlen+1)*seqlen / 2;
	bp_prob_mat_incldiag  = (float *) Scalloc( bp_prob_mat_dim, sizeof(float));


    /* create intermediate from vienna
     */
    mirror=0;
	for (vienna_idx = (seqlen+1)*seqlen/2; vienna_idx>=1; vienna_idx--) {
        if (rnafold_mat[vienna_idx] >= prob_cutoff)
            bp_prob_mat_incldiag[mirror] =  rnafold_mat[vienna_idx];
        else
            bp_prob_mat_incldiag[mirror] =  0.0;
        mirror++;
    }

    return bp_prob_mat_incldiag;
}
/***   TranformRnaFoldMatTo53   ***/





/***   Idx_2D_to_1D_inclDiag   ************************************************
 *
 *  Converts given 2D indices (ZERO-OFFSET)
 *  for a half filled matrix (i,j 0<=i<=j) to
 *  the corresponding 1d indices
 *  mat_len length of array (seq_len)
 *  Assumes Matrix contains the unneccessary diagonal (i=j)
 */
int
Idx_2D_to_1D_inclDiag (int idx1, int idx2, int mat_len)
{
    int i=0;
    int one_dim_index=0;

    if (idx2 < idx1) {
        char buf[1024];
        sprintf(buf, "Invalid Args for %s: idx2 (=%d) < idx1 (=%d)\
                               is invalid (1D!)\n", __FUNCTION__, idx1, idx2);
        Die(buf);
    }

    /* count row-offset with low index:
	 */
	while (i < idx1) {
        one_dim_index += (mat_len-1-i);
        i++;
    }
	/* now do the high index col offset
	 */
	one_dim_index += (idx2);

    return one_dim_index;
}
/***   Idx_2D_to_1D_inclDiag   ***/



/***   RnaFoldPs_CheckHeader
 *
 *
 */
int RnaFoldPs_CheckHeader(FILE_TYPE *fd)
{
    const int max_line_length = 1024;
    char  line[max_line_length];
    
    /*const char header_ident[] = "%%Creator: PS_dot.c, ViennaRNA Package";*/
    /*const char header_ident[] = "%%Title: RNA DotPlot";*/
    const char header_ident[] = "%%Creator: PS_dot.c";
    while (FILE_GETS(line, max_line_length, fd)!=NULL) {
        if (strstr(line, header_ident))
            return 1;
    }
    return 0;
}
/***  RnaFoldPs_CheckHeader  ***/



/***   RnaFoldPs_GetSequence
 *
 * read sequence string from dp file and return
 * user must free
 *
 */
char *
RnaFoldPs_GetSequence(FILE_TYPE *fd)
{
    const int max_line_length = 1024;
    char *seq;
    char  line[max_line_length];
    const char seq_start[] = "/sequence ";
    const char seq_end[] = " def";
    char *cpstr;
    int inseq;

    seq=NULL;
    inseq=0;
    while (FILE_GETS(line, max_line_length, fd)!=NULL) {
        ChompStr(line);
        if (strstr(line, seq_start)!=NULL) {
            inseq=1;
            continue;
        }
        if (inseq) {
            if (strstr(line, seq_end)!=NULL) {
                StrTolower(seq);
                break;
            }
            /* cut off trailing slash */
            #ifndef __APPLE__
              cpstr = strndup(line, strlen(line)-1);
            #else
              cpstr = strdup(line);
            #endif
            
            if (seq!=NULL) {
                seq = Srealloc(seq, strlen(seq) + 1 + strlen(cpstr));
                seq = strcat(seq, cpstr);
            } else {
                seq = strdup(cpstr);
            }
            free(cpstr);
        }
    }
    return seq;
}
/***   RnaFoldPs_GetSequence   ***/




/***   FnameIsPs
 *
 */
int FnameIsPs(char *fname)
{
    if (strstr(fname, "_dp.ps")   !=NULL ||
        strstr(fname, "_dp.ps.Z") !=NULL ||
        strstr(fname, "_dp.ps.gz")!=NULL ||
        strstr(fname, "_ti.ps")   !=NULL ||
        strstr(fname, "_ti.ps.gz")!=NULL)
        return 1;
    else
        return 0;
}
/***   FnameIsPs   ***/



/***   RnaFoldPs_GetDpMatrix
 *
 */
bpair_t ***
RnaFoldPs_GetDpMatrix(FILE_TYPE *fd, int seqlen, float cutoff)
{
    const int max_line_length = 1024;
    char box[5];
    int i, j;
    float prob;
    bpair_t ***bpmx;
    char line[max_line_length];
     
    bpmx = (bpair_t ***) mx_new (seqlen, seqlen, sizeof(bpair_t*), sizeof(bpair_t**), MX_DP);

    while (FILE_GETS(line, max_line_length, fd) != NULL) {
        ChompStr(line);
                
        if (sscanf(line, "%d %d %f %4s", &i, &j, &prob, box) != 4)
            continue;
        if (! STR_EQ(box, "ubox"))
            continue;
        if (i>=j)
            return NULL;
        
        prob=pow(prob,2);

        /* isn't this nonsense?
        bpmx[j-1][i-1] = (bpair_t*) Scalloc(1, sizeof(bpair_t));
        if (prob>cutoff) {
            bpmx[j-1][i-1]->prob = prob;
        else
            bpmx[j-1][i-1] = NULL;
        */
        if (prob>=cutoff) {
            bpmx[j-1][i-1] = (bpair_t*) Scalloc(1, sizeof(bpair_t));
            bpmx[j-1][i-1]->prob = prob;
#ifdef DEBUG
            DEBUG_P("bpmx[%d][%d] = %f\n", j, i, bpmx[j][i]->prob);
#endif
        }
    }

    return bpmx;
}
/*** RnaFoldPs_GetDpMatrix   ***/


/***   PrintTdProbMat
 *
 * Prints a Basepair-TdProbMat to console
 *
 */
void
PrintTdProbMat(char *nalseq, bpair_t ***bpmx, FILE *stream)
{
    int seqlen;
    int i=0, j=0;


    seqlen = strlen(nalseq);

    /* print sequence nt
     */
    fprintf(stream, "   ");
    for (i=0; i<seqlen ; i++)
        fprintf(stream, "  %4c  ", nalseq[i]);
    fprintf(stream, "\n");



    for (i=0; i<seqlen ; i++) {
        fprintf(stream, "%c  ", nalseq[i]);
        for (j=0; j<seqlen ; j++) {
            if (i>j)
                fprintf(stream, "  ----  ");
            else if (j==i)
                fprintf(stream, "  ++++  ");
            else
                if (bpmx[j][i]!=NULL)
                    fprintf(stream, " %0.4f ", bpmx[j][i]->prob);
                else
                    fprintf(stream, " %0.4f ", 0.0);
        }
        fprintf(stream, "\n");
    }
    fflush(stream);
}
/***   PrintTdProbMat   ***/

