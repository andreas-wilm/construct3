/******************************************************************************
* 
* tinoco.c - Create tinoco dotplots
*
* Copyright (C) 2002 Gerhard Steger <steger@biophys.uni-duesseldorf.de>
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
 *  CVS $Id: tinoco.c,v 1.7 2007-10-22 10:43:23 steger Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
  
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>
# include <unistd.h>
#include <errno.h>
#ifdef __APPLE__
    #include <sys/wait.h>
#else
    #include <wait.h>
#endif

#include <string.h>
#include <time.h>
#include <math.h>

#include "tinoco.h"

static FILE *PS_dot_common(   char *seq, char *wastlfile,                char *comment);
int          PS_dot_plot_list(char *seq, char *wastlfile, float *matrix, char *comment);

/*--------------------------------------------------------------------------*/
/* from Vienna/lib/utils.c */

void *space(unsigned size) {
  void *pointer;
  
  if ( (pointer = (void *) calloc(1, (size_t) size)) == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"SPACE: requested size: %d\n", size);
      nrerror("SPACE allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      nrerror("SPACE allocation failure -> no memory");
  }
  return  pointer;
}

#ifdef WITH_DMALLOC
#define space(S) calloc(1,(S))
#endif

#undef xrealloc
/* dmalloc.h #define's xrealloc */
void *xrealloc (void *p, unsigned size) {
  if (p == 0)
    return space(size);
  p = (void *) realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"xrealloc: requested size: %d\n", size);
      nrerror("xrealloc allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      nrerror("xrealloc allocation failure -> no memory");  
  }
  return p;
}

char *get_line(FILE *fp) /* reads lines of arbitrary length from fp */
{
  char s[512], *line, *cp;
  
  line = NULL;
  do {
    if (fgets(s, 512, fp)==NULL) break;
    cp = strchr(s, '\n');
    if (cp != NULL) *cp = '\0';
    if (line==NULL)
      line = space(strlen(s)+1);
    else
      line = (char *) xrealloc(line, strlen(s)+strlen(line)+1);
    strcat(line, s);
  } while(cp==NULL);
  
  return line;
}

/*-----------------------------------------------------------------*/

int main(int argc, char *argv[])
{
char *al_sequenz,
     *sequenz,
     *line;
char  fname[LEN_FNAME],
      format[6];
int   i, j, k, l,
      al_seqlen, seqlen,
      s, s_offset,
      helixlen = 2,
      *seqcode,
      diff,
      helix[MAXLENHELIX];
int   istty;

int     doHelixWeight = 0;
float   sum, maxsum;
float   value = 1.0;
float   *erg;

for (i=1; i<argc; i++) {
   if (argv[i][0]=='-') 
      switch( argv[i][1] ) {
         case 'l':  
            sscanf(argv[++i], "%d", &helixlen);
            break;
         case 'w':
            doHelixWeight =1;
            value = -1.0;
            break;
         case 'h':  
         default:
            usage();
      }
}

if (helixlen <= 0 || helixlen > MAXLENHELIX) {
   printf("Minimum length of helix: 0 < l <= %d\n", MAXLENHELIX);
   usage();
   exit(EXIT_SUCCESS);
}

istty = isatty(fileno(stdout));

erg = (float *)space(341*sizeof(float *));
                              /*         Xia          mfold3 */
erg[1*64+2*16+1*4+2] =  3.89; /* A:U A:U -0.93= -3.89 -0.90  */
erg[1*64+2*16+2*4+1] =  4.61; /* A:U U:A -1.10= -4.61 -1.10  */
erg[1*64+2*16+3*4+4] =  8.71; /* A:U G:C -2.08= -8.71 -2.10  */
erg[1*64+2*16+4*4+3] =  9.84; /* A:U C:G -2.35= -9.84 -2.20  */
erg[1*64+2*16+3*4+2] =  2.30; /* A:U G:U -0.55= -2.30 -0.60  */
erg[1*64+2*16+2*4+3] =  5.69; /* A:U U:G -1.36= -5.69 -1.40  */
erg[2*64+1*16+1*4+2] =  5.57; /* U:A A:U -1.33= -5.57 -1.30  */
erg[2*64+1*16+2*4+1] =  3.89; /* U:A U:A -0.93= -3.89 -0.90  */
erg[2*64+1*16+3*4+4] =  8.83; /* U:A G:C -2.11= -8.83 -2.10  */
erg[2*64+1*16+4*4+3] =  9.84; /* U:A C:G -2.35= -9.84 -2.40  */
erg[2*64+1*16+3*4+2] =  4.19; /* U:A G:U -1.00= -4.19 -1.00  */
erg[2*64+1*16+2*4+3] =  5.32; /* U:A U:G -1.27= -5.32 -1.30  */
erg[3*64+4*16+1*4+2] =  9.84; /* G:C A:U -2.35= -9.84 -2.40  */
erg[3*64+4*16+2*4+1] =  9.38; /* G:C U:A -2.24= -9.38 -2.20  */
erg[3*64+4*16+3*4+4] = 13.65; /* G:C G:C -3.26=-13.65 -3.30  */
erg[3*64+4*16+4*4+3] = 14.32; /* G:C C:G -3.42=-14.32 -3.40  */
erg[3*64+4*16+3*4+2] =  6.41; /* G:C G:U -1.53= -6.41 -1.50  */
erg[3*64+4*16+2*4+3] = 10.51; /* G:C U:G -2.51=-10.51 -2.50  */
erg[4*64+3*16+1*4+2] =  8.83; /* C:G A:U -2.11= -8.83 -2.10  */
erg[4*64+3*16+2*4+1] =  8.71; /* C:G U:A -2.08= -8.71 -2.10  */
erg[4*64+3*16+3*4+4] =  9.88; /* C:G G:C -2.36= -9.88 -2.40  */
erg[4*64+3*16+4*4+3] = 13.65; /* C:G C:G -3.26=-13.65 -3.30  */
erg[4*64+3*16+3*4+2] =  5.90; /* C:G G:U -1.41= -5.90 -1.40  */
erg[4*64+3*16+2*4+3] =  8.83; /* C:G U:G -2.11= -8.83 -2.10  */
erg[3*64+2*16+1*4+2] =  5.32; /* G:U A:U -1.27= -5.32 -1.30  */
erg[3*64+2*16+2*4+1] =  5.69; /* G:U U:A -1.36= -5.69 -1.40  */
erg[3*64+2*16+3*4+4] =  8.83; /* G:U G:C -2.11= -8.83 -2.10  */
erg[3*64+2*16+4*4+3] = 10.51; /* G:U C:G -2.51=-10.51 -2.50  */
erg[3*64+2*16+3*4+2] =  2.09; /* G:U G:U -0.50= -2.09 -0.50  */
erg[3*64+2*16+2*4+3] = -5.40; /* G:U U:G +1.29= +5.40  1.30  */
erg[2*64+3*16+1*4+2] =  4.19; /* U:G A:U -1.00= -4.19 -1.00  */
erg[2*64+3*16+2*4+1] =  2.30; /* U:G U:A -0.55= -2.30 -0.60  */
erg[2*64+3*16+3*4+4] =  5.90; /* U:G G:C -1.41= -5.90 -1.40  */
erg[2*64+3*16+4*4+3] =  6.41; /* U:G C:G -1.53= -6.41 -1.50  */
erg[2*64+3*16+3*4+2] = -1.26; /* U:G G:U +0.30= +1.26  0.30  */
erg[2*64+3*16+2*4+3] =  2.09; /* U:G U:G -0.50= -2.09 -0.50  */


sprintf(format, ">%%%ds", LEN_FNAME);
do {
    fname[0]    = '\0';
    if (istty) {
        printf("\n\nInput string (Upper or lower case); @ to quit\n");
        printf("%s%s\n", scale1, scale2);
    }
    if ((line = get_line(stdin))==NULL) break;

    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
        if (*line=='>') (void) sscanf(line, format, fname);
        printf("%s\n", line); fflush(stdout);
        free(line);
        if ((line = get_line(stdin))==NULL) break;
    }

    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;

    al_sequenz = (char *) space(strlen(line)+1);
       sequenz = (char *) space(strlen(line)+1);
    (void) sscanf(line, "%s", al_sequenz);
    free(line);
    al_seqlen = (int) strlen(al_sequenz);
    if (al_seqlen==0) break;

    for (l = 0; l < al_seqlen; l++) al_sequenz[l] = toupper(al_sequenz[l]);
    printf("\n%s\nal_seqlen = %d\n", al_sequenz, al_seqlen); fflush(stdout);

    s_offset = 0;
    for (s=0; s<=al_seqlen; s++) { /* copies nucleotids and '\n' */
        if (al_sequenz[s] == '-') {
            s_offset++;
            continue;
        }
        sequenz[s-s_offset] = al_sequenz[s];
    }

   seqlen = strlen(sequenz);
   printf("\n%s\n   seqlen = %d\n",    sequenz,    seqlen); fflush(stdout);
   
   matrix  = (float *)space((size_t)(sizeof(float)*(seqlen+1)*(seqlen+1)));
   seqcode = (int   *)space((size_t)(sizeof(int  )*(seqlen+1)));

   for (i=0; i<=seqlen-1; i++) {
      switch(sequenz[i]) {
	 case 'A':
	    seqcode[i] = 1;
	    break;
	 case 'U':
	    seqcode[i] = 2;
	    break;
	 case 'T':
	    seqcode[i] = 2;
	    break;
	 case 'G':
	    seqcode[i] = 3;
	    break;
	 case 'C':
	    seqcode[i] = 4;
	    break;
	 case '-':
	    seqcode[i] = 8;
	    break;
	 default : 
	    seqcode[i] = 10;
      }
   }
    
   for (j=0; j<=seqlen-1; j++) {
        for (i=seqlen-1; i>=j+2*helixlen-2; i--) {
            for (k=0; k<=helixlen-1; k++) {
                diff = seqcode[j+k] - seqcode[i-k];
                if (abs(diff) == 1)
                    helix[k] = seqcode[j+k] + seqcode[i-k];
                else
                    goto en;
            }
	        for (k=0; k<=helixlen-1; k++)
                matrix[(j+k)*seqlen+(i-k)] = value;
	    en:
            continue;
        }
   }
    if (doHelixWeight == 1) {
	    maxsum = 0.0;
	    for (j=0; j<=seqlen-1; j++) {
	        for (i=seqlen-1; i>=j+2*helixlen-2; i--) {
	            k   = 0;
	            sum = 0.0;
	            while (matrix[(j+k)*seqlen+(i-k)] == value) {
	                if (matrix[(j+k+1)*seqlen+(i-k-1)] == value) {
	                    sum += erg[seqcode[j+k]  *64 +
                                   seqcode[i-k]  *16 +
                                   seqcode[j+k+1]* 4 +
                                   seqcode[i-k-1]];
	                } 
                    k++;
	            }
	            /* sum = sum/k */;
                sum = sum > 0. ? sum : 0.;  /* neg. values give bad helices */
                sum = sum*sum;
	            maxsum = sum > maxsum ? sum : maxsum;
	            k = 0;
	            while (matrix[(j+k)*seqlen+(i-k)] == value) {
	                matrix[(j+k)*seqlen+(i-k)] = sum;
                    k++;
	            }
	        }
	    }
        printf("MaxSum = %f\n", maxsum);
	    for (j=0; j<=seqlen-1; j++) {
	        for (i=seqlen-1; i>=j; i--) {
	            matrix[(j)*seqlen+i] /= maxsum;
	        }
	    }
    }

   argstr[0] = '\0';
   for(i=0; i<argc; i++) {
      strcat(argstr, argv[i]);
      strcat(argstr, " ");
   } 

   if (fname[0] == '\0')
      strcpy(fname, "dot_ti");
   else
      strcat(fname, "_ti");

   i = PS_dot_plot_list(sequenz, fname, matrix, argstr); /* returns 0 for failure */
   free(matrix);
   free(seqcode);
   free(al_sequenz);
   free(sequenz);
} while (1);
free(erg);
exit(0);
}


/*--------------------------------------------------------------------------*/


void usage(void) {
   nrerror("usage: tinoco [-l #helixlength] [-w] < vie-file");
}


/*-------------------------------------------------------------------------*/

void nrerror(char *message)       /* output message upon error */
{
fprintf(stderr, "\n%s\n", message);
exit(EXIT_SUCCESS);
}

/*--------------------------------------------------------------------------*/

const char *RNAdp_prolog =
"%This file contains the square roots of the base pair probabilities in the form\n"
"% i  j  sqrt(p(i,j)) ubox\n\n"
"%%BeginProlog\n"
"/DPdict 100 dict def\n"
"DPdict begin\n"
"/logscale false def\n"
"/lpmin 1e-05 log def\n\n"
"/box { %size x y box - draws box centered on x,y\n"
"   2 index 0.5 mul sub            % x -= 0.5\n"
"   exch 2 index 0.5 mul sub exch  % y -= 0.5\n"
"   3 -1 roll dup rectfill\n"
"} bind def\n\n"
"/ubox {\n"
"   logscale {\n"
"      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n"
"   } if\n"
"   3 1 roll\n"
"   exch len exch sub 1 add box\n"
"} bind def\n\n"
"/lbox {\n"
"   3 1 roll\n"
"   len exch sub 1 add box\n"
"} bind def\n\n"
"/drawseq {\n"
"% print sequence along all 4 sides\n"
"[ [0.7 -0.3 0 ]\n"
"  [0.7 0.7 len add 0]\n"
"  [-0.3 len sub -0.4 -90]\n"
"  [-0.3 len sub 0.7 len add -90]\n"
"] {\n"
"   gsave\n"
"    aload pop rotate translate\n"
"    0 1 len 1 sub {\n"
"     dup 0 moveto\n"
"     sequence exch 1 getinterval\n"
"     show\n"
"    } for\n"
"   grestore\n"
"  } forall\n"
"} bind def\n\n"
"/drawgrid{\n"
"  0.01 setlinewidth\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"  dup 1 gt {\n"
"     dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"  } { [0.3 0.7] 0.1 setdash } ifelse\n"
"  0 exch len {\n"
"     dup dup\n"
"     0 moveto\n"
"     len lineto \n"
"     dup\n"
"     len exch sub 0 exch moveto\n"
"     len exch len exch sub lineto\n"
"     stroke\n"
"  } for\n"
"  [] 0 setdash\n"
"  0.04 setlinewidth \n"
"  currentdict /cutpoint known {\n"
"    cutpoint 1 sub\n"
"    dup dup -1 moveto len 1 add lineto\n"
"    len exch sub dup\n"
"    -1 exch moveto len 1 add exch lineto\n"
"    stroke\n"
"  } if\n"
"  0.5 neg dup translate\n"
"} bind def\n\n"
"end\n"
"%%EndProlog\n";
/*-----------------------------------------------------------------*/

char *time_stamp(void)
{
  time_t  cal_time;
  
  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}

/*-----------------------------------------------------------------*/

static FILE * PS_dot_common(char *seq, char *fname, char *comment) {
  /* write PS header etc for all dot plot variants */
  char wastlfile[LEN_DP_FNAME];
  FILE *wastl;
  char name[LEN_DP_FNAME], *c;
  int  i, length;

  strcpy(wastlfile, fname);
  strcat(wastlfile, ".ps");
  
  length= strlen(seq);
  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    fprintf(stderr, "can't open %s for dot plot\n", wastlfile);
    return NULL; /* return 0 for failure */
  }
  strncpy(name, wastlfile, LEN_DP_FNAME-1);
  if ((c=strrchr(name, '_'))!=0) *c='\0';

  fprintf(wastl,
	  "%%!PS-Adobe-3.0 EPSF-3.0\n"
	  "%%%%Title: RNA Dot Plot\n"
	  "%%%%Creator: PS_dot.c tinoco\n"
	  "%%%%CreationDate: %s", time_stamp());
  fprintf(wastl, "%%%%BoundingBox: 66 211 518 662\n");
  fprintf(wastl,
	  "%%%%DocumentFonts: Helvetica\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n"
	  "%%Options: %s\n", comment);

  if (comment) fprintf(wastl,"%% %s\n",comment);

  fprintf(wastl,"%s", RNAdp_prolog);

  fprintf(wastl,"DPdict begin\n"
	  "%%delete next line to get rid of title\n"
	  "270 665 moveto /Helvetica findfont 14 scalefont setfont "
	  "(%s) show\n\n", name);

  fprintf(wastl,"/sequence { (\\\n");
  for (i=0; i<strlen(seq); i+=255)
    fprintf(wastl, "%.255s\\\n", seq+i);
  fprintf(wastl,") } def\n");
  fprintf(wastl,"/len { sequence length } bind def\n\n");


  fprintf(wastl,"72 216 translate\n"
	  "72 6 mul len 1 add div dup scale\n");
  fprintf(wastl, "/Helvetica findfont 0.95 scalefont setfont\n\n");

  fprintf(wastl,"drawseq\n"
	    "0.5 dup translate\n"
	    "%% draw diagonal\n"
	    "0.04 setlinewidth\n"
	    "0 len moveto len 0 lineto stroke \n\n"
	    "drawgrid\n");
  return(wastl);
}

int PS_dot_plot_list(char *seq, char *wastlfile, float *matrix, char *comment) {
  FILE *wastl;
  int length;
  int i,j;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, comment);
  if (wastl==NULL) return 0; /* return 0 for failure */

  fprintf(wastl,"%%data starts here\n");
  /* print boxes in upper right half*/
  for (i=1; i<length; i++) {
     for (j=1; j<=length; j++) {
           if (j <= i) continue;
           if (matrix[(i-1)*length +j-1]>0.) {
                fprintf(wastl,"%d %d %1.9f ubox\n", i, j, sqrt(matrix[(i-1)*length +j-1]));
           }
       }
  }

  fprintf(wastl,"showpage\n"
	  "end\n"
	  "%%%%EOF\n");
  fclose(wastl);
  return 1; /* success */
}

