/******************************************************************************
* 
* CS_DrawStruct.c - 
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
 *  CVS $Id: CS_DrawStruct.c,v 1.7 2007-05-11 13:10:56 steger Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "drawstruct.h"
#include "CS_DrawStruct.h"
#include "generic_utils.h"


/*#define MAXLENGTH 2000*/          /* bis zu 2200 nt laufen auf iris (7/93) */
#define PI           3.141592654  
#define PIHALF    PI/2.

#define BRANCH     100	/* should be identical in size to MAX_NumBif in CS_OptStruct.c */

#define INIT_ANGLE     0.     /* initial bending angle */
#define INIT_X      1000.     /* coordinate of first digit */
#define INIT_Y      1000.     /* see above */
#define RADIUS       150.


struct bond {         /* bonding list */
   int i;
   int j;
};

struct bond *base_pair;

int	straight;

float angle[MAXLENGTH+5], x[MAXLENGTH+5], y[MAXLENGTH+5];

int this_loop_size[MAXLENGTH/6], this_stack_size[MAXLENGTH/6], lp, stk;
/* name:stack_size interfers with libc? */


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

static int
pairs(int i)
{

    struct bond *n;

    for (n = base_pair+1; n <= base_pair+(*base_pair).i; n++)
    {
        if (i == n->i)
            return(n->j);
    }
    return(0);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/

static void loop(int i, int j, int cc)

             /* i, j are the positions AFTER the last pair of a stack; i.e
		i-1 and j+1 are paired. */
	     /* cc corrects count of the first loop */
{
    int    count = 2;   /* counts the VERTICES of a loop polygon; that's
			   NOT necessarily the number of unpaired bases!
			   Upon entry the loop has already 2 vertices, namely
			   the pair i-1/j+1.  */

    int    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

    int    i_old, partner, k, l, start_k, start_l, fill, ladder;
    int    begin, v, diff, diff1, diff2, remember[2*BRANCH];
    float  polygon;


    /* printf("loop start: %4d %4d %4d\n",i,j,cc); */   /* sm */

    count-=cc;
    i_old = i-1, j++;         /* j has now been set to the partner of the
			       previous pair for correct while-loop
			       termination.  */
    while (i != j)
    {
    	partner = pairs(i);
    	if (!partner)
        {
    	    i++, count++, bubble++;
                /* printf("\tif:   i=%d, j=%d, pairs(i)=%d, count=%d, bubble=%d\n",i, j, pairs(i), count, bubble); */
    	}
        else
        {
                /* printf("\telse: i=%d, j=%d, pairs(i)=%d, count=%d, bubble=%d\n",i, j, pairs(i), count, bubble); */
    	    count += 2;
    	    k = i, l = partner;		/* beginning of stack */
    	    remember[++r] = k;
    	    remember[++r] = l;
    	    i = partner+1;         	/* next i for the current loop */
    
    	    start_k = k, start_l = l;
    	    ladder = 0;
    	    do
            {
    		    k++, l--, ladder++;	/* go along the stack region */
    	    } while (pairs(k) == l);
    
    	    fill = ladder-2;
    	    if (ladder >= 2)
            {
        		angle[start_k+1+fill] += PIHALF;   /*  Loop entries and    */
        		angle[start_l-1-fill] += PIHALF;   /*  exits get an        */
        		angle[start_k]        += PIHALF;   /*  additional PI/2.    */
        		angle[start_l]        += PIHALF;   /*  Why ? (exercise)    */
        		if (ladder > 2)
                {
        		    for (; fill >= 1; fill--)
                    {
            			angle[start_k+fill] = PI;  /*  fill in the angles  */
            			angle[start_l-fill] = PI;  /*  for the backbone    */
        		    }
        		}
    	    }
    	    this_stack_size[++stk] = ladder;
    	    loop(k, l, 0);
    	}
    }
    remember[++r] = j;
    begin = i_old < 0 ? 0 : i_old;
    if (straight==TRUE) {
       if ((r==3) && (begin>0) && ((remember[1]-begin)!=(remember[3]-remember[2]))) {
       								/* printf("asymm. interner loop/bulge\n\tS : count = %d; ", count); */
          count   = remember[1]-begin+remember[3]-remember[2];	/* i.e. loop size +2				*/
          polygon = PI*(count-2)/(float)count;			/* Innenwinkel					*/
       								/* printf("diff1+diff2 = %d\n", count); 	*/
          diff1   = remember[1]-begin;
       								/* printf("\t1.: diff = %d = remember[1] - begin = %d - %d\n", diff1, remember[1], begin); */
       								/* printf("\t    polygon = %f = %f; count = %d\n", polygon, polygon*180./PI, diff1-1); */
          for (fill = 1; fill < diff1; fill++)
             angle[begin+fill] += polygon;
          angle[begin      ] += PIHALF+(PI-polygon)*(diff1-1)/2.;	/* Korrektur fuer Eintrittswinkel der Helix	*/
          angle[begin+diff1] += PIHALF+(PI-polygon)*(diff1-1)/2.;	/* Korrektur fuer Austrittswinkel der Helix	*/
       								/* for (fill = 0; fill <= diff1; fill++)	*/
       								/*    printf("\t    angle(%d) = %f\n", begin+fill, angle[begin+fill]*180./PI); */
          begin   = remember[2];
          diff2   = remember[3]-begin;
       								/* printf("\t2.: diff = %d = remember[3] - begin = %d - %d\n", diff2, remember[3], begin); */
       								/* printf("\t    polygon = %f = %f; count = %d\n", polygon, polygon*180./PI, diff2-1); */
          for (fill = 1; fill < diff2; fill++)
             angle[begin+fill] += polygon;
          angle[begin      ] += PIHALF+(PI-polygon)*(diff2-1)/2.;
          angle[begin+diff2] += PIHALF+(PI-polygon)*(diff2-1)/2.;
       								/* for (fill = 0; fill <= diff2; fill++) */
       								/*    printf("\t    angle(%d) = %f\n", begin+fill, angle[begin+fill]*180./PI); */
       } else {
       								/* if ((r==3) && (remember[1]-begin==remember[3]-remember[2])) {	*/
       								/* 	 printf("symm. interner loop\n");				*/
       								/* } else {								*/
       								/*    if (r==1) {							*/
       								/* 	 printf("hairpin\n");						*/
       								/*    } else {								*/
								/*   	 printf("bifurcation\n");					*/
								/*    }									*/
								/* }									*/
          polygon = PI*(count-2)/(float)count; 			/* bending angle in loop polygon */
          for (v = 1; v <= r; v++) {
	     diff = remember[v]-begin;
    	  							/* printf("\tdiff = %d = remember[%d] - begin = %d - %d\n", diff, v, remember[v], begin); */
    	  							/* printf("\tpolygon = %f = %f; count = %d\n", polygon, polygon*180./PI, count); */
	     for (fill = 0; fill <= diff; fill++)
	        angle[begin+fill] += polygon;
          							/* for (fill = 0; fill <= diff; fill++) */
          							/*    printf("\tangle(%d) = %f\n", begin+fill, angle[begin+fill]*180./PI); */
	     if (v > r)
	        break;
	     begin = remember[++v];
          }
       }
    } else {					/* par->straight==FALSE */
       polygon = PI*(count-2)/(float)count;
       for (v = 1; v <= r; v++) {
	  diff = remember[v]-begin;
	  for (fill = 0; fill <= diff; fill++)
	     angle[begin+fill] += polygon;
	  if (v > r)
	     break;
	  begin = remember[++v];
       }
    }
    this_loop_size[++lp] = bubble;
								/* printf("loop end : %4d %4d %4d\n",i,j,bubble); */	/* sm */
}

/*---------------------------------------------------------------------------*/
#define PRINT_SIZES 0
static void
parse(char *string)
{
    int k;

    for (k = 0; k < MAXLENGTH+5; k++)
	angle[k] = 0., x[k] = 0., y[k] = 0.;

    lp = stk = 0;

    /* upon exit from loop:

		  lp-1                number of loops in the structure
       this_loop_size[i], i=1,...,lp-1     number of unpaired digits in the
				      i-th loop
       this_loop_size[lp]-2                number of external digits (free ends
				      and joins)
		  stk                 number of stacks in the structure
       this_stack_size[i], i=1,...,stk     number of pairs in the i-th stack
    */

    loop(0, strlen(string)+1, 4);

#if PRINT_SIZES
    this_loop_size[lp] -= 2;     /* correct for cheating with function loop */
    for (k = 1; k <= lp; k++)
	printf("loop %d has size %d\n", k, this_loop_size[k]);
    for (k = 1; k <= stk; k++)
	printf("stack %d has length %d\n", k, this_stack_size[k]);
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


int CS_DrawStruct(struct pair_prob *pair, char *string, PAR *par, float  coord[], struct xyplot *my_xyplot)
{
   
  float  alpha, xmin, xmax, ymin, ymax, alpha_i, delta, charsize;
  int    i, k, l, lmax, length, tag, incr, first;
  int	 par0, par9;
  float	 para, rfactor;

  length = strlen(string);
  if (length>MAXLENGTH) {
     printf("Sequence too long, not doing xy_plot\n");
     return(FALSE);
  }


/*   printf("\nCheck auf Korrektheit der Bp ..."); */
  for (i=1; i<=length; i++) {
     if ((pair[i].pair != 0) && (pair[pair[i].pair].pair != i)) {
        printf("Inconsistent pairing list: %d:%d != %d:%d\n", i, pair[i].pair, pair[i].pair, pair[pair[i].pair].pair);
        /*return(FALSE);*/
     }
     if ((pair[i].pair != 0) && (i == pair[i].pair)) {
        printf("Inconsistent pairing list: %d:%d !\n", i, pair[i].pair);
        /*return(FALSE);*/
     }
     if (pair[i].pair != 0) {
        l = 0;
        for (k=1; k<=length; k++)
           if (pair[k].pair==pair[i].pair) l++;
        if (l>1) {
           printf("Inconsistent pairing list: %d:%d && %d:%d!\n", i, pair[i].pair, k, pair[k].pair);
           /*return(FALSE);*/
        }
     }
  }
/*   printf(" ok\n"); */

/*   printf("\nEinzelne Bp loeschen ...");
  if (par->single==FALSE) {
     for (i=2; i<length; i++) {
        if (pair[i].pair != 0) {
           k = pair[i].pair;
           if ((pair[i-1].pair != k+1) && (pair[i+1].pair != k-1)) {
              * printf("\nBp(%d:%d) and Bp(%d:%d) deleted ...", i, k, k, pair[k].pair); *
              pair[k].pair = 0;
              pair[i].pair = 0;
           }
        }
     }
  }
printf(" ok\n"); */


/*   printf("\nSeq ueberweisen ..."); */
  for (i=0; i<length; i++) {
     my_xyplot[i].pair = pair[i+1].pair;
     my_xyplot[i].nuc  = *(string+i);
  }
/*   printf(" ok\n"); */


/*   printf("\nUebergabe pair -> bond ..."); */

  base_pair    = (struct bond *) Scalloc(1+length/2, sizeof(struct bond));
  base_pair->j = 0;
  base_pair++;
  l            = 0;
  tag          = 1;
  incr         = FALSE;
  first        = TRUE;
  for ( i = 1; i <= length; i++ ) {
      if ( pair[i].pair != 0 ) {
      	  incr         = TRUE;
	  base_pair->i = i;
	  base_pair->j = pair[i].pair;
	  if ( first || pair[i-1].pair == pair[i].pair+1 ) {
	     if ( first ) {
	        my_xyplot[     i      -1].tag = -1*tag;
	        my_xyplot[pair[i].pair-1].tag = -1*tag;
	        first = FALSE;
	     } else {
	        my_xyplot[     i      -1].tag =    tag;
	        my_xyplot[pair[i].pair-1].tag =    tag;
	     }
	  } else {
	     tag++;
	        my_xyplot[     i      -1].tag = -1*tag;
	        my_xyplot[pair[i].pair-1].tag = -1*tag;
      	  }
	  base_pair++;
	  pair[pair[i].pair].pair=0;
	  l++;
      } else {
      	  if ( incr ) {
      	     tag++;
      	     incr  = FALSE;
      	     first = TRUE;
      	  }
      }
  }

  base_pair-=(l+1);
  base_pair->i=l;
/*
 *   for ( i = 0; i < length; i++ ) {
 *      printf("i:j = %d:%d, tag = %d\n", i+1, my_xyplot[i].pair, my_xyplot[i].tag);
 *   }
 */
/*   printf(" ok\n"); */


/*   printf("\nparse(string) ..."); */
  straight = par->straight;
  parse(string);
/*   printf(" ok\n"); */

  
  alpha = INIT_ANGLE;
  x[0]  = INIT_X;
  y[0]  = INIT_Y;
  xmin  = xmax = x[0];
  ymin  = ymax = y[0];


/*  printf("\nPivot (%d) ...", par->pivot); */
  if (par->pivot==TRUE && par->numangle>0 ) {
    for (l=0; l<par->numangle; l++) {
	if ( par->angle0[l]<1       ) par->angle0[l]=1;		/* erstes  Nukleotid */
	if ( par->angle9[l]>=length ) par->angle9[l]=length;	/* letztes Nukleotid */
	if ( par->angle9[l]==0      ) par->angle9[l]=length;
	par->angle[l]*=PI/180.;					/* Drehwinkel */
    }

    for ( l=1; l<par->numangle; l++ ) {		/* Sortieren durch Einfuegen */
    						/* => par->angle0[i] < par->angle0[i+1] */
	par0 = par->angle0[l];
	par9 = par->angle9[l];
	para = par->angle[ l];
	i    = l-1;
	while ( i>=0 && par0 < par->angle0[i] ) {
	    par->angle0[i+1] = par->angle0[i];
	    par->angle9[i+1] = par->angle9[i];
	    par->angle[ i+1] = par->angle[ i];
	    i--;
	}
	par->angle0[i+1] = par0;
	par->angle9[i+1] = par9;
	par->angle[ i+1] = para;
    }

    l=0;
    for (i = 1; i <= length; i++) {
	if ( l<par->numangle && i==par->angle0[l] ) {		/* i=[0,length-1] !!! d.h. jetzt 5' Helixrandnukleotid+1 */
	    delta   = alpha+par->angle[l];
	    x[i]    = x[i-1]+RADIUS*cos(-delta);		/* 2. Stapel knicken */
	    y[i]    = y[i-1]+RADIUS*sin(-delta);

	    k       = par->angle9[l];				/* 3' Helixrandnukleotid+1 */

/*
 * 	    printf("l=%d\ti=%d\tk=%d\t",l,i,k);
 * 	    printf("cos(-alpha-pi/2 )=cos(%f)=%f\t",-alpha-PIHALF,  cos(-alpha-PIHALF)  );
 * 	    printf("cos(-alpha+angle)=cos(%f)=",    -alpha+angle[k]                     );
 * 	    printf(                          "%f\n",                cos(-alpha+angle[k]));
 */
/*
 *                    | 5' Helixrandnukleotid
 * 		      |			 | + 90 Grad = 3' Helixrandnukleotid
 * 		      |			 |		    | + naechster Winkel = 3' Helixrandnukleotid+1
 * 		      |			 |		    |
 */
	    x[k]    = x[i-1]+RADIUS*(cos(-alpha-PIHALF)+cos(-alpha+angle[k]));
	    y[k]    = y[i-1]+RADIUS*(sin(-alpha-PIHALF)+sin(-alpha+angle[k]));
	
	    alpha_i = alpha;
	    alpha   = delta+PI-angle[i+1];			/* mit Knick weiter */
	    par->angle[l] = alpha_i-angle[k];			/* ohne Knick ab 3' Helixrandnukleotid+1 */
	    l++;
	} else {
	    if ( x[i]==0. ) {					/* wie ueblich */
	        x[i]  = x[i-1]+RADIUS*cos(-alpha);
		y[i]  = y[i-1]+RADIUS*sin(-alpha);
		alpha+= PI-angle[i+1];
	    } else {						/* 3' Helixrandnukleotid+1 */
	        for (k=0; k<par->numangle; k++) {		/* siehe oben: x[k], y[k] */
		    if ( i==par->angle9[k] ) {
			    alpha = par->angle[k];		/* Winkel ohne Knick rueckspeichern */
			    alpha+= PI-angle[i+1];		/* und normal weiter */
			    break;				/* raus aus FOR-Schleife */
		    }
		}
	    }
	}
    }
  } else {					/* no pivot */
    for (i = 1; i <= length; i++) {
	x[i]  = x[i-1]+RADIUS*cos(-alpha);	/* korr. 5'-> 3' orient.: alpha -> -alpha */
	y[i]  = y[i-1]+RADIUS*sin(-alpha);	/* korr. 5'-> 3' orient.: alpha -> -alpha */
	alpha+= PI-angle[i+1];
    }
  }
/*   printf(" ok\n"); */


/*   printf("\nmin/max ..."); */
  for (i = 1; i < length; i++) {
     xmin = x[i] < xmin ? x[i] : xmin;
     xmax = x[i] > xmax ? x[i] : xmax;
     ymin = y[i] < ymin ? y[i] : ymin;
     ymax = y[i] > ymax ? y[i] : ymax;
  }
/*   printf(" ok\n"); */


/*   printf("\nscale factor ..."); */
  if (par->scf != 0.0) {				/* abs. scaling */  /* sm */
     xmax = xmin + par->scf;
     ymax = ymin + par->scf;
  } else {
/*
 *      printf( "MaxScalingFactor %.1f \n", MAX((fabs(xmin)+fabs(xmax)),(fabs(ymin)+fabs(ymax))) );
 */
  }
/*   printf(" ok\n"); */


/*   printf("\ntext (%d) ...", par->text); */
  charsize = .75*RADIUS;				/* char size = .5*bond length */
  xmin    -= 2.*1.5*RADIUS + 4.*charsize;		/* add space for numbering */
  xmax    += 2.*1.5*RADIUS + 4.*charsize;
  ymin    -= 2.*1.5*RADIUS + 4.*charsize;
  ymax    += 2.*1.5*RADIUS + 4.*charsize;
  if ( par->text==TRUE ) {				/* space for 4 top lines (4*0.8 cm) */
     ymax+=6.*charsize;
     coord[4] = (xmin+xmax)/2.;				/* PLOTSTRUCTURE of: squfile; month day, year; hour:min */
     coord[5] = ymax-charsize;
     coord[6] = (xmin+xmax)/2.;
     coord[7] = ymax-3.*charsize;			/* FOLD of : ... */
     coord[8] = (xmin+xmax)/2.;
     coord[9] = ymax-5.*charsize;			/* Length: ... */
  } else {
     for ( i=4; i<=9; i++ ) coord[i]=0.;
  }
  coord[ 0] = xmin;					/* size of graph */
  coord[ 1] = ymin;
  coord[ 2] = xmax;
  coord[ 3] = ymax;
  coord[10] = charsize;					/* char size = .5*bond length */
  coord[11] = charsize/2.;				/* char width */
/*   printf(" ok\n"); */


/*   printf("\nbackbone (%d) ...", par->backbone); */
  if (par->backbone==TRUE) {				/* draw backbone */
     /* printf(" draw backbone ..."); */
     for ( i = 0; i < length; i++ ) {
	my_xyplot[i].back_x1 = x[i  ];
	my_xyplot[i].back_y1 = y[i  ];
	my_xyplot[i].back_x2 = x[i+1];				/* ??????????????????? */
	my_xyplot[i].back_y2 = y[i+1];				/* ??????????????????? */
     }
  }
/*   printf(" ok\n"); */


/*   printf("\nnuc (%d) ...", par->nuc); */
  if (par->nuc==TRUE) {					/* draw nucleotides */
    if ( par->numnuc>0 ) {
	for (k=0; k<par->numnuc; k++) {
	    if ( par->nuc0[k]<1       ) par->nuc0[k]=1;
	    if ( par->nuc9[k]>=length ) par->nuc9[k]=length;
	    if ( par->nuc9[k]==0      ) par->nuc9[k]=length;
	    for (i=-1+par->nuc0[k]; i<par->nuc9[k]; i+=10 ) {
		lmax = (i+10)<par->nuc9[k]?(i+10):par->nuc9[k];
		for (l = i;   l < lmax; l++ ) {
/*
 * 		    if (*(string+l)!='X' && *(string)!='x') {
 */
			if ( l==0 ) {
			   my_xyplot[l].nuc_x = x[l]+(y[l]  -y[l+1])/2.;
                           my_xyplot[l].nuc_y = y[l]+(x[l+1]-x[l]  )/2.;
                        } else {
			   my_xyplot[l].nuc_x = x[l]+(y[l-1]-y[l+1])/4.;
                           my_xyplot[l].nuc_y = y[l]+(x[l+1]-x[l-1])/4.;
		        }
/*
 * 		    }
 */
		}
	    }
	}
    } else {
	for (i = 0; i < length; i+=10 ) {
	    lmax = (i+10)<length?(i+10):length;
	    for (l = i;   l < lmax; l++ ) {
		if (*(string+l)!='X' && *(string)!='x') {
		    my_xyplot[l].nuc_x = x[l]+(y[l]  -y[l+1])/2.;
                    my_xyplot[l].nuc_y = y[l]+(x[l+1]-x[l]  )/2.;
		}
	    }
	}
    }
  }
/*   printf(" ok\n"); */


/*   printf("\nbonds (%d) ...", par->bonds); */
  if ( par->bonds==TRUE ) {				/* draw bonds */
    for (i = 0; i < length; i++ ) {
	l = pairs(i);
	if (l) {
           my_xyplot[i-1].bond_x1 = x[i-1];
           my_xyplot[i-1].bond_y1 = y[i-1];
           my_xyplot[i-1].bond_x2 = x[l-1];
           my_xyplot[i-1].bond_y2 = y[l-1];
	}
    }
  }
/*   printf(" ok\n"); */


/*   printf("\nnumbers (%d) ...", par->numbers); */
  if (par->numbers>0) {					/* draw numbering */
    if (par->nuc==TRUE)	rfactor=0.5;
    else		rfactor=0.0;
    if (*(string)!='X' && *(string)!='x') {
          my_xyplot[0].line_x1 = x[0]+(y[length-1]-y[       1])*rfactor;
          my_xyplot[0].line_y1 = y[0]+(x[       1]-x[length-1])*rfactor;
          my_xyplot[0].line_x2 = x[0]+ y[length-1]-y[       1];
          my_xyplot[0].line_y2 = y[0]+ x[       1]-x[length-1];
          my_xyplot[0].lpos_x  = x[0]+(y[length-1]-y[       1])*1.5;
          my_xyplot[0].lpos_y  = y[0]+(x[       1]-x[length-1])*1.5;
          my_xyplot[0].label   = 1;
    }
    if (par->numbers==1) k=-1+2*par->numbers;
    else		 k=-1+  par->numbers;
    for (i = k; i < length; i+=par->numbers )
       if (*(string+i)!='X' && *(string+i)!='x') {
          my_xyplot[i].line_x1 = x[i]+(y[i-1]-y[i+1])*rfactor;
          my_xyplot[i].line_y1 = y[i]+(x[i+1]-x[i-1])*rfactor;
          my_xyplot[i].line_x2 = x[i]+ y[i-1]-y[i+1];
          my_xyplot[i].line_y2 = y[i]+ x[i+1]-x[i-1];
          my_xyplot[i].lpos_x  = x[i]+(y[i-1]-y[i+1])*1.5;
          my_xyplot[i].lpos_y  = y[i]+(x[i+1]-x[i-1])*1.5;
          my_xyplot[i].label   = i+1;
	}
  }
/*   printf(" ok\n"); */


/*   printf("\nrange (%d) ...", par->range); */
  if (par->range==TRUE) {				/* draw modified numbering */
    if (par->nuc==TRUE)	rfactor=0.5;
    else		rfactor=0.0;
    if (par->range0==1) k=-1+2*par->range0;
    else		k=-1+  par->range0;
    for (i = k; i < length; i+=par->rangestep )
	if (*(string+i)!='X' && *(string)!='x')  {
          my_xyplot[i].line_x1 = x[i]+(y[i-1]-y[i+1])*rfactor;
          my_xyplot[i].line_y1 = y[i]+(x[i+1]-x[i-1])*rfactor;
          my_xyplot[i].line_x2 = x[i]+ y[i-1]-y[i+1];
          my_xyplot[i].line_y2 = y[i]+ x[i+1]-x[i-1];
          my_xyplot[i].lpos_x  = x[i]+(y[i-1]-y[i+1])*1.5;
          my_xyplot[i].lpos_y  = y[i]+(x[i+1]-x[i-1])*1.5;
          my_xyplot[i].label   = i + 1 + par->range1 - par->range0;
	}
  }
/*   printf(" ok\n"); */


   free((struct bond *)base_pair);
/*   printf("RETURN\n"); */
   return(TRUE);
}
