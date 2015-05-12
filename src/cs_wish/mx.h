/******************************************************************************
* 
* mx.h - generic (dotplot) matrix routines
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
 *  CVS $Id: mx.h,v 1.5 2004-05-25 13:29:16 wilm Exp $    
 */



/******************************************************************************
 ***
 *** Routines for creating matrices
 ***
 ***
 *** matrices are are constructed with mx_new()
 *** dependent on int flag type valid matrix indices are :
 *** default matrix        MX_DEFAULT     :   [i][j]; 0<=i<num_row, 0<=j<num_col
 *** dotplot               MX_DP          :   [i][j]; 0<=i<num_row, i>j   (num_col will be ignored)
 *** dotplot incl diagonal MX_DPINCLDIAG  :   [i][j]; 0<=i<num_row, i>=j  (num_col will be ignored)
 ***
 *** otherwise the matrix represents a dotplot :
 ***                       [j][i]; 0<=j<nr, i<j
 ***
 *** on creation you have to cast
 *** matrices are freed with mx_free()
 ***
 ***
 *** EXAMPLE: initializing a dotplot with -1:
 ***    int **dotplot;
 ***    int i, j;
 ***    dotplot = (int**) mx_new(seqlen, 0, sizeof(int), sizeof(int*), MX_DP);
 ***    for (j=0; j<seqlen; j++)
 ***    {
 ***        for (i=j+1; i<seqlen; i++)
 ***        {
 ***            dotplot[i][j] = -1;
 ***    ...
 ***   mx_free((void**)dotplot, seqlen, MX_DP);
 ***
 ***
 *****************************************************************************/


#ifndef MATRIX_H_INCL
#define MATRIX_H_INCL

#define MX_DEFAULT  0
#define MX_DP       1
#define MX_DPDIAG   2

void **
mx_new (int num_row, int num_col, size_t size_elem, size_t size_pointer2elem, int type);

extern void
mx_free (void **mx, int nr, int type);

#endif
