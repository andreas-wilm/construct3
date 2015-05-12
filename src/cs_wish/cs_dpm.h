#ifndef CS_DPM_H_INCL
#define CS_DPM_H_INCL

/******************************************************************************
* 
* cs_dpm.h - Procedures to construct a merged consensus matrix
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
 *  CVS $Id: cs_dpm.h,v 1.9 2004-05-25 13:29:14 wilm Exp $    
 */


/***   CsDpm
 *     Returns float matrix, or NULL on error
 *     data = (prob_fac * prob, if prob > prob_thr) +
 *            (mic_fac + mic,   if mic > mic_thr)
 */
extern float **
CsDpm(float prob_fac, float mic_fac, float prob_thr, float mic_thr);


extern void
FreeCsDpm(float **csdpm);

extern  void
PrintCsDpm(char *csseq, float **csdpm, FILE *stream);

#endif

