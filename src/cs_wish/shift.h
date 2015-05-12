/******************************************************************************
* 
* shift.h - compute cs-shift between two alignments
*
* Copyright (C) 2001-2004 Institute of Physical Biology,
*                    Heinrich-Heine Universit�t D�sseldorf, Germany
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
 *  CVS $Id: shift.h,v 1.5 2004-05-25 13:29:17 wilm Exp $    
 */


#ifndef SHIFT_H_INCL
#define SHIFT_H_INCL

extern int *
CsShiftScore(aln_t *aln_trust, aln_t *aln_pred, int *pairmask);

#endif
