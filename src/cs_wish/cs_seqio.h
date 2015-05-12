/******************************************************************************
* 
* cs_seqio.h - ConStruct interface to seqio.c
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
 *  CVS $Id: cs_seqio.h,v 1.9 2004-05-25 13:29:14 wilm Exp $    
 */


#ifndef CS_SEQIO_H_INCL
#define CS_SEQIO_H_INCL




/* reads given sequence file and returns
   the aligned sequences as aln_t
 */
extern aln_t *
ReadAlnFile(char *path2seqfile);

extern aln_t *
NewAln(void);

extern void
KillAln(aln_t *aln);

extern void
CheckSeqNt(char *seqnt, char *seqid);

extern void
CheckSeqId(char *seqid);

extern int
ConvertAndUpdateAln(aln_t *store_aln, char *f_ori, char *f_out, char *outformat_ext);

#endif
