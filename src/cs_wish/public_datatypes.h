/******************************************************************************
* 
* public_datatypes.h - structs everyone needs
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
 *  CVS $Id: public_datatypes.h,v 1.3 2004-05-25 13:29:17 wilm Exp $    
 */



#include "dlist.h"



/***   bpair_t
 */
typedef struct
{
    float prob;
} bpair_t;



/***   seq_t
 */
typedef struct 
{
	char        *id;
	char        *nt;
    int         *gapshift;
    bpair_t   ***bp;        /* MX_DP: never access directly, use GetBpProbFromAlnSeq,
                               since indices are those of unaligned seq                */
} seq_t;


/***   cons_bpair_t
 */
typedef struct
{
    float     prob;       /* float DP: raw cons bp probability
                           * SUM_s[(pow(prob(s), 1.0/(double)POW_A)* weight[s])]
                           * Use GetConsBpProb to get real value */
    dlist    *contribs;   /* list of seq idxs */
} cons_bpair_t;



/***   cons_seq_t
 */
typedef struct 
{
	char           *id;
	char           *nt;
    int            *num_gaps;
    cons_bpair_t  **cons_bp;   /* cons_bpair_t DP: use GetConsBp to access */ 
} cons_seq_t;


/***   aln_t
 */
typedef struct
{
    int         num_seq;
    int         len;
    seq_t      *seq;
    cons_seq_t *cons_seq;
} aln_t;



/***   pair_t
 */
typedef struct 
{
	int ntj;
	int nti;
} pair_t;



/***   proj_entry_t
 */
typedef struct
{
	char   *seq_id;
	char   *f_bp_prob_mat;
	char   *fold_cmd;
 	float  weight;   
} proj_entry_t;



/***   proj_t
 */
typedef struct
{
	char          *name;
    char          *file;
	char          *version;

	int            num_entries;
    proj_entry_t  *entry;
    
	char          *aln_file;
    
} proj_t;



/***   tk_color_t
 */
typedef struct
{
    char *dp_bg;
    char *sel_seq;
    char *unsel_seq;
    char *act_nt;
    char *unact_nt;
    char *sel_bp;
    char *unsel_bp;
} tk_color_t;



/***   cs_opts_t
 */
typedef struct 
{
	int verbose;
	int debug;
    int do_timing;
} cs_opts_t;





