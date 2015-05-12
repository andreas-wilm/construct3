/******************************************************************************
* 
* cs_seqio.c - ConStruct interface to seqio.c
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
 *  CVS $Id: cs_seqio.c,v 1.39 2007-10-22 10:43:23 steger Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "public_datatypes.h"
#include "main.h"
#include "generic_utils.h"
#include "stopwatch.h"
#include "utils.h"
#include "mx.h"

#include "consensus.h"
#include "bp_prob_mat.h"
#include "seqio.h"
#include "cs_seqio.h"



/***   private
 */
 
#if 0
    #define DEBUG
#endif

#define   VAL_SEQNT_CHAR       "aynur-cg"
#define INVAL_SEQID_CHAR       "!\"§$%&/()=?{[]}\\'`+*~'#;:<>|"

/* supported output formats:
 * some missing: if you got time fix it and change FILE_TYPES(seqout) accordingly
 * extensions should be the same as in tcl except missing vie (written by tcl) 
 * 1. idx: exact format string as in seqio.c:file_table
 * 2. idx: corresponding extension
 */
static int   seqio_nformats = 5;
static char *seqio_outformats[][2] = {
                   {"FASTAold",  ".fa"}, {"Clustal",".aln"},
                   {"MSF",      ".msf"}, {"PHYLIP",   ".phy"},
                   {"ASN.1",    ".asn"}};


static void
KillSeq(seq_t *seq);

static void
KillConsSeq(cons_seq_t *cs_seq);
 
static void
KillBp(seq_t *seq);
/*
 ***/




/***   KillConsSeq   **********************************************************
 *
 */
void
KillConsSeq(cons_seq_t *cs_seq)
{
    FreeConsBpMx(cs_seq->cons_bp, strlen(cs_seq->nt));
    
    if (cs_seq->id!=NULL)
        free(cs_seq->id);
    if (cs_seq->nt!=NULL)
        free(cs_seq->nt);    
    if (cs_seq->num_gaps!=NULL)
        free(cs_seq->num_gaps);
}
/***   KillConsSeq   ***/



/***   KillBp   ***************************************************************
 *
 */
void
KillBp(seq_t *seq)
{
    char *degapped_seq;
    
    degapped_seq = (char*) Scalloc(strlen(seq->nt)+1, sizeof(char));    
    Degap(seq->nt, degapped_seq);
    
    FreeBpProbMx(seq->bp, strlen(degapped_seq));
    
    seq->bp=NULL;
    
    free(degapped_seq);
}
/***   KillBp   ***/



/***   KillSeq   **************************************************************
 *
 */
void
KillSeq(seq_t *seq)
{
    if (seq->bp!=NULL)
        KillBp(seq);

    free(seq->nt);        
    free(seq->id);


}
/***   KillSeq   ***/




/***   NewAln   ***************************************************************
 *
 */
aln_t *
NewAln(void)
{
    aln_t *ret_aln = NULL;

    ret_aln = (aln_t *) Scalloc(1, sizeof(aln_t));

    ret_aln->seq      = NULL;
    ret_aln->cons_seq = NULL;
    
    return ret_aln;
}
/***   NewAln   ***/




/***   KillAln   **************************************************************
 *
 */
void
KillAln(aln_t *aln)
{
    int i=0;

    if (aln->cons_seq!=NULL)
        KillConsSeq(aln->cons_seq);
    for (i=0; i<aln->num_seq; i++)
        KillSeq(&aln->seq[i]);
    free(aln->seq);
    free(aln);
    aln = NULL;
}
/***   KillAln   ***/




/***   ReadAlnFile   **********************************************************
 *
 * Reads given alignment using seqio and returns entries as aln_t
 * or NULL on failure
 *
 */
aln_t *
ReadAlnFile(char *path2seqfile)
{

    int    n_seqs;                /* number of read sequences       */
    aln_t *ret_aln   = NULL;      /* returned alignment             */
    int  aln_len = 0;             /* alignment length               */
    int  old_len = 0;             /* saved len of previous read seq
                                     used for comparison            */
    int i;
    
    SEQFILE  *sfp       = NULL;   /* seq_io internal data structure */
    SEQINFO  *sip       = NULL;   /* seq_io internal data structure */
    



    VERBOSE_P("Reading Sequences from Alignment \"%s\n", path2seqfile);

    if (FileExists(path2seqfile)==FALSE)
        return NULL;

    
    if ((sfp = seqfopen2(path2seqfile)) == NULL)
        return NULL;


    ret_aln = NewAln();
    

    /* read all entries
     *
     */
    n_seqs = 0;
    while ((sip = seqfgetinfo(sfp, 0))!=NULL)
        n_seqs++;
    seqfclose(sfp);


    ret_aln->seq = (seq_t *) Scalloc(n_seqs, sizeof(seq_t));    
    
    
    
    /***   ids and SEQINFO
     */ 
    sfp = seqfopen2(path2seqfile);
    for (i=0; i<n_seqs; i++)
    {
        int idonly = 1;
        int idlen  = 32;

        /* ret_aln->seq[nof_read_seqs-1].id = seqfdescription(sfp, 1);
          or
          ret_aln->seq[nof_read_seqs-1].id = seqfmainid(sfp, 1);
          or
        */
        ret_aln->seq[i].id = (char *) Scalloc(idlen, sizeof(char));
        sip = seqfgetinfo(sfp, 0);

        seqfoneline(sip, ret_aln->seq[i].id, idlen, idonly);
        
        DEBUG_P("id round %d: got id \"%s\"\n", i, ret_aln->seq[i].id);
    }
    seqfclose(sfp);
    

    
    /***   nts
     */ 
    sfp = seqfopen2(path2seqfile);
    for (i=0; i<n_seqs; i++)
    {
        ret_aln->seq[i].nt = seqfgetrawseq(sfp, &aln_len, 1);
        StrTolower(ret_aln->seq[i].nt);

        /* check that all sequence lengths are the same
         */        
        if ( (aln_len!=old_len) && (i!=0) )
        {
            WARN_P("Sequences are not aligned (at least %s and %s differ (%d!=%d))",
                     	ret_aln->seq[i-1].id, ret_aln->seq[i].id, old_len, aln_len);
            KillAln(ret_aln);
            return NULL;
        }
        DEBUG_P("nt round %d: got nt \"%s\"\n", i, ret_aln->seq[i].nt);
        old_len = aln_len;
    }
    seqfclose(sfp);


    /* cs needs to check for invalid characters
     * causing problems while touching bp-prob files etc.
     */
    for (i=0; i<n_seqs; i++)
    {
        /* check for invalid characters
         */
        CheckSeqNt(ret_aln->seq[i].nt, ret_aln->seq[i].id);
        CheckSeqId(ret_aln->seq[i].id);
    }        
        

    DEBUG_P("n_seqs \"%d\"\n", n_seqs);
    DEBUG_P("aln_len \"%d\"\n", aln_len);
    
    ret_aln->len      = aln_len;
    ret_aln->num_seq  = n_seqs;
    
    return ret_aln;
}
/***   ReadAlnFile   ***/



        

/***   CheckSeqNt   ***********************************************************
 *
 * Check sequence nt for invalid chars (VAL_SEQNT_CHAR) and replace them
 *
 */
void
CheckSeqNt(char *seqnt, char *seqid)
{
    const char repl_t_char   = 'u';
    const char repl_any_char = 'n';
    char val_seq_char[1024];
    int invalid;
    int i, j;

    strcpy(val_seq_char, VAL_SEQNT_CHAR);
    
    /* exit if everything's fine
     */
    if (strspn(seqnt, VAL_SEQNT_CHAR) == strlen(seqnt))
        return;

    for (i=0; i<strlen(seqnt); i++)
    {
        invalid=1;
        for (j=0; j<strlen(val_seq_char); j++)
        {
            if ((char)tolower(seqnt[i])==val_seq_char[j])
            {
                invalid=0;
                break;
            }
        }
        if (invalid)
        {
            if ((char)tolower(seqnt[i])=='t')
            {
                /* WARN_P("Replacing invalid char %c (pos %d) in sequence %s with %c\n", seqnt[i], i+1, seqid, repl_t_char); */
                seqnt[i] = repl_t_char;
            }
            else
            {
                /* WARN_P("Replacing invalid char %c (pos %d) in sequence %s with %c\n", seqnt[i], i+1, seqid, repl_any_char); */
                seqnt[i] = repl_any_char;
            }
        }
    }
}
/***   CheckSeqNt   ***/



/***   CheckSeqId   ***********************************************************
 *
 * Search sequence id for invalid chars (INVAL_SEQID_CHAR) and replace them
 */
void
CheckSeqId(char *seqid)
{
    const char repl_char   = '_';
    char inval_id_char[1024];
    int i, j;
    
    
    /* exit if everything's fine
     */
    if (strcspn(seqid, INVAL_SEQID_CHAR) == strlen(seqid))
        return;
    
    strcpy(inval_id_char, INVAL_SEQID_CHAR);
    
    for (i=0; i<strlen(seqid); i++)
    {
        for (j=0; j<strlen(inval_id_char); j++)
        {
            if (seqid[i]==inval_id_char[j])
            {
                /* WARN_P("Replacing invalid char %c in sequence-id %s with %c\n", seqid[i], seqid, repl_char); */
                seqid[i] = repl_char;
            }
        }
    }
}
/***   CheckSeqId   ***/







/***   ConvertAndUpdateAln   **************************************************
 *
 * Read input alignment to get SEQIO internals
 * then update sequence
 * and save to requested format
 * 
 * internals don't have to be updated, since truelen and
 * rawlen must be the same (overall number of gaps hasn't changed)
 * f_ori has to exist
 * sequence order in <store_aln> must be the same as in <f_ori>
 */
int
ConvertAndUpdateAln(aln_t *store_aln, char *f_ori, char *f_out, char *outformat_ext)
{
    int i, entry, len;
    char format_out[1024];
    SEQINFO *info;
    SEQFILE *insfp, *outsfp;
    char *seq;
        

    /* find requested format in supported types
     */
    for (i=0; i<seqio_nformats; i++)
    {
        if (STR_EQ(outformat_ext, seqio_outformats[i][1]))
        {
            strcpy(format_out, seqio_outformats[i][0]);
            break;
        }
    }
    if (i==seqio_nformats)
    {
        ERROR_P("unsupported sequence file extension \"%s\"\n", outformat_ext);
        return ERROR;
    }
    
    
    /* paranoia
     */
    if (!seqfisaformat(format_out))
    {
        ERROR_P("ups...seqio does not seem to support this format (\"%s\")\n", format_out);
        return ERROR;    
    }
    
    
    /* opem files
     */
    if ((insfp = seqfopen2(f_ori)) == NULL)
    {
        seqfperror(NULL);
        ERROR_P("%s\n", "couldn't open input file");
        return ERROR;
    }
    if ((outsfp = seqfopen(f_out, "w", format_out)) == NULL)
    {
        seqfperror(NULL);
        ERROR_P("%s\n", "couldn't open output file");
        seqfclose(insfp);
        return ERROR;    
    }
    
    entry=0;
    while ((seq = seqfgetseq(insfp, &len, 0)) != NULL)
    {
        if ((info = seqfinfo(insfp, 0)) != NULL)
        {
            char buf[1024];
        
            /* FIXME: sanity check oneline the same ? */
            seqfoneline(info, buf, 32, 1);

            /* FIXME: unaligned seq the same ? */
            
            /* not needed since overall number of gaps
               cannot be changed inside construct */  
            seqfwrite(outsfp, aln->seq[entry].nt, strlen(aln->seq[entry].nt), info);
            entry++;
        }
        else
        {
            WARN_P("Ups..couldn't parse info for entry %d\n", entry+1);
            WARN_P("%s\n", "outputfile may be corrupted");
        }
    }
    
    
    seqfclose(insfp);
    seqfclose(outsfp);
    
        
    return OK;
}
/***   ConvertAndUpdateAln  ***/
