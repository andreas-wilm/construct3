/******************************************************************************
* 
* generic_utils.c - generic routines
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
 *  CVS $Id: generic_utils.c,v 1.7 2004-05-25 13:29:14 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "generic_utils.h"




/***   Die   ******************************************************************
 *
 *
 */
void
Die(char *msg)
{
    fprintf(stderr, "SEVERE ERROR: %s\n", msg);
	exit(EXIT_FAILURE);
}
/***   Die   ***/



/***   FileExists   ***********************************************************
 *
 * Returns TRUE if file exists
 */
int
FileExists (char *p2f)
{
    FILE *f;
    
    if (p2f==NULL)
        return FALSE;
        
    f = fopen(p2f, "r");
	
    if (f==NULL)
    {
        return FALSE;
    }
    else 
    {
        fclose(f);
        return TRUE;
    }
}
/***   FileExists   ***/




/***   StrTolower   ***********************************************************
 *
 */
void
StrTolower(char *s)
{
    int i=0;
    while (s[i]!='\0') {
        s[i]=tolower(s[i]);
        i++;
    }
}
/***   StrTolower   ***/




/***   Scalloc   ***********************************************************
 *
 * Safe calloc
 * Die if no memory avaiable
 * 
 */
void *
Scalloc(size_t nmemb, size_t size)
{
    void *ptr=NULL;
    if ( (ptr=calloc(nmemb, size)) ==NULL)
        Die("memory allocation failure");
    return ptr;
}
/***   Scalloc   ***/



/***   Smalloc   ***********************************************************
 *
 * Safe malloc
 * Die if no memory avaiable
 * 
 */
void *
Smalloc(size_t size)
{
    void *ptr=NULL;
    if ( (ptr=malloc(size)) ==NULL)
        Die("memory allocation failure");
    return ptr;
}
/***   Smalloc   ***/



/***   Srealloc   ***********************************************************
 *
 * Safe malloc
 * Die if no memory avaiable
 * 
 */
void *
Srealloc(void *ptr, size_t size)
{
    if ( (ptr=realloc(ptr, size)) ==NULL)
        Die("memory reallocation failure");
    return ptr;
}
/***   Srealloc   ***/


/***   ChompNTrimStr   *********************************************************
 *
 * Remove all leading and trailing 'spaces' (isspace) from string
 * and return a cleaned copy
 * Caller must free
 */
char *
ChompNTrimStr(const char *in_str)
{
    char *ret_str;
    int   i=0;
    
    if (in_str==NULL)
        return NULL;
    while (isspace(in_str[i]))
        i++;
    ret_str = strdup(&in_str[i]);
    while (isspace(ret_str[strlen(ret_str)-1]))
        ret_str[strlen(ret_str)-1]='\0';
    return ret_str;
}
/***   ChompNTrimStr   ***/




/***   ChompStr   ************************************************************
 *
 * remove trailing \n from str if non-NULL
 */
void
ChompStr(char *str)
{
    if (str==NULL)
        return;
    if (str[strlen(str)-1]=='\n')
        str[strlen(str)-1]='\0';
}
/***   ChompStr   ***/





/***   printk   ***************************************************************
 *
 * Taken from the Linux kernel source
 * and slightly modified
 *
 * flag: boolean: print or don't
 *       
 */
int
printk(FILE *stream, int bool_flag, const char *fmt, ...)
{                
    va_list args;
    static char printk_buf[8192];
    int printed_len;
 
    /* Emit the output into the temporary buffer */
    va_start(args, fmt);
    printed_len = vsnprintf(printk_buf, sizeof(printk_buf), fmt, args);
    va_end(args);
    
    if (bool_flag)
    {
        fprintf(stream, "%s", printk_buf);
        fflush(stream); 
        
    }
    return printed_len;
}
 /***   printk   ***/
 

