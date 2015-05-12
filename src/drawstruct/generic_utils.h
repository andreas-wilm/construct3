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
 *  CVS $Id: generic_utils.h,v 1.4 2004-05-25 13:29:22 wilm Exp $    
 */

#ifndef GENERIC_UTILS_H_INCL
#define GENERIC_UTILS_H_INCL

/*** console io: append trailing "\n"
 */
/* print only if cs_opts.debug is true*/
#define DEBUG_P(fmt, args...)     printk(stdout, 1, "DEBUG(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* print only if cs_opts.verbose is true*/
#define VERBOSE_P(fmt, args...)   printk(stdout, 1, fmt, ## args)
/* always warn to stderr */
#define WARN_P(fmt, args...)      printk(stdout, 1, "WARNING(%s): " fmt, __FILE__, ## args)
/* always print errors to stderr*/
#define ERROR_P(fmt, args...)     printk(stderr, 1, "ERROR(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print temp debugging */
#define TMPDEBUG_P(fmt, args...)  printk(stdout, 1, "TMPDEBUG(%s|%s|%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)


/*****   string handling
 *
 */
/* INT_TO_STRLEN: how many characters has an int (excl. trailing '\0' !) */
#define INT_TO_STRLEN(a)    ((int)(ceil(log10((double)(a))))+1)
/* STR_EQ: strings are equal , case sensitive */
#define STR_EQ(a,b)         (strcmp((a),(b)) == 0)
/*  STR_NC_EQ:  strings are equal , ignoring case */
#define STR_NC_EQ(a,b)      (strcasecmp((a),(b)) == 0)


/*****   math
 *
 */
#define MAX(x,y)            (((x)>(y)) ? (x) : (y))
#define MIN(x,y)            (((x)<(y)) ? (x) : (y))
#define ABS(x)              (((x)<0.) ? (x*-1.) : (x))



#define IS_GAP(c)           ((c) == '-')


enum bool {FALSE=0, TRUE=1, ERROR=-1, OK=1, NO=0, YES=1};




/* Print given error message
   and exit
  */
extern void
Die(char *msg);


/* Console IO
 * Taken from the Linux kernel source
 */
extern int
printk(FILE *stream, int flag, const char *fmt, ...);          


/* Returns TRUE if file exists
 */
extern int
FileExists (char *p2f);

/* convert string to lowercase
 */
extern void
StrTolower(char *s);

/* safe calloc
 */
extern void *
Scalloc(size_t nmemb, size_t size);

/* safe malloc
 */
extern void *
Smalloc(size_t size);

#endif
