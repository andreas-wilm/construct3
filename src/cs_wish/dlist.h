/******************************************************************************
* 
* dlist.c - procedures for double linked lists
*           Based upon
*           Mastering Algorithms with C
*           By Kyle Loudon
*           1st Edition August 1999 
*           ISBN: 1-56592-453-3
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
 *  CVS $Id: dlist.h,v 1.9 2004-05-25 13:29:14 wilm Exp $    
 */


/* FIXME: include comp function for sorted insertion */


#ifndef DLIST_H_INCL
#define DLIST_H_INCL


/* a list element
 */
typedef struct dlist_elem_t
{
    void                 *data;
    struct dlist_elem_t  *prev;
    struct dlist_elem_t  *next;
} dlist_elem;

/* the list itself
 */
typedef struct dlist_t
{
    int           size;
    size_t        elem_data_size;
    void          (*destroy)(void *data);
    dlist_elem    *head;
    dlist_elem    *tail;

} dlist;


/* return codes */
#define DLIST_ERROR -1
#define DLIST_OK     0


/* create a dlist */
dlist *
dlist_new(size_t elem_data_size, void (*destroy)(void *data));

/* destroy a dlist */
void
dlist_destroy(dlist *list);

/* insert data after an element */
int
dlist_ins_next(dlist *list, dlist_elem *element, const void *data);

/* insert data before an element */
int
dlist_ins_prev(dlist *list, dlist_elem *element, const void *data);

/* remove one element */
int
dlist_remove(dlist *list, dlist_elem *element, void **data);

/* append data to list */
int
dlist_append(dlist *list, const void *data);

/* prepend data to list */
int
dlist_prepend(dlist *list, const void *data);



/* macros may return NULL , so check */

#define dlist_size(list)       ((list)->size)

#define dlist_head(list)       ((list)->head)

#define dlist_tail(list)       ((list)->tail)

#define dlist_is_head(element) ((element)->prev == NULL ? 1 : 0)

#define dlist_is_tail(element) ((element)->next == NULL ? 1 : 0)

#define dlist_data(element)    ((element)->data)

#define dlist_next(element)    ((element)->next)

#define dlist_prev(element)    ((element)->prev)




#endif
