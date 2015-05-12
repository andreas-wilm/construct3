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
 *  CVS $Id: dlist.c,v 1.11 2004-05-25 13:29:14 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* memcmp */

#include "dlist.h"



/***   dlist_new   ***********************************************************
 *
 * construct a dlist
 * return NULL on error
 * user must free (dlist_destroy)
 *
 */
dlist *
dlist_new(size_t elem_data_size, void (*destroy)(void *data))
{
    dlist *list = (dlist *)calloc(1, sizeof(dlist));
    if (list==NULL)
    {
        fprintf(stderr, "ERROR(%s|%s): Memory allocation failed\n",
                                             __FILE__, __FUNCTION__);
        return NULL;
    }
    list->size           = 0;
    list->elem_data_size = elem_data_size;
    list->destroy        = destroy;
    list->head           = NULL;
    list->tail           = NULL;
    
    return list;
}
/***   dlist_new   ***/




/***   dlist_destroy   ********************************************************
 * 
 * destroy a dlist,
 * using the user defined (dlist_init) destroy function for freeing the data
 *
 */
void
dlist_destroy(dlist *list)
{
    void *data;
    
    while (dlist_size(list) > 0)
    {
        if (dlist_remove(list, dlist_tail(list), (void **)&data) == DLIST_OK
            &&
            list->destroy != NULL)
        {
            /*  Call a user-defined function to
               free dynamically allocated data.
            */
            list->destroy(data);
        }
    }
     
    free(list);
    
    return;
}
/***   dlist_destroy   ***/




/***   dlist_ins_next   *******************************************************
 *
 * insert data next to <dlist_elem *element>
 * return DLIST_OK on success, DLIST_ERROR otherwise
 * actually allocates new mem and memcpy's the data
 */
int
dlist_ins_next(dlist *list, dlist_elem *element, const void *data)
{
    dlist_elem *new_element;

    /* Do not allow a NULL element unless the list is empty
    */
    if (element == NULL && dlist_size(list) != 0)
        return DLIST_ERROR;


   /* Allocate storage for the element.
    */
    if ((new_element = (dlist_elem *)malloc(sizeof(dlist_elem))) == NULL)
        return DLIST_ERROR;

    /* Insert the new element into the list.
     */
    new_element->data = (void *)data;

    if (dlist_size(list) == 0)
    {
        /*  Handle insertion when the list is empty.
         */
       list->head       = new_element;
       list->head->prev = NULL;
       list->head->next = NULL;
       list->tail       = new_element;
    }
    else
    {
        /*  Handle insertion when the list is not empty.
        */
       new_element->next = element->next;
       new_element->prev = element;
    
       if (element->next == NULL)
          list->tail = new_element;
       else
          element->next->prev = new_element;
    
       element->next = new_element;
    }

    /*  Adjust the size of the list to account for the inserted element.
     */
    list->size++;

    return DLIST_OK;
}
/***   dlist_ins_next   ***/




/***   dlist_ins_prev   *******************************************************
 *
 * insert data prior to <dlist_elem *element>
 * return DLIST_OK on success, DLIST_ERROR otherwise
 *
 *
 */
int
dlist_ins_prev(dlist *list, dlist_elem *element, const void *data)
{
    dlist_elem *new_element;

    /*  Do not allow a NULL element unless the list is empty.
     */

    if (element == NULL && dlist_size(list) != 0)
        return DLIST_ERROR;

    /*  Allocate storage to be managed by the abstract data type.
     */

    if ((new_element = (dlist_elem *)malloc(sizeof(dlist_elem))) == NULL)
        return DLIST_ERROR;

    /*  Insert the new element into the list.
     */
    new_element->data = (void *)data;


    if (dlist_size(list) == 0)
    {
        /*  Handle insertion when the list is empty.
         */
       list->head = new_element;
       list->head->prev = NULL;
       list->head->next = NULL;
       list->tail = new_element;
    }
    else
    {
        /*  Handle insertion when the list is not empty.
        */
        new_element->next = element; 
        new_element->prev = element->prev;

        if (element->prev == NULL)
            list->head = new_element;
        else
            element->prev->next = new_element;

        element->prev = new_element;
    }

    /*  Adjust the size of the list to account for the new element.
     */
    list->size++;

    return DLIST_OK;
}
/***   dlist_ins_prev   ***/




/***   dlist_remove   *********************************************************
 *
 * Remove the element pointed to by <dlist_elem *elemen> 
 * return DLIST_OK on success, DLIST_ERROR otherwise
 *
 */
int
dlist_remove(dlist *list, dlist_elem *element, void **data)
{
    /*  Do not allow a NULL element or removal from an empty list.
     */
    if (element == NULL || dlist_size(list) == 0)
        return DLIST_ERROR;

    /*  Remove the element from the list.
     */
    *data = element->data;

    if (element == list->head)
    {
        /*  Handle removal from the head of the list.
         */
       list->head = element->next;

        if (list->head == NULL)
            list->tail = NULL;
        else
            element->next->prev = NULL;
    }
    else
    {
        /*  Handle removal from other than the head of the list.
         */
        element->prev->next = element->next;

        if (element->next == NULL)
            list->tail = element->prev;
        else
            element->next->prev = element->prev;
    }

    /*  Free the data and storage allocated by the abstract data type.
     */
    free(element);


    /*  Adjust the size of the list to account for the removed element.
     */
    list->size--;

    return DLIST_OK;
}
/***   dlist_remove   ***/




/***   dlist_append   *********************************************************
 *
 * append data to the list's tail
 * return DLIST_OK on success, DLIST_ERROR otherwise
 */
int
dlist_append(dlist *list, const void *data)
{
    if (dlist_ins_next(list, dlist_tail(list), data) != 0)
        return DLIST_ERROR;
    else
        return DLIST_OK;
}
/***   dlist_append   ***/




/***   dlist_prepend   ********************************************************
 *
 * insert data at list's head
 * return DLIST_OK on success, DLIST_ERROR otherwise
 */
int
dlist_prepend(dlist *list, const void *data)
{
    if (dlist_ins_prev(list, dlist_head(list), data) != 0)
        return DLIST_ERROR;
    else
        return DLIST_OK;
}
/***   dlist_prepend   ***/

