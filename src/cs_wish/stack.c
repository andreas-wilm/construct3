/******************************************************************************
* 
* stack.c - generic stack routines
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
 *  CVS $Id: stack.c,v 1.3 2004-05-25 13:29:17 wilm Exp $    
 */


#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>

#include "stack.h"



#define STACK_ALLOC_ERROR fprintf(stderr, "ERROR(%s|%s): no more memory\n", __FILE__, __FUNCTION__)
#define STACK_WARN(msg)   fprintf(stdout, "WARNING(%s|%s): %s\n", __FILE__, __FUNCTION__, msg)



/***   stack_new   ************************************************************
 *
 */
stack_t *
stack_new(size_t elem_size, void (*free_func)(void *data))
{
    stack_t *s;
    
    if ((s=(stack_t *)calloc(1, sizeof(stack_t)))==NULL)
    {
        STACK_ALLOC_ERROR;
        return NULL;
    }
    
    s->elem_size = elem_size;
    s->size = 0;
    s->free_func = free_func;
     
    return s;
}
/***   stack_new   ***/



/***   stack_free   ***********************************************************
 *
 */
int
stack_free(stack_t *s)
{
    if (s==NULL)
    {
        STACK_WARN("stack is NULL");
        return -1;
    }
    while (stack_size(s))
        stack_pop(s, 1);
    free(s);
    s = NULL;
    return 0;
}
/***   stack_free   ***/



/***   stack_push   ***********************************************************
 *
 */
int
stack_push(stack_t *s, void *data)
{
    if (s==NULL)
    {
        STACK_WARN("stack is NULL");
        return -1;
    }
    if (data==NULL)
    {
        STACK_WARN("data is NULL");
        return -1;
    }

    if ((s->elems = (void**) realloc(s->elems, s->elem_size * (s->size+1)))==NULL)
    {
        STACK_ALLOC_ERROR;
        return 0;    
    }
    s->elems[s->size] = data;
    s->size++;
    
    return 0;
}
/***   stack_push   ***/



/***   stack_pop   ************************************************************
 *
 */
void *
stack_pop(stack_t *s, int free)
{
    void *elem;
    
    if (s==NULL)
    {
        STACK_WARN("stack is NULL");
        return NULL;
    }    
    if (stack_size(s)==0)
    {
        STACK_WARN("stack is empty");
        return NULL;
    }
    
    

    s->size--;
    elem = s->elems[s->size];
    if (free)
    {
        s->free_func(elem);
        return NULL;
    }
    else
    {
        return elem;
    }
}
/***   stack_pop   ***/



/***   stack_peek   ***********************************************************
 *
 */
void *
stack_peek(stack_t *s)
{
    if (s==NULL)
    {
        STACK_WARN("stack is NULL");
        return NULL;
    }
    if (stack_size(s)==0)
    {
        STACK_WARN("stack is empty");
        return NULL;
    }
    
    return (void*) s->elems[s->size-1];
}
/***   stack_peek   ***/



/***   stack_size   ***********************************************************
 *
 */
int
stack_size(stack_t *s)
{
    if (s==NULL)
    {
        STACK_WARN("stack is NULL");
        return 0;
    }
    return s->size;
}
/***   stack_size   ***/


