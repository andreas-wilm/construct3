/******************************************************************************
* 
* stack.h - generic stack routines
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
 *  CVS $Id: stack.h,v 1.3 2004-05-25 13:29:18 wilm Exp $    
 */



typedef struct {
    void   **elems;
    size_t   elem_size;
    int      size;
    void  (*free_func)(void *data);
} stack_t;



/*** stack_new:
 *   return new stack, or NULL on error
 */
extern stack_t *
stack_new(size_t elemsize, void (*free_func)(void *data));

/*** stack_free:
 *   free a stack
 *   return 0 on success, -1 otherwise
 */
extern int
stack_free(stack_t *s);

/*** stack_push:
 *   push data onto stack
 *   return 0 on success, -1 otherwise
 */
extern int
stack_push(stack_t *s, void *data);


/*** stack_pop:
 *   remove top element
 *   if free is 1, free element and return NULL
 *   otherwise return pointer to element or NULL on error
 */
extern void *
stack_pop(stack_t *s, int free);


/*** stack_peek:
 *   return pointer to top element or NULL on error
 */
extern void *
stack_peek(stack_t *s);

/*** stack_size:
 *   returns stack size
 */
extern int
stack_size(stack_t *s);
