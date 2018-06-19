
/*
 * Software License
 *
 * Copyright (c) 2001 Joshua M. Deutsch. All rights 
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution, if
 *    any, must include the following acknowlegement:  
 *       "This product includes software developed Joshua M. Deutsch
 *        Dept of Physics University of California, Santa Cruz"
 *    Alternately, this acknowlegement may appear in the software itself,
 *    if and wherever such third-party acknowlegements normally appear.
 *
 * 4. The name "GESSES"
 *    must not be used to endorse or promote products derived
 *    from this software without prior written permission. For written 
 *    permission, please contact josh@physics.ucsc.edu.
 *
 * 5. Products derived from this software may not be called "gesses"
 *    nor may "gesses" appear in their names without prior written
 *    permission of Joshua M. Deutsch.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL JOSHUA M. DEUTSCH OR THE UNIVERSITY
 * OF CALIFORNIA BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 */

#ifndef HTABLE_H
#define HTABLE_H

typedef struct {
   void * key;
   int * value;
} KEY_VALUE;


typedef struct {
   KEY_VALUE ** array;
   int * array_length;
   int tot_array_size;
   int key_length;
}  HTABLE;



void htable_free (HTABLE * ht);
HTABLE * htable_create(unsigned int tot_array_size,int  key_length);
int * htable_get (HTABLE * ht, int key);
void _htable_set (HTABLE * ht, void * key, int value);
int * _htable_get (HTABLE * ht, void * key);
void htable_set (HTABLE * ht, int  key, int value);
int total_num_elts (HTABLE * ht);
#endif // HTABLE_H
