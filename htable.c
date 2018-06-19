
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


#include "defs.h"
#include "htable.h"


unsigned int hash_function (HTABLE * htable, uchar * key)
//algorithm found on the web, seems pretty common, see http://burtleburtle.net/bob/hash/doobs.html
{
   int i;
   uint h;
   for (h=0, i=0; i < htable->key_length; i++) 
   {
      h += key[i];  
      h += (h<<10); 
      h ^= (h>>6); 
   } 
   h += (h<<3); 
   h ^= (h>>11); 
   h += (h<<15);
   return h % (htable->tot_array_size);
}


HTABLE * htable_create(unsigned int tot_array_size,int  key_length)
{
   HTABLE *ht;
   calloc_ (ht,1);
   calloc_ (ht->array,tot_array_size);
   calloc_ (ht->array_length,tot_array_size);
   ht->tot_array_size = tot_array_size;
   ht->key_length = key_length;
   return ht;
}


void htable_free (HTABLE * ht)
{
   int i;
   for(i=0; i < ht->tot_array_size; i++)
   {
      int j, length = ht->array_length[i];
      KEY_VALUE * elt = ht->array[i];
      for(j=0; j < length; j++)
      {
	 free (elt[j].value);
	 free (elt[j].key);
      }

      if (length > 0) free(ht->array[i]);
   }
   free(ht->array);
   free(ht->array_length);
   free(ht);
}



int * _htable_get (HTABLE * ht, void * key)
{
   int index = hash_function(ht,key);
   KEY_VALUE * elt = ht->array[index];
   int i, length = ht->array_length[index];

   if (elt)
   for(i=0; i < length; i++)
   {
      if (!memcmp (key, elt[i].key, ht->key_length)) 
      {
	 return elt[i].value;
      }
   }
   
   return 0;
}


int * htable_get (HTABLE * ht, int key)
{
    return _htable_get (ht, &key);
}



void _htable_set (HTABLE * ht, void * key, int value)
{
   int index = hash_function(ht,key);
   KEY_VALUE * elt = ht->array[index];
   int i, length = ht->array_length[index];

   
   i = length;

   if (elt)
   for(i=0; i < length; i++)
   {
      if (!memcmp (key, elt[i].key, ht->key_length)) {
	 elt[i].value[0] = value;
	 break;
      }
   }
   
   if (i == length) {
      realloc_ (ht->array[index], ++length);
      elt = ht->array[index];
      ht->array_length[index] = length;
      calloc_ (elt[length-1].value, 1);
      elt[length-1].value[0] = value;
      calloc_ (elt[length-1].key, ht->key_length);
      memcpy(elt[length-1].key, key, sizeof(elt[0].key[0]) * ht->key_length);
   }
}

inline void htable_set (HTABLE * ht, int  key, int value)
{
    _htable_set (ht,  &key, value);
}

int total_num_elts (HTABLE * ht)
{
   int index, total_num=0;
   for(index=0; index < ht->tot_array_size; index++)
   {
	 total_num += ht->array_length[index];
   }
   return total_num;
}

//#define TEST

#ifdef TEST
#include <stdio.h>

#define MEM_DEBUG

void string_test (int length)
{
    int i;
    unsigned char * key;
    HTABLE * ht =  htable_create(length/4 , 4);

   for (i=1; i < length/2; i++)
   {
      key = (uchar *) &i;
      _htable_set(ht, key, i);
   }
 
   for (i=1; i < length; i++)
   {
      int * value;
      key = (uchar *) &i;
      if (value = _htable_get(ht, key))
            printf("key = %d, value = %d\n",i, value[0]);
   }

   htable_free(ht);
}

void int_test (int length)
{
    int i;
    unsigned char * key;
    HTABLE * ht =  htable_create(length/4 , 4);

   for (i=1; i < 100*length/2; i += 100)
   {
      htable_set(ht, i, i*1000);
   }
 
   for (i=1; i < 100*length; i += 100)
   {
      int  *value;
      if (value = htable_get(ht, i))
            printf("key = %d, value = %d\n",i, value[0]);
   }

   htable_free(ht);
}

main ()
{
//   string_test(256);

#ifdef MEM_DEBUG
   mtrace();
#endif
   int_test(256);
}

#endif


