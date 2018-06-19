
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
//#define TEST_OVERLAP 
//#define TEST_SORT 
//#define TEST_SHRINK 

int sizearr(void * arr)
{
   return size_(arr);
}

int all_ones(int * arr)
{
   int i;
   int true = 1;
   for(i=0; i < size_(arr);i++)
      true = true & arr[i];
   return true;
}

int all_zeros(int * arr)
{
   int i;
   int true = 1;
   for(i=0; i < size_(arr);i++)
      true = true & !(arr[i]);
   return true;
}

int * d2i_arr(double * arr)
{
   int i;
   int * iarr=0;
   newarr_(iarr,size_(arr));
   for(i=0; i < size_(arr);i++)
      iarr[i] = (int) arr[i];

   return iarr;
}

double avearr(double * arr)
{
   double ave=0;
   avearr_(arr,ave);
   return ave;
}

double sumarr(double * arr)
{
   int i;
   double sum=0;
   for(i=0; i < size_(arr);i++)
      sum += arr[i];
   return sum;
}

double sum2darr(double ** arr)
{
   int i;
   double sum=0;
   for(i=0;i < size_(arr);i++)
      sum += sumarr(arr[i]);
   return sum;
}

void printarr_(FILE * file, double * arr)
{
   int i;
   fprintf(file, "%d:\n",i);
   for(i=0;i < size_(arr);i++)
      fprintf(file, "  %lf\n",arr[i]);
}
void printarr(double * arr)
{
  printarr_(stdout, arr);
}
 
void print2darr_(FILE * file, double ** arr)
{
   int i;
   for(i=0;i < size_(arr);i++)
      fprintf(file, "%d:\n",i);
      printarr_(file,arr[i]);
}
 
void print2darr(double ** arr)
{
  print2darr_(stdout,  arr);
}
 
void reverse_int_array(int * a, int length)
{
   int i = 0;
   int j = length-1;

   if (length <= 1) return;

   while (i < j)
   {
      SWAP(a[i],a[j]);
      i++; j--;
   }
}

void * concat(void * a1, size_t n1, void * a2, size_t n2)

{
   void * result = 0;

   calloc_(result, n1+n2);
   memcpy(result, a1, n1);

   if (n2)
      memcpy(result+n1,a2,n2);

   return result;
}


//return 1 if overlap, 0 otherwise
int overlap(void * array_1,int n1,void *  array_2,int  n2, size_t size,
      int (*array_comp)(const void *, const void * ))
{
   void * concat_array;

   if (n2 == 1) // for efficiency handle these two cases separately
   {
      int i;
      for(i=0; i < n1*size;i += size)
      {
	 if (!array_comp(array_2,array_1+i)) 
	    return 1;
      }
      return 0;
   }

   if (n1 == 1)
   {
      int i;
      for(i=0; i < n2*size;i += size)
      {
	 if (!array_comp(array_1,array_2+i)) 
	    return 1;
      }
      return 0;
   }

   concat_array = concat(array_1,n1*size,array_2,n2*size);
   qsort(concat_array, n1+n2, size * sizeof(*concat_array), array_comp);
   {
      int i;
      for(i=0; i < (n1+n2-1)*size;i += size)
      {
	 if (!array_comp(concat_array+i,concat_array+i+size)) 
	 {
	    free(concat_array);
	    return 1;
	 }
      }
    free(concat_array);
    return 0;

   }
}

int * shrink_arr_int(int *arry, int indx)//remove elt at indx
{
   int * shrunk_array;
   int n = size_(arry);
   int size = sizeof(arry[0]);
   newarr_(shrunk_array, (n-1));

   if (indx > 0) 
      memcpy(shrunk_array, arry, size*indx);

   if (indx < (n-1))
      memcpy(shrunk_array + size * indx, arry + size * (indx+1), size* (n-1-indx));

   return shrunk_array;
}


void * shrinkarr(void *arry, int indx, size_t size)//remove elt at indx
{
  int n = size_(arry);
  void * shrunk_array;
  newarr_(shrunk_array, size*(n-1));
  size_(shrunk_array) = n-1;

   if (indx > 0) 
      memcpy(shrunk_array, arry, size*indx);

   if (indx < (n-1))
      memcpy(shrunk_array + size * indx, arry + size * (indx+1), size* (n-1-indx));

   return shrunk_array;
}



void * shrink(void *arry, int n, int indx, size_t size)//remove elt at indx
{
  void * shrunk_array;
   calloc_(shrunk_array, size*(n-1));
   //shrunk_array = (void *) calloc(size*(n-1),sizeof(void));

   if (indx > 0) 
      memcpy(shrunk_array, arry, size*indx);

   if (indx < (n-1))
      memcpy(shrunk_array + size * indx, arry + size * (indx+1), size* (n-1-indx));

   return shrunk_array;
}

int * shrink_int(int * arry, int n, int indx) 
{
   return (int *) shrink(arry, n, indx, sizeof(arry[0]));
}

void * consolidate(void * array_1,int n1,void *  array_2,int  n2, size_t size,
      int (*array_comp)(const void *, const void * ), int * output_length)
{
   void * concat_array;
   void * consolidated_array=0;


   calloc_(consolidated_array, size*(n1+n2));

   concat_array = concat(array_1,n1*size,array_2,n2*size);
   qsort(concat_array, n1+n2, size * sizeof(*concat_array), array_comp);
   {
      int i,j=0;
      for(i=0; i < (n1+n2-1)*size;i += size)
      {
	 if (array_comp(concat_array+i,concat_array+i+size)) 
	 {
	    memcpy(consolidated_array + j*size, concat_array + i, size);
	    j++;
	 }
      }
      memcpy(consolidated_array + j*size, concat_array + i, size);
      j++;
      realloc_(consolidated_array, j*size);
      free(concat_array);
      output_length[0] = j;
      return consolidated_array;

   }
}


int int_comp(const void * i, const void * j)
{
   return (((int *) i)[0] - ((int *) j)[0]);
}

int * consolidate_int(int * array_1,int n1,int *  array_2,int  n2, int * output_length)
{
   return (int *) consolidate(array_1,n1,array_2,n2, sizeof(array_1[0]), int_comp, output_length); 
}

int overlap_int(int * array_1,int n_1, int *  array_2, int  n2) 
{
  return overlap( array_1,n_1,  array_2,n2, sizeof(array_1[0]), int_comp);
}
#ifdef TEST_OVERLAP
int  main()
{
   {
      int a[]={1,5,7};
      int b[]={2,3,8};
      int c[]={2,3,7};
      int right_ab[] = {1,2,3,5,7,8};
      int right_bc[] = {2,3,7,8};
      int length_ab=0;
      int length_bc=0;
      int *  consolidated_array_ab = consolidate_int(a,3,b,3, &length_ab);
      int *  consolidated_array_bc= consolidate_int(b,3,c,3, &length_bc);

      if (length_ab != 6 || memcmp_(consolidated_array_ab, right_ab, 6)) 
      {
	 int i;
	 printf("error, considation of a and b failed\n");
	 printf("length =  %d, should be 6\n", length_ab);
      }
      else printf("passed consolidation of a and b\n");

      if (length_bc != 4 || memcmp_(consolidated_array_bc, right_bc, 4)) 
      {
	 int i;
	 printf("error, considation of b and c failed\n");
	 printf("length =  %d, should be 4\n", length_bc);
      }
      else printf("passed consolidation of b and c\n");


      if (overlap_int(a,3,b,3)) printf("error, a and b don't overlap\n");
      else if (!overlap_int(a,3,c,3)) printf("error, a and c do overlap\n");
      else printf("passed the test 1st test with small array\n");


   }

   {
      int * a = 0;
      int * b = 0;
      int i;
#define TEST_SIZE 10000
      calloc_(a,TEST_SIZE);
      calloc_(b,TEST_SIZE);
      for(i=0;i < TEST_SIZE;i++)
      {
	 a[i] = 2*i;
	 b[i] = 2*i+1;
      }

      if (overlap_int(a,TEST_SIZE,b,TEST_SIZE)) 
      {
	 printf("error, a and b don't overlap\n");
	 exit(2);
      }
      b[TEST_SIZE-1] = 0;
      if (!overlap_int(a,TEST_SIZE,b,TEST_SIZE))
      {
	 printf("error, a and b do overlap\n");
	 exit(2);
      }

      else printf("passed the test 2nd test with array of %d elements\n",TEST_SIZE);

   }
   
}


#endif

void  partial_sort(void * array,size_t n,size_t size, int (* comp)(const void *, const void *), int num_sorted)
{
   int i,j, min_index=0;
   void * tmp; 
   alloca_(tmp, n);

   for (i=0; i < num_sorted;i++)
   {
      for(j=i;j < n; j++)
      {
	 if (comp(array + size*j, array + size*i) < 0) { min_index = j; }
      }
      //swap i and min_index values
      memcpy(tmp, array +  size*min_index, size);
      memcpy(array +  size*min_index, array +  size*i, size);
      memcpy(array +  size*i, tmp, size);
   }
}
 
#ifdef TEST_SORT
int  main()
{
#define length 100
   int a[length];
   int i;

   for(i=0; i < length; i++)
      a[i] = length -1 - i;
   for(i=0; i < 3; i++)
      printf("a[%d] = %d\n",i,a[i]);

   partial_sort(a,length, sizeof(a[0]), int_comp,3);
   for(i=0; i < 3; i++)
      printf("a[%d] = %d\n",i,a[i]);

   qsort(a,length, sizeof(a[0]), int_comp);
   for(i=0; i < 3; i++)
      printf("a[%d] = %d\n",i,a[i]);
}

#endif

#ifdef TEST_SHRINK
int  main()
{
#define length 10
   int a[length];
   int i;
   int * arry;

   for(i=0; i < length; i++)
      a[i] = i;


   for(i=1; i < length-1; i++)
   {
      int tmp;
      arry = shrink_int(a,length,i);
      if ((tmp = arry[i]-arry[i-1]) != 2) 
      {
	 printf("at i=%d: error in shrink_int, %d \n",i,tmp);
	 exit(2);
      }
      free(arry);
   }
}
#endif
