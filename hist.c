#include "hist.h"


HIST * new_hist(float x_min, float x_max, int num_bins)
{
   HIST * h;

   calloc_(h,1);
   h->dx = (x_max-x_min)/num_bins;
   h->min = floor(x_min/h->dx);
   h->max = ceil(x_max/h->dx);
   h->num_bins = h->max- h->min;
   h->x_min = h->min * h->dx;
   h->x_max = h->max * h->dx;
   calloc_(h->p,num_bins+1);

   return h;
}

void delete_hist(HIST * h)
{
   free(h->p);
   free(h);
}

void renew_hist(HIST * h, float new_x_min, float new_x_max)
{
   float * new_p;
   int new_num_bins = h->num_bins;
   int new_min, new_max;
   

   new_x_min = MN(h->x_min, new_x_min);
   new_x_max = MX(h->x_max, new_x_max);

   new_min = floor(new_x_min/h->dx);
   new_max = ceil(new_x_max/h->dx);

   new_num_bins = new_max - new_min;

   calloc_(new_p, new_num_bins+1);

   memmove_(new_p + h->min - new_min, h->p, h->num_bins);

   
   free(h->p);

   h->x_max = new_max * h->dx;//make max and min commensurate with dx, still a bit sloppy
   h->x_min = new_min * h->dx;
   h->max = new_max;
   h->min = new_min;
   h->num_bins = new_num_bins;
   h->p = new_p;
}

void  add_hist_weighted(HIST * h, float x, float weight)
{
   int bin;
   if (x < h->x_min)
   {
      renew_hist(h,x, h->x_max);
   }
   else if (x > h->x_max)
   {
      renew_hist(h, h->x_min, x);
   }

   bin = (int)(x/h->dx - h->min);
   assert(bin >= 0);
   assert(bin < h->num_bins+1);
   h->p[(int)(x/h->dx - h->min)] += weight;
   h->total_weight += weight;
   h->num_entries++;
}


void  add_hist(HIST * h, float x)
{
   add_hist_weighted(h,x,1.0);
}

float lookup(HIST * h, float x)
{
      int bin = (int)(x/h->dx - h->min);
      if (bin < 0 || bin > h->num_bins)
	 return 0;
      return h->p[bin];
}

void cumulate_hist(HIST * h_tot, HIST * h1, HIST * h2)
{
   int i;
   float new_min = MN(h1->min, h2->min);
   float new_max = MX(h1->max, h2->max);

   assert(h1->dx == h2->dx);

   for(i = new_min; i < new_max; i++)
   {
      float x = h1->dx*i;
      add_hist_weighted(h_tot, x, lookup(h1,x));
      add_hist_weighted(h_tot, x, lookup(h2,x));
   }
}

void print_hist(HIST * h, FILE * f)
{
   int i;
   for(i=0; i < h->num_bins; i++)
   {
      fprintf(f,"%f %f\n", (i+h->min)*h->dx, h->p[i]);
   }
}

#if 0//test
int main(int argc, char ** argv)
{
   int i;
   HIST * h = new_hist(-0.2, 0.2, 100);

   for(i=0; i < 1000; i++)
   {
      add_hist(h, ran()-0.5);
   }
   for(i=0; i < h->num_bins; i++)
   {
      printf("%f %f\n", (i+h->min)*h->dx, h->p[i]);
   }
   return 0;

}

#endif
