//#define REAL_IS_FLOAT
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "defs.h"


typedef struct {
   float dx, x_max, x_min;
   int num_entries;
   float total_weight;
   int min, max, num_bins;
   float * p;
} HIST;



HIST * new_hist(float x_min, float x_max, int num_bins);
void  add_hist(HIST * h, float x);
void cumulate_hist(HIST * h_tot, HIST * h1, HIST * h2);
void  add_hist_weighted(HIST * h, float x, float weight);
void print_hist(HIST * h, FILE * f);
