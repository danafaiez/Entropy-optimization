#include "math.h"
#include "defs.h"
#include "basis_states.h"
#include <time.h>

double gaussian_random()
{

   double rsq;
   double x,y;
   do 
   {
      x = 2*ran()-1;
      y = 2*ran()-1;
      rsq = x*x+y*y;
   } 
   while (rsq >= 1.0 || rsq == 0);

   return x * sqrt(-2*log(rsq)/rsq);
}

complx c_gaussian_random()
{
srand(time(0));
   double x,y;
   double rsq;
   do 
   {
      x = 2*ran()-1;
      y = 2*ran()-1;
      rsq = x*x+y*y;
   } 
   while (rsq >= 1.0 || rsq == 0);

   return (x+y*1.0j) * sqrt(-2*log(rsq)/rsq);
}

