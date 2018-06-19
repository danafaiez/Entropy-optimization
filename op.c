#include "fock.h"

//create: create_destroy=0, destroy: create_destroy=1
ull a_adag(int create_destroy, ull config, int pos, int length, int *sign)
{
   if (sign[0]==0)
     return 0;
   if (create_destroy^bit_val(config,pos))
   {
     *sign = 0;
     return 0;
   }
   flip_bit(config,pos);
   //__builtin_popcountll() counts the number of ones
   *sign *=  __builtin_popcountll(config&mask(pos))&1?-1:1;
   return config;
}

ull adag(ull config, int pos, int length, int *sign)
{
    return a_adag(0, config, pos, length, sign);
}

ull a(ull config, int pos, int length, int *sign)
{
    return a_adag(1, config, pos, length, sign);
}

