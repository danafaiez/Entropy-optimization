#include "io.h"
#include <sys/stat.h>
#include <string.h>
#include "defs.h"
#include "fock.h"

const char  magic_string_real[] = "EVECTOR_REAL_ONLY_DATA";
const char  magic_string_complex[] = "EVECTOR_COMPLEX_ONLY_DATA";

FILE * open_dir_file(char * dirname, char * filename, char * r_w)
{
   if (dirname)
   {
      char name[1024];
      mkdir(dirname,S_IRWXU);
      sprintf(name,"%s/%s",dirname,filename);
      return fopen(name,r_w);
   }
   else
   {
      return fopen(filename,r_w);
   }
}

#define _fwrite(p,num) fwrite((p), num*sizeof((p)[0]),1,evector_out)

FILE * write_evector(FILE * evector_out, double * evector, PARAMS * pm, char * evector_out_name, char * dirname)
{
   int i;
   int n = pm->numstates;

   if (evector == 0)
   {
      fprintf(stdout,"in write_evector:  evector is a null pointer. can't write to %s\n", evector_out_name);
      exit(4);
   }



   if (evector_out == 0)
   {
      evector_out = open_dir_file(dirname, evector_out_name,"w");

      fputs(magic_string_real,evector_out);
      _fwrite(pm,1);
   }

   _fwrite(evector,n);
   return evector_out;
}


FILE * write_all_evectors(FILE * evector_out, double * evectors, PARAMS * pm, char * evector_out_name, char * dirname)
{
   int i;
   int n = pm->numstates;

   if (evectors == 0)
   {
      fprintf(stdout,"in write_evector:  evectors is a null pointer. can't write to %s\n", evector_out_name);
      exit(4);
   }



   if (evector_out == 0)
   {
      evector_out = open_dir_file(dirname, evector_out_name,"w");

      fputs(magic_string_real,evector_out);
      _fwrite(pm,1);
   }

   _fwrite(evectors,(n*n));
   return evector_out;
}
#undef _fwrite

#define _fwrite(p,num) fwrite((p), num*sizeof((p)[0]),1,energy_out)
FILE * write_energies(FILE * energy_out, double * energies, PARAMS * pm, char * energy_out_name, char * dirname)
{
   int i;
   int n = pm->numstates;

   if (energies == 0)
   {
      fprintf(stdout,"in write_energies:  energies is a null pointer. can't write to %s\n", energy_out_name);
      exit(4);
   }



   if (energy_out == 0)
   {
      energy_out = open_dir_file(dirname, energy_out_name,"w");

      fputs(magic_string_real,energy_out);
      _fwrite(pm,1);
   }

   _fwrite(energies,(n));
   return energy_out;
}
#undef _fwrite

#define _fread(p,length) fread((p),length*sizeof((p)[0]),1,evector_in)

double * read_evectors(char * evector_in_name, PARAMS * pm)
{
   int i,n,nr;
   char *magic;
   FILE * evector_in = fopen(evector_in_name,"r");
   int magic_len = strlen(magic_string_real);
   double * evector;

   calloc_(magic,magic_len);

   if (fread(magic, sizeof(magic[0]), magic_len, evector_in) != magic_len
	 || strncmp(magic_string_real, magic, magic_len))
   {
      fprintf(stderr,"can't read %s, doesn't seem to be the right type. exiting...\n", evector_in_name);
      exit(4);
   }

   _fread(pm,1);
   
   n = pm->numstates;
   nr = pm->rbasis;

   calloc_(evector, (n*nr));
   _fread(evector, (n*nr));

   fclose(evector_in);

   return evector;
}

FILE  * read_header(char * evector_in_name, PARAMS * pm)
{
   int i,n,nr;
   char *magic;
   FILE * evector_in = fopen(evector_in_name,"r");
   int magic_len = strlen(magic_string_real);

   calloc_(magic,magic_len);

   if (fread(magic, sizeof(magic[0]), magic_len, evector_in) != magic_len
	 || strncmp(magic_string_real, magic, magic_len))
   {
      fprintf(stderr,"can't read %s, doesn't seem to be the right type. exiting...\n", evector_in_name);
      exit(4);
   }

   _fread(pm,1);
   
   free(magic);

   return evector_in;
}

double * read_evector(FILE * evector_in, PARAMS * pm)
{
   int i,n,nr;
   double * evector;

   n = pm->numstates;
   nr = pm->rbasis;

   newarr_(evector, n);
   _fread(evector, n);


   return evector;
}

#undef _fread
#define _fread(p,length) fread((p),length*sizeof((p)[0]),1,energy_in)

double * read_energies(char * energy_in_name, PARAMS * pm)
{
   int i,n,nr;
   char *magic;
   FILE * energy_in = fopen(energy_in_name,"r");
   int magic_len = strlen(magic_string_real);
   double * energies;

   calloc_(magic,magic_len);

   if (fread(magic, sizeof(magic[0]), magic_len, energy_in) != magic_len
	 || strncmp(magic_string_real, magic, magic_len))
   {
      fprintf(stderr,"can't read %s, doesn't seem to be the right type. exiting...\n", energy_in_name);
      exit(4);
   }

   _fread(pm,1);
   
   n = pm->numstates;
   nr = pm->rbasis;

   newarr_(energies, n);
   _fread(energies, n);

   fclose(energy_in);
   free(magic);
   return energies;
}


#undef _fread

void write_ascii(double * arr,char * file_name, char * dir_name )
{
   int n = size_(arr);
   int i;
   FILE * fl = open_dir_file(dir_name, file_name, "w");
   for(i=0; i < n; i++)
   {
      fprintf(fl,"%lf\n",arr[i]);
   }

   fclose(fl);
}

