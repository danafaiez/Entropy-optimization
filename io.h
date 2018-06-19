#ifndef IO_H
#define IO_H
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "data_structures.h"


FILE * open_dir_file(char * dirname, char * filename, char * r_w);
FILE * write_evector(FILE * evector_out, double * evector, PARAMS * pm, char * evector_out_name, char * dirname);
FILE * write_all_evectors(FILE * evector_out, double * evectors, PARAMS * pm, char * evector_out_name, char * dirname);
double * read_evectors(char * evector_in_name, PARAMS * pm);
FILE  * read_header(char * evector_in_name, PARAMS * pm);
double * read_evector(FILE * evector_in, PARAMS * pm);
void write_ascii(double * arr,char * file_name, char * dir_name );

#endif// IO_H
