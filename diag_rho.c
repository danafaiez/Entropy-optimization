#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "clapack.h"


char jobz[]="N";
char uplo[] = "U";

int diag_rho(int n, double * matrix_double, double * eigenvalues)
{
   integer N = n;
   double * M = (double *) matrix_double;
   double * z__;//not referenced because we don't want evectors
   integer ldz=1;
   double * work;
   //integer lwork = N+1;
   integer lwork = -1;
//   integer lrwork = N+1;
//   double * rwork;
   integer * iwork;
   integer liwork =1;
   integer info;

   calloc_(work, 1);
//   calloc_(rwork, 1);
   calloc_(iwork, 1);

   dspevd_(jobz, uplo, &N, M, eigenvalues, z__, &ldz, work, &lwork, iwork, &liwork, &info);
   lwork = work[0];
//   lrwork = rwork[0];
   liwork = iwork[0];
   free(work);
//   free(rwork);
   free(iwork);
   calloc_(work, lwork);
//   calloc_(rwork, lrwork);
   calloc_(iwork, liwork);

   dspevd_(jobz, uplo, &N, M, eigenvalues, z__, &ldz, work, &lwork, iwork, &liwork, &info);
   free(work);
//   free(rwork);
   free(iwork);

   return (int) info;

}

