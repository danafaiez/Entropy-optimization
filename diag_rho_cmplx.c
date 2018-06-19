#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "clapack.h"


char jobzc[]="N";
char uploc[] = "U";

int diag_rho_cmplx(int n, double * matrix_double, double * eigenvalues)
{
   integer N = n;
   doublecomplex * M = (doublecomplex *) matrix_double;
   doublecomplex * z__;//not referenced because we don't want evectors
   integer ldz=1;
   doublecomplex * work;
   //integer lwork = N+1;
   integer lwork = -1;
   integer lrwork = N+1;
   double * rwork;
   integer * iwork;
   integer liwork =1;
   integer info;

   calloc_(work, 1);
   calloc_(rwork, 1);
   calloc_(iwork, 1);

   zhpevd_(jobzc, uploc, &N, M, eigenvalues, z__, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
   lwork = work[0].r;
   lrwork = rwork[0];
   liwork = iwork[0];
   free(work);
   free(rwork);
   free(iwork);
   calloc_(work, lwork);
   calloc_(rwork, lrwork);
   calloc_(iwork, liwork);

   zhpevd_(jobzc, uploc, &N, M, eigenvalues, z__, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
   free(work);
   free(rwork);
   free(iwork);

   return (int) info;

}

