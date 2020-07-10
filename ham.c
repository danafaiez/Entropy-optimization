
#include <assert.h>
#include "defs.h"
#include "data_structures.h"
#include "basis_states.h"
#include "fock.h"
#include "funcs.h"
#include "io.h"
#include "ham.h"

double * energy(ull * states, PARAMS * pm)
{
   int i, L = pm->L;
   double *t = pm->t;
   double U = pm->U;
   double Up = pm->Up;
   int nn[2][L][2];
   int numstates = size_(states);

   double *h;
   newarr_(h,(numstates*numstates));
   HTABLE * hash = create_hash(states);
   int bath = pm->num_bath_sites;
   int neigh;
   for(neigh=0;neigh < 2; neigh++)
   {
      for(i=0;i < L;i++)
      {
         if (0)
         {
            nn[neigh][i][RIGHT] = (i+(neigh+1))%L;
            nn[neigh][i][LEFT] = (i-(neigh+1)+L)%L;
         }
         else
         {
            nn[neigh][i][RIGHT] = i+(neigh+1);
            nn[neigh][i][LEFT] = i-(neigh+1);
         }
      }
   }

   //Kinetic:

   for(neigh=0;neigh < 2; neigh++)
   {
      int bra_ind, ket_ind;
      //      for(bra_ind=0; bra_ind < numstates;bra_ind++)
      {
         for(ket_ind=0; ket_ind < numstates;ket_ind++)
         {
            int m,n;
            for(m=0;m < L;m++)
            {
               for(n=0;n < 2;n++)
               {
                  int m_neigh = nn[neigh][m][n];
                  if (m_neigh >= 0 && m_neigh < L)
                  {
                     
//                 if ((m>=0&&m<(L-bath)&&((nn[neigh][m][n])>=0)&&((nn[neigh][m][n])<(L-bath))) || (m>=(L-bath)&&m<(L)&&((nn[neigh][m][n])>=(L-bath))&&((nn[neigh][m][n])<(L))) )                    
  //                  {
                     int sign = 1;
                     ull ket = states[ket_ind];
                     ull ket_anh = a(ket,m_neigh,L,&sign);
                     ull ket_create_anh = adag(ket_anh,m,L,&sign);

                     int * indx = _htable_get(hash,&ket_create_anh);
                     if (!indx)
                        continue;
                     assert(ket_create_anh == states[indx[0]]);
                     H(indx[0],ket_ind) += -t[neigh]*sign;
    //               } 
                 }
               }
            }
         }
      }
   }
   //Potential
   int ket_ind;
   for(ket_ind=0; ket_ind < numstates;ket_ind++)
   {
      ull ket = states[ket_ind];
      int m,j;
      for(m=0;m < L;m++)
      {
         if (!bit_val(ket,m))
         {
            continue;
         }
         for(j=-1; j < 2; j += 2)
         {
            int neigh = m+ j;
            if (neigh >= L || neigh < 0)
               continue;
            if (bit_val(ket,neigh))
               H(ket_ind,ket_ind) += U/2;
         }
         for(j=-1; j < 2; j += 2)
         {
            int neigh = m+ 2*j;
            if (neigh >= L || neigh < 0)
               continue;
            if (bit_val(ket,neigh))
               H(ket_ind,ket_ind) += Up/2;
         }


      }
   }
   htable_free(hash);
   return h;
}

int check_hermitian(double * h,int numstates)
{
  int i,j;
  for(i=0; i < numstates;i++)
  {
     for(j=0; j < i;j++)
     {
        assert(fabs(H(i,j) - H(j,i)) < 1e-6);
     }
  }
  return 1;
}
   
double * diag(double * h, int  numstates, PARAMS * pm)
{
   char JOBZ,UPLO;
   integer M=numstates;
   integer LDA=M;
   double * ener;
   integer LWORK=2*M*M+6*M+1;
   double * WORK;
   integer INFO;
   newarr_(WORK,LWORK);
   newarr_(ener,numstates);

   JOBZ='V'; //V for eigenvectors also. N for getting only eigenvalues.
   UPLO='U';

   dsyev_(&JOBZ,&UPLO,&M,h,&LDA,ener,WORK,&LWORK,&INFO);

   if(INFO!=0)
   {
      printf("INFO=%d. Diagonalization issues.\n",(int)INFO);
   }

   printf("\nDiagonalization completed.\n");  
   //h contains the eigenvectors and ener contains the eigenvalues.
   freearr_(WORK);
   return ener;
}

void write_efiles(PARAMS * pm, double * h, double * ener)
{
   char vector[100];
   FILE *data;
   sprintf(vector, "ener_t%1.1f_tprime%1.2f_V%1.2f_Vp%1.2f_L%d_Ne%d.dat", pm->t[0], pm->t[1],pm->U,pm->Up,pm->L,pm->num_particles);
   data=fopen(vector, "w");
   int i;
   for(i=0;i<pm->numstates;i++)
   {
      fprintf(data, "%lf\n", ener[i]);
   }

   fclose(data);

   sprintf(vector, "ener_t%1.1f_tprime%1.2f_V%1.2f_Vp%1.2f_L%d_Ne%d.bin", pm->t[0], pm->t[1],pm->U,pm->Up,pm->L,pm->num_particles);
   fclose(write_energies(0, ener, pm, vector, 0));

   char vec2[100];
   sprintf(vec2, "coeff_t%1.1f_tprime%1.2f_V%1.2f_Vp%1.2f_N%d_Ne%d.dat", pm->t[0], pm->t[1],pm->U,pm->Up,pm->L,pm->num_particles);
   {
      //int tot_evs_processed=0;
      FILE * evector_out=0;
      sprintf(vec2, "real_coeff_t%1.1f_tprime%1.2f_V%1.2f_Vp%1.2f_N%d_Ne%d.bin", pm->t[0], pm->t[1],pm->U,pm->Up,pm->L,pm->num_particles);

      evector_out = write_all_evectors(evector_out, h, pm, vec2, 0);//directory arg = 0 here
      fclose(evector_out);
   }

}

