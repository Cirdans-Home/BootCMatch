/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra    IAC-CNR, IT
                       Panayot S. Vassilevski Portland State University, OR USA

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the BootCMatch group or the names of its contributors may
       not be used to endorse or promote products derived from this
       software without specific written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BootCMatch GROUP OR ITS CONTRIBUTORS
  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "bcm.h"
#include "ioutil.h"

/* #define DUMP_HIER */
void dump_on_file(const char *prefix, bcm_BootAMG *boot_amg)
{
  int j,k;
  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);
  bcm_CSRMatrix    **A_array;
  bcm_CSRMatrix    **P_array;
  bcm_CSRMatrix    **L_array;
  bcm_CSRMatrix    **U_array;
  bcm_Vector       **D_array;

  for (k=0; k<1; k++) {
    A_array           = bcm_AMGHierarchyAArray(Harray[k]);
    P_array           = bcm_AMGHierarchyPArray(Harray[k]);
    L_array           = bcm_AMGHierarchyLArray(Harray[k]);
    U_array           = bcm_AMGHierarchyUArray(Harray[k]);
    D_array           = bcm_AMGHierarchyDArray(Harray[k]);
    int num_levels    = bcm_AMGHierarchyNumLevels(Harray[k]);
    char filename[81];
    for (j=0; j<num_levels; j++) {
      sprintf(filename,"%s-P-l%3.3d.mtx",prefix,j);
      if (P_array[j]!=NULL) bcm_CSRMatrixPrintMM(P_array[j],filename);
      sprintf(filename,"%s-AC-l%3.3d.mtx",prefix,j);
      bcm_CSRMatrixPrintMM(A_array[j],filename);
    }

  }
}
typedef struct {
  int meshsize;
	double anisotropy;
	int direction;
  int solver_type;
  int max_hrc;
  double conv_ratio;
  int matchtype;
	double lambda;
  int aggrsweeps;
  int CR_it;
  int aggrtype;
  int max_levels;
  int cycle_type;
  int coarse_solver;
  int coarserelax_number;
  int relax_type;
  int prerelax_sweeps;
  int postrelax_sweeps;
  int itnlim;
  double rtol;
} parms_t;

int get_inp_data(FILE *fp, parms_t *inparms)
{
    inparms->meshsize          = int_get_fp(fp);
		inparms->anisotropy				 = double_get_fp(fp);
		inparms->direction				 = int_get_fp(fp);
    inparms->solver_type       = int_get_fp(fp);
    inparms->max_hrc           = int_get_fp(fp);
    inparms->conv_ratio        = double_get_fp(fp);
    inparms->matchtype         = int_get_fp(fp);
		inparms->lambda            = double_get_fp(fp);
    inparms->aggrsweeps        = int_get_fp(fp);
    inparms->CR_it             = int_get_fp(fp);
    inparms->aggrtype          = int_get_fp(fp);
    inparms->max_levels        = int_get_fp(fp);
    inparms->cycle_type        = int_get_fp(fp);
    inparms->coarse_solver          = int_get_fp(fp);
    inparms->coarserelax_number     = int_get_fp(fp);
    inparms->relax_type        = int_get_fp(fp);
    inparms->prerelax_sweeps   = int_get_fp(fp);
    inparms->postrelax_sweeps  = int_get_fp(fp);
    inparms->itnlim            = int_get_fp(fp);
    inparms->rtol              = double_get_fp(fp);
    return(0);
}

/* This program builds a Bootstrap AMG with a desired convergence rate
   and applies it as preconditioner for Flexible CG (F1CG) for solving
   the Laplace eq. in 3D on a unitary cube with Dirichlet boundary conditions*/

int  main(int argc, char *argv[])
{
   bcm_BootAMGBuildData *bootamg_data;
   bcm_AMGApplyData *amg_cycle;
   bcm_CSRMatrix *A;
   bcm_Vector *w;
   bcm_Vector *rhs, *Sol;
   bcm_Vector *soltrue;
   int i, *num_grid_sweeps, j, k;
   int size_rhs;
   bcm_Vector *res;
   double ratio, normnew, normold;
   parms_t inparms;
   FILE    *fp=NULL;

   fprintf(stdout,
	   "Welcome to BootCMatch version %s\n This is the testkbootsolvers program\n",
	   BCM_VERSION_STRING);

   if (argc > 1) {
     fp = fopen(argv[1], "r");
   } else {
     fp = stdin;
   }

   if (get_inp_data(fp,&inparms) != 0) {
     fprintf(stderr,"Error getting input parms \n");
     exit(1);
   }
   /* generate Laplace matrix */
    int n=inparms.meshsize;
    int N=n*n*n;
    int nnzval=0;
    int *cols, *rows, xi,yi,zi;
    double *values;
		double anisotropy=inparms.anisotropy;
    values= (double *) calloc(7*N, sizeof(double));
    cols= (int *) calloc(7*N, sizeof(int));
    rows= (int *) calloc(7*N, sizeof(int));
    for(i=0; i<N; i++){
    	 k=floor(i/(n*n));
         /* The left identity block:position i-n*n */
         if ((i-n*n)>=0){
            rows[nnzval] = i;
            cols[nnzval] = i-n*n;
						if (inparms.direction == 3){
            	values[nnzval] = -anisotropy;
						} else {
							values[nnzval] = -1.0;
						}
            nnzval=nnzval+1;
         }
         /* The left identity block:position i-n */
         if (i>=n+k*n*n && i < (k+1)*n*n) {
            rows[nnzval] = i;
            cols[nnzval] = i-n;
						if (inparms.direction == 2){
							values[nnzval] = -anisotropy;
						} else {
							values[nnzval] = -1.0;
						}
            nnzval=nnzval+1;
         }

         /* The left -1: position i-1 */
         if (i%n!=0) {
            rows[nnzval] = i;
            cols[nnzval] = i-1;
						if (inparms.direction == 1){
							values[nnzval] = -anisotropy;
						} else {
							values[nnzval] = -1.0;
						}
            nnzval=nnzval+1;
         }

         /* Set the diagonal: position i */
         rows[nnzval] = i;
         cols[nnzval] = i;
         values[nnzval] = 4.0 + 2*anisotropy;
         nnzval=nnzval+1;

         /* The right -1: position i+1 */
         if ((i+1)%n!=0) {
            rows[nnzval] = i;
            cols[nnzval] = i+1;
						if (inparms.direction == 1){
							values[nnzval] = -anisotropy;
						} else {
							values[nnzval] = -1.0;
						}
            nnzval=nnzval+1;
         }

         /* The right identity block:position i+n */

         if(i>=k*n*n && i < (k+1)*n*n-n){
            rows[nnzval] = i;
            cols[nnzval] = i+n;
						if (inparms.direction == 2){
							values[nnzval] = -anisotropy;
						} else {
							values[nnzval] = -1.0;
						}
            nnzval=nnzval+1;
         }

         /* The right identity block:position i+n*n */
         if ((i+n*n)< N){
            rows[nnzval] = i;
            cols[nnzval] = i+n*n;
						if (inparms.direction == 3){
							values[nnzval] = -anisotropy;
						} else {
							values[nnzval] = -1.0;
						}
            nnzval=nnzval+1;
         }
   }
  /*----------------------------------------------------------
   * Transform matrix from COO to CSR format
   *----------------------------------------------------------*/

  int *matrix_i, *matrix_j, k0, iad;
  double *matrix_data, x;

  matrix_i = (int *) calloc(N+1, sizeof(int));

  /* determine row lenght */
  for (j=0; j<nnzval; j++) matrix_i[rows[j]]=matrix_i[rows[j]]+1;

  /* starting position of each row */
  k=0;
  for(j=0; j<= N; j++)
    {
      k0=matrix_i[j];
      matrix_i[j]=k;
      k=k+k0;
    }
  matrix_j = (int *) calloc(nnzval, sizeof(int));
  matrix_data = (double *) calloc(nnzval, sizeof(double));

  /* go through the structure once more. Fill in matrix */
  for(k=0; k<nnzval; k++)
    {
      i=rows[k];
      j=cols[k];
      x=values[k];
      iad=matrix_i[i];
      matrix_data[iad]=x;
      matrix_j[iad]=j;
      matrix_i[i]=iad+1;
    }
  /* shift back matrix_i */
  for(j=N-1; j>=0; j--) matrix_i[j+1]=matrix_i[j];
  matrix_i[0]=0;

  A = bcm_CSRMatrixCreate(N, N, nnzval);
  bcm_CSRMatrixI(A) = matrix_i;
  bcm_CSRMatrixJ(A) = matrix_j;
  bcm_CSRMatrixData(A) = matrix_data;
  bcm_CSRMatrixNumNonzeros(A) = nnzval;

  free(rows);
  free(cols);
  free(values);

  /* bcm_CSRMatrixPrintMM(A,"lap3d.mtx"); */
  int num_rows = bcm_CSRMatrixNumRows( A );
  int num_cols = bcm_CSRMatrixNumCols( A );
  int nnz = bcm_CSRMatrixNumNonzeros( A );
  printf("numrows: %d\n",num_rows);
  printf("numcols: %d\n",num_cols);
  printf("nnz: %d\n",nnz);

     /* generate rhs vector */
     rhs=bcm_VectorCreate(num_rows);
     bcm_VectorInitialize(rhs);
     size_rhs=bcm_VectorSize(rhs);
     printf("size_rhs: %d\n",size_rhs);
     bcm_VectorSetConstantValues(rhs,1.0);
     //bcm_VectorSetRandomValues(rhs,1.0);

   assert(num_rows == size_rhs);

   /* initialize data structure for AMG building. See the setup
      routines for changing default values */
   bootamg_data = bcm_BootAMGBuildDataInitialize();
   bcm_AMGBuildData *amg_data;
   amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);
   amg_cycle= bcm_AMGCycleInitialize();

   bcm_AMGBuildDataCSRMatrix(amg_data)=A;

   /* reading AMG algorithmic parameters */

   /* composition type  and maximum number of
      components for the bootstrap AMG */
   /* Single Hierarchy Type: matchingtype, number of aggregation sweeps
      for aggressive coarsening and aggregation type (smoothed or unsmoothed) */

   bcm_BootAMGBuildSetSolverType(bootamg_data, inparms.solver_type);
   bcm_BootAMGBuildSetMaxHrc(bootamg_data, inparms.max_hrc);
   bcm_BootAMGBuildSetDesiredRatio(bootamg_data, inparms.conv_ratio);
   bcm_AMGBuildSetAggMatchType(amg_data, inparms.matchtype);
	 bcm_AMGBuildSetLambda(amg_data, inparms.lambda);
	 bcm_AMGBuildSetSweepNumber(amg_data, inparms.aggrsweeps);
   bcm_AMGBuildSetCRIterations(amg_data, inparms.CR_it);
   bcm_AMGBuildSetAggInterpType(amg_data, inparms.aggrtype);
   bcm_AMGBuildSetMaxLevels(amg_data, inparms.max_levels);
   bcm_AMGBuildSetCoarseSolver(amg_data, inparms.coarse_solver);
   bcm_AMGBuildSetCoarseRelaxNum(amg_data, inparms.coarserelax_number);
   bcm_AMGSetCycleType(amg_cycle, inparms.cycle_type);
   bcm_AMGSetRelaxType(amg_cycle, inparms.relax_type);
   bcm_AMGSetPreRelaxSteps(amg_cycle, inparms.prerelax_sweeps);
   bcm_AMGSetPostRelaxSteps(amg_cycle, inparms.postrelax_sweeps);


   printf("Composite Solver Type: %d\n",bcm_BootAMGBuildDataCompType(bootamg_data));
   printf("Max number of components: %d\n",bcm_BootAMGBuildDataMaxHrc(bootamg_data));
   printf("Desired Conv. Ratio: %e\n",bcm_BootAMGBuildDataDesRatio(bootamg_data));
   printf("Matching Type: %d\n",bcm_AMGBuildDataAggMatchType(amg_data));
   printf("Aggregation sweeps: %d\n",bcm_AMGBuildDataSweepNumber(amg_data));
   printf("CR_it: %d\n",bcm_AMGBuildDataCRIterations(amg_data));
   printf("Aggregation Type: %d\n",bcm_AMGBuildDataAggInterpType(amg_data));
   printf("max_levels: %d\n",bcm_AMGBuildDataMaxLevels(amg_data));
   printf("coarse_solver: %d\n",bcm_AMGBuildDataCoarseSolver(amg_data));
   printf("coarserelax_number: %d\n",bcm_AMGBuildDataCoarseRelaxNum(amg_data));
   printf("cycle_type: %d\n",bcm_AMGApplyDataCycleType(amg_cycle));
   printf("relax_type: %d\n",bcm_AMGApplyDataRelaxType(amg_cycle));
   printf("prerelax_sweeps: %d\n",bcm_AMGApplyDataPreRelax(amg_cycle));
   printf("postrelax_sweeps: %d\n",bcm_AMGApplyDataPostRelax(amg_cycle));

   printf("itnlim: %d\n",inparms.itnlim);
   printf("rtol: %e\n",inparms.rtol);

   if (argc > 1) fclose(fp);


   /* set maximum coarse size */

   int maxcoarse=40*pow((double)num_rows,(double)1/3);
   bcm_AMGBuildSetMaxCoarseSize(amg_data, maxcoarse);
   printf("maxcoarsesize: %d\n",bcm_AMGBuildDataMaxCoarseSize(amg_data));


   /* initialize num_grid_sweeps parameter w.r.t. the number of levels
      and the chosen cycle.
      In the final release we have to manage this in a setup routine after hierarchy building */

   num_grid_sweeps = (int *) calloc(inparms.max_levels-1, sizeof(int));
   for(i=0; i<inparms.max_levels-1; i++) num_grid_sweeps[i]=1;
   switch(inparms.cycle_type)  {
   case 1: /* H-cycle */
     {
       for(i=0; i<inparms.max_levels-1; i++) {
	 j=i%2; /*step is fixed to 2; it can be also different */
	 if(j==0) num_grid_sweeps[i]=2; /* if num_grid_sweeps is 2, we are using a hybrid V-W cycle */
       }
     }
     break;
   case 2: /* W-cycle */
     {
       for(i=0; i<inparms.max_levels-2; i++) num_grid_sweeps[i]=2;
     }
     break;
   }
   bcm_AMGApplyDataGridSweeps(amg_cycle)=num_grid_sweeps;


   /* set arbitrary initial (smooth) vector: generally
      we use unitary vector or random vector */

   w=bcm_VectorCreate(num_rows);
   bcm_VectorInitialize(w);
   int w_size=bcm_VectorSize(w);
   printf("wsize: %d\n",w_size);
   /* bcm_VectorSetRandomValues(w,1.0); */
   bcm_VectorSetConstantValues(w,1.0);

   bcm_AMGBuildDataSmoothVector(amg_data)= w;

   /* start bootstrap process */

   printf("Bootstrap starting \n");

   bcm_BootAMG *boot_amg;
   double time1=time_getWallclockSeconds();
   boot_amg=bcm_Bootstrap(bootamg_data,amg_cycle);
   double time2=time_getWallclockSeconds()-time1;

   printf("Bootstrap ended\n");

   bcm_AMGHierarchy **Harray;

   printf("Number of components:  %d\n", bcm_BootAMGNHrc(boot_amg));
   printf("Estimated convergence %e \n", bcm_BootAMGEstRatio(boot_amg));
   printf("Information on the Components\n");
   Harray=bcm_BootAMGHarray(boot_amg);
   double avgcmpx=0;
   double avgwcmpx=0;
   double avgnumlev=0;
   for(k=0;k<bcm_BootAMGNHrc(boot_amg); k++) {
     printf("Component:  %d\n", k);
     printf("Number of levels:  %d\n", bcm_AMGHierarchyNumLevels(Harray[k]));
     printf("Operator cmplx for V-cycle %e \n", bcm_AMGHierarchyOpCmplx(Harray[k]));
     printf("Operator cmplx for W-cycle %e \n", bcm_AMGHierarchyOpCmplxW(Harray[k]));
     avgcmpx=avgcmpx+ bcm_AMGHierarchyOpCmplx(Harray[k]);
     avgwcmpx=avgwcmpx+ bcm_AMGHierarchyOpCmplxW(Harray[k]);
     avgnumlev=avgnumlev+ bcm_AMGHierarchyNumLevels(Harray[k]);
   }
   printf("Wall Clock Time for Building:  %e\n", time2);
   avgcmpx=avgcmpx/bcm_BootAMGNHrc(boot_amg);
   avgwcmpx=avgwcmpx/bcm_BootAMGNHrc(boot_amg);
   avgnumlev=avgnumlev/bcm_BootAMGNHrc(boot_amg);
   printf("AVG cmpx for V-cycle %e\n", avgcmpx);
   printf("AVG cmpx for W-cycle %e\n", avgwcmpx);
   printf("AVG numlev  %e\n", avgnumlev);

   /* generate initial vector */
   Sol=bcm_VectorCreate(num_rows);
   bcm_VectorInitialize(Sol);

   int precon=1, itn, istop;
   double timetot;

   time1=time_getWallclockSeconds();
   istop=bcm_bootpcg(Sol, rhs, bootamg_data, boot_amg, amg_cycle, precon,
			 inparms.itnlim, inparms.rtol, &itn,&timetot);
   time2=time_getWallclockSeconds()-time1;
   printf("Wall Clock Time for Applying Preconditioner:  %e\n", timetot);
   printf("Wall Clock Time for Solving:  %e\n", time2);
   printf("istop=  %d\n", istop);
   printf("itn=  %d\n", itn);

   char *file_out="sol.mtx";
   bcm_VectorPrint(Sol,file_out);
/*
#ifdef DUMP_HIER
   dump_on_file("BCM-",boot_amg);
#endif
*/
   free(num_grid_sweeps);
   bcm_AMGApplyDataDestroy(amg_cycle);
   bcm_BootAMGBuildDataDestroy(bootamg_data);
   bcm_BootAMGDestroy(boot_amg);
   bcm_CSRMatrixDestroy(A);
   bcm_VectorDestroy(w);
   bcm_VectorDestroy(rhs);
   bcm_VectorDestroy(Sol);
   return (0);
}
