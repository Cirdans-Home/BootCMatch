/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra    IAC-CNR
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
#ifndef BCM_AMG_H_
#define BCM_AMG_H_
#include "bcm_matvec.h"
#include "bcm_util.h"
#ifdef HAVE_SUPERLU
#include "slu_ddefs.h"

typedef struct {
  SuperMatrix *L;
  SuperMatrix *U;
  int *perm_c;
  int *perm_r;
} factors_t;
#endif


/*--------------------------------------------------------------------------
 * bcm_AMGHierarchy includes data defining a given hierarchy
 *--------------------------------------------------------------------------*/

typedef struct
{
  /* data generated in the setup phase */

  bcm_CSRMatrix **A_array; /*array of coarse matrices including the fine ones*/
  bcm_CSRMatrix **P_array; /*array of prolongators */
  bcm_CSRMatrix **L_array; /* lower triangular part of coarse matrices */
  bcm_CSRMatrix **U_array; /* upper triangular part of coarse matrices */
  bcm_Vector    **D_array; /* diagonal of coarse matrices */
#ifdef HAVE_SUPERLU
  factors_t     *SLUfactors;
#endif
  int           num_levels; /* number of levels of the hierachy */
  double        op_cmplx;  /* operator complexity of the hierarchy for V-cycle*/
  double        op_wcmplx;  /* operator complexity of the hierarchy for W-cycle*/
  double        avg_cratio;  /* average of coarsening ratio of the hierarchy */

} bcm_AMGHierarchy;

/*--------------------------------------------------------------------------
 * Accessor functions for the bcm_AMGHierarchy structure
 *--------------------------------------------------------------------------*/

/* data generated by the setup phase */
#define bcm_AMGHierarchyAArray(amg_hierarchy) ((amg_hierarchy)->A_array)
#define bcm_AMGHierarchyPArray(amg_hierarchy) ((amg_hierarchy)->P_array)
#define bcm_AMGHierarchyLArray(amg_hierarchy) ((amg_hierarchy)->L_array)
#define bcm_AMGHierarchyUArray(amg_hierarchy) ((amg_hierarchy)->U_array)
#define bcm_AMGHierarchyDArray(amg_hierarchy) ((amg_hierarchy)->D_array)
#define bcm_AMGHierarchyNumLevels(amg_hierarchy) ((amg_hierarchy)->num_levels)
#define bcm_AMGHierarchyOpCmplx(amg_hierarchy) ((amg_hierarchy)->op_cmplx)
#define bcm_AMGHierarchyOpCmplxW(amg_hierarchy) ((amg_hierarchy)->op_wcmplx)
#define bcm_AMGHierarchyAvgCratio(amg_hierarchy) ((amg_hierarchy)->avg_cratio)
#ifdef HAVE_SUPERLU
#define bcm_AMGHierarchyLUfactors(amg_hierarchy) ((amg_hierarchy)->SLUfactors)
#endif

/*--------------------------------------------------------------------------
 * bcm_AMGBuildData for building hierarchies
 *--------------------------------------------------------------------------*/

typedef struct
{

   /* setup params */
   int      maxlevels; /* max number of levels per hierarchy */
   int      maxcoarsesize; /* maximum size of the coarsest matrix */
   int      sweepnumber; /* number of pairwise aggregation steps. Currently 0 for pairwise and 1 for double pairwise */
   int      agg_interp_type; /* 1 for smoothed aggregation, 0 for pure aggregation */
   int      agg_match_type;  /* 4 for 2/3 lambda matching (lambda-matchTwoThirdeps) 3 for exact lambda matching (lambda-matchOpt) 1 for exact matching (HSL-mc64), 2 for auction based matching, 0 for Preis approximate matching */
   int        coarse_solver; /* solver to be used on the coarsest level */
   /*     relax_type/coarse_solver = 0 ->  1 sweep of Jacobi
    *     relax_type/coarse_solver = 1 ->  1 sweep of Gauss-Seidel
    *     relax_type/coarse_solver = 2 ->  1 sweep of symm. Gauss-Seidel
    *     relax_type/coarse_solver = 9 -> Direct Solve */

   /* CR params */
   int      CRrelax_type; /* to choose relaxation scheme for Compatible Relaxation */
   double   CRrelax_weight; /* weight for weighted Jacobi in CR */
   int      CRit;  /* number of iterations for Compatible Relaxation */
   double   CRratio; /* optimal convergence ratio in Compatible Relaxation to stop coarsening*/

   /* problem data */
   bcm_CSRMatrix  *A;   /* problem matrix */
   bcm_Vector  *w;   /* current smooth vector for building new hierarchy: NB In the bootstrap process
			we update it at each new step */

} bcm_AMGBuildData;

/*--------------------------------------------------------------------------
 * Accessor functions for the bcm_AMGBuildData structure
 *--------------------------------------------------------------------------*/

#define bcm_AMGBuildDataMaxLevels(amg_data) ((amg_data)->maxlevels)
#define bcm_AMGBuildDataAggInterpType(amg_data) ((amg_data)->agg_interp_type)
#define bcm_AMGBuildDataAggMatchType(amg_data) ((amg_data)->agg_match_type)
#define bcm_AMGBuildDataMaxCoarseSize(amg_data) ((amg_data)->maxcoarsesize)
#define bcm_AMGBuildDataSweepNumber(amg_data) ((amg_data)->sweepnumber)
#define bcm_AMGBuildDataCoarseSolver(amg_data) ((amg_data)->coarse_solver)
#define bcm_AMGBuildDataCRRelaxType(amg_data) ((amg_data)->CRrelax_type)
#define bcm_AMGBuildDataCRRelaxWeight(amg_data) ((amg_data)->CRrelax_weight)
#define bcm_AMGBuildDataCRIterations(amg_data) ((amg_data)->CRit)
#define bcm_AMGBuildDataCRratio(amg_data) ((amg_data)->CRratio)
#define bcm_AMGBuildDataCSRMatrix(amg_data) ((amg_data)->A)
#define bcm_AMGBuildDataSmoothVector(amg_data) ((amg_data)->w)


/*--------------------------------------------------------------------------
 * bcm_AMGApplyData for applying built hierarchies
 *--------------------------------------------------------------------------*/

typedef struct
{

   /* cycle type params */
   int        cycle_type; /* 0 for V-cycle, 2 for W-cycle, 1 for G-cycle */
   int        *num_grid_sweeps; /* number of sweeps on a fixed level in case of G-cycle */

   int        relax_type; /* type of pre/post relaxation/smoothing */
   int        prerelax_number; /* number of pre smoothing steps */
   int        postrelax_number; /* number of post smoothing steps */

   double     relax_weight; /* weight for Jacobi relaxation */

} bcm_AMGApplyData;

/*--------------------------------------------------------------------------
 * Accessor functions for the bcm_AMGApplyData structure
 *--------------------------------------------------------------------------*/

#define bcm_AMGApplyDataCycleType(amg_cycle) ((amg_cycle)->cycle_type)
#define bcm_AMGApplyDataGridSweeps(amg_cycle) ((amg_cycle)->num_grid_sweeps)
#define bcm_AMGApplyDataRelaxType(amg_cycle) ((amg_cycle)->relax_type)
#define bcm_AMGApplyDataPreRelax(amg_cycle) ((amg_cycle)->prerelax_number)
#define bcm_AMGApplyDataPostRelax(amg_cycle) ((amg_cycle)->postrelax_number)
#define bcm_AMGApplyDataRelaxWeight(amg_cycle) ((amg_cycle)->relax_weight)

/* bcm_amg.c */
void * bcm_AMGInitialize();
int bcm_AMGBuildSetMaxLevels(void *data, int  max_levels);
int bcm_AMGBuildSetAggInterpType(void *data, int interp_type);
int bcm_AMGBuildSetAggMatchType(void *data, int match_type);
int bcm_AMGBuildSetMaxCoarseSize(void  *data, int maxcoarse_size);
int bcm_AMGBuildSetSweepNumber(void *data, int sweepnumber);
int bcm_AMGBuildSetCoarseSolver(void *data, int coarse_solver);
int bcm_AMGBuildSetCRRelaxType(void *data, int crrelax_type);
int bcm_AMGBuildSetCRRelaxWeight(void *data, double crrelax_weight);
int bcm_AMGBuildSetCRIterations(void *data, int criterations);
int bcm_AMGBuildSetCRRatio(void *data, double crratio);
int bcm_AMGBuildDataDestroy(void *data);
bcm_AMGHierarchy *bcm_AMGHierarchyCreate(int maxlevels);
int bcm_AMGHierarchyInitialize(bcm_AMGHierarchy *amg_hierarchy);
int bcm_AMGHierarchyDestroy(bcm_AMGHierarchy *amg_hierarchy);
void * bcm_AMGCycleInitialize();
int bcm_AMGSetCycleType(void *data, int cycle_type);
int bcm_AMGSetNumGridSweeps(void *data, int *num_grid_sweeps);
int bcm_AMGSetRelaxType(void *data, int relax_type);
int bcm_AMGSetPreRelaxSteps(void *data, int prerelax_number);
int bcm_AMGSetPostRelaxSteps(void *data, int postrelax_number);
int bcm_AMGSetRelaxWeight(void *data, double relax_weight);
int bcm_AMGApplyDataDestroy(void *data);

/* bcm_weightmat.c */
bcm_CSRMatrix *bcm_CSRMatrixAhat( bcm_CSRMatrix *A, bcm_Vector *w, int match_type );

/* bcm_csr_linmatch.c */
int *bcm_CSRMatrixHMatch( bcm_CSRMatrix *A );
int bcm_trymatch(int rno, int cno, bcm_CSRMatrix *W, int *jrowindex,
		 int ljrowindex, int *jjrowindex, int ljjrowindex, int *rmatch);

/* bcm_adaptivecoarsening.c */
bcm_CSRMatrix * bcm_CSRMatchingAgg(bcm_CSRMatrix *A, bcm_Vector **w,
				   bcm_CSRMatrix **P, int match_type, int num_sweeps,
				   int max_sizecoarse, int max_levels, int *ftcoarse,
				   int cr_it, int cr_relax_type, double cr_relax_weight);
int bcm_CSRMatchingPairAgg(bcm_CSRMatrix *A, bcm_Vector *w, bcm_CSRMatrix **P, int match_type);
bcm_AMGHierarchy * bcm_AdaptiveCoarsening(bcm_AMGBuildData *amg_data);

#ifdef HAVE_AMGMATCH
void c_matchLambdaOpt(int n, int nnz, int *s, int *t, double *edgeWght, double lambda, int *mateNode);
void c_matchLambdaTwoThirdeps(int n, int nnz, int *s, int *t, double *edgeWght, double lambda, int *mateNode);
#endif

#endif
