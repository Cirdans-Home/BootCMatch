#ifdef HAVE_AMGMATCH
#include <iostream>
#include <numeric>
#include <string>
#include <fstream>
#include "AMGMatch/amgmatch.h"

extern "C" void c_matchLambdaOpt(int n, int nnz, int *s, int *t, double *edgeWght, double lambda, int *mateNode){

  vector<int> vs,vt;
  vector<double> vedgeWght;
  vector<int> vmateNode;
  Amgmatch_stat pstat;

  // We need some C++ to copy the input Array into a C++ Vector
  vector<int> vs2(&s[0], &s[nnz]);
  vs.insert(vs.end(), vs2.begin(), vs2.end());
  vector<int> vt2(&t[0], &t[nnz]);
  vt.insert(vt.end(), vt2.begin(), vt2.end());
  vector<double> vedgeWght2(&edgeWght[0],&edgeWght[nnz]);
  vedgeWght.insert(vedgeWght.end(), vedgeWght2.begin(), vedgeWght2.end());

  matchLambdaOpt(n,nnz,vs,vt,vedgeWght,lambda,vmateNode,pstat);
  std::copy(vmateNode.begin(), vmateNode.end(),mateNode);

  matchingStat(n,nnz,pstat);
}

extern "C" void c_matchLambdaTwoThirdeps(int n, int nnz, int *s, int *t, double *edgeWght, double lambda, int *mateNode){

  vector<int> vs,vt;
  vector<double> vedgeWght;
  vector<int> vmateNode;
  Amgmatch_stat pstat;

  // We need some C++ to copy the input Array into a C++ Vector
  vector<int> vs2(&s[0], &s[nnz]);
  vs.insert(vs.end(), vs2.begin(), vs2.end());
  vector<int> vt2(&t[0], &t[nnz]);
  vt.insert(vt.end(), vt2.begin(), vt2.end());
  vector<double> vedgeWght2(&edgeWght[0],&edgeWght[nnz]);
  vedgeWght.insert(vedgeWght.end(), vedgeWght2.begin(), vedgeWght2.end());

  clearStat(pstat);

  fprintf(stderr, "2/3-ε matrix size %d x %d, nnz = %d, λ = %1.2f \n",n,n,nnz,lambda);
  matchLambdaTwoThirdeps(n,nnz,vs,vt,vedgeWght,lambda,vmateNode,pstat);
  fprintf(stderr, "2/3-ε Matching Complete\n");
  std::copy(vmateNode.begin(), vmateNode.end(),mateNode);

  matchingStat(n,nnz,pstat);

}

#endif
