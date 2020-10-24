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
  unsigned sSize = sizeof(s) / sizeof(int);
  unsigned tSize = sizeof(t) / sizeof(int);
  unsigned edgeWghtSize = sizeof(edgeWght) / sizeof(double);

  vector<int> vs2(&s[0], &s[sSize]);
  vs.insert(vs.end(), vs2.begin(), vs2.end());
  vector<int> vt2(&t[0], &t[tSize]);
  vt.insert(vt.end(), vt2.begin(), vt2.end());
  vector<int> vedgeWght2(&edgeWght[0],&edgeWght[edgeWghtSize]);
  vedgeWght.insert(vedgeWght.end(), vedgeWght2.begin(), vedgeWght2.end());

  matchLambdaOpt(n,nnz,vs,vt,vedgeWght,lambda,vmateNode,pstat);
  mateNode = &vmateNode[0];

  matchingStat(n,nnz,pstat);
}

extern "C" void c_matchLambdaTwoThirdeps(int n, int nnz, int *s, int *t, double *edgeWght, double lambda, int *mateNode){

  vector<int> vs,vt;
  vector<double> vedgeWght;
  vector<int> vmateNode;
  Amgmatch_stat pstat;

  // We need some C++ to copy the input Array into a C++ Vector
  unsigned sSize = sizeof(s) / sizeof(int);
  unsigned tSize = sizeof(t) / sizeof(int);
  unsigned edgeWghtSize = sizeof(edgeWght) / sizeof(double);

  vector<int> vs2(&s[0], &s[sSize]);
  vs.insert(vs.end(), vs2.begin(), vs2.end());
  vector<int> vt2(&t[0], &t[tSize]);
  vt.insert(vt.end(), vt2.begin(), vt2.end());
  vector<int> vedgeWght2(&edgeWght[0],&edgeWght[edgeWghtSize]);
  vedgeWght.insert(vedgeWght.end(), vedgeWght2.begin(), vedgeWght2.end());

  clearStat(pstat);
  
  matchLambdaTwoThirdeps(n,nnz,vs,vt,vedgeWght,lambda,vmateNode,pstat);
  mateNode = &vmateNode[0];

  matchingStat(n,nnz,pstat);

}

#endif
