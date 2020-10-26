#pragma once

#include <vector>
#include "Types.h"

using  std::vector;
using  Matchbox::Size;
using  Matchbox::Val;

struct Edge
{
  Size r;
  Size c;
  Val w;
};

struct Amgmatch_stat
{
  Val matchWght;
  Val matchTime;
  Size matchCard;
  Val lambda;

};

void matchLambdaTwoThirdeps(Size,  Size, vector<Size> &, vector<Size> &, vector<Val> &,  Val, vector<Size> &,  Amgmatch_stat &);
void matchLambdaOpt(Size,  Size, vector<Size> &, vector<Size> &, vector<Val> &,  Val, vector<Size> &,  Amgmatch_stat &);
void matchHalf(Size,  Size, vector<Size> &, vector<Size> &, vector<Val> &,  Size*,  Amgmatch_stat &);

void clearStat(Amgmatch_stat &);
void trnslWght(vector<Size> &, vector<Size> &, vector<Val> &, Val );
void readMtx(std::string , Size &n,  Size &nnz, Val *&, Size *&, Size *&, Val threshold);
void readMtxEdgeList(std::string , Size &,Size &, vector<Size> &, vector<Size> &, vector<Val> &, Val );
int writeMatching(std::string , Size , vector<Size> &);
void matchingStat(Size , Size , Amgmatch_stat &);
