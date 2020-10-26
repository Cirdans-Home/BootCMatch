
#include "Types.h"
//This is a graph from matchbox
#include "Graph.h"
#include "ConvertEngine.h"
#include "AdjGraphEngine.h"
#include "MatchingEngine.h"

#include "AMGMatch/amgmatch.h"
//#include "AMGMatch/amgmatchUtil.h"
using namespace Matchbox;
using namespace std;

/*This function takes a weighted edge list and compute 2/3-eps maximum weight matching
 * n: Number of Vertices
 * nnz: Number of Edges
 * s : vector of one end-points
 * t : vector of the other end-points.
 * edgeWgthVec: weight vector of the edges. So for any i, (s[i],t[i],edgeWgthVec[i]) is an edge with weight
 * lambda: A real value represent the weight translation
 * mateNode: output vector of the mates. -1 means the node is not matched.
 * stat: Statistics of the matching; time, cardinality and weight
 */
void matchLambdaTwoThirdeps(Size n, Size nnz, vector<Size> &s, vector<Size> &t, vector<Val> &edgWghtVec, Val lambda, vector<Size> &mateNode, Amgmatch_stat &stat)
{

    //If the lambda is non-zero translate the weights
    if(fabs(lambda) > 0)
        trnslWght(s,t,edgWghtVec,lambda);
    //create a non-bipartite graph out of source and destination vertex
    Graph g;

    Err err = g.Set(n,nnz,s,t,edgWghtVec);
    if(err != eErrNone)
    {
        cout<<"error in setting graph"<<endl;
        exit(1);
    }

    //Matchbox matchingEngine for the two-third epsilon matching algorithm
    MatchingEngine matchingEngine;


    //The matching statistics
    Size card;
    Val wght;
    Size numedgs;

    if(mateNode.empty())
        mateNode.resize(n);

    double startTime = omp_get_wtime();
    //The number of phases
    Size phases = 5;
    matchingEngine.RandomOrderTwoThirdMinusEpsilonEdgWght(g, &mateNode, &card, &wght, phases);

    //matchingEngine.ComputeHalfEdgWghtMatching(g, &mateVec, &card, &wght );
    //std::copy(mateVec.begin(),mateVec.end(),mateNode);

    double endTime = omp_get_wtime();
    bool valid = matchingEngine.CheckMatching(g, mateNode, card);
    assert(valid==true);

    //Iterate through the nodes and replace the mate of unmatched nodes by -1
    for(int i=0;i<n;i++)
    {
        if(mateNode[i] == cNullItm)
            mateNode[i] = -1;
    }
    //gather the matching statistics
    stat.matchCard = card;
    stat.matchWght = wght - stat.matchCard*lambda;
    stat.matchTime = endTime - startTime;
    stat.lambda  = lambda;
}
