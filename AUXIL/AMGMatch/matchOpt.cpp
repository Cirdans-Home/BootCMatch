#include "Types.h"

#include <omp.h>

#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/math.h>

#include "AMGMatch/amgmatch.h"

using namespace std;
using lemon::SmartGraph;
using lemon::MaxWeightedMatching;
using lemon::INVALID;

/*
 * This function takes a weighted edge list and compute optimal matching using LEMON graph library
 * n: Number of Vertices
 * nnz: Number of Edges
 * s : vector of one end-points
 * t : vector of the other end-points.
 * edgeWgthVec: weight vector of the edges. So for any i, (s[i],t[i],edgeWgthVec[i]) is an edge with weight
 * lambda: a real value represent the weight translation
 * mateNode: output vector of the mates. -1 means the node is not matched.
 * stat: Statistics of the matching; time, cardinality and weight
 */
void matchLambdaOpt(Size n, Size nnz, vector<Size> &s, vector<Size> &t, vector<Val> &edgWghtVec, Val lambda, vector<Size> &mateNode, Amgmatch_stat &stat)
{

    //If the lambda is non-zero translate the weights
    if(fabs(lambda) > 0)
        trnslWght(s,t,edgWghtVec,lambda);
    //We are using the SmartGraph class. It is th fastest but do not allow deletion
    SmartGraph G;
    //creating a set of lemon nodes
    vector<SmartGraph::Node> nodes(n);
    //typedefining the two maps for easy access
    typedef SmartGraph::EdgeMap<Val> WeightMap;
    typedef SmartGraph::NodeMap<Size> IdMap;

    //Map that contain the weights on the edge
    WeightMap edgeMap(G);
    //Map that containts the id of nodes. NOTE: There is a id in LEMON for each node but I do not know whether they will be consistent to our numbering of the nodes. To be safe I declared a node map that maps the nodes to their ids.
    IdMap idNode(G);

    //add the nodes and the ids in the map
    for(Size i=0;i<n;i++)
    {
        nodes[i] = G.addNode(); 
        idNode[nodes[i]] = i;
    }
    
    //add the edges and the correspoding weights on the map
    for(Size  i=0;i<s.size();i++)
    {
        SmartGraph::Edge e = G.addEdge(nodes[s[i]],nodes[t[i]]); 
        edgeMap[e] = edgWghtVec[i];
    }
    
    //The main part of the algorithm
    double startTime = omp_get_wtime();
    //create a object. NOTE: The example in LEMON is flawed. You have to send two template argument; the types of the inputs.
    MaxWeightedMatching<SmartGraph,WeightMap> mwm(G,edgeMap);
    //fractional initialization is the best one.
    mwm.fractionalInit();
    //run the algorithm
    mwm.run();
    double endTime = omp_get_wtime();
    //gather the statistics
    stat.matchCard = mwm.matchingSize();
    stat.matchWght = mwm.matchingWeight() - stat.matchCard*lambda;
    stat.matchTime = endTime - startTime;
    stat.lambda = lambda;
    
    //check whether the vector is allocated or not. If not then resize it.
    if(mateNode.empty()) 
        mateNode.resize(n);
    //find the matching nodes
    //We iterate trhough the nodes and find the id of a node. Then we find the mate node of the node. If the node is not invalid (i.e., matched) we find the corresponding id through the node map. TODO: Is there any better way of doing this thing without using NodeMap?
    for(SmartGraph::NodeIt nodeIt(G); nodeIt != INVALID; ++nodeIt)
    {
        Size u = idNode[nodeIt];
        SmartGraph::Node nodeV = mwm.mate(nodeIt);
        if(nodeV != INVALID)
        {
            Size v = idNode[nodeV];
            mateNode[u] = v; 
        }
        else
        {
            mateNode[u] =-1; 
        }
    }
}
