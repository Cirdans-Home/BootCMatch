#include <iostream>
#include <assert.h>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Types.h"
#include "AMGMatch/amgmatch.h"
using std::cout;
using std::endl;

#define DEBUG 0

void clearStat(Amgmatch_stat &stat)
{
    stat.matchCard = 0;
}

void trnslWght(vector<Size> &s, vector<Size> &t, vector<Val> &edgWght, Val lambda)
{
    for(Size i=0;i<edgWght.size();i++)
    {
        edgWght[i] = edgWght[i]+lambda;
    }
}
/*void writeMtx(string fname, int n1,int n2, int nnz, double *nzVal, int *rowInd, int *colPtr)
{
    std::ofstream f(fname.c_str());
    f<<"%%MatrixMarket matrix coordinate real general"<<endl;
    f<<n1<<" "<<n2<<" "<<nnz<<endl;

    for(int i=0;i<n2;i++)
    {
        for(int j=colPtr[i];j<colPtr[i+1];j++)
        {
            f<<rowInd[j]+1 <<" "<<i+1<<" "<<nzVal[j]<<endl;
        }
    }
}*/

void findMatrixProperty(std::string firstLine, bool &symmetric)
{
    symmetric = false;
    std::size_t found=firstLine.find("symmetric");
    if (found!=std::string::npos)
        symmetric = true;
}
void readMtx(std::string fileName, Size &n,  Size &nnz, Val *&nzVal, Size *&colInd, Size *&rowPtr, Val threshold)
{
    std::ifstream fileread(fileName.c_str());
    if(fileread.is_open()==false)
    {
        std::cout << "No file named "<<fileName<<std::endl;
        std::exit(1);
    }
    std::string firstLine;
    //fileread.getline(firstLine);
    std::getline(fileread,firstLine);
    bool symmetric = false;
    //matrix property
    findMatrixProperty(firstLine,symmetric);
    //Ignore header and comments
    while (fileread.peek() == '%') fileread.ignore(2048, '\n');

    Size nrow,ncol;
    Size m;
    fileread >> nrow >> ncol >> m;
    if(DEBUG==1)
    {
        std::cout<<nrow<<" "<<ncol<<" "<<nnz<<std::endl;
    }

    //std::vector<int> adjList[nrow];
    //std::vector<double> adjWeight[nrow];
    std::vector<std::vector<Size> > adjList(n);
    std::vector<std::vector<Val> > adjWeight(n);

    Size u;
    Size v;
    Val weight;
    nnz = 0;
    for(Size i=0;i<m;i++)
    {
        fileread >> u >> v >> weight;
        v--;
        u--;
        if(DEBUG==2)
        {
            std::cout<<u<<" "<<v<<" "<<weight<<std::endl;
        }
        if(fabs(weight)>threshold)
        {
            if(weight < 0)
               weight = (-1)*weight;

            adjList[u].push_back(v);
            adjWeight[u].push_back(weight);
            nnz++;
            if(symmetric==true && u!=v)
            {
                adjList[v].push_back(u);
                adjWeight[v].push_back(weight);
                nnz++;
            }
        }
    }
    nzVal = new Val[nnz];
    rowPtr = new Size[n+1];
    colInd = new Size[nnz];

    //second pass. Build the three arrays
    rowPtr[0] = 0;
    Size k=0;
    for(Size i=0;i<n;i++)
    {
        for(Size j=0;j<adjList[i].size();j++)
        {

            nzVal[k] = adjWeight[i][j];
            colInd[k] = adjList[i][j];
            k++;
        }
        rowPtr[i+1]= rowPtr[i] + adjList[i].size();
    }

    //delete the vectors
    for(Size i=0;i<nrow;i++)
    {
        adjList[i].clear();
        std::vector<Size>().swap(adjList[i]);

        adjWeight[i].clear();
        std::vector<Val>().swap(adjWeight[i]);

    }

}

void readMtxEdgeList(std::string fileName, Size &n,Size &nnz, vector<Size> &s, vector<Size> &t, vector<Val> &edgeWght, Val threshold)
{
    std::ifstream fileread(fileName.c_str());
    if(fileread.is_open()==false)
    {
        std::cout << "No file named "<<fileName<<std::endl;
        std::exit(1);
    }
    std::string firstLine;
    //fileread.getline(firstLine);
    std::getline(fileread,firstLine);
    bool symmetric = false;
    //matrix property
    findMatrixProperty(firstLine,symmetric);
    //Ignore header and comments
    while (fileread.peek() == '%') fileread.ignore(2048, '\n');

    Size nrow,ncol;
    Size m;
    fileread >> nrow >> ncol >> m;
    n = nrow;
    if(DEBUG==1)
    {
        std::cout<<nrow<<" "<<ncol<<" "<<nnz<<std::endl;
    }

    //std::vector<int> adjList[nrow];
    //std::vector<double> adjWeight[nrow];
    s.reserve(m);
    t.reserve(m);
    edgeWght.reserve(m);

    Size u;
    Size v;
    Val weight;
    nnz = 0;
    for(Size i=0;i<m;i++)
    {
        fileread >> u >> v >> weight;
        v--;
        u--;
        if(DEBUG==2)
        {
            std::cout<<u<<" "<<v<<" "<<weight<<std::endl;
        }
        if(u<v && fabs(weight)>threshold)
        {
            if(weight < 0)
               weight = (-1)*weight;
            assert(weight > 0);
            s.push_back(u);
            t.push_back(v);
            edgeWght.push_back(weight);
            nnz++;
        }
    }
}
/* Given a output file the function writes the matching into it.
 * n: # of nodes
 * outputfile: the output file
 * mateNode: vector of mates of each node
 * The output would be in the first line the number of nodes. And then n lines follows, with each line (node,mate) format with zero-indexing. If a node is unmatched the mate is -1
 */


int writeMatching(std::string outputfile, Size n, vector<Size> &mateNode)
{
    std::ofstream fileout;
    fileout.open(outputfile.c_str());

    if(fileout.is_open())
    {
        fileout<<n<<endl;
        for(int i=0;i<n;i++)
        {
            fileout<<i<<" "<<mateNode[i]<<endl;
        }
        return 0;
    }
    else
    {
        return 1;
    }
}

// DEBUG FABIO D.
int writeVec(std::string outputfile, Size n, vector<double> &mateNode)
{
    std::ofstream fileout;
    fileout.open(outputfile.c_str());

    if(fileout.is_open())
    {
        fileout<<n<<endl;
        for(int i=0;i<n;i++)
        {
            fileout<<i<<" "<<mateNode[i]<<endl;
        }
        return 0;
    }
    else
    {
        return 1;
    }
}

/* This function output the matchin g statistics */
void matchingStat(Size n, Size nnz, Amgmatch_stat &stat)
{
    cout<<"Matching Statistics"<<endl;
   cout<<"# of Nodes: "<<n<<endl;
   cout<<"# of Edges: "<<nnz<<endl;
   cout<<"Cardinality: "<<stat.matchCard<<endl;
   cout<<"Weight: "<<stat.matchWght<<endl;
   cout<<"Lambda: "<<stat.lambda<<endl;
   cout<<"Time (Sec.): "<<stat.matchTime<<endl<<endl;
}
