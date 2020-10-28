#ifndef MATCHING_ENGINE_H
#define MATCHING_ENGINE_H

#define EPS 0.0000001

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> // for fill, find, sort.
#include <numeric> // for accumulate.
#include <cassert>
#include <cmath>
#include <queue>
#include <stack>
#include <omp.h>

#include "Functional.h" // for ValGreater.
#include "Types.h"
#include "Graph.h"
#include "Utility.h" // for ResizeVector.
#include "VecItmQue.h"
#include "LstItmQue.h"


#include "VecItmStk.h"
#include "LstItmStk.h"
#include "LstItmIdxdQue.h"
#include "LstItmAndValIdxdMinPriQue.h"
#include "BhpItmAndValIdxdMinPriQue.h"
#include "LstItmAndValIdxdMaxPriQue.h"
#include "BhpItmAndValIdxdMaxPriQue.h"
#include "Timer.h"

namespace Matchbox {

    class MatchingEngine {
        public:
            /// Constructor.
            MatchingEngine():
                mInlz(false), mQueAndStkType(eVecQueAndStk), mIdxdQueType(eLstIdxdQue),
                mIdxdPriQueType(eBhpIdxdPriQue), mCardGraphSrchType(eSglSrcBfs),
                mEdgWghtGraphSrchType(eSglSrcBfs), mVtxWghtAlgType(eDrct),
                mPrcmptdMaxCard(false), mPrecision(cMaxPrecision), mFullPrint(false),
                mStatsPrint(false) {}

            /// Destructor.
            ~MatchingEngine() {}

            ///
            void ComputeMxmlMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card) const;

            ///
            void ComputeMaxCardMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card) const;

            ///
            void ComputeAprxMaxEdgWghtMatching4(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* edgWght) const;
            ///
            void ComputeMaxVtxWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght, Size augPathLenBnd) const;

            ///
            void ComputeHalfVtxWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void ComputeTwoThirdVtxWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            /////
            void SuitorVtxWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void ParSuitorVtxWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void PGVtxWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void PGDPVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ////
            void LocallyDominantVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght)  const;
            ///
            void ComputeSpecTwoThirdVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void ComputeSpecHalfVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void ComputeParSpecTwoThirdVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void ComputeParSpecHalfVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght) const;
            ///
            void ScaleOneMinusEpsVtxWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* vtxWght,
                    Val ep) const;
            ////
            void GabowCardMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card) const;
            ///
            void RandomTwoThirdMinusEpsilonEdgWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* edgWght,
                    Size iterations) const;
            ///
            void RandomOrderTwoThirdMinusEpsilonEdgWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* edgWght,
                    Size phases) const;
            ///
            void GPAEdgWght(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* edgWght) const;
            ///
            void ComputeHalfEdgWghtMatching(const Graph& graph,
                    std::vector<Size>* mateVec,
                    Size* card,
                    Val* edgWght) const;
            ///
            void nComponents(const Graph& graph) const;
            //////////
            Err GetMatchedEdgWght(const Graph& graph,
                    const std::vector<Size>& mateVec,
                    Size card,
                    Val *edgWght) const;

            ///
            Err GetMatchedVtxWght(const Graph& graph,
                    const std::vector<Size>& mateVec,
                    Size card,
                    Val *vtxWght,
                    Val *edgWght) const;

            ///
            void SetInlz(bool inlz) { mInlz = inlz; }

            ///
            Err SetQueAndStkType(QueAndStkType queAndStkType) {
                switch (queAndStkType) {
                    case eVecQueAndStk:
                    case eLstQueAndStk:
                        mQueAndStkType = queAndStkType;
                        return eErrNone;
                    default:
                        return eErrInvType;
                }
            }

            ///
            Err SetIdxdQueType(IdxdQueType idxdQueType) {
                switch (idxdQueType) {
                    case eLstIdxdQue:
                        mIdxdQueType = idxdQueType;
                        return eErrNone;
                    default:
                        return eErrInvType;
                }
            }

            ///
            Err SetIdxdPriQueType(IdxdPriQueType idxdPriQueType) {
                switch (idxdPriQueType) {
                    case eLstIdxdPriQue:
                    case eBhpIdxdPriQue:
                        mIdxdPriQueType = idxdPriQueType;
                        return eErrNone;
                    case eFhpIdxdPriQue:
                        // TODO: add fibonacci heap indexed priority queue data structure.
                    default:
                        return eErrInvType;
                }
            }

            ///
            Err SetCardGraphSrchType(GraphSrchType cardGraphSrchType) {
                switch (cardGraphSrchType) {
                    case eSglSrcBfs:
                        mCardGraphSrchType = cardGraphSrchType;
                        return eErrNone;
                    default:
                        // TODO: add other cardinality computations, especially multiple path.
                        return eErrInvType;
                }
            }

            ///
            Err SetEdgWghtGraphSrchType(GraphSrchType edgWghtGraphSrchType) {
                switch (edgWghtGraphSrchType) {
                    case eSglSrcBfs:
                        mEdgWghtGraphSrchType = edgWghtGraphSrchType;
                        return eErrNone;
                    case eMplSrcSglPath:
                        // TODO: add multiple source edge-weight computations.
                    default:
                        // TODO: reformulate the comment below in the context of more algorithms.
                        // multiple path edge-weight computations are not implemented since they
                        // are generally not expected to improve the performance, given the fact
                        // that it is generally unlikely to find a significant number of
                        // augmenting paths in the tight graph during a particular iteration.
                        return eErrInvType;
                }
            }

            ///
            Err SetVtxWghtAlgType(VtxWghtAlgType vtxWghtAlgType) {
                switch (vtxWghtAlgType) {
                    case eDrct:
                    case eSpcl:
                        mVtxWghtAlgType = vtxWghtAlgType;
                        return eErrNone;
                    default:
                        return eErrInvType;
                }
            }

            ///
            void SetPrcmptdMaxCard(bool prcmptdMaxCard) {
                mPrcmptdMaxCard = prcmptdMaxCard;
            }

            ///
            Err SetPrecision(Size precision) {
                if (precision > cMaxPrecision) {
                    return eErrInvPrecision;
                }
                mPrecision = precision;
                return eErrNone;
            }

            ///
            void SetFullPrint(bool fullPrint) { mFullPrint = fullPrint; }

            ///
            void SetStatsPrint(bool statsPrint) { mStatsPrint = statsPrint; }

            ///
            bool CheckMatching(const Graph& graph,
                    const std::vector<Size>& mateVec,
                    Size card) const;

            ///
            void PrintMatching(const std::vector<Size>& mateVec,
                    Size card) const;

            ///
            void PrintMatching(const std::vector<Size>& mateVec,
                    Size card,
                    Val wght) const;

        private:
            // controls the greedy initialization of a matching: if false then do not
            // perform greedy initialization, if true then perform greedy initialization.
            bool mInlz;
            // type of queue and stack used for vertex processing.
            QueAndStkType mQueAndStkType;
            // type of indexed queue used for vertex and weight processing.
            IdxdQueType mIdxdQueType;
            // type of indexed priority queue used for vertex and weight processing.
            IdxdPriQueType mIdxdPriQueType;
            // type of graph search used for the computation of a maximum cardinality
            // matching.
            GraphSrchType mCardGraphSrchType;
            // type of graph search used for the computation of a maximum edge-weight
            // matching.
            GraphSrchType mEdgWghtGraphSrchType;
            // type of algorithm used for the computation of a maximum vertex-weight
            // matching.
            VtxWghtAlgType mVtxWghtAlgType;
            // TODO: add comment that explains why mutable.
            // controls the precomputation of a maximum cardinality matching for the
            // speculative computation of a maximum vertex-weight matching.
            mutable bool mPrcmptdMaxCard;
            // precision used for printing floating point numbers.
            Size mPrecision;
            // controls the amount of printing: if false then print less information,
            // if true then print more information.
            bool mFullPrint;
            // controls the printing of computation statistics: if false then do not
            // print computation statistics, if true then print computation statistics.
            bool mStatsPrint;

            // forbid the copy constructor and the assignment operator.
            MatchingEngine(const MatchingEngine&);
            MatchingEngine& operator=(const MatchingEngine&);

            ///
            void rInlzGrdyForCard(const Graph& graph,
                    Size* mateArr,
                    Size* card) const;
            ///
            void rAugment(Size* mateArr,
                    Size* ptr1Arr,
                    Size* ptr2Arr,
                    Size sLast,
                    Size tLast) const;
            ///
            void rAugment(Size* mateArr,
                    Size* ptr1Arr,
                    Size* ptr2Arr,
                    Size sLast,
                    Size tLast,
                    Val* auglen) const;
            ////////
            void fix(Size* mateArr,
                    Size* ptr1Arr,
                    Size* ptr2Arr,
                    Size sLast,
                    Size tLast) const;
            ///
            Size rFind(Size* linkArr,
                    Size u) const;

            ///
            void rUnion(Size* linkArr,
                    Size* rankArr,
                    Size u,
                    Size v) const;

            ///
            void rLink(Size* linkArr,
                    Size* rankArr,
                    Size u,
                    Size v) const;

            ///
            /// without disjoint set support.
            template<class ItmQue>
                void rProcessBlsm1(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Stt* sttArr,
                        ItmQue& prcsbQue,
                        Size s1,
                        Size s2,
                        Size blsm) const;

            ///
            /// witho disjoint set support.
            template<class ItmQue>
                void rProcessBlsm2(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* linkArr,
                        Size* rankArr,
                        Stt* sttArr,
                        ItmQue& prcsbQue,
                        Size s1,
                        Size s2,
                        Size blsm/*, long double* findblsmTime, Size* findblsmCount*/) const;
            ///
            ///used by speculative MVM
            template<class ItmQue>
                void rProcessBlsm6(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* linkArr,
                        Size* rankArr,
                        Stt* sttArr,
                        ItmQue& prcsbQue,
                        Size s1,
                        Size s2,
                        Size blsm,
                        Val* minVtxWght,
                        Size* sLast,
                        Size* tLast/*,
                                     long double* findblsmTime,
                                     Size* findblsmCount*/) const;
                ////
                ///
                /// implicitly breadth-first.
                /// without disjoint set support.
                ///
                /// used in:
                /// - maximum cardinality computations.
                template<class ItmQue>
                void rFindAugPathCardSglSrc1(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Stt* sttArr,
                        ItmQue& prcsbQue,
                        ItmQue& prcsdQue,
                        Size sFirst,
                        Size* sLast,
                        Size* tLast) const;

            ///
            /// implicitly breadth-first.
            /// with disjoint set support.
            ///
            /// used in:
            /// - maximum cardinality computations.
            template<class ItmQue>
                void rFindAugPathCardSglSrc2(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* linkArr,
                        Size* rankArr,
                        Stt* sttArr,
                        ItmQue& prcsbQue,
                        ItmQue& prcsdQue,
                        Size sFirst,
                        Size* sLast,
                        Size* tLast/*, long double* findblsmTime, Size* findblsmCount,Size* ecount, Size* vcount*/) const;

            ///
            ///used by MVM
            template<class ItmQue>
                void rFindAugPathCardSglSrc3(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* linkArr,
                        Size* rankArr,
                        Stt* sttArr,
                        ItmQue& prcsbQue,
                        ItmQue& prcsdQue,
                        Size sFirst,
                        Size* sLast,
                        Size* tLast,
                        Val w2,
                        Size* sedges)const;
            ///
            ///
            ///used in restricted augmenting paths MVM
            template<class ItmQue>
                void rFindAugPathCardSglSrc1Bnd(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* ptr2Arr,
                        Size* rprsArr,
                        Stt* sttArr,
                        Size* sLenArr,
                        ItmQue& prcsbQue,
                        ItmQue& prcsdQue,
                        Size sFirst,
                        Size* sLast,
                        Size* tLast,
                        Size augPathLenBnd,
                        Size* nedges)const;
            ///
            /// the core of the single source maximum cardinality matching computation.
            /// without disjoint set support.
            template<class ItmQue>
                void rComputeMaxCardMatchingSglSrc1(const Graph& graph,
                        Size* mateArr,
                        Size* card) const;

            ///
            /// the core of the single source maximum cardinality matching computation.
            /// with disjoint set support.
            template<class ItmQue>
                void rComputeMaxCardMatchingSglSrc2(const Graph& graph,
                        Size* mateArr,
                        Size* card) const;
            ///
            ///
            /// the core of the single source maximum cardinality matching computation.
            template<class ItmQue>
                void rComputeMaxCardMatchingSglSrc3(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        std::pair<Size,Val>* sortedArr) const;
            ///
            /// the core of the single source approximate maximum weighted matching computation.
            template<class Que>
                void rComputeAprxMaxEdgWghtMatching4(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* edgWght) const;
            ///
            ///Find increasing paths it used in Specluative MVM
            template<class ItmQue> void rFindIncPathCardSglSrc3(const Graph& graph,
                    Size* mateArr,
                    Size* ptr1Arr,
                    Size* ptr2Arr,
                    Size* blsmArr,
                    Size* rprsArr,
                    Size* linkArr,
                    Size* rankArr,
                    Stt* sttArr,
                    ItmQue& prcsbQue,
                    ItmQue& prcsdQue,
                    Size sFirst,
                    Size* sLast,
                    Size* tLast) const;
            ////
            /// the core of the single source maximum vertex weighted matching computation.
            template<class ItmQue, class ItmStk>
                void rComputeMaxVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght,
                        Size augPathLenBnd) const;
            ///
            /// the core of the single source greedy 1/2-approx. vertex weighted matching computation.
            template<class ItmQue>
                void rComputeHalfVtxWghtMatching(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of the single source greedy 2/3-approx. vertex weighted matching computation.
            template<class ItmQue>
                void rComputeTwoThirdVtxWghtMatching(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of the single source Suitor 1/2-approx. vertex weighted matching computation.
            template<class ItmQue>
                void rSuitorVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            template<class ItmQue>
                /// the core of the single source locally dominant 1/2-approx. vertex weighted matching computation.
                void rLocallyDominantVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of the shared memory parallel Suitor 1/2-approx. vertex weighted matching computation.
            template<class ItmQue>
                void rParSuitorVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of the single source path growing 1/2-approx. vertex weighted matching computation.
            template<class ItmQue>
                void rPGPVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ////
            /// the core of the single source path growing and DP 1/2-approx. vertex weighted matching computation.
            template<class ItmQue>
                void rPGDPVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of the single source speculative 2/3-approx. vertex weighted matching computation.
            template<class ItmQue, class ItmStk>
                void rComputeSpecTwoThirdVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of the single source speculative 1/2-approx. vertex weighted matching computation.
            template<class ItmQue, class ItmStk>
                void rComputeSpecHalfVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of shared memory parallel speculative 2/3-approx. vertex weighted matching computation.
            template<class ItmQue, class ItmStk>
                void rComputeParSpecTwoThirdVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            /// the core of shared memory parallel speculative 1/2-approx. vertex weighted matching computation.
            template<class ItmQue, class ItmStk>
                void rComputeParSpecHalfVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght) const;
            ///
            template<class ItmQue>
                void rnComponents(const Graph& graph) const;
            ///
            ///Scaling approx. algorithm functions
            /// the core of multiple source Scaling 1-eps -approx. vertex weighted matching computation.
            template<class ItmQue, class ItmStk, class ItmIdxdQue>
                void rScaleOneMinusEpsVtxWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* vtxWght,
                        Val ep) const;
            ///
            // augment augemnting path
            // used in Scaling approximation algorithm
            template<class ItmIdxdQue,class ItmQue>
                void rAugment(Size* mateArr,
                        Size* ptrArr,
                        ItmIdxdQue& expsdQue,
                        Size sLast,
                        Size tLast,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Size num,
                        Size* label,
                        Size* blsmin) const;
            ///
            // augment through blossoms
            // used in Scaling approximation algorithm
            template<class ItmQue>
                void augmentBlsm(Size b,
                        Size s,
                        Size* mateArr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Size num) const;
            ///
            // find which blossom a vertex belong to
            // used in Scaling approximation algorithm
            template<class ItmQue>
                Size rFindblsm(Size* p,
                        Size* blsm,
                        Size* rprs,
                        Size v) const;
            ///
            // dissolve a blossom
            // used in Scaling approximation algorithm
            template<class ItmQue>
                void rDissolveBlsm3(const Graph& graph,
                        Size* blsmArr,
                        Size* blsmParentArr,
                        std::vector<Size>* blsmChildren,
                        Size b) const;
            ///
            // create and process a blossom
            // used in Scaling approximation algorithm
            template<class ItmQue,class ItmStk>
                void rProcessBlsm3(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Size* label,
                        Stt* sttArr,
                        ItmStk& prcsbStk,
                        Size s1,
                        Size s2,
                        Size* blsm,
                        Val* dualArr,
                        Size* blsmin,
                        Size* treeNumArr) const;
            ///
            // find a maximal vertex disjoint augmenting paths
            // used in Scaling approximation algorithm
            template<class ItmQue, class ItmStk, class ItmIdxdQue>
                void rFindMxmlSetAugPathsCard(const Graph& graph,
                        Size* mateArr,
                        Size* ptrArr,
                        Stt* sttArr,
                        Size* idxArr,
                        ItmStk& prcsbStk,
                        ItmQue& prcsdQue,
                        ItmIdxdQue& expsdQue,
                        ItmQue& sLastQue,
                        ItmQue& tLastQue,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Size* label,
                        Val* dualArr,
                        Size* blsm,
                        Val delta,
                        Size* blsmin,
                        Size* treeNumArr,
                        Size iter,
                        Val ep,
                        Size* se) const;
            ///
            //Gabow Max card functions
            template <class ItmQue,class ItmStk>
                void rProcessBlsm5(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Stt* sttArr,
                        ItmQue& prcsbStk,
                        Size s1,
                        Size s2,
                        Size* blsm,
                        Size* trace,
                        Size* treeNumArr,
                        int* tOuter,
                        int* Delta) const;
            ///
            template<class ItmQue>
                void rDissolveBlsm3Gabow(const Graph& graph,
                        Size* blsmArr,
                        Size* blsmParentArr,
                        std::vector<Size>* blsmChildren,
                        Size b,
                        Size blsmbdual) const;
            ///

            template<class ItmQue, class ItmStk, class ItmIdxdQue>
                void rGabowCardMatching(const Graph& graph,
                        Size* mateArr,
                        Size* card) const;
            ///
            template<class ItmQue, class ItmStk, class ItmIdxdQue>
                void rFindMxmlSetAugPathsCardGabow(const Graph& graph,
                        Size* mateArr,
                        Size* ptrArr,
                        Stt* sttArr,
                        Size* idxArr,
                        ItmStk& prcsbStk,
                        ItmIdxdQue& expsdQue,
                        ItmQue& sLastQue,
                        ItmQue& tLastQue,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        int* tOuter,
                        int* tInner,
                        int* Delta ,
                        Size* blsm,
                        Size* blsmin,
                        Size* treeNumArr) const;
            ///
            template<class ItmQue, class ItmStk, class ItmIdxdQue>
                int rEdmondSearch(const Graph& graph,
                        Size* mateArr,
                        Size* ptrArr,
                        Stt* sttArr,
                        ItmQue& searchQue,
                        ItmIdxdQue& expsdQue,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        int* tOuter,
                        int* tInner,
                        int* Delta,
                        Size* blsm,
                        Size* trace,
                        Size* treeNumArr,
                        std::deque<std::pair<Size, Size> >* deltaQArr) const;

            ///
            template <class ItmQue,class ItmStk>
                void rProcessBlsm3Gabow(const Graph& graph,
                        Size* mateArr,
                        Size* ptr1Arr,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Stt* sttArr,
                        ItmStk& prcsbStk,
                        Size s1,
                        Size s2,
                        Size* blsm,
                        Size* blsmin,
                        Size* treeNumArr) const;
            ///
            template<class ItmIdxdQue,class ItmQue>
                void rAugmentGabow(Size* mateArr,
                        Size* ptrArr,
                        ItmIdxdQue& expsdQue,
                        Size sLast,
                        Size tLast,
                        Size* blsmArr,
                        Size* rprsArr,
                        Size* blsmParent,
                        std::vector<Size>* blsmChildren,
                        std::vector<Size>* blsmL,
                        Size num,
                        Size* blsmin) const;
            ///
            template<class ItmQue>
                Size rFindblsmGabow(Size* blsm,
                        Size* rprs,
                        Size v) const;
            ///
            // find max gain 2-augmentation from an unmatched vertex
            // used in random and ordered random 2/3-eps
            template<class ItmQue>
                Size maxTwoAugUnmatched(const Graph& graph,
                        Size* mateArr,
                        Val* Mwght,
                        Val* matchingWeight,
                        Size v,
                        Size* se) const;
            ////
            // find max gain 2-augmentation from a matched vertex
            // used in random and ordered random 2/3-eps
            template<class ItmQue>
                Size maxTwoAugMatched(const Graph& graph,
                        Size* mateArr,
                        Val* markArr,
                        Val* Mwght,
                        Val* tempMwght,
                        Val* matchingWeight,
                        Size v,
                        Size* se) const;
            ///
            /// the core of the single source random 2/3-eps-approx. edge weighted matching computation.
            template<class ItmQue>
                void rRandomTwoThirdMinusEpsilonEdgWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* edgWght,
                        Size iterations) const;
            ///
            /// the core of the single source random ordered 2/3-eps-approx. edge weighted matching computation.
            template<class ItmQue>
                void rRandomOrderTwoThirdMinusEpsilonEdgWght(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* edgWght,
                        Size phases) const;
            ///
            /// the core of the single source greedy 1/2-approx. edge weighted matching computation.
            template<class ItmQue>
                void rComputeHalfEdgWghtMatching(const Graph& graph,
                        Size* mateArr,
                        Size* card,
                        Val* edgWght) const;
            ///
            // construct paths and even cycles by edges in non-increasing order of the edge weights
            // used in GPA approximation algorithm
            template<class ItmQue>
                void rGreedyPreCompute(const Graph& graph,
                        std::vector<std::deque<Size> >&,
                        std::list<std::pair<std::pair<Size,Size>,Val> > sExpsdLst2) const;
            ///
            // find an optimal matching in an even cycle by Dynamic Programming
            // used in GPA approximation algorithm
            template<class ItmQue>
                void rOptimalMatchingSingleCycle(const Graph& graph,
                        Size* mateArr,
                        Val* Mwght,
                        std::deque<Size>& cycle,
                        Val* weight) const;
            ///
            // find an optimal matching in a path by Dynamic Programming
            // used in GPA approximation algorithm
            template<class ItmQue>
                void rOptimalMatchingSinglePath(const Graph& graph,
                        Size* mateArr,
                        Val* Mwght,
                        std::deque<Size>& path,
                        Val* weight) const;
            ///
            // find optimal matchings in an collection of paths and even cycles
            // used in GPA approximation algorithm
            template<class ItmQue>
                void rOptimalMatchingColPaths(const Graph& graph,
                        Size* mateArr,
                        Val* Mwght,
                        std::vector<std::deque<Size> >& ,
                        Val* weight) const;
            ///
            /// the core of the single source GPA 1/2-eps-approx. edge weighted matching computation.
            template<class ItmQue>
                void rGreedyVariant(const Graph& graph,
                        Size* mateArr,
                        Val* Mwght,
                        Size* card,
                        Val* edgWght) const;
            ///
            Size logbase2(Val x) const;
            ///////////
    }; // class MatchingEngine
    ///

    template<class ItmQue>
        void MatchingEngine::rGreedyPreCompute(const Graph& graph, std::vector<std::deque<Size> >& pathCol,
                std::list<std::pair<std::pair<Size,Size>,Val> > sExpsdLst2) const {

            Size numVtxs =graph.mNumVtxs;
            const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            //const std::vector<Val>* edgWghtVecArr = (numVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];

            std::vector<Size> degVec;
            ResizeVector<Size>(&degVec, numVtxs);
            Size* degArr = (numVtxs == 0) ? NULL : &degVec[0];
            if (degArr != NULL) {
                std::fill(&degArr[0], &degArr[numVtxs], 0);
            }

            std::vector<Size> oppendVec;
            ResizeVector<Size>(&oppendVec, numVtxs);
            Size* oppendArr = (numVtxs == 0) ? NULL : &oppendVec[0];
            if (oppendArr != NULL) {
                std::fill(&oppendArr[0], &oppendArr[numVtxs], cNullItm);
            }

            std::vector<Size> lenVec;
            ResizeVector<Size>(&lenVec, numVtxs);
            Size* lenArr = (numVtxs == 0) ? NULL : &lenVec[0];
            if (lenArr != NULL) {
                std::fill(&lenArr[0], &lenArr[numVtxs], 0);
            }

            std::vector<std::deque<Size> > pathMember;
            ResizeVector<std::deque<Size> >(&pathMember, numVtxs);

            // iterate over the edges in non-increasing order of their weights to find eligible edges for constructng paths and even cycles
            for (std::list<std::pair<std::pair<Size,Size>, Val> >:: iterator it = sExpsdLst2.begin(); it != sExpsdLst2.end();it++) {
                std::pair<Size,Size> sedge = it->first;
                Size s = sedge.first;
                Size i = sedge.second;
                Size t = vtxVecArr[s][i];

                if(degArr[s]== 0 && degArr[t]==0)
                {
                    // non the two vertices has been considered
                    pathMember[s].push_back(i);
                    Size l = vtxVecArr[t].size();
                    for(Size j =0; j<l;j++){
                        if(vtxVecArr[t][j]==s)
                        {
                            pathMember[t].push_back(j);
                            break;
                        }

                    }

                    degArr[s]=1;
                    degArr[t]=1;
                    oppendArr[s]=t;
                    oppendArr[t]=s;
                    lenArr[s]=1;
                    lenArr[t]=1;
                }
                else if(degArr[s]== 1 && degArr[t]==1)
                {
                    if(oppendArr[t]==s)
                    {
                        // s and t are in same path
                        if(lenArr[s]==0)
                            continue;
                        else
                        {
                            pathMember[s].push_back(i);
                            Size l = vtxVecArr[t].size();

                            for(Size j =0; j<l;j++){
                                if(vtxVecArr[t][j]==s)
                                {
                                    pathMember[t].push_back(j);
                                    break;
                                }
                            }

                            oppendArr[s]=cNullItm;
                            oppendArr[t]=cNullItm;
                            lenArr[s]=0;
                            lenArr[t]=0;
                            degArr[s]++;
                            degArr[t]++;
                        }
                    }
                    else
                    {
                        // s and t are in differnt paths
                        pathMember[s].push_back(i);
                        Size l = vtxVecArr[t].size();
                        for(Size j =0; j<l;j++){
                            if(vtxVecArr[t][j]==s)
                            {
                                pathMember[t].push_back(j);
                                break;
                            }
                        }
                        Size x = oppendArr[s];
                        Size y = oppendArr[t];
                        oppendArr[x] = y;
                        oppendArr[y] = x;
                        oppendArr[s] = cNullItm;
                        oppendArr[t] = cNullItm;

                        if(	lenArr[s] > lenArr[s])
                            lenArr[x] = 1 - (lenArr[s] - lenArr[t]);
                        else
                            lenArr[x] = 1 - (lenArr[t] - lenArr[s]);

                        lenArr[y] = lenArr[x];
                        lenArr[s] = 0;
                        lenArr[t] = 0;

                        degArr[s]++;
                        degArr[t]++;
                    }
                }
                else if(degArr[s]== 1 && degArr[t]==0)
                {
                    //s is in a path and t has not been considered
                    pathMember[s].push_back(i);
                    Size l = vtxVecArr[t].size();
                    for(Size j =0; j<l;j++){
                        if(vtxVecArr[t][j]==s)
                        {
                            pathMember[t].push_back(j);
                            break;
                        }

                    }
                    oppendArr[oppendArr[s]] = t;
                    oppendArr[t] = oppendArr[s];
                    oppendArr[s] = cNullItm;

                    lenArr[t] = 1 - lenArr[s];
                    lenArr[s] = 0;

                    degArr[s]++;
                    degArr[t]++;
                }
                else if(degArr[s]== 0 && degArr[t]==1)
                {
                    //s has not been considered and t is in a path
                    pathMember[s].push_back(i);
                    Size l = vtxVecArr[t].size();
                    for(Size j =0; j<l;j++){
                        if(vtxVecArr[t][j]==s)
                        {
                            pathMember[t].push_back(j);
                            break;
                        }

                    }
                    oppendArr[oppendArr[t]] = s;
                    oppendArr[s] = oppendArr[t];
                    oppendArr[t] = cNullItm;

                    lenArr[s] = 1 - lenArr[t];
                    lenArr[t] = 0;

                    degArr[s]++;
                    degArr[t]++;
                }
            }//end for loop

            // iterate over all vertices and construct paths and even cycles
            std::deque<Size> path;
            Size succ = cNullItm;
            Size pred = cNullItm;
            for (Size v = 0; v < numVtxs; ++v)
            {
                if(degArr[v] < 1 || degArr[v]==cNullItm)
                    continue;

                path.clear();
                degArr[v]=2;
                succ = v;
                pred = cNullItm;

                while(degArr[succ]==2)
                {
                    for (Size i=0; i < pathMember[succ].size();i++)
                    {
                        Size pos = pathMember[succ][i];
                        Size t = vtxVecArr[succ][pos];
                        if(t != pred)
                        {
                            path.push_back(succ);
                            path.push_back(pos);
                            degArr[succ]=cNullItm;
                            pred = succ;
                            succ=t;
                            break;
                        }
                    }
                }
                if(degArr[succ]==cNullItm)
                {
                    pathCol.push_back(path);
                    continue;
                }
                else if(degArr[succ]==1)
                {
                    if(oppendArr[succ]==v)
                    {
                        degArr[succ]=cNullItm;
                        pathCol.push_back(path);
                        continue;
                    }
                    else
                    {
                        degArr[succ]=cNullItm;
                        for (Size i=0; i < pathMember[v].size();i++)
                        {
                            Size pos = pathMember[v][i];
                            Size t = vtxVecArr[v][pos];
                            if(degArr[t] != cNullItm)
                            {
                                path.push_front(pos);
                                path.push_front(v);
                                degArr[succ]=cNullItm;
                                pred = v;
                                succ=t;
                                break;
                            }
                        }

                        while(degArr[succ]==2)
                        {
                            for (Size i=0; i < pathMember[succ].size();i++)
                            {
                                Size pos = pathMember[succ][i];
                                Size t = vtxVecArr[succ][pos];
                                if( t !=pred )
                                {
                                    path.push_front(pos);
                                    path.push_front(succ);
                                    degArr[succ]=cNullItm;
                                    pred = succ;
                                    succ=t;
                                    break;
                                }
                            }
                        }
                        degArr[succ]=cNullItm;
                        pathCol.push_back(path);
                        continue;
                    }
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rOptimalMatchingSingleCycle(const Graph& graph, Size* mateArr,Val* Mwght, std::deque<Size>& cycle,Val* weight) const {
            /* calculates an optimal weighted matching for a cycle,
               using dynamic programming.
               */
            Size numVtxs =graph.mNumVtxs;
            const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>* edgWghtVecArr = (numVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
            Size cycleLen = cycle.size()/2;
            Val   W1[cycleLen];
            Val   W2[cycleLen];
            Size M1[cycleLen*2];
            Size M2[cycleLen*2];

            W1[0] = 0;
            W2[0] = 0;
            M1[0] = cNullItm;
            M2[0] = cNullItm;
            M1[1] = cNullItm;
            M2[1] = cNullItm;

            Size s;
            Size t;
            Val wght;
            // intial solution
            if(cycle.empty()==false)
            {
                s=cycle[0];
                t= cycle[1];
                W1[1]=edgWghtVecArr[s][cycle[1]];
                M1[2] = s;
                M1[3] = t;
            }

            s=cycle[2];
            t=cycle[3];
            W2[1]=edgWghtVecArr[s][cycle[3]];
            M2[2] = s;
            M2[3] = t;
            wght =W2[1];

            // find an optimal matching using dynamic programming
            for (Size i = 2; i < cycleLen; i++)
            {
                if (wght + W1[i-2] > W1[i-1]) {
                    W1[i] = wght + W1[i-2];
                    M1[i*2] = s;
                    M1[i*2+1] = t;
                }
                else {
                    W1[i] = W1[i-1];
                    M1[i*2] = cNullItm;
                    M1[i*2+1] = cNullItm;
                }

                s=cycle[i*2];
                t=cycle[i*2+1];
                wght = edgWghtVecArr[s][cycle[i*2+1]];
                if (wght + W2[i-2] > W2[i-1]) {
                    W2[i] = wght + W2[i-2];
                    M2[i*2] = s;
                    M2[i*2+1] = t;
                }
                else {
                    W2[i] = W2[i-1];
                    M2[i*2] =  cNullItm;
                    M2[i*2+1] =  cNullItm;
                }
            }
            // construct the matched edges from either M1 or M2
            int i = cycleLen-1;
            if (W2[i] < W1[i])
            {
                *weight =W1[i];
                while (i >= 0)
                {
                    if(cNullItm == M1[i*2])
                        i--;
                    else
                    {
                        s =M1[i*2];
                        Size pos = M1[i*2+1];
                        t= vtxVecArr[s][pos];

                        mateArr[s]=t;
                        mateArr[t]=s;
                        Mwght[s]=edgWghtVecArr[s][pos];
                        Mwght[t]=edgWghtVecArr[s][pos];
                        i-=2;
                    }
                }
            }
            else
            {
                *weight =W2[i];
                while (i >= 0)
                {
                    if (cNullItm == M2[i*2])
                        i--;
                    else
                    {
                        s =M2[i*2];
                        Size pos = M2[i*2+1];
                        t= vtxVecArr[s][pos];
                        mateArr[s]=t;
                        mateArr[t]=s;
                        Mwght[s]=edgWghtVecArr[s][pos];
                        Mwght[t]=edgWghtVecArr[s][pos];
                        i-=2;
                    }
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rOptimalMatchingSinglePath(const Graph& graph, Size* mateArr, Val* Mwght, std::deque<Size>& path,Val* weight) const {

            /* calculates an optimal weighted matching for a given path,
               using dynamic programming */
            Size numVtxs =graph.mNumVtxs;
            const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>* edgWghtVecArr = (numVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
            Size pathLen = path.size()/2;
            Val   W[pathLen+1];
            Size M[pathLen*2+2];

            W[0] = 0;
            M[0] = cNullItm;
            M[1] = cNullItm;

            Size s;
            Size t;
            Val wght;

            // initial solution
            if(path.empty()==false)
            {
                s=path[0];
                t=path[1];
                W[1]=edgWghtVecArr[s][path[1]];
                M[2] = s;
                M[3] = t;
            }

            // find an optimal matching using synamic programming
            for (Size i = 1; i < pathLen; i++)
            {

                s=path[i*2];
                t=path[i*2+1];
                wght = edgWghtVecArr[s][path[i*2 +1]];
                W[i+1] = std::max(W[i], wght + W[i-1]);
                if (std::fabs(W[i+1] - W[i] ) < 0.000000119 )
                {
                    M[i*2+2] = cNullItm;
                    M[i*2+3] = cNullItm;
                }
                else
                {
                    M[i*2+2] = s;
                    M[i*2+3] = t;
                }
            }
            *weight =W[pathLen];

            //construct the matched edges
            int i = pathLen;
            while (i > 0) {
                if ( cNullItm == M[i*2])
                    i--;

                else {
                    s  = M[i*2];
                    Size pos = M[i*2+1];
                    t= vtxVecArr[s][pos];
                    mateArr[s]=t;
                    mateArr[t]=s;
                    Mwght[s]=edgWghtVecArr[s][pos];
                    Mwght[t]=edgWghtVecArr[s][pos];
                    i-=2;
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rOptimalMatchingColPaths(const Graph& graph, Size* mateArr,Val* Mwght, std::vector<std::deque<Size> >& pathCol ,Val* weight) const {
            /* find optimal matchings for each path orcycle in pathCol*/
            Size numVtxs =graph.mNumVtxs;
            const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

            Size numPaths = pathCol.size();
            std::deque<Size> current_path;

            Size current_length;
            bool cycle = false;

            for (Size i = 0; i < numPaths; i++)
            {
                current_path   = pathCol[i];
                current_length = current_path.size();
                cycle   = false;

                Size x1 = current_path[0];
                Size x2 = vtxVecArr[x1][current_path[1]];

                Size y1 = current_path[current_length-2];
                Size y2 = vtxVecArr[y1][current_path[current_length-1]];

                // find if the path is a cycle
                if (current_length > 2)
                    if (x1 == y1 || x1 == y2 || x2 == y1 || x2 == y2) cycle = true;

                Val current_weight = 0;

                if (cycle)
                    rOptimalMatchingSingleCycle<ItmQue>(graph, mateArr,Mwght, current_path,&current_weight);
                else
                    rOptimalMatchingSinglePath<ItmQue>(graph, mateArr,Mwght, current_path,&current_weight);

                // update the weight
                (*weight)+=current_weight;
            }
        }

    template<class ItmQue>
        void MatchingEngine::rGreedyVariant(const Graph& graph, Size* mateArr, Val* Mwght, Size* card, Val* edgWght) const {

            Size numVtxs =graph.mNumVtxs;
            //Size numEdgs =graph.mNumEdgs;

            if (mateArr != NULL) {
                std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
            }

            const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>* edgWghtVecArr = (numVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];

            std::vector<std::deque<Size> > pathCol;
            std::list<std::pair<std::pair<Size,Size>, Val> > sExpsdLst2;
            // a list of edges and their weights
            for (Size sFirst = 0; sFirst < numVtxs; ++sFirst) {
                Size numEdgs = vtxVecArr[sFirst].size();
                const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
                const Val* edgWghtArr = (numEdgs == 0) ? NULL : &edgWghtVecArr[sFirst][0];
                for (Size i = 0; i < numEdgs; ++i) {
                    Size t = vtxArr[i];
                    Val wgt = edgWghtArr[i];
                    if(sFirst>t)
                        continue;
                    std::pair<Size,Size> p (sFirst,i);
                    sExpsdLst2.push_back(std::pair<std::pair<Size,Size>, Val>(p, wgt));
                }
            }
            // sort the edges
            sExpsdLst2.sort(ValGreaterPairs< std::pair< std::pair<Size,Size>, Val> >());
            Val weight = 0;
            Size rounds = 0;

            std::vector<std::deque<Size> > paths;
            Val current_weight =0;

            while(sExpsdLst2.empty()==false)
            {
                rounds++;
                paths.clear();
                current_weight = 0;
                // find a collection of paths and even cycles
                rGreedyPreCompute<ItmQue>(graph,paths,sExpsdLst2);
                // find the optimal matching in the collection of paths and even cycles
                rOptimalMatchingColPaths<ItmQue>(graph,mateArr,Mwght, paths,&current_weight);
                weight += current_weight;
                // delete matched and adjacent edges
                for (std::list<std::pair<std::pair<Size,Size>, Val> >:: iterator it = sExpsdLst2.begin(); it != sExpsdLst2.end();) {
                    Size s = (it->first).first;
                    Size t = vtxVecArr[s][(it->first).second];
                    if(mateArr[s]!=cNullItm || mateArr[t]!=cNullItm)
                        it=sExpsdLst2.erase(it);
                    else
                        ++it;
                }
            }

            *edgWght = weight;
            for (Size s = 0; s < numVtxs; ++s) {
                if (mateArr[s] != cNullItm) {
                    (*card)++;
                }
            }
            *card=*card/2;
        }


    template<class ItmQue>
        void MatchingEngine::rComputeHalfEdgWghtMatching(
                const Graph& graph, Size* mateArr,Size* card, Val* edgWght) const {

            Size numVtxs =graph.mNumVtxs;
            Size numEdgs =graph.mNumEdgs;
            const std::vector<Size>* vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>* edgWghtVecArr = (numVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
            std::vector<std::pair<std::pair<Size,Size>, Val> > sExpsdLst2;
            ResizeVector<std::pair<std::pair<Size,Size>, Val> >(&sExpsdLst2, numEdgs);
            std::pair<std::pair<Size,Size>,Val> * sExpsdLst2Arr=  (numEdgs == 0) ? NULL : &sExpsdLst2[0];
            if (mateArr != NULL) {
                std::fill(&mateArr[0], &mateArr[numVtxs], cNullItm);
            }

            *card = 0;

            Size edgeCounter =0;
            // list of edges and their weights
            for (Size sFirst = 0; sFirst < numVtxs; ++sFirst) {
                Size numEdgs = vtxVecArr[sFirst].size();
                const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
                const Val* edgWghtArr = (numEdgs == 0) ? NULL : &edgWghtVecArr[sFirst][0];
                for (Size i = 0; i < numEdgs; ++i) {
                    Size t = vtxArr[i];
                    Val wgt = edgWghtArr[i];
                    if(sFirst>t)
                        continue;
                    std::pair<Size,Size> p (sFirst,t);
                    sExpsdLst2Arr[edgeCounter]=std::pair<std::pair<Size,Size>, Val>(p, wgt);
                    edgeCounter++;
                }

            }
            // sort the edges
            std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreaterPairs<std::pair<std::pair<Size,Size>, Val> >());

            std::pair<Size,Size> sedge;

            for (Size current = 0; current < numEdgs; ++current) {
                sedge =  sExpsdLst2Arr[current].first;
                Size s = sedge.first;
                Size t = sedge.second;
                if(mateArr[s]!=cNullItm || mateArr[t]!=cNullItm)
                    continue;
                else
                {
                    // match an edge and udate the cardinality the weight
                    mateArr[t]=s;
                    mateArr[s]=t;
                    ++(*card);
                    *edgWght += sExpsdLst2Arr[current].second;
                }
            }
        }

    template<class ItmQue>
        Size MatchingEngine::maxTwoAugUnmatched(const Graph& graph, Size* mateArr,Val* Mwght, Val* matchingWeight, Size v,Size* se) const {
            Val maxgain=0;
            Size n1=cNullItm;
            Val n1wght=0;
            const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>* edgWghtVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];

            Size sNumEdgs = vtxVecArr[v].size();
            const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
            const Val* edgWghtArr = (sNumEdgs == 0) ? NULL : &edgWghtVecArr[v][0];
            // find the max gain from an unmatched vertex v
            for (Size i=0; i < sNumEdgs;i++)
            {
                Size u = vtxArr[i];
                Val wght = edgWghtArr[i];
                //(*se)++;

                //augmenting path of length 1
                if(mateArr[u]==cNullItm)
                {
                    // if u is unmathed
                    if(wght > maxgain)
                    {
                        maxgain = wght;
                        n1wght= wght;
                        n1=u;
                    }
                }
                else
                {
                    // if u is mathed
                    if(wght - Mwght[u] > maxgain)
                    {
                        maxgain = wght - Mwght[u];
                        n1wght= wght;
                        n1=u;
                    }
                }
            }
            if(n1!=cNullItm)
            {

                *matchingWeight = *matchingWeight + maxgain;
                //if the found vertex is unmatched
                if(mateArr[n1]==cNullItm)
                {
                    mateArr[n1]=v;
                    mateArr[v]=n1;
                    Mwght[n1]=n1wght;
                    Mwght[v]=n1wght;
                    return 1;
                }
                else
                {
                    Size tempN =mateArr[n1];
                    mateArr[n1]=v;
                    mateArr[v]=n1;
                    mateArr[tempN]=cNullItm;
                    Mwght[n1]=n1wght;
                    Mwght[v]=n1wght;
                    Mwght[tempN]=0;
                    return 1;
                }
            }
            else
                return 0;
        }

    template<class ItmQue>
        Size MatchingEngine::maxTwoAugMatched(const Graph& graph, Size* mateArr,Val* markArr,Val* Mwght,Val* tempMwght, Val* matchingWeight, Size v,Size*se) const {
            Val maxGain=0;
            Val currentGain=0;
            Size vMate =mateArr[v];
            std::deque<Size> aug;
            std::deque<Size> arm;
            std::deque<Size> p1;
            std::deque<Size> p2;
            std::deque<Size> markedVtcs;
            Val gainP1=0;
            Val gainP2=0;
            Val gainArm=0;
            Val vuP1wght =0;
            Val vuP2wght =0;
            Size vuP1 =0;
            Size vuP2 =0;

            const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>* edgWghtVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];

            Size sNumEdgs1 = vtxVecArr[v].size();
            const Size* vtxArr1 = (sNumEdgs1 == 0) ? NULL : &vtxVecArr[v][0];
            const Val* edgWghtArr1 = (sNumEdgs1 == 0) ? NULL : &edgWghtVecArr[v][0];
            // find the max gain augmentation of a cycle of length 4
            for (Size i=0; i < sNumEdgs1;i++)
            {
                Size u = vtxArr1[i];
                //(*se)++;
                if(u == vMate)
                    continue;
                markArr[u]=edgWghtArr1[i];
                markedVtcs.push_back(u);
            }

            Size sNumEdgs2 = vtxVecArr[vMate].size();
            const Size* vtxArr2 = (sNumEdgs2 == 0) ? NULL : &vtxVecArr[vMate][0];
            const Val* edgWghtArr2 = (sNumEdgs2 == 0) ? NULL : &edgWghtVecArr[vMate][0];
            for (Size i=0; i < sNumEdgs2;i++)
            {
                Size u = vtxArr2[i];
                Val wght = edgWghtArr2[i];
                //(*se)++;
                if(u == v)
                    continue;
                if(mateArr[u]!= cNullItm && markArr[mateArr[u]]> 0.000000119)
                {
                    if(wght + markArr[mateArr[u]]- Mwght[u]- Mwght[v] > maxGain)
                    {
                        maxGain = wght+markArr[mateArr[u]]- Mwght[u]- Mwght[v];
                        aug.clear();
                        aug.push_back(v);
                        aug.push_back(vMate);
                        aug.push_back(vMate);
                        aug.push_back(u);
                        tempMwght[vMate]=wght;
                        tempMwght[u]=wght;
                        aug.push_back(u);
                        aug.push_back(mateArr[u]);
                        aug.push_back(mateArr[u]);
                        aug.push_back(v);
                        tempMwght[mateArr[u]]=markArr[mateArr[u]];
                        tempMwght[v]=markArr[mateArr[u]];
                    }
                }
            }

            //pop out queues to reset mark array
            while (markedVtcs.empty() == false)
            {
                Size u = markedVtcs.front();
                markedVtcs.pop_front();
                markArr[u]=0.0;
            }

            // find the max gain augmentation of arms P1 and P2 if any
            for (Size i=0; i < sNumEdgs1;i++)
            {
                Size u = vtxArr1[i];
                //(*se)++;
                Val wght = edgWghtArr1[i];
                if(u == vMate)
                    continue;
                gainArm = 0;
                arm.clear();
                if(mateArr[u]!=cNullItm)
                {
                    gainArm=wght-Mwght[u];
                    arm.push_back(mateArr[u]);
                    arm.push_back(u);
                    arm.push_back(u);
                    arm.push_back(v);
                }
                else
                {
                    gainArm=wght;
                    arm.push_back(u);
                    arm.push_back(v);
                }
                if(gainArm > gainP2)
                {
                    if(gainArm>gainP1)
                    {
                        if(p1.empty()==false)
                        {
                            if((mateArr[u]!=cNullItm && vuP1 != mateArr[u] )|| mateArr[u]==cNullItm)
                            {
                                gainP2=gainP1;
                                vuP2wght=vuP1wght;
                                vuP2= vuP1;
                                p2.clear();
                                p2=p1;
                            }
                        }
                        gainP1=gainArm;
                        vuP1wght=wght;
                        vuP1=u;
                        p1.clear();
                        p1=arm;
                    }
                    else
                    {
                        if((mateArr[u]!=cNullItm && vuP1 != mateArr[u]) || mateArr[u]==cNullItm)
                        {
                            gainP2=gainArm;
                            vuP2wght=wght;
                            vuP2=u;
                            p2.clear();
                            p2=arm;
                        }
                    }
                }
            }

            // Now check every arm of vMate
            if(sNumEdgs2 ==1)
            {
                currentGain = gainP1 - Mwght[v];
                if(currentGain > maxGain)
                {
                    maxGain=currentGain;
                    aug.clear();
                    aug =p1;
                    aug.push_back(v);
                    aug.push_back(vMate);
                    tempMwght[v]=vuP1wght;
                    tempMwght[vuP1]=vuP1wght;
                }
            }
            else
            {
                // degree vMate is > 1
                if(p1.empty())
                {
                    for (Size i=0; i < sNumEdgs2;i++)
                    {
                        Size u = vtxArr2[i];
                        //(*se)++;
                        Val wght = edgWghtArr2[i];
                        if(u == v)
                            continue;
                        gainArm=0;
                        arm.clear();
                        arm.push_back(v);
                        arm.push_back(vMate);
                        arm.push_back(vMate);
                        arm.push_back(u);
                        if(mateArr[u]== cNullItm)
                        {
                            gainArm=wght-Mwght[v];
                        }
                        else
                        {
                            gainArm=wght-Mwght[v]-Mwght[u];
                            arm.push_back(u);
                            arm.push_back(mateArr[u]);
                        }
                        if(gainArm > maxGain)
                        {
                            maxGain = gainArm;
                            tempMwght[vMate]=wght;
                            tempMwght[u]=wght;
                            aug.clear();
                            aug = arm;
                        }
                    }
                }
                else
                {
                    if(p1.size()==2)
                    {
                        for (Size i=0; i < sNumEdgs2;i++)
                        {
                            Size u = vtxArr2[i];
                            //(*se)++;
                            Val wght = edgWghtArr2[i];
                            if(u == v)
                                continue;
                            gainArm=0;
                            currentGain=0;
                            arm.clear();
                            Size whichP=0;
                            if (std::find(p1.begin(), p1.end(), u) != p1.end())
                            {
                                if(p2.empty()==false)
                                {
                                    currentGain = gainP2 - Mwght[v] + wght;
                                    arm = p2;
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                    whichP=2;
                                }
                                else
                                {
                                    currentGain = wght - Mwght[v];
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                    whichP=0;
                                }
                            }
                            else
                            {
                                if(mateArr[u]== cNullItm)
                                {
                                    currentGain = gainP1 - Mwght[v] + wght;
                                    arm = p1;
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                }
                                else
                                {
                                    currentGain = gainP1 - Mwght[v] + wght - Mwght[u];
                                    arm = p1;
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                    arm.push_back(u);
                                    arm.push_back(mateArr[u]);
                                }
                                whichP=1;
                            }

                            if(currentGain > maxGain)
                            {
                                if(whichP==1)
                                {
                                    tempMwght[v]=vuP1wght;
                                    tempMwght[vuP1]=vuP1wght;
                                }
                                else if(whichP==2)
                                {
                                    tempMwght[v]=vuP2wght;
                                    tempMwght[vuP2]=vuP2wght;
                                }

                                tempMwght[vMate]=wght;
                                tempMwght[u]=wght;
                                maxGain = currentGain;
                                aug.clear();
                                aug = arm;
                            }
                        }
                    }
                    else
                    {
                        for (Size i=0; i < sNumEdgs2;i++)
                        {
                            Size u = vtxArr2[i];
                            //(*se)++;
                            Val wght = edgWghtArr2[i];
                            if(u == v)
                                continue;
                            gainArm=0;
                            currentGain=0;
                            arm.clear();
                            Size whichP=0;
                            if (std::find(p1.begin(), p1.end(), u) != p1.end())
                            {
                                if(p2.empty()==false)
                                {
                                    currentGain = gainP2 - Mwght[v] + wght- Mwght[u];
                                    arm = p2;
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                    arm.push_back(u);
                                    arm.push_back(mateArr[u]);
                                    whichP=2;
                                }
                                else
                                {
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                    if(mateArr[u]== cNullItm)
                                    {
                                        currentGain = wght - Mwght[v];
                                    }
                                    else
                                    {
                                        currentGain=wght-Mwght[v]-Mwght[u];
                                        arm.push_back(u);
                                        arm.push_back(mateArr[u]);
                                    }
                                    whichP=0;
                                }
                            }
                            else
                            {
                                if(mateArr[u]== cNullItm)
                                {
                                    currentGain = gainP1 - Mwght[v] + wght;
                                    arm = p1;
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                }
                                else
                                {
                                    currentGain = gainP1 - Mwght[v] + wght - Mwght[u];
                                    arm = p1;
                                    arm.push_back(v);
                                    arm.push_back(vMate);
                                    arm.push_back(vMate);
                                    arm.push_back(u);
                                    arm.push_back(u);
                                    arm.push_back(mateArr[u]);
                                }
                                whichP=1;
                            }

                            if(currentGain>maxGain)
                            {
                                maxGain = currentGain;
                                if(whichP==1)
                                {
                                    tempMwght[v]=vuP1wght;
                                    tempMwght[vuP1]=vuP1wght;
                                }
                                else if(whichP==2)
                                {
                                    tempMwght[v]=vuP2wght;
                                    tempMwght[vuP2]=vuP2wght;
                                }
                                tempMwght[vMate]=wght;
                                tempMwght[u]=wght;
                                aug.clear();
                                aug = arm;
                            }
                        }
                    }
                }
            }
            if (0.000000119 >= maxGain) return 0;

            // update the matching weight
            *matchingWeight = *matchingWeight + maxGain;

            //update the matched edges
            Size u1 = aug[0];
            Size u2 = aug[1];
            if(mateArr[u1]== u2)
            {
                mateArr[u1]=cNullItm;
                mateArr[u2]=cNullItm;
                Mwght[u1]=0;
                Mwght[u2]=0;
                aug.pop_front();
                aug.pop_front();
            }
            while(aug.empty() == false)
            {
                u1 = aug.front();
                aug.pop_front();
                u2 = aug.front();
                aug.pop_front();
                if(aug.empty() == false)
                {
                    Size tempU1= aug.front();
                    aug.pop_front();
                    Size tempU2= aug.front();
                    aug.pop_front();
                    mateArr[tempU1]=cNullItm;
                    mateArr[tempU2]=cNullItm;
                    Mwght[tempU1]=0;
                    Mwght[tempU2]=0;
                }
                mateArr[u1]=u2;
                mateArr[u2]=u1;
                Mwght[u1]=tempMwght[u1];
                Mwght[u2]=tempMwght[u2];
            }
            return 1;
        }

    template<class ItmQue>
        void MatchingEngine::rRandomTwoThirdMinusEpsilonEdgWght(
                const Graph& graph, Size* mateArr, Size* card,
                Val* edgWght,Size iterations) const {

            Size numVtxs =graph.mNumVtxs;
            if (mateArr != NULL) {
                std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
            }

            *card = 0;
            *edgWght=0;
            Val tempedgeWght=0;
            Size se=0;
            //	const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

            // to store matched edges weight
            std::vector<Val> MwghtVec;
            ResizeVector<Val>(&MwghtVec, graph.mNumVtxs);
            Val* Mwght = (graph.mNumVtxs == 0) ? NULL : &MwghtVec[0];
            if (Mwght != NULL) {
                std::fill(&Mwght[0], &Mwght[graph.mNumVtxs], 0.0);
            }
            // to store temporary matched edges weight
            std::vector<Val> tempMwghtVec;
            ResizeVector<Val>(&tempMwghtVec, graph.mNumVtxs);
            Val* tempMwght = (graph.mNumVtxs == 0) ? NULL : &tempMwghtVec[0];
            if (tempMwght != NULL) {
                std::fill(&tempMwght[0], &tempMwght[graph.mNumVtxs], 0.0);
            }
            // mark vertices to find cycles of length 4
            std::vector<Val> markVec;
            ResizeVector<Val>(&markVec, graph.mNumVtxs);
            Val* markArr = (graph.mNumVtxs == 0) ? NULL : &markVec[0];
            if (markArr != NULL) {
                std::fill(&markArr[0], &markArr[graph.mNumVtxs], 0.0);
            }
            srand (time(NULL));
            for (Size i = 0; i < iterations; ++i) {
                for (Size j = 0; j < numVtxs; j++) {
                    Size v = rand() % numVtxs; // pick a vertex at random
                    if (mateArr[v] != cNullItm)
                        maxTwoAugMatched<ItmQue>(graph, mateArr,markArr,Mwght,tempMwght, &tempedgeWght,v,&se);
                    else
                        maxTwoAugUnmatched<ItmQue>(graph, mateArr,Mwght,&tempedgeWght, v,&se);
                }
            }
            //timer.Stop();
            // timer.Print();
            *edgWght = tempedgeWght;

            std::cout << se<<" ";
            for (Size s = 0; s < numVtxs; ++s) {
                if (mateArr[s] != cNullItm) {
                    (*card)++;
                }
            }
            *card=*card/2;
        }

    template<class ItmQue>
        void MatchingEngine::rRandomOrderTwoThirdMinusEpsilonEdgWght(
                const Graph& graph, Size* mateArr, Size* card,
                Val* edgWght,Size phases) const {
            Size numVtxs =graph.mNumVtxs;

            if (mateArr != NULL) {
                std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
            }

            *edgWght=0;
            Val tempedgeWght =0;
            Val previousWght =0;
            Size se =0;
            //const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            // to store matched edges weight
            std::vector<Val> MwghtVec;
            ResizeVector<Val>(&MwghtVec, graph.mNumVtxs);
            Val* Mwght = (graph.mNumVtxs == 0) ? NULL : &MwghtVec[0];
            if (Mwght != NULL) {
                std::fill(&Mwght[0], &Mwght[graph.mNumVtxs], 0.0);
            }
            // to store temporary matched edges weight
            std::vector<Val> tempMwghtVec;
            ResizeVector<Val>(&tempMwghtVec, graph.mNumVtxs);
            Val* tempMwght = (graph.mNumVtxs == 0) ? NULL : &tempMwghtVec[0];
            if (tempMwght != NULL) {
                std::fill(&tempMwght[0], &tempMwght[graph.mNumVtxs], 0.0);
            }
            // mark vertices to find cycles of length 4
            std::vector<Val> markVec;
            ResizeVector<Val>(&markVec, graph.mNumVtxs);
            Val* markArr = (graph.mNumVtxs == 0) ? NULL : &markVec[0];
            if (markArr != NULL) {
                std::fill(&markArr[0], &markArr[graph.mNumVtxs], 0.0);
            }

            std::vector<Size>VtxVec;
            ResizeVector<Size>(&VtxVec, graph.mNumVtxs);
            Size* vtxArr = (graph.mNumVtxs == 0) ? NULL : &VtxVec[0];
            if (vtxArr != NULL) {
                for(Size i =0; i< graph.mNumVtxs; i++) vtxArr[i]=i;
            }

            //rGreedyVariant<ItmQue>(graph, mateArr,Mwght, card, &tempedgeWght);

            *card=0;
            //Timer timer;
            // timer.Start();

            for (Size i = 0; i < phases; ++i) {
                std::random_shuffle( VtxVec.begin(), VtxVec.end()); // permute the list of vertices
                previousWght=tempedgeWght;
                for (Size j = 0; j < numVtxs; j++) {
                    Size v =vtxArr[j];

                    //if already matched
                    if (mateArr[v] != cNullItm)
                        maxTwoAugMatched<ItmQue>(graph, mateArr,markArr,Mwght,tempMwght,&tempedgeWght,v,&se);
                    else
                        maxTwoAugUnmatched<ItmQue>(graph, mateArr,Mwght,&tempedgeWght, v,&se);
                }
                if(previousWght == tempedgeWght)
                {
                    // exit the loop if there is no improvement
                    break;
                }
            }
            //timer.Stop();
            // timer.Print();
            //Val ew=0;
            *edgWght = tempedgeWght;
            for (Size s = 0; s < numVtxs; ++s) {
                if (mateArr[s] != cNullItm) {
                    //ew+=Mwght[s];
                    (*card)++;
                }
            }
            //ew=ew/2;
            /*std::cout<<ew<<std::endl;
            if(tempedgeWght==ew)
                std::cout << "equal "<<std::endl;
            else
                std::cout << "not equal "<<std::endl;
            */
            *card=*card/2;
        }

    template<class ItmQue>
        void MatchingEngine::rProcessBlsm1(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr,
                Stt* sttArr, ItmQue& prcsbQue, Size s1, Size s2, Size blsm) const {
            // find vertex rprs, the representative of the blossom. search alternatively
            // along both paths toward the base of the blossom, instead of searching
            // along one path first and then along the other one, in order for the number
            // of steps in the search to be proportional to the number of vertices in the
            // blossom rather than being bounded by the number of vertices in the graph.
            // for a graph G = (V,E) the former approach would determine a complexity of
            // O(|V|) per augmenting path search, as opposed to the O(|V|^2) complexity
            // determined by the latter. the number of blossoms discovered during an
            // augmenting path search is O(|V|). the number of vertices "shrunk" into
            // blossoms during an augmenting path search is also O(|V|).
            Size rprs(cNullItm);
            Size tmpRprsArr[2];
            tmpRprsArr[0] = rprsArr[s1];
            tmpRprsArr[1] = rprsArr[s2];
            blsmArr[tmpRprsArr[0]] = blsm;
            blsmArr[tmpRprsArr[1]] = blsm;
            bool done = false;
            do {
                for (Size k = 0; k < 2; ++k) {
                    Size r = tmpRprsArr[k];
                    Size s = ptr1Arr[r];
                    if (s != cNullItm) {
                        r = rprsArr[s];
                        if (blsmArr[r] == blsm) {
                            rprs = r;
                            done = true;
                            break;
                        }
                        blsmArr[r] = blsm;
                        tmpRprsArr[k] = r;
                    }
                }
            } while (done == false);
            // make every T-vertex in the blossom an S-vertex and make vertex rprs the
            // representative of such a T-vertex that becomes an S-vertex. this time
            // separately along both paths toward the base of the blossom since it is
            // easy now to prevent the procedure from continuing beyond the base. for a
            // graph G = (V,E), for the reasons stated above, the complexity of this task
            // is O(|V|) per augmenting path search.
            tmpRprsArr[0] = rprsArr[s1];
            tmpRprsArr[1] = rprsArr[s2];
            for (Size k = 0; k < 2; ++k) {
                Size r = tmpRprsArr[k];
                while (r != rprs) {
                    Size t = mateArr[r];
                    ptr1Arr[t] = s1;
                    ptr2Arr[t] = s2;
                    rprsArr[t] = rprs;
                    prcsbQue.Push(t);
                    sttArr[t] = eSttBfsPrcsb;
                    Size s = ptr1Arr[r];
                    r = rprsArr[s];
                }
            }
            // update the representative of every S-vertex (original S-vertex or T-vertex
            // that later became S-vertex) whose representative is supposed to become
            // rprs, the representative of the new blossom. this is the most
            // expensive part of the blossom processing task because it visits every
            // vertex in the graph. since, for a graph G = (V,E), the number of blossoms
            // discovered during an augmenting path search is O(|V|), this determines a
            // complexity of O(|V|^2) per augmenting path search.
            for (Size s = 0; s < graph.mNumVtxs; ++s) {
                Size r = rprsArr[s];
                // check if rprs is already the representative of s or if r was involved
                // in the task of determining rprs.
                if ((r == rprs) || (blsmArr[r] != blsm)) {
                    continue;
                }
                Size t = mateArr[r];
                // check if r is actually inside the new blossom. it could have been
                // involved in the task of determining rprs, but it could still be
                // outside the blossom (beyond the blossom's base).
                if ((t == cNullItm) || (rprsArr[t] != rprs)) {
                    continue;
                }
                rprsArr[s] = rprs;
            }
        }

    template<class ItmQue>
        void MatchingEngine::rProcessBlsm2(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
                Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, Size s1, Size s2, Size blsm)
        const {
            // find vertex rprs, the representative of the blossom. search alternatively
            // along both paths toward the base of the blossom, instead of searching
            // along one path first and then along the other one, in order for the number
            // of steps in the search to be proportional to the number of vertices in the
            // blossom rather than being bounded by the number of vertices in the graph.
            // for a graph G = (V,E) the former approach would determine a complexity of
            // O(|V|) per augmenting path search, as opposed to the O(|V|^2) complexity
            // determined by the latter. the number of blossoms discovered during an
            // augmenting path search is O(|V|). the number of vertices "shrunk" into
            // blossoms during an augmenting path search is also O(|V|).
            Size rprs(cNullItm);
            Size tmpRprsArr[2];
            //Timer timer;
            //timer.Start();
            tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
            tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
            //timer.Stop();
            //(*findblsmTime)+=timer.GetTime();
            //(*findblsmCount)+=2;
            blsmArr[tmpRprsArr[0]] = blsm;
            blsmArr[tmpRprsArr[1]] = blsm;
            bool done = false;
            do {
                for (Size k = 0; k < 2; ++k) {
                    Size r = tmpRprsArr[k];
                    Size s = ptr1Arr[r];
                    if (s != cNullItm) {
                        //	timer.Start();
                        r = rprsArr[rFind(linkArr, s)];
                        //  timer.Stop();
                        // (*findblsmTime)+=timer.GetTime();
                        // (*findblsmCount)++;
                        if (blsmArr[r] == blsm) {
                            rprs = r;
                            done = true;
                            break;
                        }
                        blsmArr[r] = blsm;
                        tmpRprsArr[k] = r;
                    }
                }
            } while (done == false);
            // make every T-vertex in the blossom an S-vertex and make vertex rprs the
            // representative of such a T-vertex that becomes an S-vertex. this time
            // separately along both paths toward the base of the blossom since it is
            // easy now to prevent the procedure from continuing beyond the base. for a
            // graph G = (V,E), for the reasons stated above, the complexity of this task
            // is O(|V|) per augmenting path search.
            //timer.Start();
            tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
            tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
            //timer.Stop();
            //(*findblsmTime)+=timer.GetTime();
            //(*findblsmCount)+=2;
            for (Size k = 0; k < 2; ++k) {
                Size r = tmpRprsArr[k];
                while (r != rprs) {
                    //rprsArr[rFind(linkArr, r)] = rprs;
                    rUnion(linkArr, rankArr, r, rprs);
                    Size t = mateArr[r];
                    ptr1Arr[t] = s1;
                    ptr2Arr[t] = s2;
                    rUnion(linkArr, rankArr, t, rprs);
                    prcsbQue.Push(t);
                    sttArr[t] = eSttBfsPrcsb;

                    //timer.Start();
                    rprsArr[rFind(linkArr, rprs)] = rprs;
                    // timer.Stop();
                    Size s = ptr1Arr[r];
                    //    timer.Start();
                    r = rprsArr[rFind(linkArr, s)];
                    //    timer.Stop();
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rProcessBlsm6(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
                Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, Size s1, Size s2, Size blsm, Val* minVtxWght, Size* sLast, Size* tLast)
        const {
            // find vertex rprs, the representative of the blossom. search alternatively
            // along both paths toward the base of the blossom, instead of searching
            // along one path first and then along the other one, in order for the number
            // of steps in the search to be proportional to the number of vertices in the
            // blossom rather than being bounded by the number of vertices in the graph.
            // for a graph G = (V,E) the former approach would determine a complexity of
            // O(|V|) per augmenting path search, as opposed to the O(|V|^2) complexity
            // determined by the latter. the number of blossoms discovered during an
            // augmenting path search is O(|V|). the number of vertices "shrunk" into
            // blossoms during an augmenting path search is also O(|V|).
            //Also return the lightest vertex in the blossom which is used by the speculative MVM
            const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
            Size rprs(cNullItm);
            Size tmpRprsArr[2];
            //Timer timer;
            //timer.Start();
            tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
            tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
            //timer.Stop();
            blsmArr[tmpRprsArr[0]] = blsm;
            blsmArr[tmpRprsArr[1]] = blsm;
            bool done = false;
            do {
                for (Size k = 0; k < 2; ++k) {
                    Size r = tmpRprsArr[k];
                    Size s = ptr1Arr[r];
                    if (s != cNullItm) {
                        //	timer.Start();
                        r = rprsArr[rFind(linkArr, s)];
                        //  timer.Stop();
                        if (blsmArr[r] == blsm) {
                            rprs = r;
                            done = true;
                            break;
                        }
                        blsmArr[r] = blsm;
                        tmpRprsArr[k] = r;
                    }
                }
            } while (done == false);
            // make every T-vertex in the blossom an S-vertex and make vertex rprs the
            // representative of such a T-vertex that becomes an S-vertex. this time
            // separately along both paths toward the base of the blossom since it is
            // easy now to prevent the procedure from continuing beyond the base. for a
            // graph G = (V,E), for the reasons stated above, the complexity of this task
            // is O(|V|) per augmenting path search.
            //timer.Start();
            tmpRprsArr[0] = rprsArr[rFind(linkArr, s1)];
            tmpRprsArr[1] = rprsArr[rFind(linkArr, s2)];
            //timer.Stop();
            //(*findblsmTime)+=timer.GetTime();
            //(*findblsmCount)+=2;
            for (Size k = 0; k < 2; ++k) {
                Size r = tmpRprsArr[k];
                while (r != rprs) {
                    //rprsArr[rFind(linkArr, r)] = rprs;
                    rUnion(linkArr, rankArr, r, rprs);
                    Size t = mateArr[r];
                    ptr1Arr[t] = s1;
                    ptr2Arr[t] = s2;
                    rUnion(linkArr, rankArr, t, rprs);
                    prcsbQue.Push(t);
                    sttArr[t] = eSttBfsPrcsb;
                    if (*minVtxWght > vtxWghtArr[t]) {
                        *sLast = r;
                        *tLast = r;
                        *minVtxWght = vtxWghtArr[t];
                    }

                    //timer.Start();
                    rprsArr[rFind(linkArr, rprs)] = rprs;
                    // timer.Stop();
                    Size s = ptr1Arr[r];
                    //    timer.Start();
                    r = rprsArr[rFind(linkArr, s)];
                    //    timer.Stop();
                }
            }
        }



    template<class ItmQue>
        void MatchingEngine::rFindAugPathCardSglSrc1(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Stt* sttArr,
                ItmQue& prcsbQue, ItmQue& prcsdQue, Size sFirst, Size* sLast, Size* tLast)
        const {
            *sLast = cNullItm;
            *tLast = cNullItm;
            const std::vector<Size>*
                vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            Size blsm = 0;
            prcsbQue.Push(sFirst);
            sttArr[sFirst] = eSttBfsPrcsb;
            while (prcsbQue.Empty() == false) {
                Size s = prcsbQue.Front();
                prcsbQue.Pop();
                prcsdQue.Push(s);
                sttArr[s] = eSttBfsPrcsd;
                Size sNumEdgs = vtxVecArr[s].size();
                const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                for (Size i = 0; i < sNumEdgs; ++i) {
                    Size t = sVtxArr[i];
                    assert(s != t);
                    // if t is an outer vertex and s and t are in different blossoms then creat a blossom
                    if (sttArr[t] == eSttBfsPrcsb) {
                        if (rprsArr[s] != rprsArr[t]) {
                            rProcessBlsm1(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                                    sttArr, prcsbQue, s, t, blsm);
                        }
                        ++blsm;
                        continue;
                    } else if (sttArr[t] == eSttBfsPrcsd) {
                        continue;
                    }
                    Size ss = mateArr[t];
                    if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
                        continue;
                    }
                    if (ss == cNullItm) {
                        // an unmatched vertex is found
                        *sLast = s;
                        *tLast = t;
                        return;
                    }
                    if (sttArr[ss] != eSttIdle) {
                        continue;
                    }
                    ptr1Arr[ss] = s;
                    prcsbQue.Push(ss);
                    sttArr[ss] = eSttBfsPrcsb;
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rFindAugPathCardSglSrc2(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
                Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, ItmQue& prcsdQue,
                Size sFirst, Size* sLast, Size* tLast) const {
            *sLast = cNullItm;
            *tLast = cNullItm;
            //Timer timer;
            const std::vector<Size>*
                vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            Size blsm = graph.mNumVtxs;
            prcsbQue.Push(sFirst);
            sttArr[sFirst] = eSttBfsPrcsb;
            while (prcsbQue.Empty() == false) {
                Size s = prcsbQue.Front();
                prcsbQue.Pop();
                prcsdQue.Push(s);
                sttArr[s] = eSttBfsPrcsd;
                Size sNumEdgs = vtxVecArr[s].size();
                const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                for (Size i = 0; i < sNumEdgs; ++i) {
                    Size t = sVtxArr[i];
                    assert(s != t);
                    // if t is an outer vertex and s and t are in different blossoms then creat a blossom
                    if (sttArr[t] != eSttIdle) {
                        //timer.Start();
                        Size b1 =rprsArr[rFind(linkArr, s)] ;
                        Size b2 = rprsArr[rFind(linkArr, t)];
                        //timer.Stop();
                        if (b1!= b2) {
                            //timer.Start();
                            rProcessBlsm2(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                                    linkArr, rankArr, sttArr, prcsbQue, s, t, blsm);
                            //timer.Stop();
                        }
                        ++blsm;
                        continue;
                    } else if (sttArr[t] == eSttBfsPrcsd) {
                        continue;
                    }
                    Size ss = mateArr[t];

                    if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
                        continue;
                    }
                    if (ss == cNullItm) {
                        // an unmatched vertex is found
                        *sLast = s;
                        *tLast = t;

                        return;
                    }
                    if (sttArr[ss] != eSttIdle) {
                        continue;
                    }
                    ptr1Arr[ss] = s;
                    prcsbQue.Push(ss);
                    sttArr[ss] = eSttBfsPrcsb;
                }
            }
        }
    //////////////////////////         ////////////////////

    template<class ItmQue>
        void MatchingEngine::rFindAugPathCardSglSrc3(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
                Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, ItmQue& prcsdQue,
                Size sFirst, Size* sLast, Size* tLast, Val w2, Size* sedges)
        const {
            *sLast = cNullItm;
            *tLast = cNullItm;
            //Size num_f =0;
            //Size num_s =0;
            //Size nedges=0;
            const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
            Size blsm = graph.mNumVtxs;
            Val maxVtxWght = cNegInfVal;
            prcsbQue.Push(sFirst);
            sttArr[sFirst] = eSttBfsPrcsb;

            while (prcsbQue.Empty() == false) {
                Size s = prcsbQue.Front();
                prcsbQue.Pop();
                prcsdQue.Push(s);
                sttArr[s] = eSttBfsPrcsd;
                Size sNumEdgs = vtxVecArr[s].size();
                const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                for (Size i = 0; i < sNumEdgs; ++i) {
                    Size t = sVtxArr[i];
                    //nedges++;
                    (*sedges)++;
                    assert(s != t);
                    // if t is an outer vertex and s and t are in different blossoms then creat a blossom
                    if (sttArr[t] == eSttBfsPrcsd || sttArr[t] == eSttBfsPrcsb) {
                        if (rprsArr[rFind(linkArr, s)] != rprsArr[rFind(linkArr, t)]) {
                            rProcessBlsm2(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                                    linkArr, rankArr, sttArr, prcsbQue, s, t, blsm);

                            ++blsm;
                        }

                        continue;
                    } else if (sttArr[t] == eSttLast) {
                        continue;
                    }
                    Size ss = mateArr[t];
                    if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
                        continue;
                    }
                    if (ss == cNullItm) {
                        if (maxVtxWght < vtxWghtArr[t]) {
                            // a heavier unmatched vertex is found
                            *sLast = s;
                            *tLast = t;
                            maxVtxWght = vtxWghtArr[t];
                            if(fabs(w2 - maxVtxWght) <= 0.00001)
                                return;
                        }
                        continue;
                    }
                    if (sttArr[ss] != eSttIdle) {
                        continue;
                    }

                    ptr1Arr[ss] = s;
                    prcsbQue.Push(ss);
                    sttArr[ss] = eSttBfsPrcsb;
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rFindIncPathCardSglSrc3(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* ptr2Arr, Size* blsmArr, Size* rprsArr, Size* linkArr,
                Size* rankArr, Stt* sttArr, ItmQue& prcsbQue, ItmQue& prcsdQue,
                Size sFirst, Size* sLast, Size* tLast)
        const {
            *sLast = cNullItm;
            *tLast = cNullItm;
            const std::vector<Size>*
                vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
            Size blsm = graph.mNumVtxs;
            Val minVtxWght= vtxWghtArr[sFirst];
            prcsbQue.Push(sFirst);
            sttArr[sFirst] = eSttBfsPrcsb;

            while (prcsbQue.Empty() == false) {
                Size s = prcsbQue.Front();
                prcsbQue.Pop();
                prcsdQue.Push(s);
                sttArr[s] = eSttBfsPrcsd;
                Size sNumEdgs = vtxVecArr[s].size();
                const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];


                for (Size i = 0; i < sNumEdgs; ++i) {
                    Size t = sVtxArr[i];
                    assert(s != t);
                    // if t is an outer vertex and s and t are in different blossoms then creat a blossom
                    if (sttArr[t] == eSttBfsPrcsd || sttArr[t] == eSttBfsPrcsb) {
                        if (rprsArr[rFind(linkArr, s)] != rprsArr[rFind(linkArr, t)]) {
                            rProcessBlsm6<ItmQue>(graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr,
                                    linkArr, rankArr, sttArr, prcsbQue, s, t, blsm, &minVtxWght,sLast,tLast);
                            ++blsm;
                        }

                        continue;
                    } else if (sttArr[t] == eSttBfsPrcsd) {
                        continue;
                    }
                    Size ss = mateArr[t];
                    if ((ss != cNullItm) && (sttArr[ss] != eSttIdle)) {
                        continue;
                    }
                    if (minVtxWght > vtxWghtArr[ss]) {
                        // a lighter matched vertex is found
                        *sLast = s;
                        *tLast = t;
                        minVtxWght = vtxWghtArr[ss];
                        ptr1Arr[ss] = s;
                        prcsbQue.Push(ss);
                        sttArr[ss] = eSttBfsPrcsb;
                        continue;
                    }

                    if (sttArr[ss] != eSttIdle) {
                        continue;
                    }
                    ptr1Arr[ss] = s;
                    prcsbQue.Push(ss);
                    sttArr[ss] = eSttBfsPrcsb;
                }
            }
        }

    /////////////// Scaling approximation Algorithm functions////
    template<class ItmQue>
        Size MatchingEngine::rFindblsm(Size* p ,Size* blsm, Size* rprs, Size v) const
        {
            Size vblsm=blsm[v];
            if(rprs[vblsm] != cNullItm)
            {
                if(p[vblsm] == cNullItm)
                {
                    return vblsm;
                }
                else
                {
                    // if the blossom tree was updated then find the root blossom by climbing the blossom tree
                    while(p[vblsm] != cNullItm)
                    {
                        vblsm = p[vblsm];
                    }
                    blsm[v]=vblsm;
                    return vblsm;
                }
            }
            else
            {
                // if the blossom has been dissolved then move up the tree until the root is found
                Size newBlsm = v;
                while(p[newBlsm] != cNullItm)
                {
                    newBlsm = p[newBlsm];
                }
                blsm[v]=newBlsm;
                return newBlsm;
            }
        }

    template<class ItmQue>
        void MatchingEngine::augmentBlsm(Size b,Size s,Size* mateArr,Size* blsmArr, Size* rprsArr,
                Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,Size num) const{

            std::stack<Size> tempS;
            tempS.push(b);
            tempS.push(s);

            while(tempS.empty()==false)
            {
                s= tempS.top();
                tempS.pop();
                b = tempS.top();
                tempS.pop();
                Size t = s;
                //find the root blossom
                while(blsmParent[t]!=b)
                {
                    t = blsmParent[t];
                }

                if(t >= num)
                {
                    tempS.push(t);
                    tempS.push(s);
                }
                Size shiftB =b-num;
                // find the postion of t in the order list of blossom children
                Size pos = std::find(blsmChildren[shiftB].begin(), blsmChildren[shiftB].end(), t) - blsmChildren[shiftB].begin();
                Size pos2 =pos*2;

                int inc=0;
                if(pos%2==0)
                    inc=-1;
                else
                    inc=1;
                while(pos != 0)
                {
                    pos+=inc;
                    pos2+=2*inc;
                    t = blsmChildren[shiftB][pos];
                    Size a1 =t;
                    if(t >= num)
                    {
                        if(pos%2==0)
                            a1 = blsmL[shiftB][pos2];
                        else
                            a1 = blsmL[shiftB][pos2-1];
                        tempS.push(t);
                        tempS.push(a1);
                    }
                    pos+=inc;
                    pos2+=2*inc;
                    if(pos >= blsmChildren[shiftB].size())
                    {
                        pos = 0;
                        pos2--;
                    }
                    t = blsmChildren[shiftB][pos];
                    Size a2 =t;
                    if(t >= num)
                    {
                        if(pos%2==0)
                            a2 = blsmL[shiftB][pos2];
                        else
                            a2 = blsmL[shiftB][pos2-1];
                        tempS.push(t);
                        tempS.push(a2);
                    }
                    mateArr[a1] = a2;
                    mateArr[a2] = a1;
                }
                //change the order of the blossom children so it begins from the base and update the base
                std::rotate(blsmChildren[shiftB].begin(),blsmChildren[shiftB].begin()+pos2,blsmChildren[shiftB].end());
                std::rotate(blsmL[shiftB].begin(),blsmL[shiftB].begin()+(pos2*2),blsmL[shiftB].end());
                rprsArr[b] = s;
            }
        }


    template<class ItmIdxdQue,class ItmQue>
        void MatchingEngine::rAugment(Size* mateArr,Size* ptrArr,
                ItmIdxdQue& expsdQue, Size sLast, Size tLast, Size* blsmArr,
                Size* rprsArr, Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,Size num, Size* label, Size* blsmin) const {

            Size s = sLast;
            Size t = tLast;
            expsdQue.Erase(t);
            while (s != cNullItm) {
                Size blsmt=t;
                if(blsmParent[t]!=cNullItm)
                    blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                label[blsmt]=0;
                if(blsmt >= num && blsmin[blsmt-num]!= cNullItm)
                {
                    Size tempT =blsmin[blsmt-num];
                    //augment the blossom
                    augmentBlsm<ItmQue>(blsmArr[t],tempT,mateArr,blsmArr, rprsArr, blsmParent, blsmChildren,blsmL,num);
                    t= tempT;
                }

                Size blsms=s;
                if(blsmParent[s]!=cNullItm)
                    blsms =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                if(blsms >= num && rprsArr[blsms]!=s)
                {
                    //if s is in a blossom
                    label[blsms]=0;
                    Size base = rprsArr[blsms];
                    Size mateBase = mateArr[base];
                    //augment the blossom
                    augmentBlsm<ItmQue>(blsms,s,mateArr,blsmArr, rprsArr, blsmParent, blsmChildren,blsmL,num);

                    label[s]=0;
                    label[t]=0;
                    mateArr[s] = t;
                    mateArr[t] = s;
                    s=base;
                    if(ptrArr[s] == cNullItm || mateBase ==cNullItm )
                    {
                        expsdQue.Erase(s);
                    }
                    s = ptrArr[s];
                    t=mateBase;
                }
                else
                {
                    //if s is not in blossom
                    Size tt = mateArr[s];
                    label[s]=0;
                    label[t]=0;
                    mateArr[s] = t;
                    mateArr[t] = s;
                    if(ptrArr[s] == cNullItm)
                    {
                        expsdQue.Erase(s);
                    }
                    s = ptrArr[s];
                    t=tt;
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rDissolveBlsm3(const Graph& graph, Size* blsmArr,Size* blsmParent , std::vector<Size>* blsmChildren, Size b) const {

            for(Size i=0; i < blsmChildren[b-graph.mNumVtxs].size();i++)
            {
                Size child = blsmChildren[b-graph.mNumVtxs][i];
                blsmParent[child] = cNullItm;
                if( child < graph.mNumVtxs)
                {
                    blsmArr[child]= child;
                }
            }
        }

    template <class ItmQue,class ItmStk>
        void MatchingEngine::rProcessBlsm3(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* blsmArr, Size* rprsArr, Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL, Size* label,
                Stt* sttArr, ItmStk& prcsbStk, Size s1, Size s2, Size* blsm, Val* dualArr, Size* blsmin, Size* treeNumArr) const {

            std::queue<Size> tempQ;
            Size numVtxs = graph.mNumVtxs;
            Size bs1=s1;
            Size bs2=s2;
            if(blsmParent[s1]!=cNullItm)
                bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
            if(blsmParent[s2]!=cNullItm)
                bs2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s2);
            Size base =rprsArr[bs2];
            // Create blossom.
            Size b = *blsm;
            Size shiftB =b-numVtxs;
            rprsArr[b] = base;
            treeNumArr[b]=treeNumArr[bs2];
            blsmParent[b] = cNullItm;
            blsmParent[bs2] = b;
            blsmChildren[shiftB].push_back(bs2);
            blsmL[shiftB].push_back(s2);
            blsmL[shiftB].push_back(s1);
            // iterate over the vertices in backward order until the base is reached
            while (bs1 != b)
            {
                Size tempS1 = s1;
                blsmParent[bs1] = b;
                blsmChildren[shiftB].push_back(bs1);
                Size mates1;
                if(bs1 < numVtxs)
                {
                    blsmL[shiftB].push_back(bs1);
                    mates1 = mateArr[bs1];
                }
                else
                {
                    Size rprsbs1 = rprsArr[bs1];
                    if(tempS1 != rprsbs1)
                    {
                        ptr1Arr[tempS1]= ptr1Arr[base];
                        s1= rprsbs1;
                        tempS1 = rprsbs1;
                    }

                    blsmL[shiftB].push_back(rprsbs1);
                    mates1 = mateArr[rprsbs1];
                }
                blsmL[shiftB].push_back(mates1);
                bs1= mates1;
                if(blsmParent[mates1]!=cNullItm)
                    bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, mates1);
                blsmParent[bs1] = b;
                blsmChildren[shiftB].push_back(bs1);

                if(bs1 < numVtxs)
                {
                    blsmL[shiftB].push_back(mates1);
                    if(sttArr[bs1] == eSttIdle || (label[bs1]== 2 || label[bs1]== 0))
                    {
                        prcsbStk.Push(bs1);
                        sttArr[bs1] = eSttDfsPrcsb;
                        ptr1Arr[bs1]= ptr1Arr[base];
                    }
                }
                else
                {
                    if(blsmin[bs1-numVtxs]==cNullItm || blsmin[bs1-numVtxs]==rprsArr[bs1])
                    {
                        blsmL[shiftB].push_back(mates1);
                        ptr1Arr[mates1]= ptr1Arr[base];
                    }
                    else
                    {
                        ptr1Arr[blsmin[bs1-numVtxs]]= ptr1Arr[base];
                        blsmL[shiftB].push_back(blsmin[bs1-numVtxs]);
                    }

                    if(sttArr[rprsArr[bs1]] == eSttIdle)
                    {
                        tempQ.push(bs1-numVtxs);
                        while (tempQ.empty() == false)
                        {
                            Size tempB = tempQ.front();
                            tempQ.pop();
                            for(Size i=0; i < blsmChildren[tempB].size();i++)
                            {
                                Size child = blsmChildren[tempB][i];
                                if( child < numVtxs)
                                {
                                    if(sttArr[child] == eSttIdle)
                                    {
                                        prcsbStk.Push(child);
                                        sttArr[child] = eSttDfsPrcsb;
                                        ptr1Arr[child] = ptr1Arr[base];
                                    }
                                }
                                else
                                {
                                    tempQ.push(child-numVtxs);
                                }
                            }
                        }
                    }
                    label[bs1] = 1;
                }
                s1 = ptr1Arr[s1];
                ptr1Arr[tempS1]= ptr1Arr[base];
                blsmL[shiftB].push_back(s1);
                bs1=s1;
                if(blsmParent[s1]!=cNullItm)
                    bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
            }
            label[b] = 1;
            dualArr[b] = 0;
            (*blsm)++;
        }


    template<class ItmQue, class ItmStk, class ItmIdxdQue>
        void MatchingEngine::rFindMxmlSetAugPathsCard(const Graph& graph,
                Size* mateArr, Size* ptrArr, Stt* sttArr,
                Size* idxArr, ItmStk& prcsbStk, ItmQue& prcsdQue, ItmIdxdQue& expsdQue,
                ItmQue& sLastQue, ItmQue& tLastQue, Size* blsmArr,Size* rprsArr, Size* blsmParent,
                std::vector<Size>* blsmChildren,std::vector<Size>* blsmL, Size* label, Val* dualArr,
                Size* blsm, Val delta, Size* blsmin, Size* treeNumArr,Size iter, Val ep,  Size*se) const {

            Size numVtxs = graph.mNumVtxs;
            Val gamma = logbase2(1.0/ep);
            const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
            std::queue<Size> tempQ;
            Size sFirst = expsdQue.First();
            while (sFirst != cNullItm)
            {
                if(sttArr[sFirst] != eSttIdle)
                {
                    sFirst = expsdQue.Next(sFirst);
                    continue;
                }
                prcsbStk.Push(sFirst);
                sttArr[sFirst] = eSttDfsPrcsb;
                label[sFirst]=1;
                treeNumArr[sFirst]=sFirst;

                Size blsmf= sFirst;
                if(blsmParent[sFirst]!=cNullItm)
                    blsmf = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, sFirst);
                // if sFirst is in a blosso then add all vertices within the blossom to the search stack
                if(blsmf  >= numVtxs)
                {
                    label[blsmf]=1;
                    treeNumArr[blsmf]=sFirst;
                    tempQ.push(blsmf-numVtxs);
                    while (tempQ.empty() == false)
                    {
                        Size tempB = tempQ.front();
                        tempQ.pop();
                        for(Size i=0; i < blsmChildren[tempB].size();i++)
                        {
                            Size child = blsmChildren[tempB][i];
                            if( child < numVtxs)
                            {
                                if(sttArr[child] == eSttIdle)
                                {
                                    prcsbStk.Push(child);
                                    sttArr[child] = eSttDfsPrcsb;
                                }
                            }
                            else
                            {
                                tempQ.push(child-numVtxs);
                            }
                        }
                    }
                }
                while (prcsbStk.Empty() == false)
                {

                    Size s = prcsbStk.Top();
                    Size sNumEdgs = vtxVecArr[s].size();
                    const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                    Size cases = 0;
                    while (idxArr[s] < sNumEdgs)
                    {
                        Size t = vtxArr[idxArr[s]];
                        (*se)++;
                        Val wght = vtxWghtArr[s]+vtxWghtArr[t];
                        Size ss;
                        Size base;
                        /////////
                        if(t == mateArr[s])
                        {
                            ++(idxArr[s]);
                            continue;
                        }

                        Size blsmt=t;
                        Size blsms =s;
                        if(blsmParent[t]!=cNullItm)
                            blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                        if(blsmParent[s]!=cNullItm)
                            blsms = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);

                        if(label[blsmt]==1 &&  treeNumArr[blsms] == treeNumArr[blsmt])
                        {
                            if(fabs((dualArr[s]+dualArr[t] )- ((int)(wght/delta)*delta - delta)) < 0.000000119 && blsms != blsmt && sttArr[t] == eSttDfsPrcsd && logbase2(wght)>=iter-gamma)
                            {
                                //create a blossom
                                rProcessBlsm3<ItmQue, ItmStk>(graph, mateArr,ptrArr,blsmArr, rprsArr, blsmParent,blsmChildren,blsmL,label,sttArr, prcsbStk,t, s, blsm, dualArr,blsmin,treeNumArr);
                                break;
                            }
                            ++(idxArr[s]);
                            continue;
                        }

                        if(blsmt >=numVtxs)
                        {
                            base = rprsArr[blsmt];
                            ss = mateArr[base];
                            cases=1;
                        }
                        else
                        {
                            ss = mateArr[t];
                            if(ss !=cNullItm)
                                cases= 2;
                            else
                                cases= 3;
                        }
                        if (((ss == cNullItm) && (sttArr[t] != eSttIdle)) || label[blsmt]==2 || (treeNumArr[blsmt] !=cNullItm && treeNumArr[blsms] != treeNumArr[blsmt] )) {

                            ++(idxArr[s]);
                        }
                        else
                        {
                            // check if the edge (s,t) is an eligible edge
                            if(fabs((dualArr[s]+dualArr[t] )- ((int)(wght/delta)*delta - delta)) <= 0.000000119 && logbase2(wght)>=iter-gamma)
                            {

                                if(ss == cNullItm)
                                    break;
                                Size tempv = t;
                                if(cases == 1 )
                                {
                                    tempv = base;
                                    if(ss == cNullItm)
                                    {
                                        break;
                                    }
                                }
                                double result = dualArr[ss]+dualArr[tempv] - (floor((vtxWghtArr[ss]+vtxWghtArr[tempv])/delta)*delta );
                                bool eligible = false;
                                if(result >= 0)
                                {    // check if (t,matet) is an eligble edge
                                    if(fabs(result - floorf(result)) <= 0.000000119  &&  fabs(((round(result/delta)*delta)- result)) <= (0.000000119) && logbase2(vtxWghtArr[ss]+vtxWghtArr[tempv])>=iter-gamma)
                                    {
                                        eligible = true;
                                    }
                                }
                                if(eligible)
                                {
                                    break;
                                }
                                else
                                {
                                    label[t]=2;
                                    ++(idxArr[s]);
                                }
                            }
                            else
                                ++(idxArr[s]);
                        }
                    }
                    if(idxArr[s] < sNumEdgs)
                    {
                        Size t = vtxArr[idxArr[s]];
                        Size blsmt=t;
                        Size blsms =s;
                        if(blsmParent[t]!=cNullItm)
                            blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                        if(blsmParent[s]!=cNullItm)
                            blsms = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                        if(blsms == blsmt)
                        {
                            ++(idxArr[s]);
                            continue;
                        }
                        Size ss = mateArr[t];

                        if(ss != cNullItm)
                        {
                            Size blsmss=ss;
                            if(blsmParent[ss]!=cNullItm)
                                blsmss = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, ss);
                            if(label[blsmt]==0)
                            {
                                //if t has not been considered
                                if(blsmt >= numVtxs)
                                {
                                    // if t is on a blossom
                                    Size base = rprsArr[blsmt];
                                    blsmin[blsmt-numVtxs]=t;
                                    label[blsmt]=2;
                                    treeNumArr[blsmt]=treeNumArr[blsms];
                                    ss = mateArr[base];
                                    if(ss == cNullItm)
                                    {
                                        // if the base of a blossom is unmatched then an augmenting path is found
                                        if(sttArr[base] == eSttIdle)
                                        {
                                            label[blsmt]=1;
                                            sLastQue.Push(s);
                                            tLastQue.Push(base);
                                            sttArr[base] = eSttLast;
                                            while (prcsbStk.Empty() == false)
                                            {
                                                Size sss = prcsbStk.Top();
                                                prcsbStk.Pop();
                                                prcsdQue.Push(sss);
                                                sttArr[sss] = eSttDfsPrcsd;
                                            }
                                            break;
                                        }
                                        else
                                        {
                                            ++(idxArr[s]);
                                            continue;
                                        }
                                    }
                                    // other wise move to the base of a blossom
                                    ptrArr[ss] = s;
                                    blsmss=ss;
                                    if(blsmParent[ss]!=cNullItm)
                                        blsmss = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, ss);
                                    treeNumArr[blsmss]=treeNumArr[blsms];
                                    if(blsmss >= numVtxs)
                                    {
                                        // if mste of t is in a blosso then add all vertices within the blossom to the search stack
                                        prcsbStk.Push(ss);
                                        sttArr[ss] = eSttDfsPrcsb;
                                        label[blsmss]=1;
                                        tempQ.push(blsmss-numVtxs);
                                        while (tempQ.empty() == false)
                                        {
                                            Size tempB = tempQ.front();
                                            tempQ.pop();
                                            for(Size i=0; i < blsmChildren[tempB].size();i++)
                                            {
                                                Size child = blsmChildren[tempB][i];
                                                if( child < numVtxs)
                                                {
                                                    if(sttArr[child] == eSttIdle)
                                                    {
                                                        prcsbStk.Push(child);
                                                        sttArr[child] = eSttDfsPrcsb;
                                                        ptrArr[child] = s;
                                                    }
                                                }
                                                else
                                                {
                                                    tempQ.push(child-numVtxs);
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        label[ss]=1;
                                        prcsbStk.Push(ss);
                                        sttArr[ss] = eSttDfsPrcsb;
                                    }
                                    ++(idxArr[s]);
                                }
                                else if(blsmss >= numVtxs)
                                {
                                    // if mate of t is in a blosso then add all vertices within the blossom to the search stack
                                    label[blsmt]=2;
                                    treeNumArr[blsmt]=treeNumArr[blsms];
                                    ptrArr[ss] = s;
                                    prcsbStk.Push(ss);
                                    sttArr[ss] = eSttDfsPrcsb;
                                    label[blsmss]=1;
                                    treeNumArr[blsmss]=treeNumArr[blsms];
                                    tempQ.push(blsmss-numVtxs);
                                    while (tempQ.empty() == false)
                                    {
                                        Size tempB = tempQ.front();
                                        tempQ.pop();
                                        for(Size i=0; i < blsmChildren[tempB].size();i++)
                                        {
                                            Size child = blsmChildren[tempB][i];
                                            if( child < numVtxs)
                                            {
                                                if(sttArr[child] == eSttIdle)
                                                {
                                                    ptrArr[child] = s;
                                                    prcsbStk.Push(child);
                                                    sttArr[child] = eSttDfsPrcsb;
                                                }
                                            }
                                            else
                                            {
                                                tempQ.push(child-numVtxs);
                                            }
                                        }
                                    }
                                    ++(idxArr[s]);
                                }
                                else
                                {
                                    // add the mate of t to the search stack
                                    label[blsmt]=2;
                                    treeNumArr[blsmt]=treeNumArr[blsms];
                                    label[ss]=1;
                                    treeNumArr[ss]=treeNumArr[blsms];
                                    ptrArr[ss] = s;
                                    prcsbStk.Push(ss);
                                    sttArr[ss] = eSttDfsPrcsb;
                                    ++(idxArr[s]);
                                }
                            }
                            else
                                ++(idxArr[s]);
                        }
                        else
                        {
                            // an ugmenting path is found
                            sLastQue.Push(s);
                            tLastQue.Push(t);
                            sttArr[t] = eSttLast;
                            label[blsmt]=1;
                            if(blsmt>=numVtxs)
                                blsmin[blsmt-numVtxs]=cNullItm;
                            while (prcsbStk.Empty() == false)
                            {
                                Size sss = prcsbStk.Top();
                                prcsbStk.Pop();
                                prcsdQue.Push(sss);
                                sttArr[sss] = eSttDfsPrcsd;
                            }
                            break;
                        }
                    }
                    else
                    {
                        prcsbStk.Pop();
                        prcsdQue.Push(s);
                        sttArr[s] = eSttDfsPrcsd;
                    }
                }
                sFirst = expsdQue.Next(sFirst);
            }
        }

    template<class ItmQue, class ItmStk, class ItmIdxdQue>
        void MatchingEngine::rScaleOneMinusEpsVtxWght(const Graph& graph, Size* mateArr, Size* card, Val* vtxWght, Val ep) const {

            Size numVtxs = graph.mNumVtxs;
            if (mInlz == false) {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[numVtxs], cNullItm);
                }
                *card = 0;
            } else {
                Size iniCard(0);
                rInlzGrdyForCard(graph, mateArr, &iniCard);

                *card = iniCard;
            }
            //to store a pointer from an outer vertex to its parent in the search tree
            std::vector<Size> ptrVec;
            ResizeVector<Size>(&ptrVec,numVtxs);
            Size* ptrArr = (numVtxs == 0) ? NULL : &ptrVec[0];
            if (ptrArr != NULL) {
                std::fill(&ptrArr[0], &ptrArr[numVtxs], cNullItm);
            }
            //to store a blossom id for a vertex
            std::vector<Size> blsmVec;
            ResizeVector<Size>(&blsmVec, numVtxs);
            Size* blsmArr = (numVtxs == 0) ? NULL : &blsmVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                blsmArr[v] = v;
            }
            //to store a base vertex for a vertex or a blossom
            std::vector<Size> rprsVec;
            ResizeVector<Size>(&rprsVec, numVtxs*2);
            Size* rprsArr = (numVtxs == 0) ? NULL : &rprsVec[0];
            for (Size v = 0; v < numVtxs; ++v) {
                rprsArr[v] = v;
            }
            std::fill(&rprsArr[numVtxs], &rprsArr[numVtxs*2], cNullItm);
            for (Size v = numVtxs; v < numVtxs*2; ++v) {
                rprsArr[v] = cNullItm;
            }
            //2
            //to store a label for a vertex (0 for unlabelled vertex 1 for inner vertex 2 for an outer vertex)
            std::vector<Size> labelVec;
            ResizeVector<Size>(&labelVec, numVtxs*2);
            Size* label = (numVtxs == 0) ? NULL : &labelVec[0];
            if (label != NULL) {
                std::fill(&label[0], &label[numVtxs], 0);
            }

            //to store a blossom parent in a blossom tree
            std::vector<Size> blsmParentVec;
            ResizeVector<Size>(&blsmParentVec, numVtxs*2);
            Size* blsmParentArr = (numVtxs == 0) ? NULL : &blsmParentVec[0];
            if (blsmParentArr != NULL) {
                std::fill(&blsmParentArr[0], &blsmParentArr[numVtxs*2], cNullItm);

            }

            //to store a blossom children in a blossom tree
            std::vector<std::vector<Size> > blsmChildrenVec;
            ResizeVector<std::vector<Size> >(&blsmChildrenVec, numVtxs);
            std::vector<Size>* blsmChildrenArr = (numVtxs == 0) ? NULL : &blsmChildrenVec[0];

            //to store a blossom edges
            std::vector<std::vector<Size> > blsmLVec;
            ResizeVector<std::vector<Size> >(&blsmLVec, numVtxs);
            std::vector<Size>* blsmL = (numVtxs == 0) ? NULL : &blsmLVec[0];

            //to store which vertex an augmenting path is using in a blossom
            std::vector<Size> blsminVec;
            ResizeVector<Size>(&blsminVec, numVtxs);
            Size* blsmin = (numVtxs == 0) ? NULL : &blsminVec[0];
            if (blsmin != NULL) {
                std::fill(&blsmin[0], &blsmin[numVtxs], cNullItm);
            }
            //2
            //to store a tree id for each vertex or a blossom
            std::vector<Size> treeNumVec;
            ResizeVector<Size>(&treeNumVec, numVtxs*2);
            Size* treeNumArr = (numVtxs == 0) ? NULL : &treeNumVec[0];
            if (treeNumArr != NULL) {
                std::fill(&treeNumArr[0], &treeNumArr[numVtxs*2], cNullItm);
            }

            const std::vector<Size>*
                vtxVecArr = (numVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            //const std::vector<Val>* edgWghtVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
            const Val* vtxWghtArr = (numVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

            Val maxweight = 0.0;
            // find the max edge weight
            for (Size v = 0; v < numVtxs; ++v) {
                Size sNumEdgs = vtxVecArr[v].size();
                const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
                for (Size i = 0; i < sNumEdgs; ++i) {
                    Size t = vtxArr[i];

                    if(vtxWghtArr[v]+vtxWghtArr[t] > maxweight)
                        maxweight = vtxWghtArr[v]+vtxWghtArr[t];
                }
            }


            Val delta = maxweight*ep;
            //2
            //to store dual values for each vertex of a blossom
            std::vector<Val> dualVec;
            ResizeVector<Val>(&dualVec, numVtxs*2);
            Val* dualArr = (numVtxs== 0) ? NULL : &dualVec[0];
            std::fill(&dualArr[0], &dualArr[numVtxs], maxweight/2.0 - delta/2.0);
            std::fill(&dualArr[numVtxs], &dualArr[numVtxs*2], 0);

            //to store a status for each vertex or a blossom
            std::vector<Stt> sttVec;
            ResizeVector<Stt>(&sttVec, numVtxs);
            Stt* sttArr = (numVtxs == 0) ? NULL : &sttVec[0];
            if (sttArr != NULL) {
                std::fill(&sttArr[0], &sttArr[numVtxs], eSttIdle);
            }
            std::vector<Size> idxVec;
            ResizeVector<Size>(&idxVec, numVtxs);
            Size* idxArr = (numVtxs == 0) ? NULL : &idxVec[0];
            if (idxArr != NULL) {
                std::fill(&idxArr[0], &idxArr[numVtxs], 0);
            }


            //Timer timer;
            ItmQue sLastQue(numVtxs);
            ItmQue tLastQue(numVtxs);

            ItmStk prcsbStk(numVtxs);
            ItmQue prcsdQue(numVtxs);
            ItmIdxdQue expsdQue(numVtxs);

            Size blsm = numVtxs;

            for (Size s = 0; s < numVtxs; ++s)
            {
                expsdQue.Push(s);
            }
            Size se=0;
            for(Size i = 0; i <= logbase2(maxweight) ; i++)
            {
                bool doneScale=true;

                do
                {
                    doneScale=true;

                    //maximal set of vertex disjoint augmenting paths
                    rFindMxmlSetAugPathsCard<ItmQue, ItmStk, ItmIdxdQue>
                        (graph, mateArr, ptrArr, sttArr, idxArr,
                         prcsbStk, prcsdQue, expsdQue, sLastQue, tLastQue,blsmArr,rprsArr,blsmParentArr,blsmChildrenArr,blsmL,label,
                         dualArr, &blsm, delta,blsmin,treeNumArr,i,ep,&se);

                    // augmenting the set of maximal  augmenting paths
                    while (tLastQue.Empty() == false)
                    {
                        Size sLast = sLastQue.Front();
                        sLastQue.Pop();
                        Size tLast = tLastQue.Front();
                        tLastQue.Pop();
                        sttArr[tLast] = eSttIdle;
                        rAugment<ItmIdxdQue,ItmQue>(mateArr, ptrArr, expsdQue, sLast, tLast, blsmArr, rprsArr, blsmParentArr,
                                blsmChildrenArr,blsmL,numVtxs,label,blsmin);

                        ++(*card);
                    }
                    // vertices dual updates
                    for (Size v = 0; v < numVtxs; ++v)
                    {
                        Size blsmv=v;
                        if(blsmParentArr[v]!=cNullItm)
                        {
                            blsmv = rFindblsm<ItmQue>(blsmParentArr,blsmArr, rprsArr, v);
                        }
                        if( blsmv>=graph.mNumVtxs && label[blsmv]==1)
                        {
                            dualArr[v] = fabs(dualArr[v] - delta/2);
                        }
                        else if(blsmv>=graph.mNumVtxs && label[ blsmv]==2)
                        {
                            dualArr[v] += delta/2;
                        }
                        else if(label[v]==1 || mateArr[v]==cNullItm)
                        {
                            dualArr[v] = fabs(dualArr[v] - delta/2);
                        }
                        else if(label[v]==2)
                        {
                            dualArr[v] += delta/2;
                        }
                        label[v]=0;
                        treeNumArr[v]=cNullItm;

                    }

                    //blossoms dual updates
                    for(Size v = numVtxs; v < blsm; ++v)
                    {
                        if(rprsArr[v]!=cNullItm && blsmParentArr[v]==cNullItm)
                        {
                            if(label[v]==1)
                            {
                                dualArr[v] += delta;
                            }
                            else if(label[v]==2)
                            {
                                dualArr[v] = fabs(dualArr[v]-delta);
                            }
                        }
                        label[v]=0;
                        treeNumArr[v]=cNullItm;
                    }

                    //dissolving inner blossoms if the dual becomes zero
                    for (Size v = numVtxs; v < blsm; ++v)
                    {
                        if(rprsArr[v]!=cNullItm && blsmParentArr[v]==cNullItm && fabs(dualArr[v]) <= 0.000000119)
                        {
                            rprsArr[v]=cNullItm;
                            rDissolveBlsm3<ItmQue>(graph, blsmArr,blsmParentArr, blsmChildrenArr, v);
                        }
                    }

                    //check if we reach the limit of iterations within a scale
                    if(i <logbase2(maxweight) )
                    {
                        Size v = expsdQue.First();
                        while (v != cNullItm)
                        {

                            if(dualArr[v] > (maxweight/(pow(2,i+2)))-delta/2)
                            {
                                doneScale=false;
                                break;
                            }
                            v = expsdQue.Next(v);
                        }
                    }
                    else
                    {
                        Size v = expsdQue.First();
                        while (v != cNullItm)
                        {

                            if(dualArr[v] > 0.000000119)
                            {
                                doneScale=false;
                                break;
                            }
                            v = expsdQue.Next(v);
                        }
                    }
                    while (prcsdQue.Empty() == false)
                    {
                        Size s = prcsdQue.Front();
                        prcsdQue.Pop();
                        ptrArr[s] = cNullItm;
                        sttArr[s] = eSttIdle;
                        idxArr[s] = 0;
                    }
                }while(!doneScale);

                // update the duals
                delta = delta/2.0;
                for (Size v = 0; v < numVtxs; ++v)
                {
                    dualArr[v] += delta;
                }
            }// main scale loop

            std::cout << se<<" ";
            *vtxWght = 0.0;
            for (Size s = 0; s < numVtxs; ++s) {
                if (mateArr[s] != cNullItm) {
                    *vtxWght += vtxWghtArr[s];

                }
            }
        }

    //////////////////////////Maximum cardinality matching Gabow's Algorithm/////////////////

    /*template<class ItmQue>
      Size MatchingEngine::rFindblsmGabow(Size* blsm, Size* rprs, Size v) const
      {

      if(rprs[v] == v)
      {

      return blsm[v];
      }
      else
      {
    //std::cout << "strat find"<<std::endl;
    Size tempR = rprs[v];
    while(rprs[tempR] != tempR)
    {
    tempR = rprs[tempR];
    }
    rprs[v] = tempR;
    blsm[v]=blsm[tempR];
    //std::cout << "finish find"<<std::endl;
    return blsm[tempR];
    }
    }*/

    template<class ItmQue>
        void MatchingEngine::rDissolveBlsm3Gabow(const Graph& graph, Size* blsmArr,Size* blsmParent , std::vector<Size>* blsmChildren, Size b,Size blsmbdual) const {
            for(Size i=0; i < blsmChildren[b-graph.mNumVtxs].size();i++)
            {
                Size child = blsmChildren[b-graph.mNumVtxs][i];
                blsmParent[child] = cNullItm;
                if( child < graph.mNumVtxs)
                {
                    blsmArr[child]= child;
                }
            }
        }

    template<class ItmIdxdQue,class ItmQue>
        void MatchingEngine::rAugmentGabow(Size* mateArr,Size* ptrArr,
                ItmIdxdQue& expsdQue, Size sLast, Size tLast, Size* blsmArr, Size* rprsArr, Size* blsmParent,
                std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,Size num, Size* blsmin) const {
            Size s = sLast;
            Size t = tLast;
            expsdQue.Erase(t);
            while (s != cNullItm) {
                Size blsmt=t;
                if(blsmParent[t]!=cNullItm)
                    blsmt= rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                if(blsmt >= num && blsmin[blsmt-num]!= cNullItm)
                {
                    Size tempT =blsmin[blsmt-num];
                    augmentBlsm<ItmQue>(blsmArr[t],tempT,mateArr,blsmArr, rprsArr, blsmParent, blsmChildren,blsmL,num);
                    t= tempT;
                }
                Size blsms=s;
                if(blsmParent[s]!=cNullItm)
                    blsms= rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                if(blsms >= num && rprsArr[blsms]!=s)
                {
                    Size base = rprsArr[blsms];
                    Size mateBase = mateArr[base];
                    augmentBlsm<ItmQue>(blsms,s,mateArr,blsmArr, rprsArr, blsmParent, blsmChildren,blsmL,num);
                    mateArr[s] = t;
                    mateArr[t] = s;
                    s=base;
                    if(ptrArr[s] == cNullItm || mateBase ==cNullItm )
                    {
                        expsdQue.Erase(s);
                    }
                    s = ptrArr[s];
                    t=mateBase;
                }
                else
                {
                    Size tt = mateArr[s];
                    mateArr[s] = t;
                    mateArr[t] = s;

                    if(ptrArr[s] == cNullItm)
                    {
                        expsdQue.Erase(s);
                    }
                    s = ptrArr[s];
                    t=tt;
                }
            }
        }


    template <class ItmQue,class ItmStk>
        void MatchingEngine::rProcessBlsm3Gabow(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* blsmArr, Size* rprsArr, Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,
                Stt* sttArr, ItmStk& prcsbStk, Size s1, Size s2, Size* blsm, Size* blsmin, Size* treeNumArr) const {

            Size numVtxs = graph.mNumVtxs;
            //Timer timer;
            Size bs2=s2;
            Size bs1=s1;
            if(blsmParent[s2]!=cNullItm)
                bs2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s2);
            if(blsmParent[s1]!=cNullItm)
                bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
            Size base =rprsArr[bs2];
            // Create blossom.

            Size b = *blsm;
            Size shiftB =b-numVtxs;
            rprsArr[b] = base;
            blsmParent[b] = cNullItm;
            blsmParent[bs2] = b;
            blsmChildren[shiftB].clear();
            blsmL[shiftB].clear();
            blsmChildren[shiftB].push_back(bs2);
            blsmL[shiftB].push_back(s2);
            blsmL[shiftB].push_back(s1);

            Size tn = treeNumArr[s2];
            while (bs1 != b)
            {
                Size tempS1 = s1;
                blsmParent[bs1] = b;
                blsmChildren[shiftB].push_back(bs1);
                Size mates1;
                if(bs1 < numVtxs)
                {
                    blsmL[shiftB].push_back(bs1);
                    blsmArr[bs1]=b;
                    mates1 = mateArr[bs1];
                }
                else
                {
                    Size rprsbs1 = rprsArr[bs1];
                    if(tempS1 != rprsbs1)
                    {
                        ptr1Arr[tempS1]= ptr1Arr[base];
                        s1= rprsbs1;
                        tempS1 = rprsbs1;

                    }
                    blsmL[shiftB].push_back(rprsbs1);
                    blsmArr[rprsbs1]=b;
                    mates1 = mateArr[rprsbs1];
                }
                blsmL[shiftB].push_back(mates1);
                bs1=mates1;
                if(blsmParent[mates1]!=cNullItm)
                    bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, mates1);
                blsmParent[bs1] = b;
                blsmChildren[shiftB].push_back(bs1);
                if(bs1 < numVtxs)
                {
                    blsmL[shiftB].push_back(bs1);
                    ptr1Arr[bs1]= ptr1Arr[base];
                    blsmArr[bs1]=b;
                    if(sttArr[bs1] == eSttIdle)
                    {
                        treeNumArr[bs1] =tn;
                        sttArr[bs1] = eSttDfsPrcsb;
                    }
                }
                else
                {
                    if(blsmin[bs1-numVtxs]==cNullItm || blsmin[bs1-numVtxs]==rprsArr[bs1])
                    {
                        blsmL[shiftB].push_back(mates1);
                        ptr1Arr[mates1]= ptr1Arr[base];
                    }

                    else
                    {
                        blsmL[shiftB].push_back(blsmin[bs1-numVtxs]);
                        ptr1Arr[blsmin[bs1-numVtxs]]= ptr1Arr[base];
                    }
                    std::queue<Size> tempQ;
                    tempQ.push(bs1-numVtxs);
                    while (tempQ.empty() == false)
                    {
                        Size tempB = tempQ.front();
                        tempQ.pop();
                        for(Size i=0; i < blsmChildren[tempB].size();i++)
                        {
                            Size child = blsmChildren[tempB][i];
                            if( child < numVtxs)
                            {
                                if(sttArr[child] == eSttIdle)
                                {
                                    prcsbStk.Push(child);
                                    sttArr[child] = eSttDfsPrcsb;
                                    ptr1Arr[child] = ptr1Arr[base];
                                    treeNumArr[child] =tn;
                                }
                                blsmArr[child]=b;
                            }
                            else
                            {
                                tempQ.push(child-numVtxs);
                            }
                        }
                    }
                }
                s1 = ptr1Arr[s1];
                ptr1Arr[tempS1]= ptr1Arr[base];
                blsmL[shiftB].push_back(s1);
                bs1=s1;
                if(blsmParent[s1]!=cNullItm)
                    bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
            }
            blsmArr[s1]=b;
            blsmArr[s2]=b;
            (*blsm)++;
        }

    template <class ItmQue,class ItmStk>
        void MatchingEngine::rProcessBlsm5(const Graph& graph, Size* mateArr,
                Size* ptr1Arr, Size* blsmArr, Size* rprsArr, Size* blsmParent, std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,
                Stt* sttArr, ItmQue& prcsbStk, Size s1, Size s2, Size* blsm, Size* trace, Size* treeNumArr,
                int* tOuter, int* Delta) const {

            Size numVtxs = graph.mNumVtxs;
            //Timer timer;

            Size bs1=s1;
            Size bs2=s2;
            if(blsmParent[s2]!=cNullItm)
                bs2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s2);
            if(blsmParent[s1]!=cNullItm)
                bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
            Size base;
            Size baseBlsm;

            Size tmpRprsArr[2];
            tmpRprsArr[0] = rprsArr[bs1];
            tmpRprsArr[1] = rprsArr[bs2];
            trace[tmpRprsArr[0]] = *blsm;
            trace[tmpRprsArr[1]] = *blsm;
            bool done = false;
            do {
                for (Size k = 0; k < 2; ++k) {
                    Size r = tmpRprsArr[k];
                    Size s = ptr1Arr[r];
                    if (s != cNullItm) {
                        baseBlsm=s;
                        if(blsmParent[s]!=cNullItm)
                            baseBlsm = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                        r = rprsArr[baseBlsm];
                        if (trace[r] == *blsm) {
                            base = r;
                            done = true;
                            break;
                        }
                        trace[r] = *blsm;
                        tmpRprsArr[k] = r;
                    }
                }
            } while (done == false);
            Size b = *blsm;
            rprsArr[b] = base;
            blsmParent[b] = cNullItm;

            Size shiftB =b-numVtxs;
            blsmChildren[shiftB].clear();
            blsmL[shiftB].clear();

            Size tempB,tempM;;
            blsmL[shiftB].push_back(s2);
            blsmL[shiftB].push_back(s1);

            //Size pathlen = length[s2]+1;
            //Size templen = length[s1]+1;
            while (baseBlsm != bs1)
            {
                tempB = rprsArr[bs1];
                tempM = mateArr[tempB];
                blsmL[shiftB].push_back(tempB);
                blsmL[shiftB].push_back(tempM);
                Size tempMblsm=tempM;
                if(blsmParent[tempM]!=cNullItm)
                    tempMblsm = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, tempM);
                blsmParent[bs1] = b;
                blsmChildren[shiftB].push_back(bs1);
                blsmParent[tempMblsm] = b;
                blsmChildren[shiftB].push_back(tempMblsm);
                blsmArr[tempB]=b;
                blsmArr[tempM]=b;
                if(tempMblsm < numVtxs)
                {
                    if(sttArr[tempMblsm] == eSttIdle)
                    {
                        //length[tempMblsm] = length[s1]-length[tempB]+pathlen+1;
                        //pathlen = length[tempMblsm] +1;
                        prcsbStk.Push(tempMblsm);
                        tOuter[tempMblsm] = *Delta;
                        sttArr[tempMblsm] = eSttBfsPrcsb;
                        treeNumArr[tempMblsm]= treeNumArr[base];
                    }
                }
                s1 = ptr1Arr[tempB];
                blsmL[shiftB].push_back(tempM);
                blsmL[shiftB].push_back(s1);
                bs1=s1;
                if(blsmParent[s1]!=cNullItm)
                    bs1 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s1);
            }
            blsmChildren[shiftB].push_back(baseBlsm);
            std::reverse(blsmChildren[shiftB].begin(),blsmChildren[shiftB].end());
            std::reverse(blsmL[shiftB].begin(),blsmL[shiftB].end());
            while (baseBlsm != bs2)
            {
                tempB =rprsArr[bs2];
                tempM = mateArr[tempB];
                blsmL[shiftB].push_back(tempB);
                blsmL[shiftB].push_back(tempM);
                Size tempMblsm=tempM;
                if(blsmParent[tempM]!=cNullItm)
                    tempMblsm = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, tempM);
                blsmParent[bs2] = b;
                blsmChildren[shiftB].push_back(bs2);
                blsmParent[tempMblsm] = b;
                blsmChildren[shiftB].push_back(tempMblsm);
                blsmArr[tempB]=b;
                blsmArr[tempM]=b;
                if(tempMblsm < numVtxs)
                {
                    if(sttArr[tempMblsm] == eSttIdle)
                    {
                        //length[tempMblsm] = length[s2]-length[tempB]+pathlen+1;
                        //pathlen = length[tempMblsm] +1;
                        prcsbStk.Push(tempMblsm);
                        tOuter[tempMblsm] = *Delta;
                        sttArr[tempMblsm] = eSttBfsPrcsb;
                        treeNumArr[tempMblsm]= treeNumArr[base];
                    }
                }
                s2 = ptr1Arr[tempB];
                blsmL[shiftB].push_back(tempM);
                blsmL[shiftB].push_back(s2);
                bs2=s2;
                if(blsmParent[s2]!=cNullItm)
                    bs2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s2);
            }
            blsmParent[baseBlsm] = b;
            blsmArr[s1]=b;
            blsmArr[s2]=b;
            (*blsm)++;
        }

    template<class ItmQue, class ItmStk, class ItmIdxdQue>
        int MatchingEngine::rEdmondSearch(const Graph& graph,
                Size* mateArr, Size* ptrArr, Stt* sttArr, ItmQue& searchQue, ItmIdxdQue& expsdQue,
                Size* blsmArr,Size* rprsArr, Size* blsmParent,
                std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,  int* tOuter,
                int* tInner,int* Delta, Size* blsm, Size* trace, Size* treeNumArr,std::deque<std::pair<Size, Size> >* deltaQArr) const {

            Size numVtxs = graph.mNumVtxs;
            const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

            /*std::vector<Size> lengthVec;
              ResizeVector<Size>(&lengthVec, numVtxs);
              Size* length = (numVtxs == 0) ? NULL : &lengthVec[0];
              std::fill(&length[0], &length[numVtxs], 0);*/
            int largestD =0;
            std::queue<Size> tempQ;
            Size blsmbdual = *blsm;
            Size sFirst = expsdQue.First();

            if(sFirst == cNullItm)
            {
                return 0;
            }

            while (sFirst != cNullItm)
            {
                searchQue.Push(sFirst);
                sttArr[sFirst] = eSttBfsPrcsb;
                treeNumArr[sFirst]=sFirst;
                tOuter[sFirst] = 0;
                sFirst = expsdQue.Next(sFirst);
            }
            while(true)
            {
                while (searchQue.Empty() == false)
                {
                    Size s = searchQue.Front();
                    searchQue.Pop();
                    sttArr[s] = eSttBfsPrcsd;
                    int duals =1;

                    if(tInner[s]>-1)
                    {
                        if(tOuter[s]==-1)
                            duals = 1 + (*Delta - tInner[s]);
                        else
                            duals = 1 + (tOuter[s] - tInner[s]) - (*Delta - tOuter[s]);
                    }
                    else
                    {
                        if(tOuter[s]>-1)
                            duals = 1 - (*Delta - tOuter[s]);
                    }

                    Size sNumEdgs = vtxVecArr[s].size();
                    const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                    for(Size i =0; i<sNumEdgs; i++)
                    {
                        Size t = vtxArr[i];
                        Size ss = mateArr[t];
                        if (ss!=cNullItm && sttArr[ss]!= eSttIdle && (sttArr[t]== eSttIdle)) {
                            continue;
                        }
                        int dualt =1;
                        if(tInner[t] > -1)
                        {
                            if(tOuter[t]==-1)
                                dualt = 1 + (*Delta - tInner[t]);
                            else
                                dualt = 1 + (tOuter[t] - tInner[t]) - (*Delta - tOuter[t]);
                        }
                        else
                        {
                            if(tOuter[t]>-1)
                                dualt = 1 - (*Delta - tOuter[t]);
                        }
                        if(dualt + duals  == 0)
                        {
                            if(sttArr[t]== eSttIdle)
                            {
                                ptrArr[ss] = s;
                                //length[ss]=length[s]+2;
                                tInner[t]=*Delta;
                                tOuter[ss]=*Delta;
                                treeNumArr[ss]=treeNumArr[s];
                                sttArr[ss] = eSttBfsPrcsb;
                                searchQue.Push(ss);
                            }
                            else
                            {
                                if(treeNumArr[s] != treeNumArr[t])
                                {
                                    if(*blsm - 1 > blsmbdual)
                                    {
                                        int bbb = *blsm - 1;
                                        do
                                        {
                                            rprsArr[bbb]=cNullItm;
                                            rDissolveBlsm3Gabow<ItmQue>(graph, blsmArr,blsmParent, blsmChildren, bbb,blsmbdual);
                                            bbb--;
                                        }while(bbb > (int) blsmbdual);
                                    }
                                    for(Size j =numVtxs; j<=blsmbdual; j++)
                                    {
                                        if(blsmParent[j]==cNullItm)
                                        {
                                            tempQ.push(j-numVtxs);
                                            while (tempQ.empty() == false)
                                            {
                                                Size tempB = tempQ.front();
                                                tempQ.pop();
                                                for(Size i=0; i < blsmChildren[tempB].size();i++)
                                                {
                                                    Size child = blsmChildren[tempB][i];
                                                    if( child < numVtxs)
                                                    {
                                                        blsmArr[child]=j;
                                                    }
                                                    else
                                                    {
                                                        tempQ.push(child-numVtxs);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    for(int i =0; i<= largestD;i++)
                                    {
                                        if(deltaQArr[i].empty()==false)
                                            deltaQArr[i].clear();
                                    }
                                    return 1;
                                }
                                else
                                {

                                    Size b1=s;
                                    Size b2=t;
                                    if(blsmParent[s]!=cNullItm && blsmParent[t]!=cNullItm)
                                    {
                                        b1 =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                                        b2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                                        if( b1==b2 )
                                        {
                                            continue;
                                        }
                                    }
                                    rProcessBlsm5<ItmQue, ItmStk>(graph, mateArr,ptrArr,blsmArr, rprsArr, blsmParent,
                                            blsmChildren,blsmL,sttArr, searchQue,s, t, blsm, trace, treeNumArr,tOuter,Delta);
                                }
                            }
                        }
                        else
                        {
                            int d = duals+dualt;

                            if(sttArr[t] != eSttIdle)
                            {
                                d=d/2;
                            }
                            if(d+ (*Delta ) >  largestD)
                            {
                                largestD = d + (*Delta) ;
                            }
                            deltaQArr[*Delta+d].push_back(std::pair<Size, Size>(s,t));
                            continue;
                        }
                    }
                }
                Size l =0;
                if(deltaQArr[*Delta].empty())
                {
                    l =1;
                    while( l < numVtxs/2 +1 - *Delta)
                    {
                        if(deltaQArr[*Delta+l].empty()==false)
                            break;
                        l++;
                    }
                }
                if(deltaQArr[*Delta+l].empty()==false)
                {
                    if(l > 0)
                    {
                        (*Delta)+=l;
                        blsmbdual = *blsm - 1;
                    }
                    Size s;
                    Size t;
                    std::pair<Size, Size> edge = deltaQArr[*Delta].front();
                    s =edge.first;
                    t =edge.second;
                    deltaQArr[*Delta].pop_front();

                    Size ss = mateArr[t];
                    if (ss!=cNullItm && sttArr[ss]!= eSttIdle && sttArr[t]== eSttIdle) {
                        continue;
                    }

                    int duals =1;
                    int dualt =1;
                    if(tInner[s]==-1)
                    {
                        if(tOuter[s]!=-1)
                            duals = 1 - (*Delta - tOuter[s]);
                    }
                    else
                    {
                        if(tOuter[s]!=-1)
                            duals = 1 + (tOuter[s] - tInner[s]) - (*Delta - tOuter[s]);
                        else
                            duals = 1 + (*Delta - tInner[s]);
                    }
                    if(tInner[t]==-1)
                    {
                        if(tOuter[t]!=-1)
                            dualt = 1 - (*Delta - tOuter[t]);
                    }
                    else
                    {
                        if(tOuter[t]!=-1)
                            dualt = 1 + (tOuter[t] - tInner[t]) - (*Delta - tOuter[t]);
                        else
                            dualt = 1 + (*Delta - tInner[t]);
                    }
                    if(dualt + duals  != 0)
                    {
                        Size b11 =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                        Size b22 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                        if( b11==b22 )
                        {
                            continue;
                        }
                    }
                    if( sttArr[t]== eSttIdle)
                    {
                        Size ss = mateArr[t];
                        ptrArr[ss] = s;
                        //length[ss]=length[s]+2;
                        tInner[t]=*Delta;
                        tOuter[ss]=*Delta;
                        treeNumArr[ss]=treeNumArr[s];
                        sttArr[ss] = eSttBfsPrcsb;
                        searchQue.Push(ss);
                    }
                    else
                    {
                        if(treeNumArr[s] != treeNumArr[t])
                        {
                            if(*blsm - 1 > blsmbdual)
                            {
                                int bbb = *blsm - 1;
                                do
                                {
                                    rprsArr[bbb]=cNullItm;
                                    rDissolveBlsm3Gabow<ItmQue>(graph, blsmArr,blsmParent, blsmChildren, bbb,blsmbdual);
                                    bbb--;
                                }while(bbb > (int) blsmbdual);
                            }
                            for(Size j =numVtxs; j<=blsmbdual; j++)
                            {
                                if(blsmParent[j]==cNullItm)
                                {
                                    tempQ.push(j-numVtxs);
                                    while (tempQ.empty() == false)
                                    {
                                        Size tempB = tempQ.front();
                                        tempQ.pop();
                                        for(Size i=0; i < blsmChildren[tempB].size();i++)
                                        {
                                            Size child = blsmChildren[tempB][i];
                                            if( child < numVtxs)
                                                blsmArr[child]=j;
                                            else
                                                tempQ.push(child-numVtxs);
                                        }
                                    }
                                }
                            }
                            for(int i =0; i<= largestD;i++)
                            {
                                if(deltaQArr[i].empty()==false)
                                    deltaQArr[i].clear();
                            }
                            return 1;
                        }
                        else
                        {
                            Size b1=s;
                            Size b2=t;
                            if(blsmParent[s]!=cNullItm && blsmParent[t]!=cNullItm)
                            {
                                b1 =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                                b2 = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                            }
                            if( b1==b2 )
                            {
                                continue;
                            }
                            rProcessBlsm5<ItmQue, ItmStk>(graph, mateArr,ptrArr,blsmArr, rprsArr, blsmParent,blsmChildren,blsmL,sttArr, searchQue,s, t, blsm, trace,treeNumArr,tOuter,Delta);
                        }
                    }
                }
                else
                {
                    return 0;
                }

            }
        }

    template<class ItmQue, class ItmStk, class ItmIdxdQue>
        void MatchingEngine::rFindMxmlSetAugPathsCardGabow(const Graph& graph,
                Size* mateArr, Size* ptrArr, Stt* sttArr,
                Size* idxArr, ItmStk& prcsbStk, ItmIdxdQue& expsdQue,
                ItmQue& sLastQue, ItmQue& tLastQue, Size* blsmArr,Size* rprsArr, Size* blsmParent,
                std::vector<Size>* blsmChildren,std::vector<Size>* blsmL,  int* tOuter,int* tInner,int* Delta,
                Size* blsm, Size* blsmin, Size* treeNumArr) const {
            //Timer timer;
            Size numVtxs = graph.mNumVtxs;
            const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            std::queue<Size> tempQ;
            Size sFirst = expsdQue.First();
            while (sFirst != cNullItm)
            {
                if(sttArr[sFirst] != eSttIdle)
                {
                    sFirst = expsdQue.Next(sFirst);
                    continue;
                }

                prcsbStk.Clear();
                prcsbStk.Push(sFirst);
                treeNumArr[sFirst]=sFirst;
                sttArr[sFirst] = eSttDfsPrcsb;
                Size blsmf=sFirst;
                if(blsmParent[sFirst]!=cNullItm)
                    blsmf = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, sFirst);
                if(blsmf>=numVtxs)
                {
                    tempQ.push(blsmf-numVtxs);
                    while (tempQ.empty() == false)
                    {
                        Size tempB = tempQ.front();
                        tempQ.pop();
                        for(Size i=0; i < blsmChildren[tempB].size();i++)
                        {
                            Size child = blsmChildren[tempB][i];
                            if( child < numVtxs)
                            {
                                if(sttArr[child] == eSttIdle)
                                {
                                    prcsbStk.Push(child);
                                    sttArr[child] = eSttDfsPrcsb;
                                    treeNumArr[child]=sFirst;
                                }
                            }
                            else
                            {
                                tempQ.push(child-numVtxs);
                            }
                        }
                    }
                }
                while (prcsbStk.Empty() == false)
                {
                    Size s = prcsbStk.Top();
                    Size blsmt;
                    Size sNumEdgs = vtxVecArr[s].size();
                    const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                    while (idxArr[s] < sNumEdgs)
                    {
                        Size t = vtxArr[idxArr[s]];
                        Size ss;
                        ss = mateArr[t];
                        int duals =1;
                        int dualt =1;
                        if(tInner[s]==-1)
                        {
                            if(tOuter[s]!=-1)
                                duals = 1 - (*Delta - tOuter[s]);
                        }
                        else
                        {
                            if(tOuter[s]!=-1)
                                duals = 1 + (tOuter[s] - tInner[s]) - (*Delta - tOuter[s]);
                            else
                                duals = 1 + (*Delta - tInner[s]);
                        }

                        if(tInner[t]==-1)
                        {
                            if(tOuter[t]!=-1)
                                dualt = 1 - (*Delta - tOuter[t]);
                        }
                        else
                        {
                            if(tOuter[t]!=-1)
                                dualt = 1 + (tOuter[t] - tInner[t]) - (*Delta - tOuter[t]);
                            else
                                dualt = 1 + (*Delta - tInner[t]);
                        }
                        if(treeNumArr[s] == treeNumArr[t]  && sttArr[t] == eSttDfsPrcsd && dualt+duals == 0)
                        {
                            Size b1=s;
                            Size b2=t;
                            if(blsmParent[s]!=cNullItm && blsmParent[t]!=cNullItm)
                            {
                                b1 =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, s);
                                b2 =rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                            }
                            if( b1 != b2)
                            {
                                rProcessBlsm3Gabow<ItmQue, ItmStk>(graph, mateArr,ptrArr,blsmArr, rprsArr, blsmParent,blsmChildren,blsmL,sttArr, prcsbStk,t, s, blsm,blsmin,treeNumArr);
                                ++(idxArr[s]);
                                goto doneI;
                            }
                            ++(idxArr[s]);
                            continue;
                        }
                        if (((ss == cNullItm) && (sttArr[t] != eSttIdle)) || ((ss != cNullItm) && (sttArr[ss] != eSttIdle) )|| (treeNumArr[t] !=cNullItm && treeNumArr[s] != treeNumArr[t]) ) {
                            ++(idxArr[s]);
                            continue;
                        }
                        if(duals+dualt  == 0)
                        {
                            blsmt=t;
                            if(blsmParent[t]!=cNullItm)
                                blsmt = rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, t);
                            break;
                        }
                        else
                        {
                            ++(idxArr[s]);
                        }
                    }
                    if (idxArr[s] < sNumEdgs)
                    {

                        Size t = vtxArr[idxArr[s]];
                        Size ss = mateArr[t];
                        if(ss != cNullItm)
                        {
                            Size blsmss=ss;
                            if(blsmParent[ss]!=cNullItm)
                                blsmss=rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, ss);
                            if(sttArr[t] == eSttIdle)
                            {
                                if(blsmt >= numVtxs)
                                {
                                    Size base = rprsArr[blsmt];
                                    ss = mateArr[base];
                                    if(((ss == cNullItm) && (sttArr[base] != eSttIdle)) || ((ss != cNullItm) && (sttArr[ss] != eSttIdle)))
                                    {
                                        ++(idxArr[s]);
                                        continue;
                                    }
                                    blsmin[blsmt-numVtxs]=t;
                                    if(ss == cNullItm)
                                    {
                                        if(sttArr[base] == eSttIdle)
                                        {
                                            sLastQue.Push(s);
                                            tLastQue.Push(base);
                                            sttArr[base] = eSttLast;
                                            while (prcsbStk.Empty() == false)
                                            {
                                                Size sss = prcsbStk.Top();
                                                prcsbStk.Pop();
                                                sttArr[sss] = eSttDfsPrcsd;
                                            }
                                            break;
                                        }
                                        else
                                        {
                                            ++(idxArr[s]);
                                            continue;
                                        }
                                    }
                                    ptrArr[ss] = s;
                                    blsmss=ss;
                                    if(blsmParent[ss]!=cNullItm)
                                        blsmss=rFindblsm<ItmQue>(blsmParent,blsmArr, rprsArr, ss);
                                    if(blsmss >= numVtxs)
                                    {
                                        prcsbStk.Push(ss);
                                        sttArr[ss] = eSttDfsPrcsb;
                                        treeNumArr[ss] =treeNumArr[s];
                                        tempQ.push(blsmss-numVtxs);
                                        while (tempQ.empty() == false)
                                        {
                                            Size tempB = tempQ.front();
                                            tempQ.pop();
                                            for(Size i=0; i < blsmChildren[tempB].size();i++)
                                            {
                                                Size child = blsmChildren[tempB][i];
                                                if( child < numVtxs)
                                                {
                                                    if(sttArr[child] == eSttIdle)
                                                    {
                                                        prcsbStk.Push(child);
                                                        sttArr[child] = eSttDfsPrcsb;
                                                        ptrArr[child] = s;
                                                        treeNumArr[child] =treeNumArr[s];
                                                    }
                                                }
                                                else
                                                {
                                                    tempQ.push(child-numVtxs);
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        treeNumArr[ss] =treeNumArr[s];
                                        prcsbStk.Push(ss);
                                        sttArr[ss] = eSttDfsPrcsb;
                                    }
                                    ++(idxArr[s]);
                                    continue;
                                }
                                else if(blsmss >= numVtxs)
                                {
                                    ptrArr[ss] = s;
                                    prcsbStk.Push(ss);
                                    sttArr[ss] = eSttDfsPrcsb;
                                    treeNumArr[ss] =treeNumArr[s];
                                    tempQ.push(blsmss-numVtxs);
                                    while (tempQ.empty() == false)
                                    {
                                        Size tempB = tempQ.front();
                                        tempQ.pop();
                                        for(Size i=0; i < blsmChildren[tempB].size();i++)
                                        {
                                            Size child = blsmChildren[tempB][i];
                                            if( child < numVtxs)
                                            {
                                                if(sttArr[child] == eSttIdle)
                                                {
                                                    prcsbStk.Push(child);
                                                    sttArr[child] = eSttDfsPrcsb;
                                                    ptrArr[child] = s;
                                                    treeNumArr[child] =treeNumArr[s];
                                                }
                                            }
                                            else
                                            {
                                                tempQ.push(child-numVtxs);
                                            }
                                        }
                                    }
                                    ++(idxArr[s]);
                                    continue;
                                }
                                else
                                {
                                    treeNumArr[ss] = treeNumArr[s];
                                    ptrArr[ss] = s;
                                    prcsbStk.Push(ss);
                                    sttArr[ss] = eSttDfsPrcsb;
                                    ++(idxArr[s]);
                                }
                            }
                            else
                            {
                                ++(idxArr[s]);
                            }
                        }
                        else
                        {
                            sLastQue.Push(s);
                            tLastQue.Push(t);
                            sttArr[t] = eSttLast;
                            if(blsmt >=numVtxs)
                            {
                                blsmin[blsmt-numVtxs]=cNullItm;
                            }
                            while (prcsbStk.Empty() == false)
                            {
                                Size sss = prcsbStk.Top();
                                prcsbStk.Pop();
                                sttArr[sss] = eSttDfsPrcsd;
                            }
                            break;
                        }
                    }
                    else
                    {

                        prcsbStk.Pop();
                        sttArr[s] = eSttDfsPrcsd;
                    }
doneI:;
                }
                sFirst = expsdQue.Next(sFirst);
            }
        }

    template<class ItmQue, class ItmStk, class ItmIdxdQue>
        void MatchingEngine::rGabowCardMatching(const Graph& graph, Size* mateArr, Size* card) const {

            Size numVtxs = graph.mNumVtxs;
            if (mInlz == false) {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[numVtxs], cNullItm);
                }
                *card = 0;
            } else {
                Size iniCard(0);
                rInlzGrdyForCard(graph, mateArr, &iniCard);
                *card = iniCard;
            }

            ItmQue sLastQue(numVtxs);
            ItmQue tLastQue(numVtxs);

            ItmStk prcsbStk(numVtxs);
            ItmQue prcsbQue(numVtxs);

            ItmIdxdQue expsdQue(numVtxs);

            std::vector<Size> ptrVec;
            ResizeVector<Size>(&ptrVec,numVtxs);
            Size* ptrArr = (numVtxs == 0) ? NULL : &ptrVec[0];
            if (ptrArr != NULL) {
                std::fill(&ptrArr[0], &ptrArr[numVtxs], cNullItm);
            }

            std::vector<Size> blsmVec;
            ResizeVector<Size>(&blsmVec, numVtxs);
            Size* blsmArr = (numVtxs == 0) ? NULL : &blsmVec[0];

            std::vector<Size> rprsVec;
            ResizeVector<Size>(&rprsVec, numVtxs*2);
            Size* rprsArr = (numVtxs == 0) ? NULL : &rprsVec[0];
            for (Size v = 0; v < numVtxs; ++v) {
                rprsArr[v] = v;
                blsmArr[v] = v;
                expsdQue.Push(v);
            }
            std::fill(&rprsArr[numVtxs], &rprsArr[numVtxs*2], cNullItm);

            std::vector<Size> blsmParentVec;
            ResizeVector<Size>(&blsmParentVec, numVtxs*2);
            Size* blsmParentArr = (numVtxs == 0) ? NULL : &blsmParentVec[0];
            if (blsmParentArr != NULL) {
                std::fill(&blsmParentArr[0], &blsmParentArr[numVtxs*2], cNullItm);
            }
            std::vector<std::vector<Size> > blsmChildrenVec;
            ResizeVector<std::vector<Size> >(&blsmChildrenVec, numVtxs);
            std::vector<Size>* blsmChildrenArr = (numVtxs == 0) ? NULL : &blsmChildrenVec[0];

            std::vector<std::vector<Size> > blsmLVec;
            ResizeVector<std::vector<Size> >(&blsmLVec, numVtxs);
            std::vector<Size>* blsmL = (numVtxs == 0) ? NULL : &blsmLVec[0];

            std::vector<Size> blsminVec;
            ResizeVector<Size>(&blsminVec, numVtxs);
            Size* blsmin = (numVtxs == 0) ? NULL : &blsminVec[0];
            if (blsmin != NULL) {
                std::fill(&blsmin[0], &blsmin[numVtxs], cNullItm);
            }

            //const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            //const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

            std::vector<int> tInnerVec;
            ResizeVector<int>(&tInnerVec, numVtxs);
            int* tInner = (numVtxs == 0) ? NULL : &tInnerVec[0];
            if (tInner != NULL) {
                std::fill(&tInner[0], &tInner[numVtxs], -1);
            }
            std::vector<int> tOuterVec;
            ResizeVector<int>(&tOuterVec, numVtxs);
            int* tOuter = (numVtxs == 0) ? NULL : &tOuterVec[0];
            if (tOuter != NULL) {
                std::fill(&tOuter[0], &tOuter[numVtxs], 0);
            }
            std::vector<Stt> sttVec;
            ResizeVector<Stt>(&sttVec, numVtxs);
            Stt* sttArr = (numVtxs == 0) ? NULL : &sttVec[0];
            if (sttArr != NULL) {
                std::fill(&sttArr[0], &sttArr[numVtxs], eSttIdle);
            }
            std::vector<Size> idxVec;
            ResizeVector<Size>(&idxVec, numVtxs);
            Size* idxArr = (numVtxs == 0) ? NULL : &idxVec[0];
            if (idxArr != NULL) {
                std::fill(&idxArr[0], &idxArr[numVtxs], 0);
            }

            std::vector<Size> traceVec;
            ResizeVector<Size>(&traceVec, numVtxs);
            Size* trace = (numVtxs == 0) ? NULL : &traceVec[0];
            if (trace != NULL) {
                std::fill(&trace[0], &trace[numVtxs], cNullItm);
            }
            std::vector<Size> treeNumVec;
            ResizeVector<Size>(&treeNumVec, numVtxs);
            Size* treeNumArr = (numVtxs == 0) ? NULL : &treeNumVec[0];
            if (treeNumArr != NULL) {
                std::fill(&treeNumArr[0], &treeNumArr[numVtxs], cNullItm);
            }
            std::vector<std::deque<std::pair<Size, Size> > > deltaQVec;
            ResizeVector<std::deque<std::pair<Size, Size> > >(&deltaQVec, numVtxs/2+2);
            std::deque<std::pair<Size, Size> >* deltaQArr = (numVtxs == 0) ? NULL : &deltaQVec[0];
            Size blsm = numVtxs;
            int tNow=0;

            bool doneIter=true;

            tNow=1;
            rFindMxmlSetAugPathsCardGabow<ItmQue, ItmStk, ItmIdxdQue>
                (graph, mateArr, ptrArr, sttArr, idxArr,
                 prcsbStk, expsdQue, sLastQue, tLastQue,blsmArr,rprsArr,blsmParentArr,blsmChildrenArr,
                 blsmL, tOuter,tInner,&tNow, &blsm, blsmin,treeNumArr);

            while (tLastQue.Empty() == false)
            {
                Size sLast = sLastQue.Front();
                sLastQue.Pop();
                Size tLast = tLastQue.Front();
                tLastQue.Pop();
                sttArr[tLast] = eSttIdle;
                rAugmentGabow<ItmIdxdQue,ItmQue>(mateArr, ptrArr, expsdQue, sLast, tLast, blsmArr, rprsArr, blsmParentArr,  blsmChildrenArr,blsmL,numVtxs, blsmin);
                ++(*card);
            }

            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                blsmArr[v] = v;
                blsmParentArr[v]=cNullItm;
                trace[v]=cNullItm;
                sttArr[v]= eSttIdle;
                ptrArr[v]= cNullItm;
                tOuter[v]=-1;
                tInner[v]=-1;
                treeNumArr[v]= cNullItm;
            }
            tNow=0;
            do
            {
                doneIter=true;
                blsm = numVtxs;
                int augment = rEdmondSearch<ItmQue, ItmStk, ItmIdxdQue>
                    (graph, mateArr, ptrArr, sttArr,
                     prcsbQue, expsdQue,blsmArr,rprsArr,blsmParentArr,blsmChildrenArr,
                     blsmL, tOuter,tInner,&tNow, &blsm, trace,treeNumArr,deltaQArr);
                prcsbQue.Clear();
                if(augment == 0)
                    break;

                for(Size v =0 ; v < graph.mNumVtxs; v++)
                {
                    ptrArr[v]=cNullItm;
                    sttArr[v]=eSttIdle;
                    idxArr[v]= 0;
                    treeNumArr[v]= cNullItm;
                }

                rFindMxmlSetAugPathsCardGabow<ItmQue, ItmStk, ItmIdxdQue>
                    (graph, mateArr, ptrArr, sttArr, idxArr,
                     prcsbStk,  expsdQue, sLastQue, tLastQue,blsmArr,rprsArr,blsmParentArr,blsmChildrenArr,blsmL, tOuter,tInner,&tNow, &blsm, blsmin,treeNumArr);

                if(tLastQue.Empty()==false)
                    doneIter=false;
                while (tLastQue.Empty() == false)
                {
                    Size sLast = sLastQue.Front();
                    sLastQue.Pop();
                    Size tLast = tLastQue.Front();
                    tLastQue.Pop();
                    sttArr[tLast] = eSttIdle;
                    rAugmentGabow<ItmIdxdQue,ItmQue>(mateArr, ptrArr, expsdQue, sLast, tLast, blsmArr, rprsArr, blsmParentArr,  blsmChildrenArr,blsmL,numVtxs, blsmin);
                    ++(*card);
                }
                for (Size v = 0; v < graph.mNumVtxs; ++v) {
                    blsmArr[v] = v;
                    blsmParentArr[v]=cNullItm;
                    trace[v]=cNullItm;
                    sttArr[v]= eSttIdle;
                    ptrArr[v]= cNullItm;
                    tOuter[v]=-1;
                    tInner[v]=-1;
                    treeNumArr[v]= cNullItm;
                }
                tNow=0;

            }while(!doneIter);
        }

    template<class ItmQue>
        void MatchingEngine::rComputeMaxCardMatchingSglSrc1(const Graph& graph,
                Size* mateArr, Size* card) const {
            if (mInlz == false) {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
            } else {
                Size iniCard(0);
                //initialize with maximal cardinality
                rInlzGrdyForCard(graph, mateArr, &iniCard);
                *card = iniCard;
            }
            //to store a pointer from an outer vertex to its parent in the search tree
            std::vector<Size> ptr1Vec;
            ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
            Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
            if (ptr1Arr != NULL) {
                std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
            std::vector<Size> ptr2Vec;
            ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
            Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
            if (ptr2Arr != NULL) {
                std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
            }

            //to store a blossom id
            std::vector<Size> blsmVec;
            ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
            Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
            if (blsmArr != NULL) {
                std::fill(&blsmArr[0], &blsmArr[graph.mNumVtxs], cNullItm);
            }
            //to store a base id of a blossom
            std::vector<Size> rprsVec;
            ResizeVector<Size>(&rprsVec, graph.mNumVtxs);
            Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                rprsArr[v] = v;
            }
            //to store a status of a vertex
            std::vector<Stt> sttVec;
            ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
            Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
            if (sttArr != NULL) {
                std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
            }
            ItmQue prcsbQue(graph.mNumVtxs);
            ItmQue prcsdQue(graph.mNumVtxs);
            for (Size sFirst = 0; sFirst < graph.mNumVtxs; ++sFirst) {
                if (mateArr[sFirst] != cNullItm) {
                    continue;
                }
                Size sLast(cNullItm);
                Size tLast(cNullItm);
                // find an augmenting path
                rFindAugPathCardSglSrc1<ItmQue>
                    (graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr, sttArr, prcsbQue,
                     prcsdQue, sFirst, &sLast, &tLast);
                if (tLast != cNullItm) {
                    //augment an augmenting path
                    rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
                    ++(*card);
                }
                //reset working variable for visited vertices
                while (prcsbQue.Empty() == false) {
                    Size v = prcsbQue.Front();
                    prcsbQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    rprsArr[v] = v;
                    sttArr[v] = eSttIdle;
                }
                while (prcsdQue.Empty() == false) {
                    Size v = prcsdQue.Front();
                    prcsdQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    rprsArr[v] = v;
                    sttArr[v] = eSttIdle;
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rComputeMaxCardMatchingSglSrc2(const Graph& graph,
                Size* mateArr, Size* card) const {
            if (mInlz == false) {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
            } else {
                Size iniCard(0);
                //initialize with maximal cardinality
                rInlzGrdyForCard(graph, mateArr, &iniCard);
                *card = iniCard;
            }
            //Timer timer;
            //timer.Start();
            //to store a pointer from an outer vertex to its parent in the search tree
            std::vector<Size> ptr1Vec;
            ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
            Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
            if (ptr1Arr != NULL) {
                std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
            std::vector<Size> ptr2Vec;
            ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
            Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
            if (ptr2Arr != NULL) {
                std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a blossom id
            std::vector<Size> blsmVec;
            ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
            Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
            if (blsmArr != NULL) {
                std::fill(&blsmArr[0], &blsmArr[graph.mNumVtxs], cNullItm);
            }
            //to store a base id of a blossom
            std::vector<Size> rprsVec;
            ResizeVector<Size>(&rprsVec, graph.mNumVtxs);
            Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                rprsArr[v] = v;
            }
            //to store links between items in union find tree
            std::vector<Size> linkVec;
            ResizeVector<Size>(&linkVec, graph.mNumVtxs);
            Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                linkArr[v] = v;
            }
            //to store a rank of an item in union find tree
            std::vector<Size> rankVec;
            ResizeVector<Size>(&rankVec, graph.mNumVtxs);
            Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
            if (rankArr != NULL) {
                std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
            }
            //to store a status of a vertex
            std::vector<Stt> sttVec;
            ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
            Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
            if (sttArr != NULL) {
                std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
            }
            ItmQue prcsbQue(graph.mNumVtxs);
            ItmQue prcsdQue(graph.mNumVtxs);
            for (Size sFirst = 0; sFirst < graph.mNumVtxs; ++sFirst) {
                if (mateArr[sFirst] != cNullItm) {
                    continue;
                }
                Size sLast(cNullItm);
                Size tLast(cNullItm);
                // find an augmenting path
                rFindAugPathCardSglSrc2<ItmQue>
                    (graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr, linkArr, rankArr,
                     sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast);
                if (tLast != cNullItm) {
                    //augment an augmenting path
                    rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
                    ++(*card);
                    //reset working variable for visited vertices
                    while (prcsbQue.Empty() == false) {
                        Size v = prcsbQue.Front();
                        prcsbQue.Pop();
                        ptr1Arr[v] = cNullItm;
                        ptr2Arr[v] = cNullItm;
                        blsmArr[v] = cNullItm;
                        rprsArr[v] = v;
                        linkArr[v] = v;
                        rankArr[v] = 0;
                        sttArr[v] = eSttIdle;
                    }
                    while (prcsdQue.Empty() == false) {
                        Size v = prcsdQue.Front();
                        prcsdQue.Pop();
                        ptr1Arr[v] = cNullItm;
                        ptr2Arr[v] = cNullItm;
                        blsmArr[v] = cNullItm;
                        rprsArr[v] = v;
                        linkArr[v] = v;
                        rankArr[v] = 0;
                        sttArr[v] = eSttIdle;
                    }
                    continue;
                }
                //mark visited vertices so they are not visited again in future searches
                while (prcsbQue.Empty() == false) {
                    Size v = prcsbQue.Front();
                    prcsbQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    rprsArr[v] = v;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    sttArr[v] = eSttLast;
                }
                while (prcsdQue.Empty() == false) {
                    Size v = prcsdQue.Front();
                    prcsdQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    rprsArr[v] = v;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    sttArr[v] = eSttLast;
                }
            }
        }


    template<class ItmQue>
        void MatchingEngine::rComputeMaxCardMatchingSglSrc3(const Graph& graph,
                Size* mateArr, Size* card,std::pair<Size,Val>* sortedArr) const {
            if (mInlz == false) {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
            } else {
                //initialize with maximal cardinality
                Size iniCard(0);
                rInlzGrdyForCard(graph, mateArr, &iniCard);
                *card = iniCard;
            }
            //Timer timer;
            //timer.Start();
            //to store a pointer from an outer vertex to its parent in the search tree
            std::vector<Size> ptr1Vec;
            ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
            Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
            if (ptr1Arr != NULL) {
                std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
            std::vector<Size> ptr2Vec;
            ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
            Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
            if (ptr2Arr != NULL) {
                std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a blossom id
            std::vector<Size> blsmVec;
            ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
            Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
            if (blsmArr != NULL) {
                std::fill(&blsmArr[0], &blsmArr[graph.mNumVtxs], cNullItm);
            }
            //to store a base id of a blossom
            std::vector<Size> rprsVec;
            ResizeVector<Size>(&rprsVec, graph.mNumVtxs);
            Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                rprsArr[v] = v;
            }
            //to store links between items in union find tree
            std::vector<Size> linkVec;
            ResizeVector<Size>(&linkVec, graph.mNumVtxs);
            Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                linkArr[v] = v;
            }
            //to store a rank of an item in union find tree
            std::vector<Size> rankVec;
            ResizeVector<Size>(&rankVec, graph.mNumVtxs);
            Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
            if (rankArr != NULL) {
                std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
            }
            //to store a status of a vertex
            std::vector<Stt> sttVec;
            ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
            Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
            if (sttArr != NULL) {
                std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
            }
            ItmQue prcsbQue(graph.mNumVtxs);
            ItmQue prcsdQue(graph.mNumVtxs);
            for (Size current = 0; current <graph.mNumVtxs; ++current) {
                Size sFirst =  sortedArr[current].first;
                if (mateArr[sFirst] != cNullItm) {
                    continue;
                }
                Size sLast(cNullItm);
                Size tLast(cNullItm);
                // find an augmenting path
                rFindAugPathCardSglSrc2<ItmQue>
                    (graph, mateArr, ptr1Arr, ptr2Arr, blsmArr, rprsArr, linkArr, rankArr,
                     sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast);
                if (tLast != cNullItm) {
                    //augment an augmenting path
                    rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
                    ++(*card);
                    //reset visited vertices working variables
                    while (prcsbQue.Empty() == false) {
                        Size v = prcsbQue.Front();
                        prcsbQue.Pop();
                        ptr1Arr[v] = cNullItm;
                        ptr2Arr[v] = cNullItm;
                        blsmArr[v] = cNullItm;
                        rprsArr[v] = v;
                        linkArr[v] = v;
                        rankArr[v] = 0;
                        sttArr[v] = eSttIdle;
                    }
                    while (prcsdQue.Empty() == false) {
                        Size v = prcsdQue.Front();
                        prcsdQue.Pop();
                        ptr1Arr[v] = cNullItm;
                        ptr2Arr[v] = cNullItm;
                        blsmArr[v] = cNullItm;
                        rprsArr[v] = v;
                        linkArr[v] = v;
                        rankArr[v] = 0;
                        sttArr[v] = eSttIdle;
                    }
                    continue;
                }
                //mark visited vertices so they are not visited again in future searches
                while (prcsbQue.Empty() == false) {
                    Size v = prcsbQue.Front();
                    prcsbQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    rprsArr[v] = v;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    sttArr[v] = eSttLast;
                }
                while (prcsdQue.Empty() == false) {
                    Size v = prcsdQue.Front();
                    prcsdQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    rprsArr[v] = v;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    sttArr[v] = eSttLast;
                }
            }
        }

    template<class ItmQue>
        void MatchingEngine::rComputeAprxMaxEdgWghtMatching4(const Graph& graph,
                Size* mateArr, Size* card, Val* edgWght) const {
            if (mateArr != NULL) {
                std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
            }
            *card = 0;
            *edgWght = 0.0;
            const std::vector<Size>*
                vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const std::vector<Val>*
                edgWghtVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mEdgWghtVecVec[0];
            // candidates.
            std::vector<Size> cndtMateVec;
            ResizeVector<Size>(&cndtMateVec, graph.mNumVtxs);
            Size* cndtMateArr = (graph.mNumVtxs == 0) ? NULL : &cndtMateVec[0];
            if (cndtMateArr != NULL) {
                std::fill(&cndtMateArr[0], &cndtMateArr[graph.mNumVtxs], cNullItm);
            }
            // queue of matched vertices to be processed.
            ItmQue matchedQue(graph.mNumVtxs);
            // initialize the queue of matched vertices.
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                Size w = cNullItm;
                Val maxEdgWght = cNegInfVal;
                Size vNumEdgs = vtxVecArr[v].size();
                const Size* vVtxArr = (vNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
                const Val* vEdgWghtArr = (vNumEdgs == 0) ? NULL : &edgWghtVecArr[v][0];
                for (Size j = 0; j < vNumEdgs; ++j) {
                    Size z = vVtxArr[j];
                    Val vzEdgWght = vEdgWghtArr[j];
                    if (z == v) {
                        continue;
                    }
                    if (maxEdgWght < vzEdgWght) {
                        w = z;
                        maxEdgWght = vzEdgWght;
                    }
                }
                cndtMateArr[v] = w;
                if ((w != cNullItm) && (cndtMateArr[w] == v)) {
                    mateArr[v] = w;
                    mateArr[w] = v;
                    matchedQue.Push(v);
                    ++(*card);
                    *edgWght += maxEdgWght;
                }
            }
            // only one vertex out of each matched pair is queued, use an array to
            // process both.
            Size matchedArr[2];
            // process the matched vertices from the queue.
            while (matchedQue.Empty() == false) {
                Size u = matchedQue.Front();
                matchedQue.Pop();
                matchedArr[0] = u;
                matchedArr[1] = mateArr[u];
                for (Size k = 0; k < 2; ++k) {
                    u = matchedArr[k];
                    Size uNumEdgs = vtxVecArr[u].size();
                    const Size* uVtxArr = (uNumEdgs == 0) ? NULL : &vtxVecArr[u][0];
                    for (Size i = 0; i < uNumEdgs; ++i) {
                        Size v = uVtxArr[i];
                        if (v == u) {
                            continue;
                        }
                        if (mateArr[v] != cNullItm) {
                            continue;
                        }
                        if (cndtMateArr[v] != u) {
                            continue;
                        }
                        Size w = cNullItm;
                        Val maxEdgWght = cNegInfVal;
                        Size vNumEdgs = vtxVecArr[v].size();
                        const Size* vVtxArr = (vNumEdgs == 0) ? NULL : &vtxVecArr[v][0];
                        const Val* vEdgWghtArr = (vNumEdgs == 0) ? NULL : &edgWghtVecArr[v][0];
                        for (Size j = 0; j < vNumEdgs; ++j) {
                            Size z = vVtxArr[j];
                            Val vzEdgWght = vEdgWghtArr[j];
                            if (z == v) {
                                continue;
                            }
                            if (mateArr[z] != cNullItm) {
                                continue;
                            }
                            if (maxEdgWght < vzEdgWght) {
                                w = z;
                                maxEdgWght = vzEdgWght;
                            }
                        }
                        cndtMateArr[v] = w;
                        if ((w != cNullItm) && (cndtMateArr[w] == v)) {
                            mateArr[v] = w;
                            mateArr[w] = v;
                            matchedQue.Push(v);
                            ++(*card);
                            *edgWght += maxEdgWght;
                        }
                    }
                }
            }
        }

    //max vertex weighted matching using augmenting paths that reach a heaviest vertex
    template<class ItmQue, class ItmStk>
        void MatchingEngine::rComputeMaxVtxWght(const Graph& graph,
                Size* mateArr, Size* card, Val* vtxWght, Size augPathLenBnd) const {

            if (mateArr != NULL) {
                std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
            }
            *card = 0;
            Size numVtxs = graph.mNumVtxs;

            //to store a pointer from an outer vertex to its parent in the search tree
            std::vector<Size> ptr1Vec;
            ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
            Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
            if (ptr1Arr != NULL) {
                std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
            std::vector<Size> ptr2Vec;
            ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
            Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
            if (ptr2Arr != NULL) {
                std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
            }
            //to store a blossom id
            std::vector<Size> blsmVec;
            ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
            Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
            if (blsmArr != NULL) {
                for (Size v = 0; v < graph.mNumVtxs; ++v) {
                    blsmArr[v] = v;
                }
            }
            //to store a base id of a blossom
            std::vector<Size> rprsVec;
            ResizeVector<Size>(&rprsVec, graph.mNumVtxs*2);
            Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                rprsArr[v] = v;
            }
            std::fill(&rprsArr[numVtxs], &rprsArr[numVtxs*2], cNullItm);
            //to store a status of a vertex
            std::vector<Stt> sttVec;
            ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
            Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
            if (sttArr != NULL) {
                std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
            }
            //to store links between items in union find tree
            std::vector<Size> linkVec;
            ResizeVector<Size>(&linkVec, graph.mNumVtxs);
            Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
            for (Size v = 0; v < graph.mNumVtxs; ++v) {
                linkArr[v] = v;
            }
            //to store a rank of an item in union find tree
            std::vector<Size> rankVec;
            ResizeVector<Size>(&rankVec, graph.mNumVtxs);
            Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
            if (rankArr != NULL) {
                std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
            }


            ItmQue prcsbQue(graph.mNumVtxs);
            ItmQue prcsdQue(graph.mNumVtxs);
            //ItmStk prcsdQue(graph.mNumVtxs);
            //const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
            const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

            //std::deque<std::pair<Size, Val> > expsdLst;
            // for (Size first = 0; first < numVtxs; ++first) {
            //  expsdLst.push_back(std::pair<Size, Val>(first, vtxWghtArr[first]));
            //	}
            //expsdLst.sort(ValGreater<std::pair<Size, Val> >());

            //Size nedges_profile =0;
            //Size nblossom_profile =0;
            //Size ntedges=0;
            Size nedges=0;
            std::vector<std::pair<Size, Val> > sExpsdLst2;
            ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
            std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0];

            for (Size first = 0; first < numVtxs; ++first) {
                sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]);
            }
            //sort vertices in non-increasing order of their weights
            std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater2<std::pair<Size, Val> >());

            //Timer timer;
            //for (std::deque<std::pair<Size, Val> >::const_iterator it = expsdLst.begin();it != expsdLst.end(); ++it) {
            Val talen =0;
            for (Size current = 0; current < numVtxs; ++current) {
                Size sFirst =  sExpsdLst2Arr[current].first;
                //Size sFirst = it->first;
                //std::cout << "this is iteration " <<current << std::endl;
                //nedges=0;
                if(*card > 0)
                {
                    std::cout << "no of visited edges " <<nedges << " no of aug path " <<*card<<" sum of lengths " <<talen<<" avg lengths " <<talen/(*card)<<std::endl;
                }
                if(mateArr[sFirst]!=cNullItm)
                    continue;
                //find the next heaviest vertex
                Val secondWght =sExpsdLst2Arr[current].second;

                if(current+1<numVtxs)
                {
                    Size c = current+1;
                    while(c<numVtxs)
                    {
                        if(mateArr[sExpsdLst2Arr[c].first]==cNullItm)
                        {
                            secondWght = sExpsdLst2Arr[c].second;
                            break;
                        }
                        c++;
                    }
                }
                else
                    break;
                Size sLast(cNullItm);
                Size tLast(cNullItm);
                Val calen =0;
                //timer.Start();
                //find an augmenting path to a heaviest unmatched vertex
                rFindAugPathCardSglSrc3<ItmQue>
                    (graph, mateArr, ptr1Arr,  ptr2Arr,blsmArr, rprsArr,linkArr, rankArr,
                     sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast,secondWght,&nedges);


                if (tLast != cNullItm) {
                    //augment an augmenting path
                    rAugment(mateArr, ptr1Arr, ptr2Arr, sLast, tLast,&calen);
                    calen=calen*2 -1;
                    talen+=calen;
                    ++(*card);
                    //reset visited vertices working variables
                    while (prcsbQue.Empty() == false) {
                        Size v = prcsbQue.Front();
                        prcsbQue.Pop();
                        ptr1Arr[v] = cNullItm;
                        ptr2Arr[v] = cNullItm;
                        blsmArr[v] = cNullItm;
                        linkArr[v] = v;
                        rankArr[v] = 0;
                        rprsArr[v] = v;
                        sttArr[v] = eSttIdle;
                    }

                    while (prcsdQue.Empty() == false) {
                        Size v = prcsdQue.Front();
                        prcsdQue.Pop();
                        ptr1Arr[v] = cNullItm;
                        ptr2Arr[v] = cNullItm;
                        blsmArr[v] = cNullItm;
                        linkArr[v] = v;
                        rankArr[v] = 0;
                        rprsArr[v] = v;
                        sttArr[v] = eSttIdle;
                    }
                    continue;
                }
                //mark visited vertices so they are not visited again in future searches
                while (prcsbQue.Empty() == false) {
                    Size v = prcsbQue.Front();
                    prcsbQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    rprsArr[v] = v;
                    sttArr[v] = eSttLast;
                }

                while (prcsdQue.Empty() == false) {
                    Size v = prcsdQue.Front();
                    prcsdQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    rprsArr[v] = v;
                    sttArr[v] = eSttLast;
                }

            }
            *vtxWght = 0.0;
            for (Size s = 0; s < numVtxs; ++s) {
                if (mateArr[s] != cNullItm) {
                    *vtxWght += vtxWghtArr[s];
                }
            }
        }

        //Max vertex weighted matching using max cardinality then find increasing paths
        /*template<class ItmQue, class ItmStk>
          void MatchingEngine::rComputeMaxVtxWght(const Graph& graph,
          Size* mateArr, Size* card, Val* vtxWght, Size augPathLenBnd) const {
          if (mInlz == false) {
          if (mateArr != NULL) {

          std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
          }
         *card = 0;
         } else {

         Size iniCard(0);
        //initialize with maximal cardinality
        rInlzGrdyForCard(graph, mateArr, &iniCard);
         *card = iniCard;
         }



         Size numVtxs = graph.mNumVtxs;
        //to store a pointer from an outer vertex to its parent in the search tree
        std::vector<Size> ptr1Vec;
        ResizeVector<Size>(&ptr1Vec, graph.mNumVtxs);
        Size* ptr1Arr = (graph.mNumVtxs == 0) ? NULL : &ptr1Vec[0];
        if (ptr1Arr != NULL) {
        std::fill(&ptr1Arr[0], &ptr1Arr[graph.mNumVtxs], cNullItm);
        }
        //to store a blossom id
        std::vector<Size> blsmVec;
        ResizeVector<Size>(&blsmVec, graph.mNumVtxs);
        Size* blsmArr = (graph.mNumVtxs == 0) ? NULL : &blsmVec[0];
        if (blsmArr != NULL) {
        for (Size v = 0; v < graph.mNumVtxs; ++v) {
        blsmArr[v] = v;
        }
        }
        //to store a base id of a blossom
        std::vector<Size> rprsVec;
        ResizeVector<Size>(&rprsVec, graph.mNumVtxs*2);
        Size* rprsArr = (graph.mNumVtxs == 0) ? NULL : &rprsVec[0];
        for (Size v = 0; v < graph.mNumVtxs; ++v) {
        rprsArr[v] = v;
        }
        std::fill(&rprsArr[numVtxs], &rprsArr[numVtxs*2], cNullItm);
        //to store a status of a vertex
        std::vector<Stt> sttVec;
        ResizeVector<Stt>(&sttVec, graph.mNumVtxs);
        Stt* sttArr = (graph.mNumVtxs == 0) ? NULL : &sttVec[0];
        if (sttArr != NULL) {
        std::fill(&sttArr[0], &sttArr[graph.mNumVtxs], eSttIdle);
        }
        //to store a pointer in from an outer vertex to its parent after its status is changed in a blossom
        std::vector<Size> ptr2Vec;
        ResizeVector<Size>(&ptr2Vec, graph.mNumVtxs);
        Size* ptr2Arr = (graph.mNumVtxs == 0) ? NULL : &ptr2Vec[0];
        if (ptr2Arr != NULL) {
        std::fill(&ptr2Arr[0], &ptr2Arr[graph.mNumVtxs], cNullItm);
        }
        //to store links between items in union find tree
        std::vector<Size> linkVec;
        ResizeVector<Size>(&linkVec, graph.mNumVtxs);
        Size* linkArr = (graph.mNumVtxs == 0) ? NULL : &linkVec[0];
        for (Size v = 0; v < graph.mNumVtxs; ++v) {
        linkArr[v] = v;
        }
        //to store a rank of an item in union find tree
        std::vector<Size> rankVec;
        ResizeVector<Size>(&rankVec, graph.mNumVtxs);
        Size* rankArr = (graph.mNumVtxs == 0) ? NULL : &rankVec[0];
        if (rankArr != NULL) {
        std::fill(&rankArr[0], &rankArr[graph.mNumVtxs], 0);
        }

        ItmQue prcsbQue(graph.mNumVtxs);
        ItmQue prcsdQue(graph.mNumVtxs);
        ItmQue unmatchedQue(graph.mNumVtxs);
        const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

        std::vector<std::pair<Size, Val> > sExpsdLst2;
        ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
        std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0];

        for (Size first = 0; first < numVtxs; ++first) {
            sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]);
        }

        //sort vertices in non-increasing order of their weights
        std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater<std::pair<Size, Val> >());
        // compute maximum cardinality
        rComputeMaxCardMatchingSglSrc3<ItmQue>(graph, mateArr, card,sExpsdLst2Arr);
        //Timer timer;
        for (Size current = 0; current < numVtxs; ++current) {
            Size sFirst =  sExpsdLst2Arr[current].first;
            if(mateArr[sFirst]!=cNullItm)
                continue;
            unmatchedQue.Push(sFirst);
        }

        while(unmatchedQue.Empty()==false)
        {

            Size sFirst = unmatchedQue.Front();
            unmatchedQue.Pop();
            Size sLast(cNullItm);
            Size tLast(cNullItm);


            //timer.Start();
            // find an increasing path that reach a lightest matched vertex
            rFindIncPathCardSglSrc3<ItmQue>
                (graph, mateArr, ptr1Arr,  ptr2Arr,blsmArr, rprsArr,linkArr, rankArr,
                 sttArr, prcsbQue, prcsdQue, sFirst, &sLast, &tLast);

            if (tLast != cNullItm) {
                //reverse an increasing path
                fix(mateArr, ptr1Arr, ptr2Arr, sLast, tLast);
                //reset visited vertices working variables
                while (prcsbQue.Empty() == false) {
                    Size v = prcsbQue.Front();
                    prcsbQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    rprsArr[v] = v;
                    sttArr[v] = eSttIdle;
                }

                while (prcsdQue.Empty() == false) {
                    Size v = prcsdQue.Front();
                    prcsdQue.Pop();
                    ptr1Arr[v] = cNullItm;
                    ptr2Arr[v] = cNullItm;
                    blsmArr[v] = cNullItm;
                    linkArr[v] = v;
                    rankArr[v] = 0;
                    rprsArr[v] = v;
                    sttArr[v] = eSttIdle;
                }

                continue;
            }
            //mark visited vertices so they are not visited again in future searches
            while (prcsbQue.Empty() == false) {
                Size v = prcsbQue.Front();
                prcsbQue.Pop();
                ptr1Arr[v] = cNullItm;
                ptr2Arr[v] = cNullItm;
                blsmArr[v] = cNullItm;
                rprsArr[v] = v;
                linkArr[v] = v;
                rankArr[v] = 0;
                sttArr[v] = eSttLast;
            }
            while (prcsdQue.Empty() == false) {
                Size v = prcsdQue.Front();
                prcsdQue.Pop();
                ptr1Arr[v] = cNullItm;
                ptr2Arr[v] = cNullItm;
                blsmArr[v] = cNullItm;
                rprsArr[v] = v;
                linkArr[v] = v;
                rankArr[v] = 0;
                sttArr[v] = eSttLast;
            }
        }

        *vtxWght = 0.0;
        for (Size s = 0; s < numVtxs; ++s) {
            if (mateArr[s] != cNullItm) {
                *vtxWght += vtxWghtArr[s];
            }
        }
        }*/

        template<class ItmQue>
            void MatchingEngine::rSuitorVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
                Size numVtxs = graph.mNumVtxs;
                Size se=0;
                std::vector<Val> sWs;
                ResizeVector<Val>(&sWs, graph.mNumVtxs);
                Val* ws = (graph.mNumVtxs == 0) ? NULL : &sWs[0];
                if (ws != NULL) {
                    std::fill(&ws[0], &ws[graph.mNumVtxs], cZeroVal);
                }

                ////
                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
                for (Size first = 0; first < numVtxs; ++first)
                {
                    Size current = first;
                    bool done = false;
                    while(!done)
                    {
                        Size partner = cNullItm;
                        Val heaviest = 0.0;

                        Size sNumEdgs = vtxVecArr[current].size();
                        const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[current][0];
                        //find a heaviest neighbor with a lower or no proposal
                        for (Size i = 0; i < sNumEdgs; ++i) {
                            Size t = currentVtxArr[i];
                            se++;
                            if (vtxWghtArr[t] + vtxWghtArr[current] >  heaviest || ((fabs(vtxWghtArr[t] + vtxWghtArr[current] -  heaviest)< 0.0000001) && partner > t)){
                                if( vtxWghtArr[t] + vtxWghtArr[current] > ws[t])
                                {
                                    partner = t;
                                    heaviest = vtxWghtArr[t] + vtxWghtArr[current];
                                }
                                else if (mateArr[t] !=cNullItm && fabs((vtxWghtArr[t] + vtxWghtArr[current] )- ws[t]) < 0.0000001 && mateArr[t]>current)
                                {
                                    partner = t;
                                    heaviest = vtxWghtArr[t] + vtxWghtArr[current];
                                }

                            }
                        }
                        done = true;
                        if(heaviest > 0.0)
                        {
                            //change partner and annul y proposal is any
                            Size y = mateArr[partner];
                            mateArr[partner]=current;
                            ws[partner]=heaviest;
                            if(y != cNullItm)
                            {
                                current = y;
                                done= false;
                            }
                        }
                    }
                }

                std::cout <<se<<" ";
                *vtxWght = 0.0;
                for (Size s = 0; s < numVtxs; ++s) {
                    if (mateArr[s] != cNullItm) {
                        *vtxWght += vtxWghtArr[s];
                        ++(*card);
                    }
                }
                (*card)=(*card)/2;
            }

        template<class ItmQue>
            void MatchingEngine::rParSuitorVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {

                Size numVtxss = graph.mNumVtxs;
                std::vector<Val> sWs;
                ResizeVector<Val>(&sWs, graph.mNumVtxs);
                Val* ws = (graph.mNumVtxs == 0) ? NULL : &sWs[0];
                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
                Size nlocks[numVtxss];
                Val tempvtxWght = 0.0;
                Size tempcard=0;
#pragma omp parallel
                {
                    Size numVtxs = numVtxss;
                    //setting locks for each vertex in parallel
#pragma omp for schedule(static,256)
                    for(Size i=0;i<numVtxs;i++)
                    {
                        nlocks[i]=0;
                        ws[i]=cZeroVal;
                        mateArr[i]=cNullItm;
                    }
                    // in parallel find a heaviest partner with a lower or no proposal
#pragma omp for schedule(static,256)
                    for (Size first = 0; first < numVtxs; ++first)
                    {
                        Size current = first;
                        bool done = false;
                        while(!done)
                        {
                            Size partner = cNullItm;
                            Val heaviest = 0;

                            Size sNumEdgs = vtxVecArr[current].size();
                            const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[current][0];
                            for (Size i = 0; i < sNumEdgs; ++i) {
                                Size t = currentVtxArr[i];
                                //while(__sync_lock_test_and_set(&nlocks[t], 1));
                                if (vtxWghtArr[t] + vtxWghtArr[current] >  heaviest || ((fabs(vtxWghtArr[t] + vtxWghtArr[current] -  heaviest)< 0.0000001) && partner > t)){
                                    if( vtxWghtArr[t] + vtxWghtArr[current] > ws[t])
                                    {
                                        partner = t;
                                        heaviest = vtxWghtArr[t] + vtxWghtArr[current];
                                    }
                                    else if (mateArr[t] !=cNullItm && fabs((vtxWghtArr[t] + vtxWghtArr[current] )- ws[t]) < 0.0000001 && mateArr[t]>current)
                                    {
                                        partner = t;
                                        heaviest = vtxWghtArr[t] + vtxWghtArr[current];
                                    }

                                }
                                //__sync_lock_release(&nlocks[t]);
                            }
                            done = true;
                            if(heaviest > 0.0)
                            {
                                //lock partner
                                while(__sync_lock_test_and_set(&nlocks[partner], 1))
                                    ;
                                Size y = mateArr[partner];
                                if(y != cNullItm)
                                {
                                    if(vtxWghtArr[y] < vtxWghtArr[current] || (vtxWghtArr[y] == vtxWghtArr[current] && current < y))
                                    {
                                        mateArr[partner]=current;
                                        ws[partner]=heaviest;
                                        current = y;
                                        done= false;
                                    }
                                    else
                                    {
                                        done =false;
                                    }
                                }
                                else
                                {
                                    mateArr[partner]=current;
                                    ws[partner]=heaviest;
                                }
                                //release the lock
                                __sync_lock_release(&nlocks[partner]);
                            }
                        }
                    }

#pragma omp for schedule(static,256) nowait reduction(+:tempvtxWght,tempcard)
                    for (Size s = 0; s < numVtxs; ++s)
                    {
                        if (mateArr[s] != cNullItm) {
                            tempvtxWght += vtxWghtArr[s];
                            ++(tempcard);
                        }
                    }
                }
                *vtxWght =tempvtxWght;
                *card=tempcard/2;
            }

        template<class ItmQue>
            void MatchingEngine::rPGPVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {

                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
                Val vtxWght1=0.0;
                Val vtxWght2=0.0;
                Size card1=0;
                Size card2=0;
                Size numVtxs = graph.mNumVtxs;
                //paths counter
                Size pc=0;
                //mates of vertices for M1
                std::vector<Size> mateVec1;
                ResizeVector<Size>(&mateVec1, graph.mNumVtxs);
                Size* mateArr1 = (graph.mNumVtxs == 0) ? NULL : &mateVec1[0];
                if (mateArr1 != NULL) {
                    std::fill(&mateArr1[0], &mateArr1[graph.mNumVtxs], cNullItm);
                }
                //mate of vertices for M2
                std::vector<Size> mateVec2;
                ResizeVector<Size>(&mateVec2, graph.mNumVtxs);
                Size* mateArr2 = (graph.mNumVtxs == 0) ? NULL : &mateVec2[0];
                if (mateArr2 != NULL) {
                    std::fill(&mateArr2[0], &mateArr2[graph.mNumVtxs], cNullItm);
                }
                //marking the visitied vertices
                std::vector<Size> markVec;
                ResizeVector<Size>(&markVec, graph.mNumVtxs);
                Size* mark = (graph.mNumVtxs == 0) ? NULL : &markVec[0];
                if (mark != NULL) {
                    std::fill(&mark[0], &mark[graph.mNumVtxs], 0);
                }
                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
                std::vector<std::vector<Size> > tmpVtxVecVec;
                ResizeVector<std::vector<Size> >(&tmpVtxVecVec, numVtxs);
                std::vector<Size>* tmpVtxVecArr = (numVtxs == 0) ? NULL : &tmpVtxVecVec[0];
                for (Size first = 0; first < numVtxs; ++first){
                    Size sNumEdgs = vtxVecArr[first].size();
                    for (Size i = 0; i < sNumEdgs; ++i)
                    {
                        tmpVtxVecArr[first].push_back(0);
                    }
                }
                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
                //to alternate between M1 and M2
                bool M1turn =true;
                bool cycle =false;
                //length of the path
                Size lp=0;
                Size nedges_profile=0;
                for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
                    Size first = vfirst;
                    cycle=false;
                    pc++;
                    mark[first]=pc;
                    lp=0;
                    M1turn=true;
                    bool done = false;
                    // construct a path from first vertex
                    while(!done)
                    {
                        // the path could not be extended
                        if(M1turn && mateArr1[first]!=cNullItm)
                            break;
                        if(!M1turn && mateArr2[first]!=cNullItm)
                            break;
                        Size sNumEdgs = vtxVecArr[first].size();
                        Size second=cNullItm;
                        Val heaviest = 0.0;
                        Size eIndex=cNullItm;
                        const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[first][0];
                        // find a heaviest neighbor to extend a path
                        for (Size i = 0; i < sNumEdgs; ++i) {
                            nedges_profile++;
                            Size tempSecond = currentVtxArr[i];

                            assert(first != tempSecond);
                            // if a vertex is already visited
                            if(mark[tempSecond]==pc)
                            {
                                if(second==vfirst && lp %2 ==1)
                                    cycle=true;
                                else
                                    continue;
                            }
                            // already constructed path
                            if(mark[tempSecond]<pc && mark[tempSecond]>0)
                                continue;
                            // save a neighbor if it is heavier than the previous one
                            if (vtxWghtArr[tempSecond]> heaviest && tmpVtxVecArr[first][i] ==0 ) {
                                second= tempSecond;
                                heaviest = vtxWghtArr[tempSecond];
                                eIndex=i;
                            }
                        }

                        done =true;
                        if(heaviest >0.0)
                        {
                            mark[second]=pc;
                            // increment the length
                            lp++;
                            if(M1turn)
                            {
                                if (mateArr1[second] == cNullItm)
                                {
                                    mateArr1[second]=first;
                                    mateArr1[first]=second;
                                    vtxWght1+=vtxWghtArr[first];
                                    vtxWght1+=vtxWghtArr[second];
                                    card1++;
                                    tmpVtxVecArr[first][eIndex]=1;
                                    Size sNumEdgs1 = vtxVecArr[second].size();
                                    const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
                                    for (Size j = 0; j < sNumEdgs1; ++j) {
                                        nedges_profile++;
                                        if(currentVtxArr1[j]==first)
                                        {
                                            tmpVtxVecArr[second][j]=1;
                                            break;
                                        }
                                    }
                                    M1turn=false;
                                }
                            }
                            else
                            {
                                if (mateArr2[second] == cNullItm)
                                {
                                    mateArr2[second]=first;
                                    mateArr2[first]=second;
                                    vtxWght2+=vtxWghtArr[first];
                                    vtxWght2+=vtxWghtArr[second];
                                    card2++;
                                    tmpVtxVecArr[first][eIndex]=1;
                                    Size sNumEdgs1 = vtxVecArr[second].size();
                                    const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
                                    for (Size j = 0; j < sNumEdgs1; ++j) {
                                        nedges_profile++;
                                        if(currentVtxArr1[j]==first)
                                        {
                                            tmpVtxVecArr[second][j]=1;
                                            break;
                                        }
                                    }
                                    M1turn=true;
                                }
                            }
                            done= false;
                            first=second;
                            if(cycle)
                                done= true;
                        }
                    }
                }
                std::cout << "scanned " <<nedges_profile << " edges " <<std::endl;
                if(vtxWght1 >= vtxWght2)
                {
                    *vtxWght =vtxWght1;
                    *card=card1;
                    for (Size s = 0; s < numVtxs; ++s) {
                        mateArr[s] = mateArr1[s];
                    }
                }
                else
                {
                    *vtxWght =vtxWght2;
                    *card=card2;
                    for (Size s = 0; s < numVtxs; ++s) {
                        mateArr[s] = mateArr2[s];
                    }
                }
            }

        template<class ItmQue>
            void MatchingEngine::rPGDPVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {

                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;

                Size numVtxs = graph.mNumVtxs;

                Size midx1=0;
                Size midx2=0;
                Size pidx=0;
                Val W1=0;
                Val W2=0;
                Val temp;
                Size* tp;
                Size pc=0;
                //mates of vertices for M1
                std::vector<Size> M1Vec;
                ResizeVector<Size>(&M1Vec, 3*graph.mNumVtxs);
                Size* M1 = (graph.mNumVtxs == 0) ? NULL : &M1Vec[0];
                if (M1 != NULL) {
                    std::fill(&M1[0], &M1[3*graph.mNumVtxs], 0);
                }
                //mates of vertices for M2
                std::vector<Size> M2Vec;
                ResizeVector<Size>(&M2Vec, 3*graph.mNumVtxs);
                Size* M2 = (graph.mNumVtxs == 0) ? NULL : &M2Vec[0];
                if (M2 != NULL) {
                    std::fill(&M2[0], &M2[3*graph.mNumVtxs], 0);
                }
                //mark visited vertices
                std::vector<Size> markVec;
                ResizeVector<Size>(&markVec, graph.mNumVtxs);
                Size* mark = (graph.mNumVtxs == 0) ? NULL : &markVec[0];
                if (mark != NULL) {
                    std::fill(&mark[0], &mark[graph.mNumVtxs], 0);
                }

                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

                std::vector<std::vector<Size> > tmpVtxVecVec;
                ResizeVector<std::vector<Size> >(&tmpVtxVecVec, numVtxs);
                std::vector<Size>* tmpVtxVecArr = (numVtxs == 0) ? NULL : &tmpVtxVecVec[0];

                for (Size first = 0; first < numVtxs; ++first){
                    Size sNumEdgs = vtxVecArr[first].size();
                    for (Size i = 0; i < sNumEdgs; ++i)
                    {
                        tmpVtxVecArr[first].push_back(0);
                    }

                }

                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

                bool cycle =true;
                //edge number in apth
                Size eNumPath=0;
                //path ength
                Size lp=0;
                Size nedges_profile=0;

                for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
                    Size first = vfirst;

                    bool done = false;
                    midx1=0;
                    midx2=0;
                    W1=0;
                    W2=0;
                    eNumPath=0;
                    pc++;
                    cycle=false;
                    mark[first]=pc;
                    lp=0;
                    while(!done)
                    {
                        if(mateArr[first]!=cNullItm)
                            break;

                        Size sNumEdgs = vtxVecArr[first].size();
                        Size second=cNullItm;
                        Val heaviest = 0.0;
                        Size eIndex=cNullItm;
                        const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[first][0];
                        for (Size i = 0; i < sNumEdgs; ++i) {
                            nedges_profile++;
                            Size tempSecond = currentVtxArr[i];

                            assert(first != tempSecond);
                            // check if the vertex is already visited
                            if(mark[tempSecond]==pc)
                            {
                                if(second==vfirst && lp %2 ==1)
                                    cycle=true;
                                else
                                    continue;
                            }
                            if(mateArr[tempSecond]!=cNullItm)
                                continue;
                            //the vertex belong to different path
                            if(mark[tempSecond]<pc && mark[tempSecond]>0)
                                continue;

                            if (vtxWghtArr[tempSecond]> heaviest && tmpVtxVecArr[first][i] ==0 ) {
                                second= tempSecond;
                                heaviest = vtxWghtArr[tempSecond];
                                eIndex=i;
                            }
                        }

                        done =true;
                        if(heaviest >0.0)
                        {
                            mark[second]=pc;
                            //increment the path length
                            eNumPath++;
                            if(eNumPath==1)
                            {
                                // the first edge to be included in a path
                                W1=heaviest+vtxWghtArr[first];
                                M1[midx1]=first;
                                M1[midx1+1]=second;
                                M1[midx1+2]=heaviest;
                                midx1+=3;
                                tmpVtxVecArr[first][eIndex]=1;
                                Size sNumEdgs1 = vtxVecArr[second].size();
                                const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
                                for (Size j = 0; j < sNumEdgs1; ++j) {
                                    nedges_profile++;
                                    if(currentVtxArr1[j]==first)
                                    {
                                        tmpVtxVecArr[second][j]=1;
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                // using dynamic programming to select the heaviest matching edges in a path
                                if(W2 + heaviest+vtxWghtArr[first] > W1)
                                {
                                    temp=W2;
                                    W2=W1;
                                    W1=temp+heaviest+vtxWghtArr[first];

                                    M2[midx2]=first;
                                    M2[midx2+1]=second;
                                    M2[midx2+2]=heaviest+vtxWghtArr[first];
                                    midx2+=3;

                                    tp=M2;
                                    M2=M1;
                                    M1=tp;

                                    pidx=midx2;
                                    midx2=midx1;
                                    midx1=pidx;

                                }
                                else
                                {
                                    W2=W1;
                                    for(Size i=0;i<midx1;i++)
                                        M2[i]=M1[i];
                                    midx2=midx1;
                                }

                                tmpVtxVecArr[first][eIndex]=1;
                                Size sNumEdgs1 = vtxVecArr[second].size();
                                const Size* currentVtxArr1 = (sNumEdgs == 0) ? NULL : &vtxVecArr[second][0];
                                for (Size j = 0; j < sNumEdgs1; ++j) {
                                    nedges_profile++;
                                    if(currentVtxArr1[j]==first)
                                    {
                                        tmpVtxVecArr[second][j]=1;
                                        break;
                                    }
                                }
                            }
                            eNumPath++;
                            done= false;
                            first=second;
                            if(cycle)
                                done=true;
                        }
                    }

                    //matched vertices found by using dynamic programming
                    Size u=0,v=0;
                    if(midx1>0)
                    {
                        *vtxWght+=W1;
                        for(Size j=0;j< midx1;j+=3)
                        {
                            u=M1[j];
                            v=M1[j+1];
                            mateArr[u]=v;
                            mateArr[v]=u;
                            *card=*card+1;
                        }
                    }
                }
                std::cout << "checked " <<nedges_profile << " edges " <<std::endl;
            }

        template<class ItmQue>
            void MatchingEngine::rLocallyDominantVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;

                Size numVtxs = graph.mNumVtxs;
                // live vertices are considered (they are not matched)
                std::vector<Size> aliveVec;
                ResizeVector<Size>(&aliveVec, graph.mNumVtxs);
                Size* alive = (graph.mNumVtxs == 0) ? NULL : &aliveVec[0];
                if (alive != NULL) {
                    std::fill(&alive[0], &alive[graph.mNumVtxs], 1);
                }

                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
                bool done = false;
                Size nedges_profile=0;

                while(!done)
                {

                    for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){

                        Val heaviest=0.0;
                        //Size partner=-1,id,pidx;
                        //skip if the vertex is already matched
                        if(alive[ vfirst]==0 || mateArr[vfirst]!=cNullItm)
                            continue;

                        Size sNumEdgs = vtxVecArr[vfirst].size();
                        Size second = cNullItm;
                        const Size* currentVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                        for (Size i = 0; i < sNumEdgs; ++i) {
                            nedges_profile++;
                            Size tempSecond = currentVtxArr[i];

                            assert(vfirst != tempSecond);
                            if(alive[tempSecond]==0)
                                continue;
                            if (vtxWghtArr[tempSecond] > heaviest ||(vtxWghtArr[tempSecond]==heaviest && tempSecond>second)) {
                                second= tempSecond;
                                heaviest = vtxWghtArr[tempSecond];
                            }
                        }

                        done =true;
                        if(heaviest >0.0)
                        {
                            mateArr[vfirst]=second;
                        }
                        else
                        {
                            alive[vfirst]=0;
                        }

                    }

                    for (Size vfirst = 0; vfirst < numVtxs; ++vfirst){
                        if(alive[vfirst]==0)
                            continue;

                        Size curr =mateArr[vfirst];
                        if(mateArr[curr]!= vfirst)
                        {
                            //if the edge is not locally dominant then reset mate to be considered in the next iteration
                            mateArr[vfirst]=cNullItm;
                            done=false;
                        }
                        else
                        {
                            // match an edge
                            *vtxWght+=(vtxWghtArr[vfirst]+vtxWghtArr[curr]);
                            *card=*card+1;
                            mateArr[vfirst]=curr;
                            mateArr[curr]=vfirst;
                            alive[vfirst]=0;
                            alive[curr]=0;
                        }
                    }
                }
                std::cout << "checked " <<nedges_profile << " edges " <<std::endl;
            }

        //////////////////////////new//////////////
        template<class ItmQue, class ItmStk>
            void MatchingEngine::rComputeSpecTwoThirdVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {
                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
                //Size se=0;
                Size numVtxs = graph.mNumVtxs;

                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
                bool change = false;
                // a queue to add matched neighbors
                std::deque<Size> tempQ;
                for (Size vfirst = 0; vfirst < numVtxs; vfirst++) {

                    if(mateArr[vfirst]!=cNullItm)
                    {
                        continue;
                    }
                    Size sNumEdgs = vtxVecArr[vfirst].size();
                    const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                    //find an unmatched neighbor ( an ugmenting path of length 1)
                    for(Size i =0; i<sNumEdgs; i++)
                    {
                        Size t = vtxArr[i];//vtxArr[idxArr[vfirst]];
                        //se++;
                        if( mateArr[t] == cNullItm)
                        {
                            mateArr[vfirst]=t;
                            mateArr[t]=vfirst;
                            tempQ.clear();
                            break;
                            //goto endloop;
                        }
                        else
                        {
                            tempQ.push_back(t);
                        }
                    }
                    while (tempQ.empty() == false) {
                        Size s2 = tempQ.front();
                        Size s = mateArr[s2];
                        tempQ.pop_front();
                        Size numEdgs = vtxVecArr[s].size();
                        const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[s][0];
                        //find an ugmenting path of length 3
                        for (Size i = 0; i < numEdgs; ++i) {
                            Size t = vtxArr[i];
                            //se++;
                            if(t != vfirst && t!= s2)
                            {
                                if(mateArr[t]== cNullItm)
                                {
                                    mateArr[vfirst] = s2;
                                    mateArr[s2] = vfirst;
                                    mateArr[s] = t;
                                    mateArr[t]=s;
                                    tempQ.clear();
                                    break;
                                    //goto endloop;
                                }
                            }
                        }
                    }
                    //endloop:;
                }

                // repeate until there is no an augmenting or increasing path
                do
                {
                    change = false;
                    for (Size vfirst = 0; vfirst < numVtxs; vfirst++)
                    {
                        if(mateArr[vfirst]!=cNullItm)
                        {
                            continue;
                        }
                        Val l =vtxWghtArr[vfirst];
                        Val hea = 0.0;
                        Size n1=cNullItm;
                        Size n2=cNullItm;
                        Size n3=cNullItm;
                        Size n4=cNullItm;
                        Size sNumEdgs = vtxVecArr[vfirst].size();
                        const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                        tempQ.clear();
                        //find an augmenting path to heaviest unmatched vertex of length 1 or lightest matched vertex at distance 2
                        for(Size i =0; i<sNumEdgs; i++)
                        {
                            Size t = vtxArr[i];
                            //se++;
                            Size mateT= mateArr[t];
                            if(mateT == cNullItm)
                            {
                                if(vtxWghtArr[t] > hea)
                                {
                                    hea = vtxWghtArr[t];
                                    n1=t;
                                    n2=mateT;
                                    n3=cNullItm;
                                }
                            }
                            else
                            {
                                if(hea<=0.0)
                                {
                                    if(vtxWghtArr[mateT] < l)
                                    {
                                        l=vtxWghtArr[mateT];
                                        n1=t;
                                        n2=mateT;
                                        n3=cNullItm;
                                    }

                                }
                                tempQ.push_back(t);
                            }
                        }
                        while (tempQ.empty() == false) {
                            Size s2 = tempQ.front();
                            tempQ.pop_front();
                            Size s = mateArr[s2];
                            Size numEdgs = vtxVecArr[s].size();
                            const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[s][0];
                            //find an augmenting path to heaviest unmatched vertex of length 3 or lightest matched vertex at distance 4
                            for (Size i = 0; i < numEdgs; ++i) {
                                //se++;
                                Size t = vtxArr[i];
                                if(t != vfirst && t != s2)
                                {
                                    Size mateT = mateArr[t];
                                    if(mateT == cNullItm)
                                    {
                                        if(vtxWghtArr[t] > hea)
                                        {
                                            hea = vtxWghtArr[t];
                                            n1=s2;
                                            n2=s;
                                            n3=t;
                                        }
                                    }
                                    else if(hea<=0.0)
                                    {
                                        if(vtxWghtArr[mateT] < l)
                                        {
                                            l= vtxWghtArr[mateT];
                                            n1=s2;
                                            n2=s;
                                            n3=t;
                                            n4=mateT;
                                        }
                                    }
                                }
                            }
                        }
                        if(hea > 0.0)
                        {	// augment an augmenting path
                            change =true;
                            if(n3==cNullItm)
                            {
                                mateArr[vfirst]=n1;
                                mateArr[n1]=vfirst;
                            }
                            else
                            {
                                mateArr[vfirst]=n1;
                                mateArr[n1] = vfirst;
                                mateArr[n3]=n2;
                                mateArr[n2] = n3;
                            }
                        }
                        else if(n1!=cNullItm)
                        {
                            //reverse an increasing path
                            change = true;
                            if(n3==cNullItm)
                            {
                                mateArr[vfirst]=n1;
                                mateArr[n1]=vfirst;
                                mateArr[n2] = cNullItm;
                            }
                            else
                            {
                                mateArr[vfirst]=n1;
                                mateArr[n1] = vfirst;
                                mateArr[n3]=n2;
                                mateArr[n2] = n3;
                                mateArr[n4] = cNullItm;
                            }
                        }
                    }
                }while(change);

                //timer.Stop();
                //timer.Print();
                //std::cout <<se<<" ";
                *vtxWght = 0.0;
                *card=0;
                for (Size s = 0; s < numVtxs; ++s) {
                    if (mateArr[s] != cNullItm) {
                        *vtxWght += vtxWghtArr[s];
                        ++(*card);
                    }
                }
                (*card)=(*card)/2;
            }
        /////////
        template<class ItmQue, class ItmStk>
            void MatchingEngine::rComputeSpecHalfVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {

                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }
                *card = 0;
                Size numVtxs = graph.mNumVtxs;
                const std::vector<Size>*
                    vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
                //Size se=0;
                bool change = false;
                // find an unmatched neighbor
                for (Size vfirst = 0; vfirst < numVtxs; vfirst++) {
                    if(mateArr[vfirst]!=cNullItm)
                    {
                        continue;
                    }
                    Size sNumEdgs = vtxVecArr[vfirst].size();
                    const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                    for(Size i =0; i<sNumEdgs; i++)
                    {
                        Size t = vtxArr[i];
                        //se++;
                        if( mateArr[t] == cNullItm)
                        {
                            mateArr[vfirst]=t;
                            mateArr[t]=vfirst;
                            break;
                        }
                    }
                }
                // repeat until there is no an augmenting path or an increasing path
                do
                {
                    change = false;
                    for (Size vfirst = 0; vfirst < numVtxs; vfirst++)
                    {
                        if(mateArr[vfirst]!=cNullItm)
                        {
                            continue;
                        }
                        Val l =vtxWghtArr[vfirst];
                        Val hea = 0.0;
                        Size n1=cNullItm;
                        Size n2=cNullItm;
                        Size sNumEdgs = vtxVecArr[vfirst].size();
                        const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                        //find a heaviest unmatched meighbor or lightest matched neighbor
                        for(Size i =0; i<sNumEdgs; i++)
                        {
                            Size t = vtxArr[i];//vtxArr[idxArr[vfirst]];
                            //se++;
                            Size mateT= mateArr[t];
                            if(mateT == cNullItm)
                            {
                                if(vtxWghtArr[t] > hea)
                                {
                                    hea = vtxWghtArr[t];
                                    n1=t;
                                    n2=mateT;
                                }
                            }
                            else
                            {
                                if(hea<=0.0)
                                {
                                    if(vtxWghtArr[mateT] < l)
                                    {
                                        l=vtxWghtArr[mateT];
                                        n1=t;
                                        n2=mateT;
                                    }
                                }
                            }
                        }
                        if(hea > 0.0)
                        {
                            change =true;
                            mateArr[vfirst]=n1;
                            mateArr[n1]=vfirst;
                        }
                        else if(n1!=cNullItm)
                        {
                            change =true;
                            mateArr[vfirst]=n1;
                            mateArr[n1]=vfirst;
                            mateArr[n2] = cNullItm;
                        }
                    }
                }while(change);

                //timer.Stop();
                //timer.Print();
                //std::cout <<se<<" ";
                *vtxWght = 0.0;
                *card=0;
                for (Size s = 0; s < numVtxs; ++s) {
                    if (mateArr[s] != cNullItm) {
                        *vtxWght += vtxWghtArr[s];
                        ++(*card);
                    }
                }
                (*card)=(*card)/2;
            }

        /////////

        template<class ItmQue, class ItmStk>
            void MatchingEngine::rComputeParSpecTwoThirdVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {
                Size numVtxss = graph.mNumVtxs;
                volatile int nlocks[numVtxss];
                int change = 0;
                const std::vector<Size>*
                    vtxVecArr = (numVtxss == 0) ? NULL : &graph.mVtxVecVec[0];
                const Val* vtxWghtArr = (numVtxss == 0) ? NULL : &graph.mVtxWghtVec[0];
                Val tempvtxWght = 0.0;
                Size tempcard=0;
#pragma omp parallel
                {

                    Size numVtxs = numVtxss;
                    // set locks for each vertex in parallel
#pragma omp for schedule(static,256)
                    for(Size i=0;i<numVtxs;i++)
                    {
                        nlocks[i]=0;
                        mateArr[i]=cNullItm;
                    }
                    // a queue to save matched neighbor
                    std::deque<Size> tempQ;
                    //timeCard_start = omp_get_wtime();
                    // in parallel find 2/3 cardinality matching
#pragma omp for schedule(static,256)
                    for (Size vfirst = 0; vfirst < numVtxs; vfirst++) {
                        //startloop:
                        int done =0;

                        while(!done)
                        {

                            Size n1=cNullItm;
                            Size n2=cNullItm;
                            Size n3=cNullItm;

                            if(mateArr[vfirst]!=cNullItm)
                            {
                                done =1;
                                continue;
                            }
                            Size sNumEdgs = vtxVecArr[vfirst].size();
                            const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                            for(Size i =0; i<sNumEdgs; i++)
                            {
                                Size t = vtxArr[i];
                                if( mateArr[t] == cNullItm)
                                {
                                    if(t > vfirst)
                                    {
                                        //done =1;
                                        n1=t;
                                        n2=cNullItm;
                                        n3=cNullItm;
                                        tempQ.clear();
                                        break;
                                        //goto endloop;
                                    }
                                }
                                else
                                {
                                    tempQ.push_back(t);
                                }
                            }
                            while (tempQ.empty() == false) {
                                Size s2 = tempQ.front();
                                tempQ.pop_front();
                                Size s = mateArr[s2];
                                if(s ==vfirst)
                                {
                                    tempQ.clear();
                                    n1=cNullItm;
                                    break;
                                }
                                Size numEdgs = vtxVecArr[s].size();
                                const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[s][0];
                                for (Size i = 0; i < numEdgs; ++i) {
                                    Size t = vtxArr[i];
                                    if(t != vfirst && t!= s2)
                                    {
                                        if(mateArr[t]== cNullItm && t > vfirst)
                                        {
                                            //done =1;
                                            n1=s2;
                                            n2=s;
                                            n3=t;
                                            tempQ.clear();
                                            break;
                                            //goto endloop;
                                        }
                                    }
                                }

                            }
                            //endloop:;
                            //augment an augmenting path after locking vertices along the path
                            if(n1!=cNullItm)
                            {
                                if(n3==cNullItm)
                                {
                                    if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                    {
                                        done=1;
                                        continue;
                                    }
                                    if(mateArr[vfirst]!=cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                        continue;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    if(mateArr[n1] != cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    mateArr[vfirst]=n1;
                                    mateArr[n1]=vfirst;
                                    __sync_lock_release(&nlocks[n1]);
                                    __sync_lock_release(&nlocks[vfirst]);
                                    done =1;
                                }
                                else
                                {
                                    if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                    {
                                        done=1;
                                        continue;
                                    }
                                    if(mateArr[vfirst]!=cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                        continue;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n3], 1))
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    if(mateArr[n3]!=cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[n3]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                    {
                                        __sync_lock_release(&nlocks[n3]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    if(mateArr[n1]!=n2)
                                    {
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[n3]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n2], 1))
                                    {
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[n3]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop;
                                    }
                                    mateArr[vfirst]=n1;
                                    mateArr[n1] = vfirst;
                                    mateArr[n3]=n2;
                                    mateArr[n2] = n3;
                                    __sync_lock_release(&nlocks[n2]);
                                    __sync_lock_release(&nlocks[n1]);
                                    __sync_lock_release(&nlocks[n3]);
                                    __sync_lock_release(&nlocks[vfirst]);
                                    done =1;
                                }
                            }
                            else
                                done=1;
                        }
                    }

                    //timeCard_end  = omp_get_wtime() - timeCard_start;

                    //std::cout  <<" card time "<<timeCard_end<<" ";

                    //timeVWM_start = omp_get_wtime();
                    int tempChange = 0;
                    //repeat until there are no augmenting paths or increasing paths
                    do
                    {
#pragma omp barrier
#pragma omp single nowait
                        {
                            change = 0;
                        }
                        // in parallel find an augmenting path of length at most 3 that reach a heaviest vertex
                        // or an increasing path of length at most 4 that reach lightest vertex
#pragma omp for schedule(static,256)
                        for (Size vfirst = 0; vfirst < numVtxs; vfirst++)
                        {
                            //startloop2:
                            int done =0;
                            while(!done)
                            {

                                if(mateArr[vfirst]!=cNullItm)
                                {
                                    done =1;
                                    continue;
                                }
                                Val l =vtxWghtArr[vfirst];
                                Val hea = 0.0;
                                Size n1=cNullItm;
                                Size n2=cNullItm;
                                Size n3=cNullItm;
                                Size n4=cNullItm;
                                Size sNumEdgs = vtxVecArr[vfirst].size();
                                const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                                tempQ.clear();
                                for(Size i =0; i<sNumEdgs; i++)
                                {
                                    Size t = vtxArr[i];
                                    Size mateT= mateArr[t];
                                    if(mateT == cNullItm)
                                    {
                                        Val Tweigh = vtxWghtArr[t];
                                        if(Tweigh > hea && t > vfirst)
                                        {
                                            hea = Tweigh;
                                            n1=t;
                                            n2=mateT;
                                            n3=cNullItm;
                                        }
                                    }
                                    else
                                    {
                                        if(hea<=0.0)
                                        {
                                            Val Tmweigh = vtxWghtArr[mateT];
                                            if(Tmweigh < l)
                                            {
                                                l=Tmweigh;
                                                n1=t;
                                                n2=mateT;
                                                n3=cNullItm;
                                            }
                                        }

                                        tempQ.push_back(t);
                                    }
                                }
                                while (tempQ.empty() == false) {
                                    Size s2 = tempQ.front();
                                    tempQ.pop_front();
                                    Size s = mateArr[s2];
                                    if(s == cNullItm)
                                    {
                                        Val s2Weight = vtxWghtArr[s2];
                                        if( s2Weight > hea)
                                        {
                                            hea = s2Weight;
                                            n1=s2;
                                            n2=cNullItm;
                                            n3=cNullItm;
                                        }
                                        continue;
                                    }
                                    else if(s == vfirst)
                                    {
                                        n1=cNullItm;
                                        hea=0.0;
                                        break;

                                        //goto startloop2;
                                    }
                                    Size numEdgs = vtxVecArr[s].size();
                                    const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[s][0];
                                    for (Size i = 0; i < numEdgs; ++i) {
                                        Size t = vtxArr[i];
                                        if(t != vfirst && t != s2)
                                        {
                                            Size mateT = mateArr[t];
                                            if(mateT == cNullItm)
                                            {
                                                Val tWeight = vtxWghtArr[t];
                                                if( tWeight > hea && t > vfirst)
                                                {
                                                    hea = tWeight;
                                                    n1=s2;
                                                    n2=s;
                                                    n3=t;
                                                }
                                            }
                                            else if(hea<=0.0)
                                            {
                                                Val tmWeight = vtxWghtArr[mateT];
                                                if(tmWeight < l)
                                                {
                                                    l= tmWeight;
                                                    n1=s2;
                                                    n2=s;
                                                    n3=t;
                                                    n4=mateT;
                                                }
                                            }
                                        }
                                    }
                                }
                                // augment an augmenting path or reverse an increasing path after locking vertices along a path
                                if(hea > 0.0)
                                {
                                    if(tempChange == 0)
                                    {
                                        if(change == 0)
                                            change = 1;
                                        tempChange = 1;
                                    }
                                    if(n3==cNullItm)
                                    {
                                        if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                        {
                                            done =1;
                                            continue;
                                        }
                                        if(mateArr[vfirst]!=cNullItm)
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            done=1;
                                            continue;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(mateArr[n1] != cNullItm)
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        mateArr[vfirst]=n1;
                                        mateArr[n1]=vfirst;
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done =1;
                                    }
                                    else
                                    {
                                        if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                        {
                                            done=1;
                                            continue;
                                        }
                                        if(mateArr[vfirst]!=cNullItm)
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            done=1;
                                            continue;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n3], 1))
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(mateArr[n3]!=cNullItm)
                                        {
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                        {
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(mateArr[n1]!=n2)
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n2], 1))
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        mateArr[vfirst]=n1;
                                        mateArr[n1] = vfirst;
                                        mateArr[n3]=n2;
                                        mateArr[n2] = n3;
                                        __sync_lock_release(&nlocks[n2]);
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[n3]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                    }
                                }
                                else if(n1!=cNullItm)
                                {
                                    if(tempChange == 0)
                                    {
                                        if(change == 0)
                                            change = 1;
                                        tempChange = 1;
                                    }
                                    if(n3==cNullItm)
                                    {
                                        if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                        {
                                            done=1;
                                            continue;
                                        }
                                        if(mateArr[vfirst]!=cNullItm)
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            done=1;
                                            continue;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(mateArr[n1] != n2)
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n2], 1))
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            //goto startloop2;
                                            continue;
                                        }
                                        mateArr[vfirst]=n1;
                                        mateArr[n1]=vfirst;
                                        mateArr[n2] = cNullItm;
                                        __sync_lock_release(&nlocks[n2]);
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                    }
                                    else
                                    {
                                        if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                        {
                                            done=1;
                                            continue;
                                        }
                                        if(mateArr[vfirst]!=cNullItm)
                                        {
                                            done=1;
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n3], 1))
                                        {
                                            __sync_lock_release(&nlocks[vfirst]);
                                            //goto startloop2;
                                            continue;
                                        }
                                        if(mateArr[n3]!=n4)
                                        {
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            //goto startloop2;
                                            continue;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n4], 1))
                                        {
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                        {
                                            __sync_lock_release(&nlocks[n4]);
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            continue;
                                            //goto startloop2;
                                        }
                                        if(mateArr[n1]!=n2)
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[n4]);
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            //goto startloop2;
                                            continue;
                                        }
                                        if(__sync_lock_test_and_set(&nlocks[n2], 1))
                                        {
                                            __sync_lock_release(&nlocks[n1]);
                                            __sync_lock_release(&nlocks[n4]);
                                            __sync_lock_release(&nlocks[n3]);
                                            __sync_lock_release(&nlocks[vfirst]);
                                            //goto startloop2;
                                            continue;
                                        }
                                        mateArr[vfirst]=n1;
                                        mateArr[n1] = vfirst;
                                        mateArr[n3]=n2;
                                        mateArr[n2] = n3;
                                        mateArr[n4] = cNullItm;
                                        __sync_lock_release(&nlocks[n2]);
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[n4]);
                                        __sync_lock_release(&nlocks[n3]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                    }
                                }
                                else
                                    done=1;

                            }
                        }

                    }while(change);

#pragma omp for schedule(static,256) nowait reduction(+:tempvtxWght,tempcard)
                    for (Size s = 0; s < numVtxs; ++s)
                    {
                        if (mateArr[s] != cNullItm) {
                            tempvtxWght += vtxWghtArr[s];
                            ++(tempcard);
                        }
                    }
                }
                *vtxWght=tempvtxWght;
                *card=(tempcard)/2;
            }
        /////////////////////
        template<class ItmQue, class ItmStk>
            void MatchingEngine::rComputeParSpecHalfVtxWght(const Graph& graph,
                    Size* mateArr, Size* card, Val* vtxWght) const {
                Size numVtxss = graph.mNumVtxs;
                int nlocks[numVtxss];
                int change = 0;
                const std::vector<Size>*
                    vtxVecArr = (numVtxss == 0) ? NULL : &graph.mVtxVecVec[0];
                const Val* vtxWghtArr = (numVtxss == 0) ? NULL : &graph.mVtxWghtVec[0];
                Val tempvtxWght = 0.0;
                Size tempcard=0;
#pragma omp parallel
                {
                    Size numVtxs = numVtxss;
                    // set locks for each vertex in parallel
#pragma omp for schedule(static,256)
                    for(Size i=0;i<numVtxs;i++)
                    {
                        nlocks[i]=0;
                        mateArr[i]=cNullItm;
                    }
                    // in parallel find 2/3 cardinality matching
#pragma omp for schedule(static,256) /*firstprivate(tempQ) shared(st)*/
                    for (Size vfirst = 0; vfirst < numVtxs; vfirst++) {
                        //startloop:
                        int done =0;

                        while(!done)
                        {
                            Size n1=cNullItm;
                            if(mateArr[vfirst]!=cNullItm)
                            {
                                done=1;
                                continue;
                            }
                            Size sNumEdgs = vtxVecArr[vfirst].size();
                            const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                            for(Size i =0; i<sNumEdgs; i++)
                            {
                                Size t = vtxArr[i];
                                if( mateArr[t] == cNullItm)
                                {
                                    if(t > vfirst)
                                    {
                                        //done =1;
                                        n1=t;
                                        break;
                                        //goto endloop;
                                    }
                                }
                            }
                            //endloop:;
                            //augment an augmenting path after locking vertices along the path
                            if(n1!=cNullItm)
                            {
                                if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                {
                                    //done=1;
                                    continue;
                                }
                                if(mateArr[vfirst]!=cNullItm)
                                {
                                    __sync_lock_release(&nlocks[vfirst]);
                                    done=1;
                                    continue;
                                }
                                if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                {
                                    __sync_lock_release(&nlocks[vfirst]);
                                    continue;
                                    //goto startloop;
                                }
                                if(mateArr[n1] != cNullItm)
                                {
                                    __sync_lock_release(&nlocks[n1]);
                                    __sync_lock_release(&nlocks[vfirst]);
                                    continue;
                                    //goto startloop;
                                }
                                mateArr[vfirst]=n1;
                                mateArr[n1]=vfirst;
                                __sync_lock_release(&nlocks[n1]);
                                __sync_lock_release(&nlocks[vfirst]);
                                done=1;
                            }
                            else
                                done=1;
                        }
                    }
                    int tempChange = 0;
                    //repeat until there are no augmenting paths or increasing paths
                    do
                    {
#pragma omp barrier
#pragma omp single nowait
                        {
                            change = 0;
                        }
                        // in parallel find an augmenting path of length 1 that reach a heaviest vertex
                        // or an increasing path of length 2 that reach lightest vertex
#pragma omp for schedule(static,256)
                        for (Size vfirst = 0; vfirst < numVtxs; vfirst++)
                        {

                            int done =0;
                            //startloop2:
                            while(!done)
                            {

                                if(mateArr[vfirst]!=cNullItm)
                                {
                                    done=1;
                                    continue;
                                }
                                Val l =vtxWghtArr[vfirst];
                                Val hea = 0.0;
                                Size n1=cNullItm;
                                Size n2=cNullItm;
                                Size sNumEdgs = vtxVecArr[vfirst].size();
                                const Size* vtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[vfirst][0];
                                for(Size i =0; i<sNumEdgs; i++)
                                {
                                    Size t = vtxArr[i];
                                    Size mateT= mateArr[t];
                                    if(mateT == cNullItm)
                                    {
                                        Val Tweigh = vtxWghtArr[t];
                                        if(Tweigh > hea && t > vfirst)
                                        {
                                            hea = Tweigh;
                                            n1=t;
                                            n2=mateT;
                                        }
                                    }
                                    else
                                    {
                                        if(hea<=0.0)
                                        {
                                            Val Tmweigh = vtxWghtArr[mateT];
                                            if(Tmweigh < l)
                                            {
                                                l=Tmweigh;
                                                n1=t;
                                                n2=mateT;
                                            }
                                        }
                                    }
                                }
                                // augment an augmenting path or reverse an increasing path after locking vertices along a path
                                if(hea > 0.0)
                                {
                                    if(tempChange == 0)
                                    {
                                        if(change == 0)
                                            change = 1;
                                        tempChange = 1;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                    {
                                        //done=1;
                                        continue;
                                    }
                                    if(mateArr[vfirst]!=cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                        continue;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop2;
                                    }
                                    if(mateArr[n1] != cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop2;
                                    }
                                    mateArr[vfirst]=n1;
                                    mateArr[n1]=vfirst;
                                    __sync_lock_release(&nlocks[n1]);
                                    __sync_lock_release(&nlocks[vfirst]);
                                    done=1;
                                }
                                else if(n1!=cNullItm)
                                {
                                    if(tempChange == 0)
                                    {
                                        if(change == 0)
                                            change = 1;
                                        tempChange = 1;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[vfirst], 1))
                                    {
                                        //done=1;
                                        continue;
                                    }
                                    if(mateArr[vfirst]!=cNullItm)
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        done=1;
                                        continue;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n1], 1))
                                    {
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop2;
                                    }
                                    if(mateArr[n1] != n2)
                                    {
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop2;
                                    }
                                    if(__sync_lock_test_and_set(&nlocks[n2], 1))
                                    {
                                        __sync_lock_release(&nlocks[n1]);
                                        __sync_lock_release(&nlocks[vfirst]);
                                        continue;
                                        //goto startloop2;
                                    }
                                    mateArr[vfirst]=n1;
                                    mateArr[n1]=vfirst;
                                    mateArr[n2] = cNullItm;
                                    __sync_lock_release(&nlocks[n2]);
                                    __sync_lock_release(&nlocks[n1]);
                                    __sync_lock_release(&nlocks[vfirst]);
                                    done=1;
                                }
                                else
                                    done=1;
                            }
                        }
                    }while(change);

#pragma omp for schedule(static,256) nowait reduction(+:tempvtxWght,tempcard)
                    for (Size s = 0; s < numVtxs; ++s)
                    {
                        if (mateArr[s] != cNullItm) {
                            tempvtxWght += vtxWghtArr[s];
                            ++(tempcard);
                        }
                    }
                }

                *vtxWght=tempvtxWght;
                *card=(tempcard)/2;
            }


        ///////////////////////////////
        template<class ItmQue>
            void MatchingEngine::rnComponents(const Graph& graph) const
            {
                const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

                ItmQue prcsbQue(graph.mNumVtxs);
                ItmQue dQue(graph.mNumVtxs);

                std::vector<Size> visitVec;
                ResizeVector<Size>(&visitVec, graph.mNumVtxs);
                Size* visit = (graph.mNumVtxs == 0) ? NULL : &visitVec[0];
                if (visit != NULL) {
                    std::fill(&visit[0], &visit[graph.mNumVtxs], cZeroVal);
                }

                // const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

                Size verCount=0;
                Size compCount=0;
                //Size maxComponentSize=0;
                Size verOfMaxCom=0;
                Size edgOfMaxCom=0;
                Size tmpverOfMaxCom=0;
                Size tmpedgOfMaxCom=0;

                Size degreeSum=0;
                Size nonIsolVer=0;
                Size maxD=0;
                Size minD=cNullItm;
                //Size countMin=0;
                Size maxDepth=0;
                for (Size i=0;i<graph.mNumVtxs;i++)
                {
                    Size sNumEdgs = vtxVecArr[i].size();
                    if(sNumEdgs>0)
                    {
                        if(sNumEdgs>maxD)
                            maxD=sNumEdgs;

                        if(sNumEdgs< minD)
                            minD= sNumEdgs;
                    }
                }
                Size mecom=0,mvcom=0;
                for (Size i=0;i<graph.mNumVtxs;i++)
                {
                    if(visit[i]==0)
                    {
                        compCount++;
                        prcsbQue.Push(i);
                        dQue.Push(0);
                        while (prcsbQue.Empty() == false)
                        {
                            Size s = prcsbQue.Front();
                            Size d = dQue.Front();
                            if(d> maxDepth)
                            {
                                maxDepth =d;
                            }
                            visit[s]=1;
                            prcsbQue.Pop();
                            dQue.Pop();
                            verCount++;
                            tmpverOfMaxCom++;
                            Size sNumEdgs = vtxVecArr[s].size();
                            tmpedgOfMaxCom+=sNumEdgs;
                            if(sNumEdgs>0)
                            {
                                nonIsolVer++;
                                degreeSum+=sNumEdgs;
                            }

                            const Size* sVtxArr = (sNumEdgs == 0) ? NULL : &vtxVecArr[s][0];
                            for (Size j = 0; j < sNumEdgs; ++j) {
                                Size t = sVtxArr[j];
                                assert(s != t);
                                if(visit[t]==0)
                                {

                                    prcsbQue.Push(t);
                                    d++;
                                    dQue.Push(d);
                                    visit[t]=1;
                                }

                            }
                        }
                    }
                    if(tmpverOfMaxCom>verOfMaxCom)
                    {
                        verOfMaxCom=tmpverOfMaxCom;
                        mvcom=compCount;
                    }
                    if(tmpedgOfMaxCom>edgOfMaxCom)
                    {
                        edgOfMaxCom=tmpedgOfMaxCom;
                        mecom=compCount;
                    }

                    if(verCount >= graph.mNumVtxs)
                        break;

                    tmpverOfMaxCom=0;
                    tmpedgOfMaxCom=0;
                }
                Val dMean = degreeSum/(double)nonIsolVer;
                Val sumSD = 0.0;
                for (Size i=0; i<graph.mNumVtxs; i++)
                {
                    Size sNumEdgs = vtxVecArr[i].size();
                    sumSD+=std::pow(dMean-sNumEdgs,2.0);
                }
                Val sdev = std::sqrt(sumSD/nonIsolVer);
                std::cout << "Number of vertices: " <<graph.mNumVtxs << std::endl;
                std::cout << "Number of edges: " <<graph.mNumEdgs << std::endl;
                std::cout << "Number of components: " <<compCount << std::endl;
                std::cout << "BFS depth: " <<maxDepth << std::endl;
                std::cout << "Largest number of vertices: " <<verOfMaxCom <<" in component "<< mvcom<< std::endl;
                std::cout << "Largest number of edges: " <<edgOfMaxCom/2 <<" in component "<< mecom<< std::endl;
                std::cout << "average degrees: " <<degreeSum/(double)nonIsolVer << std::endl;
                std::cout << "standard dev: " <<degreeSum/(double)nonIsolVer << std::endl;
                std::cout << "number of non isolated vertices: " <<nonIsolVer << std::endl;
                std::cout << "max degree: " <<maxD << " min degree: " <<minD << std::endl;
                std::cout << dMean <<" " << sdev <<std::endl;

            }

        template<class ItmQue>
            void MatchingEngine::rComputeHalfVtxWghtMatching(
                    const Graph& graph, Size* mateArr, Size* card,
                    Val* vtxWght) const {

                Size numVtxs =graph.mNumVtxs;

                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }

                *card = 0;
                const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];
                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];
                std::vector<std::pair<Size, Val> > sExpsdLst2;
                ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
                std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0];
                //Timer timer;

                for (Size first = 0; first < numVtxs; ++first) {
                    sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]);
                }
                //timer.Start();
                //sort vertices in non-increasing order of their weights
                std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater<std::pair<Size, Val> >());
                //timer.Stop();
                //timer.Print();
                for (Size current = 0; current < numVtxs; ++current) {
                    Size sFirst =  sExpsdLst2[current].first;
                    if(mateArr[sFirst]!=cNullItm)
                        continue;
                    Val heaviest = 0.0;
                    Size mateT = cNullItm;
                    Size numEdgs = vtxVecArr[sFirst].size();
                    const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
                    //find  heaviest neighbor
                    for (Size i = 0; i < numEdgs; ++i) {
                        Size t = vtxArr[i];
                        if(mateArr[t]== cNullItm && vtxWghtArr[t] > heaviest)
                        {
                            heaviest = 	vtxWghtArr[t];
                            mateT=t;
                        }
                    }
                    //label:
                    if(heaviest > 0.0)
                    {
                        //match an edge
                        mateArr[mateT]=sFirst;
                        mateArr[sFirst]=mateT;
                        ++(*card);
                    }

                }
                *vtxWght = 0.0;
                for (Size s = 0; s < numVtxs; ++s) {
                    if (mateArr[s] != cNullItm) {
                        *vtxWght += vtxWghtArr[s];
                    }
                }
            }


        template<class ItmQue>
            void MatchingEngine::rComputeTwoThirdVtxWghtMatching(
                    const Graph& graph, Size* mateArr, Size* card,
                    Val* vtxWght) const {
                Size numVtxs =graph.mNumVtxs;

                if (mateArr != NULL) {
                    std::fill(&mateArr[0], &mateArr[graph.mNumVtxs], cNullItm);
                }


                *card = 0;

                const std::vector<Size>* vtxVecArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxVecVec[0];

                const Val* vtxWghtArr = (graph.mNumVtxs == 0) ? NULL : &graph.mVtxWghtVec[0];

                /* std::vector<Size> deadVec;
                   ResizeVector<Size>(&deadVec, graph.mNumVtxs);
                   Size* deadArr = (graph.mNumVtxs == 0) ? NULL : &deadVec[0];
                   if (deadArr != NULL) {
                   std::fill(&deadArr[0], &deadArr[graph.mNumVtxs], 0);
                   }*/
                // a queue to save matched neighbors
                std::deque<Size> sPrcsbQue;
                //int g=0;
                /*for (Size v = 0; v < numVtxs; ++v) {
                  Size numEdgs = vtxVecArr[v].size();
                  const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[v][0];
                  for (Size i = 1; i < numEdgs; ++i) {
                  if(vtxWghtArr[vtxArr[i-1]] < vtxWghtArr[vtxArr[i]])
                  {
                  std::cout << "smaller" << std::endl;
                  g=1;
                  break;
                  }

                  }
                  if(g)
                  break;
                  }*/
                //do test of ordered adj list
                std::vector<std::pair<Size, Val> > sExpsdLst2;
                ResizeVector<std::pair<Size, Val> >(&sExpsdLst2, numVtxs);
                std::pair<Size,Val> * sExpsdLst2Arr=  (numVtxs == 0) ? NULL : &sExpsdLst2[0];

                for (Size first = 0; first < numVtxs; ++first) {
                    sExpsdLst2Arr[first]=std::pair<Size, Val>(first, vtxWghtArr[first]);
                }
                //sort vertices in non-increasing order
                std::sort(sExpsdLst2.begin(),sExpsdLst2.end(),ValGreater<std::pair<Size, Val> >());
                Val lap=0;
                Size nap=0;
                Size nver=0;
                // for each vertex find an augmenting path of length at most 3 that reach a heaviest unmatched vertex
                for (Size current = 0; current < numVtxs; ++current) {
                    Size sFirst =  sExpsdLst2Arr[current].first;
                    if(mateArr[sFirst]!=cNullItm)
                        continue;
                    Size numEdgs = vtxVecArr[sFirst].size();
                    const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[sFirst][0];
                    Val heaviest =0;
                    Size mateT =cNullItm;
                    Size tempMate=cNullItm;
                    sPrcsbQue.clear();
                    for (Size i = 0; i < numEdgs; ++i) {
                        nver++;
                        Size t = vtxArr[i];
                        if(mateArr[t]== cNullItm)
                        {
                            if(vtxWghtArr[t] > heaviest)
                            {
                                heaviest = vtxWghtArr[t];
                                tempMate=cNullItm;
                                mateT = t;
                            }
                        }
                        else //if(!deadArr[mateArr[t]])
                        {
                            sPrcsbQue.push_back(mateArr[t]);
                        }
                    }
                    while (sPrcsbQue.empty() == false) {
                        Size s = sPrcsbQue.front();
                        sPrcsbQue.pop_front();
                        Size numEdgs = vtxVecArr[s].size();
                        const Size* vtxArr = (numEdgs == 0) ? NULL : &vtxVecArr[s][0];
                        for (Size i = 0; i < numEdgs; ++i) {
                            nver++;
                            Size t = vtxArr[i];
                            if(mateArr[t] == cNullItm && t != sFirst)
                            {
                                if(vtxWghtArr[t] > heaviest)
                                {
                                    heaviest = vtxWghtArr[t];
                                    mateT = t;
                                    tempMate=mateArr[s];
                                }
                                //break;
                            }
                        }
                        //if(i == numEdgs)
                        //deadArr[s]=1;
                    }
                    if(heaviest > 0)
                    {
                        //augment an augmenting path
                        if(tempMate == cNullItm )
                        {
                            mateArr[mateT]=sFirst;
                            mateArr[sFirst]=mateT;
                            lap=lap+1;
                            nap=nap+1;
                            ++(*card);
                        }
                        else
                        {
                            Size tempS = mateArr[tempMate];
                            mateArr[mateT]=tempS;
                            mateArr[tempS]=mateT;
                            mateArr[tempMate]=sFirst;
                            mateArr[sFirst]=tempMate;
                            lap=lap+3;
                            nap=nap+1;
                            ++(*card);
                        }
                    }

                }
                std::cout << "no of visited vertices " <<nver << " no of aug path " <<nap<<" sum of lengths " <<lap<<" avg lengths " <<lap/nap<<std::endl;
                *vtxWght = 0.0;
                for (Size s = 0; s < numVtxs; ++s) {
                    if (mateArr[s] != cNullItm) {
                        *vtxWght += vtxWghtArr[s];
                    }
                }
            }

        } // namespace Matchbox

#endif // MATCHING_ENGINE_H
