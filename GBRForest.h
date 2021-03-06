#ifndef GBRForest_HPP
#define GBRForest_HPP

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBRForest                                                            //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// Designed to be built from TMVA-trained trees, but could also be      //
// generalized to otherwise-trained trees, classification,              //
//  or other boosting methods in the future                             //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include "GBRTree.h"
#include <math.h>
#include <stdio.h>

namespace TMVA {
  class MethodBDT;
}

namespace CORE_GBR {
  class GBRForest;
}

class CORE_GBR::GBRForest {

    public:

       GBRForest();
       explicit GBRForest(const TMVA::MethodBDT *bdt);
       virtual ~GBRForest();
       
       double GetResponse(const float* vector) const;
       double GetGradBoostClassifier(const float* vector) const;
       double GetAdaBoostClassifier(const float* vector) const { return GetResponse(vector); }
       
       //for backwards-compatibility
       double GetClassifier(const float* vector) const { return GetGradBoostClassifier(vector); }
       
       void SetInitialResponse(double response) { fInitialResponse = response; }
       
       std::vector<CORE_GBR::GBRTree> &Trees() { return fTrees; }
       const std::vector<CORE_GBR::GBRTree> &Trees() const { return fTrees; }
       
    protected:
      double               fInitialResponse;
      std::vector<CORE_GBR::GBRTree> fTrees;  
      
  
};

//_______________________________________________________________________
inline double CORE_GBR::GBRForest::GetResponse(const float* vector) const {
  double response = fInitialResponse;
  for (std::vector<CORE_GBR::GBRTree>::const_iterator it=fTrees.begin(); it!=fTrees.end(); ++it) {
    response += it->GetResponse(vector);
  }
  return response;
}

//_______________________________________________________________________
inline double CORE_GBR::GBRForest::GetGradBoostClassifier(const float* vector) const {
  double response = GetResponse(vector);
  return 2.0/(1.0+exp(-2.0*response))-1; //MVA output between -1 and 1
}

#endif
