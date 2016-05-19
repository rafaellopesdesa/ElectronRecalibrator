#ifndef GBRTREE_HPP
#define GBRTREE_HPP

#include <vector>
#include <map>

namespace TMVA {
  class DecisionTree;
  class DecisionTreeNode;
}

namespace CORE_GBR {
  class GBRTree;
}

class CORE_GBR::GBRTree {

    public:

       GBRTree();
       explicit GBRTree(const TMVA::DecisionTree *tree, double scale, bool useyesnoleaf, bool adjustboundary);
       virtual ~GBRTree();
       
       double GetResponse(const float* vector) const;
       int TerminalIndex(const float *vector) const;
       
       std::vector<float> &Responses() { return fResponses; }       
       const std::vector<float> &Responses() const { return fResponses; }
       
       std::vector<unsigned char> &CutIndices() { return fCutIndices; }
       const std::vector<unsigned char> &CutIndices() const { return fCutIndices; }
       
       std::vector<float> &CutVals() { return fCutVals; }
       const std::vector<float> &CutVals() const { return fCutVals; }
       
       std::vector<int> &LeftIndices() { return fLeftIndices; }
       const std::vector<int> &LeftIndices() const { return fLeftIndices; } 
       
       std::vector<int> &RightIndices() { return fRightIndices; }
       const std::vector<int> &RightIndices() const { return fRightIndices; }
       

       
    protected:      
        unsigned int CountIntermediateNodes(const TMVA::DecisionTreeNode *node);
        unsigned int CountTerminalNodes(const TMVA::DecisionTreeNode *node);
      
        void AddNode(const TMVA::DecisionTreeNode *node, double scale, bool isregression, bool useyesnoleaf, bool adjustboundary);
        
	std::vector<unsigned char> fCutIndices;
	std::vector<float> fCutVals;
	std::vector<int> fLeftIndices;
	std::vector<int> fRightIndices;
	std::vector<float> fResponses;  
        
  
};


//_______________________________________________________________________
inline double CORE_GBR::GBRTree::GetResponse(const float* vector) const {
  
  int index = 0;
  
  unsigned char cutindex = fCutIndices[0];
  float cutval = fCutVals[0];
  
  while (true) {
     
    if (vector[cutindex] > cutval) {
      index = fRightIndices[index];
    }
    else {
      index = fLeftIndices[index];
    }
    
    if (index>0) {
      cutindex = fCutIndices[index];
      cutval = fCutVals[index];
    }
    else {
      return fResponses[-index];
    }
    
  }
  

}

//_______________________________________________________________________
inline int CORE_GBR::GBRTree::TerminalIndex(const float* vector) const {
  
  int index = 0;
  
  unsigned char cutindex = fCutIndices[0];
  float cutval = fCutVals[0];
  
  while (true) {
    if (vector[cutindex] > cutval) {
      index = fRightIndices[index];
    }
    else {
      index = fLeftIndices[index];
    }
    
    if (index>0) {
      cutindex = fCutIndices[index];
      cutval = fCutVals[index];
    }
    else {
      return (-index);
    }
    
  }
  

}
  
#endif
