#ifndef HPP_ELE_RECALIBRATOR
#define HPP_ELE_RECALIBRATOR

#include "EpCombinationTool.h"
#include "../../CMS3.h"
#include "TH2D.h"

struct resolData { 
  int runMin;
  int runMax;
  double etaMin;
  double etaMax;
  double r9Min;
  double r9Max;
  double rho;
  double err_rho;
  double phi;
  double err_phi;
};

struct scaleData {
  int runMin;
  int runMax;
  double etaMin;
  double etaMax;
  double r9Min;
  double r9Max;
  double EtMin;
  double EtMax;
  double deltaP;
  double err_deltaP;
  double err_deltaP_stat;
  double err_deltaP_syst;
};

class Recalibrator {

 protected:
  vector<resolData> smearing;
  vector<scaleData> scale;
  EpCombinationTool combine;
  
 public:
  Recalibrator(TString resolfile, TString scalefile, TString combfile, TString combtree);
  void Recalibrate(int elsIdx, bool isMC, int syst, float& recalib_energy);
};

#endif 

  

//  ep.init("data/GBRForest_data_25ns.root", "gedelectron_p4combination_25ns");
  
