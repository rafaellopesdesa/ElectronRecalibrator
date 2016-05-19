#include "Recalibrator.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"

#include <fstream>
using namespace std;

Recalibrator::Recalibrator(TString resolfile, TString scalefile, TString combfile, TString combtree) {
  
  int runMin = 0, runMax = 900000;
  TString category, dummy, phi_string, err_phi_string;
  double rho, phi, Emean, err_rho, err_phi, err_Emean;
  double deltaP, err_deltaP, err_deltaP_stat, err_deltaP_syst;
  double r9Min = -1., r9Max = 1.;

  ifstream f_resol(resolfile);
  while(!f_resol.eof()) {
    
    f_resol >> category >> 
      Emean >> err_Emean >>
      rho >> err_rho >> phi_string >> err_phi_string;
    
    if (phi_string.Contains("M_PI_2")) phi = TMath::Pi()/2.;
    else phi = phi_string.Atof();
    
    if (err_phi_string.Contains("M_PI_2")) err_phi = TMath::Pi()/2.;
    else err_phi = err_phi_string.Atof();

    if (category.Contains("-highR9")) {
      category.ReplaceAll("-highR9", "");
      r9Min = 0.94;
      r9Max = 1.00;
    } else if (category.Contains("-lowR9")) {
      category.ReplaceAll("-lowR9", "");
      r9Min = -1.00;
      r9Max = 0.94;
    }
    
    resolData data_resol;
    data_resol.runMin  = 0;
    data_resol.runMax  = 9000000;
    data_resol.rho     = rho;
    data_resol.err_rho = err_rho;
    data_resol.phi     = phi;
    data_resol.err_phi = err_phi;
    data_resol.r9Min   = r9Min;
    data_resol.r9Max   = r9Max;
    data_resol.etaMin  = static_cast<TObjString*>(category.Tokenize("_")->At(1))->String().Atof();
    data_resol.etaMax  = static_cast<TObjString*>(category.Tokenize("_")->At(2))->String().Atof();
    smearing.push_back(data_resol);    
    
  }
  f_resol.close();


  ifstream f_scale(scalefile);
  while (!f_scale.eof()) {
    f_scale >> category >> dummy
	    >> runMin >> runMax
	    >> deltaP >> err_deltaP >> err_deltaP_stat >> err_deltaP_syst;

    bool et_binning = true;
    if (category.Contains("-highR9-")) {
      category.ReplaceAll("-highR9-", "_");
      r9Min = 0.94;
      r9Max = 1.00;
    } else if (category.Contains("-lowR9-")) {
      category.ReplaceAll("-lowR9-", "_");
      r9Min = -1.00;
      r9Max = 0.94;
    } else if (category.Contains("-highR9")) {
      category.ReplaceAll("-highR9", "");
      r9Min = 0.94;
      r9Max = 1.00;
      et_binning = false;
    } else if (category.Contains("-lowR9")) {
      category.ReplaceAll("-lowR9-", "");
      r9Min = -1.00;
      r9Max = 0.94;
      et_binning = false;
    }

    scaleData data_scale;
    data_scale.runMin          = runMin;
    data_scale.runMax          = runMax;
    data_scale.deltaP          = deltaP;
    data_scale.err_deltaP      = err_deltaP;
    data_scale.err_deltaP_stat = err_deltaP_stat;
    data_scale.err_deltaP_syst = err_deltaP_syst;

    data_scale.etaMin          = (static_cast<TObjString*>(category.Tokenize("_")->At(1))->String()).Atof();
    data_scale.etaMax          = (static_cast<TObjString*>(category.Tokenize("_")->At(2))->String()).Atof();
    if (et_binning) {
      data_scale.EtMin           = (static_cast<TObjString*>(category.Tokenize("_")->At(4))->String()).Atof();
      data_scale.EtMax           = (static_cast<TObjString*>(category.Tokenize("_")->At(5))->String()).Atof();
    } else {
      data_scale.EtMin = 0.;
      data_scale.EtMin = 14000.;
    }
    scale.push_back(data_scale);    
  }
  f_scale.close();

  combine.init(combfile.Data(), combtree.Data());

}

void Recalibrator::Recalibrate(int elsIdx, bool isMC, int syst, float& recalib_energy) { 

  double electron_eta = TMath::Abs(tas::els_etaSC().at(elsIdx));
  double electron_r9  = tas::els_r9_full5x5().at(elsIdx);
  double electron_et  = tas::els_ecalEnergy().at(elsIdx)/cosh(electron_eta);
  int runnum = tas::evt_run();

  if (electron_eta >= 1.4442 && electron_eta <= 1.556) {
    recalib_energy = tas::els_p4().at(elsIdx).e();
    return;
  }

  double scale_corr = 1.;
  double scale_corr_err_stat = 0.;
  double scale_corr_err_syst = 0.;
  for (auto scale_data : scale) {
    if ((electron_eta > scale_data.etaMin) && 
	(electron_eta < scale_data.etaMax) &&
        (electron_et > scale_data.EtMin) && 
	(electron_et < scale_data.EtMax) &&
	(electron_r9 > scale_data.r9Min) &&
	(electron_r9 < scale_data.r9Max) &&
	(runnum > scale_data.runMin) &&
	(runnum < scale_data.runMax)) {
      scale_corr = scale_data.deltaP;
      scale_corr_err_stat = scale_data.err_deltaP_stat;
      scale_corr_err_syst = scale_data.err_deltaP_syst;
      scale_corr += syst*sqrt(scale_corr_err_stat*scale_corr_err_stat+scale_corr_err_syst*scale_corr_err_syst);
      break;
    }
  }

  double smearing_corr = 0.;
  for (auto smearing_data : smearing) {
    if ((electron_eta > smearing_data.etaMin) && 
	(electron_eta < smearing_data.etaMax) && 
	(electron_r9 > smearing_data.r9Min) &&
	(electron_r9 < smearing_data.r9Max) &&
	(runnum > smearing_data.runMin) &&
	(runnum < smearing_data.runMax)) {
      
      double rho = smearing_data.rho + smearing_data.err_rho * syst;
      double phi = smearing_data.phi;
      
      double constTerm =  rho*sin(phi);
      double alpha =  rho*cos(phi);

      smearing_corr = sqrt(constTerm * constTerm + alpha * alpha / electron_et);
      break;
    }
  }
	
  double new_ecal_energy, new_ecal_energy_error;
  if (isMC) {
    double corr = gRandom->Gaus(1., smearing_corr);
    new_ecal_energy       = tas::els_ecalEnergy().at(elsIdx) * corr;
    new_ecal_energy_error = std::hypot(tas::els_ecalEnergyError().at(elsIdx) * corr, smearing_corr * new_ecal_energy);
  } else {
    new_ecal_energy       = tas::els_ecalEnergy().at(elsIdx) * scale_corr;
    new_ecal_energy_error = std::hypot(tas::els_ecalEnergyError().at(elsIdx) * scale_corr, smearing_corr * new_ecal_energy);
  }
  combine.combine(elsIdx, new_ecal_energy, new_ecal_energy_error, recalib_energy);
  
}
