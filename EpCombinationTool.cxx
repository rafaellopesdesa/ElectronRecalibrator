#include "EpCombinationTool.h"
#include "GBRForest.h"
#include <TFile.h>
#include <TSystem.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;
using namespace CORE_GBR;


/*****************************************************************/
EpCombinationTool::EpCombinationTool():
    m_forest(NULL), m_ownForest(false)
/*****************************************************************/
{
}



/*****************************************************************/
EpCombinationTool::~EpCombinationTool()
/*****************************************************************/
{
    if(m_ownForest) delete m_forest;
}


/*****************************************************************/
bool EpCombinationTool::init(const string& regressionFileName, const string& bdtName)
/*****************************************************************/
{
    TFile* regressionFile = TFile::Open(regressionFileName.c_str());
    if(!regressionFile)
    {
        cout<<"ERROR: Cannot open regression file "<<regressionFileName<<"\n";
        return false;
    }
    if(m_ownForest) delete m_forest;
    m_forest = (GBRForest*) regressionFile->Get(bdtName.c_str());
    m_ownForest = true;
    //regressionFile->GetObject(bdtName.c_str(), m_forest); 
    if(!m_forest)
    {
        cout<<"ERROR: Cannot find forest "<<bdtName<<" in "<<regressionFileName<<"\n";
        regressionFile->Close();
        return false;
    }
    regressionFile->Close();
    return true;
}

bool EpCombinationTool::init(const GBRForest *forest) 
{
    if(m_ownForest) delete m_forest;
    m_forest = forest;
    m_ownForest = false;
    return true;
}


/*****************************************************************/
void EpCombinationTool::combine(int elsIdx, float old_energy, float old_energy_error, float& recalib_energy) const
/*****************************************************************/
{
    if(!m_forest)
    {
        cout<<"ERROR: The combination tool is not initialized\n";
        return;
    }

    float energy = old_energy; // tas::els_ecalEnergy().at(elsIdx);
    float energyError = old_energy_error; // tas::els_ecalEnergyError().at(elsIdx);
    float momentum = tas::els_p4In().at(elsIdx).R();
    float momentumError = tas::els_trackMomentumError().at(elsIdx);
    int  electronClass = tas::els_class().at(elsIdx);
    bool isEcalDriven = tas::els_isEcalDriven().at(elsIdx);
    bool isTrackerDriven =  tas::els_isTrackerDriven().at(elsIdx);
    bool isEB = tas::els_isEB().at(elsIdx);

    // compute relative errors and ratio of errors
    float energyRelError = energyError / energy;
    float momentumRelError = momentumError / momentum;
    float errorRatio = energyRelError / momentumRelError;

    // calculate E/p and corresponding error
    float eOverP = energy / momentum;
    float eOverPerror = sqrt(
            (energyError/momentum)*(energyError/momentum) +
            (energy*momentumError/momentum/momentum)*
            (energy*momentumError/momentum/momentum));

    // fill input variables
    float regressionInputs[11];
    regressionInputs[0]  = static_cast<float>(energy);
    regressionInputs[1]  = static_cast<float>(energyRelError);
    regressionInputs[2]  = static_cast<float>(momentum);
    regressionInputs[3]  = static_cast<float>(momentumRelError);
    regressionInputs[4]  = static_cast<float>(errorRatio);
    regressionInputs[5]  = static_cast<float>(eOverP);
    regressionInputs[6]  = static_cast<float>(eOverPerror);
    regressionInputs[7]  = static_cast<float>(isEcalDriven);
    regressionInputs[8]  = static_cast<float>(isTrackerDriven);
    regressionInputs[9]  = static_cast<float>(electronClass);
    regressionInputs[10] = static_cast<float>(isEB);
    
    // retrieve combination weight
    recalib_energy = tas::els_p4().at(elsIdx).e();
    
    float weight = 0.;
    if(eOverP>0.025 && fabs(momentum-energy)<15.*sqrt(momentumError*momentumError + energyError*energyError)) {
      weight = m_forest->GetResponse(regressionInputs);
      if(weight>1.) weight = 1.;
      else if(weight<0.) weight = 0.;
    }

    if(momentumError!=999. || weight==0.) {
      recalib_energy = weight*momentum + (1.-weight)*energy;
    }
}
