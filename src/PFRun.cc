#include "PFRun.h"

#include "TFile.h"
#include "TTree.h"
#include <iostream>
//#include <stdio.h>

PFRun::PFRun()
{
  fPulsarNumber=0;
  fRunNumber=0;
  fIsADC=false;

  fNFreq=512;
  fFreq0=109.584;
  fFreq511=112.084;
  fL511=267.6564006;
  fL0=273.7625931;
  fDL=0.0119494;

  fTau=2.4576;
  fNBins=1;
  fNPeriods=1;
  fNBinsPerPeriod=1;
  fYear=2013;
  fMonth=1;
  fDay=1;
  fHour=1;
  fMinute=1;
  fSecond=1.0;

  fReadRAWFile=NULL;
}

PFRun::~PFRun() 
{}

int PFRun::ReadRunHeader(std::string rootfile)
{
  fReadRAWFile=new TFile(rootfile.c_str());
  
  if (fReadRAWFile->IsZombie()) {
    std::cout<<"PFRun::ReadRunHeader;  root file "<<rootfile.c_str()<<" not found, STOP"<<std::endl;
    return 1;
  }

  //read Run Header
  TTree* RunHeader=(TTree*)fReadRAWFile->Get("RunHeader");
  if (RunHeader->IsZombie()) {
    std::cout<<"PFRun::ReadRunHeader;  RunHeader tree not found in "<<rootfile.c_str()<<", STOP"<<std::endl;
    return 1;
  }

  RunHeader->SetBranchAddress("fNPeriods",&fNPeriodsInput);
  RunHeader->SetBranchAddress("fNBinsPerPeriod",&fNBinsPerPeriod);

  RunHeader->SetBranchAddress("fTau",&fTau);
  RunHeader->SetBranchAddress("fPeriod",&fPeriod);

  RunHeader->SetBranchAddress("fYear",&fYear);
  RunHeader->SetBranchAddress("fMonth",&fMonth);
  RunHeader->SetBranchAddress("fDay",&fDay);
  RunHeader->SetBranchAddress("fHour",&fHour);
  RunHeader->SetBranchAddress("fMinute",&fMinute);
  RunHeader->SetBranchAddress("fSecond",&fSecond);
  
  RunHeader->SetBranchAddress("fFreq0",&fFreq0);
  RunHeader->SetBranchAddress("fFreq511",&fFreq511);

  RunHeader->GetEntry(0);

  RunHeader->Delete();

  fL511=3e10*pow(fFreq511*1e6,-1);
  fL0=3e10*pow(fFreq0*1e6,-1);
  fDL=(fL0-fL511)*pow((fNFreq-1),-1);
  
  fNBinsInput=fNPeriodsInput*fNBinsPerPeriod;
  
  fNPeriods=fNPeriodsInput-3;
  fNBins=fNPeriods*fNBinsPerPeriod;

  return 0;
}


int PFRun::ReadRAWSignal()
{
  if (!fReadRAWFile->IsOpen()) {
    std::cout<<"PFRun::ReadRAWSignal;  run is not initialised, add PFRun::ReadRunHeader function"<<std::endl;
    return 1;
  }

  // std::cout<<"assasa"<<std::endl;

  std::vector<float> fMeans;
  char tmp[100];
  for (int y=0; y<fNFreq; y++){
    sprintf(tmp,"sigTimeProfile_freqID_%d",y);
    fHPerBandSignal.push_back((TH1F*)fReadRAWFile->Get(tmp));
    fMeans.push_back(0);
    
    // std::cout<<y<<std::endl;

    for (int i=0; i<fNBinsInput; i++){
      fMeans[y]+=fHPerBandSignal[y]->GetBinContent(i);
    }
    fMeans[y]=fMeans[y]/fNBins;
    if (fMeans[y]!=0) fHPerBandSignal[y]->Scale(pow(fMeans[y],-1));
  }

  return 0;
}


int PFRun::CloseRun()
{
  fHPerBandSignal.clear();
  if (fReadRAWFile->IsOpen()) fReadRAWFile->Close();
  return 0;
}
