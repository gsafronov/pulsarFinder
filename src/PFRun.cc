#include "PFRun.h"

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <stdio>

PFRun::PFRun()
{
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
  fScanOutFile=NULL;
  fAnalysisOutFile=NULL;
}

PFRun::PFRun() 
{}

PFRun::ReadRun(std::string rootfile)
{
  fReadRAWFile=new TFile(rootfile.c_str());
  
  if (fInputFile->IsZombie()) {
    std::cout<<"PFRun::ReadRun root file "<<rootfile.c_str()<<" not found"<<std::endl;
    return 1;
  }

  //read Run Header
  // std::cout<<"blah"<<std::endl;
  TTree* RunHeader=(TTree*)fInputFile->Get("RunHeader");
  RunHeader->SetBranchAddress("fNPeriods",&fNPeriods);
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
  
  fNBins=fNPeriods*fNBinsPerPeriod;
  
}
