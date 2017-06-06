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

  //  std::cout<<"run PFRun constructor"<<std::endl;
  
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

  int nBinsPerPeriodRead;
  RunHeader->SetBranchAddress("fNPeriods",&fNPeriodsInput);
  RunHeader->SetBranchAddress("fNBinsPerPeriod",&nBinsPerPeriodRead);

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

  fNBinsPerPeriod=nBinsPerPeriodRead;
  fNBinsInput=fNPeriodsInput*fNBinsPerPeriod;
  //  std::cout<<fNBinsInput<<"    "<<fNBinsPerPeriod<<"  "<<fNPeriods<<std::endl;

  //  fNBinsInput=floor((float)(fNPeriodsInput*fNBinsPerPeriod)/(float)rebinFactor);
  //  fNBinsPerPeriod=(float)fNBinsInput/(float)fNPeriodsInput;

  //  std::cout<<"AFTER REBIN   N bins input: "<<fNBinsInput<<"    N bins per period: "<<fNBinsPerPeriod<<"  N periods: "<<fNPeriodsInput<<std::endl;
  
  fNPeriods=fNPeriodsInput-3;
  fNBins=floor(fNPeriods*fNBinsPerPeriod);

  fHFrequencyResponse=new TH1F("fHFrequencyResponse",
			       "fHFrequencyResponse",
			       fNFreq,0,fNFreq);
  
  return 0;
}


int PFRun::ReadRAWSignal()
{
  if (!fReadRAWFile->IsOpen()) {
    std::cout<<"PFRun::ReadRAWSignal;  run is not initialised, add PFRun::ReadRunHeader function"<<std::endl;
    return 1;
  }

  char tmp[100];
  for (int i=0; i<fNFreq; i++){
    sprintf(tmp,"sigTimeProfileInput_freqID_%d",i);
    TH1F*  fHPerBand=new TH1F(tmp,tmp,1,0,1);
    fHPerBandSignal.push_back(fHPerBand);
  }
  
  std::vector<float> fMeans;
  for (int y=0; y<fNFreq; y++){
    sprintf(tmp,"sigTimeProfile_freqID_%d",y);
    //get hist first
    fHPerBandSignal[y]=(TH1F*)fReadRAWFile->Get(tmp);
    //    fHPerBandSignal.push_back((TH1F*)fReadRAWFile->Get(tmp));

    //fill the frequency response
    TH1F hFreqSig("hFreqSig","hFreqSig",1000,1e6,1e9);
    for (int j=0; j<fNBins; j++){
      hFreqSig.Fill(fHPerBandSignal[y]->GetBinContent(j),1);
      //      std::cout<<fHPerBandSignal[y]->GetBinContent(j)<<std::endl;
    }
    //    std::cout<<hFreqSig.GetMean()<<std::endl;
    fHFrequencyResponse->SetBinContent(y+1,hFreqSig.GetMean());
    fHFrequencyResponse->SetBinError(y+1,hFreqSig.GetRMS()/sqrt(hFreqSig.GetEntries()));

    //rebin, remove some bins first
    
    fMeans.push_back(0);
    
    //   std::cout<<"yo:  "<<fHPerBandSignal[y]->GetNbinsX()<<std::endl;

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
  //  for (int i=0; i<fNFreq; i++) fHPerBandSignal[i]->Delete();
  fHPerBandSignal.clear();
  if (fReadRAWFile->IsOpen()) fReadRAWFile->Close();
  return 0;
}

int PFRun::RecalculateNumbers(int rebinFactor)
{
  fNBinsInput=floor((float)(fNPeriodsInput*fNBinsPerPeriod)/(float)rebinFactor);
  fNBinsPerPeriod=(float)fNBinsInput/(float)fNPeriodsInput;
  
  fNPeriods=fNPeriodsInput-3;
  fNBins=floor(fNPeriods*fNBinsPerPeriod);
  
  fTau=fTau*rebinFactor;
  return 1;
}
