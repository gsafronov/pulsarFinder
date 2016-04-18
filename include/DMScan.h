#ifndef DMSCAN_H
#define DMSCAN_H

#include <cmath>
#include <bitset>
#include "TROOT.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include <string.h>
#include <iostream>
#include <fstream> 
#include "TApplication.h"
#include "TStopwatch.h"
#include <stdio.h>      
#include <stdlib.h>
#include <cufft.h>

class DMScan
{
 public:

  //functions:
  //  int initialize();
  DMScan(std::string outputFileName, int nThreads, int nPointsToScan, float DM0, float step);
  ~DMScan();

  int initScan(std::string rootfile);
  int write(TFile* outFile, int iStep, std::string dir);

  int reset();

  int rejectSpikes(float sigmaCut);

  float getTau() {return fTau;}

  float getDM(int iStep) {return fDM0+fScanStep*iStep;}

  int getNPer() {return fNPeriods;}
  int getBPP() {return fNBinsPerPeriod;}
 
  float getScanStart() {return fDM0;}
  float getScanEnd() {return fDM0+fScanStep*(fNPointsToScan);}

  int loadDataToGPU();
  int closeGPU();
  
  int sumFrequencies_CPU(int iThread, int iStep);
  int sumFrequencies_GPU(int iStep);

  std::vector<TH1F*> sumFreq;
  std::vector<TH1F*> sumFreq_fImage;
  //  TFile* outputFile;
  std::vector<TH1F*> sigTimeProfile;

  std::vector<TH1F*> foldedProfile;
  std::vector<TH1F*> background;

  TProfile* gpuCodeTiming;
  
 private:

  int fNFreq;
  int fNThreads;
  int fNPeriods;
  int fNBins;
  int fNBinsPerPeriod;
  int fRebinFactor;
  int fYear;
  int fMonth;
  int fDay;
  int fHour;
  int fMinute;
 
  double fFreq0;
  double fFreq511;
  double fL511;
  double fL0;
  double fDL;

  double fSecond;

  std::vector<float> fMeans;

  double fPeriod;
  double fTau;
  
  int fNPointsToScan;
  float fScanStep;
  float fDM0;
  
  TFile* fInputFile;
  float* fSigArray;
  float* fDev_SigArray;
  //  float* sigSum;

  //  TFile* outputFile;
};

#endif /* DMSCAN_H */
