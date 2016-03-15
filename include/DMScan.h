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
  int readInfo(std::string txtfile);
  int initScan(std::string rootfile);
  int write(TFile* outFile, int iStep, std::string dir);
  int reset();
  //  int getInput(std::string rootfile);

  int rejectSpikes(float sigmaCut);

  float getTau() {return tau;}

  float getDM(int iStep) {return DM0+scanStep*iStep;}

  int getNPer() {return nPeriodsGlobal;}
  int getBPP() {return nBinsPerPeriod;}
 
  float getScanStart() {return DM0;}
  float getScanEnd() {return DM0+scanStep*(nPointsToScan);}

  int loadDataToGPU();
  int closeGPU();
  
  int sumFrequencies_CPU(int iThread, int iStep);
  int sumFrequencies_GPU(int iStep);

  int nPointsToScan;

  std::vector<TH1F*> sumFreq;
  std::vector<TH1F*> sumFreq_fImage;
  //  TFile* outputFile;
  std::vector<TH1F*> sigTimeProfile;

  std::vector<TH1F*> foldedProfile;
  std::vector<TH1F*> background;

  TProfile* gpuCodeTiming;
  
 private:

  int nFreq;
  int nThreads;
  int nPeriods;
  int nPeriodsGlobal;
  int nBins;
  int nBinsGlobal;
  int nBinsPerPeriod;
  int rebinFactor;
  int year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  float fsec;

  std::vector<float> means;

  float l511;
  float period;
  float DM0;
  float dL;
  float tau;
  float scanStep;
  TFile* inputFile;
  float* sigArray;
  float* d_sigArray;
  //  float* sigSum;

  //  TFile* outputFile;
};

#endif /* DMSCAN_H */
