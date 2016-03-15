#ifndef DMANALYSIS_H
#define DMANALYSIS_H

#include <cmath>
#include <bitset>
#include "TH1F.h"
#include "TH2F.h"
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
#include <sstream>
#include "TMath.h"
#include "TVirtualFFT.h"

class DMAnalysis
{
 public:
  DMAnalysis(TFile* inputfile, std::vector<std::string> rid, int nPer, int nBinsPerPeriod, int nScanPoints, float a, float b, float tau);
  ~DMAnalysis() {;}
  int write(TFile* ff, std::string dir);
  int getNScanPoints() {return _nPointsToScan;}
  int foldByPeriod(int iStep);
  std::vector<float> searchForPeaks(int iStep, int window, int rebinFactor);
  std::vector<float> localSgfScan(int iStep, int window, int rebinFactor);
  float searchPeaksInFFT(int iStep, float value, float windowSignal, float windowNoise);
    
  void fImageMerge(int iStep);
  int doFourier(int iStep);
  //  int llRatioTest(int iStep, int iFreq);
 
  std::vector<TH1F*> hFImageSum;
  std::vector<TH1F*> hFImageSumFreq;
  std::vector<TH1F*> hFImageSumPeriod;
  std::vector<TH1F*> hFoldedProfile;
  std::vector<TH1F*> hBackground;
  std::vector<TH1F*> hFullBackground;
  std::vector<TH1F*> hBckgPerPeriod;
  std::vector<TH1F*> hIndPulses;
  std::vector<TH2F*> hSNR2D;
  std::vector<TH1F*> hFFTRes;
  std::vector<TH1F*> hFFTResScale;
  std::vector<TH1F*> hFFTResPeriod;
  
  TH1F* fullFFT;
  TH2F* maxSNRperRun;
  TH2F* bckgRMSperRun;
  TH1F* hFFTPeriod_SNRvsDM;
  
 private:
  TFile* _infile;
  std::vector<std::string> _runID;
  int _nFiles;
  int _nPointsToScan;
  int _nPeriods;
  int _nBins;
  int _nBinsPerPeriod;
  float _tau;
  float _dm0;
  float _scanStep;

};

#endif
