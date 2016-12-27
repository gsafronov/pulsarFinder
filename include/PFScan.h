#ifndef PFSCAN_H
#define PFSCAN_H

#include "PFRun.h"

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

class PFScan : PFRun
{
 public:
  PFScan(int nScanPoints, float DM0, float nScanStep, int rebinFactor);
  // ~PFScan();

  int DoScan_CPU(int nThreads);
  int DoScan_GPU(int nThreadsPerBlock);
  
  int CloseGPU();

  int DoCompensation_CPU(int iThread, int iStep);
  
  int CleanSignal(float cutSpike, float cutBadRun);

  float GetDM(int iStep) {return fDM0+iStep*fScanStep;};
  float GetNScanPoints() {return fNScanPoints;};
  float GetDM0() {return fDM0;};
  int GetRebinFactor() {return fRebinFactor;}

  int InitScan(std::string, bool doFFT);
  int SaveOutput(std::string);
  int CloseScan();

 private:
  int DoCompensation_GPU(int iStep, int nThreadsPerBlock);

  int DoCuFFT(int iStep);
  int DoRooFFT(int iStep);

  int ReadRun(std::string);

  int Rebin(int rebinFactor);
  
 protected:
  //PFRun fRun;
  int fNScanPoints;
  float fDM0;
  float fScanStep;

  int fNThreads;
  
  bool fReadRun;
  bool fCleanSignal;
  bool fDoScan;
  bool fUseGPU;
  bool fDoFFT;
  bool fSaveResults;

  int fRebinFactor;
  bool fIsRebin;
  
  float* fSigArray;
  float* fDev_SigArray;

  float* fCompSigArray;
  
  std::vector<TH1F*> fHCompSig;
  //  std::vector<TH1F*> fHCompSig_CPU;
  std::vector<TH1F*> fHCompSig_FFTImage;

  TFile* fScanOutFile;

  TH1F* fHCompTiming;
  TH1F* fHFFTTiming;

};

#endif
