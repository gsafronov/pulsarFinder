#ifndef PFRUN_H
#define PFRUN_H

#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

class PFRun
{
 public:
  PFRun();
  ~PFRun();

  int ReadRunHeader(std::string rootfile);
  int ReadRAWSignal();
  int CloseRun();

 protected:
  
  int RecalculateNumbers(int rebinFactor);
  //  int SetReadRAWOutFile();
  //  int SetScanOutFile();
  //  int SetAnalysisOutFile();

  TFile* GetReadRAWOutFile() {return fReadRAWFile;};
  //  TFile* GetScanOutFile() {return fScanOutFile};
  //  TFile* GetAnalysisOutFile() {return fAnalysisOutFile};
  
  TFile* fReadRAWFile;
  //  TFile* fScanOutFile;
  //  TFile* fAnalysisOutFile;

  int fNFreq;
  double fFreq0;
  double fFreq511;
  double fL511;
  double fL0;
  double fDL;
  
  double fTau;
  double fPeriod;
  int fNPeriods;
  int fNPeriodsInput;
  int fNBins;
  int fNBinsInput;
  float fNBinsPerPeriod;
  int fYear;
  int fMonth;
  int fDay;
  int fHour;
  int fMinute;
  double fSecond;

  int fPulsarNumber;
  int fRunNumber;  

  bool fIsADC;

  //signal, compensated signal and FFT images

  std::vector<TH1F*> fHPerBandSignal;
  /*
  TH1F* fHCompSig;
  TH1F* fHCompSig_cuFFT;
  TH1F* fHCompSig_rooFFT;
  */
  /*
  float* fSigArray;
  float* fDev_SigArray;
  */
};

#endif /* PFRUN_H */
