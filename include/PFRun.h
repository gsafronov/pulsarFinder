#ifndef PFRUN_H
#define PFRUN_H

class TFile

class PFRun
{
 public:
  PFRun();
  ~PFRun();

  int ReadRun(std::string rootfile);
  
 private:
  
  int SetReadRAWOutFile();
  int SetScanOutFile();
  int SetAnalysisOutFile();

  TFile* GetReadRAWOutFile() {return fReadRAWFile};
  TFile* GetScanOutFile() {return fScanOutFile};
  TFile* GetAnalysisOutFile() {return fAnalysisOutFile};
  
  TFile* fReadRAWFile;
  TFile* fScanOutFile;
  TFile* fAnalysisOutFile;

  int fNFreq;
  double fFreq0;
  double fFreq511;
  double fL511;
  double fL0;
  double fDL;
  
  double fTau;
  int fNPeriods;
  int fNBins;
  int fNBinsPerPeriod;
  int fYear;
  int fMonth;
  int fDay;
  int fHour;
  int fMinute;
  double fSecond;

  bool fIsADC;


  //signal, compensated signal and FFT images

  TH1F* fHSignal;
  TH1F* fHCompSig;
  TH1F* fHCompSig_cuFFT;
  TH1F* fHCompSig_rooFFT;

  /*
  float* fSigArray;
  float* fDev_SigArray;
  */
}
