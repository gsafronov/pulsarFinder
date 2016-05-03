#ifndef PFANALYSIS_H
#define PFANALYSIS_H

#include "PFRun.h"
#include "PFScan.h"
#include "TH2F.h"

class PFAnalysis : PFRun
{
 public:
  PFAnalysis(std::string anOutputName, 
	     PFScan* scan,
	     bool doFFT);
  
  //  int InitAnalysis(std::vector<std::string> runs, int nScanPoints);

  int AddRun(std::string runID, 
	     std::string fnameRAW, 
	     std::string fname);

  //  std::vector<float> FindPeaks(TH1F* histogram);
  
  //  int SumPeriodsInRun();

  int SumRuns();
    
 private:
  float ProcessCompSig(int iStep, 
		       std::string runID, 
		       TFile* scanFile, 
		       TH2F* fHMaxSNRperPeriod);

  std::vector<float> FindPeaks(TH1F* histogram);
  //  int ProcessFFT(int iStep);
  //  int ReadScanOutput(int iStep);

  int FoldPeriods(TH1F* hIn, 
		  TH1F* hOut);
  
  int ProcessFFT(int iStep, 
		 std::string runID, 
		 TFile* scanFile);

  int ConvertFFTScales(TH1F* initHist, TH1F* freqHist, TH1F* periodHist);
  float SumHarmonics(TH1F* initHist, TH1F* sumHist);

  int RemoveSecond();
  
  TFile* fOutputFile;
  // TFile fScanOutFile;

  PFScan* fScan;
  int fNScanPoints;
  int fDM0;

  bool fDoFFT;

  std::vector<std::string> fRunID;
  
  //vector of 1D max SNR per run each DM point 
  std::vector<TH1F> fVHMaxSNRPerRun;
  
  //vector of 1D max SNR per run int FFT for each DM point 
  std::vector<TH1F> fVHMaxFFTSNRPerRun;

  //vector of folded periods by DM
  std::vector<TH1F> fVHFoldedPeriod;

  //  TH2F* fHMaxSNRperRun;
  
  TH1F fHFoldedPeriod;

  std::string fCurrentRunID;

  int fNumberOfRuns;

  //add various status flags for safe execution

  //FFT results
  std::vector<TH1F> fVHFFTSumRuns;
  std::vector<TH1F> fVHFFTSumRunsFreq;
  std::vector<TH1F> fVHFFTSumRunsPeriod;
  std::vector<TH1F> fVHFFTSumRunsSumHarm;


  std::vector<TH1F> fVHFFTNoise;
  std::vector<TH1F> fVHFFTLowFreq;
  
};
  
#endif  
