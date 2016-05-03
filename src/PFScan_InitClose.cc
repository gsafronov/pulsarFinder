#include "PFScan.h"
#include "TF1.h"

PFScan::PFScan(int nScanPoints, float DM0, float scanStep)
{
  fNScanPoints=nScanPoints;
  fDM0=DM0;
  fScanStep=scanStep;
  
  fReadRun=false;
  fCleanSignal=false;
  fDoScan=false;
  fDoFFT=false;
  fSaveResults=false;
  
  //  fSigArray=NULL;
  
  //  fSaveResults=true;
  //  if (fScanOutFile->IsZombie()||

}

//PFScan::~PFScan()
//{}

int PFScan::SaveOutput(std::string rootfile)
{
  if (!fDoScan) {
    std::cout<<"PFScan::SaveOutput; scan is not initialised, STOP"<<std::endl;
    return 1;
  }
  fSaveResults=true;
  fScanOutFile=new TFile(rootfile.c_str(),"RECREATE");

  if (fDoScan){
    fHCompTiming->Write();
    fHFFTTiming->Write();
    fScanOutFile->mkdir("DMCompHist");
    fScanOutFile->cd("DMCompHist");
    for (int j=0; j<fNScanPoints; j++){
      fHCompSig[j]->Write();
      if (fDoFFT) fHCompSig_FFTImage[j]->Write();
    }
    
  }
  fScanOutFile->Close();
  return 0;
}

int PFScan::InitScan(std::string rootfile, bool doFFT)
{
  std::cout<<"PFScan::InitScan"<<std::endl;
  
  //  fHCompSig.clear();
  
  fDoScan=true;

  fDoFFT=doFFT;
  
  fHCompTiming=new TH1F("fHCompTiming",
			"DM compensation time per point per run (s)",
			10000,0,20);
  fHFFTTiming=new TH1F("fHFFTTiming",
		       "FFT time per point per run (s)",
		       10000,0,20);

  ReadRun(rootfile);
  
  return 0;
}

int PFScan::ReadRun(std::string rootfile)
{
  fReadRun=true;
  ReadRunHeader(rootfile);
  ReadRAWSignal();
  
  //read file contents into the memory:
  //input array must be of sie of input histograms
  size_t size_input=(fNBinsInput*fNFreq)*sizeof(float);
  fSigArray=(float*)malloc(size_input);
  for (int i = 0; i < fNFreq; ++i){
    for (int j = 0; j < fNBinsInput; ++j){
      fSigArray[i*fNBinsInput+j] = fHPerBandSignal[i]->GetBinContent(j+1);
    }
  }
  //allocate compensated signal array:
  //it can be of length several periods less
  size_t size_output=(fNBins)*sizeof(float);
  fCompSigArray=(float*)malloc(size_output);
  for (int i = 0; i<fNBins; i++) fCompSigArray[i]=0;

  char tmp[100];
  for (int i = 0; i<fNScanPoints; i++){
    sprintf(tmp,"fHCompSig_%d",i);
    fHCompSig.push_back(new TH1F(tmp,tmp,fNBins,0,fNBins));
  }

  if (fDoFFT){
    for (int i = 0; i<fNScanPoints; i++){
      sprintf(tmp,"fHCompSig_FFTImage_%d",i);
      fHCompSig_FFTImage.push_back(new TH1F(tmp,tmp,fNBins,0,fNBins));
    }
  }
  //CloseRun();
}

int PFScan::CloseScan()
{
  std::cout<<"PFScan::CloseScan"<<std::endl;
  
  CloseRun();

  fHCompSig.clear();
  fHCompSig_FFTImage.clear();
  free(fCompSigArray);
  free(fSigArray);
  
  //  fHCompSig.clear();
  return 0;
}
