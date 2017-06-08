#include "PFScan.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include <iomanip>

PFScan::PFScan(int nScanPoints, float DM0, float scanStep, int rebinFactor)
{
  fNScanPoints=nScanPoints;
  fDM0=DM0;
  fScanStep=scanStep;
  
  fReadRun=false;
  fCleanSignal=false;
  fDoScan=false;
  fDoFFT=false;
  fSaveResults=false;
  fRebinFactor=rebinFactor;
  fIsRebin=false;
  
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
    fHFRprofile->Write();
    fHFrequencyResponse->Write();
    fHMaskedFreqResp->Write();
    fScanOutFile->mkdir("DMCompHist");
    fScanOutFile->cd("DMCompHist");

    std::ofstream fftInfo;
    char tmp[100];
    //    sprintf(tmp,"fftOutput_%s.dat",fFileName.c_str());
    fftInfo.open("fftOutput.dat");
    
    for (int j=0; j<fNScanPoints; j++){
      fHCompSig[j]->GetXaxis()->SetTitle("seconds");
      fHCompSig[j]->Write();
      if (fDoFFT) {
	fftInfo<<"DM: "<<fDM0+j*fScanStep<<"\n";
	fftInfo<<"period: "<<fPeriod<<"\n";
	fftInfo<<"tau: "<<fTau<<"\n";
	fftInfo<<"rebinFactor: "<<fRebinFactor<<"\n";
	fftInfo<<"number of periods: "<<fNPeriods<<"\n"<<"\n";
	fftInfo<<"fast fourier image [bin center (Hz)     power]:"<<"\n";
	fHCompSig_FFTImage[j]->Write();
	for (int ib=0; ib<fHCompSig_FFTImage[j]->GetNbinsX(); ib++){
	  fftInfo<<std::setprecision(5)<<fHCompSig_FFTImage[j]->GetBinCenter(ib)<<"      "<<fHCompSig_FFTImage[j]->GetBinContent(ib)<<"\n";
	}
      }
    }
    fftInfo.close();
    
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

  fHFRprofile=new TH1F("fHFRprofile","fHFRprofile",1000,1e6,1e9);

  fHMaskedFreqResp=new TH1F("fHMaskedFreqResp","fHMaskedFreqResp",
			    fNFreq,0,fNFreq);
    
  ReadRun(rootfile);
  
  return 0;
}

int PFScan::ReadRun(std::string rootfile)
{
  fReadRun=true;
  fIsRebin=false;
  ReadRunHeader(rootfile);
  ReadRAWSignal();


  //allocate frequency mask
  size_t size_fmask=(fNFreq)*sizeof(int);
  fFreqMask=(int*)malloc(size_fmask);
  for (int i=0; i<fNFreq; i++){
    fFreqMask[i]=1;
  }
  //do rebin
  
  //read file contents into the memory:
  //input array must be of size of input histograms
  size_t size_input=(fNBinsInput*fNFreq)*sizeof(float);
  fSigArray=(float*)malloc(size_input);
  for (int i = 0; i < fNFreq; ++i){
    //    fHPerBandSignal[i]->Rebin(fRebinFactor);
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
    fHCompSig.push_back(new TH1F(tmp,tmp,fNBins,0,fNBins*0.001*fTau));
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

int PFScan::Rebin(int rebinFactor)
{
  fRebinFactor=rebinFactor;
  
  RecalculateNumbers(rebinFactor);

  std::vector<float> fMeans;
  //  fMeans.clear();
  for (int i=0; i<fNFreq; i++){
    int nBinsNew=floor((float)fHPerBandSignal[i]->GetNbinsX()/(float)rebinFactor)*rebinFactor;
    TH1F* hbuf=(TH1F*)fHPerBandSignal[i]->Clone();
    //    fHPerBandSignal[i]->Clone(&hbuf);
    fHPerBandSignal[i]->SetBins(nBinsNew,0,nBinsNew);
    for (int j=0; j<nBinsNew; j++)
      fHPerBandSignal[i]->SetBinContent(i+1,hbuf->GetBinContent(j+1));
    hbuf->Delete();
    fHPerBandSignal[i]->Rebin(rebinFactor);
    fMeans.push_back(0);
    for (int j=0; j<fNBinsInput; j++){
      fMeans[i]+=fHPerBandSignal[i]->GetBinContent(j+1);
    }
    fMeans[i]=fMeans[i]/fNBins;
    if (fMeans[i]!=0) fHPerBandSignal[i]->Scale(pow(fMeans[i],-1));
  }
  
  for (int i = 0; i<fNScanPoints; i++){
    fHCompSig[i]->SetBins(fNBins,0,fNBins*0.001*fTau);
    fHCompSig_FFTImage[i]->SetBins(fNBins,0,fNBins);
  }

  
  //read file contents into the memory:
  //input array must be of size of input histograms
  size_t size_input=(fNBinsInput*fNFreq)*sizeof(float);
  void* r1=std::realloc(fSigArray,size_input);
  for (int i = 0; i < fNFreq; ++i){
    //    fHPerBandSignal[i]->Rebin(fRebinFactor);
    for (int j = 0; j < fNBinsInput; ++j){
      fSigArray[i*fNBinsInput+j] = fHPerBandSignal[i]->GetBinContent(j+1);
    }
  }
  //allocate compensated signal array:
  //it can be of length several periods less
  size_t size_output=(fNBins)*sizeof(float);
  void* r2= std::realloc(fCompSigArray, size_output);
  for (int i = 0; i<fNBins; i++) fCompSigArray[i]=0;
  
  return 1;
}
