#include "PFAnalysis.h"
#include "PFRun.h"
#include "TMath.h"
#include "Math/DistFunc.h"
#include "TF1.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

double CorrectForLEE(double significance, 
		     int nTrials)
{
  double retSNR;
  if (significance<8) {
    double Cmod;
    if (significance>0) Cmod=ROOT::Math::gaussian_cdf(significance);
    else {
      Cmod=1;
    }
    double pval=2*(1-Cmod);
    double plee=1-pow(1-pval,nTrials);
    //find new quantile
    double SNR=ROOT::Math::normal_quantile(1-0.5*plee,1);
    retSNR=SNR;
  }
  else retSNR=significance;
  return retSNR;
}

PFAnalysis::PFAnalysis(std::string anOutputFname, 
		       PFScan* scan,
		       bool doFFT)
{
  fOutputFile=new TFile(anOutputFname.c_str(), "RECREATE");
  
  fScan=scan;
  fNScanPoints=scan->GetNScanPoints();
  fDM0=scan->GetDM0();
  fDoFFT=doFFT;
  fNumberOfRuns=0;
  fRebinFactor=scan->GetRebinFactor();

  TH1::SetDefaultSumw2(kFALSE);
  for (int i=0; i<fNScanPoints; i++){
    char tmp[100];
    sprintf(tmp, "fVHFoldedPeriod_%d", i);
    fVHFoldedPeriod.push_back(TH1F(tmp,tmp,1,0,1));
    fVHFoldedPeriod[i].TH1F::Sumw2(kFALSE);
    sprintf(tmp, "fVHFFTSumRuns_%d", i);
    
    fVHFFTSumRuns.push_back(TH1F(tmp,tmp,1,0,1));
    sprintf(tmp, "fVHFFTSumRunsFreq_%d", i);
    fVHFFTSumRunsFreq.push_back(TH1F(tmp,tmp,1,0,1));
    sprintf(tmp, "fVHFFTSumRunsPeriod_%d", i);
    fVHFFTSumRunsPeriod.push_back(TH1F(tmp,tmp,1,0,1));
    sprintf(tmp, "fVHFFTSumRunsSumHarm_%d", i);
    fVHFFTSumRunsSumHarm.push_back(TH1F(tmp,tmp,1,0,1));
  }
}

int PFAnalysis::ConvertFFTScales(TH1F* initHisto, TH1F* freqHisto, TH1F* periodHisto)
{
  //convert FFT to understandable scales
  //scale x-axis by 1/T for frequencyfloat periodBins[_nBins+1];
  double periodBins[fNBins+1];
    for (int j=0; j<=fNBins; j++){
      periodBins[fNBins-j]=pow((j+1)/(fNBins*fTau*0.001),-1);
    }
  freqHisto->SetBins(fNBins,pow((fNBins)*fTau*0.001,-1), pow(fTau*0.001,-1));
  periodHisto->SetBins(fNBins,periodBins);

  for (int i=0; i<initHisto->GetNbinsX(); i++){
    freqHisto->SetBinContent(i+1,initHisto->GetBinContent(i+1));
    periodHisto->SetBinContent(i+1,initHisto->GetBinContent(initHisto->GetNbinsX()-i));
  }

  return 0;  
}

int PFAnalysis::SumRuns()
{
  std::cout<<"PFAnalysis::SumRuns"<<std::endl;
  if (fNumberOfRuns==0) return 1;

  //array of run names
  std::vector<std::string> labels;
  for (int i=0; i<fRunID.size(); i++) {
    labels.push_back(fRunID[i].substr(0,6));
  }
  
  TH1::SetDefaultSumw2(kFALSE);

  //maximum SRN per run for each DM
  TH2F hMaxSNRperRun("hMaxSNRperRun",
		      "Maximum SNR per run vs DM",
		      fNScanPoints,fDM0,fScan->GetDM(fNScanPoints),
		      fRunID.size(),0,fRunID.size());
  

  for (int i=0; i<fVHMaxSNRPerRun.size(); i++){
    for (int j=0; j<fNScanPoints; j++){
      //     std::cout<<"run: "<<fRunID[i]<<"    DM: "<<fScan->GetDM(j)<<"    "
      //  	       <<fVHMaxSNRPerRun[i].GetBinContent(j+1)<<std::endl;
      hMaxSNRperRun.Fill(fScan->GetDM(j),i,fVHMaxSNRPerRun[i].GetBinContent(j+1));
    }
  }

  //  std::cout<<"popopo"<<std::endl;

  
  fVHMaxSNRPerRun.clear();
  
  //normalize folded periods
  
  for (int i=0; i<fVHFoldedPeriod.size(); i++){
    fVHFoldedPeriod[i].Scale(pow(fRunID.size(),-1));
  }


  if (fDoFFT){
    
    //convert FFT to understandable scales
    //scale x-axis by 1/T for frequencyfloat periodBins[_nBins+1];
    for (int iStep=0; iStep<fNScanPoints; iStep++){
      ConvertFFTScales(&fVHFFTSumRuns[iStep],&fVHFFTSumRunsFreq[iStep], &fVHFFTSumRunsPeriod[iStep]);
      //sum harmonics
      SumHarmonics(&fVHFFTSumRunsFreq[iStep], &fVHFFTSumRunsSumHarm[iStep]);
    }
  }
  TH1F hMaxFFTSNRTotal("hMaxFFTSNRTotal",
		       "Maximum sum FFT SNR vs DM",
		       fNScanPoints,fDM0,fScan->GetDM(fNScanPoints));

  //  std::cout<<"momomo"<<std::endl;
  
  TH2F hMaxFFTSNRperRun("hMaxFFTSNRperRun",
			"Maximum FFT SNR per run vs DM",
			fNScanPoints,fDM0,fScan->GetDM(fNScanPoints),
			fRunID.size(),0,fRunID.size());
  if (fDoFFT){
    for (int j=0; j<fNScanPoints; j++){
      float maxSNR=0;
      for (int i=0; i<fVHMaxFFTSNRPerRun.size(); i++) {
	hMaxFFTSNRperRun.Fill(fScan->GetDM(j),i,fVHMaxFFTSNRPerRun[i].GetBinContent(j+1));
	if (fVHMaxFFTSNRPerRun[i].GetBinContent(j+1)>maxSNR) maxSNR=fVHMaxFFTSNRPerRun[i].GetBinContent(j+1);
      }
      hMaxFFTSNRTotal.Fill(fScan->GetDM(j),maxSNR);
    }
    fVHMaxFFTSNRPerRun.clear();
  }
  
  char tmp[100];
  if (fOutputFile->IsOpen()) {
    fOutputFile->cd();

    for(int k=0; k<fRunID.size(); k++){
      hMaxSNRperRun.GetYaxis()->SetBinLabel(k+1,labels[k].c_str());
      if (fDoFFT) hMaxFFTSNRperRun.GetYaxis()->SetBinLabel(k+1,labels[k].c_str());
    }
    //    std::cout<<"zozozo"<<std::endl;
    hMaxSNRperRun.SetStats(kFALSE);
    hMaxSNRperRun.SetDrawOption("COLZ");
    hMaxSNRperRun.Write();
    if (fDoFFT){
      hMaxFFTSNRperRun.SetStats(kFALSE);
      hMaxFFTSNRperRun.SetDrawOption("COLZ");
      hMaxFFTSNRperRun.Write();
      hMaxFFTSNRTotal.GetSumw2()->Set(0);
      hMaxFFTSNRTotal.Write();
    }

    fOutputFile->mkdir("FoldedPeriod");
    fOutputFile->cd("FoldedPeriod");
    for (int j=0; j<fVHFoldedPeriod.size(); j++){
      sprintf(tmp,"Sum of Periods in all runs for DM = %f",fScan->GetDM(j));
      fVHFoldedPeriod[j].SetTitle(tmp);
      fVHFoldedPeriod[j].GetSumw2()->Set(0);
      fVHFoldedPeriod[j].Write();
    }
    //std::cout<<"youououiuiu"<<std::endl;
    if (fDoFFT){
      fOutputFile->mkdir("FFT");
      fOutputFile->cd("FFT");
      for (int j=0; j<fNScanPoints; j++){
	fVHFFTSumRunsFreq[j].Write();
	fVHFFTSumRunsPeriod[j].Write();
	fVHFFTSumRunsSumHarm[j].Write();
	fVHFFTNoise[j].Write();
	fVHFFTLowFreq[j].Write();
      }
    }
    else 
    fOutputFile->Close();
  }
  //  std::cout<<"yououou"<<std::endl;
  return 0;
}

int PFAnalysis::AddRun(std::string runID, 
		       std::string fnameRAW, 
		       std::string fnameScanOut)
		       
{
  std::cout<<"PFAnalysis::AddRun"<<std::endl;
  
  TH1::SetDefaultSumw2(kFALSE);

  fRunID.push_back(runID);
  fNumberOfRuns++;

  ReadRunHeader(fnameRAW);
  
  RecalculateNumbers(fRebinFactor);

  // std::cout<<fNBinsPerPeriod<<" "<<fNPeriods<<std::endl;
  
  //max SNR per period histogram
  std::string hname="fHMaxSNRperPeriod_"+runID;
  std::string htitle="Maximum SNR per period for run "+runID;
  TH2F hMaxSNRperPeriod(hname.c_str(),htitle.c_str(),
			fNScanPoints, 0, fNScanPoints, 
			fNPeriods,0,fNPeriods);

  //rebin some histograms when 1'st run is added

  if (fRunID.size()==1){
    for (int i=0; i<fNScanPoints; i++){
      fVHFoldedPeriod[i].SetBins(floor(fNBinsPerPeriod),0,0.001*fTau*(floor(fNBinsPerPeriod)));
      if (fDoFFT) fVHFFTSumRuns[i].SetBins(fNBins,0,fNBins);
      
      /*
      fVHFFTSumRunsFreq[i].SetBins(fNBins,pow((fNBins)*fTau*0.001,-1), pow(fTau*0.001,-1));
      fVHFFTSumRunsPeriod[i].SetBins(fNBins,periodBins);
      fVHFFTSumRunsSumHarm[i].SetBins(fNBins,pow((fNBins)*fTau*0.001,-1), pow(fTau*0.001,-1));
      */
    }
  }   

  //add max per iStep 1D histogram
  hname="fHMaxSNR_"+runID;
  fVHMaxSNRPerRun.push_back(TH1F(hname.c_str(),hname.c_str(),
				 fNScanPoints, 0, fNScanPoints));
  if (fDoFFT){
    hname="fHMaxFFTSNR_"+runID;
    fVHMaxFFTSNRPerRun.push_back(TH1F(hname.c_str(),hname.c_str(),
				      fNScanPoints, 0, fNScanPoints));
  }
  //std::cout<<fNScanPoints<<"   "<<fDM0<<"   "<<fScan->GetDM(fNScanPoints)<<"   "<<fRunID.size()-1<<std::endl;

  if (fOutputFile->IsOpen()){
    std::string directory=runID+"/foldedPeriod";
    fOutputFile->cd();
    fOutputFile->mkdir(runID.c_str());
    fOutputFile->mkdir(directory.c_str());
  }

  TFile scanOutFile(fnameScanOut.c_str());
  for (int i=0; i<fNScanPoints; i++){
    float maxSNR=ProcessCompSig(i, runID, &scanOutFile, &hMaxSNRperPeriod);
    //    std::cout<<"maxSNR: "<<maxSNR<<std::endl;
    fVHMaxSNRPerRun[fNumberOfRuns-1].Fill(i,maxSNR);

    if (fDoFFT) {
      ProcessFFT(i,runID,&scanOutFile);
    }
      
  }
  scanOutFile.Close();
  
  if (fOutputFile->IsOpen()){
    fOutputFile->cd();
    fOutputFile->cd(runID.c_str());
    hMaxSNRperPeriod.GetSumw2()->Set(0);
    hMaxSNRperPeriod.SetStats(kFALSE);
    hMaxSNRperPeriod.SetDrawOption("COLZ");
    hMaxSNRperPeriod.Write();
  }
  
  return 0;
}

float PFAnalysis::ProcessCompSig(int iStep, 
				 std::string runID, 
				 TFile* scanOutFile, 
				 TH2F* hMaxSNRperPeriod)
{
  float returnValue=-1;
  //  std::cout<<"PFAnalysis::ProcessCompSig   step: "<<iStep<<std::endl;
 
  //get iStep histogram
  char tmp[100];
  sprintf(tmp, "DMCompHist/fHCompSig_%d", iStep);
  TH1F* hStepCompSig=(TH1F*)scanOutFile->Get(tmp);

  hStepCompSig->TH1F::Sumw2(kFALSE);

  //run FindPeaks
  std::vector<float> vMaxSNRperPeriod;
  vMaxSNRperPeriod=FindPeaks(hStepCompSig);

  //add to foldPeriod
  //write into runID/foldedPeriod directory
  if (fOutputFile->IsOpen()){
    std::string directory=runID+"/foldedPeriod";
    fOutputFile->cd(directory.c_str());
    
    TH1F foldedHist("foldedHist","foldedHist",floor(fNBinsPerPeriod), 0, 0.001*fTau*(floor(fNBinsPerPeriod)));
    //    std::cout<<floor(fNBinsPerPeriod)<<"  "<<fTau<<"  "<<fTau*floor(fNBinsPerPeriod)<<std::endl;
    foldedHist.TH1F::Sumw2(kFALSE);
    FoldPeriods(hStepCompSig,&foldedHist);
    fVHFoldedPeriod[iStep].Add(&foldedHist);
    
    sprintf(tmp,"hFoldedPeriod_%d",iStep);
    foldedHist.SetName(tmp);
    std::string title="Sum of all periods in run "+runID;
    foldedHist.SetTitle(title.c_str());
    foldedHist.GetSumw2()->Set(0);
    foldedHist.Write();
  }

  //hStepCompSig->Delete();
  
  //fill fHMaxSNRperPeriod, find maximum per run
  float maxPerRun=0;
  for (int i=0; i<vMaxSNRperPeriod.size(); i++){
    //    std::cout<<"max per period "<<i<<" "<<vMaxSNRperPeriod[i]<<std::endl;
    hMaxSNRperPeriod->Fill(iStep, i, vMaxSNRperPeriod[i]);
    if (vMaxSNRperPeriod[i]>maxPerRun) maxPerRun=vMaxSNRperPeriod[i];
  }

  //  std::cout<<"before lee corr: "<<maxPerRun<<std::endl;
  
  maxPerRun=CorrectForLEE(maxPerRun,fNPeriods);

  //  std::cout<<"after lee corr: "<<maxPerRun<<std::endl;
  
  returnValue=maxPerRun;

  //resize
  //fHMaxSNRperRun->GetYaxis()->Set(
  //  fHMaxSNRperRun->Fill(iStep, fRunID.size(), maxPerRun);
  
  //  std::cout<<"or here?"<<std::endl;

  //  std::cout<<"Leave compsig"<<std::endl;
  //resize and fill fHMaxSNRperRun
  return returnValue;
}


int PFAnalysis::ProcessFFT(int iStep, 
			     std::string runID, 
			     TFile* scanOutFile)
{
  //get iStep histogram
  char tmp[100];
  sprintf(tmp, "DMCompHist/fHCompSig_FFTImage_%d", iStep);
  TH1F* hStepCompSigFFT=(TH1F*)scanOutFile->Get(tmp);

  if (hStepCompSigFFT->GetNbinsX()==fNBins) fVHFFTSumRuns[iStep].Add(hStepCompSigFFT,1);
  else std::cout<<"merge fImage:  skipping wrong binning  "<<"  "<<fRunID[fRunID.size()-1]<<"    Nbins: "<<hStepCompSigFFT->GetNbinsX()<<"     NbinsOrig: "<<fNBins<<std::endl;

  TH1F hFreq("hFreq","hFreq",1,0,1);
  TH1F hPeriod("hPeriod","hPeriod",1,0,1);
  ConvertFFTScales(hStepCompSigFFT,&hFreq,&hPeriod);

  TH1F hSum("hSum","hSum",1,0,1);
  float maxFFTSNR=SumHarmonics(&hFreq, &hSum);
  fVHMaxFFTSNRPerRun[fVHMaxFFTSNRPerRun.size()-1].Fill(iStep,maxFFTSNR);

  return 0;
}


float PFAnalysis::SumHarmonics(TH1F* hFFTFreq, TH1F* hSum)
{
  //returnValue;
  float maxSNR=0;

  //FFT processing parameters: lowest frequency, number of harmonics to sum, highest frequency
  float fLowestFreq=0.2;
  int fNumSumHarm=4;
  float fHighestFreq=20;
  int fNPeaksToRemove=20;

  //fit the background
  TF1 bckgFunc("bckgFunc","[0]+[1]/x+[2]/(x*x)+[3]/(x*x*x)+[4]/(x*x*x*x)",0.05,100);
  bckgFunc.SetParameter(0,10000);
  bckgFunc.SetParameter(1,10000);
  bckgFunc.SetParameter(2,10000);
  bckgFunc.SetParameter(3,10000);
  bckgFunc.SetParameter(4,10000);
  hFFTFreq->Fit("bckgFunc","QN","",0.05,fHighestFreq);
  
  //find bckg sigma
  int iStartBin=hFFTFreq->FindBin(fLowestFreq);
  int iEndBin=hFFTFreq->FindBin(fHighestFreq);
  //  std::cout<<iStartBin<<"      "<<iEndBin<<std::endl;
  char tmp[100];
  sprintf(tmp,"fVHFFTNoise_%d",(int)fVHFFTNoise.size());
  fVHFFTNoise.push_back(TH1F(tmp,tmp,1000,-100000,100000));
  
  for (int i=iStartBin; i<iEndBin; i++){
    float subtrNoise=hFFTFreq->GetBinContent(i)-bckgFunc.Eval(hFFTFreq->GetBinCenter(i));
    fVHFFTNoise[fVHFFTNoise.size()-1].Fill(subtrNoise,1);
  }
  float sigmaNoiseRaw=fVHFFTNoise[fVHFFTNoise.size()-1].GetRMS();
  
  //fill new small histogram with only few first harmonics
  sprintf(tmp,"fVHFFTLowFreq_%d",(int)fVHFFTLowFreq.size());
  fVHFFTLowFreq.push_back(TH1F(tmp,tmp,iEndBin,0,
			       hFFTFreq->GetBinLowEdge(iEndBin)+hFFTFreq->GetBinWidth(iEndBin)));
  for (int j=1; j<=fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetNbinsX(); j++){
    float bicoBckg=bckgFunc.Eval(fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetBinCenter(j));
    if (j>iStartBin) fVHFFTLowFreq[fVHFFTLowFreq.size()-1].SetBinContent(j,hFFTFreq->GetBinContent(j)-bicoBckg);
  }

  
  /*
    //SUBTRACTION OF 1s HARMONICS SHOWED BAD PERFORMANCE - PDF IS NOT EXACTLY GAUSS, CONSEQUENTLY 
    //LARGE DEVIATIONS PERSIST ON THE HISTOGRAM AFTER SUBTRACTION. INSTEAD CUTTING THE 1S REGION AWAY SHOULD BE USED

  TF1 peakFunc("peakFunc","[0]*exp(-(x-[1])*(x-[1])/(2*([2]*[2])))",0.5,100);
  
  hSum->Reset();
  hSum->Add(hFFTFreq);

  std::cout<<"second peaks:  ";

  std::vector<TF1> vPeakFunc;

  for (int i=0; i<3; i++){
    peakFunc.SetParameter(0,100000);
    peakFunc.SetParameter(1,i+1+0.001);
    peakFunc.SetParameter(2,0.005);
    fVHFFTLowFreq[fVHFFTLowFreq.size()-1].Fit("peakFunc","Q","",i+1-0.04,i+1+0.05);
    vPeakFunc.push_back(peakFunc);

    std::cout<<peakFunc.GetParameter(1)<<"/"<<peakFunc.GetParameter(2)<<"    ";
  }
  
  //subtract peaks from FFT
  for (int i=0; i<vPeakFunc.size(); i++) {
    //std::cout<<"peakfunc: "<<vPeakFunc[i].GetParameter(1)<<std::endl;
    int peakBin=fVHFFTLowFreq[fVHFFTLowFreq.size()-1].FindBin(vPeakFunc[i].GetParameter(1));
    //go 5 bins left and 5 bins right and subtract the fit function
    for (int j=peakBin-5; j<=peakBin+5; j++) {
      float funcVal=vPeakFunc[i].Eval(fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetBinCenter(j));
      float bico=fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetBinContent(j);
      fVHFFTLowFreq[fVHFFTLowFreq.size()-1].SetBinContent(j,bico-funcVal);
    }
  }
    
  */

  //replace 1s peaks with const
  for (int i=0; i<fNPeaksToRemove; i++){
    int iBin=fVHFFTLowFreq[fVHFFTLowFreq.size()-1].FindBin(i+1.01);
    for (int j=iBin-5; j<=iBin+5; j++){
      fVHFFTLowFreq[fVHFFTLowFreq.size()-1].SetBinContent(j,0);
    }
  }

  //sum harmonics
  hSum->Reset();
  hSum->SetBins(fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetNbinsX(),0,
		fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetBinLowEdge(iEndBin)+fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetBinWidth(iEndBin));
  hSum->Add(&fVHFFTLowFreq[fVHFFTLowFreq.size()-1]);

  for (int iHarm=2; iHarm <= fNumSumHarm; iHarm++){
    //B
    //    TH1F squeezedFFT=*hFFTFreq;
    //    squeezedFFT.Rebin(iHarm);
    for (int iBin=1; iBin<floor(fNBins/iHarm); iBin++){
      float squeezedCo=0;
      for (int k=iBin*iHarm; k<(iBin+1)*iHarm; k++) squeezedCo+=fVHFFTLowFreq[fVHFFTLowFreq.size()-1].GetBinContent(k);
      float sumco=hSum->GetBinContent(iBin)+squeezedCo;
      hSum->SetBinContent(iBin, sumco);
    }   
  }

  float maxFreq=0;
  for (int i=1; i<=hSum->FindBin(fHighestFreq/fNumSumHarm); i++){
    if (hSum->GetBinContent(i)>maxSNR) {
      maxSNR=hSum->GetBinContent(i);
      maxFreq=hSum->GetBinCenter(i);
    }
  }
  
  maxSNR=maxSNR/(sigmaNoiseRaw*sqrt(fNumSumHarm));
  
  //  std::cout<<"FFT SNR: "<<maxSNR<<"     corr. freq: "<<maxFreq<<std::endl;

  return maxSNR;
  
  //    hSum->Reset();
  //    hSum->Add(hFFTFreq);

  //squeeze the frequency representation by x2, x3, x4, x5 etc and add to initial freq histogram
  //A) split the bins in the initial histogram by 1/2, 1/3, 1/4, etc and divide contents between new bins
  //B) merge bins in added histogram
  
  //  TH1F summedFFT=*hFFTFreq;
  //  TH1F summedFFT("summedFFT","summedFFT",fNBins, pow((fNBins)*fTau*0.001,-1), pow(fTau*0.001,-1));

  //  hSum->Reset();
  //  hSum->Add(hFFTFreq);
  //  std::cout<<"FIRST FIT: "<<std::endl;
  //  hSum->Fit("multiPeakFunc","","",0.5,10);

  //clean the initial histogram

  //  std::cout<<"REFIT: "<<std::endl;
  //  hSum->Fit("multiPeakFunc","","",0.5,10);
  /*
  for (int iHarm=2; iHarm < 10 ; iHarm++){
    //B
    //    TH1F squeezedFFT=*hFFTFreq;
    //    squeezedFFT.Rebin(iHarm);
    for (int iBin=1; iBin<floor(fNBins/iHarm); iBin++){
      float sqco=0;
      for (int k=iBin*iHarm; k<(iBin+1)*iHarm; k++) sqco+=hFFTFreq->GetBinContent(k);
      float sumco=hSum->GetBinContent(iBin)+sqco;
      hSum->SetBinContent(iBin, sumco);
    }   
  }
  */
  //  return 0;  
}

std::vector<float> PFAnalysis::FindPeaks(TH1F* hCompSig)
{
  //  std::cout<<"PFAnalysis::FindPeaks"<<std::endl;
  std::vector<float> vectorToReturn;
  
  //find background
  std::vector<TH1F> hBckgPerPeriod;
  char tmp[100];
  int iPerBuf=0;
  hBckgPerPeriod.push_back(TH1F(tmp,tmp,300,400,700));
  for (int i=0; i<hCompSig->GetNbinsX(); i++){
    sprintf(tmp,"hBckgPer_%d",i);
    int iPer=floor((float)i/(float)fNBinsPerPeriod);
    if (iPer!=iPerBuf){
      hBckgPerPeriod.push_back(TH1F(tmp,tmp,300,400,700));
      iPerBuf=iPer;
    }
    //    std::cout<<iPer<<"    "<<i<<"   "<<hCompSig->GetBinContent(i+1)<<std::endl;
    hBckgPerPeriod[iPer].Fill(hCompSig->GetBinContent(i+1),1);
  }

  //background is gauss with a good precision
  //fit core to determine PDF
  TF1 ff("ff","gaus");
  float mean[fNPeriods];
  float rms[fNPeriods];
  float fitMean[fNPeriods];
  float fitSigma[fNPeriods];
  float fitChiNDOF[fNPeriods];
  
  //  std::cout<<"yopop"<<std::endl;
  for (int p=0; p<fNPeriods; p++){
    fitMean[p]=-1;
    fitSigma[p]=-1;
    //    std::cout<<"blabla   "<<p<<"  "<<fNPeriods<<"  "<<hBckgPerPeriod.size()<<std::endl;
    //check if there are many hits not falling into 400-700 amplitude range
    int nMiss=hBckgPerPeriod[p].GetBinContent(0)+
      hBckgPerPeriod[p].GetBinContent(hBckgPerPeriod[p].GetNbinsX()+1);
    if (nMiss>=fNBinsPerPeriod-10||hBckgPerPeriod[p].GetEntries()==0){
      //      std::cout<<"misisisisi"<<std::endl;
      continue;
    }
    
    mean[p]=hBckgPerPeriod[p].GetMean();
    rms[p]=hBckgPerPeriod[p].GetRMS();
    
    ff.SetParameter(0,hBckgPerPeriod[p].GetEntries());
    ff.SetParameter(1,mean[p]);
    ff.SetParameter(2,rms[p]);

    //   std::cout<<"heere   "<<p<<"  "<<hBckgPerPeriod.size()<<std::endl;
    
    hBckgPerPeriod[p].Fit(&ff,"QN","",mean[p]-2*rms[p],mean[p]+2*rms[p]);
    
    fitMean[p]=ff.GetParameter(1);
    fitSigma[p]=ff.GetParameter(2);
        
    float chi2=ff.GetChisquare();
    int ndof=ff.GetNDF();

    fitChiNDOF[p]=(float)chi2/(float)ndof;
    
    //safety agains failed fit
    
    if ((float)chi2/(float)ndof>5||fabs((fitMean[p]-mean[p])/rms[p])>2||hBckgPerPeriod[p].GetEntries()<50) {
      fitMean[p]=mean[p];
      fitSigma[p]=rms[p];
    }
  }
  
  //fill hIndPulses histograms
  float SNRmax=-1;
  float SNRmin=1000;
  iPerBuf=0;
  for (int i=0; i<hCompSig->GetNbinsX(); i++){
    int iPer=floor((float)i/(float)fNBinsPerPeriod);
    float binco;
    int bID;
    if (iPer!=iPerBuf&&fitSigma[iPerBuf]!=0){
      //      std::cout<<" bibi:    "<<iPer<<"  "<<iPerBuf<<"   "<<i<<"  "<<bID<<"  "<<binco<<"   "<<fitMean[iPerBuf]<<"  "<<mean[iPerBuf]<<"   "<<fitSigma[iPerBuf]<<"   "<<rms[iPerBuf]<<"   "<<fitChiNDOF[iPerBuf]<<std::endl;
      iPerBuf=iPer;

      //if (SNRmax<0) std::cout<<" bibi:    "<<binco<<"   "<<fitMean[iPerBuf]<<std::endl;
      //<<"  "<<<<std::endl;
      //std::cout<<"before corr: "<<fNBinsPerPeriod<<"  "<<SNRmax<<std::endl;
      SNRmax=CorrectForLEE(SNRmax,floor(fNBinsPerPeriod));
      //std::cout<<"after corr: "<<SNRmax<<std::endl;
      
      
      vectorToReturn.push_back(SNRmax);
      
      SNRmax=-1;
      SNRmin=1000;
    }
    //calculate SNR
    float bico=hCompSig->GetBinContent(i+1);
    if (bico<=0) bico=fitMean[iPer];
    if ((bico-fitMean[iPer])/fitSigma[iPer]>SNRmax) {
      SNRmax=(bico-fitMean[iPer])/fitSigma[iPer];
      binco=bico;
      bID=i;
    }
    if ((bico-fitMean[iPer])/fitSigma[iPer]<SNRmin) SNRmin=(bico-fitMean[iPer])/fitSigma[iPer];
    //if (fabs(SBmaxPerRun)>100000) std::cout<<"mamin: "<<_runID[ifi]<<"  "<<i<<"  "<<iPer<<"   content: "<<bico<<"  fmean: "<<fitMean[iPer]<<"  fsigma: "<<fitSigma[iPer]<<"  "<<SBmax<<"   "<<SBmin<<std::endl;
  }
  //  std::cout<<"rerere"<<std::endl;
  return vectorToReturn;

}

//////

int PFAnalysis::FoldPeriods(TH1F* histIn, 
			    TH1F* histOut)
{
  //TH1F foldedHist("foldedHist","foldedHist",fNBinsPerPeriod, 0, fNBinsPerPeriod); 
  //  fVHFoldedPeriod.push_back(TH1F("fVHFP","fVHFP",fNBinsPerPeriod, 0, fNBinsPerPeriod));
  if (histIn->IsZombie()||histOut->IsZombie()) return 1;

  for (int i=0; i<histIn->GetNbinsX(); i++){
    int iPeriod=floor((float)(i)/(float)fNBinsPerPeriod);
    float cbin=histIn->GetBinCenter(i+1);
    //    std::cout<<i<<"  "<<iPeriod<<"  "<<fNBinsPerPeriod<<"  "<<cbin<<"  "<<fTau*0.001*fTau*fNBinsPerPeriod*iPeriod<<"  "<<cbin-fTau*0.001*fNBinsPerPeriod*iPeriod<<"  "<<histIn->GetBinContent(i+1)<<std::endl;
    histOut->Fill(cbin-fTau*0.001*fNBinsPerPeriod*iPeriod,histIn->GetBinContent(i+1));
  }
  histOut->Scale(pow(fNPeriods,-1));

  return 0;
}


