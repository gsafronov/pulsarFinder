#include "DMAnalysis.h"
#include "TMath.h"
#include "Math/DistFunc.h"
#include "TCanvas.h"
#include <limits>

double CorrectForLEE(double significance, int nTrials)
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

DMAnalysis::DMAnalysis(TFile* inf, std::vector<std::string> rid, int nper, int nBinsPerPeriod, int npo, float dm, float dmL, float tauIn)
{
  _infile=inf;
  _runID=rid;
  _nFiles=_runID.size();
  _nPointsToScan=npo;
  _nPeriods=nper;
  _nBinsPerPeriod=nBinsPerPeriod;
  _dm0=dm;
  _scanStep=(dmL-_dm0)/_nPointsToScan;
  _nBins=_nPeriods*_nBinsPerPeriod;
  _tau=tauIn;

  std::cout<<"nBins: "<<_nBins<<std::endl;
  std::cout<<"nPeriods: "<<_nPeriods<<std::endl;
  

  std::vector<std::string> labels; //[_nFiles];
  for (int i=0; i<_nFiles; i++) {
    labels.push_back(_runID[i].substr(0,6));
  }
  
  maxSNRperRun=new TH2F("maxSNR_PerRun_vs_DM","maxSNR_PerRun_vs_DM",_nPointsToScan,_dm0, dmL, _nFiles, 0, _nFiles);
  
  bckgRMSperRun=new TH2F("bckgRMS_PerRun","bckgRMS_PerRun_vs_DM",_nPointsToScan,_dm0,dmL,_nFiles, 0, _nFiles);

  hFFTPeriod_SNRvsDM=new TH1F("hFFTPeriod_SNRvsDM","hFFTPeriod_SNRvsDM",_nPointsToScan, _dm0, _dm0+_nPointsToScan*_scanStep);

  for(int k=0; k<_nFiles; k++){
    maxSNRperRun->GetYaxis()->SetBinLabel(k+1,labels[k].c_str());
    bckgRMSperRun->GetYaxis()->SetBinLabel(k+1,labels[k].c_str());
  }
  
  char tmp[100];
  for (int i=0; i<_nPointsToScan; i++){
    sprintf(tmp,"hFoldedProfile_%d",i);
    hFoldedProfile.push_back(new TH1F(tmp,tmp,1000,0,1000));
    
    sprintf(tmp,"hBackground_%d",i);
    hBackground.push_back(new TH1F(tmp,tmp,300,400,700));
    hBackground[i]->TH1F::Sumw2();
    
    //    sprintf(tmp,"hFullBackground_%d",i);
    //    hFullBackground.push_back(new TH1F(tmp,tmp,300,400,700));
    //    hFullBackground[i]->TH1F::Sumw2();
    
    sprintf(tmp,"fImageSum_%d",i);
    hFImageSum.push_back(new TH1F(tmp,tmp,_nBins,0,_nBins));
    
    sprintf(tmp,"fImageSumFreq_%d",i);
    hFImageSumFreq.push_back(new TH1F(tmp,tmp,_nBins,pow((_nBins)*_tau*0.001,-1), pow(_tau*0.001,-1)));
    
    sprintf(tmp,"fImageSumPeriod_%d",i);
    float periodBins[_nBins+1];
    for (int j=0; j<=_nBins; j++){
      periodBins[_nBins-j]=pow((j+1)/(_nBins*_tau*0.001),-1);
    }
    hFImageSumPeriod.push_back(new TH1F(tmp,tmp,_nBins,periodBins));
    
    sprintf(tmp,"hIndPulses_%d",i);
    hIndPulses.push_back(new TH1F(tmp,tmp,100,0,100));
    
    sprintf(tmp,"SNR_vs_PeriodNumber_iDM=%d",i);
    hSNR2D.push_back(new TH2F(tmp,tmp,_nPeriods,0,_nPeriods,_nFiles,0,_nFiles));
    for(int k=0; k<_nFiles; k++){
      hSNR2D[i]->GetYaxis()->SetBinLabel(k+1,labels[k].c_str());
    }
    
    sprintf(tmp,"hFFTRes_%d",i);
    hFFTRes.push_back(new TH1F(tmp,tmp,_nBins,0,_nBins));
    
    sprintf(tmp,"hFFTResScale_%d",i);
    hFFTResScale.push_back(new TH1F(tmp,tmp,_nBins,0,1));
  }

  std::cout<<"_nPointsToScan histograms created"<<std::endl;
  
  for (int i=0; i<_nPeriods+200; i++){
    sprintf(tmp,"hBckgPerPeriod_%d",i);
    hBckgPerPeriod.push_back(new TH1F(tmp,tmp,300,400,700));
    hBckgPerPeriod[i]->TH1F::Sumw2();
  }
  
  //FFT 
  //period representation
  //define binning:
  double blim[_nBins+1];
  for (int k=0; k<=_nBins; k++){
    if (k==0) blim[_nBins-k]=pow((k+1)*pow(_nBins,-1),-1)+100;
    else blim[_nBins-k]=pow(k*pow(_nBins,-1),-1);
  }
  
  for (int i=0; i<_nPointsToScan; i++){
    sprintf(tmp,"hFFTResPeriod_%d",i);
    hFFTResPeriod.push_back(new TH1F(tmp,tmp,_nBins,blim));
  }

  std::cout<<"end of constructor"<<std::endl;
  
}

int DMAnalysis::write(TFile* ff, std::string path)
{
  ff->cd(path.c_str());
  //TCanvas* canv[_nPointsToScan];
  char tmp[100];
  for (int i=0; i<_nPointsToScan; i++){
    hBackground[i]->Write();
    hFoldedProfile[i]->Write();
    //    hFullBackground[i]->Write();
    if (i==10) hBckgPerPeriod[10]->Write();
    //sprintf(tmp,"peaks_DM=%d",i);
    //canv[i]=new TCanvas(tmp,tmp);
    //canv[i]->SetBatch(kTRUE);
    //canv[i]->Divide(2,1);
    //TVirtualPad* p1=canv[i]->cd(1);
    //p1->SetLogy();
    hIndPulses[i]->SetLineWidth(2);
    hIndPulses[i]->Write();
    //canv[i]->cd(2);
    hSNR2D[i]->SetStats(kFALSE);
    hSNR2D[i]->SetDrawOption("COLZ");
    hSNR2D[i]->Write();
    //hSNR2D[i]->SetDrawOption("COLZ");
    //canv[i]->Write();
    hFFTRes[i]->Write();
    hFFTResScale[i]->Write();
    hFFTResPeriod[i]->Write();
    
    hFImageSum[i]->Write();
    hFImageSumFreq[i]->Write();
    hFImageSumPeriod[i]->Write();
  }
  ff->cd();
  //  TCanvas* maxSNR=new TCanvas("maxSNR","maxSNR");
  // maxSNR->SetBatch(kTRUE);
  // maxSNR->cd();
  maxSNRperRun->SetStats(kFALSE);
  maxSNRperRun->SetDrawOption("COLZ");
  maxSNRperRun->Write();
  //  maxSNRperRun->SetDrawOption("COLZ");
  //  maxSNR->Write();
  //  maxSNRperRun->Write();
  //  TCanvas* bckgQual=new TCanvas("bckgQual","bckgQUal");
  //  bckgQual->SetBatch(kTRUE);
  //  bckgQual->cd();
  bckgRMSperRun->SetStats(kFALSE);
  bckgRMSperRun->SetDrawOption("COLZ");
  bckgRMSperRun->Write();
  //  bckgRMSperRun->SetDrawOption("COLZ");
  //  bckgRMSperRun->Draw("Q");
  //  bckgQual->Write();
  bckgRMSperRun->Write();
  hFFTPeriod_SNRvsDM->Write();
}

int DMAnalysis::foldByPeriod(int iStep)
{
  //may also lauch this function into CUDA it would be:
  //each thread will fold 1 value of period
  //total number of fold operations: nValues*nScanPoints
  int summedPeriods=0;
  for (int i=0; i<_nFiles; i++){
    std::stringstream ss;
    ss << iStep;
    std::string step = ss.str();
    std::string hname=_runID[i]+"/"+"sumFreq_"+step;
    TH1F* histoToFold=(TH1F*)_infile->Get(hname.c_str());  
    if (histoToFold==NULL) {
      std::cout<<"DMAnalysis::foldByPeriod: histogram "<<hname<<" not found"<<std::endl;
      return 1;
    }
    
    hFoldedProfile[iStep]->SetBins(_nBinsPerPeriod,0,_nBinsPerPeriod);
    
    summedPeriods=0;
    for (int i=1; i<=floor(_nBins/_nBinsPerPeriod)*_nBinsPerPeriod; i++){
      if (i%_nBinsPerPeriod==0&&histoToFold->GetBinContent(i)>0) summedPeriods++;
      if (histoToFold->GetBinContent(i)>0) {
	hFoldedProfile[iStep]->Fill((i-1)%_nBinsPerPeriod,histoToFold->GetBinContent(i));
	hBackground[iStep]->Fill(histoToFold->GetBinContent(i),1);
      }
    }
  }
  return summedPeriods;
}



std::vector<float> DMAnalysis::searchForPeaks(int iStep, int window, int rebinFactor)
{      
  std::vector<float> output;
  float minPval=100;
  float minPos=-1;
  int _nFiles=_runID.size();
  
  float SBmaxTot=-1;
  
  for (int ifi=0; ifi<_nFiles; ifi++){ 
    float SBmaxPerRun=-1;

    std::stringstream ss;
    ss << iStep;
    std::string step = ss.str();
    std::string hname=_runID[ifi]+"/"+"sumFreq_"+step;
    
    TH1F* hAllPeriods=(TH1F*)_infile->Get(hname.c_str());  
    hAllPeriods->Rebin(rebinFactor);
    
    //find background
    hBackground[iStep]->Reset();
    for (int p=0; p<_nPeriods; p++)
      {
	hBckgPerPeriod[p]->Reset();
      }
    
    int nPeriodsCur=hAllPeriods->GetNbinsX()/_nBinsPerPeriod;
    float nPeriodCur1=hAllPeriods->GetNbinsX()/_nBinsPerPeriod;
    
    int nPeriodsToRun=std::min(nPeriodsCur,_nPeriods);
    
    for (int i=0; i<nPeriodsToRun*_nBinsPerPeriod; i++){
      int iPer=floor(i/_nBinsPerPeriod);
      hBckgPerPeriod[iPer]->Fill(hAllPeriods->GetBinContent(i),1);
      //      hFullBackground[iStep]->Fill(hAllPeriods->GetBinContent(i),1);
    }
    
    //background is gauss with a good precision
    //fit core to determine PDF
    TF1 ff("ff","gaus");
    float mean[nPeriodsToRun];
    float rms[nPeriodsToRun];
    float fitMean[nPeriodsToRun];
    float fitSigma[nPeriodsToRun];
    for (int p=0; p<nPeriodsToRun; p++){
      fitMean[p]=-1;
      fitSigma[p]=-1;
      
      int nMiss=hBckgPerPeriod[p]->GetBinContent(0)+hBckgPerPeriod[p]->GetBinContent(hBckgPerPeriod[p]->GetNbinsX()+1);
      if (nMiss>=_nBinsPerPeriod-10||hBckgPerPeriod[p]->GetEntries()==0) continue;
      
      mean[p]=hBckgPerPeriod[p]->GetMean();
      rms[p]=hBckgPerPeriod[p]->GetRMS();
      
      ff.SetParameter(0,hBckgPerPeriod[p]->GetEntries());
      ff.SetParameter(1,mean[p]);
      ff.SetParameter(2,rms[p]);
      
      hBckgPerPeriod[p]->Fit(&ff,"QN","",mean[p]-2*rms[p],mean[p]+2*rms[p]);
      
      fitMean[p]=ff.GetParameter(1);
      fitSigma[p]=ff.GetParameter(2);
      
      float chi2=ff.GetChisquare();
      int ndof=ff.GetNDF();
      
      //safety agains failed fit
      
      if (chi2/ndof>5||fabs((fitMean[p]-mean[p])/rms[p])>2) {
	fitMean[p]=mean[p];
	fitSigma[p]=rms[p];
      }
      
    }
    
    //fill hIndPulses histograms
    float SBmax=1;
    float SBmin=1000;
    int perBuf=0;
    for (int i=1; i<=hAllPeriods->GetNbinsX(); i++){
      //skip first period
      int iPer=floor((i*rebinFactor)/_nBinsPerPeriod);
      if (iPer==0) continue;
      if (iPer>=nPeriodsToRun-1) continue;
      if (iPer!=perBuf&&fitSigma[perBuf]!=0){
	hSNR2D[iStep]->Fill(perBuf, ifi, SBmax);
	
	perBuf=iPer;
	
	hIndPulses[iStep]->Fill(SBmax,1);
	
	SBmax=CorrectForLEE(SBmax,_nBinsPerPeriod);
	
	//if (SBmax >= DBL_MAX || SBmax <= -DBL_MAX) std::cout<<"iStep: "<<iStep<<"   SBmax:     "<<SBmaxPerRun<<"   "<<SBmax<<"   run: "<<_runID[ifi]<<"  period: "<<iPer<<std::endl;
	if (SBmax>SBmaxTot) SBmaxTot=SBmax;
	if (SBmax>SBmaxPerRun) {
	  SBmaxPerRun=SBmax;
	}
	SBmax=-1;
	SBmin=1000;
      }
      //calculate SNR
      float bico=hAllPeriods->GetBinContent(i);
      if (bico<=0) bico=fitMean[iPer];
      if ((bico-fitMean[iPer])/fitSigma[iPer]>SBmax) SBmax=(bico-fitMean[iPer])/fitSigma[iPer];
      if ((bico-fitMean[iPer])/fitSigma[iPer]<SBmin) SBmin=(bico-fitMean[iPer])/fitSigma[iPer];
      //if (fabs(SBmaxPerRun)>100000) std::cout<<"mamin: "<<_runID[ifi]<<"  "<<i<<"  "<<iPer<<"   content: "<<bico<<"  fmean: "<<fitMean[iPer]<<"  fsigma: "<<fitSigma[iPer]<<"  "<<SBmax<<"   "<<SBmin<<std::endl;
    }
    
    SBmaxPerRun=CorrectForLEE(SBmaxPerRun,_nPeriods);
    
    maxSNRperRun->Fill(_dm0+iStep*_scanStep,ifi,SBmaxPerRun);
    hAllPeriods->Delete();
  }
  
  output.push_back(SBmaxTot);
  output.push_back(1);
  return output;
}

void DMAnalysis::fImageMerge(int iStep)
{
  int _nFiles=_runID.size();
  for(int ifi=0; ifi<_nFiles; ifi++){
    std::stringstream ss;
    ss << iStep;
    std::string step = ss.str();
    std::string hname=_runID[ifi]+"/"+"sumFreq_fImage_"+step;
    
    TH1F* hFImage=(TH1F*)_infile->Get(hname.c_str());
    
    if(hFImage->GetNbinsX()==_nBins) hFImageSum[iStep]->Add(hFImage,1);
    else std::cout<<"merge fImage:  skipping wrong binning  "<<ifi<<"  "<<hname<<"    Nbins: "<<hFImage->GetNbinsX()<<"     NbinsOrig: "<<hFImageSum[iStep]->GetNbinsX()<<std::endl;
  }
  
  //convert to understandable scales
  //scale x-axis by 1/T for frequency
  for (int i=0; i< hFImageSum[iStep]->GetNbinsX(); i++){
    hFImageSumFreq[iStep]->SetBinContent(i+1,hFImageSum[iStep]->GetBinContent(i+1));
    hFImageSumPeriod[iStep]->SetBinContent(i+1,hFImageSum[iStep]->GetBinContent(hFImageSum[iStep]->GetNbinsX()-i));
  }

  hFImageSumPeriod[iStep]->Scale(pow(_nFiles,-1));
}


std::vector<float> DMAnalysis::localSgfScan(int iStep, int window, int rebinFactor)
{
  int _nFiles=_runID.size();
  hFoldedProfile[iStep]->Rebin(rebinFactor);
  
  if (window>_nBinsPerPeriod) std::cout<<"localSgfScan: warning - window larger than number of bins"<<std::endl;
  
  //define output: minimal P-value and its position
  std::vector<float> output;
  
  //distributions of signal in integration interval within period  
  char tmp[100];
  std::vector<TH1F> sigInBin;
  std::vector<TH1F> sigInBinFull; //[_nBinsPerPeriod];
  for (int h=0; h<_nBinsPerPeriod; h++){
    sprintf(tmp,"sigInBin_%d",h);
    sigInBin.push_back(TH1F(tmp,tmp,10000,300,1000)); //folded by DM
    sprintf(tmp,"sigInBinFull_%d",h);
    sigInBinFull.push_back(TH1F(tmp,tmp,1000,-1,10)); //full
  }
  
  //background is gauss with a good precision
  //fit core to determine PDF
  TF1 ff("ff","gaus");
  float mean=hBackground[iStep]->GetMean();
  float rms=hBackground[iStep]->GetRMS();
  
  hBackground[iStep]->Fit(&ff,"QN","",mean-5*rms,mean+5*rms);
  float fitMean=ff.GetParameter(1);
  float fitSigma=ff.GetParameter(2);
  int nFitEvents=hBackground[iStep]->GetEntries();
    
  //p-value of the fit is small due to low statistics in several bins, need to rebin
  //estimate peak significance for summed histogram
  
  //set errors
  int iPeriod=0;
  for (int j=0; j<_nFiles; j++){
    std::stringstream ss;
    ss << iStep;
    std::string step = ss.str();
    std::string hname=_runID[j]+"/"+"sumFreq_"+step;
    TH1F* histo=(TH1F*)_infile->Get(hname.c_str());
    for (int i=1; i<=_nBins; i++){
      if ((i-1)%_nBinsPerPeriod==0) iPeriod++;
      sigInBin[(i-1)%_nBinsPerPeriod].Fill(histo->GetBinContent(i),1);
    }
  }
  
  for (int i=1; i<=_nBinsPerPeriod; i++){
    float rms=sigInBin[i-1].GetRMS();
    float mean=sigInBin[i-1].GetMean();
    hFoldedProfile[iStep]->SetBinError(i,rms*pow(mean,-1)*hFoldedProfile[iStep]->GetBinContent(i)*pow(sqrt(_nPeriods),-1));
  }
  
  //local fit to define local significance: fit with const backg in <windowWidth> bins
  TF1 ffH0("ffH0","[0]");
  ffH0.SetParameter(0,fitMean);
  
  float minPval=100;
  float minPos=-1;
  
  for (int i=1; i<=_nBinsPerPeriod-window; i++){
    float loBord=(i-1);
    float upBord=(i-1)+window;
    hFoldedProfile[iStep]->Fit(&ffH0,"QN","",loBord,upBord);
    float chi=ffH0.GetChisquare();
    int ndf=ffH0.GetNDF();
    if (TMath::Prob(chi,ndf)<minPval){
      minPval=TMath::Prob(chi,ndf);
      minPos=(i-1)+window/2;
    }
  }
  
  output.push_back(minPval);
  output.push_back(minPos);
  
  std::cout<<"iStep:  "<<iStep<<"    max excess (p-val, significance):  "<<minPval<<"  "<<ROOT::Math::normal_quantile(minPval,1)<<std::endl;

  return output;
}

int DMAnalysis::doFourier(int iStep)
{
  for (int j=0; j<_nFiles; j++){
    std::stringstream ss;
    ss << iStep;
    std::string step = ss.str();
    std::string hname=_runID[j]+"/"+"sumFreq_"+step;
    TH1F* histo=(TH1F*)_infile->Get(hname.c_str());
    TH1F* histoForFour=new TH1F("histoForFour","histoForFour",_nBins-_nBinsPerPeriod*2,0,_nBins-_nBinsPerPeriod*2); //histogram with non-empty bins
    int iEntr=0;
    for (int i=1; i<_nBins+1; i++){
      float bico=histo->GetBinContent(i);
      if (bico>0){
	iEntr++;
	histoForFour->SetBinContent(iEntr,bico);
      }
    }
    
    TH1* hm=new TH1F("hm","hm",_nBins,0,_nBins);
    TVirtualFFT::SetTransform(0);
    hm = histoForFour->FFT(hm, "MAG");
    hFFTRes[iStep]->Add(hm,1);
    hm->Delete();
    histoForFour->Delete();
  }
  
  for (int j=1; j<=hFFTRes[iStep]->GetNbinsX(); j++){
    hFFTResScale[iStep]->SetBinContent(j,hFFTRes[iStep]->GetBinContent(j));
    hFFTResPeriod[iStep]->SetBinContent(hFFTRes[iStep]->GetNbinsX()-j,hFFTRes[iStep]->GetBinContent(j));
  }      
}


float DMAnalysis::searchPeaksInFFT(int iStep, float value, float windowSignal, float windowNoise)
{
  //estimate background
  TH1F bckg("bckg","bckg",10000,0,10000);
  int iStart=hFImageSumPeriod[iStep]->FindBin(value-windowNoise);
  int iFinish=hFImageSumPeriod[iStep]->FindBin(value+windowNoise);
  std::cout<<iStart<<"    "<<iFinish<<std::endl;
  for (int i=iStart; i<iFinish; i++){
    //skip 1s harmonics:
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<0.255&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.245) continue;
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<0.34&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.325) continue;
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<0.51&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.48) continue;
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<1.02&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.96) continue;
    
    //    if (fImageSumPeriod->GetBinCenter(i)<0.&&fImageSumPeriod->GetBinCenter(i)>0.49) continue;
    bckg.Fill(hFImageSumPeriod[iStep]->GetBinContent(i),1);
  }
  //find mean and sigma of noise 
  float bckgMean=bckg.GetMean();
  float bckgSigma=bckg.GetRMS();
  std::cout<<"FFT period:   bckgMean "<<bckgMean<<"  "<<bckgSigma<<std::endl;
  
  //check SNR in the region of period
  float maxSNR=-10;
  iStart=hFImageSumPeriod[iStep]->FindBin(value-windowSignal);
  iFinish=hFImageSumPeriod[iStep]->FindBin(value+windowSignal);
  for (int i=iStart; i<iFinish; i++) {
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<0.255&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.245) continue;
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<0.34&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.325) continue;
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<0.51&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.48) continue;
    if (hFImageSumPeriod[iStep]->GetBinCenter(i)<1.02&&hFImageSumPeriod[iStep]->GetBinCenter(i)>0.96) continue;
    float snr=(hFImageSumPeriod[iStep]->GetBinContent(i)-bckgMean)*pow(bckgSigma,-1);
    if (snr>maxSNR) maxSNR=snr;
  }
  std::cout<<"FFT period:    maxSNR "<<maxSNR<<std::endl;
  hFFTPeriod_SNRvsDM->SetBinContent(iStep+1,maxSNR);
  return maxSNR;  
}
