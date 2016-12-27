#include "PFScan.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"

//interface to CPU multithreading 
struct CompThreadID{
  PFScan* pfscan;
  int iThread;
  int iStep;
};

void* DoCompensation_CPU_interface(void *threadarg)
{
  struct CompThreadID *ctID;
  
  ctID = (struct CompThreadID *) threadarg;
  
  int rv = (ctID->pfscan)->DoCompensation_CPU(ctID->iThread, ctID->iStep);

  pthread_exit(NULL);
}
//////////////////////////////////////////

int PFScan::CleanSignal(float sigmaCut, 
			float cutBadRun)
{
  std::cout<<"PFScan::CleanSignal"<<std::endl;

  fCleanSignal=true;
  
  TH1F sumSigRef("sumSigRef","sumSigRef",fHPerBandSignal[0]->GetNbinsX(),0,fHPerBandSignal[0]->GetXaxis()->GetXmax());

  TH1F backgRef("backgRef","backgRef",5010,-10,5000);
  TH1F backgPerPeriod("backgPerPeriod","backgPerPeriod",1000,0,1000);

  for (int y=0; y<fNFreq; y++) sumSigRef.Add(fHPerBandSignal[y],1);
      
  for (int y=0; y<fNBins; y++){
    if (sumSigRef.GetBinContent(y+1)>0) backgRef.Fill(sumSigRef.GetBinContent(y+1),1);
//      if (sumSigRef.GetBinContent(y+1)==0) std::cout<<"ZeroBin"<<std::endl;
  }

  //  double refMean=backgRef.GetMean();
  //  double refRMS=backgRef.GetRMS();

  TF1 gfit("gfit","gaus");
  backgRef.Fit(&gfit,"QN","",450,650);
  double refMean=gfit.GetParameter(1);
  double refRMS=gfit.GetParameter(2);
  
  std::cout<<"removing spikes     "<<"reference mean: "<<refMean<<"     RMS: "<<refRMS<<"    number of bins: "<<fNBins<<std::endl;
  
  TH1F backgAvg("backgAvg","backgAvg",1000,0,1000);
  
  for (int i=0; i<fNPeriodsInput; i++){
    for (int j=0; j<floor(fNBinsPerPeriod); j++){
      backgPerPeriod.Fill(sumSigRef.GetBinContent(floor(i*fNBinsPerPeriod)+j),1);
    }
    backgAvg.Fill(backgPerPeriod.GetMean(),1);
    //      std::cout<<"period: "<<i<<"   bin: "<<i*fNBinsPerPeriod<<"    average: "<<backgPerPeriod.GetMean()<<std::endl;
    backgPerPeriod.Reset();
  }
  
  //  std::cout<<"RMS of period average: "<<backgAvg.GetRMS()<<"     Mean of period average: "<<backgAvg.GetMean()<<std::endl;
  //  std::cout<<"if RMS > 50 do not take the run"<<std::endl;
  
  bool isBadRun=false;
  if (backgAvg.GetRMS()>cutBadRun) isBadRun=true; 
  if (backgAvg.GetMean()>fNFreq*1.2||backgAvg.GetMean()<fNFreq*0.8) isBadRun=true;
  
  for (int q=1; q<fNBinsInput+1; q++){
    float bico=sumSigRef.GetBinContent(q);
    if (bico>refMean+sigmaCut*refRMS) std::cout<<"SPIKE "<<"     value: "<<bico<<"  bin: "<<q<<std::endl;
    for (int y=0; y<fNFreq; y++){
      if (isBadRun){
	refMean=fNFreq;
	fHPerBandSignal[y]->SetBinContent(q,refMean/fNFreq);
	fSigArray[y*fNBinsInput+q-1]=refMean/fNFreq;
      }
      if (bico>refMean+sigmaCut*refRMS){
	fHPerBandSignal[y]->SetBinContent(q,refMean/fNFreq);
	fSigArray[y*fNBinsInput+q-1]=refMean/fNFreq;
      }
    }
  }
  //  std::cout<<"finish"<<std::endl;
  return 0;
}

int PFScan::DoCompensation_CPU(int iThread, 
			       int iStep)
{
  //std::cout<<"DOCOMPENSATIONthread#"<<iThread<<"  "<<std::endl;
  //  std::cout<<fNBins<<"  "<<fNBins/fNBinsPerPeriod<<std::endl;
  float dm=GetDM(iStep);
  //define which bins current thread would sum
  int startPeriod=iThread*floor((fNPeriods+10)/(fNThreads));//add 10 for to cover all bins
  int endPeriod=(iThread+1)*floor((fNPeriods+10)/(fNThreads));
  if (iThread==fNThreads-1) endPeriod=fNPeriods;
  //  std::cout<<fTau<<"  "<<fPeriod<<" "<<fNPeriods<<" "<<fNThreads<<"   start: "<<startPeriod<<"  "<<endPeriod<<std::endl;
  for (int i=floor(startPeriod*fNBinsPerPeriod)+1; i<floor(endPeriod*fNBinsPerPeriod)+1; i++){
      //(i<=fNBinsPerPeriod)||i>fNBins-fNBinsPerPeriod) continue;

  //for (int i=1; i<fNBins+1; i++){
    float bico=0;
    for (int y=0; y<fNFreq; y++) {
      //take sigTimeProfile[511-y] as 511-th is shorter wavelength
      //calculate delay wrt 511th for particular freq[511-y]
      float dT=4.6*(-fL511*fL511+pow(fL511+y*fDL,2))*dm*0.001; //covert to ms
      //calculate residual difference to nearest positive side pulse
      float dTnearest=dT-fPeriod*floor(dT*pow(fPeriod,-1));
      //move frequency band by -dTnearest bins, add lower bins to "upper side"
      float delta=dTnearest*pow(fTau,-1);
      //float bico1=fSigArray[((fNFreq-1)-y)*fNBinsInput+int(floor(i+delta))];
      //float bico2=fSigArray[((fNFreq-1)-y)*fNBinsInput+int(floor(i+delta)+1)];
      //float bico1=sigArray[((511-y)*nBinsInput+int(floor(i+delta)))%(nBinsInput*nFreq)];
      //float bico2=sigArray[((511-y)*nBinsInput+int((floor(i+delta)+1)))%(nBinsInput*nFreq)];

      float bico1=fHPerBandSignal[(fNFreq-1)-y]
	->GetBinContent(int(floor(i+delta))%fNBinsInput);
      float bico2=fHPerBandSignal[(fNFreq-1)-y]
	->GetBinContent(int((floor(i+delta)+1))%fNBinsInput);
      float lowerBinFrac=1-((i+delta)-floor(i+delta));
      float upperBinFrac=1-lowerBinFrac;
      if (floor(i+delta)+1<fNBinsInput) bico+=lowerBinFrac*bico1+upperBinFrac*bico2;
    }
    fCompSigArray[i-1]=bico;
    //    std::cout<<"thread: "<<iThread<<"   bin index: "<<i-1<<"   value: "<<bico<<std::endl;
    //    if (bico==bico) fHCompSig[iStep]->Fill(i-1,bico);
  }
    return 0;
}

int PFScan::DoScan_CPU(int nThreads)
{
  std::cout<<"PFScan::DoScan_CPU"<<std::endl;
  fNThreads=nThreads;

  //do rebin here
  if (!fIsRebin) {
    Rebin(fRebinFactor);
    fIsRebin=true;
  }
  
  void *status;
  pthread_t threads[fNThreads];
  
  struct CompThreadID thrData[fNThreads];

  for (int i=0; i<fNScanPoints; i++){
    TStopwatch stwch;
    if ((i+1)%10==0||i+1==1) std::cout<<"PFScan::DoScan_CPU; process point "<<i+1<<std::endl;
    for(int iTh=0; iTh < fNThreads; iTh++ ){
      thrData[iTh].iThread=iTh;
      thrData[iTh].iStep=i;
      //problem with accessing protected members in mutithreading mode
      //problematic call to the function though "this" pointer
      thrData[iTh].pfscan=this;
      
      int rc = pthread_create(&threads[iTh], NULL, DoCompensation_CPU_interface, (void *)&thrData[iTh]);
      if (rc){
	std::cout << "Error:unable to create thread," << rc << std::endl;
	exit(-1);
      }
      //else std::cout<<"thread #"<<iTh<<" created"<<std::endl;
    }
    for ( int iTh=0; iTh < fNThreads; iTh++ ){	    
      int rc = pthread_join(threads[iTh], &status);
      if (rc){
	std::cout << "Error:unable to join," << rc << std::endl;
	exit(-1);
      }
    }

    for (int k=0; k<fNBins; k++) {
      fHCompSig[i]->SetBinContent(k+1,fCompSigArray[k]);
      fHCompSig[i]->SetBinError(k+1,0);
    }

    fHCompTiming->Fill(stwch.RealTime(),1);

    if (fDoFFT) DoRooFFT(i);
  }

  return 0;
}

int PFScan::DoRooFFT(int iStep)
{
  TStopwatch stwch;
  TH1F hNonZeroCompSig("hNonZeroCompSig",
		       "hNonZeroCompSig",
		       fNBins,
		       //		       -fNBinsPerPeriod*2,
		       0,
		       fNBins);
		       //-fNBinsPerPeriod*2); //histogram with non-empty bins
  int iEntr=0;
  for (int i=1; i<fNBins+1; i++){
    //   std::cout<<"run on bin "<<i<<std::endl;
    float bico=fHCompSig[iStep]->GetBinContent(i);
    if (bico>0){
      iEntr++;
      hNonZeroCompSig.SetBinContent(i,bico);
    }
    else  hNonZeroCompSig.SetBinContent(i,fNFreq);
  }
  
  //  std::cout<<"filled histo"<<std::endl;

  TH1* hm=new TH1F("hm","hm",fNBins,0,fNBins);
  TVirtualFFT::SetTransform(0);
  hm = hNonZeroCompSig.FFT(hm, "MAG");
  fHCompSig_FFTImage[iStep]->Add(hm,1);
  hm->Delete();

  fHFFTTiming->Fill(stwch.RealTime(),1);
  return 0;
  /*
  for (int j=1; j<=fHCompSig_rooFFT[iStep]->GetNbinsX(); j++){
    hFFTResScale[iStep]->SetBinContent(j,hFFTRes[iStep]->GetBinContent(j));
    hFFTResPeriod[iStep]->SetBinContent(hFFTRes[iStep]->GetNbinsX()-j,hFFTRes[iStep]->GetBinContent(j));
  } 
  */
}
