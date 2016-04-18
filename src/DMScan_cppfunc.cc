#include "DMScan.h"
#include "TMath.h"
#include "Math/DistFunc.h"

//UNDERSTAND NORMALIZATION OF BANDS AND ERRORS IN COMPENSATED AND FOLDED PROFILE

DMScan::DMScan(std::string configFile, int nThr, int nPo, float sdm, float step)
{
  fNThreads=nThr;
  fNPointsToScan=nPo;
  fDM0=sdm;
  fScanStep=step;  
  gpuCodeTiming=new TProfile("hGpuCodeTiming","hGpuCodeTiming",20,0,20,0,10000); 
}

DMScan::~DMScan()
{
}

int DMScan::sumFrequencies_CPU(int iThread, int iStep)
{
  //  std::cout<<fNBins<<"  "<<fNBins/fNBinsPerPeriod<<std::endl;
  float DM=getDM(iStep);
  int startPeriod=iThread*floor((fNPeriods+10)/fNThreads); //add 10 for safety
  //  std::cout<<"start: "<<startPeriod<<std::endl;
  int endPeriod=(iThread+1)*floor((fNPeriods+10)/fNThreads);
  for (int i=startPeriod*fNBinsPerPeriod; i<endPeriod*fNBinsPerPeriod; i++)
    {
      if (i<=fNBinsPerPeriod||i>fNBins-fNBinsPerPeriod) continue;
      float bico=0;
      for (int y=0; y<512; y++) //debug 
	//512; y++)
	{
	  //take sigTimeProfile[511-y] as 511-th is shorter wavelength
	  //calculate delay wrt 511th for particular freq[511-y]
	  float dT=4.6*(-fL511*fL511+pow(fL511+y*fDL,2))*DM*0.001; //covert to ms
	  //calculate residual difference to nearest positive side pulse
	  float dTnearest=dT-fPeriod*floor(dT*pow(fPeriod,-1));
	  //move frequency band by -dTnearest bins, add lower bins to "upper side"
	  float delta=dTnearest*pow(fTau,-1);
	  float bico1=sigTimeProfile[511-y]->GetBinContent(int(floor(i+delta))%fNBins);
	  float bico2=sigTimeProfile[511-y]->GetBinContent(int((floor(i+delta)+1))%fNBins);
	  float lowerBinFrac=1-((i+delta)-floor(i+delta));
	  float upperBinFrac=1-lowerBinFrac;
	  if (floor(i+delta)+1<fNBins) bico+=lowerBinFrac*bico1+upperBinFrac*bico2;
	}
      sumFreq[iStep]->Fill(i-1,bico);
    }
    return 0;
}

int DMScan::write(TFile* outFile, int iStep, std::string dir)
{
  outFile->cd(dir.c_str());
  sumFreq[iStep]->Write();
  sumFreq_fImage[iStep]->Write();
  return 0;
}

int DMScan::reset()
{
  //here is the segmentation violation
  sumFreq.clear();
  sumFreq_fImage.clear();
  sigTimeProfile.clear();
  fInputFile->Close();
  return 0;
}

int DMScan::rejectSpikes(float sigmaCut)
{
  TH1F sumSigRef("sumSigRef","sumSigRef",fNBins,0,fNBins);
  TH1F backgRef("backgRef","backgRef",5010,-10,5000);
  TH1F backgPerPeriod("backgPerPeriod","backgPerPeriod",1000,0,1000);

  for (int y=0; y<512; y++)
    {
      sumSigRef.Add(sigTimeProfile[y],1);
    }
  
  for (int y=0; y<fNBins; y++)
    {
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
  
  for (int i=0; i<fNPeriods; i++)
    {
      for (int j=0; j<fNBinsPerPeriod; j++)
	{
	  backgPerPeriod.Fill(sumSigRef.GetBinContent(i*fNBinsPerPeriod+j),1);
	}
      backgAvg.Fill(backgPerPeriod.GetMean(),1);
//      std::cout<<"period: "<<i<<"   bin: "<<i*fNBinsPerPeriod<<"    average: "<<backgPerPeriod.GetMean()<<std::endl;
      backgPerPeriod.Reset();
    }
  
  std::cout<<"RMS of period average: "<<backgAvg.GetRMS()<<"     Mean of period average: "<<backgAvg.GetMean()<<std::endl;
  std::cout<<"if RMS > 50 do not take the run"<<std::endl;
    
  bool isBadRun=false;
  if (backgAvg.GetRMS()>50) isBadRun=true; 
  if (backgAvg.GetMean()>612||backgAvg.GetMean()<412) isBadRun=true;

  for (int q=1; q<fNBins+1; q++)
    {
      float bico=sumSigRef.GetBinContent(q);
      if (bico>refMean+sigmaCut*refRMS) std::cout<<"SPIKE "<<"     value: "<<bico<<"  bin: "<<q<<std::endl;
      for (int y=0; y<512; y++)
	{
	  if (isBadRun)
	    {
	      refMean=512;
	      sigTimeProfile[y]->SetBinContent(q,refMean/512);
	      fSigArray[y*fNBins+q-1]=refMean/512;
	    }
	  if (bico>refMean+sigmaCut*refRMS)
	    {
	      sigTimeProfile[y]->SetBinContent(q,refMean/512);
	      fSigArray[y*fNBins+q-1]=refMean/512;
	    }
	}
    }
  
  /*
  for (int q=0; q<fNBins; q++)
    {
      float bicoMid=sumSigRef.GetBinContent(q);
      float bicoFw;
      if (q<fNBins-2) bicoFw=sumSigRef.GetBinContent(q+2);
      else bicoFw=refMean;
      float bicoBckw;
      if (q>2) bicoBckw=sumSigRef.GetBinContent(q-2);
      else bicoBckw=refMean;
      if (bicoMid>refMean+sigmaCut*refRMS&&bicoFw<refMean+sigmaCut*refRMS&&bicoBckw<refMean+sigmaCut*refRMS)
	{
	  for (int y=0; y<512; y++)
	    {
	      sigTimeProfile[y]->SetBinContent(q,refMean/512);
	    }
	}
    }
  */
  //  sigTimeProfile[y]->Scale(pow(norm[y],1));
  return 0;
}


int DMScan::initScan(std::string rootfile)
{
  fInputFile=new TFile(rootfile.c_str());

  if (fInputFile->IsZombie()) {
    std::cout<<"DMScan::initScan root file "<<rootfile.c_str()<<" not found"<<std::endl;
    return 1;
  }

  //read Run Header
  std::cout<<"blah"<<std::endl;
  TTree* RunHeader=(TTree*)fInputFile->Get("RunHeader");
  RunHeader->SetBranchAddress("fNPeriods",&fNPeriods);
  RunHeader->SetBranchAddress("fNBinsPerPeriod",&fNBinsPerPeriod);

  RunHeader->SetBranchAddress("fTau",&fTau);
  RunHeader->SetBranchAddress("fPeriod",&fPeriod);

  RunHeader->SetBranchAddress("fYear",&fYear);
  RunHeader->SetBranchAddress("fMonth",&fMonth);
  RunHeader->SetBranchAddress("fDay",&fDay);
  RunHeader->SetBranchAddress("fHour",&fHour);
  RunHeader->SetBranchAddress("fMinute",&fMinute);
  RunHeader->SetBranchAddress("fSecond",&fSecond);
  
  RunHeader->SetBranchAddress("fFreq0",&fFreq0);
  RunHeader->SetBranchAddress("fFreq511",&fFreq511);

  RunHeader->GetEntry(0);

  RunHeader->Delete();

  fL511=3e10*pow(fFreq511*1e6,-1);
  fL0=3e10*pow(fFreq0*1e6,-1);
  fDL=(fL0-fL511)*pow(512,-1);
  
  fNBins=fNPeriods*fNBinsPerPeriod;

  fNFreq=512;

  //std::cout<<fNFreq<<"  "<<fNBins<<"  "<<fNBinsPerPeriod<<"  "<<fL511<<"  "<<fDL<<"  "<<fTau<<"   "<<fPeriod<<std::endl;

  char tmp[100];
  for (int y=0; y<512; y++){
    sprintf(tmp,"sigTimeProfile_freqID_%d",y);
    sigTimeProfile.push_back((TH1F*)fInputFile->Get(tmp));
    //sigTimeProfile[y]=(TH1F*)inputFile->Get(tmp);
    fMeans.push_back(0);
    //    fNBins=sigTimeProfile[y]->GetNbinsX();
    
    for (int i=0; i<fNBins; i++){
      fMeans[y]+=sigTimeProfile[y]->GetBinContent(i);
    }
    fMeans[y]=fMeans[y]/fNBins;
    //      std::cout<<"means:   "<<means[y]<<std::endl;
    if (fMeans[y]!=0) sigTimeProfile[y]->Scale(pow(fMeans[y],-1));
  }
  
  for (int i=0; i<fNPointsToScan; i++){
    sprintf(tmp,"sumFreq_%d",i);
    sumFreq.push_back(new TH1F(tmp,tmp,fNBins,0,fNBins));
  }
  
  for (int i=0; i<fNPointsToScan; i++){
    sprintf(tmp,"sumFreq_fImage_%d",i);
    sumFreq_fImage.push_back(new TH1F(tmp,tmp,fNBins,0,fNBins));
  }
  
  //read file contents into the memory:
  size_t size_input=(fNBins*512)*sizeof(float);
  fSigArray=(float*)malloc(size_input);
  for (int i = 0; i < 512; ++i){
    for (int j=0; j < fNBins; ++j){
	  fSigArray[i*fNBins+j] = sigTimeProfile[i]->GetBinContent(j+1);
	}
    }
  
  //Device input vector:
  fDev_SigArray = NULL;
 
  return 0;
}
