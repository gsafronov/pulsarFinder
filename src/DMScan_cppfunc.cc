#include "DMScan.h"
#include "TMath.h"
#include "Math/DistFunc.h"

//UNDERSTAND NORMALIZATION OF BANDS AND ERRORS IN COMPENSATED AND FOLDED PROFILE

DMScan::DMScan(std::string configFile, int nThr, int nPo, float sdm, float step)
{
  nThreads=nThr;
  nPointsToScan=nPo;
  DM0=sdm;
  scanStep=step;  
  std::ifstream cff;
  cff.open(configFile.c_str());
  std::string param;
  char tmp[100];
  while (!(cff.eof()))
    {
      cff>>param;
      if (cff.eof()) break;
      if (param=="tau") cff>>tau;
      if (param=="period") cff>>period;
      if (param=="l511") cff>>l511;
      if (param=="dL") cff>>dL;
      if (param=="nFrequencies") cff>>nFreq;
      if (param=="nPeriods") cff>>nPeriodsGlobal;
      if (param=="nBinsPerPeriod") cff>>nBinsPerPeriod;
      cff.getline(tmp,100,'\n');
      //std::cout<<"dd "<<cff.end<<std::endl;
    }
  nBinsGlobal=nPeriodsGlobal*nBinsPerPeriod;
  gpuCodeTiming=new TProfile("gpuCodeTiming","gpuCodeTiming",20,0,20,0,10000); 

}

DMScan::~DMScan()
{
}
/*
int DMScan::initScan(std::string rootfile)
{
  inputFile=new TFile(rootfile.c_str());

  
  int nPer, nBPer;
  
  TTree* tr=(TTree*)inputFile->Get("trInfo");
  tr->SetBranchAddress("nPeriods",&nPer);
  tr->SetBranchAddress("nBinsPerPeriod",&nBPer);
  //  tr->SetBranchAddress("year",&yy);
  //  tr->SetBranchAddress("month",&mm);
  //  tr->SetBranchAddress("day",&dd);
  //  tr->SetBranchAddress("hour",&hr);
  //  tr->SetBranchAddress("minute",&min);
  //  tr->SetBranchAddress("second",&sec);
  //  tr->SetBranchAddress("fsec",&fsec);
  tr->GetEntry(0);
  
  nPeriods=nPer;
  nBinsPerPeriod=nBPer;
  nBins=nBinsPerPeriod*nPeriods;
  

  char tmp[100];
  for (int y=0; y<512; y++)
    {
      sprintf(tmp,"sigTimeProfile_freqID_%d",y);
      sigTimeProfile.push_back((TH1F*)inputFile->Get(tmp));
	//sigTimeProfile[y]=(TH1F*)inputFile->Get(tmp);
      means.push_back(0);
      nBins=sigTimeProfile[y]->GetNbinsX();

      for (int i=0; i<nBins; i++)
	{
	  means[y]+=sigTimeProfile[y]->GetBinContent(i);
	}
      means[y]=means[y]/nBins;
      //      std::cout<<"means:   "<<means[y]<<std::endl;
      if (means[y]!=0) sigTimeProfile[y]->Scale(pow(means[y],-1));
    }

  nPeriods=nBins/nBinsPerPeriod;

  std::cout<<"number of periods: "<<nPeriods<<"    number of bins: "<<nBins<<"   nBinsPerPeriod: "<<nBinsPerPeriod<<std::endl;
  
  //estimate signal value;
  //TH1F* signalEst=(TH1F*)inputFile->Get("sigTimeProfile_freqID_0");

  for (int i=0; i<nPointsToScan; i++)
    {
      sprintf(tmp,"sumFreq_%d",i);
      sumFreq.push_back(new TH1F(tmp,tmp,nBins,0,nBins));
    }

  //read file contents into the memory:
  size_t size_input=(nBins*512)*sizeof(float);
  sigArray=(float*)malloc(size_input);
  for (int i = 0; i < 512; ++i)
    {
      for (int j=0; j < nBins; ++j)
	{
	  sigArray[i*nBins+j] = sigTimeProfile[i]->GetBinContent(j+1);
	}
    }

  //Device input vector:
  d_sigArray = NULL;
  
  
  //Allocate the host output vector C
  size_t size_output=nBins*sizeof(float);
  sigSum = (float *)malloc(size_output);
  for (int k=0; k<nBins; k++)
    {
      sigSum[k]=0;
    }
  
  //Device input vector:
  d_sigSum = NULL;
  
  
  for (int i=0; i<nPointsToScan; i++)
    {
      sprintf(tmp,"foldedProfile_%d",i);
      foldedProfile.push_back(new TH1F(tmp,tmp,nBinsPerPeriod,0,nBinsPerPeriod));
    }

  for (int i=0; i<nPointsToScan; i++)
    {
      sprintf(tmp,"background_%d",i);
      background.push_back(new TH1F(tmp,tmp,1000,0,10));
    }
  
  return 0;
}
*/
int DMScan::sumFrequencies_CPU(int iThread, int iStep)
{
  //  std::cout<<nBins<<"  "<<nBins/nBinsPerPeriod<<std::endl;
  float DM=getDM(iStep);
  int startPeriod=iThread*floor((nPeriods+10)/nThreads); //add 10 for safety
  //  std::cout<<"start: "<<startPeriod<<std::endl;
  int endPeriod=(iThread+1)*floor((nPeriods+10)/nThreads);
  for (int i=startPeriod*nBinsPerPeriod; i<endPeriod*nBinsPerPeriod; i++)
    {
      if (i<=nBinsPerPeriod||i>nBins-nBinsPerPeriod) continue;
      float bico=0;
      for (int y=0; y<512; y++) //debug 
	//512; y++)
	{
	  //take sigTimeProfile[511-y] as 511-th is shorter wavelength
	  //calculate delay wrt 511th for particular freq[511-y]
	  float dT=4.6*(-l511*l511+pow(l511+y*dL,2))*DM*0.001; //covert to ms
	  //calculate residual difference to nearest positive side pulse
	  float dTnearest=dT-period*floor(dT*pow(period,-1));
	  //move frequency band by -dTnearest bins, add lower bins to "upper side"
	  float delta=dTnearest*pow(tau,-1);
	  float bico1=sigTimeProfile[511-y]->GetBinContent(int(floor(i+delta))%nBins);
	  float bico2=sigTimeProfile[511-y]->GetBinContent(int((floor(i+delta)+1))%nBins);
	  float loFrac=1-((i+delta)-floor(i+delta));
	  float upFrac=1-loFrac;
	  if (floor(i+delta)+1<nBins) bico+=loFrac*bico1+upFrac*bico2;
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
  inputFile->Close();
  return 0;
}

int DMScan::rejectSpikes(float sigmaCut)
{
  TH1F sumSigRef("sumSigRef","sumSigRef",nBins,0,nBins);
  TH1F backgRef("backgRef","backgRef",5010,-10,5000);
  TH1F backgPerPeriod("backgPerPeriod","backgPerPeriod",1000,0,1000);

  for (int y=0; y<512; y++)
    {
      sumSigRef.Add(sigTimeProfile[y],1);
    }
  
  for (int y=0; y<nBins; y++)
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
  
  std::cout<<"removing spikes     "<<"reference mean: "<<refMean<<"     RMS: "<<refRMS<<"    number of bins: "<<nBins<<std::endl;
  
  TH1F backgAvg("backgAvg","backgAvg",1000,0,1000);
  
  for (int i=0; i<nPeriods; i++)
    {
      for (int j=0; j<nBinsPerPeriod; j++)
	{
	  backgPerPeriod.Fill(sumSigRef.GetBinContent(i*nBinsPerPeriod+j),1);
	}
      backgAvg.Fill(backgPerPeriod.GetMean(),1);
//      std::cout<<"period: "<<i<<"   bin: "<<i*nBinsPerPeriod<<"    average: "<<backgPerPeriod.GetMean()<<std::endl;
      backgPerPeriod.Reset();
    }
  
  std::cout<<"RMS of period average: "<<backgAvg.GetRMS()<<"     Mean of period average: "<<backgAvg.GetMean()<<std::endl;
  std::cout<<"if RMS > 50 do not take the run"<<std::endl;
    
  bool isBadRun=false;
  if (backgAvg.GetRMS()>50) isBadRun=true; 
  if (backgAvg.GetMean()>612||backgAvg.GetMean()<412) isBadRun=true;

  for (int q=1; q<nBins+1; q++)
    {
      float bico=sumSigRef.GetBinContent(q);
      if (bico>refMean+sigmaCut*refRMS) std::cout<<"SPIKE "<<"     value: "<<bico<<"  bin: "<<q<<std::endl;
      for (int y=0; y<512; y++)
	{
	  if (isBadRun)
	    {
	      refMean=512;
	      sigTimeProfile[y]->SetBinContent(q,refMean/512);
	      sigArray[y*nBins+q-1]=refMean/512;
	    }
	  if (bico>refMean+sigmaCut*refRMS)
	    {
	      sigTimeProfile[y]->SetBinContent(q,refMean/512);
	      sigArray[y*nBins+q-1]=refMean/512;
	    }
	}
    }
  
  /*
  for (int q=0; q<nBins; q++)
    {
      float bicoMid=sumSigRef.GetBinContent(q);
      float bicoFw;
      if (q<nBins-2) bicoFw=sumSigRef.GetBinContent(q+2);
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


