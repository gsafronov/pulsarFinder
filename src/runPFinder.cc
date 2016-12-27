#include "PFScan.h"
#include "PFAnalysis.h"
//#include "TMath.h"
#include "Math/DistFunc.h"
#include "TApplication.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>

int convertStringParam(std::string IN, bool* OUT)
{
  int returnValue=0;
  if (IN=="yes") *OUT=true;
  else if (IN=="no") *OUT=false;
  else returnValue=1;
  return returnValue;
}

int main(int argc, char *argv[])
{
  //config file name
  std::string configName="config/mainConfig.cff";
  //read input with options
  if (argc > 3) 
    {
      std::cout<<"too many parameters"<<std::endl;
      return 1;
    }
  if (argc==3)
    {
      if (strcmp("-f",argv[1])==0) configName=std::string(argv[2]);
      else 
	{
	  std::cout<<argv[1]<<"can use only -f option"<<std::endl;
	  return 1;
	}
    }
  if (argc==2)
    {
      std::cout<<argv[1]<<": do not understand this parameter"<<std::endl;
      return 1;
    }
  
  if (argc==1)
    {
      configName="config/mainConfig.cff";
    }
  
  int nCPUThreads=1;
  int nGPUThreadsPerBlock=1;
  bool useGPU=false;
  int startFileNumber=0;
  int nFiles=0;
  int nPointsToScan=1;
  float DM0=0;
  float scanStep=1;
  int rebinFactor=1;
  bool removeSpikes=true;

  bool doDMScan=false;
  bool doFFT=false;
  bool doAnalysis;
  bool saveScanOutput=false;

  //string inputs
  std::string sRemoveSpikes="yes";
  std::string sUseGPU="no";
  std::string sDoDMScan="no";
  std::string sDoFFT="no";
  std::string sSaveScanOutput="no";
  std::string sDoAnalysis="no";

  std::string anOutputFileName="out/anOutput.root";
  std::string scanOutputLocation=".";
  std::string dirWithData=".";

  TStopwatch stwch;
    
  TApplication theApp("App", &argc, argv);

  std::cout<<"reading config"<<std::endl;

  char tmp[200];

  std::ifstream flist;
  flist.open(configName.c_str());

  if (flist.is_open()){
    std::string confParam="default";
    flist>>confParam;
    while(confParam!="runs:"){		 
      //std::cout<<confParam<<std::endl;
      if (confParam=="CPU_NumberOfThreads") flist>>nCPUThreads;
      else if (confParam=="GPU_NumberOfThreadsPerBlock") flist>>nGPUThreadsPerBlock;
      else if (confParam=="useGPU") flist>>sUseGPU;
      else if (confParam=="doDMScan") flist>>sDoDMScan;
      else if (confParam=="doFFT") flist>>sDoFFT;
      else if (confParam=="saveScanOutput") flist>>sSaveScanOutput;
      else if (confParam=="doAnalysis") flist>>sDoAnalysis;
      else if (confParam=="firstFileNumber") flist>>startFileNumber;
      else if (confParam=="nFilesToProcess") flist>>nFiles;
      else if (confParam=="inputDataLocation") flist>>dirWithData;
      else if (confParam=="scanOutputLocation") flist>>scanOutputLocation;
      else if (confParam=="anOutputFileName") flist>>anOutputFileName;
      else if (confParam=="DM0") flist>>DM0;
      else if (confParam=="scanStep") flist>>scanStep;
      else if (confParam=="rebinFactor") flist>>rebinFactor;
      else if (confParam=="removeSpikes") flist>>sRemoveSpikes;
      //      else if (confParam=="scanOutputFile") flist>>scanOutputFname;
      //else if (confParam=="rebinFactor") flist>>rebinFactor;
      //else if (confParam=="fitWindow") flist>>fitWindow;
      else if (confParam=="nPointsToScan") flist>>nPointsToScan;
      //std::cout<<"ddd"<<std::endl;
      /*
      else 
	{
	  std::cout<<"run list: unrecognized parameter, break"<<std::endl;
	  break;
	}
      */
      flist.getline(tmp,200,'\n');
      flist>>confParam;
    }
  }
  else {
    std::cout<<"config file "<<configName<<" was not found"<<std::endl;
    return 1;
  }

  int err=0;
  err+=convertStringParam(sUseGPU, &useGPU);
  err+=convertStringParam(sDoDMScan, &doDMScan);
  err+=convertStringParam(sDoFFT, &doFFT);
  err+=convertStringParam(sSaveScanOutput, &saveScanOutput);
  err+=convertStringParam(sDoAnalysis, &doAnalysis);
  err+=convertStringParam(sRemoveSpikes, &removeSpikes);

  if (err>0) {
    std::cout<<"misprint in yes/no parameters"<<std::endl;
    return 1;
  }

  //loop over runs
  std::vector<std::string> runID;
  int r=0;
  std::string rID="001";
  
  while (flist>>rID){
    r++;
    flist.getline(tmp,100,'\n');
    if (r<startFileNumber) continue;
    if (r>=startFileNumber+nFiles) break;
    runID.push_back(rID);
  }

  std::cout<<"create PFScan"<<"    rebin: "<<rebinFactor<<std::endl;
  PFScan scan(nPointsToScan, DM0, scanStep, rebinFactor);
  std::cout<<"create PFAnalysis"<<std::endl;
  PFAnalysis analysis(anOutputFileName, &scan, doFFT);

  for (int i=0; i<runID.size(); i++){
    flist.getline(tmp,100,'\n');
    // if (r<startFileNumber) continue;
    // if (r>=startFileNumber+nFiles) break;
    // runID.push_back(rID);
    std::string fnameRAW=dirWithData+"/readRAW_"+runID[i]+".root";
    std::string fnameScanOut=scanOutputLocation+"/scanOut_"+runID[i]+".root";
    
    std::cout<<"Initalize run #"<<i+1<<"   "<<runID[i]<<";  "<<fnameRAW<<std::endl;
    //    std::cout<<"timer: "<<stwch.CpuTime()<<"   "<<stwch.RealTime()<<std::endl;
    if (doDMScan){
      scan.InitScan(fnameRAW, doFFT);
      if (removeSpikes) scan.CleanSignal(5, 50);

      if (useGPU) scan.DoScan_GPU(nGPUThreadsPerBlock);
      else scan.DoScan_CPU(nCPUThreads);
      
      if (saveScanOutput){
	scan.SaveOutput(fnameScanOut);
      }
      
      scan.CloseScan();
    }
    if (doAnalysis){
      analysis.AddRun(runID[i], fnameRAW, fnameScanOut);
    }
    
  }
  if (doAnalysis) {
    analysis.SumRuns();
  }
  std::cout<<"done"<<std::endl;
}


