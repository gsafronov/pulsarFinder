#include "DMScan.h"
#include "DMAnalysis.h"
#include "thread_iface.h"
#include "TMath.h"
#include "Math/DistFunc.h"

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

  int nThreads=1;
  int useGPU=0;
  int startFileNumber=0;
  int nFiles=0;
  int nPointsToScan=1;
  int DM0=0;
  int scanStep=1;
  int rebinFactor=1;
  int fitWindow=100;
  
  std::string dirWithData=".";

  TStopwatch stwch;
    
  TApplication theApp("App", &argc, argv);

  std::cout<<"read config"<<std::endl;

  std::string fname;
  char tmp[200];

  //std::cout<<"sdsdsd"<<std::endl;

  std::string scanConf="config/scanParam.cff";
  std::string scanOutputFname="scanOut/scanOutput.root";

  //std::cout<<"bb"<<std::endl;
  //set default values
    
  std::ifstream flist;
  flist.open(configName.c_str());

  if (flist.is_open()){
    std::string confParam="default";
    flist>>confParam;
    while(confParam!="runs:"){		 
      //std::cout<<confParam<<std::endl;
      if (confParam=="numberOfCPUThreads") flist>>nThreads;
      else if (confParam=="useGPU") flist>>useGPU;
      else if (confParam=="firstFileNumber") flist>>startFileNumber;
      else if (confParam=="nFilesToProcess") flist>>nFiles;
      else if (confParam=="dataLocation") flist>>dirWithData;
      else if (confParam=="runConf") flist>>scanConf;
      else if (confParam=="DM0") flist>>DM0;
      else if (confParam=="scanStep") flist>>scanStep;
      else if (confParam=="scanOutputFile") flist>>scanOutputFname;
      else if (confParam=="rebinFactor") flist>>rebinFactor;
      else if (confParam=="fitWindow") flist>>fitWindow;
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
  
  flist.getline(tmp,100,'\n');

  std::cout<<"create dmscan"<<std::endl;

  std::vector<std::string> runID;

  std::string ff=scanConf;
  
  if (useGPU==1) nThreads=1;
  DMScan dmscan(ff, nThreads, nPointsToScan, DM0, scanStep);

  std::cout<<"number of threads:  "<<nThreads<<std::endl; 
  pthread_t threads[nThreads];

  TFile scoutf(scanOutputFname.c_str(),"RECREATE");

  TH1F* timePerDM=new TH1F("timePerDM","timePerDM",10000,0,10);

  int r=0;
  std::string rID="001";
  while (flist>>rID)
    //for (int r=0; r<nFiles; r++)
    {
      r++;
      flist.getline(tmp,100,'\n');
      if (r<startFileNumber) continue;
      if (r>=startFileNumber+nFiles) break;
      //flist>>rID
      runID.push_back(rID);
      fname=dirWithData+"/readRAW_"+runID[r-startFileNumber]+".root";

      std::cout<<"NEXT RUN:   initalize scan of run #"<<r<<"   "<<rID<<";  "<<fname<<"; timer: "<<stwch.CpuTime()<<"   "<<stwch.RealTime()<<std::endl;
      stwch.Continue();

      scoutf.mkdir(runID[r-startFileNumber].c_str());
      
      int err = dmscan.initScan(fname.c_str());
      if (err!=0) {
	std::cout<<"Cannot initialize scan for run "<<runID[r-startFileNumber]<<". Skipping the run"<<std::endl;
	continue;
      }
      
      float sigmaCut=6;
      std::cout<<"reject DM0 spikes, cut at SNR "<<sigmaCut<<std::endl; 

      dmscan.rejectSpikes(sigmaCut);
      //      dmscan.sumFrequencies(0,0);

      //      dmscan.rejectSpikes(sigmaCut);

      if (useGPU==1) dmscan.loadDataToGPU();
      std::cout<<"start scan, timer: "<<stwch.CpuTime()<<"   "<<stwch.RealTime()<<std::endl;
      stwch.Continue();
      
      //      float prevTime=stwtch.RealTime();
      //      stwtch.Continue();
      for (int k=0; k<nPointsToScan; k++)
	{								
	  float prevTime=stwch.RealTime();
	  stwch.Continue();
	  if (k%10==0||k==nPointsToScan-1) std::cout<<"point #"<<k<<"    timer: "<<stwch.CpuTime()<<"   "<<stwch.RealTime()<<std::endl;
	  stwch.Continue();
	  void *status;
	  struct sumFreqID thrData[nThreads];
	  
	  for(int iTh=0; iTh < nThreads; iTh++ )
	    {
	      thrData[iTh].iThread=iTh;
	      thrData[iTh].iStep=k;
	      thrData[iTh].dmscan=&dmscan;
	      
	      int rc;
	      if (useGPU==0) rc = pthread_create(&threads[iTh], NULL, sumFreq_CPU_interface, (void *)&thrData[iTh]);
	      else rc = pthread_create(&threads[iTh], NULL, sumFreq_GPU_interface, (void *)&thrData[iTh]);
	      
	      if (rc){
		
		std::cout << "Error:unable to create thread," << rc << std::endl;
		exit(-1);
	      }
	    }
	  for ( int ii=0; ii < nThreads; ii++ ){	    
	    int rc = pthread_join(threads[ii], &status);
	    if (rc){
	      std::cout << "Error:unable to join," << rc << std::endl;
	      exit(-1);
	    }
	      //	      std::cout << "Main: completed thread id :" << ii ;
	      //std::cout << "  exiting with status :" << status << std::endl;
	    }

	  dmscan.write(&scoutf, k, runID[r-startFileNumber]);

	  float fullTime=stwch.RealTime()-prevTime;
	  stwch.Continue();
	  timePerDM->Fill(fullTime,1);
	}
      dmscan.closeGPU();
      dmscan.reset();
    }
  scoutf.cd();
  timePerDM->Write();
  dmscan.gpuCodeTiming->Write();
  scoutf.Close();

  //  timePerDM->Delete();
  
  return 0;
}
