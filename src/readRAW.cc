#include <iostream>
#include <fstream> 
#include <stdint.h>
#include <cmath>
#include <bitset>
#include "TH1F.h"
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"


float pulseToFloat(unsigned int pulse, float tau)
{
  float exp, spectr;
  spectr=(pulse&0xFFFFFF);
  exp=int(pulse&0x7F000000) >> 24;
  exp=exp-64-24;
  float ratio=tau/0.2048;
  spectr=spectr*std::pow(2,exp)/ratio;
  return spectr;
}

int readRAW(std::string runID, std::string rawdata_dir, std::string output_dir)
{
  std::string fname=output_dir+"/readRAW_"+runID+".root";
  TFile tfout(fname.c_str(),"RECREATE");

  fname = rawdata_dir+"/"+runID;
  std::ifstream data(fname.c_str(),std::ios::binary|std::ios::in);
  
  std::cout<<"READING RUN "<<runID<<"    path: "<<fname<<std::endl;

  //read header
  int length = 40; 
  char * buffer = new char [length];
  int sizeHeader;
  
  int iper, nBinsPerPeriod, nPeriods;
  int year, month, day, hour, minute, second;
  double fsec;
  
  //      std::cout<<"READING HEADER"<<std::endl;
  for (int k=0; k<13; k++)
    {
      // std::cout << "Reading " << length << " characters... "<<std::endl;
      // read data as a block:
      data.read(buffer,length);
      //  int size = data.tellg();
      sizeHeader = data.tellg();
      std::cout<<k<<"   "<<buffer<<std::endl;
      for (int q=0; q<length; q++)
	{
	  //	  std::cout<<q<<":"<<buffer[q]<<" ";
	}
      
      //covert buffer to a number
      int number;
      if (buffer[17]==' ') 
	{
	  number=1000*(int(buffer[13])-int('0'))+100*(int(buffer[14])-int('0'))+10*(int(buffer[15])-int('0'))+(int(buffer[16])-int('0'));
	}
      if (buffer[16]==' ') 
	{
	  number=100*(int(buffer[13])-int('0'))+10*(int(buffer[14])-int('0'))+(int(buffer[15])-int('0'));
	}
      if (buffer[15]==' ') 
	{
	  number=10*(int(buffer[13])-int('0'))+(int(buffer[14])-int('0'));
	}
      
      if (k==5) nPeriods=number;
      if (k==7) nBinsPerPeriod=number;
      
      //decode time
      if (k==12)
	{
	  day=(int(buffer[13])-int('0'))*10+(int(buffer[14])-int('0'));
	  month=(int(buffer[16])-int('0'))*10+(int(buffer[17])-int('0'));
	  year=(int(buffer[19])-int('0'))*10+(int(buffer[20])-int('0'));
	  hour=(int(buffer[22])-int('0'))*10+(int(buffer[23])-int('0'));
	  minute=(int(buffer[25])-int('0'))*10+(int(buffer[26])-int('0'));
	  second=(int(buffer[28])-int('0'))*10+(int(buffer[29])-int('0'));
	  fsec=(int(buffer[31])-int('0'))*0.1+(int(buffer[32])-int('0'))*0.01+(int(buffer[33])-int('0'))*0.001+(int(buffer[34])-int('0'))*0.0001+(int(buffer[35])-int('0'))*0.00001+(int(buffer[36])-int('0'))*0.000001+(int(buffer[37])-int('0'))*0.0000001;
	  //	      std::cout<<day<<" "<<month<<" "<<year<<" "<<hour<<" "<<minute<<" "<<second<<" "<<fsec<<std::endl;
	}
      
      
      //      std::cout<<"position: "<<size<<std::endl;
    }   

  float signal[512][nBinsPerPeriod];
  TTree trPeriod("trPeriod","trPeriod");
  trPeriod.Branch("iper",&iper,"iper/I");
  trPeriod.Branch("nBinsPerPeriod",&nBinsPerPeriod,"nBinsPerPeriod/I");
  trPeriod.Branch("nPeriods",&nPeriods,"nPeriods/I");
  trPeriod.Branch("signal",signal,"signal[512][nBinsPerPeriod]/F");
  
  TTree trInfo("trInfo","trInfo");
  trInfo.Branch("nPeriods",&nPeriods,"nPeriods/I");
  trInfo.Branch("nBinsPerPeriod",&nBinsPerPeriod,"nBinsPerPeriod/I");
  trInfo.Branch("year",&year,"year/I");
  trInfo.Branch("month",&month,"month/I");
  trInfo.Branch("day",&day,"day/I");
  trInfo.Branch("hour",&hour,"hour/I");
  trInfo.Branch("minute",&minute,"minute/I");
  trInfo.Branch("second",&second,"second/I");
  trInfo.Branch("fsec",&fsec,"fsec/D");
  
  //  trInfo->Branch(""
  trInfo.Fill();
  
  int nBins=nPeriods*nBinsPerPeriod;
  TH1F loFreq("loFreq","loFreq",nBins,0,nBins);
  TH1F hiFreq("hiFreq","hiFreq",nBins,0,nBins);
  TH1F foldPeriod("foldPeriod","foldPeriod",nBinsPerPeriod,0,nBinsPerPeriod);
  
  char tmp[100];
  TH1F sigTimeProfile[512];
  for (int k=0; k<512; k++)
    {
      sprintf(tmp,"sigTimeProfile_freqID_%d",k);
      sigTimeProfile[k]=TH1F(tmp,tmp,nBins,0,nBins);
    }
  
  std::cout<<"READING DATA    numPeriods: "<<nPeriods<<"   binsPerPeriod: "<<nBinsPerPeriod<<std::endl;
  // int lengthData= 8*sizeof(uint32_t);
  int lengthData = sizeof(uint32_t);
  // lengthData=4;
  char* fileContents;
  fileContents = new char[lengthData];
  //      std::cout<<"length: "<<lengthData<<std::endl;
  unsigned int number=0;
  //  char * bufferData = new char [lengthData];
  float tau=1.2288;
  int ipos;
  int iPoint=0;
  int iPointAbs=0;
  int iFreq=0;
  int iPeriod=0;
  while(data.good())
    {
      data.read((char *) &number,lengthData);
      int ipos = data.tellg();
      if (!data.good()) continue;
      //	  if ((ipos/4)%512==0&&((ipos/4)/512)%10000==0) std::cout<<(ipos/4)/512<<std::endl;
      float ampl=pulseToFloat(number,tau);
      int iFreq=(((ipos-sizeHeader)/lengthData-1)%512);
      signal[iFreq][iPoint]=ampl;
      if (iFreq!=513) 
	{
	  signal[iFreq][iPoint]=ampl;
	  sigTimeProfile[iFreq].Fill(iPointAbs,ampl);
	  if (iPointAbs>4000) foldPeriod.Fill(iPoint,ampl);
	}
      if (iFreq==511) 
	{
	  iPoint++;
	  if (iPoint%nBinsPerPeriod==0) 
	    {
	      iPoint=0;
	      iper=iPeriod;
	      trPeriod.Fill();
	      iPeriod++;
	    }	      
	  iPointAbs++;
	}	
      //break;
    }
  
  data.close();
  delete fileContents;
  delete buffer;
  
  tfout.cd();
  //      std::cout<<"blabla"<<std::endl;
  foldPeriod.Write();
  trPeriod.Write();
  for (int q=0; q<512; q++)
    {
      sigTimeProfile[q].Write(); 
    }
  trInfo.Write();
  tfout.Close();
  //char c = data.get();
  // while (data.good()) {
  //  std::cout << c;
  //  c = data.get();
  // }
  std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<std::endl;
  return 0;
}

struct readRunID{
  std::string runID;
  std::string rawdata_dir;
  std::string output_dir;
};


void* readRAW_CPU_interface(void *threadarg)
{
  struct readRunID *rrID;
  
  rrID = (struct readRunID *) threadarg;
  
  int rv = readRAW(rrID->runID, rrID->rawdata_dir, rrID->output_dir);

  pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
      //read input with options
  if (argc > 3) 
    {
      std::cout<<"too many params"<<std::endl;
      return 1;
    }
  
  std::string configName;

  if (argc==3)
    {
      if (strcmp("-f",argv[1])==0) configName=std::string(argv[2]);
      else 
	{
	  std::cout<<argv[1]<<": wrong parameter option"<<std::endl;
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
  
  std::ifstream flist;
  flist.open(configName.c_str());
  char tmp[100];
  int nFiles, nThreads;
  int startFileNumber, endFileNumber;
  std::string confParam;
  std::string rawdata_dir, output_dir;
  //  std::string fname;
  std::vector<std::string> runID;
  flist>>confParam;
  while(confParam!="runs:")
    {		 
      std::cout<<confParam<<std::endl;
      if (confParam=="inputDir") flist>>rawdata_dir;
      else if (confParam=="outputDir") flist>>output_dir;
      else if (confParam=="startFileNumber") flist>>startFileNumber;
      else if (confParam=="nRuns") flist>>nFiles;
      flist.getline(tmp,100,'\n');
      flist>>confParam;
      std::cout<<"linex"<<std::endl;
    }
  std::string rID;
  int runCounter=0;
  while (flist>>rID)
    {
      std::cout<<"bubu"<<std::endl;
      runCounter++;
      flist.getline(tmp,100,'\n');
      std::cout<<runCounter<<std::endl;
      if (runCounter<startFileNumber) continue;
      std::cout<<"rrr: "<<runCounter<<std::endl;
      if (runCounter>=startFileNumber+nFiles) break;
      runID.push_back(rID);
      std::cout<<"ddd"<<std::endl;
    }

  for (int iPack=0; iPack<floor(nFiles); iPack++)
    {
      readRAW(runID[iPack], rawdata_dir, output_dir);
    }
  return 0;
}
