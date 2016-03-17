#include "DMScan.h"
#include "DMAnalysis.h"
#include "thread_iface.h"
#include "TMath.h"
#include "Math/DistFunc.h"

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
      if (strcmp("-f",argv[1])==0) std::cout<<"missing file name"<<std::endl;
      std::cout<<argv[1]<<": wrong parameter option"<<std::endl;
      return 1;
    }
  
  if (argc==1)
    {
      configName="config/mainConfig.cff";
    }
  //  std::cout<<"input: "<<argc<<std::endl;
  //  std::cout<<argv[0]<<std::endl;
  //  std::cout<<argv[1]<<std::endl;
  int nThreads=1;
  int nFiles=0;
  int useGPU=0;
  std::string dirWithData=".";
  
  TStopwatch stwch;
  stwch.Start();
  
  TApplication theApp("App", &argc, argv);

  std::cout<<"read config: "<<configName<<std::endl;

  std::ifstream flist;
  flist.open(configName.c_str());
  //  flist.open("config/mainConfig.cff")
  if (!flist.is_open())
    {
      std::cout<<configName<<": failed to open config file"<<std::endl;
      return 1;
    }

  std::string fname;
  char tmp[100];
  
  std::string confParam;
  std::string scanConf;
  std::string anConf;
  std::string scanOutputFname;
  std::string analysisOutputFname;

  int startFileNumber;
  int rebinFactor, fitWindow, nPointsToScan;
  float DM0, scanStep;
  flist>>confParam;
  while(confParam!="runs:")
    {		 
      //      std::cout<<"ddd   "<<nThreads<<std::endl;
      if (confParam=="numberOfCPUThreads") flist>>nThreads;
      else if (confParam=="useGPU") flist>>useGPU;
      else if (confParam=="firstFileNumber") flist>>startFileNumber;
      else if (confParam=="nFilesToProcess") flist>>nFiles;
      else if (confParam=="dataLocation") flist>>dirWithData;
      else if (confParam=="runConf") flist>>scanConf;
      else if (confParam=="DM0") flist>>DM0;
      else if (confParam=="scanStep") flist>>scanStep;
      else if (confParam=="scanOutputFile") flist>>scanOutputFname;
      else if (confParam=="analysisOutputFile") flist>>analysisOutputFname;
      else if (confParam=="rebinFactor") flist>>rebinFactor;
      else if (confParam=="fitWindow") flist>>fitWindow;
      else if (confParam=="nPointsToScan") flist>>nPointsToScan;
      //      else 
      //     	   {
      //		std::cout<<"run list: unrecognized parameter, break"<<std::endl;
      //	        break;
      //	   }
       flist.getline(tmp,100,'\n');
       flist>>confParam;
     }
  flist.getline(tmp,100,'\n');

  std::cout<<"create dmscan"<<std::endl;

  std::vector<std::string> runID;

  int r=0;
  std::string rID;
  while (flist>>rID)
    //for (int r=0; r<nFiles; r++)
    {
      r++;
      flist.getline(tmp,100,'\n');
      if (r<startFileNumber) continue;
      if (r>=startFileNumber+nFiles) break;
      //flist>>rID
      runID.push_back(rID);
    }
  
  std::string ff=scanConf;
  DMScan dmscan(ff, nThreads, nPointsToScan, DM0, scanStep);

  std::cout<<"create dmanalysis"<<std::endl;

  TFile* outf=new TFile(analysisOutputFname.c_str(),"RECREATE");

  TH1F* sigVSdm=new TH1F("sigVSdm","sigVSdm",dmscan.nPointsToScan,dmscan.getScanStart(),dmscan.getScanEnd());
  TH1F* minPosVSdm=new TH1F("minPosVSdm","minPosVSdm",dmscan.nPointsToScan,dmscan.getScanStart(),dmscan.getScanEnd()); 

  TH1F* hFFTPeriod_SNRvsDM=new TH1F("hFFTPeriod_SNRvsDM","hFFTPeriod_SNRvsDM",nPointsToScan, DM0, DM0+nPointsToScan*scanStep);
  
  TFile* anreadf=new TFile(scanOutputFname.c_str());

  //  for (int i=0; i<nFiles; i++)
  //   {

  int nper=dmscan.getNPer();
  int bpp=dmscan.getBPP();
  float dm0=dmscan.getScanStart();
  float dmL=dmscan.getScanEnd();
  float tau=dmscan.getTau();
  DMAnalysis dma(anreadf, runID, nper, bpp, nPointsToScan, dm0, dmL, tau);
  
  std::cout<<"number of scan points to process:: "<<dma.getNScanPoints()<<"   binsPerPeriod: "<<bpp<<std::endl;

  float fftPeak=0;
  
  for (int j=0; j<dma.getNScanPoints(); j++)
    {
      std::cout<<"DM: "<<dmscan.getDM(j);
      //fold by period
      dma.foldByPeriod(j);
      
      //run significance scan
      //std::cout<<"run significance scan"<<std::endl;
      //      std::vector<float> sign=dma.localSgfScan(j,fitWindow,rebinFactor);
      //std::cout<<"search for giant peaks"<<std::endl;
      std::vector<float> sign=dma.searchForPeaks(j,fitWindow,rebinFactor);
      std::cout<<"   max significance:   "<<sign[0]<<"      "<<sign[1]<<std::endl;
      //      std::cout<<"do fourier"<<std::endl;
      //      int fft = dma.doFourier(j);
      // std::cout<<"done"<<std::endl;
      dma.fImageMerge(j);
      //fftPeak=dma.searchPeaksInFFT(j, 1.3, 0.2, 0.5);   //for 0301
      //fftPeak=dma.searchPeaksInFFT(j, 0.44, 0.2, 0.2);   //for 0357
      fftPeak=dma.searchPeaksInFFT(j, 0.319, 0.01, 0.1);   //for 2055
      hFFTPeriod_SNRvsDM->Fill(dmscan.getDM(j)+0.01, fftPeak);
      sigVSdm->Fill(dmscan.getDM(j)+0.01, sign[0]); 
      //-ROOT::Math::normal_quantile(sign[0],1));
      minPosVSdm->Fill(dmscan.getDM(j), sign[1]);
    }

  outf->mkdir("foldedHistograms");
  dma.write(outf,"foldedHistograms");

  std::cout<<"timer: "<<stwch.CpuTime()<<"   "<<stwch.RealTime()<<std::endl;
  
  /*
  TH1* ftest=0;
  TF1 *fsin = new TF1("fsin", "sin(x)",0, 4*TMath::Pi());
  //+sin(2*x)+sin(0.5*x)+1", 0, 4*TMath::Pi());
  Int_t n=25;
  TH1D *hsin = new TH1D("hsin", "hsin", n+1, 0, 4*TMath::Pi());
  Double_t x;
  
  //Fill the histogram with function values
  for (Int_t i=0; i<=n; i++){
    x = (Double_t(i)/n)*(4*TMath::Pi());
    hsin->SetBinContent(i+1, fsin->Eval(x));
  
  }
  */
  /*  
  TVirtualFFT::SetTransform(0);
  ftest = hsin->FFT(ftest, "MAG");
  int ftestNB=ftest->GetNbinsX();
  double ftestSTART=ftest->GetBinLowEdge(1);
  double ftestFINISH=ftest->GetBinLowEdge(ftestNB)+ftest->GetBinWidth(ftestNB);
  std::cout<<"axis limits: "<<ftestSTART/(4*TMath::Pi())<<"  "<<ftestFINISH/(4*TMath::Pi())<<"  "<<ftestFINISH<<std::endl;
  TH1D* ftestScale=new TH1D("ftestScale","ftestScale",ftestNB,ftestSTART/(4*TMath::Pi()),ftestFINISH/(4*TMath::Pi()));
  for (int i=1; i<ftestNB; i++)
    {
      ftestScale->SetBinContent(i,ftest->GetBinContent(i));
    }
  */
  outf->cd();
  //  std::cout<<"fff"<<std::endl;
  //  hsin->Write();
  //  ftest->Write();
  //  ftestScale->Write();
  hFFTPeriod_SNRvsDM->Write();
  sigVSdm->Write();
  minPosVSdm->Write();
  //  std::cout<<"sss"<<std::endl;
  dma.hFoldedProfile[0]->Write();
  dma.hBackground[0]->Write();
  //dma.hFFTPeriod_SNRvsDM->Write();
  outf->Close();
  // std::cout<<"aaa"<<std::endl;
  return 1;
}
