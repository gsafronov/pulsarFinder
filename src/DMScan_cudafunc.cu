#include "include/DMScan.h"
#include <cuda_runtime.h>

//before class definition define CUDA kernel func
//see comments on actions in DMScan::sumFrequencies() function
//treat array as 1dim: in thread i sum_j(511-j+i), fullID - ID on a given freq band
__global__ void sumFreq_kernel(float *sigArray, float *sigSum, int nFreq, int nBins, int nBinsPerPeriod, float DM, float l511, float dL, float tau, float period)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  float bico=0;
  int fullID=i;
  for (int iFreq=0; iFreq<nFreq; iFreq++){
    float dT=4.6*(-l511*l511+(l511+iFreq*dL)*(l511+iFreq*dL))*DM*0.001;
    float dTnearest=dT-period*floor(dT/period);
    float delta=dTnearest/tau;
    float bico1=sigArray[((511-iFreq)*nBins+int(floor(fullID+delta)))%(nBins*nFreq)];
    float bico2=sigArray[((511-iFreq)*nBins+int((floor(fullID+delta)+1)))%(nBins*nFreq)];
    float loFrac=1-((fullID+delta)-floor(fullID+delta));
    float upFrac=1-loFrac;
    bico+=loFrac*bico1+upFrac*bico2;
  }
  sigSum[fullID]=bico;
}

int DMScan::loadDataToGPU()
{
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
  // Print the vector length to be used, and compute its size
  size_t size_input=(nBins*nFreq)*sizeof(float);
  size_t size_output=(nBins)*sizeof(float);	
  
  // Allocate the device input signal vector 
  err = cudaMalloc((void **)&d_sigArray, size_input);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to allocate the device input signal vector (error code %s)!\n", cudaGetErrorString(err));
    
    free(sigArray);
    
    return 0;
    //exit(EXIT_FAILURE);
  }
  
  // Copy the host input signal vector to the device signal vector in device memory
  err = cudaMemcpy(d_sigArray, sigArray, size_input, cudaMemcpyHostToDevice);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to copy input signal vector from host to device (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  free(sigArray);
  return 0;
}

int DMScan::sumFrequencies_GPU(int iStep)
{
  float DM=DM0+scanStep*iStep;
  
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
  size_t size_output=nBins*sizeof(float);

  //host compensated DM
  float* sigSum=(float *)malloc(size_output);
  
  // device compensated DM
  float* d_sigSum;
  err = cudaMalloc((void **)&d_sigSum, size_output);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to allocate compensated signal vector on the device (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  
  // Launch the Vector Add CUDA Kernel
  int threadsPerBlock = 32;
  int blocksPerGrid = floor(nBins/threadsPerBlock)+1;
  //  printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
  //       std::cout<<"parameters: "<<nFreq<<"  "<<nBins<<"  "<<DM<<"  "<<l511<<"  "<<dL<<"  "<<tau<<"  "<<period<<"   blocksPerGrid "<<blocksPerGrid<<std::endl;
  sumFreq_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_sigArray, d_sigSum, nFreq, nBins, nBinsPerPeriod, DM, l511, dL, tau, period);
  err = cudaGetLastError();
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to launch sumFreq kernel (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  
  // Copy the device result vector in device memory to the host result vector in host memory.
  //  printf("Copy output data from the CUDA device to the host memory:");
  //     std::cout<<"address: "<<sigSum<<"  "<<&sigSum[0]<<"  "<<d_sigSum<<"  "<<size_output<<"  "<<cudaMemcpyDeviceToHost<<std::endl;
  err = cudaMemcpy(sigSum, d_sigSum, size_output, cudaMemcpyDeviceToHost);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to copy compensated signal vector from device to host (error code: %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  
  for (int i=0; i<nBinsGlobal; i++){
    if (sigSum[i]==sigSum[i]&&i<nBins) sumFreq[iStep]->SetBinContent(i+1,sigSum[i]);
    else sumFreq[iStep]->SetBinContent(i+1,512);
  }
  
  //allocate host vector with non-zero elements and image
  //use common number of bins for all runs (nBinsGlobal)
  
  int signal_memSize=sizeof(float)*nBinsGlobal;
  float* h_sigSum_NZ = (float*)malloc(signal_memSize);
  //fill the vector
  for (int i=0; i<nBinsGlobal; i++){
    if (sigSum[i]!=0&&sigSum[i]==sigSum[i]&&i<nBins) h_sigSum_NZ[i]=sigSum[i];
    else h_sigSum_NZ[i]=512;
  }
  
  //allocate device vector with non-zero elements:
  int image_memSize = sizeof(float2) * (floor(nBinsGlobal/2)+1);
  float* d_sigSum_NZ;
  err=cudaMalloc((void **)&d_sigSum_NZ, image_memSize);
  if (err!=cudaSuccess){			
    std::cout<<"Failed to allocate device signal NZ vector"<<std::endl; 
    return 1;
  }
  
  //copy host to device
  err = cudaMemcpy(d_sigSum_NZ, h_sigSum_NZ, image_memSize, cudaMemcpyHostToDevice);

  if (err != cudaSuccess){
    fprintf(stderr, "Failed to copy signal NZ vector from host to device (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  /////////
  
  //allocate host and device images
  float2* h_fImage = (float2*)malloc(image_memSize);

  //Allocate device output vector for fImage
  float2* d_fImage;
  err = cudaMalloc((void **)&d_fImage, image_memSize);
  
  if (err!=cudaSuccess){			
    std::cout<<"Failed to allocate device fImage:"<<std::endl; 
    return 1;
  }
  
  // create CUFFT plan
  cufftHandle plan;
  if (cufftPlan1d(&plan, nBinsGlobal, CUFFT_R2C, 1)!=CUFFT_SUCCESS){
    std::cout<<"FFT plan creation failed"<<std::endl;
    return 1;
  }
  
  //DO THE FOURIER TRANSFORM OF d_sigSum;
  
  //run transform
  if (cufftExecR2C(plan,d_sigSum_NZ, d_fImage) != CUFFT_SUCCESS){
    std::cout<<"FFT: ExecR2C failed"<<std::endl;
    return 1;
  }
  //  else 
  
  //copy the output to the host
  
  //  int image_memSize = sizeof(float2) * nBins;
  err= cudaMemcpy(h_fImage, d_fImage, image_memSize, cudaMemcpyDeviceToHost);	
  
  if (err!=cudaSuccess){
    std::cout<<"Failed to copy fImage from device to host"<<std::endl;
    return 1;
  }   
  
  //fill the histogram
  for (int i=0; i<nBinsGlobal; i++){
    float amplitude;
    if (i<floor(nBinsGlobal/2)+1)	{  
      amplitude=sqrt(pow(h_fImage[i].x,2)+pow(h_fImage[i].y,2));
      if (amplitude==amplitude) sumFreq_fImage[iStep]->SetBinContent(i+1,amplitude);
      else std::cout<<"DMScan::sumFrequencies_GPU: nan in FFT bin "<<i<<", iStep: "<<iStep<<std::endl; 
	     //sumFreq_fImage[iStep]->SetBinContent(i+1,amplitude);
    }  
    else{
      int ind=nBinsGlobal-i;
      amplitude=sqrt(pow(h_fImage[ind].x,2)+pow(h_fImage[ind].y,2));
      if (amplitude==amplitude) sumFreq_fImage[iStep]->SetBinContent(i+1,amplitude);
      else std::cout<<"DMScan::sumFrequencies_GPU: nan in FFT bin "<<i<<", iStep: "<<iStep<<std::endl;
    }	
  }

  free(sigSum);
  free(h_sigSum_NZ);
  free(h_fImage);
    
  cudaFree(d_fImage);
  cudaFree(d_sigSum_NZ);
  cudaFree(d_sigSum);
  cufftDestroy(plan);
  
  return 0;
}

int DMScan::closeGPU()
{
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
  // Free device global memory
  err = cudaFree(d_sigArray);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to free device input signal vector (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  
  err = cudaDeviceReset();
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  return 0;
  //STEP 12
}

int DMScan::initScan(std::string rootfile)
{
  inputFile=new TFile(rootfile.c_str());

  if (inputFile->IsZombie()) {
    std::cout<<"DMScan::initScan root file "<<rootfile.c_str()<<" not found"<<std::endl;
    return 1;
  }
  
  char tmp[100];
  for (int y=0; y<512; y++){
    sprintf(tmp,"sigTimeProfile_freqID_%d",y);
    sigTimeProfile.push_back((TH1F*)inputFile->Get(tmp));
    //sigTimeProfile[y]=(TH1F*)inputFile->Get(tmp);
    means.push_back(0);
    nBins=sigTimeProfile[y]->GetNbinsX();
    
    for (int i=0; i<nBins; i++){
      means[y]+=sigTimeProfile[y]->GetBinContent(i);
    }
    means[y]=means[y]/nBins;
    //      std::cout<<"means:   "<<means[y]<<std::endl;
    if (means[y]!=0) sigTimeProfile[y]->Scale(pow(means[y],-1));
  }
  
  nPeriods=nBins/nBinsPerPeriod;
  
  for (int i=0; i<nPointsToScan; i++){
    sprintf(tmp,"sumFreq_%d",i);
    sumFreq.push_back(new TH1F(tmp,tmp,nBinsGlobal,0,nBinsGlobal));
  }
  
  for (int i=0; i<nPointsToScan; i++){
    sprintf(tmp,"sumFreq_fImage_%d",i);
    sumFreq_fImage.push_back(new TH1F(tmp,tmp,nBinsGlobal,0,nBinsGlobal));
  }
  
  //read file contents into the memory:
  size_t size_input=(nBins*512)*sizeof(float);
  sigArray=(float*)malloc(size_input);
  for (int i = 0; i < 512; ++i){
    for (int j=0; j < nBins; ++j){
	  sigArray[i*nBins+j] = sigTimeProfile[i]->GetBinContent(j+1);
	}
    }
  
  //Device input vector:
  d_sigArray = NULL;
  
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;

  if (err != cudaSuccess){
    fprintf(stderr, "Failed to allocate pinned host memory (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  } 
  
  return 0;
}

