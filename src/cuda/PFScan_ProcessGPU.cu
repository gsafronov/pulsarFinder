#include "PFScan.h"
#include "TF1.h"
#include "TStopwatch.h"
#include <cuda_runtime.h>
#include <cufft.h>

//before class definition define CUDA kernel func
//see comments on actions in DMScan::sumFrequencies() function
//treat array as 1dim: in thread i sum_j(511-j+i), fullID - ID on a given freq band
__global__ void sumFreq_kernel(float *sigArray, 
			       float *sigSum, 
			       int nFreq, 
			       int nBinsInput, 
			       int nBinsPerPeriod, 
			       float DM, 
			       float l511, 
			       float dL, 
			       float tau, 
			       float period)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  float bico=0;
  int fullID=i;
  for (int iFreq=0; iFreq<nFreq; iFreq++){
    float dT=4.6*(-l511*l511+(l511+iFreq*dL)*(l511+iFreq*dL))*DM*0.001;
    float dTnearest=dT-period*floor(dT/period);
    float delta=dTnearest/tau;
    //QUESTION: WHY %nBins*nFreq???, keep for safety omit few last periods later
    float bico1=sigArray[(((nFreq-1)-iFreq)*nBinsInput
			  +int(floor(fullID+delta)))%(nBinsInput*nFreq)];
    float bico2=sigArray[(((nFreq-1)-iFreq)*nBinsInput
			  +int((floor(fullID+delta)+1)))%(nBinsInput*nFreq)];
    float loFrac=1-((fullID+delta)-floor(fullID+delta));
    float upFrac=1-loFrac;
    bico+=loFrac*bico1+upFrac*bico2;
  }
  sigSum[fullID]=bico;
}

int PFScan::DoScan_GPU(int nThreadsPerBlock)
{
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
  // Print the vector length to be used, and compute its size
  size_t size_input=(fNBinsInput*fNFreq)*sizeof(float);
  //  size_t size_output=(fNBins)*sizeof(float);	
  
  // Allocate the device input signal vector 
  err = cudaMalloc((void **)&fDev_SigArray, size_input);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to allocate the device input signal vector (error code %s)!\n", cudaGetErrorString(err));
    
    //    free(fSigArray);
    
    return 0;
    //exit(EXIT_FAILURE);
  }
  
  // Copy the host input signal vector to the device signal vector in device memory
  err = cudaMemcpy(fDev_SigArray, fSigArray, size_input, cudaMemcpyHostToDevice);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to copy input signal vector from host to device (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  //  free(fSigArray);

  for (int i=0; i<fNScanPoints; i++){
    if ((i+1)%10==0||i+1==1) std::cout<<"PFScan::DoScan_GPU;  process point "<<i+1<<std::endl;
    DoCompensation_GPU(i, nThreadsPerBlock);
    if (fDoFFT) DoCuFFT(i);
  }

  CloseGPU();

  return 0; 
}


int PFScan::DoCompensation_GPU(int iStep, int nThreadsPerBlock)
{
  TStopwatch stwch;
  float DM=fDM0+fScanStep*iStep;
  
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
  size_t size_output=fNBinsInput*sizeof(float);

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
  int nBlocksPerGrid = floor(fNBinsInput/nThreadsPerBlock)+1;
  //  printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
  //  std::cout<<"parameters: "<<fNFreq<<"  "<<fNBins<<"  "<<DM<<"  "<<fL511<<"  "<<fDL<<"  "<<fTau<<"  "<<fPeriod<<"   blocksPerGrid "<<nBlocksPerGrid<<std::endl;
  sumFreq_kernel<<<nBlocksPerGrid, nThreadsPerBlock>>>
    (fDev_SigArray, d_sigSum, fNFreq, fNBinsInput, 
     fNBinsPerPeriod, DM, fL511, fDL, fTau, fPeriod);
  err = cudaGetLastError();
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to launch sumFreq kernel (error code %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  
  // Copy the device result vector in device memory to the host result vector in host memory.
  //  printf("Copy output data from the CUDA device to the host memory:");
  //  std::cout<<"address: "<<sigSum<<"  "<<sigSum[0]<<"  "<<d_sigSum<<"  "<<size_output<<"  "<<cudaMemcpyDeviceToHost<<std::endl;
  err = cudaMemcpy(sigSum, d_sigSum, size_output, cudaMemcpyDeviceToHost);
  
  if (err != cudaSuccess){
    fprintf(stderr, "Failed to copy compensated signal vector from device to host (error code: %s)!\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  
  for (int i=0; i<fNBins; i++){
    if (sigSum[i]==sigSum[i]&&i<fNBins) fHCompSig[iStep]->SetBinContent(i+1,sigSum[i]);
    else fHCompSig[iStep]->SetBinContent(i+1,fNFreq);
  }
  
  free(sigSum);
  cudaFree(d_sigSum);

  fHCompTiming->Fill(stwch.RealTime(),1);

  return 0;  
}

int PFScan::DoCuFFT(int iStep)
{
  TStopwatch stwch;
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
    //allocate host vector with non-zero elements and image
  //use common number of bins for all runs (fNBins)
  
  // std::cout<<"NBINS: "<<fNBins<<std::endl;

  int signal_memSize=sizeof(float)*fNBins;
  float* h_sigSum_NZ = (float*)malloc(signal_memSize);
  //fill the vector
  for (int i=0; i<fNBins; i++){
    float bico=fHCompSig[iStep]->GetBinContent(i+1);
    if (bico!=0&&bico==bico&&i<fNBins) h_sigSum_NZ[i]=bico;
    else h_sigSum_NZ[i]=fNFreq;
  }
  
  //allocate device vector with non-zero elements:
  int image_memSize = sizeof(float2) * (floor(fNBins/2)+1);
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
  if (cufftPlan1d(&plan, fNBins, CUFFT_R2C, 1)!=CUFFT_SUCCESS){
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
  
  //  int image_memSize = sizeof(float2) * fNBins;
  err= cudaMemcpy(h_fImage, d_fImage, image_memSize, cudaMemcpyDeviceToHost);	
  
  if (err!=cudaSuccess){
    std::cout<<"Failed to copy fImage from device to host"<<std::endl;
    return 1;
  }   
  
  //fill the histogram
  for (int i=0; i<fNBins; i++){
    if (i<floor(fNBins/2)+1)	{  
      float amplitude=sqrt(pow(h_fImage[i].x,2)+pow(h_fImage[i].y,2));
      if (amplitude==amplitude) fHCompSig_FFTImage[iStep]->SetBinContent(i+1,amplitude);
      else std::cout<<"DMScan::sumFrequencies_GPU: nan in FFT bin "
		    <<i<<", iStep: "<<iStep<<std::endl; 
	     //sumFreq_fImage[iStep]->SetBinContent(i+1,amplitude);
    }  
    else{
      int ind=fNBins-i;
      float amplitude=sqrt(pow(h_fImage[ind].x,2)+pow(h_fImage[ind].y,2));
      if (amplitude==amplitude) fHCompSig_FFTImage[iStep]->SetBinContent(i+1,amplitude);
      else std::cout<<"DMScan::sumFrequencies_GPU: nan in FFT bin "
		    <<i<<", iStep: "<<iStep<<std::endl;
    }	
  }

  free(h_sigSum_NZ);
  free(h_fImage);
    
  cudaFree(d_fImage);
  cudaFree(d_sigSum_NZ);
  cufftDestroy(plan);

  fHFFTTiming->Fill(stwch.RealTime(),1);

  return 0;
}


int PFScan::CloseGPU()
{
  // Error code to check return values for CUDA calls
  cudaError_t err = cudaSuccess;
  
  // Free device global memory
  err = cudaFree(fDev_SigArray);
  
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
