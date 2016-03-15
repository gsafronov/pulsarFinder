#include "DMScan.h"

//interface to multithreading

struct sumFreqID{
  DMScan* dmscan;
  int iThread;
  int iStep;
};


void* sumFreq_CPU_interface(void *threadarg)
{
  struct sumFreqID *sfID;
  
  sfID = (struct sumFreqID *) threadarg;
  
  int rv = (sfID->dmscan)->sumFrequencies_CPU(sfID->iThread, sfID->iStep);

  pthread_exit(NULL);
}


void* sumFreq_GPU_interface(void *threadarg)
{
  struct sumFreqID *sfID;
  
  sfID = (struct sumFreqID *) threadarg;
  
  int rv = (sfID->dmscan)->sumFrequencies_GPU(sfID->iStep);

  pthread_exit(NULL);
}

