PulsarFinder: program for PRAO LPI BSA data processing

INSTRUCTIONS:

1). Install fftw 3.3.4 (http://www.fftw.org/download.html).

2). Install root v6.06 (https://root.cern.ch/downloading-root).
Before the compilation root should be configured to use fftw libs:

>./configure --prefix=<...> --with-fftw3-incdir="path to fftw headers" --with-fftw3-libdir="path to fftw libs"

3). Add root environment to your shell startup script (example for ~/.bashrc):

>export ROOTSYS="full path to your root installation directory"

>export PATH=${PATH}:${ROOTSYS}/bin

>export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib/root

4). Install CUDA v7.5 (https://developer.nvidia.com/cuda-downloads).

5). Set the CUDAHOME variable in your shell startup script (example for ~/.bashrc):

>export CUDAHOME="full path to your cuda installation dir"

6). Install the latest version of git (https://git-scm.com/downloads).

7). Download and install the PulsarFinder program:

>git clone https://github.com/gsafronov/pulsarFinder.git

>cd pulsarFinder

>./configure --prefix="path where you want to install PF (default is /usr/local/)"

>make

>make install

8). Add PF installation directory to your PATH (bash example):

export PATH=${PATH}:"path used during configure step"/bin



Examples of program configure files are in the ./examples directory.

At first RAW data should be converted to root histograms, this is performed by "readRAW" program. An example configure file is ./examples/readRAW/readRAW_0301.cfg. To run the readRAW program use this command:

>readRAW -f "config file"


After some number of runs is converted to the root format one can run de-dispersion for different DM values, search for peaks in de-dispersed data and search for peaks in fourier images of de-dispersed data. All this is performed by runPFinder program. An example config file is located in ./examples/PFinder/runPF_0301.cfg. To run the runPFinder program use this command:

>runPFinder -f "config file"