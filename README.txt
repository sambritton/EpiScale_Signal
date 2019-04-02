CUDA code for developmental biology using Subcellular Element Method

Hardware requirement: 
Nvidia video card that supports SM 2.0+ and CUDA 4.0 

Software environment requirement: 
CMAKE ----- Build system.
CUDA  ----- Provide runtime support for parallel GPU computation.
Thrust ---- Build-in library of cuda, similar to STL of C++
Paraview -- (Optional) Visualization software for animation purpose. 


Location of configuration files:
 ./resources
*******************************************
To run simulation on CRC clusters "acms.crc.nd.edu" which are based on SGE (Sun Grid Engine) cluster software :
   (1) git clone https://github.com/AliNemat/EpiScale_Signal.git
   (2) module load cmake gcc/7.1.0 cuda/10.0 bertini matlab  
   (3) In the directory  ~/SceCells write the command "ccmake . 
   (4) In the directory ~/SceCells write the command "make"
   (5)Submit your simulation with the command "qsub EpiScale_run2.sh"  # Note: Other .sh files in ~/script are not active anymore#


