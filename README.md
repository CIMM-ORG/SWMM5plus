# SWMMengine
Prototype Fortran 2008 engine for the EPA Storm Water Managemetn Model (SWMM)

## How to compile
To compile the software using the GNU gfortran compiler simply run ./Allmake.sh
at your terminal. Make sure the file is executable. 

Once it is compiled, run ./SWMM to run the default case. Running ./SWMM ./[inputfile].inp 
will enable the program to pull information from a SWMM5 input file.

./Allclean is also provided to clean the
case after simulation. 

- Once Lib_FSWMM is downloaded, it is necessary to comment the necessary lines in interface.f08 to stop debugging
## HAPPY SWMMing!
