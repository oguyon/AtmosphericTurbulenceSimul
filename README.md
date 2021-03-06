# Turbulence Simulator  {#mainpage} 


Simulates effect of atmospheric turbulence on wavefronts at visible and near-IR wavelength.
Includes chromatic diffractive propagation between layers, inner, outer scales, wind speed and direction for each layer, CN2 profile. 

OPD and scintillation can be computed for multiple wavelengths and multiple directions.

Output can be FITS files (cubes), or shared memory image for real-time use by other programs. In shared memory mode, physical time and computer clock time can be synced, with a one-to-one match or a slowing down / speeding up factor.


## Downloading source code
Latest distribution is on [github](
https://github.com/oguyon/).
You can clone the repository, or download the latest .tar.gz distribution.


## Compilation
The source code follows the standard GNU build process:

./configure

make

make install


## Documentation 
Please consult the [online documentation]{http://oguyon.github.io/AtmosphericTurbulenceSimul/index.html} (generated by doxygen).


## Libraries
The following libraries are used:
- readline, for reading the command line input
- flex, for parsing the command line input
- bison, to interpret the command line input
- fftw, for performing Fourier Transforms
- gsl, for math functions and tools
- fitsio, for reading and writing FITS image files

## Source Code Architecture 
Written in C.
The main is a command line interface (CLI). Source code is in CLIcore.c and CLIcore.h.
Key data structures (such as the image data structure) are declared in CLIcore.h.

## How to run the turbulence simulator
Copy scripts from ./src/AtmosphericTurbulence/scripts directory to working directory.
Edit and execute the main script "runturb"


## Credits
This software was developed with support from the National Science Foundation (award #1006063) 
