% Atmospheric Turbulence Simulation
% Olivier Guyon
% Oct 2016

# Overview

## Scope

Atmospheric Turbulence Simulation

- Multi-layer simulation
- Uses Fresnel propagation engine between optical elements
- Computes phase, amplitude
- Accurate treatment of chromatic effects 
- FITS files output


## Usage

Scripts to run the software are located within the source code directory:

	./AtmosphericTurbulenceSimul/src/AtmosphericTurbulence/scripts/

The scripts can be linked to your working directory by executing the following command:

	ln -s $PWD/syncscripts /myworkdirectory/syncscripts

Then, execute in your work directory:

	./syncscripts

This will install all required scripts in workdirectory.

Code is composed of a several layers (from high to low) :

-------------------- -----------------------------------------------------------
Script                   Description
-------------------- -----------------------------------------------------------
**runturb**          Top level script

**mkHVturbprof**     Create Hufnager-Valley turbulence profile
-------------------- -----------------------------------------------------------


# Quick Start

Execute main script with -h option for help:

~~~
./runturb -h
~~~

The script takes a single argument, the wavelength [micron] :

~~~
./runturb 1.65
~~~


# Refractive index

By default, Sellmeier equations are used to compute atmosphere refractive index. You can also edit the runturb script to point to the location of "RIA file" (Refractive Index Absorption). See file for details.

