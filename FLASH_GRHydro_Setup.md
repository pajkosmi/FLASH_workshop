### Setting up FLASH Starting from Github clone

* Clone github repo from `ECP-Astro/FLASH5`
* Switch to branch `964SparkGRMHD`
* Creating a makefile for your system
  * Enter directory `FLASH5/sites/`
  * Create a new directory with the name of your system
  * Copy `Makefile.h` from `nagini.uchicago.edu` into your new directory
  * Edit `Makefile.h` in your new directory, specifying relevant paths & libraries 
* Return to `FLASH5/` directory
* Run `./setup TOVStar -auto -1d +spherical -debug -nxb=64 +sparkGR +amrex +serialio -unit=IO/IOMain/hdf5/serial/AM -objdir=obj_TOVStar_evolve_AMReX useGRMHD=True`
  * If successful, the term `SUCCESS` should appear on the last line
* Enter the newly built directory:  `obj_TOVStar_evolve_AMReX`
  * This directory contains all the relevant linked files for our GRHydro problem
* Execute `./make` (or `./make -j4` if you have 4 cores)
  * If successful, the term `SUCCESS` should appear on the last line
* To run the TOV Star test problem, exceute the newly built binary: `./flash5`
