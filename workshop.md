# FLASH Workshop
## O'Shea Group Meeting 9 July 2020
### Mike Pajkos

Athena reference url: https://hackmd.io/qpV_Vo07SSWta5Ha6m9TkA?view

* [History & Documentation](#flash-history-and-documentation)
* [Features](#flash-features)
* [Architecture](#flash-architecture)
* [Get the Code](#downloading-flash)
* [Test Problem](#building-a-test-problem)
* [Run on HPCC](#running-on-hpcc)
* [New Test Problem](#creating-a-new-test-problem)
* [Visualization](#visualization-and-data-analysis)

### FLASH History and Documentation
* FLASH is a legacy code having been developed for over 2 decades
* FLASH Center website: http://flash.uchicago.edu/site/index.shtml
* Code overview: http://flash.uchicago.edu/site/flashcode/
* User guide: http://flash.uchicago.edu/site/flashcode/user_support/flash4_ug_4p62.pdf

### FLASH Features
* physics/grid features from website 
* what it's lacking compared to other codes (see Claire/forrest's athena/enzo comments)

### FLASH Architecture
* FLASH is composable; at runtime it picks and chooses the necessary 'components'
  * Cuts down ~2e6 total lines of code to ~2e4 for a typical simulation 
* FLASH is separated into 6 main units:
  * Physics: contains relevant 'science' for a given simulation
  * Simulation: contains the initial conditions to start a simulation
  * Driver: evolves the simulation
  * IO: manages input & output
  * Infrastructure: handles things like how the grid is setup
  * Monitor: tracks the simulation progress with log files and timing measures
* API and inheritance structure

### Downloading FLASH
* Please clone code from the master branch found at: https://github.com/snaphu-msu/BANG
  * By downloading and working with FLASH, you agree to the following terms: http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node3.html

### Building a Test Problem
* Proper makefile (taken care of already for hpcc)
* the setup process
* the make process
* ready to launch executable
* the flash.par file

### Running on HPCC
* Outline batch submission
* different scaling runs
* highlight where to look in log file
* Scaling plots?

### Creating a New Test Problem
* (Mike this you will walk through)
* New directory in simulation
* Config: why it's there, what it is, what's needed
* Simulation_init(Block)
* Customizing the flash.par file

### Visualization and Data Analysis
* VizIt for quick diagnostic
* yt For pretty pictures and robust analysis


-----------------------------------------------------

#FLASH Workshop OShea Group
##Mike Pajkos

## Before the workshop 
* Getting the code and compiling the code
  * Load the necessary modules
   
   ```
   Module list here
   asd
   ```	
  * Clone the repository from the website into your hpcc home directory (mike ensure makefile = gnu works) using the https:/ link and enter the new `BANG/` directory
  * On the command line run `./setup Sedov -auto -2d -debug -nxb=18 -nyb=18 +spark +pm4dev -gridinterpolation=native -parfile=test_paramesh_2d.par -objdir=obj_Sedov_2D -makefile=gnu`
  * This will create an 'object directory' called `obj_Sedov_2D`; enter the `obj_Sedov_2D` directory
  * Run `make -j` to compile the code
  * Within `obj_Sedov_2D`, there will be an executable and parameter (`.par`) file: `flash4` and `flash.par`, respectively
  * It is standard to run FLASH simulations in a directory outside of the object directory (I'll explain why in the workshop), so return to your HPCC home directory and create a directory for running this simulation: `run_Sedov_FLASH/`
  * Copy your executable (`flash4`) and `.par` file (`flash.par`) from your object directory `~/BANG/obj_Sedov_2D/` into your new run directory `~/run_Sedov_FLASH/`
  * Mike step for downloading slm.sub script
    
  ```
  slm stuff
  asdf
  asdf
  ```
  
  * You're now all set up for the workshop.  If you'd like, feel free to see the output by running `./flash4`

## Features:
### Capabilites
* Strengths
  * Multiphysics: (Relativistic) MHD, Gravity, Radiation Transport, Equations of State, Nuclear reactions, Particles, Cosmology
  * Adaptive mesh refinement
  * OpenMP/MPI compatible
  * 

* Weaknesses
  * Uniform timestepping (no time subcycling)
  * 

### Architecture
* FLASH source tree, API, inheritance structure


### Units 
#### Infrastructure

#### Driver

#### IO

#### physics

#### Simulation

## Setting up a Simulation

### Explaining the setup process in detail
* `./setup` is used to 'pick and choose' the lines of code FLASH needs for a simulation
  * Formally it links the needed source code into the object directory
  * We use a variety of setup *options* (ex. `-2d`, `-nxb=18`) & setup *shortcuts* (ex. `+spark`, `+pm4dev`) to specify everything from grid geometry to what kind of hydro solvers to use (Full list available in chapter 5.1 of user guide)
* `make -j` compiles the selected code (in parallel) within the object directory
  * This is how we get the files we're interested in, namely the executable and parameter file
  * The code is compiled based on the `Makefile` for dev-intel16.  If you look in the `BANG/sites/` directory, there are a variety of directories for different 'sites' or hardware platforms FLASH runs on.  
  * If you were to download FLASH on your laptop, you would add a directory with the same name as the `hostname` of your laptop and add a Makefile specifying the file paths needed for various libraries & compilers
* `flash4` is the binary executable file
  * Remember this is built *for your system* (in our case hpcc).  You could not copy this over to a machine with a different configuration (ex. your laptop).  Instead you would have to setup & make again on the different machine.
  * The binary can be run in serial `./flash4` or parallel `mpirun -np NPROCS_HERE ./flash4`, given the proper parallel libraries were specified in the Makefile
* `flash.par` is the parameter file
  * The parameter file contains runtime parameters that allow the user to partially control the runtime environment (This must be in the same directory as the binary executable)
  * Some examples: the energy of the sedov explosion, physical size of the domain, CFL parameter, or how often to output a file (full list of possible parameters is generated in the `setup_params` file within the object directory

## Output Files
* talk about `.dat` `.log` `plt` `chk` files

## Scaling Study
* directions for different procs.  Maybe sean's plots here?
 
## Visualization
  
* yt & VizIt capabilites.  When to use which.
  
  
## Creating a New Problem
  
* New directory in SimulationMain
* Config file
* flash.par file
* Simulation Init
* Simulation Init block


## New physics?
* New directory in physics unit.  Account for REQUIRES in Config? 
