# FLASH Workshop
## O'Shea Group Meeting 9 July 2020
### Mike Pajkos

Athena reference url: https://hackmd.io/qpV_Vo07SSWta5Ha6m9TkA?view

* [History & Documentation](#flash-history-and-documentation)
* [Before the Workshop](#before-the-workshop)
* [Features](#flash-features)
* [Architecture](#flash-architecture)
* [Test Problem](#building-a-test-problem)
* [Run on HPCC](#running-on-hpcc)
* [Output Files](#output-files)
* [New Test Problem](#creating-a-new-test-problem)
* [Visualization](#visualization-and-data-analysis)
* [New Physics](#new-physics)

### FLASH History and Documentation
* FLASH is a legacy code having been developed for over 2 decades
* FLASH Center website: http://flash.uchicago.edu/site/index.shtml
* Code overview: http://flash.uchicago.edu/site/flashcode/
* User guide: http://flash.uchicago.edu/site/flashcode/user_support/flash4_ug_4p62.pdf

### Before the Workshop
 * Log onto HPCC, dev-intel16
 * Load the necessary modules

  ```
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load HDF5/1.8.20
  ```	
 * Please clone code from the master branch found at: https://github.com/snaphu-msu/BANG and enter the new `BANG/` directory (mike ensure makefile = gnu works)
   * By downloading and working with FLASH, you agree to the following terms: http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node3.html
 * On the command line run `./setup Sedov -auto -2d -debug -nxb=18 -nyb=18 +spark +pm4dev -gridinterpolation=native -parfile=test_paramesh_2d.par -objdir=obj_Sedov_2D -makefile=gnu`
 * This will create an 'object directory' called `obj_Sedov_2D`; enter the `obj_Sedov_2D` directory
 * Run `make -j` to compile the code
 * Within `obj_Sedov_2D`, there will be an executable and parameter (`.par`) file: `flash4` and `flash.par`, respectively
 * It is standard to run FLASH simulations in a directory outside of the object directory (I'll explain why in the workshop), so return to your HPCC home directory and create a directory for running this simulation: `run_Sedov_FLASH/`
 * Copy your executable (`flash4`) and `.par` file (`flash.par`) from your object directory `~/BANG/obj_Sedov_2D/` into your new run directory `~/run_Sedov_FLASH/`
 * Mike step for downloading slm.sub script

 ```
#!/bin/bash -login
#SBATCH --time=00:20:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -J sedov_Xcores
#SBATCH --mem-per-cpu=2G
#SBATCH -e %J.e
#SBATCH -o %J.o
#SBATCH -C NOAUTO:intel16

### load necessary modules, e.g.
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load HDF5/1.8.20


### change to the working directory where your code is located
cd /mnt/home/your_username/path_to_run_directory/run_Sedov_FLASH/

## call your executable
srun -n $SLURM_NTASKS ./flash4 -par_file flash.par
 
### output the resource usage of the job
qstat -f $PBS_JOBID
 ```

 * You're now all set up for the workshop.  If you'd like, feel free to see the output by running `./flash4`


### FLASH Features
* Strengths
  * Multiphysics: (Relativistic) MHD, Gravity, Radiation Transport, Equations of State, Nuclear reactions, Particles, Cosmology
  * Adaptive mesh refinement (AMR) with flux correction on fine/coarse boundaries
  * OpenMP/MPI compatible
  * IO compatibility with HDF5
  * 

* Weaknesses
  * Uniform timestepping (no time subcycling)
  * what it's lacking compared to other codes (see Claire/forrest's athena/enzo comments)
  *


### FLASH Architecture
* FLASH is written mostly in FORTRAN90
* FLASH is composable; at runtime it picks and chooses the necessary 'components' it needs to model a physical situation
  * Cuts down ~2e6 total lines of code to ~2e4 for a typical simulation 
* It has a 'tree like' structure beginning with broad, parent 'Unit directories' (always beginning with a capital letter) and branch to more specific functions and routines
* In the high level directories there are often 'stub implementations'
  * These are essentially 'no ops' that allow FLASH to compile easily even if a certain features of the code are left out (MIKE EXAMPLE grid refine)
* Kernels are typically found at the finest (leaf) level of the FLASH architecture 
  * This design helps modularity, by allowing new functionalities/solvers to be integrated into FLASH without rewriting existing code 

* FLASH is separated into 6 main units:
#### Physics
* contains relevant 'science' for a given simulation
#### Simulation
* contains the initial conditions to start a simulation

#### Driver
* evolves the simulation

#### IO
* manages input & output

#### Infrastructure
* handles things like how the grid is setup

#### Monitor
* tracks the simulation progress with log files and timing measures


* API and inheritance structure


### Building a Test Problem
* `./setup` is used to 'pick and choose' the lines of code FLASH needs for a simulation
  * Formally it links the needed source code into the object directory
  * We use a variety of setup *options* (ex. `-2d`, `-nxb=18`) & setup *shortcuts* (ex. `+spark`, `+pm4dev`) to specify everything from grid geometry to what kind of hydro solvers to use (Full list available in chapter 5.1 of user guide or in `BANG/bin/Readme.SetupVars`)
* `make -j` compiles the selected code (in parallel) within the object directory
  * This is how we get the files we're interested in, namely the executable and parameter file
  * The code is compiled based on the `Makefile.h` for dev-intel16.  If you look in the `BANG/sites/` directory, there are a variety of directories for different 'sites' or hardware platforms FLASH runs on.  
  * If you were to download FLASH on your laptop, you would add a directory with the same name as the `hostname` of your laptop and add a Makefile specifying the file paths needed for various libraries & compilers
* `flash4` is the binary executable file
  * Remember this is built *for your system* (in our case hpcc).  You could not copy this over to a machine with a different configuration (ex. your laptop).  Instead you would have to setup & make again on the different machine.
  * The binary can be run in serial `./flash4` or parallel `mpirun -np NPROCS_HERE ./flash4`, given the proper parallel libraries were specified in the Makefile
* `flash.par` is the parameter file
  * The parameter file contains runtime parameters that allow the user to partially control the runtime environment (This must be in the same directory as the binary executable)
  * Some examples: the energy of the sedov explosion, physical size of the domain, CFL parameter, or how often to output a file (full list of possible parameters is generated in the `setup_params` file within the object directory


### Running on HPCC
* Outline batch submission
* different scaling runs
* highlight where to look in log file
* directions for different procs.  Maybe sean's plots here?

### Output Files
* The data file (`.dat`) file
  * Contains globally summed qunatities: ex total mass, y momentum, magnetic energy
  * User defined quantities can be defined in `IO_writeIntegralQuantities.F90`
* Checkpoint (`_chk_`) files are used to restart a simulation
  * Output at user defined frequency
  * Used as preventative measure in case a system or code crashes
  * Contains variables, species, grid reconstruction data, scalar values, and meta-data
  * Note these files can fill directory space fast if output too often
* Plot (`plt`) files contain info to interpret grid info
  * Stores user defined variables (ex density, pressure) for analyzing/post-processing
  * Not used for restarting simulations so these typically take up less data than `chk` files
* Particle files contain info needed for analyzing/processing particle data
* Log (`.log`) file contains various useful (meta)data , warnings, and errors--if present
  * Can contain relevant git version, `setup` call, and parallel info
  * Holds runtime parameters (`.par` info) and physical constants
  * Also tracks when grid refinement occurs and number of blocks in domain
  * Useful for performance measures by containing timing information in each unit (ex Hydro, Grid, Gravity, ...)
  

### Creating a New Test Problem
* (Mike this you will walk through)
* For a new problem, you will need certain 'ingredients'.  Here we will setup Rayleigh-Taylor instability.
  * First a we setup a directory in `/BANG/source/Simulation/SimulationMain/RayleighTaylor`
  * `Config` files outline the necessary Units, physics, and variables (as well as their default values) needed to model the simulation
    * These files are what is parsed by `setup` and specify exactly which lines of code are to be linked to in the object directory
  * `Makefile` lists which source `.F90` files will be compiled
  * `flash.par` specifies certain runtime parameters (ex. density of top fluid or number of refinement level)
    * Note this is what is copied into the object directory when we run setup 
  * `Simulation_data.F90` stores the variables specified in `flash.par`
  * `Simulation_init.F90` initializes the values in `Simulation_data.F90` based on the values specified in `flash.par`
  * `Simulation_initBlock.F90` assigns these newly initialized values onto the computational domain
    * For this Rayliegh-Taylor example, it initially defines the denser fluid on top and less dense fluid below
  * Run `./setup RayleighTaylor -auto -2d -debug -nxb=18 -nyb=18 +spark +pm4dev -gridinterpolation=native -objdir=obj_RTinstab_2D -makefile=gnu`
  

### Visualization and Data Analysis
* VizIt for quick diagnostic
* yt For pretty pictures and robust analysis


### New physics?
* New directory in physics unit.  Account for REQUIRES in Config? 
