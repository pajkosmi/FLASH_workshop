# FLASH Workshop
## O'Shea Group Meeting 9 July 2020
### Mike Pajkos

* [History & Documentation](#flash-history-and-documentation)
* [Before the Workshop](#before-the-workshop)
* [Features](#flash-features-at-a-glance)
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
 * Log onto HPCC, dev-intel18
 * Load the necessary modules

  ```
module purge
module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load HDF5/1.8.20
  ```	
 * Please clone code using the https protocol found at: https://github.com/snaphu-msu/BANG and enter the new `BANG/` directory 
   * By downloading and working with FLASH, you agree to the following terms: http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node3.html
 * Switch branches using `git checkout Workshop_OShea`
 * On the command line run `./setup Sedov -auto -2d -debug -nxb=18 -nyb=18 +spark +pm4dev -gridinterpolation=native -parfile=workshop_flash.par -objdir=obj_Sedov_2D -makefile=gnu`
 * This will create an 'object directory' called `obj_Sedov_2D`; enter the `obj_Sedov_2D` directory
 * Run `make -j` to compile the code
 * Within `obj_Sedov_2D`, there will be an executable and parameter (`.par`) file: `flash4` and `flash.par`, respectively
 * It is standard to run FLASH simulations in a directory outside of the object directory (I'll explain why in the workshop), so return to your HPCC home directory and create a directory for running this simulation: `run_Sedov_FLASH/`
 * Copy your executable (`flash4`) and `.par` file (`flash.par`) from your object directory `~/BANG/obj_Sedov_2D/` into your new run directory `~/run_Sedov_FLASH/`


 * You're now all set up for the workshop.  If you'd like, feel free to see the output by running `./flash4`


### FLASH Features at a Glance
* Strengths
  * Flexibility to construct new test problems
  * Multiphysics: (Relativistic) MHD, Gravity, Radiation Transport, Equations of State, Nuclear reactions, Particles, Cosmology
  * Adaptive mesh refinement (AMR) & multiple geometries
  * OpenMP/MPI compatible
  * Parallel IO compatibility with HDF5

* Weaknesses
  * Only uniform timestepping (no time subcycling)
  * Limited compatibility with GPUs (Jared is currently working on this)
  * AMR is restricted to octree
    * No 3D spherical geometry
    * No logarithmically sized cells in spherical geometry
  

### FLASH Architecture
* FLASH is written mostly in FORTRAN90
* FLASH is composable; at runtime it picks and chooses the necessary 'components' it needs to model a physical situation
  * Cuts down ~2e6 total lines of code to ~2e4 for a typical simulation 
* It has a 'tree like' structure beginning with broad, parent 'Unit directories' (always beginning with a capital letter) and branch to more specific functions and routines
* In the high level directories there are often 'stub implementations'
  * These are essentially 'no ops' that allow FLASH to compile easily even if a certain features of the code are left out (MIKE EXAMPLE grid refine)
* Kernels are typically found at the finest (leaf) level of the FLASH architecture 
  * This design helps modularity, by allowing new functionalities/solvers to be integrated into FLASH without rewriting existing code 

* Below outlines 5 of the main units in FLASH:
#### Driver
* Controls initialization and evolution of simulations
* Organizes interaction between units
* Responsible for launching the parallel environment

#### Infrastructure
* Grid
  * Stores the data the physics unit will access
  * Manages Eulerian grid (spatial domain) and movement of LaGrangian tracer particles
    * Cartesian, spherical, cylindrical, & polar geometries
  * Handles Boundary conditions and solving partial differential equations on the grid (ex Multipole solver)
  * Can use uniform or adaptive grid (Octree based)
* Input/Output (IO)
  * Serial and parallel HDF5 output
  * Creates files with various simulation (meta)data (outlined in **Output Files** section below)
* Runtime parameters
  * Stores and maintains runtime parameters used in the simulation
  * Handles default parameters (found in `Config` files) and those user defined (found in `flash.par` file)
* Multispecies
  * Tracks multiple kinds of fluids (ex. air vs water or Ni56 vs C12)
  * Typically treated as mass scalars advected along with the fluid
* Physical Constants
  * Stores commonly used physical constants in various units
  * New constants can be added by the user in `PhysicalConstants init`
  * ex: pi, Newton's gravitational constant, mass of an electron, Avogadro's number,...

#### Physics
* (Relativistic) MHD
  * Manages fluid flow
  * Directionally split (PPM) & unsplit solvers
  * `Spark` is the current GRMHD solver currently in development (today we will use its Newtonian version)
    * Reconstruction techniques: TVD, FOG, WENO5, PPM, MP5
    * Directionally unsplit, finite-volume, compressible
    * Riemann solvers: HLLC, HLLD (MHD), HLLE (GRMHD)
    * Divergence cleaning (glm) method for maintaining divergence(B) = 0
    
* Equation of State (EOS)
  * Ensures thermodynamic consistency between variables for hydro and nuclear burning
  * Different EOSs available: perfect gas (with multiple adiabatic indices), Helmholtz (accounts for relativistic/degenerate matter), & nuclear (tabulated)
  
* Local Source Terms
  * Contain terms responsible for emitting/absorbing energy 
    * Nuclear burning (7, 13, & 19 isotope reaction networks)
    * Ionization in plasmas: He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni, H, & electrons
    * Stirring or driving terms in hydro simulations
    * Energy deposited by lasers incident on the computational domain
    * Heat exchange between ion, electron, & radiation components
    * Advection-diffusion-reactions (commonly called flames or deflagration fronts)

* Diffusive Terms
  * Implements diffusive effects: head conduction, viscosity, & mass diffusivity

* Gravity 
  * Computes gravitational source terms
    * Constant, plane-parallel, point mass, user defined, & self gravity
      * Multipole and multigrid methods to solve Poisson's equation

* Particles
  * Follows active (massive, charged, sink) or passive (tracer) particles
  
* Cosmology
  * Solves Friedmann equation for scale factor in expanding universe
  * Applies cosmological redshift to hydro quantities
    * Calculations assumed to take place in comoving coordinates
  * See Chapter 21 for more details
  
* Material Properties
  * Computes thermal conductivity, magnetitic resistivity, viscocity, & opacity values


* Radiative Transfer
  * Solves the radiative transfer equation
    * Multigroup diffusion
    * Neutrino leakage (exclusive to this development branch)
    * M1 neutrino treatment (exclusive to this development branch)
 
* High Energy Density Physics Experiments
  * Multi-temperature treatment of plasmas: ion, electron, & radiation (see Chapter 13 for relevant equations)
  

#### Monitor
* Logfile
  * Manages output log (more specifics on log file in **Output Files** section)
* Timer
  * Manages timing routines to monitor performance
  * Users can manually place labelled timers to measure performance for specific parts of the code or find load imbalances
* Profiler
  * Interface provided for third-party profilers or tracing libraries

#### Simulation
* Contains (M)HD, gravity, particle, burn, and radiative transfer test problems
* This is where we would create a new simulation

#### Other Units
* Note there are other units avaialable beyond the scope of this workshop.  Below they are mentioned
  * Proton Imaging & Proton Emission
  * Thomson Scattering
  * Xray Imaging
  * Piecewise cubic interpolation
  * Quadratic, cubic, quartic root finders
  * RungeKutta integration
  

### Building a Test Problem
* `./setup` is used to 'pick and choose' the lines of code FLASH needs for a simulation
  * Formally it links the needed source code into the object directory
  * We use a variety of setup *options* (ex. `-2d`, `-nxb=18`) & setup *shortcuts* (ex. `+spark`, `+pm4dev`) to specify everything from grid geometry to what kind of hydro solvers to use (Full list available in chapter 5.1 of user guide or in `BANG/bin/Readme.SetupVars`)
* `make -j` compiles the selected code (in parallel) within the object directory
  * This is how we get the files we're interested in, namely the executable and parameter file
  * The code is compiled based on the `Makefile.h` for dev-intel18.  If you look in the `BANG/sites/` directory, there are a variety of directories for different 'sites' or hardware platforms FLASH runs on.  
  * If you were to download FLASH on your laptop, you would add a directory with the same name as the `hostname` of your laptop and add a `Makefile.h` specifying the file paths needed for various libraries & compilers
* `flash4` is the binary executable file
  * Remember this is built *for your system* (in our case hpcc).  You could not copy this over to a machine with a different configuration (ex. your laptop).  Instead you would have to setup & make again on the different machine.
  * The binary can be run in serial `./flash4` or parallel `mpirun -np NPROCS_HERE ./flash4`, given the proper parallel libraries were specified in the `Makefile.h`
* `flash.par` is the parameter file
  * The parameter file contains runtime parameters that allow the user to partially control the runtime environment (This must be in the same directory as the binary executable)
  * Some examples: the energy of the sedov explosion, physical size of the domain, CFL parameter, or how often to output a file (full list of possible parameters is generated in the `setup_params` file within the object directory


### Running on HPCC
* Today we will perform scaling tests for a sedov explosion
* Enter the node reservation following Brian's instructions
* Submit an interactive job with 16 cores using `salloc -N 1 -c 16 --time=0:20:00`
* Strong Scaling (Fixed problem size, different number of processors)
  * Enter the `~/run_Sedov_FLASH/` directory 
  * Edit the `flash.par` parameter `log_file = "sedov.log"` to `log_file = "sedov_strong_1P.log"`
  * Run `mpirun -np 1 ./flash4`
  * Edit the `flash.par` parameter `log_file = "sedov.log"` to `log_file = "sedov_strong_2P.log"`
  * Run `mpirun -np 2 ./flash4`
  * Repeat the above steps for 4, 8, & 16 processors
  * Type `less sedov_strong_1P.log` & go to the end of the file by typing `G`
  * In the row labeled `accounting unit` is a column labeled `avg/proc (s)` (means average time per processor)
    * Follow the `avg/proc (s)` column down to the `evolution` row.  Record this number for all of your `.log` files
  * Plot these times vs number of processors (1, 2, 4, 8, & 16) to get a strong scaling plot
* Weak scaling (Increasing problem size with number of processors)
  * Enter the `~/run_Sedov_FLASH/` directory 
  * Edit the `flash.par` parameter to set `log_file = "sedov_weak_1P.log"`
  * In `flash.par` ensure `nblockx = 1` & `nblocky = 1`
  * Run `mpirun -np 1 ./flash4`
  * Edit the `flash.par` parameter to set `log_file = "sedov_weak_4P.log"`
  * In `flash.par` ensure `nblockx = 2` & `nblocky = 2`
  * Run `mpirun -np 4 ./flash4`
  * Edit the `flash.par` parameter to set `log_file = "sedov_weak_16P.log"`
  * In `flash.par` ensure `nblockx = 4` & `nblocky = 4`
  * Run `mpirun -np 16 ./flash4`
  * Once again at the bottom of each `.log` file, record the entry in the `avg/proc (s)` column and `evolution` row
    * Plot these times vs number of processors (1, 4, & 16) to get a weak scaling plot
* Maybe Sean's production plots here?

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
#### VizIt
* Useful for quick diagnostics or 'order of magnitude estimates'
* Select `FLASH` filetype when loading in `chk` files
* Found at: https://wci.llnl.gov/simulation/computer-codes/visit/

### yt 
* Useful for pretty pictures and robust analysis
* Load in HDF5 `chk` files with `ds = yt.load('FLASH_chk_file_here')`
* Found at: https://yt-project.org/

### New physics
* Incorporating new physics into FLASH comes in two flavors
  * Updating existing physics (ex: Hydro --> MHD)
    * Add new features to existing routines, like in `physics/Hydro/HydroMain/Spark/` 
  * Creating new physics modules alltogether (ex: including lattice QCD)
    * Write new routines from scratch and give them their own directory (ex: `physics/LatticeQCD/`)
* Make sure to include the new physics in the `Config` file in the Simulation directory 
