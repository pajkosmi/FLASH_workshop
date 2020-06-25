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


