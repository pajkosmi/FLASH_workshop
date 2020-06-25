# FLASH Workshop
## O'Shea Group Meeting 9 July 2020
### Mike Pajkos

Athena reference url: https://hackmd.io/qpV_Vo07SSWta5Ha6m9TkA?view

* [History & Documentation](#flash-history-&-documentation)

### FLASH History & Documentation
* FLASH Center website: http://flash.uchicago.edu/site/index.shtml
* Code overview: http://flash.uchicago.edu/site/flashcode/
* User guide: http://flash.uchicago.edu/site/flashcode/user_support/flash4_ug_4p62.pdf

### FLASH Features
* physics/grid features from website 
* what it's lacking compared to other codes (see Claire/forrest's athena/enzo comments)

### FLASH Architecture
* Outline here or in slides?
* Monitor (Logfile/timer), Driver, IO, physics, Infrastructure (Grid), Simulation
* API and inheritance structure

### Downloading FLASH
* cloning instructions 
* 

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

### Building a New Test Problem
* (Mike this you will walk through)
* New directory in simulation
* Config: why it's there, what it is, what's needed
* Simulation_init(Block)
* Customizing the flash.par file

### Visualization and Data Analysis
* VizIt for quick diagnostic
* yt For pretty pictures and robust analysis


