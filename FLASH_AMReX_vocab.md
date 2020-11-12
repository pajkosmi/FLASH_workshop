## FLASH & AMReX Data Structures (ExaStar)

This document is meant to guide the reader from the inner workings of FLASH (specificially the MHD solver `Spark`) to the AMReX multifab data structures.  For a more holistic overview of `FLASH`, take a look at the **Features** and **Architecture** sections of the `workshop.md` guide.  For the purposes of ExaStar, the changes that incorporate general relativity reside in the `964SparkGRMHD` branch on the ECP-Astro/FLASH5 github.

**FLASH**:

* In FLASH (an Eulerian code), basic subsets of the domain are called 'blocks' or 'tiles'
  * There is a distiction between blocks and tiles but it is beyond the scope of this writeup.  For now we'll only worry about blocks.
* The hydro solver that evolves the fluids on these blocks is called `Spark` 
  * Once in the `FLASH5` directory, `Spark` can be found in the `source/physics/Hydro/HydroMain/Spark/general_relativity/` directory
* Commonly 
