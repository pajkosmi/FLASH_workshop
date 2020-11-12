## FLASH & AMReX Data Structures (ExaStar)

This document is meant to guide the reader from the inner workings of FLASH (specificially the MHD solver `Spark`) to the AMReX multifab data structures.  For a more holistic overview of `FLASH`, take a look at the **Features** and **Architecture** sections of the `workshop.md` guide.  For the purposes of ExaStar, the changes that incorporate general relativity reside in the `964SparkGRMHD` branch on the `ECP-Astro/FLASH5` github.

**From FLASH...**

* In FLASH (an Eulerian code), basic subsets of the domain are called 'blocks' or 'tiles'
  * There is a distiction between blocks and tiles but it is beyond the scope of this writeup.  For now we'll only worry about blocks
* The hydro solver that evolves the fluids on these blocks is called `Spark` 
  * Once in the `FLASH5` directory, `Spark` can be found in the `source/physics/Hydro/HydroMain/Spark/general_relativity/` directory
* Accessing/modifying the data in a block is done through `call blockDesc%getDataPtr(solnData,CENTER)` (ex. `hy_rk_updateSoln.F90`)
  * Here `blockDesc` is the 'block descriptor' that specifies which block is being operated on
  * `%getDataPtr` gets a pointer whose target is the block of data 
  * `(solnData,CENTER)` specifies to point to the cell centered data.  The pointer targetting these data is called `solnData`
  * `solnData(:,i,j,k)` is a 4D array indexed by physical quantity (ex density, pressure, ...), x-coordinate, y-coordinate, & z-coordinate.  It is what will be used to construct the stress energy tensor
* In general, the cell centered data of `blockDesc` is a 'physical multifab'.  When using AMReX, these physical multifabs are of the same type as `amrex_multifab`
  * Formally defined in `source/Grid/GridMain/AMR/Amrex/gr_physicalMultifabs.F90` where `solnData()` points to `unk(:)`
  * The `amrex_multifab` data type is defined in `amrex_amr_module`, which resides in the `AMReX` library

**...to AMReX**
