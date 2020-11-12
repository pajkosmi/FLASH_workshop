## FLASH & AMReX Data Structures (ExaStar)

This document is meant to guide the reader from the inner workings of FLASH to the AMReX multifab data structures.  For a more holistic overview of `FLASH`, take a look at the XXX section of the `workshop.md` guide.

**FLASH**:

* In FLASH (an Eulerian code), basic subsets of the domain are called 'blocks'  
* The hydro solver that evolves the fluids on these blocks is called `Spark` 
