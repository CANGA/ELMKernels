CLM/ELM indexes column variables in a manner that cannot be reproduced using C++. In this project, 
we use basic C-style arrays with an index starting at 0 for storage of column data. Described below 
are the conventions we use to map column data indices between CLM/ELM and this project.  

SUBSURFACE VARIABLES (CLM example: (1:nlevgrnd) or (1:nlevsoi))
SHIFT CLM INDEX BY -1
CLM/ELM variables that exist only in the subsurface are indexed top to bottom from 1:nlevgrnd (typically 15), 
where var(1) is the topmost layer, and var(nlevgrnd) is the bottom layer. This project's equivalent 
data structure will be shifted by -1, so the top layer is var[0] and the bottom layer is var[nlevgrnd-1].

ABOVE SURFACE VARIABLES (CLM example: (-nlevsno+1:0) or (-nlevsno+1:1))
SHIFT CLM INDEX BY  nlevsno - 1
CLM/ELM variables that exist on and above the surface (snow, canopy) are indexed from top to bottom,
like (-nlevcan+1:0) (only above ground component) or (-nlevsno+1:1) (ground surface and above ground component).
This project's equivalent data structure will be shifted by nlevsno - 1, (or nlevcan - 1) so the top layer is 
var[0] and the bottom layer is var[nlevsno] (if ground surface included) or var[nlevsno] (if ground surface excluded). 

SURFACE SUBSURFACE VARIABLES (CLM example: (-nlevsno+1:nlevgrnd))
SHIFT CLM INDEX BY  nlevsno - 1
CLM/ELM variables that exist both in the snowpack/canopy and through the subsurface are indexed from top to bottom,
following logically from the convention discussed above. This project's equivalent data structure will be shifted by 
nlevsno - 1, similar to the surface data. The surface/subsurface data


  IMPORTANT INDEX CONVERSIONS:
  SUBSURFACE VARIABLES
  top layer - CLM: var(1), THIS PROJECT: var[0]
  bottom layer - CLM: var(nlevgrnd), THIS PROJECT: var[nlevgrnd-1]

  ABOVE SURFACE DATA AND ABOVE AND BELOW SURFACE DATA:
  snow or canopy bottom layer (adjacent to ground surface) - CLM: var(0), THIS PROJECT: var[nlevsno-1]
  surface layer in above surface data (or top subsurface) - CLM: var(1), THIS PROJECT: var[nlevsno]
  top active snow layer -  CLM: var(snl+1), THIS PROJECT: var[nlevsno-snl] -- snl is neg in CLM, pos here
  bottom subsurface layer - CLM: var(nlevgrnd), THIS PROJECT: var[nlevgrnd+nlevsno-1]

  CONCURRENTLY WORKING WITH SUBSURFACE ONLY AND ABOVE SURFACE DATA:
  subsurface data's index is offset from above surface data by nlevsno  
    AboveAndBelow[i + nlevsno] is same layer as BelowOnly[i]  


THOUGHTS:
CLM's data structures conveniently result in the top subsurface layer's index = 1. This library will not be 
able to reproduce this feature. This shouldn't cause problems, but will require index conversion to occur 
throughout the code, and special consideration when concurrently operating on subsurface only and above ground data.  

ABOVE GROUND DATA
CLM variable indices      indices
|-4|   top snow layer     |0|
 --                       ---
|-3|   snow layer         |1|
 --                       ---
|-2|   snow layer         |2|
 --                       ---
|-1|   snow layer         |3|
 --                       ---
| 0|   bottom snow layer  |4|  [nlevsno - 1]
-----  ground interface  ------
| 1|   top soil layer     |5|  [nlevsno]


ABOVE/BELOW GROUND DATA
CLM variable indices      indices
|-4|   top snow layer     |0|
 --                       ---
|-3|   snow layer         |1|
 --                       ---
|-2|   snow layer         |2|
 --                       ---
|-1|   snow layer         |3|
 --                       ---
| 0|   bottom snow layer  |4|  [nlevsno - 1]
-----  ground interface  ------
|1|  top subsurface layer |5|  [nlevsno]
 --                       ---
|2|    soil layer         |6|
 --                       ---
|3|    soil layer         |7|
 --                       ---
|4|    soil layer         |8|
.
.
--                        ---
|nlevgrnd|   botom layer  |nlevgrnd + nlevsno - 1|


BELOW GROUND DATA
CLM variable indices      indices
-----  ground interface  ------
| 1| top subsurface layer |0| 
 --                       ---
|2|    soil layer         |1|
 --                       ---
|3|    soil layer         |2|
 --                       ---
|4|    soil layer         |3|
.
.
--                        ---
|nlevgrnd|   botom layer  |nlevgrnd - 1|