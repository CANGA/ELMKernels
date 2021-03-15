ELM column variables are indexed in a manner that result in the top subsurface layer's 
index = 1 in every type of column-based data structure. This indexing cannot be reproduced 
using the basic subset of C++ we have constrained ourselves to without introducing unacceptable 
overhead. 

In this project,we use basic C-style arrays with an index starting at 0 for storage of all column 
data. This shouldn't cause problems, but will require index conversion to occur throughout the 
code, and special consideration when concurrently operating on subsurface only and above ground 
data.  Described below are the conventions we use to map column data indices between ELM and 
this project.  

SUBSURFACE ONLY VARIABLES (ELM example: (1:nlevgrnd) or (1:nlevsoi))
SHIFT ELM INDEX BY -1
ELM variables that exist only in the subsurface are indexed top to bottom from 1:nlevgrnd (typically 15), 
where var(1) is the topmost layer, and var(nlevgrnd) is the bottom layer. This project's equivalent 
data structure will be shifted by -1, so the top layer is var[0] and the bottom layer is var[nlevgrnd-1].

ABOVE SURFACE VARIABLES (ELM example: (-nlevsno+1:0) or (-nlevsno+1:1))
SHIFT ELM INDEX BY  nlevsno - 1
ELM variables that exist on and above the surface (snow, canopy) are indexed from top to bottom,
like (-nlevcan+1:0) (only above ground component) or (-nlevsno+1:1) (ground surface and above ground component).
This project's equivalent data structure will be shifted by (nlevsno - 1), (or nlevcan - 1) so the top 
layer is var[0] and the bottom layer is var[nlevsno] (if ground surface included) or var[nlevsno - 1] 
(if ground surface excluded). 

SURFACE SUBSURFACE VARIABLES (ELM example: (-nlevsno+1:nlevgrnd))
SHIFT ELM INDEX BY  nlevsno - 1
ELM variables that exist both in the snowpack/canopy and through the subsurface are indexed from top 
to bottom, following logically from the convention discussed above. This project's equivalent data structure 
will be shifted by (nlevsno - 1), similar to the surface data.


  IMPORTANT INDEX CONVERSIONS:
  SUBSURFACE VARIABLES
  top layer - ELM: var(1), THIS PROJECT: var[0]
  bottom layer - ELM: var(nlevgrnd), THIS PROJECT: var[nlevgrnd-1]

  ABOVE SURFACE DATA AND ABOVE AND BELOW SURFACE DATA:
  snow or canopy bottom layer (adjacent to ground surface) - ELM: var(0), THIS PROJECT: var[nlevsno-1]
  surface/top subsurface layer - ELM: var(1), THIS PROJECT: var[nlevsno]
  bottom subsurface layer - ELM: var(nlevgrnd), THIS PROJECT: var[nlevgrnd+nlevsno-1]
  top active snow layer -  ELM: var(snl+1), THIS PROJECT: var[nlevsno-snl] -- snl is negative in ELM, positive here
  second from top active snow layer -  ELM: var(snl+2), THIS PROJECT: var[nlevsno-snl+1]

  CONCURRENTLY WORKING WITH SUBSURFACE ONLY AND ABOVE SURFACE DATA:
  subsurface data's index is offset from above surface data by nlevsno  
    AboveAndBelow[i + nlevsno] is same layer as BelowOnly[i],
  top subsurface layer - AboveAndBelow[nlevsno] & BelowOnly[0]
  bottom layer - AboveAndBelow[nlevsno + nlevgrnd - 1] & BelowOnly[nlevgrnd - 1]


VISUAL INDEX MAPPING:

ABOVE GROUND DATA
ELM variable indices    new indices
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
ELM variable indices    new indices
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
|nlevgrnd|   bottom layer  |nlevgrnd + nlevsno - 1|


BELOW GROUND DATA
ELM variable indices    new indices
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
|nlevgrnd|   bottom layer  |nlevgrnd - 1|



Special zi mapping
ELM                      new indices
      zi                 zi
grid                        grid
 --   -5                  0  ---
|-4|                         |0|
 --   -4                  1  ---
|-3|                         |1|
 --   -3                  2  ---
|-2|                         |2|
 --   -2                  3  ---
|-1|                         |3|
 --   -1                  4  ---  [nlevsno - 1]
| 0|                         |4|  
-----  0 ground interface 5 ----- [nlevsno]
| 1|                         |5|  
 --    1                  6  ---  [nlevsno + 1]
 