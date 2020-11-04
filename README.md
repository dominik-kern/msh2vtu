# msh2vtu

Note, currently restricted to 2D meshes!

**Usage:**
```
usage: msh2vtu.py [-h] [--renumber] [--ogs] [-a] [-o OUTPUT] filename

Prepare Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-meshes and save them in vtu-format. Note that all mesh entities must belong to a physical group!

positional arguments:
  filename              Gmsh mesh file (*.msh) as input data.

optional arguments:
  -h, --help            show this help message and exit
  --renumber            Renumber physical IDs of domains starting by zero (boundary IDs are ignored).
  --ogs                 Rename "gmsh:physical" to "MaterialIDs" for domains and change type of corresponding cell data to INT32.
  -a, --ascii           Save output files (*.vtu) in ascii format.
  -o OUTPUT, --output OUTPUT
                        Base name of output files; if not given, then it defaults to basename of inputfile.

No cell data are written for boundaries (lines).
```

**Example:**
A geological model of a sediment basin by Christian Silbermann

From the input file *geolayers_2d.msh* are generated:    
-*geolayers_2d_boundary.vtu*,
-*geolayers_2d_domain.vtu*,                 
-*geolayers_2d_physical_group_RockBed.vtu*,
-*geolayers_2d_physical_group_SedimentLayer1.vtu*,
-*geolayers_2d_physical_group_SedimentLayer2.vtu*,
-*geolayers_2d_physical_group_Bottom.vtu*,  
-*geolayers_2d_physical_group_SedimentLayer3.vtu*,
-*geolayers_2d_physical_group_Left.vtu*,    
-*geolayers_2d_physical_group_Top.vtu*,
-*geolayers_2d_physical_group_Right.vtu*.

