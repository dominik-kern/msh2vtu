# msh2vtu

This script depends on [meshio](https://github.com/nschloe/meshio).
It was tested with meshio 4.3.6 [Python 3.8.5] and gmsh 4.4.1.

Note that msh2vtu is still restricted to linear simplex meshes (3-point triangle in 2D, 4-point tetraeder in 3D)!

**Usage:**
```
usage: msh2vtu.py [-h] [-g] [-r] [-a] [-d DIM] [-o OUTPUT] [-z] filename

Prepares a Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-submeshes and saves them in vtu-format. Note that all mesh entities must belong to physical groups.

positional arguments:
  filename              Gmsh mesh file (*.msh) as input data

optional arguments:
  -h, --help            show this help message and exit
  -g, --ogs             rename "gmsh:physical" to "MaterialIDs" for domains and change type of corresponding cell data to INT32
  -r, --rdcd            renumber domain cell data, physical IDs (cell data) of domains get numbered beginning with zero
  -a, --ascii           save output files (*.vtu) in ascii format
  -d DIM, --dim DIM     spatial dimension (2 or 3), trying automatic detection, if not given
  -o OUTPUT, --output OUTPUT
                        basename of output files; if not given, then it defaults to basename of inputfile
  -z, --delz            deleting z-coordinate, for 2D-meshes with z=0

Tested with Meshio 4.3.6 and Gmsh 4.4.1. Check for changes between versions, if there are errors.
```

**Example:**
A geological model of a sediment basin by Christian Silbermann

``python3 msh2vtu example/geolayers_2d.msh`` generates from the input file *geolayers_2d.msh* (gmsh 4.4.1):

- *geolayers_2d_boundary.vtu*,
- *geolayers_2d_domain.vtu*,                 
- *geolayers_2d_physical_group_RockBed.vtu*,
- *geolayers_2d_physical_group_SedimentLayer1.vtu*,
- *geolayers_2d_physical_group_SedimentLayer2.vtu*,
- *geolayers_2d_physical_group_SedimentLayer3.vtu*,
- *geolayers_2d_physical_group_Bottom.vtu*,  
- *geolayers_2d_physical_group_Left.vtu*,    
- *geolayers_2d_physical_group_Right.vtu*.
- *geolayers_2d_physical_group_Top.vtu*,
