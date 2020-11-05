# msh2vtu

This script depends on [meshio](https://github.com/nschloe/meshio).
It was tested with meshio 4.3.3 [Python 3.8.5]

Note also, currently restricted to 2D meshes!

**Usage:**
```
usage: msh2vtu.py [-h] [--renumber] [--ogs] [-a] [-o OUTPUT] filename

Prepares a Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-meshes and saves them in
vtu-format. Cell data are only written for domains but not for boundaries. Note that all mesh entities must belong
to physical groups!

positional arguments:
  filename              Gmsh mesh file (*.msh) as input data

optional arguments:
  -h, --help            show this help message and exit
  --renumber            renumber physical IDs of domains starting by zero (boundaries are ignored)
  --ogs                 rename "gmsh:physical" to "MaterialIDs" for domains and change type of corresponding cell
                        data to INT32
  -a, --ascii           save output files (*.vtu) in ascii format
  -o OUTPUT, --output OUTPUT
                        base name of output files; if not given, then it defaults to basename of inputfile

Tested with Meshio 4.3.3 and Gmsh 4.4.1. Check for changes between versions, if there are errors!
```

**Example:**
A geological model of a sediment basin by Christian Silbermann

From the input file *geolayers_2d.msh* (gmsh 4.4.1) are generated
``python3 msh2vtu example/geolayers_2d.msh``:

- *geolayers_2d_boundary.vtu*,
- *geolayers_2d_domain.vtu*,                 
- *geolayers_2d_physical_group_RockBed.vtu*,
- *geolayers_2d_physical_group_SedimentLayer1.vtu*,
- *geolayers_2d_physical_group_SedimentLayer2.vtu*,
- *geolayers_2d_physical_group_Bottom.vtu*,  
- *geolayers_2d_physical_group_SedimentLayer3.vtu*,
- *geolayers_2d_physical_group_Left.vtu*,    
- *geolayers_2d_physical_group_Top.vtu*,
- *geolayers_2d_physical_group_Right.vtu*.

