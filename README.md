# msh2vtu

Note, currently restricted to 2D meshes!

*Usage:*
python3 msh2vtu.py [-h] [--renumber] [--rename] [-o OUTPUT] filename

Prepare Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-meshes and
save them in vtu-format.

positional arguments:
  filename              Gmsh mesh (*.msh)

optional arguments:
  -h, --help            show this help message and exit
  --renumber            Renumber physical IDs of domains starting by zero (boundary IDs are not
                        changed)
  --rename              Rename "gmsh:physical" to "MaterialIDs"
  -o OUTPUT, --output OUTPUT
                        Base name of output files; if not given, then it defaults to basename of
                        inputfile.

*Example:*
A geological model of a sediment basin by Christian Silbermann
