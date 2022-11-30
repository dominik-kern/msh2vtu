[![No Maintenance Intended](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)

Development continues at OGStools(https://gitlab.opengeosys.org/ogs/tools/ogstools) since November 2022!

# msh2vtu

This script depends on [meshio](https://github.com/nschloe/meshio).
It was tested with meshio 5.3.4 [Python 3.10.6] and gmsh 4.10.5.

Supported element types

- lines (linear and quadratic) in 1D
- triangles and quadrilaterals (linear and quadratic) in 2D
- tetra- and hexahedrons (linear and quadratic) in 3D

## Usage

```
pip install msh2vtu
msh2vtu [-h] [-g] [-r] [-a] [-d DIM] [-o OUTPUT] [-z] [-s] [-v] filename

Prepares a Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-submeshes
and saves them in vtu-format. Note that all mesh entities should belong to physical groups.

positional arguments:
  filename              Gmsh mesh file (*.msh) as input data

optional arguments:
  -h, --help            show this help message and exit
  -g, --ogs             rename "gmsh:physical" to "MaterialIDs" for domains and change type of
                        corresponding cell data to INT32
  -r, --rdcd            renumber domain cell data, physical IDs (cell data) of domains get
                        numbered beginning with zero
  -a, --ascii           save output files (*.vtu) in ascii format
  -d DIM, --dim DIM     spatial dimension (1, 2 or 3), trying automatic detection, if not given
  -o OUTPUT, --output OUTPUT
                        basename of output files; if not given, then it defaults to basename of
                        inputfile
  -z, --delz            deleting z-coordinate, for 2D-meshes with z=0 (dimension must be dim=2)
  -s, --swapxy          swap x and y coordinate
  -v, --version         show program's version number and exit

```
In addition it may be used as Python module with an emulated command line call
```
from msh2vtu import run   # to run mesh conversion
import sys   # to emulate command line call
import argparse   # to parse emulated command line call
parser = argparse.ArgumentParser()

# generate a mesh, e.g. my_mesh.msh with Gmsh

args = argparse.Namespace(filename='my_mesh.msh', output='', dim=0, delz=False, swapxy=False, rdcd=True, ogs=True, ascii=False)   # filename, output="", dim=0, delz, swapxy, rdcd, ogs, ascii
run(args)
```

## Examples

A geological model (2D) of a sediment basin by Christian Silbermann and a terrain model (3D) from the official Gmsh tutorials (x2).

``msh2vtu example/geolayers_2d.msh`` generates from the input file *geolayers_2d.msh* (gmsh 4.4.1):

- *geolayers_2d_boundary.vtu*
- *geolayers_2d_domain.vtu*
- *geolayers_2d_physical_group_RockBed.vtu*
- *geolayers_2d_physical_group_SedimentLayer1.vtu*
- *geolayers_2d_physical_group_SedimentLayer2.vtu*
- *geolayers_2d_physical_group_SedimentLayer3.vtu*
- *geolayers_2d_physical_group_Bottom.vtu*  
- *geolayers_2d_physical_group_Left.vtu*
- *geolayers_2d_physical_group_Right.vtu*
- *geolayers_2d_physical_group_Top.vtu*

## Build instructions

Local development:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
msh2vtu --help
```

Distribution:

```bash
# python -m venv .venv
# source .venv/bin/activate
pip install --upgrade build twine
python -m build
python -m twine upload [--repository testpypi] dist/*
```
