# mesh unit square with triangle elements and run msh2vtu directly
import numpy   # for numerics
import gmsh   # for meshing
from msh2vtu import run   # for mesh conversion
import sys   # to emulate command line call
import argparse   # to parse emulated command line call
parser = argparse.ArgumentParser()

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("square")

# Dimensions
dim1=1
dim2=2

lc=0.5 # characteristic length for mesh size

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
gmsh.model.geo.addPoint(0, 1, 0, lc, 4)

# Lines connecting points
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

# The third elementary entity is the surface. In order to define a surface 
# from the curves defined above, a curve loop has first to be defined.
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

# Add plane surfaces defined by one or more curve loops.
gmsh.model.geo.addPlaneSurface([1], 1)

# Here we define physical curves that groups
Bottom = gmsh.model.addPhysicalGroup(dim1, [1])
gmsh.model.setPhysicalName(dim1, Bottom, "Bottom")

Right = gmsh.model.addPhysicalGroup(dim1, [2])
gmsh.model.setPhysicalName(dim1, Right, "Right")

Top = gmsh.model.addPhysicalGroup(dim1, [3])
gmsh.model.setPhysicalName(dim1, Top, "Top")

Left = gmsh.model.addPhysicalGroup(dim1, [4])
gmsh.model.setPhysicalName(dim1, Left, "Left")

Rectangle = gmsh.model.addPhysicalGroup(dim2, [1])
gmsh.model.setPhysicalName(dim2, Rectangle, "UnitSquare")

# Before it can be meshed, the internal CAD representation must be synchronized
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(dim2)
#gmsh.model.mesh.setOrder(2)   # higher order, for simplex elements there is no difference between Lagrange and Serendipity
gmsh.write("square_tri.msh")   # if meshio could directly access a gmsh object then this intermediate file could be skipped
gmsh.finalize()

# emulate command line and run msh2vtu
args = argparse.Namespace(filename='square_tri.msh', output='', dim=0, delz=False, swapxy=False, rdcd=True, ogs=True, ascii=False)   # filename, output="", dim=0, delz, swapxy, rdcd, ogs, ascii
run(args)
