# mesh unit square with triangle elements (higher order)
import numpy
import gmsh

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("line")

# Dimensions
dim0=0
dim1=1

lc=0.1 # characteristic length for mesh size

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 0, 0, lc, 2)

# Lines connecting points
gmsh.model.geo.addLine(1, 2, 1)

# Here we define physical curves that groups
Left = gmsh.model.addPhysicalGroup(dim0, [1])
gmsh.model.setPhysicalName(dim0, Left, "left")

Right = gmsh.model.addPhysicalGroup(dim0, [2])
gmsh.model.setPhysicalName(dim0, Right, "right")

Domain = gmsh.model.addPhysicalGroup(dim1, [1])
gmsh.model.setPhysicalName(dim1, Domain, "domain")

# Before it can be meshed, the internal CAD representation must be synchronized
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(dim1)
#gmsh.model.mesh.setOrder(2)   # higher order, for simplex elements there is no difference between Lagrange and Serendipity
gmsh.write("line.msh")

gmsh.finalize()
