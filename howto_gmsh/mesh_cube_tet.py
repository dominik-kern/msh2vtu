# mesh unit cube with tetraeders
import numpy
import gmsh

# init
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("cube")

dim1=1
dim2=2
dim3=3
lc = 1.0    # mesh size

# vertices
gmsh.model.geo.addPoint(1, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 1, 0, lc, 2)
gmsh.model.geo.addPoint(0, 1, 0, lc, 3)
gmsh.model.geo.addPoint(0, 0, 1, lc, 4)
gmsh.model.geo.addPoint(1, 0, 1, lc, 5)
gmsh.model.geo.addPoint(1, 1, 1, lc, 6)
gmsh.model.geo.addPoint(0, 1, 1, lc, 7)
gmsh.model.geo.addPoint(0, 0, 0, lc, 8)

# edges
gmsh.model.geo.addLine(7, 6, 1)
gmsh.model.geo.addLine(6, 5, 2)
gmsh.model.geo.addLine(5, 1, 3)
gmsh.model.geo.addLine(1, 8, 4)
gmsh.model.geo.addLine(8, 3, 5)
gmsh.model.geo.addLine(3, 7, 6)
gmsh.model.geo.addLine(7, 4, 7)
gmsh.model.geo.addLine(4, 8, 8)
gmsh.model.geo.addLine(4, 5, 9)
gmsh.model.geo.addLine(2, 1, 10)
gmsh.model.geo.addLine(2, 6, 11)
gmsh.model.geo.addLine(2, 3, 12)

# faces
gmsh.model.geo.addCurveLoop([6, 1, -11, 12], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.addCurveLoop([11, 2, 3, -10], 2)
gmsh.model.geo.addPlaneSurface([2], 2)

gmsh.model.geo.addCurveLoop([2, -9, -7, 1], 3)
gmsh.model.geo.addPlaneSurface([3], -3)

gmsh.model.geo.addCurveLoop([6, 7, 8, 5], 4)
gmsh.model.geo.addPlaneSurface([4], -4)

gmsh.model.geo.addCurveLoop([8, -4, -3, -9], 5)
gmsh.model.geo.addPlaneSurface([5], 5)

gmsh.model.geo.addCurveLoop([10, 4, 5, -12], 6)
gmsh.model.geo.addPlaneSurface([6], 6)

# volume
gmsh.model.geo.addSurfaceLoop([6,2,1,4,3,5], 1)
gmsh.model.geo.addVolume([1], 1)

# physical groups
DCBF = gmsh.model.addPhysicalGroup(dim2, [4,3,2,6])
gmsh.model.setPhysicalName(dim2, DCBF, "Tunnel")

A = gmsh.model.addPhysicalGroup(dim2, [1])
gmsh.model.setPhysicalName(dim2, A, "Ausgang")

E = gmsh.model.addPhysicalGroup(dim2, [5])
gmsh.model.setPhysicalName(dim2, E, "Eingang")

W = gmsh.model.addPhysicalGroup(dim3, [1])
gmsh.model.setPhysicalName(dim3, E, "Wuerfel")

# mesh
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(dim3)
gmsh.model.mesh.setOrder(2)   # higher order, no difference between Lagrange and Serendipity elements

gmsh.write("cube_tet.msh")
gmsh.finalize()
