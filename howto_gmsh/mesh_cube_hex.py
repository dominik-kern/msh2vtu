# mesh unit cube with hexaeders
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

# vertices (dim=0)
gmsh.model.geo.addPoint(1, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 1, 0, lc, 2)
gmsh.model.geo.addPoint(0, 1, 0, lc, 3)
gmsh.model.geo.addPoint(0, 0, 0, lc, 4)
# edges (dim=1)
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

# faces (dim=2)
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# prepare structured grid
nx = 3   # number of mesh divisions in x-direction
ny = 3   # .. and y
gmsh.model.geo.mesh.setTransfiniteCurve(1, ny)
gmsh.model.geo.mesh.setTransfiniteCurve(2, nx)
gmsh.model.geo.mesh.setTransfiniteCurve(3, ny)
gmsh.model.geo.mesh.setTransfiniteCurve(4, nx)
gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.geo.mesh.setRecombine(dim2, 1)

# volume (dim=3)
# extrusion of a surface creates: top surface, volume, side surfaces
# parameters: (dim, tag), x,y,z, divisions per layer, heights per layer, Recombine
newEntities=gmsh.model.geo.extrude([(dim2, 1)], 0, 0, 0.5, [2], [], True)
top_tag=newEntities[0][1]
vol_tag=newEntities[1][1]
side_tags=[nE[1] for nE in newEntities[2:] ]

Bottom = gmsh.model.addPhysicalGroup(dim2, [1])
gmsh.model.setPhysicalName(dim2, Bottom, "bottom")

Top = gmsh.model.addPhysicalGroup(dim2, [top_tag])
gmsh.model.setPhysicalName(dim2, Top, "top")

Sides = gmsh.model.addPhysicalGroup(dim2, side_tags)
gmsh.model.setPhysicalName(dim2, Sides, "sides")

Vol = gmsh.model.addPhysicalGroup(dim3, [vol_tag])
gmsh.model.setPhysicalName(dim3, Vol, "volume")


# mesh
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(dim3)
gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 1) # serendipity elements
gmsh.model.mesh.setOrder(2)   # higher order elements (quadratic approximation)

gmsh.write("cube_hex.msh")
gmsh.finalize()
