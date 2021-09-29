#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 11:31:50 2021

@author: dominik

generate structured grid

TODO
    add physical group as cell data
"""

import pyvista as pv
from PVGeo.grids import SurferGridReader
import numpy as np

#   read file
filename = "relief.grd"
dem = SurferGridReader().apply(filename)

relief = dem.warp_by_scalar(scale_factor=1)
relief.plot(cmap='terrain')   # 

# create meshgrid
z_spacing = np.array([1, 2, 3])   
 
xx = np.repeat(relief.x, len(z_spacing), axis=-1)
yy = np.repeat(relief.y, len(z_spacing), axis=-1)
z_offset = np.cumsum(z_spacing).reshape((1, 1, -1))
# since the top z-coordinates (relief) repeat, we must subtract the z-offset for each layer
zz = np.repeat(relief.z, len(z_spacing), axis=-1) - z_offset

mesh = pv.StructuredGrid(xx, yy, zz)
mesh["Elevation"] = zz.ravel(order="F")   # flatten nested array to 1D-array

#List of camera position, focal point, and view up, either vector or string, e.g. "xy"
#cpos = [(1826736.796308761, 5655837.275274233, 4676.8405505181745),
# (1821066.1790519988, 5649248.765538796, 943.0995128226014),
# (-0.2797856225380979, -0.27966946337594883, 0.9184252809434081)]

# as "Elevation is the only data field it gets plotted
mesh.plot(show_edges=True)   # lighting=False

#mesh.save("terrain.vtk")   
pv.save_meshio("pv_terrain.vtu", mesh)