#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:54:03 2021

@author: dominik
"""
import pyvista as pv
from PVGeo.grids import SurferGridReader

filename = "relief.grd"
dem = SurferGridReader().apply(filename)   # GRD reader
relief = dem.warp_by_scalar(scale_factor=1.)
relief.plot()
pv.save_meshio("relief.stl", relief.triangulate())   # STL writer