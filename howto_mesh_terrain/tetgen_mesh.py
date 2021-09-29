#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 09:29:51 2021

@author: dominik

TetGen starts from a closed triangulated surface mesh.
So one needs to find a way to create such a surface (pyvista).

FAVORIZED STRATEGY:
    - create a pv.Box, adjust its top surface z=z_top to DEM z=z(x,y)

ALTERNATIVES 
    - create complete Polydata (vertices, faces) for DEM and remaining sides 
    - generate DEM for all sides and merge them
    - clip a volume
"""
import pyvista as pv
import tetgen
import numpy as np
pv.set_plot_theme('document')

box = pv.Box(bounds=(-1., 1., -1., 1., -1., 1.), level=1, quads=False)
#box.plot(show_edges=True, opacity=0.4)

tet = tetgen.TetGen(box)
tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
grid = tet.grid
#grid.plot(show_edges=True, opacity=0.4)

# plot mesh quality
cell_qual = grid.compute_cell_quality()['CellQuality']
grid.plot(scalars=cell_qual, scalar_bar_args={'title': 'Quality'},
              cmap='bwr', clim=[0, 1], flip_scalars=True,
              show_edges=True,)

