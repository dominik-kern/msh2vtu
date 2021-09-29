#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 15:56:56 2021
@author: dominik

Changed from Sphere (https://tetgen.pyvista.org basic example) to Box
"""
import pyvista as pv
import tetgen
import numpy as np
pv.set_plot_theme('document')

box = pv.Box(bounds=(-1., 1., -1., 1., -1., 1.), level=2, quads=False)
tet = tetgen.TetGen(box)
#box.plot(show_edges=True, opacity=0.4)
tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
grid = tet.grid
grid.plot(show_edges=True, opacity=0.4)


# get cell centroids
cell_center = grid.cell_centers().points

# extract cells below the xy plane (z=0)
mask = cell_center[:, 2] < 0
cell_ind = mask.nonzero()[0]
subgrid = grid.extract_cells(cell_ind)

# advanced plotting
plotter = pv.Plotter()
plotter.add_mesh(subgrid, 'lightgrey', lighting=True, show_edges=True)
plotter.add_mesh(box, 'r', 'wireframe')
plotter.add_legend([[' Input Mesh ', 'r'],
                    [' Tessellated Mesh ', 'black']])
plotter.show()

cell_qual = subgrid.compute_cell_quality()['CellQuality']
print(f'Mean cell quality: {cell_qual.mean():.3}')

# plot quality
subgrid.plot(scalars=cell_qual, scalar_bar_args={'title': 'Quality'},
              cmap='bwr', clim=[0, 1], flip_scalars=True,
              show_edges=True,)
