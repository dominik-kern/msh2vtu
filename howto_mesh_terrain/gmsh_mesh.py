#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 10:22:22 2021

@author: dominik

This script is derived from the gmsh-demo   terrain_stl.py
Requirements: relief.stl (possibly converted from relief.grd)
It is assumed that curves and points of read-in boundary are ordered.

TODO 
    - N-E-S-W instead of front/left/...
    - use transfinite curves/surfaces/volumes to create a hex-mesh
    - make class for SIDES
    - try to detect side points automatically from stl
"""
import gmsh
import math
import os
import sys
import numpy as np

dim1 = 1
dim2 = 2
dim3 = 3
EPS = 1e-6   # for collinearity check, 
MIN_SIZE = 0.1 # minimum element size
MAX_SIZE = 1.0 # maximum elemenz size

# side definitions (straight lines), points chosen outside to prevent collocation
X0 =  0 
X1 = 10
DX =  1
Y0 =  0
Y1 = 10
DY =  1
Z = 7   # bottom_level

side_points = {"front": {"p1": np.array([X1, Y0-DY]), "p2": np.array([X1, Y1+DY])}, 
               "right": {"p1": np.array([X0-DX, Y1]), "p2": np.array([X1+DX, Y1])}, 
               "back":  {"p1": np.array([X0, Y0-DY]), "p2": np.array([X0, Y1+DY])}, 
               "left":  {"p1": np.array([X0-DX, Y0]), "p2": np.array([X1+DX, Y0])}}   # a line is defined by two points: p1 (x,y), p2 (x,y)

side_surface_ids = {"front": [], 
                    "right": [], 
                    "back":  [], 
                    "left":  []}   


def collinear2D(p0, p1, p2):   #
    x1, y1 = p0[0] - p1[0], p0[1] - p1[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    CROSS_PRODUCT = x1 * y2 - x2 * y1
    print("p0-p1=[{}, {}]".format(x1,y1))
    print("p2-p0=[{}, {}]".format(x2,y2))
    print(CROSS_PRODUCT)
    return abs(CROSS_PRODUCT) < EPS

def on_line2D(xyz, guess):
    p0 = np.array([xyz[0], xyz[1]])
    print(p0)
    points = side_points[guess]
    if collinear2D(p0, points["p1"], points["p2"]):
        return guess
    else:
        for side, points in side_points.items():
            print(side)
            if collinear2D(p0, points["p1"], points["p2"]):
                return side
    print("point " + str(p0) + " not on a given side")            
    return "NO SIDE FOUND" 


gmsh.initialize(sys.argv) # use finalize to unload from memory

# load an STL surface
path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, 'relief.stl'))

# classify the surface mesh according to given angle, and create discrete model
# entities (surfaces, curves and points) accordingly; curveAngle forces bounding
# curves to be split on sharp corners
gmsh.model.mesh.classifySurfaces(math.pi, curveAngle = 0.0*math.pi)   # angle, boundary = True, forReparametrization = False, curveAngle = pi
# angle=pi, selects the surface as one, no matter what angles are between the STL-patches
# curveAngle=0 selects each STL line segment as curve, even if they continue in the same direction

# create a geometry for the discrete curves and surfaces
gmsh.model.mesh.createGeometry()   

# retrieve the surface, its boundary curves and corner points
top_surfaces = gmsh.model.getEntities(2)   # discrete surface 
if len(top_surfaces) == 1:
    top_surface = top_surfaces[0]
else:
    print("More than one top surface detected.")
    gmsh.finalize()
    sys.exit()

top_curves = gmsh.model.getEntities(1)   # discrete surface 
top_points = gmsh.model.getEntities(0)   # discrete surface 

## create geometric entities to form one volume below the terrain surface
bottom_point_ids = []
side_curve_ids = []   # vertical lines
for top_point in top_points:
    xyz = gmsh.model.getValue(0, top_point[1], [])   # get x,y,z coordinates   
    bottom_point_id = gmsh.model.geo.addPoint(xyz[0], xyz[1], Z)
    bottom_point_ids.append(bottom_point_id)
    side_curve_id = gmsh.model.geo.addLine(bottom_point_id, top_point[1])
    side_curve_ids.append(side_curve_id)
gmsh.model.geo.synchronize()

Nc = len(top_curves)
bottom_curve_ids = []   # horizontal lines

guessed_side = "left"
for ci, top_curve in enumerate(top_curves):  
    cip1 = (ci+1) % Nc   # switch to next and from last to first (cycle)
    xyz_i = gmsh.model.getValue(0, bottom_point_ids[ci], [])   # get x,y,z coordinates 
    xyz_ip1 = gmsh.model.getValue(0, bottom_point_ids[cip1], [])   # get x,y,z coordinates 
    xyz = 0.5*(xyz_i+xyz_ip1)   # midpoint
    detected_side = on_line2D(xyz, guessed_side)
    if detected_side not in side_points:
        print("Geometry error")
        gmsh.finalize()
        sys.exit()
    else:
        guessed_side =  detected_side   # next guess
    
    bottom_curve_id = gmsh.model.geo.addLine(bottom_point_ids[ci], bottom_point_ids[cip1])
    bottom_curve_ids.append(bottom_curve_id)
    side_ll = gmsh.model.geo.addCurveLoop([bottom_curve_id, side_curve_ids[cip1], -top_curve[1], -side_curve_ids[ci]])
    side_surface_id = gmsh.model.geo.addPlaneSurface([side_ll])
    side_surface_ids[detected_side].append(side_surface_id)
    
bottom_ll = gmsh.model.geo.addCurveLoop(bottom_curve_ids)
bottom_surface_id = gmsh.model.geo.addPlaneSurface([bottom_ll])

all_side_ids=[]
for surface_ids in side_surface_ids.values():
    all_side_ids += surface_ids

surface_loop1 = gmsh.model.geo.addSurfaceLoop( [top_surface[1]]+ all_side_ids + [bottom_surface_id] )
volume1 = gmsh.model.geo.addVolume([surface_loop1])

gmsh.model.geo.synchronize()

# physical groups
Top_Surface = gmsh.model.addPhysicalGroup(dim2, [top_surface[1]])
gmsh.model.setPhysicalName(dim2, Top_Surface, "top")

Bottom_Surface = gmsh.model.addPhysicalGroup(dim2, [bottom_surface_id])
gmsh.model.setPhysicalName(dim2, Bottom_Surface, "bottom")

Side_Surfaces = []
for side in side_surface_ids:
     Side_Surface = gmsh.model.addPhysicalGroup(dim2, side_surface_ids[side])
     gmsh.model.setPhysicalName(dim2, Side_Surface, side)

Volume1 = gmsh.model.addPhysicalGroup(dim3, [volume1])
gmsh.model.setPhysicalName(dim3, Volume1, "volume")

# meshing

gmsh.option.setNumber('Mesh.MeshSizeMin', MIN_SIZE)
gmsh.option.setNumber('Mesh.MeshSizeMax', MAX_SIZE)
gmsh.model.mesh.generate(3)
gmsh.write('gmsh_terrain.msh')

#gmsh.fltk.run()   # GUI
gmsh.finalize()   # to remove all objects from memory
