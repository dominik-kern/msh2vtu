#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 10:22:22 2021

@author: dominik

This script is derived from the gmsh-demo   terrain_stl.py
Requirements: relief.stl (possibly converted from relief.grd)
It is assumed that curves and points of read-in boundary are ordered.

TODO 
    - use transfinite curves/surfaces/volumes to create a hex-mesh
    - make class for SIDES
    - try to detect side points automatically from stl
"""
import gmsh
import math
import os
import sys
import numpy as np

dim1=1
dim2=2
dim3=3
EPS = 1e-12   # for collinearity check, 

# side definitions (straight lines), points chosen outside to prevent collocation
side_points = {"front": {"p1": np.array([10, -1]), "p2": np.array([10, -2])}, 
               "right": {"p1": np.array([-1, 10]), "p2": np.array([-2, 10])}, 
               "back":  {"p1": np.array([ 0, -1]), "p2": np.array([ 0, -2])}, 
               "left":  {"p1": np.array([-1,  0]), "p2": np.array([-2,  0])}}   # a line is defined by two points: p1 (x,y), p2 (x,y)

side_surface_ids = {"front": [], 
                    "right": [], 
                    "back":  [], 
                    "left":  []}   


def collinear2D(p0, p1, p2):   
    x1, y1 = p1[0] - p0[0], p1[1] - p0[1]
    x2, y2 = p2[0] - p0[0], p2[1] - p0[1]
    return abs(x1 * y2 - x2 * y1) < EPS

def on_line2D(xyz, guess):
    p0 = np.array([xyz[0], xyz[1]])
    points = side_points[guess]
    if collinear2D(p0, points["p1"], points["p2"]):
        return guess
    else:
        for side, points in side_points.items():
            if collinear2D(p0, points["p1"], points["p2"]):
                return side
    return "point not on one of the give sides"


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
z = 7   # bottom_level
bottom_point_ids = []
side_curve_ids = []   # vertical lines
for top_point in top_points:
    xyz = gmsh.model.getValue(0, top_point[1], [])   # get x,y,z coordinates   
    bottom_point_id = gmsh.model.geo.addPoint(xyz[0], xyz[1], z)
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

gmsh.option.setNumber('Mesh.MeshSizeMin', 0.1)
gmsh.option.setNumber('Mesh.MeshSizeMax', 1.0)
gmsh.model.mesh.generate(3)
gmsh.write('gmsh_terrain.msh')

#gmsh.fltk.run()   # GUI
gmsh.finalize()   # to remove all objects from memory
