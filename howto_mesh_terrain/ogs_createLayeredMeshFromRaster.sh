#!/usr/bin/bash
generateStructuredMesh -e tri --lx 10 --nx 20 --ly 10 --ny 20 -o tri2d_mesh.vtu
createLayeredMeshFromRasters -i tri2d_mesh.vtu -o wedge3d_mesh.vtu -r layer_file_list
Vtu2Grid -i wedge3d_mesh.vtu -o hex3d_mesh.vtu -x 0.4 
