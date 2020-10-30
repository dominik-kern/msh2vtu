'''
read gmsh file and write vtu-files for domain, boundaries and possibly subdomains
according to physical groups in gmsh 
(all mesh entities must belong to a physical group!).
'''

# TODO extension to 3D

import meshio
import os
from numpy import array
import argparse
import warnings


# parsing command line arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Prepare Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-meshes and save them in vtu-format.')
parser.add_argument('filename', help='Gmsh mesh (*.msh)')
parser.add_argument('--renumber', action='store_true',help='Renumber physical IDs of domains starting by zero (boundary IDs are not changed)')
parser.add_argument('--rename', action='store_true', help='Rename "gmsh:physical" to "MaterialIDs"')
parser.add_argument('-o','--output', default='', help='Base name of output files; if not given, then it defaults to basename of inputfile.')
args = parser.parse_args()


# check if input file exists and is in gmsh-format
if os.path.isfile(args.filename):
	filename_without_extension=os.path.splitext(args.filename)[0]
	file_extension=os.path.splitext(args.filename)[1]
	if file_extension != '.msh':
		warnings.warn("Warning, input file seems not to be in gmsh-format (*.msh)")
else:
	raise FileNotFoundError


# derive output filenames
if args.output=='':	# no parameter given, use same as input file
	output_basename=filename_without_extension
else:
	output_basename=args.output


# read in mesh 
mesh = meshio.read(args.filename)
# points: 2Darray( point, xyz)
# cells: list( cellblocks( type='line'/'triangle', data=2Darray(element, points) ))
# cell_data: dict('gmsh:physical', 2Darray(geo_tag, ph_tag))
# field_data: dict('top'/'bottom'/PH_NAME, 1Darray (ph_tag, geo_type)	
points, cells = mesh.points, mesh.cells
cell_data, field_data = mesh.cell_data, mesh.field_data


# some variable declarations
type_index=0	# to access type ('line' or 'triangle') of a cellblock
data_index=1	# to access data (nodes of elements) of a cellblock
gmshdict={1:'line', 2: 'triangle'}	# gmsh convention
gmsh_string='gmsh:physical'
ogs_string='MaterialIDs'
# an index field to access the cells corresponding to a physical group
physical_cell_data=cell_data[gmsh_string]	# array of ph_no for each element 

# if user wants to change 'gmsh:physical" to MaterialIDs or not
if args.rename:
	boundary_data_string=ogs_string
	domain_data_string=ogs_string
	selection_data_string=ogs_string
else:
	boundary_data_string=gmsh_string
	domain_data_string=gmsh_string
	selection_data_string=gmsh_string

# if user wants to set physical group numbering to beginning with zero, 
# because OGS wants MaterialIDs to start by 0
id_offset=0	# initial value, zero will not change anything
if args.renumber:
# find minimum physical_id of domains (triangles)
	id_list_domains=[]
	for dataset in field_data.values():	# go through all physical groups
		if dataset[1]==2:	# only for domains (triangles), ignore boundary (lines)
			id_list_domains.append(dataset[0])  # append physical id
	if len(id_list_domains):
		id_offset=min(id_list_domains)
print(id_offset)	

# A mesh consists of cellblocks, now we go through them
# first for domain and boundary, easily recognizable by celltype
domain_cells=[] 	# start with empty list
domain_cell_data=[]
boundary_cells=[]	# start with empty list
boundary_cell_data=[]
for cellblock_data, cellblock in zip(physical_cell_data, cells):

	if cellblock[type_index]=='line':
		boundary_cells.append(cellblock)
		boundary_cell_data.append(cellblock_data)

	if cellblock[type_index]=='triangle':
		domain_cells.append(cellblock)
		domain_cell_data.append(cellblock_data-id_offset)

# write results to file
if len(domain_cell_data):  # only if there are data to write
	meshio.write(output_basename+"_domain.vtu", 
	meshio.Mesh(points=points, cells=domain_cells, cell_data={domain_data_string: domain_cell_data})) 

if len(boundary_cell_data): # only if there are data to write
	meshio.write(output_basename+"_boundary.vtu", 
	meshio.Mesh(points=points, cells=boundary_cells, cell_data={boundary_data_string: boundary_cell_data}) )


# now we want to extract subdomains given by physical groups in gmsh
# so we need an additional loop 

# name=user-defined name of physical group, data=[physical_id, geometry_type]
for name, data in field_data.items():
	selected_cells=[] 	
	selected_cell_data=[]
	ph_id=data[0]	# selection by physical id
	cell_type=gmshdict[data[1]]	# 'line' or 'triangle'

	for cellblock_data, cellblock in zip(physical_cell_data, cells):

		# access matching data by an index field
		selected_data=array(cellblock[data_index][cellblock_data==ph_id]) #
		if len(selected_data):	# append only nonzero-data
			selected_cells.append(meshio.CellBlock(cell_type, selected_data))
			selected_cell_data.append(cellblock_data)
	
	if len(selected_cell_data):
		outputfilename=output_basename+"_physical_group_"+name+".vtu"	
		meshio.write(outputfilename, 
		meshio.Mesh(points=points, cells=selected_cells, 
		cell_data={selection_data_string: selected_cell_data} )) 

