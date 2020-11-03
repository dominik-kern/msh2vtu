# TODO extension to 3D

import meshio
import os
import numpy 
import argparse
import warnings


# parsing command line arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Prepare Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-meshes and save them in vtu-format. Note that all mesh entities must belong to a physical group!', epilog='No cell data are written for boundaries (lines).')
parser.add_argument('filename', help='Gmsh mesh file (*.msh) as input data.')
parser.add_argument('--renumber', action='store_true',help='Renumber physical IDs of domains starting by zero (boundary IDs are not changed).')
parser.add_argument('--rename', action='store_true', help='Rename "gmsh:physical" to "MaterialIDs for domains".')
parser.add_argument('-a','--ascii', action='store_true', help='Save output files (*.vtu) in ascii format.')
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
if args.output=='':	# no parameter given, use same basename as input file
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
type_index=0	# to access type ('line' or 'triangle') of a cellblock in cells
data_index=1	# to access data (nodes of elements) of a cellblock in cells
ph_index=0	# to access physical id in field data
geo_index=1	# to access geometrical id in field data
line_id=1	# geometry type 
triangle_id=2	# geometry type
gmshdict={line_id: 'line', triangle_id: 'triangle'}	# gmsh convention
gmsh_string='gmsh:physical'
ogs_string='MaterialIDs'
# an index field to access the cells corresponding to a physical group
physical_cell_data=cell_data[gmsh_string]	# array of ph_no for all elements in all cellblocks

# if user wants to change 'gmsh:physical" to "MaterialIDs" or not
if args.rename:
	#boundary_data_string=ogs_string # do not write MaterialID for boundaries
	domain_data_string=ogs_string
	selection_data_string=ogs_string	# only used for domains
else:
	#boundary_data_string=gmsh_string
	domain_data_string=gmsh_string
	selection_data_string=gmsh_string	# only used for domains

# if user wants physical group numbering beginning with zero, 
# because OGS wants MaterialIDs to start by zero
id_offset=0	# initial value, zero will not change anything
if args.renumber:
# find minimum physical_id of domains (triangles)
	id_list_domains=[]
	for dataset in field_data.values():	# go through all physical groups
		if dataset[geo_index]==triangle_id:	# only for domains (triangles), ignore boundary (lines)
			id_list_domains.append(dataset[0])  # append physical id
	if len(id_list_domains):	# if there are some domain..
		id_offset=min(id_list_domains)	# ..then find minimal physical id
	

# A mesh consists of cellblocks, now we go through them
# first for domain and boundary, easily recognizable by celltype
domain_cells=[] 	# start with empty list
domain_cell_data=[]	# start with empty list
boundary_cells=[]	# start with empty list
#boundary_cell_data=[] 	# not written (may cause trouble)
for cellblock_data, cellblock in zip(physical_cell_data, cells):

	if cellblock[type_index]==gmshdict[line_id]:
		boundary_cells.append(cellblock)
		#boundary_cell_data.append(cellblock_data)

	if cellblock[type_index]==gmshdict[triangle_id]:
		domain_cells.append(cellblock)
		domain_cell_data.append(cellblock_data-id_offset)

# write results to file
if len(domain_cells):  # only if there are data to write
	domain_mesh=meshio.Mesh(points=points, cells=domain_cells, cell_data={domain_data_string: domain_cell_data}) 
	#domain_mesh.prune()	# get rid of out-of-mesh nodes RESULTS IN MESH DEFECTS
	meshio.write(output_basename+"_domain.vtu", domain_mesh, binary=not args.ascii)

if len(boundary_cells): # only if there are data to write
	boundary_mesh=meshio.Mesh(points=points, cells=boundary_cells)
	#boundary_mesh.prune() # get rid of out-of-mesh nodes GIVES AN ERROR
	meshio.write(output_basename+"_boundary.vtu", boundary_mesh, binary=not args.ascii)	# TODO out-of-mesh nodes


# Now we want to extract subdomains given by physical groups in gmsh,
# so we need an additional loop. 

# name=user-defined name of physical group, data=[physical_id, geometry_type]
for name, data in field_data.items():
	
	selected_cells=[] 	
	selected_cell_data=[]
	ph_id=data[ph_index]	# selection by physical id
	cell_type=gmshdict[data[geo_index]]	# 'line' or 'triangle'

	for cellblock_data, cellblock in zip(physical_cell_data, cells):

		# access matching data by an index field
		selected_data=numpy.array(cellblock[data_index][cellblock_data==ph_id]) # select cells by physical id
		selected_cellblock_data=cellblock_data[cellblock_data==ph_id] # also select corresponding cell data
		if len(selected_data):	# append only nonzero-data
			selected_cells.append(meshio.CellBlock(cell_type, selected_data))
			selected_cell_data.append(selected_cellblock_data)
	
	if len(selected_cells):
		outputfilename=output_basename+"_physical_group_"+name+".vtu"	
		if data[geo_index]==line_id: 	# manual pruning of orphaned nodes
			old_cells=numpy.concatenate([selected_cells[k][1] for k in range(len(selected_cells)) ])	# join active cells from all cellblocks 
			shape2d=old_cells.shape
			old_points=old_cells.flatten()	# "old" means from the input mesh and "new" appearing in the actual physical group			
			unique_points, unique_inverse = numpy.unique(old_points, return_inverse=True)
			new_points=points[unique_points]	# extract only used nodes
			new_cells=unique_inverse.reshape(shape2d)+1	# cell connectivity corresponding to new node numbering (starting with 1)				
			physical_submesh=meshio.Mesh(points=new_points, cells=[meshio.CellBlock(cell_type, new_cells)])
		else:
			physical_submesh=meshio.Mesh( points=points, cells=selected_cells, 
			cell_data={selection_data_string: selected_cell_data} ) 
			physical_submesh.prune()  # somehow here prune works correctly (if it is only one cellblock?)
		meshio.write(outputfilename, physical_submesh, binary=not args.ascii)

