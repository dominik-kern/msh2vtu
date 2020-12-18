#!/usr/bin/env python3
# Author: Dominik Kern (TU Bergakademie Freiberg)
# Python version 3.8.5
import meshio  # https://github.com/nschloe/meshio
import os
import sys
import numpy
import argparse
import warnings


# 	obsolete auxiliary function needed for meshio version 4.0.16
# def line_mesh_prune(points, input_cells):  # remove orphaned points
#    # "old" means from the input mesh and "new" the mesh of connected points only
#    original_shape = input_cells.shape  # for reconstrution after flatten
#    old_points = input_cells.flatten()  # 1d-array needed
#    unique_points, unique_inverse = numpy.unique(old_points, return_inverse=True)
#    new_points = points[unique_points]  # extract only used nodes
#    new_cells = unique_inverse.reshape(original_shape)  # update cell connectivity
#    return new_points, new_cells


# print info for mesh: statistics and data field names
def print_info(mesh):
    N, D = mesh.points.shape
    print(str(N) + " points in " + str(D) + " dimensions", end='; ')
    cell_info = "cells: "
    for cell_type, cell_values in mesh.cells_dict.items():
        cell_info += str(len(cell_values)) + " " + cell_type + ", "
    print(cell_info[0:-2], end='; ')
    print("point_data=" + str(list(mesh.point_data)), end='; ')
    print("cell_data=" + str(list(mesh.cell_data)), end='; ')
    print("cell_sets=" + str(list(mesh.cell_sets)))
    print("##")


# function to create node connectivity list, i.e. store for each node to which domain elements it belongs
def cells_at_nodes(cells, node_count, cell_start_index):	# depending on the numbering of mixed meshes in OGS one may think of an object-oriented way to add elements (of different type) to node connectivity
    node_connectivity = [ set() for _ in range(node_count) ]  # initialize list of sets
    cell_index = cell_start_index
    for cell in cells:
        for node in cell:
            node_connectivity[node].add(cell_index)
        cell_index += 1
    return node_connectivity


# function to find out to which domain elements a boundary element belongs 
def connected_domain_cells(boundary_cells_values, node_connectivity):
    number_of_boundary_cells=len(boundary_cells_values)
    domain_cells_array=numpy.zeros(number_of_boundary_cells) 	# to return unique common connected domain cell to be stored as cell_data ("bulk_element_id"), if there are more than one do not store anything
    domain_cells_number=numpy.zeros(number_of_boundary_cells)	# number of connected domain_cells
    for cell_index, cell_values in enumerate(boundary_cells_values): 	# cell lists node of which it is comprised
        connected_domain_cells=[]
        for node in cell_values:
            connected_domain_cells.append(node_connectivity[node])
        common_domain_cells=set.intersection(*connected_domain_cells) 
        number_of_connected_domain_cells=len(common_domain_cells)
        domain_cells_number[cell_index]=number_of_connected_domain_cells
        if number_of_connected_domain_cells == 1:	# there should be one domain cell for each boundary cell, however cells of boundary dimension may be in the domain (e.g. as sources)
            domain_cells_array[cell_index]=common_domain_cells.pop()  # assign only one (unique) connected dmain cell
        elif number_of_connected_domain_cells <1:
            warnings.warn( "connected domain cells: cell " + str(cell_index)  + " of boundary dimension does not belong to any domain cell")
    return domain_cells_array, domain_cells_number


# some variable declarations
ph_index = 0	# to access physical group id in field data
geo_index = 1	# to access geometrical dimension in field data
dim0 = 0
dim1 = 1
dim2 = 2
dim3 = 3
ogs_point_data_key = "bulk_node_ids"	# for all points, as the selection goes via the cells and subsequent trim
available_cell_types = { dim0: {"vertex"}, dim1: {"line"}, dim2: {"triangle", "quad"}, dim3: {"tetra", "wedge"} }  
gmsh_physical_cell_data_key = "gmsh:physical"	
ogs_domain_cell_data_key = "MaterialIDs"
ogs_boundary_cell_data_key = "bulk_elem_ids"

tested_meshio_version = "4.3.6"
tested_gmsh_version = "4.4.1"

if meshio.__version__ < tested_meshio_version:
    warnings.warn(
        "Warning, out-dated Meshio version. In case of errors watch for commented code fragments from previous versions in this script (msh2vtu)."
    )
elif meshio.__version__ > tested_meshio_version:
    print(
        "Newer version of Meshio than supported. Backward compatibility may be missing!"
    )

# parsing command line arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
    description="Prepares a Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-submeshes and saves them in vtu-format. Note that all mesh entities should belong to physical groups.",
    epilog="Tested with Meshio "
    + tested_meshio_version
    + " and Gmsh "
    + tested_gmsh_version
    + ". Check for changes between versions, if there are errors.",
)
parser.add_argument("filename", help="Gmsh mesh file (*.msh) as input data")
parser.add_argument(
    "-g",
    "--ogs",
    action="store_true",
    help='rename "gmsh:physical" to "MaterialIDs" for domains and change type of corresponding cell data to INT32',
)
parser.add_argument(
    "-r",
    "--rdcd",
    action="store_true",
    help="renumber domain cell data, physical IDs (cell data) of domains get numbered beginning with zero",
)
parser.add_argument(
    "-a",
    "--ascii",
    action="store_true",	
    help="save output files (*.vtu) in ascii format",
)
parser.add_argument(
    "-d",
    "--dim",
    type=int,
    default=0,
    help="spatial dimension (2 or 3), trying automatic detection, if not given",
)
parser.add_argument(
    "-o",
    "--output",
    default="",
    help="basename of output files; if not given, then it defaults to basename of inputfile",
)
parser.add_argument(
    "-z",
    "--delz",
    action="store_true",
    help="deleting z-coordinate, for 2D-meshes with z=0 (dimension must be dim=2)",
)

args = parser.parse_args()

# check if input file exists and is in gmsh-format
if os.path.isfile(args.filename):
    filename_without_extension = os.path.splitext(args.filename)[0]
    file_extension = os.path.splitext(args.filename)[1]
    if file_extension != ".msh":
        warnings.warn("Warning, input file seems not to be in gmsh-format (*.msh)")
else:
    warnings.warn("No input file (mesh) found.")
    raise FileNotFoundError

# derive output filenames
if args.output == "":  # no parameter given, use same basename as input file
    output_basename = filename_without_extension
else:
    output_basename = args.output

# read in mesh (be aware of shallow copies, i.e. by reference)
mesh = meshio.read(args.filename)
points, point_data = mesh.points, mesh.point_data
cells, cells_dict, cell_data, cell_data_dict = (
    mesh.cells,
    mesh.cells_dict,
    mesh.cell_data,
    mesh.cell_data_dict,
)
field_data = mesh.field_data
number_of_original_points = len(points)
existing_cell_types=set(mesh.cells_dict.keys())   # elements found in mesh

print("Original mesh (read)")
print_info(mesh)
# write original mesh (only file conversion)
# original_mesh=meshio.Mesh(points=points, cells=cells, cell_data=cell_data, field_data=field_data)
meshio.write(output_basename + "_original.vtu", mesh, binary=not args.ascii)


# check if element types are supported in current version of this script
all_available_cell_types=set()	# initial value
for cell_types in available_cell_types.values():	
    all_available_cell_types=all_available_cell_types.union(cell_types)
for cell_type in existing_cell_types:
    if cell_type not in all_available_cell_types:
        warnings.warn('Unsupported cell type found')

# set spatial dimension of mesh
if args.dim==0:
    # automatically detect spatial dimension of mesh
    dim=dim0	# initial value
    for test_dim, test_cell_types in available_cell_types.items():
        if len(test_cell_types.intersection(existing_cell_types)) and test_dim>dim:
            dim=test_dim

    print('Detected mesh dimension: ' + str(dim))
    print('##')
else:
    dim=args.dim	# trust the user

# delete third dimension if wanted  
if args.delz and dim==dim2:
	mesh.prune_z_0()

# boundary and domain cell types depend on dimension
if dim1<=dim and dim<=dim3:
    boundary_dim=dim-1
    domain_dim=dim
    boundary_cell_types=existing_cell_types.intersection(available_cell_types[boundary_dim])
    domain_cell_types=existing_cell_types.intersection(available_cell_types[domain_dim])
else:
    warnings.warn("Error, invalid dimension dim=" + str(dim) + "!")
    sys.exit()

# Check for existence of physical groups 
if gmsh_physical_cell_data_key in cell_data_dict:
    physical_groups_found=True

    # if user wants physical group numbering of domains beginning with zero
    id_offset = 0  # initial value, zero will not change anything
    if args.rdcd:  # prepare renumber-domain-cell-data (rdcd)
        # find minimum physical_id of domains (dim)
        id_list_domains = []
        for dataset in field_data.values():  # go through all physical groups
            if (dataset[geo_index] == domain_dim):  # only for domains, ignore lower dimensional entities
                id_list_domains.append(dataset[ph_index])  # append physical id
        if len(id_list_domains):  # if there are some domains..
            id_offset = min(id_list_domains)  # ..then find minimal physical id

else:
    print("No physical groups found.") 
    physical_groups_found=False


###############################################################################
# Extract domain mesh, note that meshio 4.3.3. offers remove_lower_dimensional_cells(), but we want to keep a uniform style for domain and subdomains. Make sure to use domain_mesh=deepcopy(mesh) in this case!
###############################################################################
domain_cells=[]	# list
if args.ogs:
    domain_cell_data_key = ogs_domain_cell_data_key
else:
    domain_cell_data_key = gmsh_physical_cell_data_key
domain_cell_data={} # dict ..
domain_cell_data[domain_cell_data_key]=[]	# .. with values: list of arrays

for domain_cell_type in domain_cell_types:
    # cells
    domain_cells_values = cells_dict[domain_cell_type]
    number_of_domain_cells=len(domain_cells_values)
    domain_cells_block=(domain_cell_type, domain_cells_values)
    domain_cells.append(domain_cells_block)

    # cell_data
    if physical_groups_found:
        if domain_cell_type in cell_data_dict[gmsh_physical_cell_data_key]: 
            domain_in_physical_group=True
        else:
            domain_in_physical_group=False
    else:
        domain_in_physical_group=False
	        
    if domain_in_physical_group:
        if args.ogs:
            domain_cell_data_values = cell_data_dict[gmsh_physical_cell_data_key][domain_cell_type]
            domain_cell_data_values = numpy.int32(domain_cell_data_values - id_offset)  # ogs needs MaterialIDs as int32, possibly beginning with zero (by id_offset)
        else:
            domain_cell_data_values = cell_data_dict[gmsh_physical_cell_data_key][domain_cell_type]
    else:
        domain_cell_data_values = numpy.zeros((number_of_domain_cells), dtype=int)
        print("Some domain cells are not in a physical group, their PhysicalID/MaterialID is set to zero.")
    domain_cell_data[domain_cell_data_key].append(domain_cell_data_values)

domain_mesh = meshio.Mesh( points=points, cells=domain_cells, cell_data=domain_cell_data)
# domain_mesh.prune()	# for older meshio version (4.0.16)
domain_mesh.remove_orphaned_nodes()


if len(domain_mesh.points) == number_of_original_points: 
    meshio.write( output_basename + "_domain.vtu", domain_mesh, binary=not args.ascii )
    print("Domain mesh (written)")
    print_info(domain_mesh)

    # store node numbers for use as bulk_node_id (point_data)
    domain_mesh_node_numbers = numpy.arange(number_of_original_points)
    
    # prepare data needed for bulk_elem_id (cell_data)
    cell_start_index=0    # node connectivity for a mixed mesh (double check element numbering!), needed with and without ogs option to identify boundary cells
    node_connectivity = [ set() for _ in range(number_of_original_points) ]  # initialize list of sets
    # make a list for each type of domain cells
    for cell_block in domain_cells:
        block_node_connectivity = cells_at_nodes(cell_block[1], number_of_original_points, cell_start_index)  # later used for boundary mesh and submeshes
        cell_start_index+=len(cell_block[1]) # assume consective cell numbering (as it is written to vtu)
        # add connectivities of current cell type to entries (sets) of total connectivity (list of sets)
        for total_list, block_list in zip(node_connectivity, block_node_connectivity):
            total_list.update(block_list)

else:
    warnings.warn( "There are nodes outside domain, this may lead to ambiguities, no domain-mesh written.")# 
    if args.ogs:
        sys.exit()	# cannot associate bulk data without node_connectivity 


###############################################################################
# Extract boundary mesh
###############################################################################

# points, process full list (all points), later trimmed according to cell selection, deep copy needed because removed_orphaned_nodes() operates on shallow copy of point_data
if args.ogs:	
    boundary_point_data={}	# dict
    boundary_point_data[ogs_point_data_key]= numpy.uint64(domain_mesh_node_numbers) 
else:
    boundary_point_data = {key: value[:] for key, value in point_data.items()}  # deep copy

# cells and cell data
boundary_cells=[]	# list
boundary_cell_data={}	# dict
if args.ogs:
    boundary_cell_data_key = ogs_boundary_cell_data_key
else:
    boundary_cell_data_key = gmsh_physical_cell_data_key
boundary_cell_data[boundary_cell_data_key]=[]	# list

for boundary_cell_type in boundary_cell_types:
        
    # cells 
    boundary_cells_values = cells_dict[boundary_cell_type] 	# preliminary, as there may be cells of boundary dimension inside domain (i.e. which are no boundary cells)
    connected_cells, connected_cells_count  = numpy.uint64( connected_domain_cells(boundary_cells_values, node_connectivity) )
    boundary_index = connected_cells_count == 1	# a boundary cell is connected with exactly one domain cell
    boundary_cells_values=boundary_cells_values[boundary_index]  # final boundary cells
    boundary_cells_block = (boundary_cell_type, boundary_cells_values)
    boundary_cells.append(boundary_cells_block)
    if not boundary_index.all():
        print("For information, there are cells of boundary dimension not on the boundary (e.g. inside domain).")

    # cell_data
    if physical_groups_found:
        if boundary_cell_type in cell_data_dict[gmsh_physical_cell_data_key]: 
            boundary_in_physical_group=True
        else:
            boundary_in_physical_group=False
    else:
        boundary_in_physical_group=False

    if args.ogs:
        boundary_cell_data_values = connected_cells[boundary_index]
    else:
        if boundary_in_physical_group:
            boundary_cell_data_values = cell_data_dict[gmsh_physical_cell_data_key][boundary_cell_type] 
        else:   
            number_of_boundary_cells=len(boundary_cells_values) # cells of specific type
            boundary_cell_data_values = numpy.zeros((number_of_boundary_cells), dtype=int)
            print("Some boundary cells are not in a physical group, their PhysicalID is set to zero.")
    boundary_cell_data[boundary_cell_data_key].append(boundary_cell_data_values)

boundary_mesh = meshio.Mesh( points=points, point_data=boundary_point_data, cells=boundary_cells, cell_data=boundary_cell_data )
boundary_mesh.remove_orphaned_nodes()
meshio.write( output_basename + "_boundary.vtu", boundary_mesh, binary=not args.ascii )
print("Boundary mesh (written)")
print_info(boundary_mesh)


###############################################################################
# Now we want to extract subdomains given by physical groups in gmsh
# name=user-defined name of physical group, data=[physical_id, geometry_id]
###############################################################################
if not physical_groups_found:
    sys.exit()

for name, data in field_data.items():
    ph_id = data[ph_index]  # selection by physical id (user defined)
    subdomain_dim = data[geo_index]  # 0 or 1 or 2 or 3
    if dim0<=subdomain_dim and subdomain_dim<=dim3:
        subdomain_cell_types=existing_cell_types.intersection(available_cell_types[subdomain_dim])
    else:
        warnings.warn("Invalid dimension found in physical groups.")
        continue

    # point data
    if args.ogs:	
        subdomain_point_data={}	# dict
        subdomain_point_data[ogs_point_data_key]= numpy.uint64(domain_mesh_node_numbers) 
    else:
        subdomain_point_data = {key: value[:] for key, value in point_data.items()}  # deep copy

    # cells, cell_data
    subdomain_cells=[]		# list
    subdomain_cell_data={}	# dict
    if args.ogs:
        if subdomain_dim==domain_dim:
            subdomain_cell_data_key = ogs_domain_cell_data_key
        elif subdomain_dim==boundary_dim:
            subdomain_cell_data_key = ogs_boundary_cell_data_key
        else:
            subdomain_cell_data_key = gmsh_physical_cell_data_key	# use gmsh, as the requirements from OGS
    else:
        subdomain_cell_data_key = gmsh_physical_cell_data_key	# same for all dimensions
    subdomain_cell_data[subdomain_cell_data_key]=[]	# list
    subdomain_cell_data_trouble=False	# flag to indicate invalid bulk_element_ids, then no cell data will be written    

    for cell_type in subdomain_cell_types:

        # cells
        selection_index = cell_data_dict[gmsh_physical_cell_data_key][cell_type] == ph_id
        selection_cells_values = cells_dict[cell_type][selection_index]
        if len(selection_cells_values):  # if there are some data
            selection_cells_block=(cell_type, selection_cells_values)
            subdomain_cells.append(selection_cells_block)

            # cell data
            if args.ogs:
            
                if subdomain_dim == boundary_dim: 
                    connected_cells, connected_cells_count = connected_domain_cells(selection_cells_values, node_connectivity)
                    boundary_index = connected_cells_count == 1 	# a boundary cell is connected with one domain cell, needed to write bulk_elem_id
                    selection_cell_data_values = numpy.uint64(connected_cells)
                    if not boundary_index.all():  
                        print("In physical group " + name + " are bulk_elem_ids not uniquely defined, e.g. for cells of boundary dimension inside the domain, and thus not written. If bulk_elem_ids should be written for a physical group, then make sure all its cells of boundary dimension are located at the boundary.")
                        subdomain_cell_data_trouble=True
                elif subdomain_dim == domain_dim: 
                    selection_cell_data_values = numpy.int32( cell_data_dict[gmsh_physical_cell_data_key][cell_type][selection_index] -id_offset )
                
                else:  # any cells of lower dimension than boundary
                    selection_cell_data_values = numpy.int32( cell_data_dict[gmsh_physical_cell_data_key][cell_type][selection_index] )

            else:
                selection_cell_data_values = cell_data_dict[gmsh_physical_cell_data_key][cell_type][selection_index]
                
            subdomain_cell_data[subdomain_cell_data_key].append(selection_cell_data_values)
    
    outputfilename = output_basename + "_physical_group_" + name + ".vtu"
    if subdomain_cell_data_trouble:
        submesh = meshio.Mesh(points=points, point_data=subdomain_point_data, cells=subdomain_cells) 	# do not write invalid cell_data 
    else:
        submesh = meshio.Mesh(points=points, point_data=subdomain_point_data, cells=subdomain_cells, cell_data=subdomain_cell_data)  

    submesh.remove_orphaned_nodes() # submesh.prune() for meshio_version 4.0.16
    outputfilename = output_basename + "_physical_group_" + name + ".vtu"
    meshio.write(outputfilename, submesh, binary=not args.ascii)
    print("Submesh " + name + " (written)")
    print_info(submesh)

