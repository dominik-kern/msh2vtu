#!/usr/bin/env python
# Author: Dominik Kern (TU Bergakademie Freiberg)
# Python version 3.8.5
import meshio  # https://github.com/nschloe/meshio
import os
import sys
import numpy
import argparse
import warnings


# 	auxiliary function needed for meshio version 4.0.16
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
    for cell_type, cell_array in mesh.cells_dict.items():
        cell_info += str(len(cell_array)) + " " + cell_type + ", "
    print(cell_info[0:-2], end='; ')
    print("point_data=" + str(list(mesh.point_data)), end='; ')
    print("cell_data=" + str(list(mesh.cell_data)), end='; ')
    print("cell_sets=" + str(list(mesh.cell_sets)))
    print("##")


# function to create node connectivity list
def cells_at_nodes(cells, node_count):
    empty_list = []
    node_connectivity = [
        empty_list[:] for _ in range(node_count)
    ]  # initialize list of lists
    cell_number = 0
    for cell in cells:
        for node in cell:
            node_connectivity[node].append(cell_number)
        cell_number += 1
    return node_connectivity


# function to find out to which domain elements a line element belongs
def connected_domain_cells(cells_array, node_connectivity, dim):
    boundary_cell_data_list = []
    for cell in cells_array:
        node1 = cell[0]
        node2 = cell[1]
        if dim==2:
            common_domain_cell = [
                a_cell
                for a_cell in node_connectivity[node1]
                if a_cell in node_connectivity[node2]
            ]  # find intersection of 2 lists
        elif dim==3:
            node3 = cell[2]        
            common_domain_cell = [
                a_cell
                for a_cell in node_connectivity[node1]
                if (a_cell in node_connectivity[node2]) and (a_cell in node_connectivity[node3])
            ] 	# find intersection of 3 lists
        else:
            warnings.warn("domain_cells: invalid dimension")
 
        if len(common_domain_cell) == 1:
            boundary_cell_data_list.append(
                common_domain_cell
            )  # to be stored as cell data
        else:
            warnings.warn(
                "domain cells: boundary cell does not or not uniquely belong to a domain cell"
            )

    return numpy.array(boundary_cell_data_list)


# some variable declarations
ph_index = 0	# to access physical group id in field data
geo_index = 1	# to access geometrical dimension in field data
vertex_id = 0	# dimension
line_id = 1	# dimension
triangle_id = 2	# dimension
tetra_id = 3 	# dimension
gmshdict = {
    vertex_id: "vertex",
    line_id: "line",
    triangle_id: "triangle",
    tetra_id: "tetra"
}  # gmsh convention
gmsh_cell_physical = "gmsh:physical"
gmsh_point = "gmsh:dim_tags"
ogs_domain_cell = "MaterialIDs"
ogs_boundary_point = "bulk_node_ids"
ogs_boundary_cell = "bulk_elem_ids"

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
    description="Prepares a Gmsh-mesh for use in OGS by extracting domain-, boundary- and physical group-submeshes and saves them in vtu-format. Note that all mesh entities must belong to physical groups.",
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
    help="deleting z-coordinate, for 2D-meshes with z=0",
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

if args.dim==0:
    # automatically detect spatial dimension of mesh
    element_types=mesh.cells_dict.keys()
    if gmshdict[tetra_id] in element_types:
        dim=3
    elif gmshdict[triangle_id] in element_types:
        dim=2
    elif gmshdict[line_id] in element_types:
        dim=1
    else:
        warnings.warn("Neither 3D nor lower dimensional elements (2D, 1D) found.")
        dim=0
    print('Detected mesh dimension: ' + str(dim))
else:
    dim=args.dim

print("Original mesh (read)")
print_info(mesh)
# write original mesh (only file conversion)
# original_mesh=meshio.Mesh(points=points, cells=cells, cell_data=cell_data, field_data=field_data)
meshio.write(output_basename + "_original.vtu", mesh, binary=not args.ascii)

if args.delz:
	mesh.prune_z_0()


# check if element types are supported in current version of this script
for element_type in mesh.cells_dict.keys():
    if not element_type in gmshdict.values():
        warnings.warn('Unsupported element type found')

# boundary and domain cell types depend on dimension
if dim==1:
    boundary_id=vertex_id
    domain_id=line_id
    boundary_cell_type=gmshdict[vertex_id]
    domain_cell_type=gmshdict[line_id]
elif dim==2:
    boundary_id=line_id	# dimension
    domain_id=triangle_id	# dimension
    boundary_cell_type=gmshdict[line_id]
    domain_cell_type=gmshdict[triangle_id]
elif dim==3:
    boundary_id=triangle_id	# dimension
    domain_id=tetra_id		# dimension
    boundary_cell_type=gmshdict[triangle_id]
    domain_cell_type=gmshdict[tetra_id]
else:
    warnings.warn("Error, invalid dimension")
    sys.exit()

# Check for existence of physical groups 
if gmsh_cell_physical in cell_data_dict:
    physical_groups_found=True

    # check whether there are physical groups for domain and boundary cells
    if domain_cell_type in cell_data_dict[gmsh_cell_physical]:
        domain_physical_group_found=True
        # if user wants physical group numbering of domains beginning with zero
        id_offset = 0  # initial value, zero will not change anything
        if args.rdcd:  # prepare renumber-domain-cell-data (rdcd)
            # find minimum physical_id of domains (triangles)
            id_list_domains = []
            for dataset in field_data.values():  # go through all physical groups
                if (dataset[geo_index] == triangle_id):  # only for domains (triangles), ignore boundaries (lines) and vertices (points)
                    id_list_domains.append(dataset[ph_index])  # append physical id
            if len(id_list_domains):  # if there are some domains..
                id_offset = min(id_list_domains)  # ..then find minimal physical id

    else:
        domain_physical_group_found=False

    if boundary_cell_type in cell_data_dict[gmsh_cell_physical]:
        boundary_physical_group_found=True
    else:
        boundary_physical_group_found=False 

else:
    print("No physical groups found.")
    physical_groups_found=False
    domain_physical_group_found=False
    boundary_physical_group_found=False


###############################################################################
# Extract domain mesh, note that meshio 4.3.3. offers remove_lower_dimensional_cells(), but we want to keep a uniform style for domain and subdomains. Make sure to use domain_mesh=deepcopy(mesh) in this case!
if domain_cell_type in cells_dict:
    domain_cells_array = cells_dict[domain_cell_type]
    number_of_domain_cells=len(domain_cells_array)
    domain_cells=(domain_cell_type, domain_cells_array)
	        
    if args.ogs:
        domain_cell_data_string = ogs_domain_cell
        if domain_physical_group_found:  
            domain_cell_data_array = cell_data_dict[gmsh_cell_physical][domain_cell_type]
            domain_cell_data_array = numpy.int32(domain_cell_data_array - id_offset)  # ogs needs MaterialIDs as int32, possibly beginning with zero (by id_offset)
        else:
            domain_cell_data_array = numpy.zeros((number_of_domain_cells), dtype=int)
        
        domain_cell_data = {domain_cell_data_string: [domain_cell_data_array]}
        node_connectivity = cells_at_nodes(domain_cells_array, number_of_original_points)  # later used for boundary mesh and submeshes
    else:
        domain_cell_data_string = gmsh_cell_physical
        if domain_physical_group_found:
            domain_cell_data_array = cell_data_dict[gmsh_cell_physical][domain_cell_type]
            domain_cell_data = {domain_cell_data_string: [domain_cell_data_array]}
        else:
            domain_cell_data = {}
    
    domain_mesh = meshio.Mesh( points=points, cells=[domain_cells], cell_data=domain_cell_data)
    # domain_mesh.prune()	# for older meshio version (4.0.16)
    domain_mesh.remove_orphaned_nodes()
    if len(domain_mesh.points) == number_of_original_points: 
        meshio.write( output_basename + "_domain.vtu", domain_mesh, binary=not args.ascii)
        print("Domain mesh (written)")
        print_info(domain_mesh)
    else:
        print( "There are nodes outside domain, this may lead to ambiguities, no domain-mesh written.")
else:
    print("No domain-cells found, no domain-mesh written.")
    node_connectivity=[]


###############################################################################
# Extract boundary mesh
if boundary_cell_type in cells_dict:
    boundary_cells_array = cells_dict[boundary_cell_type]
    boundary_cells = (boundary_cell_type, boundary_cells_array)

    if args.ogs:
        domain_mesh_node_numbers = numpy.arange(number_of_original_points)
        boundary_point_data_string = ogs_boundary_point
        boundary_point_data_array = numpy.uint64(domain_mesh_node_numbers)
    else:
        boundary_point_data_string = gmsh_point
        boundary_point_data_array = point_data[gmsh_point]
    boundary_point_data = {boundary_point_data_string: boundary_point_data_array} 

    if args.ogs:
        boundary_cell_data_string = ogs_boundary_cell
        boundary_cell_data_array = numpy.uint64( connected_domain_cells(boundary_cells_array, node_connectivity, dim) )
        boundary_cell_data = {boundary_cell_data_string: [boundary_cell_data_array]}
    else:
        if boundary_physical_group_found: 
            boundary_cell_data_string = gmsh_cell_physical
            boundary_cell_data_array = cell_data_dict[gmsh_cell_physical][boundary_cell_type] 
            boundary_cell_data = {boundary_cell_data_string: [boundary_cell_data_array]}
        else:
            boundary_cell_data = {}

    boundary_mesh = meshio.Mesh( points=points, point_data=boundary_point_data, cells=[boundary_cells], cell_data=boundary_cell_data )

    boundary_mesh.remove_orphaned_nodes()
    meshio.write( output_basename + "_boundary.vtu", boundary_mesh, binary=not args.ascii)
    print("Boundary mesh (written)")
    print_info(boundary_mesh)

else:
    print("No boundary-cells found")


###############################################################################
# Now we want to extract subdomains given by physical groups in gmsh
# name=user-defined name of physical group, data=[physical_id, geometry_id]
if not physical_groups_found:
    sys.exit()

for name, data in field_data.items():
    ph_id = data[ph_index]  # selection by physical id (user defined)
    geo_id = data[geo_index]  # 0 or 1 or 2 or 3
    cell_type = gmshdict[geo_id]  # 'vertex' or 'line' or 'triangle' or 'tetra'
    selection_index = cell_data_dict[gmsh_cell_physical][cell_type] == ph_id
    selection_cells_array = cells_dict[cell_type][selection_index]
    if len(selection_cells_array):  # if there are some data
        if data[geo_index] == boundary_id: 
            if args.ogs:
                selection_point_data_string = ogs_boundary_point
                selection_point_data_array = numpy.uint64(
                    domain_mesh_node_numbers  # all points, will be trimmed later
                )
                selection_cell_data_string = ogs_boundary_cell
                selection_cell_data_array = numpy.uint64(
                    connected_domain_cells(selection_cells_array, node_connectivity, dim)
                )
            else:
                selection_point_data_string = gmsh_point
                selection_point_data_array = point_data[
                    gmsh_point
                ]  # all points, will be trimmed later
                selection_cell_data_string = gmsh_cell_physical
                selection_cell_data_array = cell_data_dict[gmsh_cell_physical][
                    cell_type
                ][selection_index]

            submesh = meshio.Mesh(
                points=points,
                point_data={selection_point_data_string: selection_point_data_array},
                cells=[(cell_type, selection_cells_array)],
                cell_data={selection_cell_data_string: [selection_cell_data_array]},
            )
            submesh.remove_orphaned_nodes()  # trim mesh
            # workaround for meshio version 4.0.16
            # new_points, new_cells = line_mesh_prune(points, selected_cells)
            # submesh = meshio.Mesh(points=new_points, cells=[(cell_type, new_cells)])
        elif data[geo_index] == domain_id: 
            selection_cell_data_array = (
                cell_data_dict[gmsh_cell_physical][cell_type][selection_index]
                - id_offset
            )
            if args.ogs:
                selection_cell_data_string = ogs_domain_cell
                selection_cell_data_array = numpy.int32(selection_cell_data_array)
            else:
                selection_cell_data_string = gmsh_cell_physical

            submesh = meshio.Mesh(
                points=points,
                cells=[(cell_type, selection_cells_array)],
                cell_data={selection_cell_data_string: [selection_cell_data_array]},
            )  # point_data not needed
            # submesh.prune()	# for meshio_version 4.0.16
            submesh.remove_orphaned_nodes()
        elif 0<= data[geo_index] and data[geo_index]<=3:  
            selection_cell_data_string = gmsh_cell_physical	# keep gmsh string, since there are no requirements from OGS 
            selection_cell_data_array = cell_data_dict[gmsh_cell_physical][cell_type][
                selection_index
            ]
            if args.ogs:
                selection_cell_data_array = numpy.int32(selection_cell_data_array)

            submesh = meshio.Mesh(
                points=points,
                cells=[(cell_type, selection_cells_array)],
                cell_data={selection_cell_data_string: [selection_cell_data_array]},
            )  # point_data not needed
            # submesh.prune()   # for meshio_version 4.0.16
            submesh.remove_orphaned_nodes()
        else:
            print("Unknown geometry id encountered, empty submesh.")
            submesh = meshio.Mesh(points=[], cells=[])
        outputfilename = output_basename + "_physical_group_" + name + ".vtu"
        meshio.write(outputfilename, submesh, binary=not args.ascii)
        print("Submesh " + name + " (written)")
        print_info(submesh)

    else:
        print("No cells found for physical group " + name + ", no submesh written.")

