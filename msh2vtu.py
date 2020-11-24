#!/usr/bin/env python
# Author: Dominik Kern (TU Bergakademie Freiberg)
import meshio # https://github.com/nschloe/meshio
import os
import sys
import numpy
import argparse
import warnings

tested_meshio_version = "4.3.3"
tested_gmsh_version = "4.4.1"

if meshio.__version__ < tested_meshio_version:
    warnings.warn(
        "Warning, out-dated Meshio version. In case of errors watch for commented code fragments from previous versions in this script (msh2vtu)."
    )
elif meshio.__version__ > tested_meshio_version:
    print(
        "Newer version of Meshio than supported. Backward compatibility may be missing!"
    )


# 	auxiliary function needed for meshio version 4.0.16
# def line_mesh_prune(points, input_cells):  # remove orphaned points
#    # "old" means from the input mesh and "new" the mesh of connected points only
#    original_shape = input_cells.shape  # for reconstrution after flatten
#    old_points = input_cells.flatten()  # 1d-array needed
#    unique_points, unique_inverse = numpy.unique(old_points, return_inverse=True)
#    new_points = points[unique_points]  # extract only used nodes
#    new_cells = unique_inverse.reshape(original_shape)  # update cell connectivity
#    return new_points, new_cells


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
def domain_cells(cells_array, node_connectivity):
    boundary_cell_data_list = []
    for cell in cells_array:
        node1 = cell[0]
        node2 = cell[1]
        common_domain_cell = [
            a_cell
            for a_cell in node_connectivity[node1]
            if a_cell in node_connectivity[node2]
        ]  # find intersection of lists
        if len(common_domain_cell) == 1:
            boundary_cell_data_list.append(
                common_domain_cell
            )  # to be stored as cell data
        else:
            warnings.warn(
                "Boundary cell does not or not uniquely belong to a domain cell"
            )
    return numpy.array(boundary_cell_data_list)


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
    "--z-del",
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

# write original mesh (only file conversion)
# original_mesh=meshio.Mesh(points=points, cells=cells, cell_data=cell_data, field_data=field_data)
meshio.write(output_basename + "_original.vtu", mesh, binary=not args.ascii)


# some variable declarations
ph_index = 0  # to access physical id in field data
geo_index = 1  # to access geometrical id in field data
line_id = 1  # geometry type
triangle_id = 2  # geometry type
gmshdict = {line_id: "line", triangle_id: "triangle"}  # gmsh convention
gmsh_cell_physical = "gmsh:physical"
gmsh_point = "gmsh:dim_tags"
ogs_domain_cell = "MaterialIDs"
ogs_boundary_point = "bulk_node_ids"
ogs_boundary_cell = "bulk_elem_ids"
number_of_original_points = len(points)
domain_cell_type = gmshdict[triangle_id]
boundary_cell_type = gmshdict[line_id]

# if user wants physical group numbering of domains beginning with zero
id_offset = 0  # initial value, zero will not change anything
if args.rdcd:  # prepare renumber-domain-cell-data (rdcd)
    # find minimum physical_id of domains (triangles)
    id_list_domains = []
    for dataset in field_data.values():  # go through all physical groups
        if (
            dataset[geo_index] == triangle_id
        ):  # only for domains (triangles), ignore boundaries (lines)
            id_list_domains.append(dataset[ph_index])  # append physical id
    if len(id_list_domains):  # if there are some domains..
        id_offset = min(id_list_domains)  # ..then find minimal physical id


# Extract domain mesh, note that meshio 4.3.3. offers remove_lower_dimensional_cells(), but we want to keep a uniform style for domain and subdomains. Make sure to use domain_mesh=deepcopy(mesh) in this case!
if gmsh_cell_physical not in cell_data_dict:
    warnings.warn("Warning, no physical groups found.")
    sys.exit()

if domain_cell_type in cells_dict:
    domain_cells_array = cells_dict[domain_cell_type]
    domain_cell_data_array = cell_data_dict[gmsh_cell_physical][domain_cell_type]
    if args.ogs:
        domain_cell_data_string = ogs_domain_cell
        domain_cell_data_array = numpy.int32(
            domain_cell_data_array
        )  # ogs needs MaterialIDs as int32
    else:
        domain_cell_data_string = gmsh_cell_physical
    
    if len(domain_cells_array) and len(domain_cells_array) == len(
        domain_cell_data_array
    ):  # only if there are (consistent) data to write
        domain_mesh = meshio.Mesh(
            points=points,
            cells=[(domain_cell_type, domain_cells_array)],
            cell_data={domain_cell_data_string: [domain_cell_data_array - id_offset]},
        )
        # domain_mesh.prune()	# for older meshio version (4.0.16)
        domain_mesh.remove_orphaned_nodes()
        if len(domain_mesh.points) == number_of_original_points:
            meshio.write(
                output_basename + "_domain.vtu", domain_mesh, binary=not args.ascii
            )
        else:
            print(
                "There are nodes outside domain, this may lead to ambiguities, no domain-mesh written."
            )
    else:
        print("Inconsistent domain-cells, no domain-mesh written.")
else:
    print("No domain-cells found")

# Extract boundary mesh
if boundary_cell_type in cells_dict:
    boundary_cells_array = cells_dict[boundary_cell_type]
    if args.ogs:
        node_connectivity = cells_at_nodes(domain_cells_array, number_of_original_points)
        domain_mesh_node_numbers = numpy.arange(number_of_original_points)
        boundary_point_data_string = ogs_boundary_point
        boundary_point_data_array = numpy.uint64(domain_mesh_node_numbers)
        boundary_cell_data_string = ogs_boundary_cell
        boundary_cell_data_array = numpy.uint64(
            domain_cells(boundary_cells_array, node_connectivity)
        )
    else:
        boundary_point_data_string = gmsh_point
        boundary_point_data_array = point_data[gmsh_point]
        boundary_cell_data_string = gmsh_cell_physical
        boundary_cell_data_array = cell_data_dict[gmsh_cell_physical][boundary_cell_type]
    
    if len(boundary_cells_array) and len(boundary_cells_array) == len(
        boundary_cell_data_array
    ):
        boundary_mesh = meshio.Mesh(
            points=points,
            point_data={boundary_point_data_string: boundary_point_data_array},
            cells=[(boundary_cell_type, boundary_cells_array)],
            cell_data={boundary_cell_data_string: [boundary_cell_data_array]},
        )
        boundary_mesh.remove_orphaned_nodes()
        meshio.write(
            output_basename + "_boundary.vtu", boundary_mesh, binary=not args.ascii
        )
    else:
        print("Inconsistent boundary-cells, no boundary-mesh written.")
else:
    print("No boundary-cells found")


# Now we want to extract subdomains given by physical groups in gmsh
# name=user-defined name of physical group, data=[physical_id, geometry_id]
for name, data in field_data.items():

    ph_id = data[ph_index]  # selection by physical id (user defined)
    geo_id = data[geo_index]  # 1 or 2 ..
    cell_type = gmshdict[geo_id]  # .. 'line' or 'triangle'
    selection_index = cell_data_dict[gmsh_cell_physical][cell_type] == ph_id
    selection_cells_array = cells_dict[cell_type][selection_index]
    if len(selection_cells_array):  # if there are some data
        if data[geo_index] == line_id:  # boundary
            if args.ogs:
                selection_point_data_string = ogs_boundary_point
                selection_point_data_array = numpy.uint64(
                    domain_mesh_node_numbers  # all points, will be trimmed later
                )
                selection_cell_data_string = ogs_boundary_cell
                selection_cell_data_array = numpy.uint64(
                    domain_cells(selection_cells_array, node_connectivity)
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
                cells=[(boundary_cell_type, selection_cells_array)],
                cell_data={selection_cell_data_string: [selection_cell_data_array]},
            )
            submesh.remove_orphaned_nodes()  # trim mesh
            # workaround for meshio version 4.0.16
            # new_points, new_cells = line_mesh_prune(points, selected_cells)
            # submesh = meshio.Mesh(points=new_points, cells=[(cell_type, new_cells)])
        elif data[geo_index] == triangle_id:  # domain
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
                cells=[(domain_cell_type, selection_cells_array)],
                cell_data={selection_cell_data_string: [selection_cell_data_array]},
            )  # point_data not needed
            # submesh.prune()	# for meshio_version 4.0.16
            submesh.remove_orphaned_nodes()
        else:# TODO handle vertex
            print("Unknown geometry id encountered, empty submesh.")
            submesh = meshio.Mesh(points=[], cells=[])
        outputfilename = output_basename + "_physical_group_" + name + ".vtu"
        meshio.write(outputfilename, submesh, binary=not args.ascii)
    else:
        print("No cells found for physical group " + name + ", no submesh written.")
